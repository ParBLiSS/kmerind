/*
 * file.hpp
 *
 *  Created on: Feb 28, 2016
 *      Author: tpan
 */

#ifndef FILE_HPP_
#define FILE_HPP_

#include <string>
#include <cstring>      // memcpy, strerror

#include <ios>          // ios_base::failure
#include <iostream>     // ios_base::failure
#include <unistd.h>     // sysconf
#include <sys/mman.h>   // mmap
#include <fcntl.h>      // for open64 and close
#include <sstream>      // stringstream
#include <exception>    // std exception

#define USE_MPI

#if defined(USE_MPI)
#include <mpi.h>
#endif

// mxx
#include <mxx/collective.hpp>
#include <mxx/shift.hpp>

#include <io/file_loader.hpp>
#include <io/fastq_loader.hpp>
#include <io/fasta_loader.hpp>
#include <partition/range.hpp>



namespace bliss {

namespace io {

/**
 * @brief loaded file data.
 */
struct file_data {
  using iterator = unsigned char*;

	// type of ranges
	using range_type = ::bliss::partition::range<size_t>;

	// range from which the data came
	range_type parent_range_bytes;

	// range loaded in memory
	range_type in_mem_range_bytes;

	// valid range for this
	range_type valid_range_bytes;

	// storage for actual data
	std::vector<unsigned char> data;

	unsigned char * begin() {
		return data.data() + valid_range_bytes.start - in_mem_range_bytes.start;
	}
	unsigned char * end() {
		return data.data() + valid_range_bytes.end - in_mem_range_bytes.start;
	}
	range_type getRange() {
		return valid_range_bytes;
	}

};



/// base class to wrap a file object with structured data or bytes
class base_file {

public:
	/// flag to indicate read
	/// flag to indicate write

	using range_type = ::bliss::partition::range<size_t>;

protected:
	/// name of file.
	std::string filename;

	/// size of file in bytes
	range_type file_range_bytes;

	/// virtual function for computing the size of a file.
	size_t get_file_size() {

		// get the file size.
		struct stat64 filestat;
		// get the file state
		int ret = stat64(filename.c_str(), &filestat);

		// handle any error
		if (ret < 0) {
			::std::stringstream ss;
			int myerr = errno;
			ss << "ERROR in file size calc: ["  << this->filename << "] error " << myerr << ": " << strerror(myerr);
			throw new ::std::ios_base::failure(ss.str());
		}

		return static_cast<size_t>(filestat.st_size);
	}


public:

	/**
	 * @brief  bulk load all the data and return it in a newly constructed vector
	 * @param range_bytes	  range to read, in bytes
	 * @return				  vector containing data as bytes.
	 */
	virtual std::vector<unsigned char> read_range(range_type const & range_bytes) {
		::std::vector<unsigned char> out;
		this->read_range(out, range_bytes);
		return out;
	}

	/**
	 * @brief  bulk load all the data and return it in a newly constructed vector.  reuse vector
	 * @note	virtual so different file reading mechanisms can be defined
	 * @param range_bytes	range to read, in bytes
	 * @param output		vector containing data as bytes.
	 */
	virtual range_type read_range(std::vector<unsigned char> & output, range_type const & range_bytes) = 0;


	/// flag to indicate read
	/// flag to indicate write

	/**
	 * @brief initializes a file for reading/writing
	 * @note 	do not call get_file_size here because
	 * 			parallel_file inherits this class and defines a parallel get_file_size.
	 * 			we want to call this constructor from subclass.
	 *
	 * @param _filename 	name of file to open
	 */
	base_file(std::string const & _filename) :
		filename(_filename), file_range_bytes(0, 0) {};

	/// default destructor
	virtual ~base_file() {};


	/**
	 * @brief  read the whole file
	 * @return 				file_data object containing data and various ranges.
	 */
	::bliss::io::file_data read_file() {
		file_data out;
		read_file(out);
		return out;
	}

	/**
	 * @brief  read the whole file.  reuse allocated file_data object.
	 * @note	virtual so that parallel readers can compute range to read.
	 * @param output 		file_data object containing data and various ranges.
	 */
	virtual void read_file(::bliss::io::file_data & output) {

		output.in_mem_range_bytes =
			output.valid_range_bytes =
				this->read_range(output.data, file_range_bytes);

		output.parent_range_bytes = file_range_bytes;
	}

	/// get file size
	size_t size() { return file_range_bytes.end; };

	/// get file name
	::std::string const & get_filename() const { return filename; };
};



/**
 * @brief  subclass for memmapped file
 */
class mmap_file : public ::bliss::io::base_file {

	// share open and close with file class.
protected:
	/// BASE type
	using BASE = ::bliss::io::base_file;

	/// file descriptor.  -1 means not open.
	int fd;

	/// mapped data's pointer
	unsigned char * mapped_data;

	/// range of the mapped data.  in units of bytes
	typename BASE::range_type mapped_range_bytes;

	/// the page size for the os.   in units of bytes
	const size_t page_size;


	/// opens a file. side effect computes size of the file.
	void open_file() {
		// first close file.
		close_file();

		// open the file and get a handle.
		fd = open64(this->filename.c_str(), O_RDONLY);
		if (fd == -1)
		{
			this->file_range_bytes.end = 0;

			// if open failed, throw exception.
			::std::stringstream ss;
			int myerr = errno;
			ss << "ERROR in file open: ["  << this->filename << "] error " << myerr << ": " << strerror(myerr);
			throw new ::std::ios_base::failure(ss.str());
		}
	}


	/// funciton for closing a file
	void close_file() {
		if (fd >= 0) {
			close(fd);
			fd = -1;
			this->file_range_bytes.end = 0;
		}
	}
public:

	// this is needed to prevent overload name hiding.  see http://stackoverflow.com/questions/888235/overriding-a-bases-overloaded-function-in-c/888337#888337
	using BASE::read_range;

	/**
	 *  bulk load all the data and return it in pre-allocated vector
	 *  @param output   vector where the results are to be stored.
	 */
	virtual typename BASE::range_type read_range(std::vector<unsigned char> & output,
			typename BASE::range_type const & range_bytes) {

		// map
		map(range_bytes);

		if (mapped_data == nullptr) {
			std::cout << "WARNING: mapped data is null" << std::endl;
			output.clear();
		}

		// ensure the portion to copy is within the mapped region.
		typename BASE::range_type target =
				BASE::range_type::intersect(mapped_range_bytes, range_bytes);

		if (target.size() == 0) {
			// print error through exception.
			std::cout << "WARNING: read_range: requested " << range_bytes << " not in mapped " << mapped_range_bytes;
			output.clear();
			return target;
		}

		// resize output's capacity
		if (output.capacity() < target.size()) output.resize(target.size());

		// copy the data into memory.  vector is contiguous, so this is okay.
		memmove(output.data(), mapped_data + (target.start - mapped_range_bytes.start), target.size());

		// unmap
		unmap();

		return target;
	}


	unsigned char * const & get_mapped_data() const {
		return mapped_data;
	}

	range_type const & get_mapped_range() const {
		return mapped_range_bytes;
	}

	size_t const & get_page_size() const {
		return page_size;
	}


	/**
	 * @brief   map the specified portion of the file to memory.
	 * @note    AGNOSTIC of overlaps
	 * @param range_bytes    range specifying the portion of the file to map.
	 */
	void map(typename BASE::range_type const & range_bytes) {
		// first clear.
		unmap();

		// if no file
		if (this->fd == -1)
			throw ::std::ios_base::failure("ERROR: map: file is not yet open.");

		typename BASE::range_type target =
				BASE::range_type::intersect(this->file_range_bytes, range_bytes);

		if (target.size() == 0) {
			// print error through exception.
			std::cout << "WARN: map: requested " << range_bytes << " not in file " << this->file_range_bytes;

			return;
		}

		// get start and end positions that are page aligned.
	    mapped_range_bytes.start = BASE::range_type::align_to_page(target, page_size);
	    mapped_range_bytes.end = target.end;  // okay for end not to align - made 0.

	    // NOT using MAP_POPULATE. (SLOW)  no need to prefault the entire range - use read ahead from madvice.
	    // NOTE HUGETLB not supported for file mapping, only anonymous.  also, kernel has to enable it and system has to have it reserved,
	    // MAP_SHARED so that we don't have CoW (no private instance) (potential sharing between processes?)  slightly slower by 1%?
	    // MAP_NORESERVE so that swap is not allocated.
	    mapped_data = (unsigned char*)mmap64(nullptr, mapped_range_bytes.size(),
	                                 	 PROT_READ,
	                                 	 MAP_SHARED | MAP_NORESERVE, fd,
	                                 	 mapped_range_bytes.start);

	    // if mmap failed,
	    if (mapped_data == MAP_FAILED)
	    {
			// print error through exception.
			std::stringstream ss;
			int myerr = errno;
			ss << "ERROR in mmap: " << myerr << ": " << strerror(myerr);
			throw ::std::ios_base::failure(ss.str());
	    }

		// set the madvice info.  SEQUENTIAL vs RANDOM does not appear to make a difference in running time.
		int madv_result = madvise(mapped_data, mapped_range_bytes.size(),
				MADV_SEQUENTIAL | MADV_WILLNEED);
		if ( madv_result == -1 ) {
			std::stringstream ss;
			int myerr = errno;
			ss << "ERROR in madvise: " << myerr << ": " << strerror(myerr);
			throw std::ios_base::failure(ss.str());
		}
	}

	void unmap() {
		// unmap it.
		if ((mapped_data != nullptr) && (mapped_range_bytes.size() > 0)) {
			munmap(mapped_data, mapped_range_bytes.size());

			mapped_data = nullptr;
			mapped_range_bytes.start = mapped_range_bytes.end;
		}
	}

	/**
	 * initializes a file for reading/writing via memmap
	 * @param _filename 	name of file to open
	 */
	mmap_file(std::string const & _filename) :
		::bliss::io::base_file(_filename), fd(-1), mapped_data(nullptr), mapped_range_bytes(0, 0),
		page_size(sysconf(_SC_PAGE_SIZE)) {
		this->open_file();

		this->file_range_bytes.end = this->get_file_size();
	};

	/// default destructor
	virtual ~mmap_file() {
		unmap();
		this->close_file();
	};


	// this is needed to prevent overload name hiding.  see http://stackoverflow.com/questions/888235/overriding-a-bases-overloaded-function-in-c/888337#888337
	using BASE::read_file;

};

/**
 * @brief    file wrapper that uses c stdio calls.  has buffering.
 */
class stdio_file : public ::bliss::io::base_file {

protected:

	using BASE = ::bliss::io::base_file;

	// file handle;
	FILE *fp;

	/// virtual function for opening a file. side effect computes size of the file.
	void open_file() {
		this->close_file();

		fp = fopen64(this->filename.c_str(), "r");

		if (fp == nullptr) {
			this->file_range_bytes.end = 0;

			// if open failed, throw exception.
			::std::stringstream ss;
			int myerr = errno;
			ss << "ERROR in file open: ["  << this->filename << "] error " << myerr << ": " << strerror(myerr);
			throw new ::std::ios_base::failure(ss.str());

		}
	}

	/// virtual funciton for closing a file
	void close_file() {
		if (fp != nullptr) {
			fclose(fp);
			fp = nullptr;
			this->file_range_bytes.end = 0;
		}
	}
public:

	// this is needed to prevent overload name hiding.  see http://stackoverflow.com/questions/888235/overriding-a-bases-overloaded-function-in-c/888337#888337
	using BASE::read_range;

	/**
	 * @brief  bulk load all the data and return it in a newly constructed vector.  reuse vector
	 * @note	virtual so different file reading mechanisms can be defined AND parallel
	 * @param range_bytes	range to read, in bytes
	 * @param output		vector containing data as bytes.
	 */
	virtual range_type read_range(std::vector<unsigned char> & output, range_type const & range_bytes) {
		if (fp == nullptr) {
			throw std::logic_error("ERROR: read_range: file pointer is null");
		}

		// ensure the portion to copy is within the mapped region.
		typename BASE::range_type target =
				BASE::range_type::intersect(this->file_range_bytes, range_bytes);

		if (target.size() == 0) {
			// print error through exception.
			std::cout << "WARNING: read_range: requested " << range_bytes << " not in file " << file_range_bytes;
			output.clear();
			return target;
		}

		// use fseek instead lseak to preserve buffering.
		int res = fseeko64(fp, target.start, SEEK_SET);
		if ( res == -1 ) {
			std::stringstream ss;
			int myerr = errno;
			ss << "ERROR: fseeko64: file " << this->filename << " error " << myerr << ": " << strerror(myerr);
			throw std::ios_base::failure(ss.str());
		}

		// resize output's capacity
		if (output.capacity() < target.size()) output.resize(target.size());

		size_t read = fread_unlocked(output.data(), 1, target.size(), fp);

		if (read < 0) {
			std::stringstream ss;
			int myerr = errno;
			ss << "ERROR: fread: file " << this->filename << " error " << myerr << ": " << strerror(myerr);
			throw std::ios_base::failure(ss.str());
		}

		return target;
	}
	/**
	 * initializes a file for reading/writing
	 * @param _filename 	name of file to open
	 */
	stdio_file(std::string const & _filename) : ::bliss::io::base_file(_filename), fp(nullptr) {
		this->file_range_bytes.end = this->get_file_size();
		this->open_file();
	};

	/// destructor
	virtual ~stdio_file() {
		this->close_file();
	};

	// this is needed to prevent overload name hiding.  see http://stackoverflow.com/questions/888235/overriding-a-bases-overloaded-function-in-c/888337#888337
	using BASE::read_file;

};






namespace parallel {


class base_file : public ::bliss::io::base_file {
protected:
	using BASE = ::bliss::io::base_file;

	/// communicator used.  object instead of being a reference - lifetime of a src comm temp object in constructor is just that of the constructor call.
	const ::mxx::comm comm;

		/// compute file size in parallel (1 proc, then broadcast.
	size_t get_file_size() {

		size_t file_size = 0;

		// rank 0 gets the size
		if (comm.rank() == 0)
			file_size = BASE::get_file_size();

		// then broadcast.
		if (comm.size() > 1)
			MPI_Bcast(&file_size, 1, MPI_UNSIGNED_LONG, 0, comm);

		return file_size;
	}
public:

	// this is needed to prevent overload name hiding.  see http://stackoverflow.com/questions/888235/overriding-a-bases-overloaded-function-in-c/888337#888337
	using BASE::read_range;

	/**
	 * @brief constructor
	 * @note	_comm could be a temporary constructed from MPI_Comm.  in which case, the lifetime is the call of the constructor.
	 * 			to save it, we cannot just save the reference - the underlying object would go out of scope and die.
	 * 			we have to copy it.  ::mxx::comm does not have copy constructor, but has explict copy function and move constructor.  so use them.
	 * 			note that copy function is collective.  move constructor is not collective, so it could cause deadlock.
	 * @param _filename 		name of file to open
	 * @param _comm				MPI communicator to use.
	 */
	base_file(std::string const & _filename, ::mxx::comm const & _comm = ::mxx::comm()) :
		::bliss::io::base_file(_filename), comm(::std::move(_comm.copy())) {  // _comm could be a temporary constructed from MPI_Comm.
		this->file_range_bytes.end = this->get_file_size();
	};

	/// destructor
	virtual ~base_file() {};  // will call super's unmap.

	// this is needed to prevent overload name hiding.  see http://stackoverflow.com/questions/888235/overriding-a-bases-overloaded-function-in-c/888337#888337
	using BASE::read_file;

};


template <typename FileReader, typename FileParser = ::bliss::io::BaseFileParser<unsigned char*> >
class partitioned_file : public ::bliss::io::parallel::base_file {

protected:
	using BASE = ::bliss::io::parallel::base_file;

	using range_type = typename ::bliss::io::base_file::range_type;

	/// FileReader
	FileReader reader;

	/// overlap amount
	const size_t overlap;

	/// partitioner to use.
	::bliss::partition::BlockPartitioner<range_type> partitioner;


	/**
	 * @brief partitions the specified range by the number of processes in communicator
	 * @note  does not add overlap.  this is strictly for block partitioning a range.
	 * 			if overlap is needed, use the overload version.
	 * @note  assumes input parameter is in memory and inside file.
	 */
	range_type partition(range_type const & range_bytes) {
		range_type target =
				BASE::range_type::intersect(range_bytes, this->file_range_bytes);

		if (this->comm.size() == 1) {
			return target;
		}

		partitioner.configure(target, this->comm.size());

		return partitioner.getNext(this->comm.rank());
	}

	/**
	 * @brief  partition the in mem valid range by the number of processes in communicator
	 * @note	overlap is added to the end of the result ranges.  This may make the result extend
	 * 			beyong initial in_mem_range.
	 * @note	assumes that both input parameters are within file range and in memory.
	 */
	::std::pair<range_type, range_type>
	partition(range_type const & in_mem_range_bytes,
			range_type const & valid_range_bytes) {
		::std::cout << "rank " << this->comm.rank() << " partition comm_size " << this->comm.size() << ::std::endl;
		std::cout << " partition comm object at " << &(this->comm) << " for this " << this << std::endl;
        std::cout << "partition base comm is set to world ? " <<  (MPI_COMM_WORLD == comm ? "y" : "n") << std::endl;

		range_type in_mem =
				BASE::range_type::intersect(in_mem_range_bytes, this->file_range_bytes);

		// and restrict valid to in memory data.
		range_type valid =
				BASE::range_type::intersect(valid_range_bytes, in_mem);

		// single process, just return
		if (this->comm.size() == 1) {
			return ::std::make_pair(in_mem, valid);
		}
		// === multi process.

		// partition valid range
		partitioner.configure(valid, this->comm.size());
		valid = partitioner.getNext(this->comm.rank());

		// compute the in mem range.  extend by overlap
		range_type in_mem_valid = valid;
		in_mem_valid.end += overlap;

		in_mem_valid.intersect(in_mem);

		return ::std::make_pair(in_mem_valid, valid);
	}

public:
	// this is needed to prevent overload name hiding.  see http://stackoverflow.com/questions/888235/overriding-a-bases-overloaded-function-in-c/888337#888337
	using BASE::read_range;

	/**
	 * @brief  bulk load the data and return it in a newly constructed vector.  block decomposes the range.  reuse vector
	 * @note   virtual so other overloads of this function can also block decompose the range.
	 * @param range_bytes	range to read, in bytes
	 * @param output		vector containing data as bytes.
	 */
	virtual range_type read_range(std::vector<unsigned char> & output,
				range_type const & range_bytes) {
		// first get rough partition
		range_type target = partition(range_bytes);

		// then read the range via sequential version
		reader.read_range(output, target);

		return target;
	}

	/**
	 * @brief constructor
	 * @param _filename 		name of file to open
	 * @param _comm				MPI communicator to use.
	 */
	partitioned_file(std::string const & _filename, size_t const & _overlap = 0UL, ::mxx::comm const & _comm = ::mxx::comm()) :
		::bliss::io::parallel::base_file(_filename, _comm), reader(_filename), overlap(_overlap) {};

	/// destructor
	virtual ~partitioned_file() {};  // will call super's unmap.


	// this is needed to prevent overload name hiding.  see http://stackoverflow.com/questions/888235/overriding-a-bases-overloaded-function-in-c/888337#888337
	using BASE::read_file;

	/**
	 * @brief  read the whole file.  reuse allocated file_data object.
	 * @note	virtual so that parallel readers can compute range to read.
	 * @param output 		file_data object containing data and various ranges.
	 */
	virtual void read_file(::bliss::io::file_data & output) {

		range_type in_mem_partitioned;
		range_type valid_partitioned;

		::std::tie(in_mem_partitioned, valid_partitioned) =
				partition(this->file_range_bytes, this->file_range_bytes);

		// then read the range via sequential version
		output.in_mem_range_bytes = reader.read_range(output.data, in_mem_partitioned);

		output.valid_range_bytes = valid_partitioned;
		output.parent_range_bytes = file_range_bytes;
	}
};

template <typename FileReader >
class partitioned_file<FileReader, ::bliss::io::FASTQParser<unsigned char* > > :
	public ::bliss::io::parallel::base_file {

protected:
	using BASE = ::bliss::io::parallel::base_file;

	/// FileReader
	FileReader reader;

	/// overlap amount
	const size_t overlap;

	/// partitioner to use.
	::bliss::partition::BlockPartitioner<range_type> partitioner;

	/**
	 * @brief partitions the specified range by the number of processes in communicator
	 * @note  does not add overlap.  this is strictly for block partitioning a range.
	 * 			if overlap is needed, use the overload version.
	 * @note  assumes input parameter is in memory and inside file.
	 */
	range_type partition(range_type const & range_bytes) {

		// ensure that the search range is inside the file
		// restrict in_mem to within file
		range_type target =
				BASE::range_type::intersect(range_bytes, this->file_range_bytes);

		// covers whole file, so done.
		if (target.size() == 0) return target;

		// allocate a mmap file for doing the search.
		::bliss::io::mmap_file record_finder(this->filename);

		// configure the partitioner and get the local partition.
		partitioner.configure(target, this->comm.size());
		range_type hint = partitioner.getNext(this->comm.rank());

		// now define the search ranges.  for case when partition is smaller than a record,
		// we need to search more, so make the end of search range (page_size + end).
		range_type search = hint;
		// but with staggered start
		search.end = hint.end + record_finder.get_page_size();  // assume FileParser can check the previous char.

		// ensure that previous char is inspected.
		range_type load = search;
		load.start = (search.start == this->file_range_bytes.start ?
				this->file_range_bytes.start : search.start - 1);  // read extra char before


		//===  search the region for the starting position, using mmap file
		range_type found;

		// map the search region
		record_finder.map(load);

		//=== search for record starting point using FileParser.
		::bliss::io::FASTQParser<unsigned char *> parser;

		parser.init_parser(record_finder.get_mapped_data(), this->file_range_bytes,
				record_finder.get_mapped_range(), search, this->comm);

		// mark the first entry found.
		found.start = parser.find_first_record(record_finder.get_mapped_data(),
				this->file_range_bytes, record_finder.get_mapped_range(), search);

		record_finder.unmap();

		//	std::cout << "rank " << comm.rank() << " start " << final.start << std::endl;

		// each rank will find a start, which may be outside of hint.

		// now get the end from next proc that found a start.


		// not using comm split because the subcommunicator may be constructed in worse than log p time.
		// next we do reverse exclusive scan with rev communicator
		found.end = ::mxx::exscan(found.start, [](size_t const & x, size_t const & y) {
			return ::std::min(x, y);
		}, this->comm.reverse());
		if (this->comm.rank() == (this->comm.size() - 1)) found.end = target.end;   // last process.

		// finally, check final_start again, and if start is pointing to end of file, let it point to the current range's end.
		//if (found.start == target.end) found.start = found.end;
//		std::cout << "rank " << this->comm.rank() << " found = " << found <<  std::endl;

		return found;
	}



	/**
	 * @brief  partition the in mem valid range by the number of processes in communicator
	 * @note	overlap is added to the end of the result ranges.  This may make the result extend
	 * 			beyong initial in_mem_range.
	 * @note	assumes that both input parameters are within file range and in memory.
	 */
	::std::pair<range_type, range_type>
	partition(range_type const & in_mem_range_bytes,
				range_type const & valid_range_bytes) {

		// search for starting point using valid_range_bytes

		// restrict in_mem to within file
		range_type in_mem =
				range_type::intersect(in_mem_range_bytes, this->file_range_bytes);

		// and restrict valid to in memory data.
		range_type valid =
				range_type::intersect(valid_range_bytes, in_mem);


		range_type found = partition(valid);

		return std::make_pair(found, found);
	}

public:
	// this is needed to prevent overload name hiding.  see http://stackoverflow.com/questions/888235/overriding-a-bases-overloaded-function-in-c/888337#888337
	using BASE::read_range;

	/**
	 * @brief  bulk load the data and return it in a newly constructed vector.  block decomposes the range.  reuse vector
	 * @note   virtual so other overloads of this function can also block decompose the range.
	 * @param range_bytes	range to read, in bytes
	 * @param output		vector containing data as bytes.
	 */
	virtual range_type read_range(std::vector<unsigned char> & output,
				range_type const & range_bytes) {

		// first get rough partition
		range_type target = partition(range_bytes);

		// then read the range via sequential version
		reader.read_range(output, target);

		return target;
	}

	/**
	 * @brief constructor
	 * @param _filename 		name of file to open
	 * @param _comm				MPI communicator to use.
	 */
	partitioned_file(std::string const & _filename, size_t const & _overlap = 0UL, ::mxx::comm const & _comm = ::mxx::comm()) :
		::bliss::io::parallel::base_file(_filename, _comm), reader(_filename), overlap(_overlap) {};

	/// destructor
	virtual ~partitioned_file() {};  // will call super's unmap.


	// this is needed to prevent overload name hiding.  see http://stackoverflow.com/questions/888235/overriding-a-bases-overloaded-function-in-c/888337#888337
	using BASE::read_file;

	/**
	 * @brief  read the whole file.  reuse allocated file_data object.
	 * @note	virtual so that parallel readers can compute range to read.
	 * @param output 		file_data object containing data and various ranges.
	 */
	virtual void read_file(::bliss::io::file_data & output) {
		range_type in_mem_partitioned;
		range_type valid_partitioned;

		::std::tie(in_mem_partitioned, valid_partitioned) =
				partition(this->file_range_bytes, this->file_range_bytes);


		// then read the range via sequential version
		output.in_mem_range_bytes = reader.read_range(output.data, in_mem_partitioned);

		output.valid_range_bytes = valid_partitioned;
		output.parent_range_bytes = file_range_bytes;
	}
};

// TODO: 2 level parallel file io, first level at node, second level at processes in node.





class mpiio_base_file : public ::bliss::io::parallel::base_file {

	static_assert(sizeof(MPI_Offset) > 4, "ERROR: MPI_Offset is defined as an integer 4 bytes or less.  Do not use mpiio_file ");

protected:
	using BASE = ::bliss::io::parallel::base_file;

	/// overlap amount
	const size_t overlap;

	/// MPI file handle
	MPI_File fh;

	/// partitioner to use.
	::bliss::partition::BlockPartitioner<range_type> partitioner;

	std::string get_error_string(std::string const & op_name, int const & return_val) {
		char error_string[BUFSIZ];
		int length_of_error_string, error_class;
		std::stringstream ss;

		MPI_Error_class(return_val, &error_class);
		MPI_Error_string(error_class, error_string, &length_of_error_string);

		ss << "ERROR in mpiio: rank " << comm.rank() << " " << op_name << " " << this->filename << " error: " << error_string << std::endl;
		return ss.str();
	}

	std::string get_error_string(std::string const & op_name, int const & return_val, MPI_Status const & stat) {
		char error_string[BUFSIZ];
		int length_of_error_string, error_class;
		std::stringstream ss;

		MPI_Error_class(return_val, &error_class);
		MPI_Error_string(error_class, error_string, &length_of_error_string);

		ss << "ERROR in mpiio: rank " << comm.rank() << " " << op_name << " " << this->filename << " error: " << return_val << " [" << error_string << "]";

//		// status.MPI_ERROR does not appear to be decodable by error_class.  google search did not find how to decode it.
//		MPI_Error_class(stat.MPI_ERROR, &error_class);
//		MPI_Error_string(error_class, error_string, &length_of_error_string);

		ss << " MPI_Status error: [" << stat.MPI_ERROR << "]" << std::endl;

		return ss.str();
	}


	/// MPI_IO get file size
	size_t get_file_size() {
		if (fh == MPI_FILE_NULL) {
			std::stringstream ss;
			ss << "ERROR in mpiio: rank " << comm.rank() << " file " << this->filename << " not yet open " << std::endl;
			throw ::std::ios_base::failure(ss.str());
		}

		MPI_Offset s;

		int res = MPI_File_get_size(fh, &s);
		if (res != MPI_SUCCESS) {
			throw ::std::ios_base::failure(get_error_string("get_size", res));
		}

		return static_cast<size_t>(s);
	}

	/// opens a file. side effect computes size of the file.
	void open_file() {
		// first clear previously open file
		close_file();

		// open the file
		int res = MPI_File_open(this->comm, const_cast<char *>(this->filename.c_str()), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);

		if (res != MPI_SUCCESS) {
			throw ::std::ios_base::failure(get_error_string("open", res));
		}

		// ensure atomicity is turned off
		MPI_File_set_atomicity(fh, 0);
	}

	/// funciton for closing a file
	void close_file() {
		if (fh != MPI_FILE_NULL) {
			int res = MPI_File_close(&fh);
			if (res != MPI_SUCCESS) {
				throw ::std::ios_base::failure(get_error_string("close", res));
			}
			fh = MPI_FILE_NULL;
		}
	}



public:
	// this is needed to prevent overload name hiding.  see http://stackoverflow.com/questions/888235/overriding-a-bases-overloaded-function-in-c/888337#888337
	using BASE::read_range;

	/**
	 * @brief  bulk load the data and return it in a newly constructed vector.  block decomposes the range.  reuse vector
	 * @note   virtual so other overloads of this function can also block decompose the range.
	 * @param range_bytes	range to read, in bytes
	 * @param output		vector containing data as bytes.
	 */
	virtual range_type read_range(std::vector<unsigned char> & output,
				range_type const & range_bytes) {
		if (fh == MPI_FILE_NULL) {
			std::stringstream ss;
			ss << "ERROR in mpiio: rank " << comm.rank() << " file " << this->filename << " not yet open " << std::endl;
			throw ::std::ios_base::failure(ss.str());
		}

		// ensure valid range is used.
		range_type target =
				BASE::range_type::intersect(range_bytes, this->file_range_bytes);

		// do equal partition
		if (comm.size() > 1) {
			partitioner.configure(target, comm.size());
			target = partitioner.getNext(comm.rank());
		}

		// compute the size to read.
		range_type read_range = target;
		read_range.end += this->overlap;
		read_range.intersect(this->file_range_bytes);

		output.reserve(target.size() + overlap);  // reserve more

		output.resize(read_range.size());    // set size for reading.

		// then read the file using mpiio
		MPI_Status stat;

		// NOTE:  file offset is in units of byte (up to size_t).  number of elements to read has type int.

		size_t step_size = 1UL << 30 ;  // using 2^30 vs 2^31-2^15 does not make a huge performance difference (will only matter for low proc count anyways)
		size_t rem = read_range.size() % step_size;
//		size_t steps = read_range.size() / step_size;
		int count = 0;
		int res = MPI_SUCCESS;


// ======= DOES NOT WORK. number of elements has type int.
//		res =MPI_File_read_at_all(fh, read_range.start, output.data(), read_range.size(), MPI_BYTE, &stat);
//		if (res != MPI_SUCCESS) throw ::std::ios_base::failure(get_error_string("read", res, stat));


// ======= DOES NOT WORK to read rem first as BYTES, then read rest as 1GB blocks.  rem reads okay, but the first byte for big type fails.
//		res = MPI_File_read_at_all(fh, read_range.start, output.data(), rem, MPI_BYTE, &stat);
//		if (res != MPI_SUCCESS) throw ::std::ios_base::failure(get_error_string("read", res, stat));
//		// count is okay here.
//		res = MPI_Get_count(&stat, MPI_BYTE, &count);
//		if (res != MPI_SUCCESS) throw ::std::ios_base::failure(get_error_string("count", res));
//		if (static_cast<size_t>(count) != rem) {
//			std::stringstream ss;
//			ss << "ERROR in mpiio: rank " << comm.rank() << " remainder read error. request " << rem << " bytes got " << count << " bytes" << std::endl;
//			throw ::std::ios_base::failure(ss.str());
//		}
//
//		// now read the big data type (size of 1GB)
//		if (steps > 0) {
//			::mxx::datatype dt = ::mxx::get_datatype<unsigned char>().contiguous(step_size);  // make element  2^30 in size.
//		    res = MPI_File_read_at_all(fh, read_range.start + rem, output.data() + rem,
//	                         steps, dt.type(), &stat);
//			if (res != MPI_SUCCESS) throw ::std::ios_base::failure(get_error_string("read", res, stat));
//			// count comes back as -32677.
//			res = MPI_Get_count(&stat, dt.type(), &count);
//			if (res != MPI_SUCCESS) throw ::std::ios_base::failure(get_error_string("read", res));
//			if (static_cast<size_t>(count) != steps) {
//				std::stringstream ss;
//				ss << "ERROR in mpiio: rank " << comm.rank() << " remainder read error. request " << steps << " 2^30 byte blocks got " << count << " blocks" << std::endl;
//				throw ::std::ios_base::failure(ss.str());
//			}
//		}


// ============== does not work to read big data type first then rem.
//		if (steps > 0) {
//			// first blocks failed - reads a vector of 0's
//			::mxx::datatype dt = ::mxx::get_datatype<unsigned char>().contiguous(step_size);  // make element  2^30 in size.
//		    res = MPI_File_read_at_all(fh, read_range.start, output.data(),
//	                         steps, dt.type(), &stat);
//	   	    if (res != MPI_SUCCESS) throw ::std::ios_base::failure(get_error_string("read", res, stat));
//			// get count got -32766.
//		    res = MPI_Get_count(&stat, dt.type(), &count);
//        		if (res != MPI_SUCCESS) throw ::std::ios_base::failure(get_error_string("read", res));
//			if (static_cast<size_t>(count) != steps) {
//				std::stringstream ss;
//				ss << "ERROR in mpiio: rank " << comm.rank() << " remainder read error. request " << steps << " 2^30 byte blocks got " << count << " blocks" << std::endl;
//				throw ::std::ios_base::failure(ss.str());
//			}
//		}
//		res = MPI_File_read_at_all(fh, read_range.start + steps * step_size, output.data() + steps * step_size, rem, MPI_BYTE, &stat);
//		if (res != MPI_SUCCESS) throw ::std::ios_base::failure(get_error_string("read", res, stat));

//		res = MPI_Get_count(&stat, MPI_BYTE, &count);
//		if (res != MPI_SUCCESS) throw ::std::ios_base::failure(get_error_string("read", res));
//		if (static_cast<size_t>(count) != rem) {
//			std::stringstream ss;
//			ss << "ERROR in mpiio: rank " << comm.rank() << " remainder read error. request " << rem << " bytes got " << count << " bytes" << std::endl;
//			throw ::std::ios_base::failure(ss.str());
//		}


// ========= Single big data type that is the whole size does not work.  hangs when nprocs = 1 and for nprocs = 2, does not appear to populate data at all
//		::mxx::datatype dt = ::mxx::get_datatype<unsigned char>().contiguous(read_range.size());  // make element the whole range in size.
//		res = MPI_File_read_at_all(fh, read_range.start, output.data(),  1, dt.type(), &stat);
//		if (res != MPI_SUCCESS) throw ::std::ios_base::failure(get_error_string("read", res, stat));
//
// 		// getCount returns -32766, so cannot be used.
//		res = MPI_Get_count(&stat, dt.type(), &count);
//		if (res != MPI_SUCCESS) throw ::std::ios_base::failure(get_error_string("read", res));
//		if (static_cast<size_t>(count) != read_range.size()) {
//			std::stringstream ss;
//			ss << "ERROR in mpiio: rank " << comm.rank() << " remainder read error. request " << read_range.size() << " bytes got " << count << " bytes" << std::endl;
//			throw ::std::ios_base::failure(ss.str());
//		}

// ===========  iterative works
		size_t iter_step_size;
		for (size_t s = 0; s < read_range.size(); s += step_size, rem -= step_size) {
			iter_step_size = std::min(step_size, read_range.size() - s);
			res = MPI_File_read_at_all(fh, read_range.start + s, output.data() + s,
					iter_step_size, MPI_BYTE, &stat);
			if (res != MPI_SUCCESS) throw ::std::ios_base::failure(get_error_string("read", res, stat));

			res = MPI_Get_count(&stat, MPI_BYTE, &count);
			if (res != MPI_SUCCESS) throw ::std::ios_base::failure(get_error_string("read", res, stat));

			if (static_cast<size_t>(count) != iter_step_size) {
				std::stringstream ss;
				ss << "ERROR in mpiio: rank " << comm.rank() << " remainder read error. request " << iter_step_size << " bytes got " << count << " bytes" << std::endl;
				throw ::std::ios_base::failure(ss.str());
			}

		}

//		std::cout << "rank " << comm.rank() << " done reading " << read_range << std::endl;

		return target;
	}


	mpiio_base_file(::std::string const & _filename, size_t const _overlap = 0UL,  ::mxx::comm const & _comm = ::mxx::comm()) :
		::bliss::io::parallel::base_file(_filename, _comm), overlap(_overlap), fh(MPI_FILE_NULL) {
		open_file();
		this->file_range_bytes.end = get_file_size();
	};

	~mpiio_base_file() { close_file(); };

	// this is needed to prevent overload name hiding.  see http://stackoverflow.com/questions/888235/overriding-a-bases-overloaded-function-in-c/888337#888337
	using BASE::read_file;


};


/**
 * parallel file abstraction for MPIIO, block partition with overlap
 */
template <typename FileParser = ::bliss::io::BaseFileParser<unsigned char*> >
class mpiio_file : public ::bliss::io::parallel::mpiio_base_file {

	protected:
		using BASE = ::bliss::io::parallel::mpiio_base_file;


public:
		// this is needed to prevent overload name hiding.  see http://stackoverflow.com/questions/888235/overriding-a-bases-overloaded-function-in-c/888337#888337
		using BASE::read_range;


		mpiio_file(::std::string const & _filename, size_t const & _overlap = 0UL, ::mxx::comm const & _comm = ::mxx::comm()) :
			::bliss::io::parallel::mpiio_base_file(_filename, _overlap, _comm) {};

		~mpiio_file() {};

		// this is needed to prevent overload name hiding.  see http://stackoverflow.com/questions/888235/overriding-a-bases-overloaded-function-in-c/888337#888337
		using BASE::read_file;


		/**
		 * @brief  read the whole file.  reuse allocated file_data object.
		 * @note	virtual so that parallel readers can compute range to read.
		 * @param output 		file_data object containing data and various ranges.
		 */
		virtual void read_file(::bliss::io::file_data & output) {

			// then read the range via sequential version
			output.valid_range_bytes = read_range(output.data, this->file_range_bytes);

			// overlap is part of the read so there is no left shift needed.

			// set up the ranges.
			output.in_mem_range_bytes = output.valid_range_bytes;

			output.in_mem_range_bytes.end += overlap;
			output.in_mem_range_bytes.intersect(this->file_range_bytes);

			output.parent_range_bytes = file_range_bytes;
		}
};

/**
 * parallel file abstraction for MPIIO.  FASTQParser, which searches the input.
 */
template <>
class mpiio_file<::bliss::io::FASTQParser<unsigned char*> > : public ::bliss::io::parallel::mpiio_base_file {

	protected:
		using BASE = ::bliss::io::parallel::mpiio_base_file;


public:
		// this is needed to prevent overload name hiding.  see http://stackoverflow.com/questions/888235/overriding-a-bases-overloaded-function-in-c/888337#888337
		using BASE::read_range;


		mpiio_file(::std::string const & _filename, size_t const & _overlap = 0UL, ::mxx::comm const & _comm = ::mxx::comm()) :
			::bliss::io::parallel::mpiio_base_file(_filename, sysconf(_SC_PAGE_SIZE), _comm) {
			// specify 1 page worth as overlap
		};

		~mpiio_file() { };

		// this is needed to prevent overload name hiding.  see http://stackoverflow.com/questions/888235/overriding-a-bases-overloaded-function-in-c/888337#888337
		using BASE::read_file;



		/**
		 * @brief  read the whole file.  reuse allocated file_data object.
		 * @note	virtual so that parallel readers can compute range to read.
		 * @param output 		file_data object containing data and various ranges.
		 */
		virtual void read_file(::bliss::io::file_data & output) {

			// overlap is set to page size, so output will have sufficient space.
			// note that this is same behavior as the serial mmap_file

			// then read the range via sequential version
			range_type block_range = read_range(output.data, this->file_range_bytes);

//			std::cout << "rank " << this->comm.rank() << " read " << block_range << std::endl;
//			std::ostream_iterator<unsigned char> oit(std::cout, "");
//			std::copy(output.data.begin(), output.data.begin() + 100, oit);
//			std::cout << std::endl;

			range_type in_mem = block_range;
			in_mem.end = in_mem.start + output.data.size();

//			std::cout << "rank " << this->comm.rank() << " in mem " << in_mem << std::endl;


			// now search for the true start.
			::bliss::io::FASTQParser<unsigned char *> parser;
			parser.init_parser(output.data.data(), this->file_range_bytes,
					in_mem, block_range, this->comm);

			// mark the first entry found.
			size_t real_start = parser.find_first_record(output.data.data(),
					this->file_range_bytes, in_mem, in_mem);

//			std::cout << "rank " << this->comm.rank() << " real start " << real_start << std::endl;

//			std::cout << "rank " << this->comm.rank() << " data size before remove overlap " << output.data.size() << std::endl;

			// now clear the region outside of the block range  (the overlap portion)
			if (output.data.size() > block_range.size()) {
				output.data.erase(output.data.begin() + block_range.size(), output.data.end());
			}

//			std::cout << "rank " << this->comm.rank() << " data size after remove overlap " << output.data.size() << std::endl;

			// ==== shift values to the left (as vector?)
			// actually, not shift. need to do all to all - if a rank found no internal starts
			// first compute the target proc id via exscan
			bool not_found = (real_start >= block_range.end);  // if real start is outside of partition, not found
			real_start = std::min(real_start, block_range.end);
			int target_rank = not_found ? 0 : comm.rank();
			target_rank = ::mxx::exscan(target_rank, [](int const & x, int const & y) {
				return (x < y) ? y : x;
			}, comm);

//			std::cout << "rank " << this->comm.rank() << " adjusted real start " << real_start << std::endl;
//			std::cout << "rank " << this->comm.rank() << " target rank " << target_rank << std::endl;


			// MPI 2 does not have neighbor_ collective operations, so we can't create
			// graph with sparse edges.  in any case, the amount of data we send should be
			// small so alltoallv should be enough
			std::vector<size_t> send_counts(comm.size(), 0);
			if (comm.rank() > 0) send_counts[target_rank] = real_start - block_range.start;

			// copy the region to shift
			std::vector<unsigned char> shifted =
					::mxx::all2allv(output.data, send_counts, comm);
			output.data.insert(output.data.end(), shifted.begin(), shifted.end());

//			std::cout << "rank " << this->comm.rank() << " shifted " << shifted.size() << std::endl;
//			std::cout << "rank " << this->comm.rank() << " new data size " << output.data.size() << std::endl;

			// adjust the ranges.
			output.in_mem_range_bytes = block_range;
			output.in_mem_range_bytes.end = block_range.start + output.data.size();

//			std::cout << "rank " << this->comm.rank() << " final in mem " << output.in_mem_range_bytes << std::endl;


			output.valid_range_bytes.start = real_start;
			output.valid_range_bytes.end =
					not_found ? block_range.end : output.in_mem_range_bytes.end;

//			std::cout << "rank " << this->comm.rank() << " final valid " << output.valid_range_bytes << std::endl;

			output.parent_range_bytes = file_range_bytes;

//			std::cout << "rank " << this->comm.rank() << " file  " << output.parent_range_bytes << std::endl;

		}
};


}  // namespace parallel


//======================================


using RangeType = bliss::partition::range<size_t>;
using DataType = unsigned char;

struct MemData {
	using iterator = DataType*;
	DataType * data;
	RangeType mem_range;
	RangeType valid_range;
	RangeType parent_range;

	DataType * begin() {
		return data + valid_range.start - mem_range.start;
	}
	DataType * end() {
		return data + valid_range.end - mem_range.start;
	}
	RangeType getRange() {
		return valid_range;
	}
};




/**
 * @brief unmaps a file region from memory
 * @note      AGNOSTIC of overlaps
 */
void unload_data(MemData & data) {

	if (data.data != nullptr) {
		delete [] data.data;
		data.data = nullptr;
		data.mem_range.start = data.mem_range.end;
		data.valid_range.start = data.valid_range.end;
	}
}



namespace memmap {





// load a file.  single process reading entire file.
template <typename FileParser = ::bliss::io::BaseFileParser<DataType *> >
MemData load_file(::std::string const & filename, size_t overlap = 0) {

	bliss::io::mmap_file fobj(filename);
	bliss::io::file_data fdata = fobj.read_file();

	MemData mdata;

	mdata.mem_range = fdata.in_mem_range_bytes;
	mdata.valid_range = fdata.valid_range_bytes;
	mdata.parent_range = fdata.parent_range_bytes;
	mdata.data = new unsigned char[fdata.in_mem_range_bytes.size()];

	memmove(mdata.data, fdata.data.data(), fdata.in_mem_range_bytes.size());

	// return the data
	return mdata;
}





} // namespace memmap

namespace stdio {



// load a file.  single process reading entire file.
template <typename FileParser = ::bliss::io::BaseFileParser<DataType *> >
MemData load_file(::std::string const & filename, size_t overlap = 0) {

	bliss::io::stdio_file fobj(filename);
	bliss::io::file_data fdata = fobj.read_file();

	MemData mdata;

	mdata.mem_range = fdata.in_mem_range_bytes;
	mdata.valid_range = fdata.valid_range_bytes;
	mdata.parent_range = fdata.parent_range_bytes;
	mdata.data = new unsigned char[fdata.in_mem_range_bytes.size()];

	memmove(mdata.data, fdata.data.data(), fdata.in_mem_range_bytes.size());

	// return the data
	return mdata;
}


} // namespace stdio

namespace parallel {

#if defined(USE_MPI)


namespace memmap {

// for FASTQ sequences (short).  multiple processes reading different parts of a file.  use mpi communicator to
// get start and end.
template <typename FileParser = ::bliss::io::BaseFileParser<DataType *> >
typename ::std::enable_if<(::std::is_same<FileParser, ::bliss::io::FASTQParser<DataType *> >::value), MemData>::type
load_file(::std::string const & filename, ::mxx::comm const & comm, size_t overlap = 0) {


	::bliss::io::parallel::partitioned_file<::bliss::io::mmap_file, FileParser> fobj(filename, 0, comm);
	::bliss::io::file_data fdata = fobj.read_file();

	MemData mdata;

	mdata.mem_range = fdata.in_mem_range_bytes;
	mdata.valid_range = fdata.valid_range_bytes;
	mdata.parent_range = fdata.parent_range_bytes;
	mdata.data = new unsigned char[fdata.in_mem_range_bytes.size()];

	memmove(mdata.data, fdata.data.data(), fdata.in_mem_range_bytes.size());


	// return the data
	return mdata;
}


// for FASTA sequences (short).  multiple processes reading different parts of a file
template <typename FileParser = ::bliss::io::BaseFileParser<DataType *> >
typename ::std::enable_if<!(::std::is_same<FileParser, ::bliss::io::FASTQParser<DataType *> >::value), MemData>::type
load_file(::std::string const & filename, ::mxx::comm const & comm, size_t overlap = 0) {


	bliss::io::parallel::partitioned_file<bliss::io::mmap_file, FileParser> fobj(filename, 0, comm);
	bliss::io::file_data fdata = fobj.read_file();

	MemData mdata;

	mdata.mem_range = fdata.in_mem_range_bytes;
	mdata.valid_range = fdata.valid_range_bytes;
	mdata.parent_range = fdata.parent_range_bytes;
	mdata.data = new unsigned char[fdata.in_mem_range_bytes.size()];

	memmove(mdata.data, fdata.data.data(), fdata.in_mem_range_bytes.size());

	// return the data
	return mdata;
}

// TODO: load multiple files concurrently, via subcommunicator.
} // namespace memmap


namespace mpi_io {
// for FASTA sequences (short).  multiple processes reading different parts of a file
template <typename FileParser = ::bliss::io::BaseFileParser<DataType *> >
typename ::std::enable_if<!(::std::is_same<FileParser, ::bliss::io::FASTQParser<DataType *> >::value), MemData>::type
load_file(::std::string const & filename, ::mxx::comm const & comm, size_t overlap = 0) {


	bliss::io::parallel::mpiio_file<FileParser> fobj(filename, 0, comm);
	bliss::io::file_data fdata = fobj.read_file();

	MemData mdata;

	mdata.mem_range = fdata.in_mem_range_bytes;
	mdata.valid_range = fdata.valid_range_bytes;
	mdata.parent_range = fdata.parent_range_bytes;
	mdata.data = new unsigned char[fdata.in_mem_range_bytes.size()];

	memmove(mdata.data, fdata.data.data(), fdata.in_mem_range_bytes.size());

	// return the data
	return mdata;
}
}  // namespace mpi_io

#endif
} // parallel


} // io

} // bliss



#endif /* FILE_HPP_ */
