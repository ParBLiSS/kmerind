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
#include <fcntl.h>      // for open
#include <sstream>      // stringstream


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
 * @brief  opens a file and return the file handle.  Throws IOException if can't open.
 * @param fn  name of file to open
 * @return    file handle/descriptor for the file
 */
int open_file(::std::string const & fn) {
  // open the file and get a handle.
  int fd = open64(fn.c_str(), O_RDONLY);
  if (fd == -1)
  {
    // if open failed, throw exception.
    ::std::stringstream ss;
    int myerr = errno;
    ss << "ERROR in file open: ["  << fn << "] error " << myerr << ": " << strerror(myerr);
    throw new ::std::ios_base::failure(ss.str());
  }
  return fd;
}

/**
 * @brief closes a file
 * @param fd  the file descriptor of the file to close.
 */
void close_file(int &fd) {
  if (fd >= 0) {
    close(fd);
    fd = -1;
  }
}


/**
 * @brief get the file size  (supports 64bit) of a file given the file descriptor
 * @note  the method uses fstat64, which does not require opening the file and seeking.  also avoids file encoding and text/binary issues.
 * @param fd    file descriptor, from fstat or fstat64
 * @return      size in bytes, data type size_t.
 */
size_t get_file_size(const int& fd) throw (bliss::io::IOException) {

  struct stat64 filestat;
  // get the file state
  int ret = fstat64(fd, &filestat);

  // handle any error
  if (ret < 0) {
    throw ::std::ios_base::failure("ERROR in file size detection");
  }

  return static_cast<size_t>(filestat.st_size);
}

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

struct MapData {
	DataType * memdata;
	RangeType memrange;
};






/**
 * @brief     map the specified portion of the file to memory.
 * @note      AGNOSTIC of overlaps
 * @param r   range specifying the portion of the file to map.
 * @return    memory address (pointer) to where the data is mapped.
 *            TODO: should return std::pair<PointerType, page aligned offset>
 */
MapData map_file_range(int const & fd, const RangeType & file_range, const RangeType &r) {
	// if no file
	if (fd < 0) {
		throw ::std::ios_base::failure("ERROR specified invalid file handle for mmap.  was 'open_file' called?");
	}

	MapData data;
	data.memdata = nullptr;
	data.memrange = r;

	// if no range
	if (r.size() == 0) return data;

	// get the page size for the system
	const size_t page_size = sysconf(_SC_PAGE_SIZE);

	/// memory map.  get start and end positions that are page aligned.
    data.memrange.start = RangeType::align_to_page(r, page_size);
    data.memrange.end = r.end % page_size;
    data.memrange.end = (data.memrange.end == 0) ? r.end : (r.end + (page_size - data.memrange.end));

    data.memrange.end = ::std::min(data.memrange.end, file_range.end);

    // NOT using MAP_POPULATE.  it slows things done when testing on single node.  NOTE HUGETLB not supported for file mapping.
    data.memdata = (DataType*)mmap64(nullptr, data.memrange.size(),
                                 PROT_READ,
                                 MAP_PRIVATE, fd,
                                 data.memrange.start);  // TODO: range means bytes or number of elements?

    // if mmap failed,
    if (data.memdata == MAP_FAILED)
    {
		// print error through exception.
		std::stringstream ss;
		int myerr = errno;
		ss << "ERROR in mmap: " << myerr << ": " << strerror(myerr);
		throw ::std::ios_base::failure(ss.str());
    }

	// set the madvice info.
	int madv_result = madvise(data.memdata, data.memrange.size(), MADV_SEQUENTIAL | MADV_WILLNEED);
	if ( madv_result == -1) {
		std::stringstream ss;
		int myerr = errno;
		ss << "ERROR in madvise: " << myerr << ": " << strerror(myerr);
		throw std::ios_base::failure(ss.str());
	}

	return data;
}

void unmap_file_range(MapData & data) {
	// unmap it.
	if ((data.memdata != nullptr) && (data.memrange.size() > 0))
		munmap(data.memdata, data.memrange.size());
}

/**
 * @param data   	    range stored in memory currently.
 * @param mem_range	    range to keep in memory.
 * @param valid_range   range containing meaningful data (for current process)
 */
MemData load_data(MapData & data, RangeType const & parent_range, RangeType const & mem_range, RangeType const & valid_range ) {

	if (! data.memrange.contains(mem_range)) {
		std::cout << "data memrange = " << data.memrange << " mem_range " << mem_range << std::endl;
		throw std::invalid_argument("ERROR: load data specified to-keep memory range outside of in memory data.");
	}
	if (! mem_range.contains(valid_range))
		throw std::invalid_argument("ERROR: load data specified valid_range outside of to-keep memory range");


	MemData output;
	output.mem_range = mem_range;
	output.valid_range = valid_range;
	output.data = nullptr;

	if ( (data.memdata == nullptr) || (data.memrange.size() == 0) )
		return output;

	if (mem_range.size() == 0)
		return output;

	// copy the data into memory
	output.data = new DataType[mem_range.size()];
	memcpy(output.data, data.memdata + (mem_range.start - data.memrange.start), mem_range.size());

	output.parent_range = parent_range;

	return output;
}







// load a file.  single process reading entire file.
template <typename FileParser = ::bliss::io::BaseFileParser<DataType *> >
MemData load_file(::std::string const & filename, size_t overlap = 0) {

	// open file
	int fd = ::bliss::io::open_file(filename);

	// get file size
	size_t fsize = ::bliss::io::get_file_size(fd);

	// get the file range
	RangeType file_range(0, fsize);

	// map the data
	MapData mapped = ::bliss::io::memmap::map_file_range(fd, file_range, file_range);

	// copy the data
	MemData mdata = ::bliss::io::memmap::load_data(mapped, file_range, file_range, file_range);

	// unmap the data
	unmap_file_range(mapped);

	// close file
	::bliss::io::close_file(fd);

	// return the data
	return mdata;
}





} // namespace memmap


namespace parallel {

#if defined(USE_MPI)

// parallel file size - 1 proc load and broadcast.
size_t get_file_size(const int & fd, mxx::comm const & comm) {
  size_t file_size = 0;
  if (comm.rank() == 0) {
	// get size on rank 0, then broadcast.
	file_size = ::bliss::io::get_file_size(fd);
  }

  // broadcast.
  int p = comm.size();

  if (p > 1) MPI_Bcast(&file_size, 1, MPI_UNSIGNED_LONG, 0, comm);

  return file_size;
}


namespace memmap {

// for FASTQ sequences (short).  multiple processes reading different parts of a file.  use mpi communicator to
// get start and end.
template <typename FileParser = ::bliss::io::BaseFileParser<DataType *> >
typename ::std::enable_if<(::std::is_same<FileParser, ::bliss::io::FASTQParser<DataType *> >::value), MemData>::type
load_file(::std::string const & filename, mxx::comm const & comm, size_t overlap = 0) {


	// get the MPI size and rank
	int nprocs = comm.size();
	MPI_Comm_size(comm, &nprocs);

	// single proc version
	if (nprocs == 1) return ::bliss::io::memmap::load_file(filename, overlap);

	int rank = comm.rank();

	// open file
	int fd = ::bliss::io::open_file(filename);

	// get the file range
	size_t fsize = ::bliss::io::parallel::get_file_size(fd, comm);
	RangeType file_range(0, fsize);

	// get the partitioned range
	::bliss::partition::BlockPartitioner<RangeType> partitioner;
	partitioner.configure(file_range, nprocs);

	// first range
	RangeType r = partitioner.getNext(rank);

	// set search range to reach end of file.  rely on system to not load too much.
	RangeType search = r;
	search.end = file_range.end;

	//==  adjust the partitioned range

	// load the range.  each rank loads to end of file.  (for small file, this is okay.  for large file???)
	::bliss::io::memmap::MapData mapped = ::bliss::io::memmap::map_file_range(fd, file_range, search);

	// search for start
	FileParser parser;
	parser.init_parser(mapped.memdata, file_range, mapped.memrange, search, comm);

	RangeType final;
	final.start = parser.find_first_record(mapped.memdata, file_range, mapped.memrange, search);

//	std::cout << "rank " << comm.rank() << " start " << final.start << std::endl;

	//=== modify the range via mpi
	// now get the end from next proc that found a start.
	// first setup the end element - if search returned end, then set final.start to end of file.
	if (final.start == r.end) final.start = file_range.end;  // nothing found.

//	std::cout << "rank " << comm.rank() << " start " << final.start << std::endl;


	// not using comm split because the subcommunicator may be constructed in worse than log p time.
	// next we do reverse exclusive scan with rev communicator
	final.end = mxx::exscan(final.start, [](size_t const & x, size_t const & y) {
		return ::std::min(x, y);
	}, comm.reverse());
	if (comm.rank() == (comm.size() - 1)) final.end = file_range.end;   // last process.

	// finally, check final_start again, and if start is pointing to end of file, let it point to the current range's end.
	if (final.start == file_range.end) final.start = final.end;

	//=== reload the data
	::bliss::io::memmap::unmap_file_range(mapped);
	mapped = ::bliss::io::memmap::map_file_range(fd, file_range, final);

	// copy the data using the actual start/end
	MemData mdata = ::bliss::io::memmap::load_data(mapped, file_range, final, final);

	// unmap the data.
	::bliss::io::memmap::unmap_file_range(mapped);

	// close file
	::bliss::io::close_file(fd);

	// return the data
	return mdata;
}


// for FASTA sequences (short).  multiple processes reading different parts of a file
template <typename FileParser = ::bliss::io::BaseFileParser<DataType *> >
typename ::std::enable_if<!(::std::is_same<FileParser, ::bliss::io::FASTQParser<DataType *> >::value), MemData>::type
load_file(::std::string const & filename, mxx::comm const & comm, size_t overlap = 0) {


	// get the MPI size and rank
	int nprocs = comm.size();

	// single proc version
	if (nprocs == 1) return ::bliss::io::memmap::load_file(filename, overlap);

	int rank = comm.rank();

	// open file
	int fd = ::bliss::io::open_file(filename);

	// get the file range
	size_t fsize = ::bliss::io::parallel::get_file_size(fd, comm);
	RangeType file_range(0, fsize);

	// get the partitioned range
	::bliss::partition::BlockPartitioner<RangeType> partitioner;
	partitioner.configure(file_range, nprocs);

	// first range
	RangeType r = partitioner.getNext(rank);

	//==  adjust the partitioned range
	RangeType final = r;
	if (comm.rank() < (comm.size() - 1)) final.end += overlap;

	std::cout << "rank: " << comm.rank() << " final range : " << final << ::std::endl;

	// then load it.
	::bliss::io::memmap::MapData mapped = ::bliss::io::memmap::map_file_range(fd, file_range, final);

	// copy the data using the actual start/end
	MemData mdata = ::bliss::io::memmap::load_data(mapped, file_range, final, r);

	// unmap the data.
	::bliss::io::memmap::unmap_file_range(mapped);

	// close file
	::bliss::io::close_file(fd);

	// return the data
	return mdata;
}


// TODO: load multiple files concurrently, via subcommunicator.


} // namespace memmap

#endif
} // parallel


} // io

} // bliss



#endif /* FILE_HPP_ */
