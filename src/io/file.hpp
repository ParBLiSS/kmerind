/*
 * Copyright 2015 Georgia Institute of Technology
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/*
 * file.hpp
 *
 * @note  special notes for Lustre:
 *  from http://www.nas.nasa.gov/hecc/support/kb/lustre-best-practices_226.html
 *  from http://www.opensfs.org/wp-content/uploads/2014/04/D1_S8_Lustre2.5PerformanceEvaluation.pdf.  IOPS limited.
 *  1. avoid repeated stat - e.g. used for getting file size.  lustre IOS is limited.  possible symptom: some processes may not find the file
 *  2. avoid multiple processes opening same file at the same time. symptom:  some processes may not find the file. aggregating is better
 *  3. open once, close once, instead of repeated open
 *  4. stripe align IO requests (default page size is 4K for mmap, so is already aligned, but need to prefetch)
 *  5. unlocked read.
 *
 * based on this https://www.nics.tennessee.edu/computing-resources/file-systems/io-lustre-tips
 * it appears that strip size of 32MB is best.
 *
 *  Created on: Feb 28, 2016
 *      Author: tpan
 */

#ifndef FILE_HPP_
#define FILE_HPP_

#include "bliss-config.hpp"

#include <string>
#include <cstring>      // memcpy, strerror

#include <ios>          // ios_base::failure
#include <iostream>     // ios_base::failure
#include <unistd.h>     // sysconf, usleep, lseek,
#include <sys/mman.h>   // mmap
#include <fcntl.h>      // for open64 and close
#include <sstream>      // stringstream
#include <exception>    // std exception

#if defined(USE_MPI)
#include <mpi.h>
#endif

// mxx
#include <mxx/collective.hpp>
#include <mxx/shift.hpp>

#include <io/io_exception.hpp>

#include <io/file_loader.hpp>
#include <io/fastq_loader.hpp>
#include <io/fasta_loader.hpp>
#include <io/unix_domain_socket.h>
#include <partition/range.hpp>

#include <utils/exception_handling.hpp>

#include <utils/memory_usage.hpp>

// define so that pread64 can support reading larger than int bytes.
// http://www.ibm.com/support/knowledgecenter/ssw_i5_54/apis/pread64.htm
#define _LARGE_FILE_API

// DONE: directly expose the mmapped region  (possible a variant for mmap without caching.)
// DONE: large number of PARALLEL FOPEN AND FREAD is not good.  mmap is okay because of common file descriptor between processes?  can that be guaranteed?  or is it okay because of open?
//        lseek and read may be better, provide the same file descriptor can be used for processes on the same node.
// TODO: read then shuffle boundaries for partitioned file, FASTQ case.
// TODO: change file open behavior to reduce congestion. - open when using, retry until success.
// DONE: util to clear linux disk cache  http://www.linuxatemyram.com/play.html,  or O_DIRECT is supposed to bypass cache,
//        but seems to carry a lot of complications.

// mmap_file
//    mmap range:  whole file, access region, or sliding window?
//    first, mmap will adapt to physical memory available.  mmapping more than physical + swap is okay
//    second.  maxrss reporting includes mmapped region, up to physical memory, but actual application memory use is lower.
//    a. whole large file by each process - The TLB table will contain Page Table Entries for whole mapped range.
//    b. regions of interest - TLB may still be somewhat large, and is dependent on size of input and num procs.
//        http://stackoverflow.com/questions/8506366/does-mmap-with-map-noreserve-reserve-physical-memory
//        https://forums.freebsd.org/threads/48980/
//    c. sliding mmapped window may reduce number of Page Table Entries and thus invalidation for other processes.
//      overhead of map/unmap is potentially expensive.  no test to indicate this is beneficial, so don't do it.
//    CONCLUSION:  map access region.

//    concurrent access to same file:
//      open64 is used with mmap.  fds are assigned sequentially and separately by each process. if they are the same, it's coincidental.
//        confirmed via mpirun -np 4 on 1 node
//      multiple calls create multiple open file descriptions for same inode, and each can seek/etc.
//      sharing file descriptor via unix socket results in same thing as dup(2),
//          since offset is shared, need to use pread, which does not modify the file descriptor
//          http://stackoverflow.com/questions/1997622/can-i-open-a-socket-and-pass-it-to-another-process-in-linux
//          http://blog.varunajayasiri.com/passing-file-descriptors-between-processes-using-sendmsg-and-recvmsg
//          http://www.thomasstover.com/uds.html
//      testing with 4096 procs still appears to work correctly for mmap.  2 possibilities :
//        1. kernel is aggregating the "open files" - perhaps at inode level?
//        2. lustre allows more concurrent file descriptors (from open64) than FILE streams (from fopen)?  untested and unlikely.
//      mapping same region of same file will return different virtual address in process space.  if same, coincidental?
//      however, the physical memory for the mapping is shared if MAP_SHARED is used.
//        https://www.mkssoftware.com/docs/man3/mmap.3.asp
//        each proc is going to read different pages. so this does not matter.
//    CONCLUSION:  the number of concurrent open file descriptors is an issue.  See TODO below

//    MAP_SHARED can exceed physical + swap memory.  MAP_PRIVATE requires reservation.
//        on linux, MAP_PRIVATE + MAP_NORESERVE allows swapspace overcommitting.
//        http://stackoverflow.com/questions/7222164/mmap-an-entire-large-file
//    MAP_NORESERVE prevents allocating swap (swap is more relevant to write).
//        http://stackoverflow.com/questions/8506366/does-mmap-with-map-noreserve-reserve-physical-memory
//    CONCLUSION:  use MAP_SHARED | MAP_NORESERVE

//    largest file size for mmap?  length up to size_t, but is probably limited by addressable memory addresses + ulimit.
//    CONCLUSION:  change nothing here.

//  DONE: minimize the number of open files.
//    1 proc per node open file.  then
//        1. share file descriptor, or
//        2. mmap cache/fread on 1 and share (require contiguous, and idle C-1 cores, communication at worst shared mem at best)
//        3. don't open file multiple times.
//        4. serialize/space out file open over time.
//            delay,
//            open only when doing something
//        5. map remains even when file is closed.
//    chose 1.

//  DONE: other things to change:
//      1. DONE. for mmap, reuse file handle/mmap, and mremap as needed (not needed).
//      2. DONE. use open/lseek/pread instead of fopen/fseek/fread/
//

// DONE:  open/lseek/read instead of fopen/fseek/fread
// DONE:  refactored mmap_file with a mapped_data object
// DONE:  remove 1 extra mmap from FASTQParser partitioned_file
// TODO:  move file open/close to closer to actual reading
//          close right after map, before unmap
// DONE:  copy file descriptor to processes - only 1 on each node opens.
//          all closes (same behavior as when using dup)
//          require sequential constructors to reuse fd.
//          NOTE: fread will move the fd, so sharing a file descriptor does not work.
//          NOTE: pread64 can only read 2GB chunks at a time.  offset is 64bit though.



namespace bliss {

namespace io {

/**
 * @brief loaded file data.
 */
struct file_data {
	using container = std::vector<unsigned char>;
  using iterator = typename container::iterator;
  using const_iterator = typename container::const_iterator;

	// type of ranges
	using range_type = ::bliss::partition::range<size_t>;

	// range from which the data came
	range_type parent_range_bytes;

	// range loaded in memory.  INCLUDES OVERLAP
	range_type in_mem_range_bytes;

	// valid range for this.  EXCLUDES OVERLAP
	range_type valid_range_bytes;

	// storage for actual data
	container data;

	/// beginning of the valid range
	iterator begin() {
		return data.begin() + valid_range_bytes.start - in_mem_range_bytes.start;
	}
	/// end of valid range
	iterator end() {
		return data.begin() + valid_range_bytes.end - in_mem_range_bytes.start;
	}

	/// beginning of the valid range
	const_iterator cbegin() const {
		return data.cbegin() + valid_range_bytes.start - in_mem_range_bytes.start;
	}
	/// end of valid range
	const_iterator cend() const {
		return data.cbegin() + valid_range_bytes.end - in_mem_range_bytes.start;
	}


	/// start of inmem range
	iterator in_mem_begin() {
		return data.begin();
	}
	/// end of in mem range
	iterator in_mem_end() {
		return data.end();
	}


	/// start of inmem range
	const_iterator in_mem_cbegin() const {
		return data.cbegin();
	}
	/// end of in mem range
	const_iterator in_mem_cend() const {
		return data.cend();
	}

	range_type getRange() const {
		return valid_range_bytes;
	}
};



/**
 * mmapped data.  wrapper for moving it around.
 */
class mapped_data {
  protected:
    /// pointer to mapped address
    unsigned char * data;

    /// range type
    using range_type = ::bliss::partition::range<size_t>;

    /// mapped range
    range_type range_bytes;

    const size_t page_size;

  public:


    /**
     * @brief   map the specified portion of the file to memory.
     * @note    AGNOSTIC of overlaps
     * @param range_bytes    range specifying the portion of the file to map.
     */
    mapped_data(int const & _fd, range_type const & target) :
      data(nullptr), range_bytes(0, 0),
      page_size(sysconf(_SC_PAGE_SIZE))
    {

      // if no file
      if (_fd == -1) {
        throw ::bliss::utils::make_exception<::bliss::io::IOException>("ERROR: map: file is not yet open.");
      }

      if (target.size() == 0) {
        // print error through exception.
        std::cout << "WARN: map: zero bytes in range" << std::endl;

        return;
      }

      // get start and end positions that are page aligned.
      range_bytes.start = range_type::align_to_page(target, page_size);
      range_bytes.end = target.end;  // okay for end not to align - made 0.

      // NOT using MAP_POPULATE. (SLOW)  no need to prefault the entire range - use read ahead from madvice.
      // NOTE HUGETLB not supported for file mapping, only anonymous.  also, kernel has to enable it and system has to have it reserved,
      // MAP_SHARED so that we don't have CoW (no private instance) (potential sharing between processes?)  slightly slower by 1%?
      // MAP_NORESERVE so that swap is not allocated.
      data = (unsigned char*)mmap64(nullptr, range_bytes.size(),
                                   PROT_READ,
                                   MAP_SHARED | MAP_NORESERVE, _fd,
                                   range_bytes.start);

      // if mmap failed,
      if (data == MAP_FAILED)
      {
        // print error through exception.
        std::stringstream ss;
        int myerr = errno;
        ss << "ERROR in mmap: " << myerr << ": " << strerror(myerr);
        throw ::bliss::utils::make_exception<::bliss::io::IOException>(ss.str());
      }

      // set the madvice info.  SEQUENTIAL vs RANDOM does not appear to make a difference in running time.
      int madv_result = madvise(data, range_bytes.size(),
          MADV_SEQUENTIAL | MADV_WILLNEED);
      if ( madv_result == -1 ) {
        std::stringstream ss;
        int myerr = errno;
        ss << "ERROR in madvise: " << myerr << ": " << strerror(myerr);

        throw ::bliss::utils::make_exception<bliss::io::IOException>(ss.str());
      }
      // for testing
      //std::cout << "serial fd=" << _fd << " mapped region = " << range_bytes << " pointer is " << (const void*)data << ::std::endl;
    }

    ~mapped_data() {
      // unmap it.
      if ((data != nullptr) && (range_bytes.size() > 0)) {
        munmap(data, range_bytes.size());

        data = nullptr;
        range_bytes.start = range_bytes.end;
      }
    }

    // copy constructor and assignment operators are deleted.
    mapped_data(mapped_data const & other) = delete;
    mapped_data& operator=(mapped_data const & other) = delete;

    mapped_data(mapped_data && other) :
      data(other.data), range_bytes(other.range_bytes), page_size(other.page_size) {

      other.data = nullptr;
      other.range_bytes.start = other.range_bytes.end;
    }
    mapped_data& operator=(mapped_data && other) {
      // first unmap
      if ((data != nullptr) && (range_bytes.size() > 0))
        munmap(data, range_bytes.size());

      // then move other to here.
      data = other.data;      other.data = nullptr;
      range_bytes = other.range_bytes;   other.range_bytes.start = other.range_bytes.end;

      return *this;
    }



    /// accessor for mapped data
    unsigned char * const & get_data() const {
      return data;
    }

    /// explicit cast operator
    explicit operator unsigned char *() const {
      return data;
    }

    /// accessor for mapped data range
    range_type const & get_range() const {
      return range_bytes;
    }

    /// get the size of this mapping.
    size_t size() const {
      return range_bytes.size();
    }

    /// accessor for page size.
    size_t const & get_page_size() const {
      return page_size;
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

	/// file descriptor
	int fd;

	/// size of file in bytes
	range_type file_range_bytes;

	/// virtual function for computing the size of a file.
	size_t get_file_size() {

    if (this->filename.length() == 0) return 0UL;

//	  printf("base file get_file_size\n");

		// get the file size.
		struct stat64 filestat;
		// get the file state
		int ret = stat64(filename.c_str(), &filestat);

		// handle any error
		if (ret < 0) {
			::std::stringstream ss;
			int myerr = errno;
			ss << "ERROR : bliss::io::base_file::get_file_size: ["  << filename << "] " << myerr << ": " << strerror(myerr);
			throw ::bliss::utils::make_exception<::bliss::io::IOException>(ss.str());
		}

		return static_cast<size_t>(filestat.st_size);
	}

  /// opens a file. side effect computes size of the file.
  void open_file() {
//printf("open seq file\n");

    if (this->filename.length() == 0) return;

    // first close file.
    this->close_file();

    // open the file and get a handle.
    this->fd = open64(this->filename.c_str(), O_RDONLY);
    if (this->fd == -1)
    {
      this->file_range_bytes.end = 0;

      // if open failed, throw exception.
      ::std::stringstream ss;
      int myerr = errno;
      ss << "ERROR in base_file open: ["  << this->filename << "] error " << myerr << ": " << strerror(myerr);
      throw ::bliss::utils::make_exception<::bliss::io::IOException>(ss.str());
    }
  }

  /// funciton for closing a file
  void close_file() {
    if (this->fd >= 0) {
      close(this->fd);
      this->fd = -1;
      this->file_range_bytes.end = 0;
    }
  }




  /**
   * @brief initializes a file for reading/writing.  for use when the size is already obtained, as in the compositional parallel file case.
   * @note  do not call get_file_size here because
   *      parallel_file inherits this class and defines a parallel get_file_size.
   *      we want to call this constructor from subclass.
   *
   * @param _filename   name of file to open
   * @param _file_size  size of file, previously obtained.
   */
  base_file(std::string const & _filename, size_t const & _file_size, size_t const & delay_ms = 0) :
    filename(_filename), fd(-1), file_range_bytes(0, _file_size) {

    if (delay_ms > 0) usleep(delay_ms * 1000UL);

    this->open_file();
  };

  /**
   * @brief constructor that takes an existing open file descriptor and size.
   * @param _fd   open file descriptor
   * @param _file_size  size of file, previously obtained.
   */
  base_file(int const & _fd, size_t const & _file_size) :
    filename(::std::string()), fd(dup(_fd)), file_range_bytes(0, _file_size) {};

public:

	/**
	 * @brief  bulk load all the data and return it in a newly constructed vector
	 * @param range_bytes	  range to read, in bytes
	 * @return				  vector containing data as bytes.
	 */
	virtual typename ::bliss::io::file_data::container read_range(range_type const & range_bytes) {
		typename ::bliss::io::file_data::container out;
		this->read_range(out, range_bytes);
		return out;
	}

	/**
	 * @brief  bulk load all the data and return it in a newly constructed vector.  reuse vector
	 * @note	virtual so different file reading mechanisms can be defined
	 * @param range_bytes	range to read, in bytes
	 * @param output		vector containing data as bytes.
	 */
	virtual range_type read_range(typename ::bliss::io::file_data::container & output, range_type const & range_bytes) = 0;


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
		filename(_filename), fd(-1), file_range_bytes(0, 0) {
		this->file_range_bytes.end = this->get_file_size();
	  this->open_file();
	};

	/// default destructor
	virtual ~base_file() {
    this->close_file();
	};


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
	size_t size() {
		return file_range_bytes.end;
	};

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

public:

	// this is needed to prevent overload name hiding.  see http://stackoverflow.com/questions/888235/overriding-a-bases-overloaded-function-in-c/888337#888337
	using BASE::read_range;

	/**
	 *  bulk load all the data and return it in pre-allocated vector
	 *  @param output   vector where the results are to be stored.
	 */
	virtual typename BASE::range_type read_range(typename ::bliss::io::file_data::container & output,
			typename BASE::range_type const & range_bytes) {

	  typename BASE::range_type file_range = BASE::range_type::intersect(range_bytes, this->file_range_bytes);

	  output.clear();

	  if (file_range.size() == 0) {
	    return file_range;
	  }

		// map
		mapped_data md(this->fd, file_range);

		//
		unsigned char * md_data = md.get_data();
		typename BASE::range_type mapped_range = md.get_range();

    if (md.size() == 0) {
      std::cout << "WARNING: mapped data size is 0" << std::endl;
      return mapped_range;
    }
		if (md_data == nullptr) {
			throw std::logic_error("ERROR: mapped data is null, but mapped range size is larger than 0");
		}

		// ensure the portion to copy is within the mapped region.
		typename BASE::range_type target =
				BASE::range_type::intersect(mapped_range, range_bytes);

		if (target.size() == 0) {
			// print error through exception.
			std::cout << "WARNING: read_range: requested " << range_bytes << " not in mapped " << mapped_range << std::endl;
			return target;
		}

		// resize output's capacity
		if (output.capacity() < target.size()) output.resize(target.size());

		// copy the data into memory.  vector is contiguous, so this is okay.
		memmove(output.data(), md_data + (target.start - mapped_range.start), target.size());

		return target;

		// unmap when go out of scope.
	}

	inline mapped_data map(typename BASE::range_type const & range_bytes) {
	  return mapped_data(this->fd, BASE::range_type::intersect(range_bytes, this->file_range_bytes));
	}

	/**
	 * initializes a file for reading/writing via memmap
	 * @param _filename 	name of file to open
	 */
	mmap_file(std::string const & _filename) :
		::bliss::io::base_file(_filename) {

//    // for testing:  multiple processes on the same node maps to the same ptr address
//    map(this->file_range_bytes);
//    std::cout << " serial fd = " << this->fd << " mapped region = " << mapped_range_bytes << " pointer is " << (const void*)mapped_data << std::endl;
//    unmap();
	};

  /**
   * initializes a file for reading/writing via memmap.  for use by parallel file (composition)
   * @param _filename   name of file to open
   * @param _file_size  previously determined file size.
   */
  mmap_file(std::string const & _filename, size_t const & _file_size, size_t const & delay_ms = 0) :
    ::bliss::io::base_file(_filename, _file_size, delay_ms) {}

  /**
   * initializes a file for reading/writing via memmap.  for use by parallel file (composition)
   * @param _filename   name of file to open
   * @param _file_size  previously determined file size.
   */
  mmap_file(int const & _fd, size_t const & _file_size) :
    ::bliss::io::base_file(_fd, _file_size) {}

	/// default destructor
	virtual ~mmap_file() {
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
	void open_file_stream() {
//		printf("open seq file stream\n");


		fp = fdopen(dup(this->fd), "r");  // make a copy of file descriptor first.  this allows closing file pointer.

		if (fp == nullptr) {
			this->file_range_bytes.end = 0;

      // if open failed, throw exception.
      ::std::stringstream ss;
      int myerr = errno;
      ss << "ERROR in stdio_file open: this->fd "  << this->fd << " error " << myerr << ": " << strerror(myerr);

      throw ::bliss::utils::make_exception<::bliss::io::IOException>(ss.str());

		}
	}

	/// virtual funciton for closing a file
	void close_file_stream() {
		if (fp != nullptr) {
			fclose(fp);    // okay to close since fd is dup'ed
			fp = nullptr;
//			this->fd = -1;
//			this->file_range_bytes.end = 0;
		}
	}



public:

	// this is needed to prevent overload name hiding.  see http://stackoverflow.com/questions/888235/overriding-a-bases-overloaded-function-in-c/888337#888337
	using BASE::read_range;
	using BASE::size;

	/**
	 * @brief  bulk load all the data and return it in a newly constructed vector.  reuse vector
	 * @note	virtual so different file reading mechanisms can be defined AND parallel
	 * @param range_bytes	range to read, in bytes
	 * @param output		vector containing data as bytes.
	 */
	virtual range_type read_range(typename ::bliss::io::file_data::container & output, range_type const & range_bytes) {
	  this->open_file_stream();

		// ensure the portion to copy is within the mapped region.
		typename BASE::range_type target =
				BASE::range_type::intersect(this->file_range_bytes, range_bytes);

		if (target.size() == 0) {
			// print error through exception.
			std::cout << "WARNING: read_range: requested " << range_bytes << " not in file " << file_range_bytes << std::endl;
			output.clear();
			return target;
		}

		// use fseek instead lseek to preserve buffering.
		int res = fseeko64(fp, target.start, SEEK_SET);
		if ( res == -1 ) {
			std::stringstream ss;
			int myerr = errno;
			ss << "ERROR: fseeko64: file " << this->filename << " error " << myerr << ": " << strerror(myerr);

			throw ::bliss::utils::make_exception<bliss::io::IOException>(ss.str());
		}

//		std::cout << "curr pos in fd is " << ftell(this->fp) << std::endl;

		// resize output's capacity
		if (output.capacity() < target.size()) output.resize(target.size());

		size_t read = fread_unlocked(output.data(), 1, target.size(), fp);

		if ((read == 0) && (target.size() > 0)) {
			std::stringstream ss;
			int myerr = errno;
			ss << "ERROR: fread: file " << this->filename << " error " << myerr << ": " << strerror(myerr);

			throw ::bliss::utils::make_exception<bliss::io::IOException>(ss.str());
		}

//    std::cout << "curr pos in fd after read is " << ftell(this->fp) << std::endl;

		this->close_file_stream();

//		std::cout << "read " << target << std::endl;
		return target;
	}

	/**
	 * initializes a file for reading/writing
	 * @param _filename 	name of file to open
	 */
	stdio_file(std::string const & _filename) : ::bliss::io::base_file(_filename), fp(nullptr) {};


  /**
   * initializes a file for reading/writing.  for use by a parallel file (composition pattern)
   * @param _filename   name of file to open
   * @param _file_size  previously computed file size.
   */
  stdio_file(std::string const & _filename, size_t const & _file_size, size_t const & delay_ms) :
    ::bliss::io::base_file(_filename, _file_size, delay_ms), fp(nullptr) {};

  /**
   * initializes a file for reading/writing.  for use by a parallel file (composition pattern)
   * @param _fd    previously opened file descriptor.
   * @param _file_size  previously computed file size.
   */
  stdio_file(int const & _fd, size_t const & _file_size) :
    ::bliss::io::base_file(_fd, _file_size), fp(nullptr) {
  };

	/// destructor
	virtual ~stdio_file() {
		this->close_file_stream();
	};

	// this is needed to prevent overload name hiding.  see http://stackoverflow.com/questions/888235/overriding-a-bases-overloaded-function-in-c/888337#888337
	using BASE::read_file;



};


/**
 * @brief    file wrapper that uses c stdio calls.  has buffering.
 */
class posix_file : public ::bliss::io::base_file {


protected:

  using BASE = ::bliss::io::base_file;

public:

  // this is needed to prevent overload name hiding.  see http://stackoverflow.com/questions/888235/overriding-a-bases-overloaded-function-in-c/888337#888337
  using BASE::read_range;

  /**
   * @brief  bulk load all the data and return it in a newly constructed vector.  reuse vector
   * @note  virtual so different file reading mechanisms can be defined AND parallel
   * @param range_bytes range to read, in bytes
   * @param output    vector containing data as bytes.
   * @return  the range for the read data.
   */
  virtual range_type read_range(typename ::bliss::io::file_data::container & output, range_type const & range_bytes) {
    if (this->fd == -1) {
      throw ::bliss::utils::make_exception<std::logic_error>("ERROR: read_range: file pointer is null");
    }

    // ensure the portion to copy is within the mapped region.
    typename BASE::range_type target =
        BASE::range_type::intersect(this->file_range_bytes, range_bytes);

    if (target.size() == 0) {
      // print error through exception.
      std::cout << "WARNING: read_range: requested " << range_bytes << " not in file " << this->file_range_bytes << std::endl;
      output.clear();
      return target;
    }

//    // use fseek instead lseek to preserve buffering.
//    int res = lseek64(this->fd, target.start, SEEK_SET);
//    if ( res == -1 ) {
//      std::stringstream ss;
//      int myerr = errno;
//      ss << "ERROR: lseek64: file " << this->filename << " error " << myerr << ": " << strerror(myerr);
//
//      throw ::bliss::utils::make_exception<bliss::io::IOException>(ss.str());
//    }

    //std::cout << "curr pos in fd is " << lseek64(this->fd, 0, SEEK_CUR) << ::std::endl;

    // resize output's capacity
    if (output.capacity() < target.size()) output.resize(target.size());

    size_t s = 0;
    long count;

    //pread64 can only read 2GB at a time
    for (; s < target.size(); ) {
     	count = pread64(this->fd, output.data() + s, std::min(1UL << 30, target.size() - s ),
    			static_cast<__off64_t>(target.start + s));

        if (count < 0) {
          std::stringstream ss;
          int myerr = errno;
          ss << "ERROR: pread64: file " << this->filename << " error " << myerr << ": " << strerror(myerr);

          throw ::bliss::utils::make_exception<bliss::io::IOException>(ss.str());
        }

    	s += count;
    }


    if (s != target.size()) {
        std::stringstream ss;
        ss << "ERROR: pread64: file " << this->filename << " read " << s << " less than range: " << target.size();
        throw ::bliss::utils::make_exception<bliss::io::IOException>(ss.str());

    }

    return target;
  }

  /**
   * initializes a file for reading/writing
   * @param _filename   name of file to open
   */
  posix_file(std::string const & _filename) : ::bliss::io::base_file(_filename) {};


  /**
   * initializes a file for reading/writing.  for use by a parallel file (composition pattern)
   * @param _filename   name of file to open
   * @param _file_size  previously computed file size.
   */
  posix_file(std::string const & _filename, size_t const & _file_size, size_t const & delay_ms) :
    ::bliss::io::base_file(_filename, _file_size, delay_ms) {};

  /**
   * initializes a file for reading/writing.  for use by a parallel file (composition pattern)
   * @param _fd   previously opened file descriptor
   * @param _file_size  previously computed file size.
   */
  posix_file(int const & _fd, size_t const & _file_size) :
    ::bliss::io::base_file(_fd, _file_size) {};

  /// destructor
  virtual ~posix_file() {};

  // this is needed to prevent overload name hiding.  see http://stackoverflow.com/questions/888235/overriding-a-bases-overloaded-function-in-c/888337#888337
  using BASE::read_file;



};

#ifdef USE_MPI


namespace parallel {




class base_file : public ::bliss::io::base_file {
protected:
	using BASE = ::bliss::io::base_file;

	using range_type = typename BASE::range_type;


	/// communicator used.  object instead of being a reference - lifetime of a src comm temp object in constructor is just that of the constructor call.
	const ::mxx::comm comm;

		/// compute file size in parallel (1 proc, then broadcast.
	size_t get_file_size() {

		size_t file_size = 0;

		// rank 0 gets the size.
		if (comm.rank() == 0) {
//		  printf("parallel base file get_file_size\n");
		  file_size = BASE::get_file_size();
		}
		// then broadcast.
		if (comm.size() > 1)
			MPI_Bcast(&file_size, 1, MPI_UNSIGNED_LONG, 0, comm);

		return file_size;
	}

  base_file(::mxx::comm const & _comm = ::mxx::comm()) :
    ::bliss::io::base_file(static_cast<int>(-1), static_cast<size_t>(0)), // will find real size very soon
     comm(_comm.copy()) {  // _comm could be a temporary constructed from MPI_Comm.
  };


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
		::bliss::io::base_file(static_cast<int>(-1), static_cast<size_t>(0)), // will find real size very soon
		 comm(_comm.copy()) {  // _comm could be a temporary constructed from MPI_Comm.  std::move not needed.  copy elision is in effect.
	  this->filename = _filename;
		this->file_range_bytes.end = this->get_file_size();
		this->BASE::open_file();
	};



	/// destructor
	virtual ~base_file() {};  // will call super's unmap.

	// this is needed to prevent overload name hiding.  see http://stackoverflow.com/questions/888235/overriding-a-bases-overloaded-function-in-c/888337#888337
	using BASE::read_file;

};



class base_shared_fd_file : public ::bliss::io::parallel::base_file {
protected:
  using BASE = ::bliss::io::parallel::base_file;

  using range_type = typename BASE::range_type;


  /// opens a file. side effect computes size of the file.
  void open_file() {
    // first split the communicator


//	if (this->comm.rank() == 0) printf("open parallel file, grouped by node\n");

    ::mxx::comm shared = this->comm.split_shared();

    if (shared.rank() == 0) {
      this->::bliss::io::base_file::open_file();  // fd is populated on shared comm rank 0.
    }

    int id = ::mxx::allreduce(this->comm.rank(), [](int const & x, int const & y){ return ::std::min(x, y); }, shared);

    shared.barrier();

    if (shared.size() > 1) {
      // if we have more than 1, then we broadcast the file id via interprocess communicator, so that
      // all mpi processes on the same node share the same file handle.

      ::bliss::io::util::broadcast_file_descriptor(this->fd, id, shared.size(), shared.rank());

	if (this->comm.rank() == id) 
		printf("broadcasting file descriptor %d from rank G%d N%d to rest.\n", this->fd, id, shared.rank()); 
    }
    // that's it.

    //std::cout << "after copying, curr pos in fd is " << lseek64(this->fd, 0, SEEK_CUR) << std::endl;

    shared.barrier();
  }

public:

  // this is needed to prevent overload name hiding.  see http://stackoverflow.com/questions/888235/overriding-a-bases-overloaded-function-in-c/888337#888337
  using BASE::read_range;

  /**
   * @brief constructor
   * @note  _comm could be a temporary constructed from MPI_Comm.  in which case, the lifetime is the call of the constructor.
   *      to save it, we cannot just save the reference - the underlying object would go out of scope and die.
   *      we have to copy it.  ::mxx::comm does not have copy constructor, but has explict copy function and move constructor.  so use them.
   *      note that copy function is collective.  move constructor is not collective, so it could cause deadlock.
   * @param _filename     name of file to open
   * @param _comm       MPI communicator to use.
   */
  base_shared_fd_file(std::string const & _filename, ::mxx::comm const & _comm = ::mxx::comm()) :
    ::bliss::io::parallel::base_file(_comm) {  // _comm could be a temporary constructed from MPI_Comm.
    this->filename = _filename;
    this->file_range_bytes.end = this->BASE::get_file_size();
    this->open_file();
  };

  /// destructor
  virtual ~base_shared_fd_file() {
    this->close_file();
  };  // will call super's close.

  // this is needed to prevent overload name hiding.  see http://stackoverflow.com/questions/888235/overriding-a-bases-overloaded-function-in-c/888337#888337
  using BASE::read_file;

};



template <typename FileReader,
          typename FileParser = ::bliss::io::BaseFileParser<typename ::bliss::io::file_data::const_iterator>,
          typename BaseType = ::bliss::io::parallel::base_file >
class partitioned_file : public BaseType {

protected:
	using BASE = BaseType;

	using range_type = typename ::bliss::io::base_file::range_type;
	using FileParserType = FileParser;

	/// FileReader
	FileReader reader;

	/// overlap amount
	const size_t overlap;

	/// partitioner to use.
	::bliss::partition::BlockPartitioner<typename BASE::range_type> partitioner;


	/**
	 * @brief partitions the specified range by the number of processes in communicator
	 * @note  does not add overlap.  this is strictly for block partitioning a range.
	 * 			if overlap is needed, use the overload version.
	 * @note  assumes input parameter is in memory and inside file.
	 */
	typename BASE::range_type partition(typename BASE::range_type const & range_bytes) {
	  typename BASE::range_type target =
				BASE::range_type::intersect(range_bytes, this->file_range_bytes);

		if (this->comm.size() == 1) {
			return target;
		}

		partitioner.configure(target, this->comm.size());

		return partitioner.getNext(this->comm.rank());
//		typename BASE::range_type result = partitioner.getNext(this->comm.rank());
//
//		std::cout << "rank = " << this->comm.rank() << " range " << result << std::endl;
//		return result;
	}

	/**
	 * @brief  partition the in mem valid range by the number of processes in communicator
	 * @note	overlap is added to the end of the result ranges.  This may make the result extend
	 * 			beyong initial in_mem_range.
	 * @note	assumes that both input parameters are within file range and in memory.
	 * @return  pair of ranges, first range is the in memory, partitioned, second range is the valid range.
	 *
	 */
	::std::pair<typename BASE::range_type, typename BASE::range_type>
	overlapped_partition(typename BASE::range_type const & in_mem_range_bytes,
	          typename BASE::range_type const & valid_range_bytes) {
//		::std::cout << "rank " << this->comm.rank() << " partition comm_size " << this->comm.size() << ::std::endl;
//		std::cout << " partition comm object at " << &(this->comm) << " for this " << this << std::endl;
//        std::cout << " partition comm is set to world ? " <<  (MPI_COMM_WORLD == this->comm ? "y" : "n") << std::endl;

	  typename BASE::range_type in_mem =
				BASE::range_type::intersect(in_mem_range_bytes, this->file_range_bytes);

		// and restrict valid to in memory data.
		typename BASE::range_type valid =
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
		typename BASE::range_type in_mem_valid = valid;
		in_mem_valid.end += overlap;

		in_mem_valid.intersect(in_mem);

		return ::std::make_pair(in_mem_valid, valid);
	}

public:
	// this is needed to prevent overload name hiding.  see http://stackoverflow.com/questions/888235/overriding-a-bases-overloaded-function-in-c/888337#888337
	using BASE::read_range;

	/**
	 * @brief  bulk load the data and return it in a newly constructed vector.  block decomposes the range.  reuse vector.  no overlap
	 * @note   virtual so other overloads of this function can also block decompose the range.
	 * @param range_bytes	range to read, in bytes
	 * @param output		vector containing data as bytes.
	 */
	virtual typename BASE::range_type read_range(typename ::bliss::io::file_data::container & output,
	                                             typename BASE::range_type const & range_bytes) {
		// first get rough partition
	  typename BASE::range_type target = partition(range_bytes);

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
		BaseType(_filename, _comm),
		 reader(this->fd, this->file_range_bytes.end), overlap(_overlap) {};

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

	  typename BASE::range_type in_mem_partitioned;
	  typename BASE::range_type valid_partitioned;

		::std::tie(in_mem_partitioned, valid_partitioned) =
				overlapped_partition(this->file_range_bytes, this->file_range_bytes);

		// then read the range via sequential version
		output.in_mem_range_bytes = reader.read_range(output.data, in_mem_partitioned);

		output.valid_range_bytes = valid_partitioned;
		output.parent_range_bytes = this->file_range_bytes;
	}
};


template <typename FileReader, typename BaseType>
class partitioned_file<FileReader, ::bliss::io::FASTQParser<typename ::bliss::io::file_data::const_iterator >, BaseType > :
	public BaseType {

protected:
	using BASE = BaseType;
	using FileParserType = ::bliss::io::FASTQParser<typename ::bliss::io::file_data::const_iterator >;
	using range_type = typename ::bliss::io::base_file::range_type;

	/// FileReader
	FileReader reader;

	/// overlap amount
	const size_t overlap;

	/// partitioner to use.
	::bliss::partition::BlockPartitioner<typename BASE::range_type> partitioner;

	/**
	 * @brief partitions the specified range by the number of processes in communicator
	 * @note  does not add overlap.  this is strictly for block partitioning a range.
	 * 			if overlap is needed, use the overload version.
	 * @note  assumes input parameter is in memory and inside file.
	 */
	typename BASE::range_type partition(typename BASE::range_type const & range_bytes) {
	  typename BASE::range_type target =
				BASE::range_type::intersect(range_bytes, this->file_range_bytes);

		if (this->comm.size() == 1) {
			return target;
		}

		partitioner.configure(target, this->comm.size());

		return partitioner.getNext(this->comm.rank());
//		typename BASE::range_type result = partitioner.getNext(this->comm.rank());
//
//		std::cout << "rank = " << this->comm.rank() << " range " << result << std::endl;
//		return result;
	}


	/**
	 * @brief  partition the in mem valid range by the number of processes in communicator
	 * @note	overlap is added to the end of the result ranges.  This may make the result extend
	 * 			beyong initial in_mem_range.
	 * @note	assumes that both input parameters are within file range and in memory.
   * @return  pair of ranges, first range is the in memory, partitioned, second range is the valid range.
	 */
	::std::pair<typename BASE::range_type, typename BASE::range_type>
	overlapped_partition(typename BASE::range_type const & in_mem_range_bytes,
	          typename BASE::range_type const & valid_range_bytes) {
//		::std::cout << "rank " << this->comm.rank() << " partition comm_size " << this->comm.size() << ::std::endl;
//		std::cout << " partition comm object at " << &(this->comm) << " for this " << this << std::endl;
//        std::cout << " partition comm is set to world ? " <<  (MPI_COMM_WORLD == this->comm ? "y" : "n") << std::endl;

	  typename BASE::range_type in_mem =
				BASE::range_type::intersect(in_mem_range_bytes, this->file_range_bytes);

		// and restrict valid to in memory data.
		typename BASE::range_type valid =
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
		typename BASE::range_type in_mem_valid = valid;
		in_mem_valid.end += overlap;

		in_mem_valid.intersect(in_mem);

		return ::std::make_pair(in_mem_valid, valid);
	}

public:
	// this is needed to prevent overload name hiding.  see http://stackoverflow.com/questions/888235/overriding-a-bases-overloaded-function-in-c/888337#888337
	using BASE::read_range;

	/**
	 * @brief  bulk load the data and return it in a newly constructed vector.  block decomposes the range.  reuse vector.  no overlap
	 * @note   virtual so other overloads of this function can also block decompose the range.
	 * @param range_bytes	range to read, in bytes
	 * @param output		vector containing data as bytes.
	 */
	virtual typename BASE::range_type read_range(typename ::bliss::io::file_data::container & output,
	                                             typename BASE::range_type const & range_bytes) {

		// first get rough partition
	  typename BASE::range_type target = partition(range_bytes);


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
		BaseType(_filename, _comm),
		 reader(this->fd, this->file_range_bytes.end), overlap(0UL) {};

	/// destructor
	virtual ~partitioned_file() {};  // will call super's unmap.


	// this is needed to prevent overload name hiding.  see http://stackoverflow.com/questions/888235/overriding-a-bases-overloaded-function-in-c/888337#888337
	using BASE::read_file;

	/**
	 * @brief  read the whole file.  reuse allocated file_data object.
	 * @note	virtual so that parallel readers can compute range to read.
	 * 			same as mpiio file's read_file
	 * @param output 		file_data object containing data and various ranges.
	 */
	virtual void read_file(::bliss::io::file_data & output) {

		// overlap is set to page size, so output will have sufficient space.
		// note that this is same behavior as the serial mmap_file

		// then read the range via sequential version
		range_type partition_range = read_range(output.data, this->file_range_bytes);

//		std::cout << "rank " << this->comm.rank() << " read " << partition_range << std::endl;
//			std::ostream_iterator<unsigned char> oit(std::cout, "");
//			std::copy(output.data.begin(), output.data.begin() + 100, oit);
//			std::cout << std::endl;

		range_type in_mem = partition_range;
		in_mem.end = in_mem.start + output.data.size();

//		std::cout << "rank " << this->comm.rank() << " in mem " << in_mem << std::endl;
//		std::cout << "rank " << this->comm.rank() << " overlap " << this->overlap << std::endl;


		// now search for the true start.
		FileParserType parser;
  // mark the first entry found.
		size_t real_start = parser.init_parser(output.in_mem_cbegin(), this->file_range_bytes,
				in_mem, partition_range, this->comm);


//		std::cout << "rank " << this->comm.rank() << " real start " << real_start << std::endl;
//		std::cout << "rank " << this->comm.rank() << " data size before remove overlap " << output.data.size() << std::endl;

		// now clear the region outside of the block range  (the overlap portion)
		if (output.data.size() > partition_range.size()) {
			output.data.erase(output.data.begin() + partition_range.size(), output.data.end());
		}

//		std::cout << "rank " << this->comm.rank() << " data size after remove overlap " << output.data.size() << std::endl;

		// ==== shift values to the left (as vector?)
		// actually, not shift. need to do all to all - if a rank found no internal starts
		// first compute the target proc id via exscan
		bool not_found = (real_start >= partition_range.end);  // if real start is outside of partition, not found
		real_start = std::min(real_start, partition_range.end);
		int target_rank = not_found ? 0 : this->comm.rank();
		target_rank = ::mxx::exscan(target_rank, [](int const & x, int const & y) {
			return (x < y) ? y : x;
		}, this->comm);

//		std::cout << "rank " << this->comm.rank() << " adjusted real start " << real_start << std::endl;
//		std::cout << "rank " << this->comm.rank() << " target rank " << target_rank << std::endl;


		// MPI 2 does not have neighbor_ collective operations, so we can't create
		// graph with sparse edges.  in any case, the amount of data we send should be
		// small so alltoallv should be enough
		std::vector<size_t> send_counts(this->comm.size(), 0);
		if (this->comm.rank() > 0) send_counts[target_rank] = real_start - in_mem.start;

		// copy the region to shift
		typename ::bliss::io::file_data::container shifted =
				::mxx::all2allv(output.data, send_counts, this->comm);
		output.data.insert(output.data.end(), shifted.begin(), shifted.end());

//		std::cout << "rank " << this->comm.rank() << " shifted " << shifted.size() << std::endl;
//		std::cout << "rank " << this->comm.rank() << " new data size " << output.data.size() << std::endl;

		// adjust the ranges.
		output.in_mem_range_bytes = partition_range;
		output.in_mem_range_bytes.end = partition_range.start + output.data.size();

//		std::cout << "rank " << this->comm.rank() << " final in mem " << output.in_mem_range_bytes << std::endl;


		output.valid_range_bytes.start = real_start;
		output.valid_range_bytes.end =
				not_found ? partition_range.end : output.in_mem_range_bytes.end;

//		std::cout << "rank " << this->comm.rank() << "/" << this->comm.size() << " final valid " << output.valid_range_bytes << std::endl;

		output.parent_range_bytes = this->file_range_bytes;

//		std::cout << "rank " << this->comm.rank() << " file  " << output.parent_range_bytes << std::endl;

	}

};


template <typename FileReader, typename BaseType>
class partitioned_file<FileReader, ::bliss::io::FASTAParser<typename ::bliss::io::file_data::const_iterator>, BaseType > :
	public BaseType {

protected:
	using BASE = BaseType;

	using FileParserType = ::bliss::io::FASTAParser<typename ::bliss::io::file_data::const_iterator>;

	using range_type = typename ::bliss::io::base_file::range_type;

	/// FileReader
	FileReader reader;

	/// overlap amount
	const size_t overlap;

	/// partitioner to use.
	::bliss::partition::BlockPartitioner<typename BASE::range_type> partitioner;

	/**
	 * @brief partitions the specified range by the number of processes in communicator
	 * @note  does not add overlap.  this is strictly for block partitioning a range.
	 * 			if overlap is needed, use the overload version.
	 * @note  assumes input parameter is in memory and inside file.
	 */
	typename BASE::range_type partition(typename BASE::range_type const & range_bytes) {
	  typename BASE::range_type target =
				BASE::range_type::intersect(range_bytes, this->file_range_bytes);

		if (this->comm.size() == 1) {
			return target;
		}

		partitioner.configure(target, this->comm.size());

		return partitioner.getNext(this->comm.rank());
//		typename BASE::range_type result = partitioner.getNext(this->comm.rank());
//
//		std::cout << "rank = " << this->comm.rank() << " range " << result << std::endl;
//		return result;
	}


	/**
	 * @brief  partition the in mem valid range by the number of processes in communicator
	 * @note	overlap is added to the end of the result ranges.  This may make the result extend
	 * 			beyong initial in_mem_range.
	 * @note	assumes that both input parameters are within file range and in memory.
   * @return  pair of ranges, first range is the in memory, partitioned, second range is the valid range.
	 */
	::std::pair<typename BASE::range_type, typename BASE::range_type>
	overlapped_partition(typename BASE::range_type const & in_mem_range_bytes,
	          typename BASE::range_type const & valid_range_bytes) {
//		::std::cout << "rank " << this->comm.rank() << " partition comm_size " << this->comm.size() << ::std::endl;
//		std::cout << " partition comm object at " << &(this->comm) << " for this " << this << std::endl;
//        std::cout << " partition comm is set to world ? " <<  (MPI_COMM_WORLD == this->comm ? "y" : "n") << std::endl;

	  typename BASE::range_type in_mem =
				BASE::range_type::intersect(in_mem_range_bytes, this->file_range_bytes);

		// and restrict valid to in memory data.
		typename BASE::range_type valid =
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
		typename BASE::range_type in_mem_valid = valid;
		in_mem_valid.end += 2 * overlap;

		in_mem_valid.intersect(in_mem);

		return ::std::make_pair(in_mem_valid, valid);
	}

public:
	// this is needed to prevent overload name hiding.  see http://stackoverflow.com/questions/888235/overriding-a-bases-overloaded-function-in-c/888337#888337
	using BASE::read_range;

	/**
	 * @brief  bulk load the data and return it in a newly constructed vector.  block decomposes the range.  reuse vector.  no overlap
	 * @note   virtual so other overloads of this function can also block decompose the range.
	 * @param range_bytes	range to read, in bytes
	 * @param output		vector containing data as bytes.
	 */
	virtual typename BASE::range_type read_range(typename ::bliss::io::file_data::container & output,
	                                             typename BASE::range_type const & range_bytes) {

		// first get rough partition
	  typename BASE::range_type target = partition(range_bytes);


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
		BaseType(_filename, _comm),
		 reader(this->fd, this->file_range_bytes.end), overlap(_overlap) {};

	/// destructor
	virtual ~partitioned_file() {};  // will call super's unmap.


	// this is needed to prevent overload name hiding.  see http://stackoverflow.com/questions/888235/overriding-a-bases-overloaded-function-in-c/888337#888337
	using BASE::read_file;

	/**
	 * @brief  read the whole file.  reuse allocated file_data object.
	 * @note	virtual so that parallel readers can compute range to read.
	 * 			same as mpiio file's read_file
	 * @param output 		file_data object containing data and various ranges.
	 */
	virtual void read_file(::bliss::io::file_data & output) {

		// overlap is set to page size, so output will have sufficient space.
		// note that this is same behavior as the serial mmap_file

		// then read the range via sequential version
		::std::tie(output.in_mem_range_bytes, output.valid_range_bytes) =
				overlapped_partition(this->file_range_bytes, this->file_range_bytes);

		// then read the range via sequential version
		output.in_mem_range_bytes = reader.read_range(output.data, output.in_mem_range_bytes);

		output.parent_range_bytes = this->file_range_bytes;

//		std::cout << "rank " << this->comm.rank() << " read " << partition_range << std::endl;
//			std::ostream_iterator<unsigned char> oit(std::cout, "");
//			std::copy(output.data.begin(), output.data.begin() + 100, oit);
//			std::cout << std::endl;

//		std::cout << "rank " << this->comm.rank() << " in mem " << in_mem << std::endl;
//		std::cout << "rank " << this->comm.rank() << " overlap " << this->overlap << std::endl;

		// now search for the true start.
		FileParserType parser;
  // mark the first entry found.
		size_t overlap_end = parser.find_overlap_end(output.in_mem_cbegin(), output.parent_range_bytes,
				output.in_mem_range_bytes, output.valid_range_bytes.end, overlap);

		// erase the extra.
		output.in_mem_range_bytes.end = overlap_end;
		output.data.erase(output.data.begin() + output.in_mem_range_bytes.size(), output.data.end());


//		std::cout << "rank " << this->comm.rank() << "/" << this->comm.size() << " final valid " << output.valid_range_bytes << std::endl;

//		std::cout << "rank " << this->comm.rank() << " file  " << output.parent_range_bytes << std::endl;

	}

};

/// disable shared fd for stdio file.
template <typename FileParser>
class partitioned_file<::bliss::io::stdio_file, FileParser, ::bliss::io::parallel::base_shared_fd_file >
{
  private:
	template <typename dummy = FileParser>
    partitioned_file(std::string const & _filename, size_t const & _overlap = 0UL, ::mxx::comm const & _comm = ::mxx::comm()) {
		// static_assert(false, "ERROR: stdio_file is not compatible with base_shared_fd_file, as there is no fread that will maintain the old position");
    }

	template <typename dummy = FileParser>
	partitioned_file() {
		// static_assert(false, "ERROR: stdio_file is not compatible with base_shared_fd_file, as there is no fread that will maintain the old position");
	};

	~partitioned_file() {};
};



// multilevel parallel file io relies on MPIIO.
template <typename FileParser = ::bliss::io::BaseFileParser<typename ::bliss::io::file_data::const_iterator> >
class mpiio_base_file : public ::bliss::io::base_file {

	static_assert(sizeof(MPI_Offset) > 4, "ERROR: MPI_Offset is defined as an integer 4 bytes or less.  Do not use mpiio_file ");

protected:
	using BASE = ::bliss::io::base_file;

	/// overlap amount
	const size_t overlap;

  /// communicator used.  object instead of being a reference - lifetime of a src comm temp object in constructor is just that of the constructor call.
	const ::mxx::comm comm;

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

//	  if (comm.rank() == 0) printf("mpiio base file get_file_size\n");

		if (fh == MPI_FILE_NULL) {
			std::stringstream ss;
			ss << "ERROR in mpiio: rank " << comm.rank() << " file " << this->filename << " not yet open " << std::endl;

			throw ::bliss::utils::make_exception<::bliss::io::IOException>(ss.str());
		}

		MPI_Offset s;

		int res = MPI_File_get_size(fh, &s);
		if (res != MPI_SUCCESS) {
			throw ::bliss::utils::make_exception<::bliss::io::IOException>(get_error_string("get_size", res));
		}

		return static_cast<size_t>(s);
	}

	/// opens a file. side effect computes size of the file.
	void open_file() {
//	  if (this->comm.rank() == 0) printf("open mpiio file\n");
	

	// first clear previously open file
		close_file();

		// open the file
		int res = MPI_File_open(this->comm, const_cast<char *>(this->filename.c_str()), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);

		if (res != MPI_SUCCESS) {
			throw ::bliss::utils::make_exception<::bliss::io::IOException>(get_error_string("open", res));
		}

		// ensure atomicity is turned off
		MPI_File_set_atomicity(fh, 0);
	}

	/// funciton for closing a file
	void close_file() {
		if (fh != MPI_FILE_NULL) {
			int res = MPI_File_close(&fh);
			if (res != MPI_SUCCESS) {
				throw ::bliss::utils::make_exception<::bliss::io::IOException>(get_error_string("close", res));
			}
			fh = MPI_FILE_NULL;
		}
	}



public:
	// this is needed to prevent overload name hiding.  see http://stackoverflow.com/questions/888235/overriding-a-bases-overloaded-function-in-c/888337#888337
	using BASE::read_range;
	using FileParserType = FileParser;

	/**
	 * @brief  bulk load the data and return it in a newly constructed vector.  block decomposes the range.  reuse vector
	 * @note   virtual so other overloads of this function can also block decompose the range.
	 * @param range_bytes	range to read, in bytes
	 * @param output		vector containing data as bytes.
	 */
	virtual range_type read_range(typename ::bliss::io::file_data::container & output,
				range_type const & range_bytes) {
		if (fh == MPI_FILE_NULL) {
			std::stringstream ss;
			ss << "ERROR in mpiio: rank " << comm.rank() << " file " << this->filename << " not yet open " << std::endl;

			throw ::bliss::utils::make_exception<::bliss::io::IOException>(ss.str());
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
		read_range.end += std::is_same<FileParser, ::bliss::io::FASTAParser<typename bliss::io::file_data::const_iterator> >::value ? 2 * this->overlap : this->overlap;
		read_range.intersect(this->file_range_bytes);

		output.resize(read_range.size());    // set size for reading.

		// then read the file using mpiio
		MPI_Status stat;

		// NOTE:  file offset is in units of byte (up to size_t).  number of elements to read has type int.

		size_t step_size = 1UL << 30 ;  // using 2^30 vs 2^31-2^15 does not make a huge performance difference (will only matter for low proc count anyways)
		size_t rem = read_range.size() % step_size;
//	size_t steps = read_range.size() / step_size;
		int count = 0;
		int res = MPI_SUCCESS;


// ======= DOES NOT WORK. number of elements has type int.
//		res =MPI_File_read_at_all(fh, read_range.start, output.data(), read_range.size(), MPI_BYTE, &stat);
//		if (res != MPI_SUCCESS) throw ::bliss::utils::make_exception<::bliss::io::IOException>(get_error_string("read", res, stat));


// ======= DOES NOT WORK to read rem first as BYTES, then read rest as 1GB blocks.  rem reads okay, but the first byte for big type fails.
//		res = MPI_File_read_at_all(fh, read_range.start, output.data(), rem, MPI_BYTE, &stat);
//		if (res != MPI_SUCCESS) throw ::bliss::utils::make_exception<::bliss::io::IOException>(get_error_string("read", res, stat));
//		// count is okay here.
//		res = MPI_Get_count(&stat, MPI_BYTE, &count);
//		if (res != MPI_SUCCESS) throw ::bliss::utils::make_exception<::bliss::io::IOException>(get_error_string("count", res));
//		if (static_cast<size_t>(count) != rem) {
//			std::stringstream ss;
//			ss << "ERROR in mpiio: rank " << comm.rank() << " remainder read error. request " << rem << " bytes got " << count << " bytes" << std::endl;
//			throw ::bliss::utils::make_exception<::bliss::io::IOException>(ss.str());
//		}
//
//		// now read the big data type (size of 1GB)
//		if (steps > 0) {
//			::mxx::datatype dt = ::mxx::get_datatype<unsigned char>().contiguous(step_size);  // make element  2^30 in size.
//		    res = MPI_File_read_at_all(fh, read_range.start + rem, output.data() + rem,
//	                         steps, dt.type(), &stat);
//			if (res != MPI_SUCCESS) throw ::bliss::utils::make_exception<::bliss::io::IOException>(get_error_string("read", res, stat));
//			// count comes back as -32677.
//			res = MPI_Get_count(&stat, dt.type(), &count);
//			if (res != MPI_SUCCESS) throw ::bliss::utils::make_exception<::bliss::io::IOException>(get_error_string("read", res));
//			if (static_cast<size_t>(count) != steps) {
//				std::stringstream ss;
//				ss << "ERROR in mpiio: rank " << comm.rank() << " remainder read error. request " << steps << " 2^30 byte blocks got " << count << " blocks" << std::endl;
//				throw ::bliss::utils::make_exception<::bliss::io::IOException>(ss.str());
//			}
//		}


// ============== does not work to read big data type first then rem.
//		if (steps > 0) {
//			// first blocks failed - reads a vector of 0's
//			::mxx::datatype dt = ::mxx::get_datatype<unsigned char>().contiguous(step_size);  // make element  2^30 in size.
//		    res = MPI_File_read_at_all(fh, read_range.start, output.data(),
//	                         steps, dt.type(), &stat);
//	   	    if (res != MPI_SUCCESS) throw ::bliss::utils::make_exception<::bliss::io::IOException>(get_error_string("read", res, stat));
//			// get count got -32766.
//		    res = MPI_Get_count(&stat, dt.type(), &count);
//        		if (res != MPI_SUCCESS) throw ::bliss::utils::make_exception<::bliss::io::IOException>(get_error_string("read", res));
//			if (static_cast<size_t>(count) != steps) {
//				std::stringstream ss;
//				ss << "ERROR in mpiio: rank " << comm.rank() << " remainder read error. request " << steps << " 2^30 byte blocks got " << count << " blocks" << std::endl;
//				throw ::bliss::utils::make_exception<::bliss::io::IOException>(ss.str());
//			}
//		}
//		res = MPI_File_read_at_all(fh, read_range.start + steps * step_size, output.data() + steps * step_size, rem, MPI_BYTE, &stat);
//		if (res != MPI_SUCCESS) throw ::bliss::utils::make_exception<::bliss::io::IOException>(get_error_string("read", res, stat));

//		res = MPI_Get_count(&stat, MPI_BYTE, &count);
//		if (res != MPI_SUCCESS) throw ::bliss::utils::make_exception<::bliss::io::IOException>(get_error_string("read", res));
//		if (static_cast<size_t>(count) != rem) {
//			std::stringstream ss;
//			ss << "ERROR in mpiio: rank " << comm.rank() << " remainder read error. request " << rem << " bytes got " << count << " bytes" << std::endl;
//			throw ::bliss::utils::make_exception<::bliss::io::IOException>(ss.str());
//		}


// ========= Single big data type that is the whole size does not work.  hangs when nprocs = 1 and for nprocs = 2, does not appear to populate data at all
//		::mxx::datatype dt = ::mxx::get_datatype<unsigned char>().contiguous(read_range.size());  // make element the whole range in size.
//		res = MPI_File_read_at_all(fh, read_range.start, output.data(),  1, dt.type(), &stat);
//		if (res != MPI_SUCCESS) throw ::bliss::utils::make_exception<::bliss::io::IOException>(get_error_string("read", res, stat));
//
// 		// getCount returns -32766, so cannot be used.
//		res = MPI_Get_count(&stat, dt.type(), &count);
//		if (res != MPI_SUCCESS) throw ::bliss::utils::make_exception<::bliss::io::IOException>(get_error_string("read", res));
//		if (static_cast<size_t>(count) != read_range.size()) {
//			std::stringstream ss;
//			ss << "ERROR in mpiio: rank " << comm.rank() << " remainder read error. request " << read_range.size() << " bytes got " << count << " bytes" << std::endl;
//			throw ::bliss::utils::make_exception<::bliss::io::IOException>(ss.str());
//		}

// ===========  iterative works
		size_t iter_step_size;
		// since collective, first get the maximum read size.
		size_t max_read_size = read_range.size();
		max_read_size = ::mxx::allreduce(max_read_size, [](size_t const & x, size_t const & y){
		  return ::std::max(x, y);
		}, this->comm);

		// compute the steps for from the max and local_steps
		size_t steps = (max_read_size + step_size - 1) >> 30;
    size_t local_full_steps = read_range.size() >> 30;
		size_t offset = 0;

		for (size_t s = 0; s < steps; ++s) {
		  // compute the iter_step_size.
		  // below local_full_steps, full step size.
		  // at local full steps (== last step), read_range.size() % step_size;
		  // above local full steps, 0.
		  iter_step_size = (s < local_full_steps) ? step_size :
		      (s == local_full_steps) ? rem : 0;

			res = MPI_File_read_at_all(fh, read_range.start + offset, output.data() + offset,
			                           iter_step_size, MPI_BYTE, &stat);
			offset += iter_step_size;

			if (res != MPI_SUCCESS)
			  throw ::bliss::utils::make_exception<::bliss::io::IOException>(get_error_string("read", res, stat));

			res = MPI_Get_count(&stat, MPI_BYTE, &count);
			if (res != MPI_SUCCESS)
			  throw ::bliss::utils::make_exception<::bliss::io::IOException>(get_error_string("read count", res, stat));

			if (static_cast<size_t>(count) != iter_step_size) {
				std::stringstream ss;
				ss << "ERROR in mpiio: rank " << comm.rank() << " remainder read error. request " << iter_step_size << " bytes got " << count << " bytes" << std::endl;

				throw ::bliss::utils::make_exception<::bliss::io::IOException>(ss.str());
			}

		}

//		std::cout << "rank " << comm.rank() << " done reading " << read_range << std::endl;

		return target;
	}


	mpiio_base_file(::std::string const & _filename, size_t const _overlap = 0UL,  ::mxx::comm const & _comm = ::mxx::comm()) :
	  BASE(static_cast<int>(-1), static_cast<size_t>(0)),
	 	 overlap(_overlap),
				  comm(_comm.copy()), fh(MPI_FILE_NULL) {
		this->filename = _filename;
	  this->open_file();
		this->file_range_bytes.end = this->get_file_size();  // call after opening file
	};

	~mpiio_base_file() { this->close_file(); };

	// this is needed to prevent overload name hiding.  see http://stackoverflow.com/questions/888235/overriding-a-bases-overloaded-function-in-c/888337#888337
	using BASE::read_file;


};


/**
 * parallel file abstraction for MPIIO, block partition with overlap
 */
template <typename FileParser = ::bliss::io::BaseFileParser<typename ::bliss::io::file_data::const_iterator> >
class mpiio_file : public ::bliss::io::parallel::mpiio_base_file<FileParser> {

	protected:
		using BASE = ::bliss::io::parallel::mpiio_base_file<FileParser>;


public:
		// this is needed to prevent overload name hiding.  see http://stackoverflow.com/questions/888235/overriding-a-bases-overloaded-function-in-c/888337#888337
		using BASE::read_range;
		using FileParserType = typename BASE::FileParserType;


		mpiio_file(::std::string const & _filename, size_t const & _overlap = 0UL, ::mxx::comm const & _comm = ::mxx::comm()) :
			BASE(_filename, _overlap, _comm) {};

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

			output.in_mem_range_bytes.end += this->overlap;
			output.in_mem_range_bytes.intersect(this->file_range_bytes);

			output.parent_range_bytes = this->file_range_bytes;

		}
};

/**
 * parallel file abstraction for MPIIO.  FASTQParser, which searches the input.
 */
template <>
class mpiio_file<::bliss::io::FASTQParser<typename ::bliss::io::file_data::const_iterator> > :
	public ::bliss::io::parallel::mpiio_base_file<::bliss::io::FASTQParser<typename ::bliss::io::file_data::const_iterator> > {

	protected:
		using BASE = ::bliss::io::parallel::mpiio_base_file<::bliss::io::FASTQParser<typename ::bliss::io::file_data::const_iterator> > ;


public:
		// this is needed to prevent overload name hiding.  see http://stackoverflow.com/questions/888235/overriding-a-bases-overloaded-function-in-c/888337#888337
		using BASE::read_range;
		using FileParserType = typename BASE::FileParserType;


		mpiio_file(::std::string const & _filename, size_t const & _overlap = 0UL, ::mxx::comm const & _comm = ::mxx::comm()) :
			BASE(_filename, 0UL, _comm) {};      // specify 1 page worth as overlap

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
			range_type partition_range = read_range(output.data, this->file_range_bytes);

//			std::cout << "rank " << this->comm.rank() << " read " << partition_range << std::endl;
//			std::ostream_iterator<unsigned char> oit(std::cout, "");
//			std::copy(output.data.begin(), output.data.begin() + 100, oit);
//			std::cout << std::endl;

			range_type in_mem = partition_range;
			in_mem.end = in_mem.start + output.data.size();

//			std::cout << "rank " << this->comm.rank() << " in mem " << in_mem << std::endl;
//			std::cout << "rank " << this->comm.rank() << " overlap " << this->overlap << std::endl;


			// now search for the true start.
			FileParserType parser;
      // mark the first entry found.
			size_t real_start = parser.init_parser(output.in_mem_cbegin(), this->file_range_bytes,
					in_mem, partition_range, this->comm);


//			std::cout << "rank " << this->comm.rank() << " real start " << real_start << std::endl;
//			std::cout << "rank " << this->comm.rank() << " data size before remove overlap " << output.data.size() << std::endl;

			// now clear the region outside of the block range  (the overlap portion)
			if (output.data.size() > partition_range.size()) {
				output.data.erase(output.data.begin() + partition_range.size(), output.data.end());
			}

//			std::cout << "rank " << this->comm.rank() << " data size after remove overlap " << output.data.size() << std::endl;

			// ==== shift values to the left (as vector?)
			// actually, not shift. need to do all to all - if a rank found no internal starts
			// first compute the target proc id via exscan
			bool not_found = (real_start >= partition_range.end);  // if real start is outside of partition, not found
			real_start = std::min(real_start, partition_range.end);
			int target_rank = not_found ? 0 : this->comm.rank();
			target_rank = ::mxx::exscan(target_rank, [](int const & x, int const & y) {
				return (x < y) ? y : x;
			}, this->comm);

//			std::cout << "rank " << this->comm.rank() << " adjusted real start " << real_start << std::endl;
//			std::cout << "rank " << this->comm.rank() << " target rank " << target_rank << std::endl;


			// MPI 2 does not have neighbor_ collective operations, so we can't create
			// graph with sparse edges.  in any case, the amount of data we send should be
			// small so alltoallv should be enough
			std::vector<size_t> send_counts(this->comm.size(), 0);
			if (this->comm.rank() > 0) send_counts[target_rank] = real_start - in_mem.start;

			// copy the region to shift
			typename ::bliss::io::file_data::container shifted =
					::mxx::all2allv(output.data, send_counts, this->comm);
			output.data.insert(output.data.end(), shifted.begin(), shifted.end());

//			std::cout << "rank " << this->comm.rank() << " shifted " << shifted.size() << std::endl;
//			std::cout << "rank " << this->comm.rank() << " new data size " << output.data.size() << std::endl;

			// adjust the ranges.
			output.in_mem_range_bytes = partition_range;
			output.in_mem_range_bytes.end = partition_range.start + output.data.size();

//			std::cout << "rank " << this->comm.rank() << " final in mem " << output.in_mem_range_bytes << std::endl;


			output.valid_range_bytes.start = real_start;
			output.valid_range_bytes.end =
					not_found ? partition_range.end : output.in_mem_range_bytes.end;

//			std::cout << "rank " << this->comm.rank() << "/" << this->comm.size() << " final valid " << output.valid_range_bytes << std::endl;

			output.parent_range_bytes = this->file_range_bytes;

//			std::cout << "rank " << this->comm.rank() << " file  " << output.parent_range_bytes << std::endl;

		}
};

/**
 * parallel file abstraction for MPIIO.  FASTQParser, which searches the input.
 */
template <>
class mpiio_file<::bliss::io::FASTAParser<typename ::bliss::io::file_data::const_iterator> > :
	public ::bliss::io::parallel::mpiio_base_file<::bliss::io::FASTAParser<typename ::bliss::io::file_data::const_iterator> > {

	protected:
		using BASE = ::bliss::io::parallel::mpiio_base_file<::bliss::io::FASTAParser<typename ::bliss::io::file_data::const_iterator> > ;


public:
		// this is needed to prevent overload name hiding.  see http://stackoverflow.com/questions/888235/overriding-a-bases-overloaded-function-in-c/888337#888337
		using BASE::read_range;
		using FileParserType = typename BASE::FileParserType;


		mpiio_file(::std::string const & _filename, size_t const & _overlap = 0UL, ::mxx::comm const & _comm = ::mxx::comm()) :
			BASE(_filename, _overlap, _comm) {};      // specify 1 page worth as overlap

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
			// then read the range via sequential version
			output.valid_range_bytes = read_range(output.data, this->file_range_bytes);

			// overlap is part of the read so there is no left shift needed.

			// set up the ranges.
			output.in_mem_range_bytes = output.valid_range_bytes;

			output.in_mem_range_bytes.end += 2 * this->overlap;
			output.in_mem_range_bytes.intersect(this->file_range_bytes);

			output.parent_range_bytes = this->file_range_bytes;

	//		std::cout << "rank " << this->comm.rank() << " read " << partition_range << std::endl;
	//			std::ostream_iterator<unsigned char> oit(std::cout, "");
	//			std::copy(output.data.begin(), output.data.begin() + 100, oit);
	//			std::cout << std::endl;

	//		std::cout << "rank " << this->comm.rank() << " in mem " << in_mem << std::endl;
	//		std::cout << "rank " << this->comm.rank() << " overlap " << this->overlap << std::endl;

			// now search for the true start.
			FileParserType parser;
	  // mark the first entry found.
			size_t overlap_end = parser.find_overlap_end(output.in_mem_cbegin(), output.parent_range_bytes,
					output.in_mem_range_bytes, output.valid_range_bytes.end, overlap);

			// erase the extra.
			output.in_mem_range_bytes.end = overlap_end;
			output.data.erase(output.data.begin() + output.in_mem_range_bytes.size(), output.data.end());


	//		std::cout << "rank " << this->comm.rank() << "/" << this->comm.size() << " final valid " << output.valid_range_bytes << std::endl;

	//		std::cout << "rank " << this->comm.rank() << " file  " << output.parent_range_bytes << std::endl;

		}
};


}  // namespace parallel
#endif  // USE_MPI


} // io

} // bliss



#endif /* FILE_HPP_ */
