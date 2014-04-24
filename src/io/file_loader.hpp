/**
 * file_loader.hpp
 *
 *  Created on: Feb 18, 2014
 *      Author: tpan
 *
 *
 *      opens a file via mmap, and optionally loads the entire content to memory.
 *      similar to a container, can return size, begin iterator and end iterator.
 *
 *
 *  this class should do the partitioning across multiple processors using MPI.
 */

#ifndef FILE_LOADER_HPP_
#define FILE_LOADER_HPP_

#include "config.hpp"

#if defined(USE_MPI)
#include "mpi.h"
#endif

#include <string>
#include <cstring>      // memcpy
#include <exception>    // ioexception
#include <sstream>      // stringstream

#include <unistd.h>     // sysconf
#include <sys/stat.h>   // block size.
#include <sys/mman.h>   // mmap
#include <fcntl.h>      // for open

#include "sys/sysinfo.h"  // for meminfo

#include "iterators/range.hpp"
#include "utils/logging.h"

namespace bliss
{
  namespace io
  {

    class io_exception : public std::exception
    {
      protected:
        std::string message;

      public:
        io_exception(const char* _msg)
        {
          message = _msg;
        }

        io_exception(const std::string &_msg)
        {
          message = _msg;
        }

        virtual const char* what() const throw ()
        {
          return message.c_str();
        }
    };

    /**
     *
     */
    template <typename T>
    class file_loader
    {
      public:

        typedef bliss::iterator::range<size_t> RangeType;

        /**
         * opens the file and save file handle.  specifies the range to load.
         *
         * does not mmap the file.  please call load() before running
         * use the MPI_Comm object to inform all nodes about file size and to generate approximate partitions.
         *
         * Adjustment to partitions should be done outside, and then the loader's range object is re-set.
         *  e.g. copy the fullRange object to range, or change the overlap, or align to sequence boundaries.
         *
         *  NOTE: range's overlap is included in range.end.
         *
         * @param _filename       input file name
         * @param _comm           MPI Communicator (defined only if CMake is configured with MPI).
         *                        default is MPI_COMM_WORLD, which means all nodes participate.
         *                        to get just the calling node, use MPI_COMM_SELF.
         */
#if defined(USE_MPI)
        file_loader(const std::string &_filename,
                    const MPI_Comm& _comm = MPI_COMM_WORLD ) throw (io_exception)
            : filename(_filename), range(), fullRange(), file_handle(-1),
              data(nullptr), aligned_data(nullptr), loaded(false),
              nprocs(1), rank(0), comm(_comm)
        {
          // get the processor rank and nprocessors.
          MPI_Comm_rank(comm, &rank);
          MPI_Comm_size(comm, &nprocs);
#else
        file_loader(const std::string &_filename) throw (io_exception)
            : filename(_filename), range(), fullRange(), file_handle(-1),
              data(nullptr), aligned_data(nullptr), loaded(false),
              nprocs(1), rank(0)
        {
#endif
          assert(filename.length() > 0);
          page_size = sysconf(_SC_PAGE_SIZE);

          /// get the file size.
          uint64_t file_size = 0;
          if (rank == 0)
          {
            struct stat filestat;
            stat(filename.c_str(), &filestat);
            file_size = static_cast<uint64_t>(filestat.st_size);
            INFO("block size is " << filestat.st_blksize);
            INFO("sysconf block size is " << page_size);
          }

#if defined(USE_MPI)
          if (nprocs > 1) {
            /// broadcast filesize to all
            MPI_Bcast(&file_size, 1, MPI_UNSIGNED_LONG_LONG, 0, comm);
          }
#endif


          // compute the full range
          fullRange.block_start = 0;
          fullRange.start = 0;
          fullRange.end = file_size / sizeof(T);
          fullRange.overlap = 0;

          // partition the full range
          range = fullRange.block_partition(nprocs, rank);
        }


        /**
         * closes the file.
         */
        virtual ~file_loader()
        {
          unload();
        }

        /**
         * return copy of pointer to the start of the actual data.  not the page aligned one.
         *
         * use this like an iterator.
         */
        inline T* begin()
        {
          return data;
        }
        inline const T* begin() const
        {
          return this->begin();
        }
        /**
         * return copy of pointer to the end of the actual data. (after last element)
         *
         * use this like an iterator.
         */
        inline T* end()
        {
          return data + this->size();
        }
        inline const T* end() const
        {
          return this->end();
        }


        /**
         * recall that range.end includes the overlap
         *
         * @return    size of the mmap region (or preloaded data) in units of T.
         */
        inline size_t size() {
          return range.end - range.start;
        }

        /**
         * return the full range for this file mmap.  (in units of data type T)
         * @return
         */
        const RangeType& getFullRange() const {
          return fullRange;
        }

        /**
         * return the range for this file mmap.  (in units of data type T)
         * @return
         */
        const RangeType& getRange() const {
          return range;
        }


        void setRange(const RangeType& _range)
        {
          range = _range;
        }


      protected:
        std::string filename;
        RangeType range;  // offset in file from where to read
        RangeType fullRange;  // offset in file from where to read

        size_t page_size;
        int file_handle;  // file handle

        T* data;
        T* aligned_data;  // memmapped data, page aligned.  strictly internal

        bool loaded;

        int nprocs;
        int rank;

#if defined(USE_MPI)
        MPI_Comm comm;
#endif

        /**
         * TODO: test on remotely mounted file system.
         */


        /**
         *  performs memmap, and optionally preload the data into memory.
         *
         * @param _memUseFraction value between 0.0 and 0.5, representing percentage of FREE memory to be used.  0 means no preloading.
         */
        void load(const float _memUseFraction = 0.0f) throw (io_exception)
        {
          /// configure the range.
          range = range.align_to_page(page_size);
          assert(range.end > range.start);



          /// open the file and get a handle.
          file_handle = open(filename.c_str(), O_RDONLY);
          if (file_handle == -1)
          {
            std::stringstream ss;
            int myerr = errno;
            ss << "ERROR in file open: ["  << filename << "] error " << myerr << ": " << strerror(myerr);
            throw io_exception(ss.str());
          }


          /// do the mem map
          // NOT using MAP_POPULATE.  it slows things done when testing on single node.
          aligned_data = (T*)mmap(nullptr, (range.end - range.block_start ) * sizeof(T),
                                     PROT_READ,
                                     MAP_PRIVATE, file_handle,
                                     range.block_start * sizeof(T));

          if (aligned_data == MAP_FAILED)
          {

            if (file_handle != -1)
            {
              close(file_handle);
              file_handle = -1;
            }

            std::stringstream ss;
            int myerr = errno;
            ss << "ERROR in mmap: " << myerr << ": " << strerror(myerr);
            throw io_exception(ss.str());
          }


          /// check if can load the region into memory.
          float memUseFraction = _memUseFraction;

          if (_memUseFraction > 0.5f) memUseFraction = 0.5f;
          if (_memUseFraction < 0.0f) memUseFraction = 0.0f;


          loaded = false;
          if (memUseFraction > 0.0f) {
            /// check if we can preload.
            struct sysinfo memInfo;
            sysinfo (&memInfo);
            long long limit = static_cast<long long>(static_cast<float>(memInfo.freeram * memInfo.mem_unit) * memUseFraction);

            if (limit < (this->size() * sizeof(T)) ) {
              loaded = true;
            } else {
              WARNING("Insufficient memory requested during file loading.  Not preloading file.");
            }
          }

          if (loaded)
          {
            // allocate space
            data = new T[this->size() + 1];
            //copy data over
            memcpy(data, aligned_data + (range.start - range.block_start),
                   sizeof(T) * this->size());

            data[this->size()] = 0;

            // close the input
            munmap(aligned_data, (range.end - range.block_start) * sizeof(T));

            aligned_data = data;
          }
          else
          {

            data = aligned_data + (range.start - range.block_start);

          }
        }

        void unload()
        {
          if (loaded)
          {
            if (data != nullptr)
              delete [] data;
            loaded = false;
          }
          else
          {
            if (aligned_data != nullptr || aligned_data != MAP_FAILED)
              munmap(aligned_data, (range.end - range.block_start) * sizeof(T));
          }
          aligned_data = nullptr;
          data = nullptr;

          //printf("unloading complete.\n");
          if (file_handle != -1) {
            close(file_handle);
            file_handle = -1;
          }
        }


    };

  } /* namespace functional */
} /* namespace bliss */
#endif /* FILE_LOADER_HPP_ */
