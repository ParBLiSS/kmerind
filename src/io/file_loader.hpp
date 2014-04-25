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
     *  real data:  mmap is better for large files and limited memory.
     *              preloading is better for smaller files and/or large amount of memory.
     *  stream processing means data does not need to be buffered in memory - more efficient.
     *
     *
     *  file access:  better to work with a few pages at a time, or to work with large block?
     *
     *  Usage:
     *    instantiate file_loader()  - default partitions the file equally
     *    call adjustPartition()  (optionally)
     *    then call load()
     *
     *    do some work.  e.g. via begin(), end(),
     *      or via getNextChunk()
     *
     *    call unload()
     *    destroy file_loader
     *
     */
    template <typename T, typename SizeT = size_t>
    class file_loader
    {
      public:

        typedef SizeT                             SizeType;
        typedef bliss::iterator::range<SizeType>  RangeType;
        typedef T*                                IteratorType;

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
         *        range is in units of T
         *
         * TODO: test on remotely mounted file system.
         *
         * @param _filename       input file name
         * @param _comm           MPI Communicator (defined only if CMake is configured with MPI).
         *                        default is MPI_COMM_WORLD, which means all nodes participate.
         *                        to get just the calling node, use MPI_COMM_SELF.
         */
#if defined(USE_MPI)
        explicit file_loader(const std::string &_filename,
                    const MPI_Comm& _comm = MPI_COMM_WORLD ) throw (io_exception)
            : filename(_filename), range(), fullRange(), file_handle(-1),
              data(nullptr), aligned_data(nullptr), loaded(false),
              nprocs(1), rank(0), comm(_comm)
#else
        explicit file_loader(const std::string &_filename) throw (io_exception)
            : filename(_filename), range(), fullRange(), file_handle(-1),
              data(nullptr), aligned_data(nullptr), loaded(false),
              nprocs(1), rank(0)
#endif
        {
#if defined(USE_MPI)
          // get the processor rank and nprocessors.
          MPI_Comm_rank(comm, &rank);
          MPI_Comm_size(comm, &nprocs);
#endif
          //DEBUG("CONSTRUCT");

          assert(filename.length() > 0);
          page_size = sysconf(_SC_PAGE_SIZE);

          /// get the file size.
          SizeType file_size = 0;
          if (rank == 0)
          {
            struct stat filestat;
            stat(filename.c_str(), &filestat);
            file_size = static_cast<SizeT>(filestat.st_size);
//            INFO("block size is " << filestat.st_blksize);
//            INFO("sysconf block size is " << page_size);
          }

#if defined(USE_MPI)
          if (nprocs > 1) {
            /// broadcast filesize to all
            MPI_Bcast(&file_size, 1, MPI_UNSIGNED_LONG_LONG, 0, comm);
          }
#endif
          if (rank == nprocs - 1)
          {
//            INFO("file size for " << filename << " is " << file_size);
          }

          // compute the full range
          fullRange.block_start = 0;
          fullRange.start = 0;
          fullRange.end = file_size / sizeof(T);   // range is in units of T
          fullRange.overlap = 0;

          // partition the full range
          range = fullRange.block_partition(nprocs, rank);
          range = range.align_to_page(page_size);

          /// open the file and get a handle.
          file_handle = open(filename.c_str(), O_RDONLY);
          if (file_handle == -1)
          {
            std::stringstream ss;
            int myerr = errno;
            ss << "ERROR in file open: ["  << filename << "] error " << myerr << ": " << strerror(myerr);
            throw io_exception(ss.str());
          }

          //DEBUG("file_loader initialized for " << filename);

        }

        /**
         * closes the file.
         */
        virtual ~file_loader()
        {
          //DEBUG("DESTROY");

          unload();

          //printf("unloading complete.\n");
          if (file_handle != -1) {
            close(file_handle);
            file_handle = -1;
          }
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
        inline SizeType size() {
          return range.end - range.start;
        }

        /**
         * return the full range for this file mmap.  (in units of data type T)
         * @return
         */
        const RangeType& getFullRange() const {
          //DEBUG("full range: " << fullRange);
          return fullRange;
        }

        /**
         * return the range for this file mmap.  (in units of data type T)
         * @return
         */
        const RangeType& getRange() const {
          return range;
        }

        /**
         * return the range for this file mmap.  (in units of data type T)
         * @return
         */
        void setRange(const RangeType& r) {
          //DEBUG("set range to: " << r);
          range = r;
          range.align_to_page(page_size);
        }

        /**
         * search at the start and end of the block partition for some matching condition,
         * and set as new partition positions.
         */
        template <typename PartitionHelper>
        void adjustRange(PartitionHelper &helper) throw (io_exception) {

          static_assert(std::is_same<typename PartitionHelper::SizeType, SizeType>::value, "PartitionHelper for adjustRange() has a different SizeType than FileLoader\n");

          // range where to search for the new end
          RangeType next = range >> (range.end - range.start);
          // limit to the full range length
          next &= fullRange;

          // concatenate current and next ranges
          RangeType searchRange = range | next;
          searchRange.align_to_page(page_size);

          // for output
          RangeType output(range);

          // map the content
          T* searchData = this->map(searchRange) + searchRange.start - searchRange.block_start;

          if (range.start > 0)
          {
            // not first block

            // search for new start using finder
            output.start = helper(searchData, range.start, range.end);
          }

          // get the new ending offset
          if (range.end < fullRange.end )
          {
            // not the last block
            output.end = helper(searchData + (range.end - range.start), next.start, next.end);
          }

          // clean up and unmap
          this->unmap(searchData, searchRange);

          // readjuste
          output.align_to_page(page_size);
          range = output;
        }


        /**
         *  performs memmap, and optionally preload the data into memory.
         *
         * @param _memUseFraction value between 0.0 and 0.5, representing percentage of FREE memory to be used.  0 means no preloading.
         */
        void load(const float _memUseFraction = 0.0f) throw (io_exception)
        {
          // clean up any previous runs.
          //DEBUG("Loading");
          unload();

          //DEBUG("loading");

          range.align_to_page(page_size);

          /// do the mem map
          aligned_data = map(range);

          /// check if can load the region into memory.
          float memUseFraction = (_memUseFraction > 0.5f) ? 0.5f :
              (_memUseFraction < 0.0f) ? 0.0f : _memUseFraction;


          loaded = false;
          if (memUseFraction > 0.0f) {
            /// check if we can preload.
            struct sysinfo memInfo;
            sysinfo (&memInfo);
            long long limit = static_cast<long long>(static_cast<float>(memInfo.freeram * memInfo.mem_unit) * memUseFraction);

            if (limit < (this->size() * sizeof(T)) ) {
              loaded = true;
            } else {
              //WARNING("Insufficient memory requested during file loading.  Not preloading file.");
            }
          }

          //DEBUG("mapped");

          if (loaded)
          {
            // allocate space
            data = new T[this->size() + 1];
            //copy data over
            memcpy(data, aligned_data + (range.start - range.block_start),
                   sizeof(T) * this->size());

            data[this->size()] = 0;

            // close the input
            unmap(aligned_data, range);

            aligned_data = data;
          }
          else
          {

            data = aligned_data + (range.start - range.block_start);

          }
          //DEBUG("loaded");

        }

        void unload()
        {
          //DEBUG("Unloading");

          if (loaded)
          {
            if (data != nullptr)
              delete [] data;
            loaded = false;
          }
          else
          {
            if (aligned_data != nullptr && aligned_data != MAP_FAILED)
              unmap(aligned_data, range);
          }
          aligned_data = nullptr;
          data = nullptr;

          //DEBUG("unloaded");

        }



      protected:
        std::string filename;
        RangeType range;  // offset in file from where to read
        RangeType fullRange;  // offset in file from where to read

        SizeType page_size;
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
         * map a region of the file to memory.
         * @param r
         * @return
         */
        T* map(RangeType r) throw (io_exception) {

          r.align_to_page(page_size);

          // NOT using MAP_POPULATE.  it slows things done when testing on single node.
          T* result = (T*)mmap(nullptr, (r.end - r.block_start ) * sizeof(T),
                                     PROT_READ,
                                     MAP_PRIVATE, file_handle,
                                     r.block_start * sizeof(T));

          if (result == MAP_FAILED)
          {

            if (file_handle != -1)
            {
              close(file_handle);
              file_handle = -1;
            }
            //DEBUG("Range: " << r);
            //DEBUG("file handle: " << file_handle);


            std::stringstream ss;
            int myerr = errno;
            ss << "ERROR in mmap: " << myerr << ": " << strerror(myerr);
            throw io_exception(ss.str());
          }

          //DEBUG("mapped");

          return result;
        }

        void unmap(T* d, RangeType r) {

          munmap(d, (r.end - r.block_start) * sizeof(T));
          //DEBUG("unmapped");

        }

      private:
        // disallow default constructor.
        file_loader() {};

    };

  } /* namespace functional */
} /* namespace bliss */
#endif /* FILE_LOADER_HPP_ */
