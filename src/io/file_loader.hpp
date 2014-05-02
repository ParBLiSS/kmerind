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
#include <cstring>      // memcpy, strerror
#include <exception>    // ioexception
#include <sstream>      // stringstream
#include <memory>

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
         *                        to get just the calling node, use MPI_COMM_SELF (gives the right nprocs).
         */
#if defined(USE_MPI)
        explicit file_loader(const std::string &_filename,
                    const MPI_Comm& _comm = MPI_COMM_SELF) throw (io_exception)
            : filename(_filename), range(), fullRange(), file_handle(-1),
              data(nullptr), aligned_data(nullptr), loaded(false), preloaded(false),
              chunkPos(0),
              nprocs(1), rank(0), comm(_comm)
#else
        explicit file_loader(const std::string &_filename,
                             const int _nParts = 1, const int _rank = 0) throw (io_exception)
            : filename(_filename), range(), fullRange(), file_handle(-1),
              data(nullptr), aligned_data(nullptr), loaded(false),preloaded(false),
              chunkPos(0),
              nprocs(_nParts), rank(_rank)
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
//            std::cerr << "file size is " << file_size;
//            std::cerr << " block size is " << filestat.st_blksize;
//            std::cerr << " sysconf block size is " << page_size << std::endl;
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

          chunkPos = range.start;

          /// open the file and get a handle.
          file_handle = open64(filename.c_str(), O_RDONLY);
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
          assert(loaded);

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
          assert(loaded);

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
          chunkPos = range.start;
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


          // map the content
          T* searchData = this->map(searchRange) + searchRange.start - searchRange.block_start;

          // for output
          RangeType output(range);
          if (range.start > fullRange.start)
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
          chunkPos = range.start;
        }

        template <typename PartitionHelper>
        RangeType adjustChunkRange(PartitionHelper &helper, const RangeType &chunkRange) {
          // make sure the file is open
          assert(data != nullptr);

          // make sure the search range is within bounds
          RangeType next = (chunkRange >> (chunkRange.end - chunkRange.start)) & range;

          // adjust the end only - assume start is okay.
          auto e = helper(data + (next.start - range.start), next.start, next.end);

          RangeType result(chunkRange);
          result.end = e;
          return result;
        }

        void resetChunks() {
          chunkPos = range.start;
        }

        /**
         * search at the start and end of the block partition for some matching condition,
         * and set as new partition positions.
         *
         * @param[in]  helper       to find the partition boundaries
         * @param[out] start        output start pointer. at least "begin".  if copying, then caller needs to manage "start"
         * @param[out] end          output end pointer.  at most "end"
         * @param[in]  chunkSize    suggested partition size.  default to 0, which is translated to system page size.
         * @param[in]  copying      if copying, then start and end point to a memory block that is a copy of the underlying file mmap.
         * @return                  actual chunk size created
         */
        template <typename PartitionHelper>
        SizeType getNextChunk(PartitionHelper &helper, T* &startPtr, T* &endPtr, const SizeType &chunkSize = 0, const bool copying = true) {

          static_assert(std::is_same<typename PartitionHelper::SizeType, SizeType>::value, "PartitionHelper for getNextChunk() has a different SizeType than FileLoader\n");

//          assert(loaded);

          // traversed all.  so done.
          SizeType len = range.end - range.start;
          SizeType cs = (chunkSize == 0 ? page_size : chunkSize);
          SizeType s = range.start, e = range.end;
          SizeType readLen = 0;

          /// part that needs to be sequential

          //DEBUG("range " << range << " chunk " << cs << " chunkPos " << chunkPos);

#pragma omp critical (computeChunkBoundary)
          {
            if (chunkPos >= range.end) {
              endPtr = startPtr = data + len;
              readLen = 0;
            } else {

              s = chunkPos;
              RangeType next(chunkPos + cs, chunkPos + 2 * cs);
              next &= range;
              /// search end part only.  update chunkPos for next search.
              try {
                // search for end.
                e = helper(data + (next.start - range.start), next.start, next.end);
              } catch (io_exception& ex) {
                // did not find the end, so set e to next.end.
                // TODO: need to handle this scenario better - should keep search until end.
                WARNING(ex.what());
                e = next.end;
              }


              chunkPos = e;
#pragma omp flush(chunkPos)
              readLen = e - s;
              //printf("s = %lu, e = %lu, chunkPos = %lu, readLen 1 = %lu\n", s, e, chunkPos, readLen);
              //std::cout << "range: " << range << " next: " << next << std::endl;
            }
          }
          /// end part that needs to be serial

          //DEBUG(" new chunkPos " << chunkPos << " start: " << s << " end: " << e);

          // now can set in parallel
          // set the start.  guaranteed to be at the start of a record according to the partitionhelper.

          if (readLen > 0) {
             if (copying) {
              //unique_ptr is not appropriate, as ownership is either by this object (no copy) or by calling thread. (copying)
              startPtr = new T[readLen + 1];
              endPtr = startPtr + readLen;
              *endPtr = 0;
              memcpy(startPtr, data + (s - range.start), readLen * sizeof(T));
            } else {
              startPtr = data + (s - range.start);
              endPtr = data + (e - range.start);
            }
         }
//         DEBUG("read " << readLen << " elements, start at " << s << " [" << *startPtr << "] end at "<< e << " [" << *endPtr << "]");

          if (startPtr == nullptr || endPtr == nullptr) {
            std::cerr << "ERROR: file loader get chunk returning null ptrs. readlen = " << readLen << " start and end are " << s << "-" << e << " range is " << range << std::endl;
          }
        return readLen;
      }


        /**
         * search at the start and end of the block partition for some matching condition,
         * and set as new partition positions.
         *
         * @param[in]  helper       to find the partition boundaries
         * @param[out] start        output start pointer. at least "begin".  if copying, then caller needs to manage "start"
         * @param[out] end          output end pointer.  at most "end"
         * @param[in]  chunkSize    suggested partition size.  default to 0, which is translated to system page size.
         * @param[in]  copying      if copying, then start and end point to a memory block that is a copy of the underlying file mmap.
         * @return                  actual chunk size created
         */
        template <typename PartitionHelper>
        SizeType getNextChunkAtomic(PartitionHelper &helper, T* &startPtr, T* &endPtr, const SizeType &chunkSize = 0, const bool copying = true) {

          static_assert(std::is_same<typename PartitionHelper::SizeType, SizeType>::value, "PartitionHelper for getNextChunk() has a different SizeType than FileLoader\n");

//          assert(loaded);

          // traversed all.  so done.
          SizeType len = range.end - range.start;
          SizeType cs = (chunkSize == 0 ? page_size : chunkSize);
          SizeType s = range.start, e = range.end;
          SizeType readLen = 0;

          /// part that needs to be sequential

          //DEBUG("range " << range << " chunk " << cs << " chunkPos " << chunkPos);


#pragma omp atomic capture
          { s = chunkPos; chunkPos += cs; }  // get the current value and then update

          if (s > range.end) {
              s = range.end;
              readLen = 0;
          } else {
            // since we are just partitioning, need to search for start and for end.
            RangeType curr(s, s + cs);
            curr &= range;
            /// search start part.
            try {
              // search for end.
              s = helper(data + (curr.start - range.start), curr.start, curr.end);
            } catch (io_exception& ex) {
              // did not find the end, so set e to next.end.
              // TODO: need to handle this scenario better - should keep search until end.
              WARNING(ex.what());
              s = curr.start;
            }

            RangeType next = curr >> cs;
            next &= range;
            /// search end part only.  update localPos for next search.
            try {
              // search for end.
              e = helper(data + (next.start - range.start), next.start, next.end);
            } catch (io_exception& ex) {
              // did not find the end, so set e to next.end.
              // TODO: need to handle this scenario better - should keep search until end.
              WARNING(ex.what());
              e = next.end;
            }

            readLen = e - s;
            //printf("s = %lu, e = %lu, chunkPos = %lu, readLen 1 = %lu\n", s, e, chunkPos, readLen);
            //std::cout << "range: " << range << " next: " << next << std::endl;
          }

          /// end part that needs to be serial

          //DEBUG(" new chunkPos " << chunkPos << " start: " << s << " end: " << e);

          // now can set in parallel
          // set the start.  guaranteed to be at the start of a record according to the partitionhelper.

          if (readLen > 0) {
             if (copying) {
              //unique_ptr is not appropriate, as ownership is either by this object (no copy) or by calling thread. (copying)
              startPtr = new T[readLen + 1];
              endPtr = startPtr + readLen;
              *endPtr = 0;
              memcpy(startPtr, data + (s - range.start), readLen * sizeof(T));
            } else {
              startPtr = data + (s - range.start);
              endPtr = data + (e - range.start);
            }
         }
//         DEBUG("read " << readLen << " elements, start at " << s << " [" << *startPtr << "] end at "<< e << " [" << *endPtr << "]");

//          if (startPtr == nullptr || endPtr == nullptr) {
//            std::cerr << "ERROR: file loader get chunk returning null ptrs. readlen = " << readLen << " start and end are " << s << "-" << e << " range is " << range << std::endl;
//          }
        return readLen;
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


          preloaded = false;
          if (memUseFraction > 0.0f) {
            /// check if we can preload.
            struct sysinfo memInfo;
            sysinfo (&memInfo);
            long long limit = static_cast<long long>(static_cast<float>(memInfo.freeram * memInfo.mem_unit) * memUseFraction);

            if (limit < (this->size() * sizeof(T)) ) {
              preloaded = true;
            } else {
              //WARNING("Insufficient memory requested during file loading.  Not preloading file.");
            }
          }

          //DEBUG("mapped");

          if (preloaded)
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
          loaded = true;

        }

        void unload()
        {
          //DEBUG("Unloading");

          if (preloaded)
          {
            if (data != nullptr)
              delete [] data;
            preloaded = false;
          }
          else
          {
            if (aligned_data != nullptr && aligned_data != MAP_FAILED)
              unmap(aligned_data, range);
          }
          aligned_data = nullptr;
          data = nullptr;
          //DEBUG("unloaded");
          loaded = false;
        }



      protected:
        std::string filename;
        RangeType range;      // offset in file from where to read
        RangeType fullRange;  // offset in file from where to read

        SizeType page_size;
        int file_handle;      // file handle

        T* data;              // actual start of data
        T* aligned_data;      // mem-mapped data, page aligned.  strictly internal

        bool loaded;
        bool preloaded;

        SizeType chunkPos;    // for chunked iteration.  size since "data".

        int nprocs;           // for partitioning.
        int rank;             // for partitioning.

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
          T* result = (T*)mmap64(nullptr, (r.end - r.block_start ) * sizeof(T),
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
