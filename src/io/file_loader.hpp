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

#if defined(USE_OPENMP)
#include "omp.h"
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
#include "io/data_block.hpp"
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
    template <typename T, bool Buffering = true, bool Preloading = false, typename SizeT = size_t>
    class file_loader
    {

        /////// type defintions.

      public:
        typedef SizeT                             SizeType;
        typedef bliss::iterator::range<SizeType>  RangeType;
        typedef T*                                        InputIteratorType;

      protected:
        // internal DataBlock type.  uses T* as iterator type.
        typedef typename std::conditional<Preloading,
                                          bliss::io::BufferedDataBlock<InputIteratorType, RangeType>,
                                          bliss::io::UnbufferedDataBlock<InputIteratorType, RangeType> >::type     DataType;

      public:
        typedef typename std::conditional<Buffering,
                                          bliss::io::BufferedDataBlock<typename DataType::iterator, RangeType>,
                                          bliss::io::UnbufferedDataBlock<typename DataType::iterator, RangeType> >::type     DataBlockType;

        typedef typename DataType::iterator               IteratorType;               // this is for use with the Preloaded datablock
        typedef typename DataBlockType::iterator          BlockIteratorType;     // this is for use with the buffered data chunk

        ////// member variables
      protected:
        SizeType page_size;
        int file_handle;      // file handle
        std::string filename;

        RangeType fileRange;  // offset in file from where to read
        InputIteratorType aligned_data;      // mem-mapped data, page aligned.  strictly internal


        DataType srcData;
        RangeType range;      // offset in file from where to read


        bool loaded;
        bool preloaded;

        SizeType chunkPos;    // for chunked iteration.  size since "data".

        int nprocs;           // for partitioning.
        int rank;             // for partitioning.
#if defined(USE_MPI)
        MPI_Comm comm;
#endif

        DataBlockType *dataBlocks;

        SizeType seqSize;

        /////////////// Constructor and Destructor
      public:
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
            : file_handle(-1), filename(_filename), fileRange(), aligned_data(nullptr),
              range(), loaded(false), preloaded(false),
              chunkPos(0),
              nprocs(1), rank(0), comm(_comm), dataBlocks(nullptr), seqSize(0)
#else
        explicit file_loader(const std::string &_filename,
                             const int _nParts = 1, const int _rank = 0) throw (io_exception)
            : file_handle(-1), filename(_filename), fileRange(), aligned_data(nullptr),
              range(), loaded(false), preloaded(false),
              chunkPos(0),
              nprocs(_nParts), rank(_rank), dataBlocks(nullptr), seqSize(0)
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
          fileRange.block_start = 0;
          fileRange.start = 0;
          fileRange.end = file_size / sizeof(T);   // range is in units of T
          fileRange.overlap = 0;

          // partition the full range
          range = fileRange.block_partition(nprocs, rank);
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


          /// allocate per thread datablocks
          int nthreads = 1;
#if defined(USE_OPENMP)
          nthreads = omp_get_max_threads();
#endif
          dataBlocks = new DataBlockType[nthreads];
        }

        /**
         * closes the file.
         */
        virtual ~file_loader()
        {
          //DEBUG("DESTROY");
          delete [] dataBlocks;

          unload();

          //printf("unloading complete.\n");
          if (file_handle != -1) {
            close(file_handle);
            file_handle = -1;
          }
        }


      private:
        // disallow default constructor.
        file_loader() {};


        ////// PUBLIC METHODS


      public:

        /**
         *  performs memmap, and optionally preload the data into memory.
         *
         * @param _memUseFraction value between 0.0 and 0.5, representing percentage of FREE memory to be used.  0 means no preloading.
         *                        if template specified preloading=false, then this param is ignored.
         */
        void load() throw (io_exception)
        {
          // clean up any previous runs.
          //DEBUG("Loading");
          unload();

          //DEBUG("loading");

          range.align_to_page(page_size);

          /// do the mem map
          aligned_data = map(range);


//////////////// check to see if there is enough memory for preloading.
///  since Preloading is a template param, can't choose a srcData type if Preloading is true but there is not enough memory.
/// so we have to assert it.
//          preloaded = false;
//          if (Preloading) {
//            /// check if can load the region into memory.
//            float memUseFraction = (_memUseFraction > 0.5f) ? 0.5f :
//                (_memUseFraction < 0.0f) ? 0.0f : _memUseFraction;
//
//            if (memUseFraction > 0.0f) {
//              /// check if we can preload.
//              struct sysinfo memInfo;
//              sysinfo (&memInfo);
//              long long limit = static_cast<long long>(static_cast<float>(memInfo.freeram * memInfo.mem_unit) * memUseFraction);
//
//              if (limit < (this->size() * sizeof(T)) ) {
//                preloaded = true;
//              } else {
//                //WARNING("Insufficient memory requested during file loading.  Not preloading file.");
//              }
//            }
//          }
          if (Preloading) {
            /// check if can load the region into memory.

            /// check if we can preload.  use at most 1/20 of the available memory because of kmer expansion.
            struct sysinfo memInfo;
            sysinfo (&memInfo);

            if ((this->size() * sizeof(T)) > (memInfo.freeram * memInfo.mem_unit / 20)) { // use too much mem.  throw exception.
              std::stringstream ss;
              ss << "ERROR in file preloading: ["  << filename << "] mem required " << (this->size() * sizeof(T)) << " bytes; (1/20th) available: " << (memInfo.freeram * memInfo.mem_unit / 20);
              throw io_exception(ss.str());
            }
          }



          //DEBUG("mapped");
          srcData.assign(aligned_data + (range.start - range.block_start), aligned_data + (range.end - range.block_start), range);

          if (Preloading)
          {
//            // allocate space
//            data = new T[this->size() + 1];
//            //copy data over
//            memcpy(data, aligned_data + (range.start - range.block_start),
//                   sizeof(T) * this->size());
//
//            data[this->size()] = 0;
//
//            // close the input
//            unmap(aligned_data, range);
//
//            aligned_data = data;
            unmap(aligned_data, range);
//          }
//          else
//          {
//            data = aligned_data + (range.start - range.block_start);
          }
          //DEBUG("loaded");
          loaded = true;
        }

        void unload()
        {
          //DEBUG("Unloading");

          if (!Preloading)
          {
            if (aligned_data != nullptr && aligned_data != MAP_FAILED)
              unmap(aligned_data, range);
          }
//          else
//          {
//            //            if (data != nullptr)
//            //              delete [] data;
//                        preloaded = false;
//          }
          srcData.clear();
          aligned_data = nullptr;
          //DEBUG("unloaded");
          loaded = false;
        }

        DataType& getData() {
          assert(loaded);

          return srcData;
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
         * return the full range for this file.  (in units of data type T)
         * @return
         */
        const RangeType& getFileRange() const {
          //DEBUG("full range: " << fullRange);
          return fileRange;
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
        void adjustRange(const PartitionHelper &partitioner) throw (io_exception) {

          static_assert(std::is_same<typename PartitionHelper::SizeType, SizeType>::value, "PartitionHelper for adjustRange() has a different SizeType than FileLoader\n");
          static_assert(std::is_same<typename PartitionHelper::ValueType, T>::value, "PartitionHelper for adjustRange() has a different ValueType than FileLoader\n");
          static_assert(std::is_same<typename PartitionHelper::IteratorType, InputIteratorType>::value, "PartitionHelper for adjustRange() has a different IteratorType than FileLoader\n");


          // range where to search for the new end
          RangeType next = range >> (range.end - range.start);
          // limit to the full range length

          // concatenate current and next ranges
          RangeType searchRange = range | next;
          searchRange &= fileRange;
          searchRange.align_to_page(page_size);

          // map the content
          InputIteratorType searchData = this->map(searchRange) + searchRange.start - searchRange.block_start;

          // NOTE: assume that the range is large enough that we would not run into trouble with partition search not finding a starting position.

          // for output
          RangeType output(range);

          // search for new start using finder
          output.start = partitioner(searchData, fileRange, range);

          // get the new ending offset
          output.end = partitioner(searchData + (range.end - range.start), fileRange, next);

          // clean up and unmap
          this->unmap(searchData, searchRange);

          // readjust
          this->setRange(output);
        }



//
//
//        template <typename PartitionHelper>
//        RangeType adjustChunkRange(PartitionHelper &partitioner, const RangeType &chunkRange) {
//          // make sure the file is open
//          assert(loaded);
//
//          // make sure the search range is within bounds
//          RangeType next = (chunkRange >> (chunkRange.end - chunkRange.start)) & range;
//
//          // adjust the end only - assume start is okay.
//          auto e = partitioner(srcData.begin() + (next.start - range.start), range, next);
//
//          RangeType result(chunkRange);
//          result.end = e;
//          return result;
//        }

        void resetChunks() {
          chunkPos = range.start;
        }

//        /**
//         * search at the start and end of the block partition for some matching condition,
//         * and set as new partition positions.
//         *
//         * @param[in]  partitioner       to find the partition boundaries
//         * @param[out] start        output start pointer. at least "begin".  if copying, then caller needs to manage "start"
//         * @param[out] end          output end pointer.  at most "end"
//         * @param[in]  chunkSize    suggested partition size.  default to 0, which is translated to system page size.
//         * @param[in]  copying      if copying, then start and end point to a memory block that is a copy of the underlying file mmap.
//         * @return                  absolute range of for the chunk returned.
//         */
//        template <typename PartitionHelper>
//        RangeType getNextChunk(PartitionHelper &partitioner, T* &startPtr, T* &endPtr, const SizeType &chunkSize = 0, const bool copying = true) {
//
//          static_assert(std::is_same<typename PartitionHelper::SizeType, SizeType>::value, "PartitionHelper for getNextChunk() has a different SizeType than FileLoader\n");
//          assert(chunkSize >= 0);
//          assert(loaded);
//
////          assert(loaded);
//
//          // traversed all.  so done.
//          SizeType len = range.end - range.start;
//          SizeType cs = (chunkSize == 0 ? page_size : chunkSize);
//          SizeType s = range.end, e = range.end;
//          SizeType readLen = 0;
//
//          /// part that needs to be sequential
//
//          //DEBUG("range " << range << " chunk " << cs << " chunkPos " << chunkPos);
//
//#pragma omp critical (computeChunkBoundary)
//          {
//            if (chunkPos >= range.end) {
//              endPtr = startPtr = data + len;
//              readLen = 0;
//            } else {
//
//              s = chunkPos;
//              RangeType next(chunkPos + cs, chunkPos + 2 * cs);
//              next &= range;
//              /// search end part only.  update chunkPos for next search.
//              try {
//                // search for end.
//                e = partitioner(data + (next.start - range.start), range.start, next.start, next.end);
//              } catch (io_exception& ex) {
//                // did not find the end, so set e to next.end.
//                // TODO: need to handle this scenario better - should keep search until end.
//                WARNING(ex.what());
//                e = next.end;
//              }
//
//
//              chunkPos = e;
//#pragma omp flush(chunkPos)
//              readLen = e - s;
//              //printf("s = %lu, e = %lu, chunkPos = %lu, readLen 1 = %lu\n", s, e, chunkPos, readLen);
//              //std::cout << "range: " << range << " next: " << next << std::endl;
//            }
//          }
//          /// end part that needs to be serial
//
//          //DEBUG(" new chunkPos " << chunkPos << " start: " << s << " end: " << e);
//
//          // now can set in parallel
//          // set the start.  guaranteed to be at the start of a record according to the partitionhelper.
//
//          if (readLen > 0) {
//             if (copying) {
//              //unique_ptr is not appropriate, as ownership is either by this object (no copy) or by calling thread. (copying)
//              startPtr = new T[readLen + 1];
//              endPtr = startPtr + readLen;
//              *endPtr = 0;
//              memcpy(startPtr, data + (s - range.start), readLen * sizeof(T));
//            } else {
//              startPtr = data + (s - range.start);
//              endPtr = data + (e - range.start);
//            }
//         }
////         DEBUG("read " << readLen << " elements, start at " << s << " [" << *startPtr << "] end at "<< e << " [" << *endPtr << "]");
//
//          if (startPtr == nullptr || endPtr == nullptr) {
//            std::cerr << "ERROR: file loader get chunk returning null ptrs. readlen = " << readLen << " start and end are " << s << "-" << e << " range is " << range << std::endl;
//          }
//        return RangeType(s, e);
//      }


        /**
         * search at the start and end of the block partition for some matching condition,
         * and set as new partition positions.
         *
         * @param[in]  partitioner       to find the partition boundaries
         * @param[out] start        output start pointer. at least "begin".  if copying, then caller needs to manage "start"
         * @param[out] end          output end pointer.  at most "end"
         * @param[in]  chunkSize    suggested partition size.  default to 0, which is translated to system page size.
         * @param[in]  copying      if copying, then start and end point to a memory block that is a copy of the underlying file mmap.
         * @return                  actual chunk size created
         */
        template <typename PartitionHelper>
        DataBlockType& getNextChunkAtomic(const PartitionHelper &partitioner, const SizeType &chunkSize = 0) {

          getSeqSize(partitioner, 3);


          static_assert(std::is_same<typename PartitionHelper::SizeType, SizeType>::value, "PartitionHelper for getNextChunk() has a different SizeType than FileLoader\n");
          assert(loaded);
//          printf("chunkSize = %lu\n", chunkSize);
          assert(chunkSize >= (seqSize * 2));


          int tid = 0;
#if defined(USE_OPENMP)
          tid = omp_get_thread_num();
#endif


          SizeType cs = chunkSize;  // we need at least seqSize * 2
          SizeType s = range.end,
                   e = range.end;
          SizeType readLen = 0;


          //DEBUG("range " << range << " chunk " << cs << " chunkPos " << chunkPos);

          /// part that needs to be sequential
#pragma omp atomic capture
          { s = chunkPos; chunkPos += cs; }  // get the current value and then update

          if (s >= range.end) {
            printf("should not be here often\n");
              s = range.end;
              readLen = 0;
          } else {
            // since we are just partitioning, need to search for start and for end.
            RangeType curr(s, s + cs);
            curr &= srcData.getRange();
            /// search start part.

             //printf("curr: %lu %lu\n", s, s+cs);
            try {
              // search for end.
              s = partitioner(srcData.begin() + (curr.start - range.start), range, curr);
            } catch (io_exception& ex) {
              // did not find the end, so set e to next.end.
              // TODO: need to handle this scenario better - should keep search until end.
              WARNING(ex.what());

              printf("got an exception:  %s \n", ex.what());
              s = curr.end;
            }

            RangeType next = curr >> cs;
            next &= srcData.getRange();
            /// search end part only.  update localPos for next search.
            //printf("next: %lu %lu\n", next.start, next.end);
            try {
              // search for end.
              e = partitioner(srcData.begin() + (next.start - range.start), range, next);
            } catch (io_exception& ex) {
              // did not find the end, so set e to next.end.
              // TODO: need to handle this scenario better - should keep search until end.
              printf("got an exception\n");
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

          dataBlocks[tid].assign(srcData.begin() + (s - range.start), srcData.begin() + (e - range.start), RangeType(s,e) & range);
//         DEBUG("read " << readLen << " elements, start at " << s << " [" << *startPtr << "] end at "<< e << " [" << *endPtr << "]");

//          if (startPtr == nullptr || endPtr == nullptr) {
//            std::cerr << "ERROR: file loader get chunk returning null ptrs. readlen = " << readLen << " start and end are " << s << "-" << e << " range is " << range << std::endl;
//          }
        return dataBlocks[tid];
      }

      template<typename PartitionHelper>
      SizeType getSeqSize(const PartitionHelper &partitioner, int iterations) {

        //// TODO: if a different partitioner is used, seqSize may be incorrect.
        ////   seqSize is a property of the partitioner applied to the data.

        assert(loaded);

        if (seqSize == 0) {
          SizeType s, e;
          RangeType r(range);

          s = partitioner(srcData.begin(), range, r);

          for (int i = 0; i < iterations; ++i) {
            r.start = s + 1;   // advance by 1, in order to search for next entry.
            r &= range;
            e = partitioner(srcData.begin() + (r.start - range.start), range, r);
            seqSize += (e - s);
            s = e;
          }
          seqSize /= iterations;

          // DEBUG("Sequence size is " << seqSize);
        }
        return seqSize;
      }



      protected:



        /**
         * map a region of the file to memory.
         * @param r
         * @return
         */
        InputIteratorType map(RangeType r) throw (io_exception) {

          r.align_to_page(page_size);

          // NOT using MAP_POPULATE.  it slows things done when testing on single node.
          InputIteratorType result = (InputIteratorType)mmap64(nullptr, (r.end - r.block_start ) * sizeof(T),
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

        void unmap(InputIteratorType d, RangeType r) {

          munmap(d, (r.end - r.block_start) * sizeof(T));
          //DEBUG("unmapped");

        }

    };

  } /* namespace functional */
} /* namespace bliss */
#endif /* FILE_LOADER_HPP_ */
