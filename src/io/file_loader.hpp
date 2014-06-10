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

#include "partition/range.hpp"
#include "io/data_block.hpp"
#include "io/io_exception.hpp"
#include "utils/logging.h"


namespace bliss
{
  namespace io
  {


    /**
     *  real data:  mmap is better for large files and limited memory.
     *              preloading is better for smaller files and/or large amount of memory.
     *  stream processing means data does not need to be buffered in memory - more efficient.
     *
     *
     *  file access:  better to work with a few pages at a time, or to work with large block?
     *
     *
     *  To reuse chunk and range, we'd want to limit to only a fixed number of entries.
     *  however, this means we it's most convenient to say that the max number of entires is the number of threads.
     *  else we don't have a good way of matching the tid to threads.
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
    template <typename T, bool Buffering = true, bool Preloading = false,
          typename ChunkPartitioner = bliss::partition::DemandDrivenPartitioner<bliss::partition::range<size_t> >,
          typename Derived = void >
    class FileLoader
    {

        /////// type defintions.

      public:
        typedef size_t                             SizeType;
        typedef bliss::partition::range<SizeType>  RangeType;
        typedef T*                                 InputIteratorType;

        // internal DataBlock type.  uses T* as iterator type.
        typedef typename std::conditional<Preloading,
                                          bliss::io::BufferedDataBlock<InputIteratorType, RangeType>,
                                          bliss::io::UnbufferedDataBlock<InputIteratorType, RangeType> >::type              DataType;

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
        T* mmap_data;      // mem-mapped data, page aligned.  strictly internal
        RangeType mmap_range;      // offset in file from where to read

        DataType srcData;

        bool loaded;
        bool preloaded;

//        SizeType chunkPos;    // for chunked iteration.  size since "data".

        int nthreads;
        size_t chunkSize;
        int nprocs;           // for partitioning.
        int rank;             // for partitioning.
#if defined(USE_MPI)
        MPI_Comm comm;
#endif

        DataBlockType *dataBlocks;

        bliss::partition::BlockPartitioner<RangeType > partitioner;   // block partitioning is the default.
        ChunkPartitioner chunkPartitioner;                          // construct when mmap_range change, on load.


        SizeType recordSize;

        /////////////// Constructor and Destructor
      public:

        /// defining move constructor will disallow automatic copy constructor.
        //  this is to prevent copying the buffers, to keep a single handle on the input file.
        FileLoader(FileLoader<T, Buffering, Preloading, ChunkPartitioner, Derived>&& other) {
          page_size   = other.page_size;                         other.page_size = 1;
          file_handle = other.file_handle;                       other.file_handle = -1;
          filename.swap(other.filename);
          fileRange   = other.fileRange;                         other.fileRange = RangeType();
          mmap_data   = other.mmap_data;                         other.mmap_data = nullptr;
          mmap_range  = other.mmap_range;                        other.mmap_range = RangeType();
          srcData     = std::move(other.srcData);
          loaded      = other.loaded;                            other.loaded = false;
          preloaded   = other.preloaded;                         other.preloaded = false;
          nthreads    = other.nthreads;                          other.nthreads = 1;
          chunkSize   = other.chunkSize;                         other.chunkSize = 1;
          nprocs      = other.nprocs;                            other.nprocs = 1;
          rank        = other.rank;                              other.rank = 0;
#if defined(USE_MPI)
          comm        = other.comm;                              other.comm = MPI_COMM_NULL;
#endif
          dataBlocks  = other.dataBlocks;                        other.dataBlocks = nullptr;
          partitioner = other.partitioner;                       other.partitioner = bliss::partition::BlockPartitioner<RangeType>();
          chunkPartitioner = other.chunkPartitioner;             other.chunkPartitioner = ChunkPartitioner();
          recordSize  = other.recordSize;                        other.recordSize = 1;
        }

        /// defining move assignment will disallow automatic copy assignment operator
        FileLoader& operator=(FileLoader<T, Buffering, Preloading, ChunkPartitioner, Derived>&& other) {
          if (this != &other) {
            //DEBUG("DESTROY");
            if (dataBlocks != nullptr) delete [] dataBlocks;

            unload();

            //printf("unloading complete.\n");
            if (file_handle != -1) {
              close(file_handle);
              file_handle = -1;
            }

            page_size   = other.page_size;                         other.page_size = 1;
            file_handle = other.file_handle;                       other.file_handle = -1;
            filename.swap(other.filename);
            fileRange   = other.fileRange;                         other.fileRange = RangeType();
            mmap_data   = other.mmap_data;                         other.mmap_data = nullptr;
            mmap_range  = other.mmap_range;                        other.mmap_range = RangeType();
            srcData     = std::move(other.srcData);
            loaded      = other.loaded;                            other.loaded = false;
            preloaded   = other.preloaded;                         other.preloaded = false;
            nthreads    = other.nthreads;                          other.nthreads = 1;
            chunkSize   = other.chunkSize;                         other.chunkSize = 1;
            nprocs      = other.nprocs;                            other.nprocs = 1;
            rank        = other.rank;                              other.rank = 0;
  #if defined(USE_MPI)
            comm        = other.comm;                              other.comm = MPI_COMM_NULL;
  #endif
            dataBlocks  = other.dataBlocks;                        other.dataBlocks = nullptr;
            partitioner = other.partitioner;                       other.partitioner = bliss::partition::BlockPartitioner<RangeType>();
            chunkPartitioner = other.chunkPartitioner;             other.chunkPartitioner = ChunkPartitioner();
            recordSize  = other.recordSize;                        other.recordSize = 1;
          }
          return *this;
        }


        /// defining a constructor automatically disallows the default no arg constructor.


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
         * default behavior is to block partition the file based on number of processors.  to make it load only part, call getMMapRange, then modify it.
         *
         * TODO: test on remotely mounted file system.
         *
         * @param _filename       input file name
         * @param _comm           MPI Communicator (defined only if CMake is configured with MPI).
         *                        default is MPI_COMM_WORLD, which means all nodes participate.
         *                        to get just the calling node, use MPI_COMM_SELF (gives the right nprocs).
         */
#if defined(USE_MPI)
        FileLoader(const std::string &_filename, const MPI_Comm& _comm, const int _nThreads = 1, const size_t _chunkSize = 1) throw (bliss::io::IOException)
            : file_handle(-1), filename(_filename), fileRange(), mmap_data(nullptr),
              mmap_range(), loaded(false), preloaded(false),
              nthreads(_nThreads), chunkSize(_chunkSize), comm(_comm), dataBlocks(nullptr), recordSize(1)
        {
          assert(filename.length() > 0);
          assert(_nThreads > 0);
          assert(_chunkSize > 0);

          // get the processor rank and nprocessors.
          MPI_Comm_rank(comm, &rank);
          MPI_Comm_size(comm, &nprocs);

          //DEBUG("CONSTRUCT");

          /// get the file size.
          SizeType file_size = 0;
          if (rank == 0)
          {
            struct stat filestat;
            int ret = stat(filename.c_str(), &filestat);

            if (ret < 0) {

              std::stringstream ss;
              ss << "ERROR in file size detection: ["  << filename << "] error ";
              throw IOException(ss.str());
            }

            file_size = static_cast<SizeType>(filestat.st_size);
//            std::cerr << "file size is " << file_size;
//            std::cerr << " block size is " << filestat.st_blksize;
//            std::cerr << " sysconf block size is " << page_size << std::endl;
          }

          if (nprocs > 1) {
            /// broadcast filesize to all
            MPI_Bcast(&file_size, 1, MPI_UNSIGNED_LONG_LONG, 0, comm);
          }

//          if (rank == nprocs - 1)
//          {
//            INFO("file size for " << filename << " is " << file_size);
//          }

          init(file_size);
        }
#endif


        FileLoader(const std::string &_filename,
                            const int _nProcs = 1, const int _rank = 0,
                            const int _nThreads = 1, const size_t _chunkSize = 1 ) throw (bliss::io::IOException)
            : file_handle(-1), filename(_filename), fileRange(), mmap_data(nullptr),
              mmap_range(), loaded(false), preloaded(false),
              nthreads(_nThreads), chunkSize(_chunkSize), nprocs(_nProcs), rank(_rank), dataBlocks(nullptr), recordSize(1)
        {
          //DEBUG("CONSTRUCT");

          assert(filename.length() > 0);
          assert(_nThreads > 0);
          assert(_chunkSize > 0);
          assert(_rank >= 0);
          assert(_nProcs > _rank);


          /// get the file size.
          SizeType file_size = 0;
          struct stat filestat;
          int ret = stat(filename.c_str(), &filestat);

          if (ret < 0) {

            std::stringstream ss;
            ss << "ERROR in file size detection: ["  << filename << "] error ";
            throw IOException(ss.str());
          }

          file_size = static_cast<SizeType>(filestat.st_size);
//            std::cerr << "file size is " << file_size;
//            std::cerr << " block size is " << filestat.st_blksize;
//            std::cerr << " sysconf block size is " << page_size << std::endl;

//          if (rank == nprocs - 1)
//          {
//            INFO("file size for " << filename << " is " << file_size);
//          }

          init(file_size);
        }


        /**
         * closes the file.
         */
        virtual ~FileLoader()
        {
          //DEBUG("DESTROY");
          if (dataBlocks != nullptr) delete [] dataBlocks;

          unload();

          //printf("unloading complete.\n");
          if (file_handle != -1) {
            close(file_handle);
            file_handle = -1;
          }
        }



        ////// PUBLIC METHODS

        const size_t& getChunkSize() const {
          return chunkSize;
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
         * return the range for this file mmap.  (in units of data type T).  valid only after loading.
         * @return
         */
        const RangeType& getMMapRange() const {
          assert(loaded);
          return mmap_range;
        }


        void resetPartitionRange() {
          partitioner.reset();
        }

        /**
         * get the Partitioned Range.  this method allows calling Derived class' version, to further refine the partitioned range.
         * have a default.  if default is -1, then use rank.  not in inner loop so okay to check param value each call
         *
         * @param pid
         * @return
         */
        template<typename D = Derived>
        typename std::enable_if<std::is_void<D>::value, RangeType>::type getNextPartitionRange(const int pid = -1) {
          // if just FileLoader calls this.
          if (pid == -1)
            return partitioner.getNext(rank);
          else
            return partitioner.getNext(pid);
        }
        template<typename D = Derived>
        typename std::enable_if<!std::is_void<D>::value, RangeType>::type getNextPartitionRange(const int pid = -1) {
          // use the derived one.
          if (pid == -1)
            return static_cast<Derived*>(this)->getNextPartitionRangeImpl(rank);
          else
            return static_cast<Derived*>(this)->getNextPartitionRangeImpl(pid);
        }



        void resetChunkRange() {
          assert(loaded);
          chunkPartitioner.reset();
        }

        // partitioner has chunk size. and current position.  at most as many as number of threads (specified in constructor) calling this function.
        // not using a default tid because this function could be called a lot.
        /**
         *
         * @param tid
         * @return
         */
        template <typename D = Derived>
        typename std::enable_if<std::is_void<D>::value, RangeType>::type getNextChunkRange(const int tid) {
          assert(loaded);
          return chunkPartitioner.getNext(tid);
        }
        template <typename D = Derived>
        typename std::enable_if<!std::is_void<D>::value, RangeType>::type getNextChunkRange(const int tid) {
          return static_cast<Derived*>(this)->getNextChunkRangeImpl(tid);
        }




        /**
         *  performs memmap, and optionally preload the data into memory.
         *
         * @param r   range that will be mapped and loaded into memory (with buffering or without)
         */
        void load(const RangeType &range) throw (IOException)
        {
          // clean up any previous runs.
          //DEBUG("Loading");
          unload();

          //DEBUG("loading");

          mmap_range = range & fileRange;
          mmap_range.align_to_page(page_size);

          /// do the mem map
          mmap_data = map(mmap_range);

//////////////// check to see if there is enough memory for preloading.
///  since Preloading is a template param, can't choose a srcData type if Preloading is true but there is not enough memory.
/// so we have to assert it.

          if (Preloading) {
            /// check if can load the region into memory.

            /// check if we can preload.  use at most 1/20 of the available memory because of kmer expansion.
            struct sysinfo memInfo;
            sysinfo (&memInfo);

            if ((mmap_range.size() * sizeof(T)) > (memInfo.freeram * memInfo.mem_unit / 2)) { // use too much mem.  throw exception.  linux freeram is limited to 10 to 15 % of physical.
              std::stringstream ss;
              ss << "ERROR in file preloading: ["  << filename << "] mem required " << (mmap_range.size() * sizeof(T)) << " bytes; (1/20th) available: " << (memInfo.freeram * memInfo.mem_unit / 2) << " free: " << memInfo.freeram << " unit " << memInfo.mem_unit;
              throw IOException(ss.str());
            }
          }


          //DEBUG("mapped");
          // okay to use + since mmap_data is a pointer so it's a random access iterator.
          srcData.assign(mmap_data + (mmap_range.start - mmap_range.block_start), mmap_data + (mmap_range.end - mmap_range.block_start), mmap_range);

          if (Preloading)
          {
            unmap(mmap_data, mmap_range);
          }
          //DEBUG("loaded");
          loaded = true;

          recordSize = getRecordSize<Derived>(3);   // look through 3 records to see the max sizeof records.
          chunkPartitioner.configure(mmap_range, nthreads, std::max(chunkSize, 2 * recordSize));   // TODO: copy constructor works?
                // at least 2 x the record size for a chunk.

        }

        /**
         * unmap the file
         */
        void unload()
        {
          //DEBUG("Unloading");

          if (!Preloading)
          {
            if (mmap_data != nullptr && mmap_data != MAP_FAILED)
              unmap(mmap_data, mmap_range);
          }

          srcData.clear();
          mmap_data = nullptr;
          //DEBUG("unloaded");
          loaded = false;
          mmap_range = RangeType();
        }



        /**
         * get the loaded source data.
         * @return
         */
        DataType& getData() {
          assert(loaded);

          return srcData;
        }

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
        DataBlockType& getChunk(const int tid, const RangeType &chunkRange) {
          assert(loaded);
//          printf("chunkSize = %lu\n", chunkSize);
          RangeType r = chunkRange & mmap_range;

          auto s = srcData.begin();
          std::advance(s, (r.start - mmap_range.start));
          auto e = srcData.begin();
          std::advance(e, (r.end - mmap_range.start));


          dataBlocks[tid].assign(s, e, r);
//         DEBUG("read " << readLen << " elements, start at " << s << " [" << *startPtr << "] end at "<< e << " [" << *endPtr << "]");

//          if (startPtr == nullptr || endPtr == nullptr) {
//            std::cerr << "ERROR: file loader get chunk returning null ptrs. readlen = " << readLen << " start and end are " << s << "-" << e << " range is " << range << std::endl;
//          }
        return dataBlocks[tid];
      }




      protected:

        void init(const size_t &file_size) {
          page_size = sysconf(_SC_PAGE_SIZE);

          // compute the full mmap_range
          fileRange = RangeType(0, file_size / sizeof(T));   // range is in units of T

          /// open the file and get a handle.
          file_handle = open64(filename.c_str(), O_RDONLY);
          if (file_handle == -1)
          {
            std::stringstream ss;
            int myerr = errno;
            ss << "ERROR in file open: ["  << filename << "] error " << myerr << ": " << strerror(myerr);
            throw IOException(ss.str());
          }

          //DEBUG("file_loader initialized for " << filename);

          dataBlocks = new DataBlockType[nthreads];  // as many as there are threads, to allow delayed access to the datablocks.

          // partition the full mmap_range
          partitioner.configure(fileRange, nprocs, chunkSize);
        }


        // partitioner has chunk size. and current position.  at most as many as number of threads (specified in constructor) calling this function.
        template<typename D = Derived>
        typename std::enable_if<std::is_void<D>::value, SizeType>::type getRecordSize(const int iterations = 3) {
          return 1;
        }
        template<typename D = Derived>
        typename std::enable_if<!std::is_void<D>::value, SizeType>::type getRecordSize(const int iterations = 3) {
          return static_cast<Derived*>(this)->getRecordSizeImpl(iterations);
        }


        /**
         * map a region of the file to memory.
         * @param r
         * @return
         */
        InputIteratorType map(RangeType &r) throw (IOException) {

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
            throw IOException(ss.str());
          }

          //DEBUG("mapped");

          return result;
        }

        void unmap(InputIteratorType &d, RangeType &r) {

          munmap(d, (r.end - r.block_start) * sizeof(T));
          //DEBUG("unmapped");

        }

    };

  } /* namespace functional */
} /* namespace bliss */
#endif /* FILE_LOADER_HPP_ */
