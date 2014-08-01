/**
 * @file    file_loader.hpp
 * @ingroup bliss::io
 * @author  tpan
 * @brief   contains a generic FileLoader class to perform distributed and concurrent file access.
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
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

#include <cassert>

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
#include "partition/partitioner.hpp"
#include "io/data_block.hpp"
#include "io/io_exception.hpp"
#include "utils/logging.h"


namespace bliss
{
  namespace io
  {

    // TODO: be more specific about meaning of things, e.g. chunk,  and rename functions approriately.

    /**
     * @class FileLoader
     * @brief Opens a file (or a portion of it) via memmap, and optionally copies into memory.
     * @details FileLoader is similar to a container: can return size, begin iterator and end iterator.
     *
     *          Memmap maps a file to a memory address range, and allows direct, random access to the
     *          file.  it utilizes the same mechanism as virtual memory paging, so it's efficient.
     *          Performance of mmap is good generally compared to traditional open-read-close-process
     *
     *          Using iterators allows streaming processing, which may be more efficient.
     *
     *          File loading happens with a 2 level partitioning, in order to support multiprocess (mpi)
     *          and multithreaded (openmp) concurrent file access.
     *
     *          The first level partitioning generates "L1 BLOCKSs".
     *          When loading the file, a L1Block is loaded with its own mmap operation
     *          The L1Blockss are consumed by distinct reading processes, e.g. MPI processes
     *
     *          The second level partitioning generates "L2 BLOCKSs".
     *          When loading access a L2Block, the L2Block is generated from the parent L1Block, and therefore
     *            does not have its own mmap call.
     *          The L2Blocks in the a L1Block are consumed by distinct computing processes, e.g. openmp threads.
     *
     *          Both L1 and L2 Blocks are specified using Range, and the data is accessed through DataBlocks,
     *          which optionally provides buffering at the L1 and/or L2 Block level.
     *
     *          Having a 2 level partitioning provides flexibility and encourages sequential file access (for L2Blocks)
     *          It is not necessary that L2Block be used, however, as L1 Blocks can be used to traverse the data as well.
     *
     *          Any finer grain partitioning, such as partitioning by record, should be done by subclasses of FileLoader.
     *
     *  Usage:
     *    instantiate file_loader() - // default partitioning is based on number of MPI processes.
     *                                // but user can also specify the number of partitions.
     *    get range for a L1Block     // allows iteration over partition ids, mapping of process to partition other than 1 to 1.
     *    load using the range        // actually open the file and memmap
     *
     *    get the iterators (begin(),end()) and do some work.  (L1 Block traversal)
     *    or
     *    get a L2Block  via getNextChunk(), then get the iterators (begin(),end()) and do some work.  (L2 Block traversal)
     *
     *    unload                       // mem unmap
     *    destroy file_loader          // close the file.
     *
     *
     * @note:   For real data,  mmap is better for large files and limited memory, while
     *              preloading is better for smaller files and/or large amount of memory.
     *          FileLoader is the base class in a CRTP static polymorphism pattern.
     *
     *          Memmap require page aligned addresses.  internally, this class has pointers to the
     *          start of the range of interest, not to the page aligned starting position.
     *          DataBlocks therefore has at the "begin" iterator the first element OF INTEREST in the mmapped range,
     *            not the first element mapped.
     *
     */
    template <typename T, bool Buffering = true, bool Preloading = false,
          typename L2Partitioner = bliss::partition::DemandDrivenPartitioner<bliss::partition::range<size_t> >,
          typename Derived = void >
    class FileLoader
    {


      public:
        /////// type defintions.
        typedef bliss::partition::range<size_t>     RangeType;
        typedef T*                                  PointerType;

        // internal DataBlock type.  uses T* as iterator type.
        typedef typename std::conditional<Preloading,
                                          bliss::io::BufferedDataBlock<  PointerType, RangeType>,
                                          bliss::io::UnbufferedDataBlock<PointerType, RangeType> >::type                        L1BlockType;

        typedef typename std::conditional<Buffering,
                                          bliss::io::BufferedDataBlock<  typename L1BlockType::iterator, RangeType>,
                                          bliss::io::UnbufferedDataBlock<typename L1BlockType::iterator, RangeType> >::type     L2BlockType;


      protected:
        ////// member variables

        /// "constants"
        mutable size_t pageSize;
        mutable size_t recordSize;

        /// status
        bool loaded;
        bool preloaded;

        /// file loading
        mutable std::string filename;
        mutable int fileHandle;      // file handle
        mutable RangeType fileRange;  // offset in file from where to read


        /// L1 partitioner
#if defined(USE_MPI)
        MPI_Comm comm;
#endif
        int nprocs;           // for partitioning.
        int rank;             // for partitioning.
        bliss::partition::BlockPartitioner<RangeType > L1partitioner;   // block partitioning is the default.

        /// memory mapping
        RangeType mmapRange;      // offset in file from where to read
        T* mmapData;      // mem-mapped data, page aligned.  strictly internal

        /// L1 Block
        L1BlockType L1Block;


        /// L2 Blocks and partitioner
        int nthreads;
        size_t L2BlockSize;
        L2Partitioner L2partitioner;                          // construct when mmapRange change, on load.
        L2BlockType *L2Blocks;




        /////////////// Constructor and Destructor
      public:

        /// defining move constructor will disallow automatic copy constructor.
        //  this is to prevent copying the buffers, to keep a single handle on the input file.
        FileLoader(FileLoader<T, Buffering, Preloading, L2Partitioner, Derived>&& other) {
          pageSize   = other.pageSize;                         other.pageSize = 1;
          fileHandle = other.fileHandle;                       other.fileHandle = -1;
          filename.swap(other.filename);
          fileRange   = other.fileRange;                         other.fileRange = RangeType();
          mmapData   = other.mmapData;                         other.mmapData = nullptr;
          mmapRange  = other.mmapRange;                        other.mmapRange = RangeType();
          L1Block     = std::move(other.L1Block);
          loaded      = other.loaded;                            other.loaded = false;
          preloaded   = other.preloaded;                         other.preloaded = false;
          nthreads    = other.nthreads;                          other.nthreads = 1;
          L2BlockSize   = other.L2BlockSize;                         other.L2BlockSize = 1;
          nprocs      = other.nprocs;                            other.nprocs = 1;
          rank        = other.rank;                              other.rank = 0;
#if defined(USE_MPI)
          comm        = other.comm;                              other.comm = MPI_COMM_NULL;
#endif
          L2Blocks  = other.L2Blocks;                        other.L2Blocks = nullptr;
          L1partitioner = other.L1partitioner;                       other.L1partitioner = bliss::partition::BlockPartitioner<RangeType>();
          L2partitioner = other.L2partitioner;             other.L2partitioner = L2Partitioner();
          recordSize  = other.recordSize;                        other.recordSize = 1;
        }

        /// defining move assignment will disallow automatic copy assignment operator
        FileLoader& operator=(FileLoader<T, Buffering, Preloading, L2Partitioner, Derived>&& other) {
          if (this != &other) {
            //DEBUG("DESTROY");
            if (L2Blocks != nullptr) delete [] L2Blocks;

            unload();

            //printf("unloading complete.\n");
            if (fileHandle != -1) {
              close(fileHandle);
              fileHandle = -1;
            }

            pageSize   = other.pageSize;                         other.pageSize = 1;
            fileHandle = other.fileHandle;                       other.fileHandle = -1;
            filename.swap(other.filename);
            fileRange   = other.fileRange;                         other.fileRange = RangeType();
            mmapData   = other.mmapData;                         other.mmapData = nullptr;
            mmapRange  = other.mmapRange;                        other.mmapRange = RangeType();
            L1Block     = std::move(other.L1Block);
            loaded      = other.loaded;                            other.loaded = false;
            preloaded   = other.preloaded;                         other.preloaded = false;
            nthreads    = other.nthreads;                          other.nthreads = 1;
            L2BlockSize   = other.L2BlockSize;                         other.L2BlockSize = 1;
            nprocs      = other.nprocs;                            other.nprocs = 1;
            rank        = other.rank;                              other.rank = 0;
  #if defined(USE_MPI)
            comm        = other.comm;                              other.comm = MPI_COMM_NULL;
  #endif
            L2Blocks  = other.L2Blocks;                        other.L2Blocks = nullptr;
            L1partitioner = other.L1partitioner;                       other.L1partitioner = bliss::partition::BlockPartitioner<RangeType>();
            L2partitioner = other.L2partitioner;             other.L2partitioner = L2Partitioner();
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
        FileLoader(const std::string &_filename, const MPI_Comm& _comm, const int _nThreads = 1, const size_t _L2BlockSize = 1) throw (bliss::io::IOException)
            : fileHandle(-1), filename(_filename), fileRange(), mmapData(nullptr),
              mmapRange(), loaded(false), preloaded(false),
              nthreads(_nThreads), L2BlockSize(_L2BlockSize), comm(_comm), L2Blocks(nullptr), recordSize(1)
        {
          assert(filename.length() > 0);
          assert(_nThreads > 0);
          assert(_L2BlockSize > 0);

          // get the processor rank and nprocessors.
          MPI_Comm_rank(comm, &rank);
          MPI_Comm_size(comm, &nprocs);

          //DEBUG("CONSTRUCT");

          /// get the file size.
          size_t file_size = 0;
          if (rank == 0)
          {
            struct stat filestat;
            int ret = stat(filename.c_str(), &filestat);

            if (ret < 0) {

              std::stringstream ss;
              ss << "ERROR in file size detection: ["  << filename << "] error ";
              throw IOException(ss.str());
            }

            file_size = static_cast<size_t>(filestat.st_size);
//            std::cerr << "file size is " << file_size;
//            std::cerr << " block size is " << filestat.st_blksize;
//            std::cerr << " sysconf block size is " << pageSize << std::endl;
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
                            const int _nThreads = 1, const size_t _L2BlockSize = 1 ) throw (bliss::io::IOException)
            : fileHandle(-1), filename(_filename), fileRange(), mmapData(nullptr),
              mmapRange(), loaded(false), preloaded(false),
              nthreads(_nThreads), L2BlockSize(_L2BlockSize), nprocs(_nProcs), rank(_rank), L2Blocks(nullptr), recordSize(1)
        {
          //DEBUG("CONSTRUCT");

          assert(filename.length() > 0);
          assert(_nThreads > 0);
          assert(_L2BlockSize > 0);
          assert(_rank >= 0);
          assert(_nProcs > _rank);


          /// get the file size.
          size_t file_size = 0;
          struct stat filestat;
          int ret = stat(filename.c_str(), &filestat);

          if (ret < 0) {

            std::stringstream ss;
            ss << "ERROR in file size detection: ["  << filename << "] error ";
            throw IOException(ss.str());
          }

          file_size = static_cast<size_t>(filestat.st_size);
//            std::cerr << "file size is " << file_size;
//            std::cerr << " block size is " << filestat.st_blksize;
//            std::cerr << " sysconf block size is " << pageSize << std::endl;

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
          if (L2Blocks != nullptr) delete [] L2Blocks;

          unload();

          //printf("unloading complete.\n");
          if (fileHandle != -1) {
            close(fileHandle);
            fileHandle = -1;
          }
        }



        ////// PUBLIC METHODS

        const size_t& getChunkSize() const {
          return L2BlockSize;
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
          return mmapRange;
        }


        void resetPartitionRange() {
          L1partitioner.reset();
        }

        /**
         * get the Partitioned Range.  this method allows calling Derived class' version, to further refine the partitioned range.
         * have a default.  if default is -1, then use rank.  not in inner loop so okay to check param value each call
         *
         * @param pid
         * @return
         *
         * TODO: overload this to avoid "-1"
         */
        template<typename D = Derived>
        typename std::enable_if<std::is_void<D>::value, RangeType>::type getNextPartitionRange(const int pid = -1) {
          // if just FileLoader calls this.
          if (pid < 0)
            return L1partitioner.getNext(rank);
          else
            return L1partitioner.getNext(pid);
        }
        template<typename D = Derived>
        typename std::enable_if<!std::is_void<D>::value, RangeType>::type getNextPartitionRange(const int pid = -1) {
          // use the derived one.
          if (pid < 0)
            return static_cast<Derived*>(this)->getNextPartitionRangeImpl(rank);
          else
            return static_cast<Derived*>(this)->getNextPartitionRangeImpl(pid);
        }



        void resetChunkRange() {
          assert(loaded);
          L2partitioner.reset();
        }

        // L1partitioner has chunk size. and current position.  at most as many as number of threads (specified in constructor) calling this function.
        // not using a default tid because this function could be called a lot.
        /**
         *
         * @param tid
         * @return
         */
        template <typename D = Derived>
        typename std::enable_if<std::is_void<D>::value, RangeType>::type getNextChunkRange(const int tid) {
          assert(loaded);
          return L2partitioner.getNext(tid);
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

          mmapRange = RangeType::intersect(range, fileRange);
          size_t block_start = mmapRange.align_to_page(pageSize);

          /// do the mem map
          mmapData = map(mmapRange);

//////////////// check to see if there is enough memory for preloading.
///  since Preloading is a template param, can't choose a L1Block type if Preloading is true but there is not enough memory.
/// so we have to assert it.

          if (Preloading) {
            /// check if can load the region into memory.

            /// check if we can preload.  use at most 1/20 of the available memory because of kmer expansion.
            struct sysinfo memInfo;
            sysinfo (&memInfo);

            if ((mmapRange.size() * sizeof(T)) > (memInfo.freeram * memInfo.mem_unit / 2)) { // use too much mem.  throw exception.  linux freeram is limited to 10 to 15 % of physical.
              std::stringstream ss;
              ss << "ERROR in file preloading: ["  << filename << "] mem required " << (mmapRange.size() * sizeof(T)) << " bytes; (1/20th) available: " << (memInfo.freeram * memInfo.mem_unit / 2) << " free: " << memInfo.freeram << " unit " << memInfo.mem_unit;
              throw IOException(ss.str());
            }
          }


          //DEBUG("mapped");
          // okay to use + since mmapData is a pointer so it's a random access iterator.
          L1Block.assign(mmapData + (mmapRange.start - block_start), mmapData + (mmapRange.end - block_start), mmapRange);

          if (Preloading)
          {
            unmap(mmapData, mmapRange);
          }
          //DEBUG("loaded");
          loaded = true;

          recordSize = getRecordSize<Derived>(3);   // look through 3 records to see the max sizeof records.

          // update the L2BlockSize to at least 2x record size.  uses move assignment operator
          L2BlockSize = std::max(L2BlockSize, 2 * recordSize);
          // then configure the partitinoer to use the right number of threads.
          L2partitioner.configure(mmapRange, nthreads, L2BlockSize);

        }

        /**
         * @brief   unloads the memmap of (a partition of) the file from memory
         * @details if preloaded, clear the data
         *          if not preloaded, then unmemmap the file region from memory
         */
        void unload()
        {
          // if the partition was not preloaded,
          if (!Preloading)
          {
            // and mmap did not fail
            if (mmapData != nullptr && mmapData != MAP_FAILED)
              unmap(mmapData, mmapRange);
            // then we unmap

          } // else the partition was preloaded, so already unmapped

          // clean up variables.
          L1Block.clear();
          mmapData = nullptr;
          loaded = false;
          mmapRange = RangeType();
        }

        /**
         * @brief   get the Partition that's been memmapped/loaded.
         * @details the data can then be used directly, using the iterator accessors of DataBlock
         * @return  DataBlock (buffered or unbuffered) wrapping the data mapped/read from the file.
         */
        L1BlockType& getL1Data() {

          assert(loaded);

          return L1Block;
        }


        /**
         * @brief   get the DataBlock from the Partition
         * @param tid
         * @param chunkRange
         * @return
         */
        L2BlockType& getL2DataForRange(const int tid, const RangeType &L2BlockRange) {

          assert(loaded);

          RangeType r = RangeType::intersect(L2BlockRange, mmapRange);

          auto s = L1Block.begin();
          std::advance(s, (r.start - mmapRange.start));
          auto e = s;
          std::advance(e, (r.size()));


          L2Blocks[tid].assign(s, e, r);
          return L2Blocks[tid];
      }




      protected:

        /**
         * @brief   initializes the FileLoader
         * @details This function does the following:
         *            get the page size
         *            construct the fileRange from file size
         *            opening the file
         *            allocating the internal cache of DataBlocks
         *            configure the L1partitioner.
         * @param file_size   size of the entire file being read.
         */
        void init(const size_t &file_size) {
          // get the page size
          pageSize = sysconf(_SC_PAGE_SIZE);

          // compute the full range of the file
          fileRange = RangeType(0, file_size / sizeof(T));   // range is in units of T

          // open the file and get a handle.
          fileHandle = open64(filename.c_str(), O_RDONLY);
          if (fileHandle == -1)
          {
            // if open failed, throw exception.
            std::stringstream ss;
            int myerr = errno;
            ss << "ERROR in file open: ["  << filename << "] error " << myerr << ": " << strerror(myerr);
            throw IOException(ss.str());
          }

          // allocate cache for DataBlocks, for persistent access by threads, so we only need as many as there are threads.
          L2Blocks = new L2BlockType[nthreads];

          // Set up the partitioner to block partition the full file range.  Overlap is 0.
          L1partitioner.configure(fileRange, nprocs);
        }

        /**
         * @brief   compute the approximate size of a data record in the file by reading a few records.
         * @details For FileLoader, which does not treat the file as a collection of records, size returned is 1.
         *          This function delegates to the subclass' implementation, as subclass will have knowledge of how to
         *          parse a record.
         *          This method provides a hint of the size of a record, which is useful for caching
         *          and determining size of a DataBlock.
         *
         * @param count    number of records to read to compute the approximation. default = 3
         * @return  approximate size of a record.
         */
        template<typename D = Derived>
        typename std::enable_if<std::is_void<D>::value, size_t>::type getRecordSize(const int count = 3) {
          // this method is enabled if Derived is void, which indicates self.
          return 1;
        }

        /**
         * @brief   compute the approximate size of a data record in the file by reading a few records.
         * @details For subclasses of FileLoader, this function delegates to the subclass' implementation,
         *          as subclass will have knowledge of how to parse a record.
         *          This method provides a hint of the size of a record, which is useful for caching
         *          and determining size of a DataBlock.
         * @note    at most as many threads as specified in the c'tor will call this function.
         * @param count    number of records to read to compute the approximation.  default = 3.
         * @return  approximate size of a record.
         */
        template<typename D = Derived>
        typename std::enable_if<!std::is_void<D>::value, size_t>::type getRecordSize(const int count = 3) {
          return static_cast<Derived*>(this)->getRecordSizeImpl(count);
        }


        /**
         * @brief     map the specified portion of the file to memory.
         * @param r   range specifying the portion of the file to map
         * @return    memory address (pointer) to where the data is mapped.
         */
        PointerType map(RangeType &r) throw (IOException) {

          /// memory map.  requires that the starting position is block aligned.
          size_t block_start = r.align_to_page(pageSize);

          // NOT using MAP_POPULATE.  it slows things done when testing on single node.
          PointerType result = (PointerType)mmap64(nullptr, (r.end - block_start ) * sizeof(T),
                                     PROT_READ,
                                     MAP_PRIVATE, fileHandle,
                                     block_start * sizeof(T));

          // if mmap failed,
          if (result == MAP_FAILED)
          {
            // clean up.
            if (fileHandle != -1)
            {
              close(fileHandle);
              fileHandle = -1;
            }

            // print error through exception.
            std::stringstream ss;
            int myerr = errno;
            ss << "ERROR in mmap: " << myerr << ": " << strerror(myerr);
            throw IOException(ss.str());
          }

          return result;
        }

        /**
         * @brief unmaps a file region from memory
         * @param d   The pointer to the memory address
         * @param r   The range that was mapped.
         */
        void unmap(PointerType &d, RangeType &r) {
          munmap(d, (r.end - r.align_to_page(pageSize)) * sizeof(T));
        }

    };

  } /* namespace functional */
} /* namespace bliss */
#endif /* FILE_LOADER_HPP_ */
