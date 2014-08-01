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
     *          When loading the file, a L1Block is loaded with its own mmap call.
     *          The L1Blockss are consumed by distinct, independent reading processes, e.g. MPI processes
     *
     *          The second level partitioning generates "L2 BLOCKSs".
     *          When loading access a L2Block, the L2Block is generated from the parent L1Block, and therefore
     *            does not have its own mmap call.
     *          The L2Blocks in the a L1Block are consumed by distinct but cooperative computing processes, e.g. openmp threads.
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
     *          L1Partitioner is hard coded to BlockPartitioner at the moment
     *
     * @tparam  T                 type of each element read from file
     * @tparam  L1Buffering       bool indicating if L1 partition blocks should be buffered.  default to false
     * @tparam  L2Buffering       bool indicating if L2 partition blocks should be buffered.  default to true (to avoid contention between threads)
     * @tparam  L2PartitionerT    L2 partitioner, default to DemandDrivenPartitioner
     * @tparam  Derived           Derived (sub) class name.  For FileLoader, this is void.
     */
    template <typename T, bool L1Buffering = false, bool L2Buffering = true,
          typename L2PartitionerT = bliss::partition::DemandDrivenPartitioner<bliss::partition::range<size_t> >,
          typename Derived = void >
    class FileLoader
    {
      public:
        //========== type definitions.
        /**
         * @typedef RangeType
         * @brief   range object types
         */
        typedef bliss::partition::range<size_t>     RangeType;

        /**
         * @typedef PointerType
         * @brief   type for the raw data from the file.  element is of templated type T.
         */
        typedef T*                                  PointerType;
        using InputIteratorType = PointerType;


        /**
         * @typedef L1BlockType
         * @brief   DataBlock type for blocks from L1 partitioning.  blocks may be buffered or unbuffered depending on template parameter L1Buffering.
         */
        typedef typename std::conditional<L1Buffering,
                                          bliss::io::BufferedDataBlock<  PointerType, RangeType>,
                                          bliss::io::UnbufferedDataBlock<PointerType, RangeType> >::type                        L1BlockType;

        /**
         * @typedef L2BlockType
         * @brief   DataBlock type for blocks from L2 partitioning.  blocks may be buffered or unbuffered depending on template parameter L2Buffering.
         */
        typedef typename std::conditional<L2Buffering,
                                          bliss::io::BufferedDataBlock<  typename L1BlockType::iterator, RangeType>,
                                          bliss::io::UnbufferedDataBlock<typename L1BlockType::iterator, RangeType> >::type     L2BlockType;

      protected:
        /**
         * @typedef L1PartitionerT
         * @brief   Type of the Level 1 Partitioner to generate the range of the file to load from disk
         */
        typedef bliss::partition::BlockPartitioner<RangeType >  L1PartitionerT;

        //========== member variables

        //====== "constants"
        /**
         * @brief size of a disk page.
         */
        const size_t pageSize;


        /**
         * @brief size of a record in a file.
         * @note  conceptually a constant but needs to be set via a function call, so mutable.
         */
        mutable size_t recordSize;

        /// boolean indicating if file has been loaded
        bool loaded;

        //====== file loading
        /**
         * @brief name of file to load.
         * @note  FileLoader has the same filename for its lifetime.
         */
        //const std::string filename;

        /**
         * @brief file handle for the opened file
         * @note  as file loader is tied to a file, this is conceptually a constant, so mutable.
         */
        mutable int fileHandle;

        /**
         * @brief full range of the file, in units of data type  T size.  [0, file_size / sizeof(T) ).
         * @note  as file loader is tied to a file, this is conceptually a constant, so mutable.
         */
        mutable RangeType fileRange;  // full size of file

        //==== L1 Partitioning.
#if defined(USE_MPI)
        /// MPI communication object.  for use when we use 1 loader per MPI process.
        MPI_Comm comm;
#endif


        /**
         * @brief   id of loader that will be performing the memmap and reading of a subrange of the file.
         * @note    The mapping of processes/threads to loaders is flexible.  The usecase we have here is 1 to 1 mapping
         *           between MPI processes to loaders, so loaderId == mpi comm rank.
         *          Id is assigned to a process running a FileLoader, and is constant for the lifetime of the FileLoader
         */
        mutable int loaderId;

        /**
         * @brief Level 1 Partitioner, default to BlockPartitioner.
         * @details  produces L1 Block ranges that will be used to memmap the file.
         */
        bliss::partition::BlockPartitioner<RangeType > L1Partitioner;

        //==== memory mapping

        /// range in the file to be memmapped and read.  produced by L1Partitioner
        RangeType mmapRange;

        /**
         * @brief  memmapped data as raw pointer.  loaded from file using mmapRange as specification.
         * @details the pointer points to file data with position that is page aligned, as required by mmap.
         */
        T* mmapData;      // mem-mapped data, page aligned.  strictly internal

        //===   L1 Block
        /// L1Block instance that wraps mmapData to provide buffering and abstraction of interface.
        L1BlockType L1Block;


        //=== L2 Blocks and partitioner
        /**
         * @brief   The number of concurrent consuming threads that are processing the memmapped region of a file at L2 partition level
         * @note    The number of consuming threads may be different than the total number of L2Blocks, if Cyclic or DemandDriven Partitioners are used
         *           instead of BlockPartitioner
         *          The mapping of OS threads to consuming threads is flexible.  The use case we have here is 1 to 1 mapping
         *           between openmp threads to consuming threads.
         *
         *          the threads consuming the L2 blocks should all live in the same process since they share memmapped data
         *          The number of Consuming Threads is constant for the lifetime of the FileLoader
         */
        const size_t nConsumingThreads;

        /**
         * @brief   The size ofteh L2 Blocks.  User tunable.
         * @note    Constant for teh lifetime of the File Loader.
         */
        const size_t L2BlockSize;

        /**
         * @brief Level 2 Partitioner, default to DemandDrivenPartitioner.
         * @details  produces L2 Block ranges from the memmapped region.  The blocks are consumed by computing threads
         */
        L2PartitionerT L2Partitioner;

        /**
         * @brief L2 Block cache.  size of array is | nConsumingThreads |
         * @details The purpose of this cache is to allow slower computation to access the data, especially if buffered, without
         *          implicitly or explicitly copying the data to the consuming threads
         *          It also provide some thread safety as a consuming thread has its own slot in the array, and changes to that element
         *          is done by only that thread.
         */
        L2BlockType *L2Blocks;



      public:
        //================== Constructor and Destructor

        /// Removed default constructor
        FileLoader() = delete;

#if defined(USE_MPI)
        /**
         * @brief MPI enabled constructor.
         * @details   MPI Comm object is used to get nConcurrentLoaders, thus the L1 partitions, and to broadcast file size
         *
         *            The constructor first opens the file,
         *                then get the number of concurrent loaders and its loader id (comm_size, and comm_rank respectively)
         *                then configures the L1 partitioner
         *                and initializes the L2 block caching (per thread)
         *            The file is opened, not yet mmap (load).  this is done to allow sequencing of events by caller
         *            and subclasses.
         *
         * @note      Default behavior is to block-partition the file based on number of processors.
         *            To adjust the L1 block range, call getMMapRange, and modify it before calling load()
         *
         *            It is necessary to call "load" 1 time before using accessing the data.
         *
         *            Range is in units of T, overlap is included in range.end
         *
         * @param _filename       input file name
         * @param _comm           MPI Communicator (defined only if CMake is configured with MPI).
         *                        default is MPI_COMM_WORLD, which means all nodes participate.
         * @param _nThreads       number of threads to consume L2 blocks
         * @param _L2BlockSize    size of each L2Block.
         *
         * TODO: test on remotely mounted file system.
         *
         */
        FileLoader(const std::string &_filename, const MPI_Comm& _comm, const size_t _nThreads = 1, const size_t _L2BlockSize = 1) throw (bliss::io::IOException)
            : pageSize(sysconf(_SC_PAGE_SIZE)), recordSize(1), loaded(false), fileHandle(-1), fileRange(), comm(_comm),
              L1Partitioner(), mmapRange(), mmapData(nullptr), L1Block(),
              nConsumingThreads(_nThreads), L2BlockSize(_L2BlockSize), L2Partitioner(), L2Blocks(nullptr)
        {
          // TODO: remove asserts.
          assert(_filename.length() > 0);
          assert(_nThreads > 0);
          assert(_L2BlockSize > 0);

          // open the file
          fileHandle = openFile(_filename);

          // get from the communicator the number of concurrent loaders, and the id of the current loader.
          int nConcurrentLoaders = 0;
          MPI_Comm_rank(comm, &loaderId);
          MPI_Comm_size(comm, &nConcurrentLoaders);

//          // get the file size.
//          size_t file_size = 0;
//          if (loaderId == 0)
//          {
//            size_t file_size = getFileSize(fileHandle);
//          }
//          // TODO: check if we can avoid the broadcast.
//          if (nConcurrentLoaders > 1) {
//            // broadcast file_size to all
//            MPI_Bcast(&file_size, 1, MPI_UNSIGNED_LONG_LONG, 0, comm);
//          }

          // get the file size.  all processes do this, since all have to open the file and read from it anyways.
          // configure the L1 partitioner
          configL1(getFileSize(fileHandle), nConcurrentLoaders);

          // init the L2 block cache.
          initL2();
        }
#endif

        /**
         * @brief NON-MPI constructor.
         * @details   The constructor first opens the file,
         *                then configures the L1 partitioner with number of concurrent loaders (e.g. nthreads, or mpi comm_size)
         *                and its loader id (e.g. thread id, or mpi rank).
         *                and initializes the L2 block caching (per thread)
         *            The file is opened, not yet mmap (load).  this is done to allow sequencing of events by caller
         *            and subclasses.
         *
         * @note      Default behavior is to block-partition the file based on number of processors.
         *            To adjust the L1 block range, call getMMapRange, and modify it before calling load()
         *
         *            It is necessary to call "load" 1 time before using accessing the data.
         *
         *            Range is in units of T, overlap is included in range.end
         *
         * @param _filename       input file name
         * @param _comm           MPI Communicator (defined only if CMake is configured with MPI).
         *                        default is MPI_COMM_WORLD, which means all nodes participate.
         * @param _nThreads       number of threads to consume L2 blocks
         * @param _L2BlockSize    size of each L2Block.
         *
         * TODO: test on remotely mounted file system.
         *
         */
        FileLoader(const std::string &_filename,
                            const int _nConcurrentLoaders = 1, const int _loaderId = 0,
                            const size_t _nThreads = 1, const size_t _L2BlockSize = 1 ) throw (bliss::io::IOException)
            : pageSize(sysconf(_SC_PAGE_SIZE)), recordSize(1), loaded(false), fileHandle(-1), fileRange(),
              loaderId(_loaderId), L1Partitioner(), mmapRange(), mmapData(nullptr), L1Block(),
              nConsumingThreads(_nThreads), L2BlockSize(_L2BlockSize), L2Partitioner(), L2Blocks(nullptr)
        {
          assert(_filename.length() > 0);
          assert(_nThreads > 0);
          assert(_L2BlockSize > 0);
          assert(_loaderId >= 0);
          assert(_nConcurrentLoaders > _loaderId);

          // open the file
          fileHandle = openFile(_filename);

          // get the file size and configure L1 partitioner
          configL1(getFileSize(fileHandle), _nConcurrentLoaders);

          // initialize L2 block cache.
          initL2();
        }



        /// Removed default copy constructor
        FileLoader(const FileLoader<T, L1Buffering, L2Buffering, L2PartitionerT, Derived>& other) = delete;
        /// Removed default copy assignment operator
        FileLoader& operator=(const FileLoader<T, L1Buffering, L2Buffering, L2PartitionerT, Derived>& other) = delete;


        /**
         * @brief move constructor.  here to ensure that the DataBlocks are moved.
         * @param other FileLoader to move
         */
        FileLoader(FileLoader<T, L1Buffering, L2Buffering, L2PartitionerT, Derived>&& other) :
          pageSize(other.pageSize), recordSize(other.recordSize), loaded(other.loaded),
          fileHandle(other.fileHandle), fileRange(other.fileRange),
          loaderId(other.loaderId), L1Partitioner(other.L1Partitioner), mmapRange(other.mmapRange), mmapData(other.mmapData),
          L1Block(std::move(other.L1Block)), nConsumingThreads(other.nConsumingThreads), L2BlockSize(other.L2BlockSize),
          L2Partitioner(other.L2Partitioner), L2Blocks(other.L2Blocks)


        {
          other.pageSize = 1;
          other.fileHandle = -1;
          other.fileRange = RangeType();
          other.mmapData = nullptr;
          other.mmapRange = RangeType();

          other.loaded = false;
          other.nConsumingThreads = 1;
          other.L2BlockSize = 1;
          other.loaderId = 0;
          other.L2Blocks = nullptr;
          other.L1Partitioner = L1PartitionerT();
          other.L2Partitioner = L2PartitionerT();
          other.recordSize = 1;

#if defined(USE_MPI)
          comm        = other.comm;
          other.comm = MPI_COMM_NULL;
#endif
        }

        /**
         * @brief move assignement operator.  here to ensure that the DataBlocks are moved.
         *
         * @param other   FileLoader to move
         * @return        updated object with moved data
         */
        FileLoader& operator=(FileLoader<T, L1Buffering, L2Buffering, L2PartitionerT, Derived>&& other) {
          if (this != &other) {
            // remove the old
            if (L2Blocks != nullptr) delete [] L2Blocks;
            // unload and close my file before copying in the other one.
            unload();
            // close my file
            closeFile(fileHandle);

            pageSize   = other.pageSize;                         other.pageSize = 1;
            fileHandle = other.fileHandle;                       other.fileHandle = -1;
            fileRange   = other.fileRange;                         other.fileRange = RangeType();
            mmapData   = other.mmapData;                         other.mmapData = nullptr;
            mmapRange  = other.mmapRange;                        other.mmapRange = RangeType();
            L1Block     = std::move(other.L1Block);
            loaded      = other.loaded;                            other.loaded = false;
            nConsumingThreads    = other.nConsumingThreads;                          other.nConsumingThreads = 1;
            L2BlockSize   = other.L2BlockSize;                         other.L2BlockSize = 1;
            loaderId        = other.loaderId;                              other.loaderId = 0;
  #if defined(USE_MPI)
            comm        = other.comm;                              other.comm = MPI_COMM_NULL;
  #endif
            L2Blocks  = other.L2Blocks;                        other.L2Blocks = nullptr;
            L1Partitioner = other.L1Partitioner;                       other.L1Partitioner = bliss::partition::BlockPartitioner<RangeType>();
            L2Partitioner = other.L2Partitioner;             other.L2Partitioner = L2PartitionerT();
            recordSize  = other.recordSize;                        other.recordSize = 1;
          }
          return *this;
        }


        /**
         * @brief default destructor.  unloads the data and closes the file
         */
        virtual ~FileLoader()
        {
          if (L2Blocks != nullptr) delete [] L2Blocks;

          unload();

          closeFile(fileHandle);
        }

        /**
         * @brief     Get the L2 partition Block Size.
         * @return    L2 block size
         */
        const size_t& getL2BlockSize() const {
          return L2BlockSize;
        }

        /**
         * @brief return the full range for this file.  (in units of data type T)
         * @return  the range that spans the entire file.
         */
        const RangeType& getFileRange() const {
          return fileRange;
        }

        /**
         * @brief  return the range for this file loader's memory map.  (in units of data type T).  valid only after loading.
         * @return  the range object spanning the entire memory mapped region.
         */
        const RangeType& getMMappedRange() const {
          //TODO: remove assert.
          assert(loaded);
          return mmapRange;
        }




        /**
         * @brief memory maps the specified region for reading,
         * @details  The input range SHOULD come from the call to getNextL1BlockRange
         *           optionally buffer the data into memory.
         *
         * @param r   range that will be mapped and loaded into memory (with buffering or without)
         */
        void load(const RangeType &range) throw (bliss::io::IOException)
        {
          // clean up any previous runs.
          unload();

          // make sure the mmapRange is within the file range, and page align it.
          mmapRange = RangeType::intersect(range, fileRange);
          size_t block_start = RangeType::align_to_page(mmapRange.start, pageSize);

          // map the region of file to memory
          mmapData = map(mmapRange);

          // check to see if there is enough memory for preloading.
          // if not, then we can't use Buffering even if L1Block chooses Buffering, so throw exception
          if (L1Buffering) {
            // check if we can buffer.  use at most 1/20 of the available memory because of kmer data size expansion.
            struct sysinfo memInfo;
            sysinfo (&memInfo);

            // linux freeram is limited to 10 to 15 % of physical, hence we just divide by 2, and that gives us 1/20 of the available memory
            if ((mmapRange.size() * sizeof(T)) > (memInfo.freeram * memInfo.mem_unit / 2)) {
              // not enough memory available
              std::stringstream ss;
              ss << "ERROR in file buffering: requires " << (mmapRange.size() * sizeof(T)) <<
                  " bytes; but (1/20th) available: " << (memInfo.freeram * memInfo.mem_unit / 2) <<
                  " free: " << memInfo.freeram << " unit " << memInfo.mem_unit;
              throw IOException(ss.str());
            }
          }

          // now assign the data to the L1Block object. (if buffering, will copy).
          // use "+" since mmapData is a pointer so it's a random access iterator.
          L1Block.assign(mmapData + (mmapRange.start - block_start), mmapData + (mmapRange.end - block_start), mmapRange);

          // if we are buffering, then it's okay to unmap the data, as we already have a copy
          if (L1Buffering)
          {
            unmap(mmapData, mmapRange);
          }

          // mark data as loaded.
          loaded = true;

          // now configure L2 partitioner to traverse this data.
          configL2();
        }

        /**
         * @brief   unloads the memmapped region of a file from memory
         * @details if buffering, clear the data
         *          if not buffering, then unmemmap the file region from memory
         */
        void unload()
        {
          // if the partition was not buffered,
          if (!L1Buffering)
          {
            // and data has been memmapped,
            if (mmapData != nullptr && mmapData != MAP_FAILED)
              // then we unmap
              unmap(mmapData, mmapRange);
          } // else the partition was preloading, so already unmapped

          // clean up L1Block
          L1Block.clear();

          // reset the data pointer, loaded flag, and the variable storing the memmapped range.
          mmapData = nullptr;
          loaded = false;
          mmapRange = RangeType();
        }


        /**
         * @brief reset the L1 partitioner so it can be reused.
         */
        void resetL1Partitioner() {
          L1Partitioner.reset();
        }


HERE!!!


        /**
         * get the Partitioned Range.  this method allows calling Derived class' version, to further refine the partitioned range.
         * have a default.  if default is -1, then use loaderId.  not in inner loop so okay to check param value each call
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
            return L1Partitioner.getNext(loaderId);
          else
            return L1Partitioner.getNext(pid);
        }
        template<typename D = Derived>
        typename std::enable_if<!std::is_void<D>::value, RangeType>::type getNextPartitionRange(const int pid = -1) {
          // use the derived one.
          if (pid < 0)
            return static_cast<Derived*>(this)->getNextPartitionRangeImpl(loaderId);
          else
            return static_cast<Derived*>(this)->getNextPartitionRangeImpl(pid);
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


        void resetChunkRange() {
          assert(loaded);
          L2Partitioner.reset();
        }

        // L1Partitioner has chunk size. and current position.  at most as many as number of threads (specified in constructor) calling this function.
        // not using a default tid because this function could be called a lot.
        /**
         *
         * @param tid
         * @return
         */
        template <typename D = Derived>
        typename std::enable_if<std::is_void<D>::value, RangeType>::type getNextChunkRange(const int tid) {
          assert(loaded);
          return L2Partitioner.getNext(tid);
        }
        template <typename D = Derived>
        typename std::enable_if<!std::is_void<D>::value, RangeType>::type getNextChunkRange(const int tid) {
          return static_cast<Derived*>(this)->getNextChunkRangeImpl(tid);
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
         * @brief   configure the L1 Partitioning
         * @details This function does the following:
         *            construct the fileRange from file size
         *            configure the L1Partitioner.
         *
         *          file_size is passed in because for the MPI version of c'tor, the file_size was broadcast first
         * @param file_size   size of the entire file being read.
         */
        void configL1(const size_t &file_size, const int& nConcurrentLoaders) {

          // compute the full range of the file
          fileRange = RangeType(0, file_size / sizeof(T));   // range is in units of T

          // Set up the partitioner to block partition the full file range.  Overlap is 0.
          L1Partitioner.configure(fileRange, nConcurrentLoaders);
        }

        /**
         * @brief   initializes the L2 Block caching
         * @details This function does the following:
         *            allocating the internal cache of DataBlocks
         */
        void initL2() {
          // allocate cache for DataBlocks, for persistent access by threads, so we only need as many as there are threads.
          L2Blocks = new L2BlockType[nConsumingThreads];
        }

        /**
         * @brief configure the L2 Partitioning
         * @details: This function does the following
         *            get the approximate record size
         *            configure L2 partitioner for the memmapped range, number of threads, and L2Block Size.
         */
        void configL2() {
          // look through 3 records to see the approximate size of records.
          recordSize = getRecordSize<Derived>(3);

          // update the L2BlockSize to at least 2x record size.
          L2BlockSize = std::max(L2BlockSize, 2 * recordSize);

          // then configure the partitinoer to use the right number of threads.
          L2Partitioner.configure(mmapRange, nConsumingThreads, L2BlockSize);
        }


        /**
         * @brief  opens a file and return the file handle.  Throws IOException if can't open.
         * @param fn  name of file to open
         * @return    file handle/descriptor for the file
         */
        int openFile(const std::string& fn) throw (bliss::io::IOException) {
          // open the file and get a handle.
          int output = open64(fn.c_str(), O_RDONLY);
          if (output == -1)
          {
            // if open failed, throw exception.
            std::stringstream ss;
            int myerr = errno;
            ss << "ERROR in file open: ["  << fn << "] error " << myerr << ": " << strerror(myerr);
            throw IOException(ss.str());
          }
          return output;
        }

        /**
         * @brief closes a file
         * @param fd  the file descriptor of the file to close.
         */
        void closeFile(int &fd) {
          if (fd != -1) {
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
        size_t getFileSize(const int& fd) throw (bliss::io::IOException) {
          struct stat64 filestat;

          // get the file state
          int ret = fstat64(fd, &filestat);

          // handle any error
          if (ret < 0) {
            throw IOException("ERROR in file size detection");
          }

          // return file size.
          return static_cast<size_t>(filestat.st_size);
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
