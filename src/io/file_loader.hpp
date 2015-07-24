/**
 * @file    file_loader.hpp
 * @ingroup io
 * @author  Tony Pan <tpan7@gatech.edu>
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
   template<typename FileType>
   struct PositionType;

    /**
     * @class FileLoader
     * @brief Opens a file (or a portion of it) via memmap, and optionally copies(portion of) data into memory.
     * @details FileLoader is similar to an iterator but does not use the same API.  It has a set of getNext__ functions
     *            for traversing the file in BLOCKS.
     *
     *            The decision not to follow iterator API is due to thread support:  the thread id would need to be
     *            passed in to the iterator operators, thus would already create a departure from the iterator API.
     *            Also, the increment and dereference operations need to be bundled together to support atomicity, which
     *            also is different than the standard iterator approach.
     *
     *            FileLoader "getNext__" functions returns a DataBlock.  DataBlock has begin and end iterators for traversing
     *              the individual data elements, and can provide buffering.
     *
     *        2-Level Partitioning:
     *          FileLoader uses a 2 level partitioning, in order to support multiprocess (mpi)
     *          and multithreaded (openmp) concurrent file access.
     *
     *          The first level partitioning generates "L1 BLOCKs".
     *          L1Block of a file is loaded from disk/network with its own mmap call.
     *          This makes L1Blocks suited for distinct, independent FileLoader processes, e.g. MPI processes
     *
     *          The second level partitioning generates "L2 BLOCKSs".
     *          L2Blocks represent a further partitioning of a parent L1Block, and therefore does not have its own mmap call.
     *          L2Blocks are suited for distinct but cooperative computing processes, e.g. openmp threads.
     *
     *          Both L1 and L2 Blocks are specified using bliss::partition::Range, and represented by bliss::io::DataBlocks.
     *          DataBlocks optionally provides buffering at the L1 and/or L2 Block level.  See bliss::io::DataBlock for
     *          additional details.
     *
     *          Having a 2 level partitioning provides flexibility and encourages sequential file access (for L2Blocks).
     *          It is not necessary that L2Block be used, however, as L1 Blocks can be used directly for data traversal.
     *
     *        MemMap
     *          Memmap maps a file to a memory address range, and allows direct, random access to the
     *          file.  it utilizes the same mechanism as virtual memory paging, so it's efficient.
     *          Performance of mmap is good generally compared to traditional open-read-close-process
     *
     *  Simple Usage:
     *    1. instantiate file_loader - // default partitioning is based on MPI communicator
     *                                 // but user can also specify the number of partitions and partition id directly
     *    2. getNextL1Block()          // get the next L1 DataBlock for the current partition id (e.g. rank)
     *
     *    3. [get the iterators (begin(),end()) and do some work.  (L1 Block traversal)]
     *      or
     *      loop
     *        getNextL2Block(),
     *        [then get the iterators (begin(),end()) and do some work.  (L2 Block traversal)]
     *
     *    4. destroy file_loader          // close the file.
     *
     *
     *  Any modification to the L1 and L2 partitioning, specifically modifying the associated Range object, or additional partitioning,
     *          such as partitioning by record, should be done by subclasses of FileLoader.  This falls under the Advanced Usage show below.
     *
     *  Advanced Usage:  This should be relevant to subclasses only, as most of these methods are protected.
     *    1. instantiate file_loader - // default partitioning is based on number of MPI processes.
     *                                // but user can also specify the number of partitions.
     *    2. getNextL1BlockRange()       // allows iteration over partition ids, mapping of process to partition other than 1 to 1.
     *    ...                         // modify range
     *    3. getL1DataForRange(range)    // actually open the file and memmap
     *
     *
     *    4. [get the iterators (begin(),end()) and do some work.  (L1 Block traversal)]
     *      or
     *      loop
     *        getNextL2BlockRange(),
     *        ...                      // modify range
     *        getL2DataForRange(range) // get the data. and optionally buffering
     *        [then get the iterators (begin(),end()) and do some work.  (L2 Block traversal)]
     *
     *    5. unload                       // mem unmap
     *    6. destroy file_loader          // close the file.
     *
     *
     * @note:   For real data,  No Buffering is better for large files and limited memory, while
     *              buffering is better for smaller files and/or large amount of memory.
     *          FileLoader is the base class in a CRTP static polymorphism pattern.
     *
     *          Memmap require page aligned addresses.  Internally, this class has pointers to the
     *          start of the range of interest, not to the page aligned starting position.
     *          DataBlocks therefore has at the "begin" iterator the first element OF INTEREST in the mmapped range,
     *            not the first element mapped.
     *
     * @tparam  T                 type of each element read from file
     * @tparam  L1Buffering       bool indicating if L1 partition blocks should be buffered.  default to false
     * @tparam  L2Buffering       bool indicating if L2 partition blocks should be buffered.  default to true (to avoid contention between threads)
     * @tparam  L1PartitionerT    Type of the Level 1 Partitioner to generate the range of the file to load from disk
     * @tparam  L2PartitionerT    L2 partitioner, default to DemandDrivenPartitioner
     * @tparam  Derived           Derived (sub) class name.  For FileLoader, this is void.
     */
    template <typename T,
          bool L2Buffering = true,
          bool L1Buffering = false,
          typename L2PartitionerT = bliss::partition::DemandDrivenPartitioner<bliss::partition::range<size_t> >,
          typename L1PartitionerT = bliss::partition::BlockPartitioner<bliss::partition::range<size_t> >,
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
         * @typedef InputIteratorType
         * @brief   type for the raw data from the file.  element is of templated type T.
         */
        typedef T*                                  InputIteratorType;


        /**
         * @typedef L1BlockType
         * @brief   DataBlock type for blocks from L1 partitioning.  blocks may be buffered or unbuffered depending on template parameter L1Buffering.
         */
        typedef typename std::conditional<L1Buffering,
                                          bliss::io::BufferedDataBlock<  InputIteratorType, RangeType>,
                                          bliss::io::UnbufferedDataBlock<InputIteratorType, RangeType> >::type                        L1BlockType;

        /**
         * @typedef L2BlockType
         * @brief   DataBlock type for blocks from L2 partitioning.  blocks may be buffered or unbuffered depending on template parameter L2Buffering.
         */
        typedef typename std::conditional<L2Buffering,
                                          bliss::io::BufferedDataBlock<  typename L1BlockType::iterator, RangeType>,
                                          bliss::io::UnbufferedDataBlock<typename L1BlockType::iterator, RangeType> >::type     L2BlockType;

      protected:
        /**
         * @typedef InputIteratorType
         * @brief   aliased to PointerType
         */
        using PointerType = InputIteratorType;


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
         * @brief file handle for the opened file
         * @note  a FileLoader instance is tied to a file, this is conceptually a constant, so mutable.
         */
        mutable int fileHandle;

        /**
         * @brief full range of the file, in units of data type  T size.  [0, file_size / sizeof(T) ).
         * @note  as a FileLoader instance is tied to a file, this is conceptually a constant, so mutable.
         */
        mutable RangeType fileRange;  // full size of file

        //==== L1 Partitioning.
#if defined(USE_MPI)
        /// MPI communication object.  for use when we use 1 loader per MPI process.
        MPI_Comm comm;
#endif

        /**
         * @brief   id of a FileLoader instance that will be performing the memmap and reading of a subrange of the file.
         * @note    The mapping of processes/threads to loaders is flexible.  The default use case we have here is 1 to 1 mapping
         *           between MPI processes to loaders, so loaderId == mpi comm rank.
         *          Id is assigned to a process running a FileLoader, and is constant for the lifetime of the FileLoader
         */
        mutable int loaderId;

        /**
         * @brief Level 1 Partitioner, produces L1 Block ranges that will be used to memmap the file.
         * @details  defaults to BlockPartitioner
         */
        L1PartitionerT L1Partitioner;

        //==== memory mapping
        /**
         * @brief  memmapped data as raw pointer.
         * @details the pointer points to file data with position that is page aligned, as required by mmap.
         *          loaded from file using mmapRange as specification.
         */
        T* mmapData;      // mem-mapped data, page aligned.  strictly internal

        //===   L1 Block
        /// L1Block instance that wraps mmapData to provide buffering and abstraction of interface.
        L1BlockType L1Block;


        //=== L2 Blocks and partitioner
        /**
         * @brief   The number of concurrent consuming threads that are processing a file at L2 partition level
         * @details The mapping of OS threads to consuming threads is flexible.  The use case we have here is 1 to 1 mapping
         *           between openmp threads to consuming threads.
         *
         *          the threads consuming the L2 blocks should all live in the same process since they share memmapped data
         *          The number of Consuming Threads is constant for the lifetime of the FileLoader
         *
         * @note    The number of consuming threads may be different than the total number of L2Blocks, if Cyclic or DemandDriven Partitioners are used
         *           instead of BlockPartitioner
         *
         */
        const size_t nConsumingThreads;

        /**
         * @brief   The size of each L2 Blocks.  User tunable.  Also autoadjust to at least 2x record size.
         * @note    Constant for the lifetime of the File Loader.
         */
        mutable size_t L2BlockSize;

        /**
         * @brief Level 2 Partitioner, produces L2 Block ranges from the memmapped region.
         * @details  default to DemandDrivenPartitioner.  The blocks are consumed by computing threads
         */
        L2PartitionerT L2Partitioner;

        /**
         * @brief L2 Block cache to buffer access to DataBlocks during computation.  size of array is | nConsumingThreads |
         * @details The purpose of this cache is to allow slower computation to access the data, without
         *          implicitly or explicitly copying the data to the consuming threads.
         *          It also provide some thread safety as a consuming thread has its own slot in the array, and changes to that element
         *          is done by only that thread.
         */
        L2BlockType *L2Blocks;


      public:
        //================== Constructor and Destructor


#if defined(USE_MPI)
        /**
         * @brief MPI comm based constructor.  uses comm rank and size for initialization
         * @details   MPI Comm object is used to set get the current process id to use as the loader id, thus the L1 partitions, and to broadcast file size
         *
         *            The constructor first opens the file,
         *                then get the number of concurrent loaders and its loader id (comm_size, and comm_rank respectively)
         *                then configures the L1 partitioner
         *                and initializes the L2 block caching (per thread)
         *            The file is opened, not yet mmap (load).  this sequence provides flexility for caller
         *            and subclasses to customize the process.
         *
         * @note      Default behavior is to block-partition the file based on number of processors.
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
        explicit FileLoader(const MPI_Comm& _comm, const std::string &_filename, const size_t _nThreads = 1, const size_t _L2BlockSize = 1) throw (bliss::io::IOException)
            : pageSize(sysconf(_SC_PAGE_SIZE)), recordSize(1), loaded(false), fileHandle(-1), fileRange(), comm(_comm),
              L1Partitioner(), mmapData(nullptr), L1Block(),
              nConsumingThreads(_nThreads), L2BlockSize(_L2BlockSize), L2Partitioner(), L2Blocks(nullptr)
        {
          if (_filename.length() <= 0) throw std::invalid_argument("ERROR: Filename Length is less than 1");
          if (_nThreads <= 0) throw std::invalid_argument("ERROR: Number of threads is less than 1");
          if (_L2BlockSize <= 0) throw std::invalid_argument("ERROR: L2Blocksize is less than 1");

          // open the file
          fileHandle = openFile(_filename);

          // get from the communicator the number of concurrent loaders, and the id of the current loader.
          int nConcurrentLoaders = 0;
          MPI_Comm_rank(comm, &loaderId);
          MPI_Comm_size(comm, &nConcurrentLoaders);

          //====  every node can get the file size directly, so no need for MPI_Bcast.

          // get the file size.  all processes do this, since all have to open the file and read from it anyways.
          // compute the full range of the file
          fileRange = RangeType(0, getFileSize(fileHandle) / sizeof(T));   // range is in units of T

          // configure the L1 partitioner
          configL1Partitioner(L1Partitioner, fileRange, nConcurrentLoaders);

        }
#endif

        /**
         * @brief NON-MPI constructor.
         * @details   The constructor first opens the file,
         *                then configures the L1 partitioner with number of concurrent loaders (e.g. nthreads, or mpi comm_size)
         *                and its loader id (e.g. thread id, or mpi rank).
         *                and initializes the L2 block caching (per thread)
         *            The file is opened, not yet mmap (load).  this is done to support flexible sequencing of events by caller
         *            and subclasses.
         *
         * @note      Default behavior is to block-partition the file based on number of processors.
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
        FileLoader(const size_t _nConcurrentLoaders, const int _loaderId, const std::string &_filename,
                   const size_t _nThreads = 1, const size_t _L2BlockSize = 1 ) throw (bliss::io::IOException)
            : pageSize(sysconf(_SC_PAGE_SIZE)), recordSize(1), loaded(false), fileHandle(-1), fileRange(),
              loaderId(_loaderId), L1Partitioner(), mmapData(nullptr), L1Block(),
              nConsumingThreads(_nThreads), L2BlockSize(_L2BlockSize), L2Partitioner(), L2Blocks(nullptr)
        {
          if (_filename.length() <= 0) throw std::invalid_argument("ERROR: Filename Length is less than 1");
          if (_nThreads <= 0) throw std::invalid_argument("ERROR: Number of threads is less than 1");
          if (_L2BlockSize <= 0) throw std::invalid_argument("ERROR: L2Blocksize is less than 1");
          if (_loaderId < 0) throw std::invalid_argument("ERROR: Loader ID is less than 0");
          if (_nConcurrentLoaders <= _loaderId) throw std::invalid_argument("ERROR: Loader ID is greater than number of loaders");

          // open the file
          fileHandle = openFile(_filename);

          // get the file size.  all processes do this, since all have to open the file and read from it anyways.
          // compute the full range of the file
          fileRange = RangeType(0, getFileSize(fileHandle) / sizeof(T));   // range is in units of T

          // configure the L1 partitioner
          configL1Partitioner(L1Partitioner, fileRange, _nConcurrentLoaders);
        }



        /// Removed default copy constructor
        FileLoader(const FileLoader<T, L2Buffering, L1Buffering, L2PartitionerT, L1PartitionerT, Derived>& other) = delete;
        /// Removed default copy assignment operator
        FileLoader& operator=(const FileLoader<T, L2Buffering, L1Buffering, L2PartitionerT, L1PartitionerT, Derived>& other) = delete;


        /**
         * @brief move constructor.  here to ensure that the DataBlock instances are moved.
         * @param other FileLoader to move
         */
        FileLoader(FileLoader<T, L2Buffering, L1Buffering, L2PartitionerT, L1PartitionerT, Derived>&& other) :
          pageSize(other.pageSize), recordSize(other.recordSize), loaded(other.loaded),
          fileHandle(other.fileHandle), fileRange(other.fileRange),
          loaderId(other.loaderId), L1Partitioner(other.L1Partitioner), mmapData(other.mmapData),
          L1Block(std::move(other.L1Block)), nConsumingThreads(other.nConsumingThreads), L2BlockSize(other.L2BlockSize),
          L2Partitioner(other.L2Partitioner), L2Blocks(other.L2Blocks)
        {
          other.pageSize = 1;
          other.fileHandle = -1;
          other.fileRange = RangeType();
          other.mmapData = nullptr;

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
         * @brief move assignment operator.  here to ensure that the DataBlock instances are moved.
         *
         * @param other   FileLoader to move
         * @return        updated object with moved data
         */
        FileLoader& operator=(FileLoader<T, L2Buffering, L1Buffering, L2PartitionerT, L1PartitionerT, Derived>&& other) {
          if (this != &other) {
            // remove the old
            if (L2Blocks != nullptr) delete [] L2Blocks;
            // unload and close my file before copying in the other one.
            unloadL1Data();
            // close my file
            closeFile(fileHandle);

            pageSize   = other.pageSize;                         other.pageSize = 1;
            fileHandle = other.fileHandle;                       other.fileHandle = -1;
            fileRange   = other.fileRange;                         other.fileRange = RangeType();
            mmapData   = other.mmapData;                         other.mmapData = nullptr;
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

          unloadL1Data();

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
         * @brief   get the DataBlock that represents the current L1 partition block for this process.
         * @details the data can then be used directly, using the iterator accessors of DataBlock
         *          The data may be buffered depending on L1BlockType (DataBlock)
         *
         *          does NOT change the FileLoader
         *
         * @param pid           L1 Partition Id from which to get the Block
         * @return  DataBlock (buffered or unbuffered) wrapping the data mapped/read from the file.
         */
        L1BlockType& getCurrentL1Block() {
          if (!loaded) throw std::logic_error("ERROR: getting L1Block before file is loaded");

          return L1Block;
        }

        /**
         * @brief   get the DataBlock that represents the current L2 partition block for this thread.
         * @details the data can then be used directly, using the iterator accessors of DataBlock
         *          The data may be buffered depending on L2BlockType (DataBlock)

         *          does NOT change the FileLoader
         *
         * @param tid           L2 Partition Id from which to get the Block (thread id)
         * @return  DataBlock (buffered or unbuffered) wrapping the data mapped/read from the file.
         */
        L2BlockType& getCurrentL2Block(const size_t &tid) {
          if (!loaded) throw std::logic_error("ERROR: getting L2Block before file is loaded");

          return L2Blocks[tid];
        }

        /**
         * @brief   get the next L1 DataBlock from the L1 partitioner for this process.
         * @details the data can then be used directly, using the iterator accessors of DataBlock
         *          The data may be buffered depending on L1BlockType (DataBlock)
         *
         *          This method is a convenience method that gets the next L1BlockRange, and then load the specified data from disk.
         *          it does not allow modification of the range before data loading
         *
         *          this method defaults to the loaderId specified at FileLoader construction.
         *
         * @return  DataBlock (buffered or unbuffered) wrapping the data mapped/read from the file.
         */
        L1BlockType& getNextL1Block() {
          // first get the Range
          RangeType r = getNextL1BlockRange<Derived>();

          return loadL1DataForRange(r);
        }



        /**
         * @brief   get the DataBlock that represents the current L2 partition block for this thread.
         * @details the data can then be used directly, using the iterator accessors of DataBlock
         *          The data may be buffered depending on L2BlockType (DataBlock)
         *
         *          This method is a convenience method that gets the next L2BlockRange, and then gets the DataBlock with that range
         *          it does not allow modification of the range before data loading
         *
         * @param tid           L2 Partition Id from which to get the Block (thread id)
         * @return  DataBlock (buffered or unbuffered) wrapping the data from the parent L1 block.
         */
        L2BlockType& getNextL2Block(const size_t &tid) {

          RangeType r = getNextL2BlockRange<Derived>(tid);

          return getL2DataForRange(tid, r);
        }



        /**
         * @brief reset the L1 partitioner so it can be reused.
         */
        void resetL1Partitioner() {
          L1Partitioner.reset();
        }


        /**
         * @brief Reset the L2 partitioner so we can reuse it.
         */
        void resetL2Partitioner() {
          L2Partitioner.reset();
        }


      protected:
        /**
         * @brief get the next available range for the L1 partition block, given the L1 partition id
         * @details   This method uses CRTP to allow calling Derived class' implementation, which can further refine the partitioned range.
         * @note      If the default BlockPartitioner were used, then there is only 1 NextL1Block.  else there may be more than 1.
         *
         * @param pid   The L1 partition id.
         * @return      The range of the next L1 partition block.
         */
        template<typename D = Derived>
        typename std::enable_if<std::is_void<D>::value, RangeType>::type getNextL1BlockRange(const size_t pid) {
          // just FileLoader calls this.
          return L1Partitioner.getNext(pid);
        }
        /**
         * @brief get the next available range for the L1 partition block, given the L1 partition id
         * @details   This method uses CRTP to allow calling Derived class' implementation, which can further refine the partitioned range.
         * @note      If the default BlockPartitioner were used, then there is only 1 NextL1Block.  else there may be more than 1.
         *
         * @param pid   The L1 partition id.
         * @return      The range of the next L1 partition block.
         */
        template<typename D = Derived>
        typename std::enable_if<!std::is_void<D>::value, RangeType>::type getNextL1BlockRange(const size_t pid) {
          // use the derived one.
          return static_cast<Derived*>(this)->getNextL1BlockRangeImpl(pid);
        }
        /**
         * @brief get the next available range for the L1 partition block, given the L1 partition id
         * @details   This method uses CRTP to allow calling Derived class' implementation, which can further refine the partitioned range.
         *            The loaderId is set at initialization, possibly as MPI comm rank (if the comm version of the constructor were called)
         * @note      If the default BlockPartitioner were used, then there is only 1 NextL1Block.  else there may be more than 1.
         *
         * @return      The range of the next L1 partition block.
         */
        template<typename D = Derived>
        typename std::enable_if<std::is_void<D>::value, RangeType>::type getNextL1BlockRange() {
          // just FileLoader calls this.
          return L1Partitioner.getNext(loaderId);
        }
        /**
         * @brief get the next available range for the L1 partition block, using the loaderID
         * @details   This method uses CRTP to allow calling Derived class' implementation, which can further refine the partitioned range.
         *            The loaderId is set at initialization, possibly as MPI comm rank.  (if the comm version of the constructor were called)
         * @note      If the default BlockPartitioner were used, then there is only 1 NextL1Block.  else there may be more than 1.
         *
         * @return      The range of the next L1 partition block.
         */
        template<typename D = Derived>
        typename std::enable_if<!std::is_void<D>::value, RangeType>::type getNextL1BlockRange() {
          // use the derived one.
          return static_cast<Derived*>(this)->getNextL1BlockRangeImpl(loaderId);
        }

        /**
         * @brief Loads the file by memory maps the specified range for reading,
         * @details  Optionally buffers the data into memory
         *           The input range SHOULD come from the call to getNextL1BlockRange
         *
         * @param r   range that will be mapped and loaded into memory (with buffering or without)
         * @return    loaded L1 DataBlock for the specified Range
         */
        L1BlockType& loadL1DataForRange(const RangeType &range) throw (bliss::io::IOException)
        {
          // clean up any previous runs.
          unloadL1Data();

          // don't do anything if range is empty.
          if (range.size() == 0) return L1Block;

          // make sure the mmapRange is within the file range, and page align it.
          RangeType mmapRange = RangeType::intersect(range, fileRange);
          size_t block_start = RangeType::align_to_page(mmapRange.start, pageSize);

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

          // map the region of file to memory
          mmapData = map(mmapRange);

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
          configL2Partitioner(L2Partitioner, L2Blocks, mmapRange, nConsumingThreads, L2BlockSize);

          return L1Block;
        }

        /**
         * @brief   unloads the memmapped region of a file from memory
         * @details if buffering, clear the data
         *          if not buffering, then unmemmap the file region from memory
         */
        void unloadL1Data()
        {
          // if the partition was not buffered,
          if (!L1Buffering)
          {
            // and data has been memmapped,
            if (mmapData != nullptr && mmapData != MAP_FAILED)
              // then we unmap
              this->unmap(mmapData, L1Block.getRange());
          } // else the partition was preloading, so already unmapped

          // clean up L1Block
          L1Block.clear();

          // reset the data pointer, loaded flag, and the variable storing the memmapped range.
          mmapData = nullptr;
          loaded = false;

          resetL2Partitioner();
        }



        /**
         * @brief get the next L2 Partition Block given a thread id.
         * @details  for derived File Loader classes, their implementation methods are called.
         *            This call needs to properly support concurrent range computation.
         * @note  not providing a no-arg default impl because the threading is setup before these calls and after construction.,
         *        and we do not want to get the id in this method from the threading library.
         * @param tid   Thread Id.  L2 Partition is between threads within a group
         * @return  The range of the next L2 Partition Block.
         */
        template <typename D = Derived>
        typename std::enable_if<std::is_void<D>::value, RangeType>::type getNextL2BlockRange(const size_t tid) {
          if (!loaded) throw std::logic_error("ERROR: getting L2Block range before file is loaded");
          return L2Partitioner.getNext(tid);
        }

        /**
         * @brief get the next L2 Partition Block given a thread id.
         * @details  for derived File Loader classes, their implementation methods are called.
         *            This call needs to properly support concurrent range computation.
         * @note  not providing a no-arg default impl because the threading is setup before these calls and after construction.,
         *        and we do not want to get the id in this method from the threading library.
         * @param tid   Thread Id.  L2 Partition is between threads within a group
         * @return  The range of the next L2 Partition Block.
         */
        template <typename D = Derived>
        typename std::enable_if<!std::is_void<D>::value, RangeType>::type getNextL2BlockRange(const size_t tid) {
          if (!loaded) throw std::logic_error("ERROR: getting L2Block range before file is loaded");
          return static_cast<Derived*>(this)->getNextL2BlockRangeImpl(tid);
        }


        /**
         * @brief   get the DataBlock from L2 Partitioner for the given thread Id.
         * @note    Depending on data type of L2Block, the data may be buffered (copied to) memory.
         * @param tid         id of thread calling this function.
         * @param L2BlockRange   The range of the L2Block to retrieve. computed by getNextL2BlockRange.
         * @return    reference to loaded L2Block for the specified thread id and range.
         */
        L2BlockType& getL2DataForRange(const size_t tid, const RangeType &L2BlockRange) {
          if (!loaded) throw std::logic_error("ERROR: getting L2Block data before file is loaded");

          // make sure the range is within the mmaped region.
          RangeType r = RangeType::intersect(L2BlockRange, L1Block.getRange());

          // now compute the start and end iterators of this range.  use std::advance in case the iterator is not a random access iterator.
          auto s = L1Block.begin();
          std::advance(s, (r.start - L1Block.getRange().start));
          auto e = s;
          std::advance(e, r.size());

          // get (optionally buffer) the data into a L2Block, and cache it for repeated access
          L2Blocks[tid].assign(s, e, r);

          return L2Blocks[tid];
        }



        // these methods are used by subclasses, possibly overridden by subclasses.
        /**
         * @brief   configure the L1 Partitioning
         * @details This function does the following:
         *            construct the fileRange from file size
         *            configure the L1Partitioner.
         *
         *          file_size is passed in because for the MPI version of c'tor, the file_size was broadcast first
         * @param[in/out] partitioner L1 Partitioner to configure
         * @param[in] fileRange   portion of file to be partitioned by partitioner (0 to file size by defautl).
         * @param[in] nPartitions  number of partitions to create.
         * @param[in] chunkSize   the size of each partitioned chunk (== L1Block size).  default to 0, which means computed.
         */
        void configL1Partitioner(L1PartitionerT& partitioner, const RangeType& fileRange, const size_t& nPartitions, const size_t& chunkSize = 0) {

          // Set up the partitioner to partition the full file range with nPartitions.  Overlap is 0.
          partitioner.configure(fileRange, nPartitions, chunkSize);
        }




        /**
         * @brief configure the L2 Partitioning
         * @details: This function does the following
         *            get the approximate record size and update the chunkSize for the L2Partitioner
         *            configure L2 partitioner for the memmapped range, number of threads, and L2Block Size.
         *            initialize the L2Block cache, one per thread.
         * @param[in/out] partitioner L2 Partitioner to configure
         * @param[in/out] cache   L2 Block cache, 1 per thread.
         * @param[in] mmapRange   memory mapped region's range, to be partitioned by partitioner.
         * @param[in] nThreads    number of threads to create.
         * @param[in/out] chunkSize   the size of each partitioned chunk (== L1Block size).
         */
        void configL2Partitioner(L2PartitionerT& partitioner, L2BlockType* &cache, const RangeType& mmapRange, const size_t& nThreads, size_t& chunkSize) {
          //=====  adjust the L2BlockSize adaptively.
          // look through 3 records to see the approximate size of records.
          // update the chunkSize  to at least 2x record size, (passes back to called).
          chunkSize = std::max(chunkSize, 2 * getRecordSize<Derived>(10));

          // then configure the partitinoer to use the right number of threads.
          partitioner.configure(L1Block.getRange(), nThreads, chunkSize);

          // also configure the cache.
          if (cache != nullptr) delete [] cache;
          cache = new L2BlockType[nThreads];
        }

      public:
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
        typename std::enable_if<std::is_void<D>::value, size_t>::type getRecordSize(const int count = 10) {
          // this method is enabled if Derived is void, which indicates self.
          return 1;
        }

      public:

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
        typename std::enable_if<!std::is_void<D>::value, size_t>::type getRecordSize(const int count = 10) {
          return static_cast<Derived*>(this)->getRecordSizeImpl(count);
        }

        /// get number of estimated kmers, given k.
        template<typename D = Derived>
        typename std::enable_if<std::is_void<D>::value, size_t>::type getKmerCountEstimate(const int k) {
          return this->getFileRange().size() - k + 1;
        }

        /// get number of estimated kmers, given k.
        template<typename D = Derived>
        typename std::enable_if<!std::is_void<D>::value, size_t>::type getKmerCountEstimate(const int k) {
          return static_cast<Derived*>(this)->getKmerCountEstimateImpl(k);
        }

      private:
        // these methods are not meant to be overridden by subclasses.

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


      protected:

        /**
         * @brief     map the specified portion of the file to memory.
         * @param r   range specifying the portion of the file to map
         * @return    memory address (pointer) to where the data is mapped.
         */
        PointerType map(const RangeType &r) throw (IOException) {

          /// memory map.  requires that the starting position is block aligned.
          size_t block_start = RangeType::align_to_page(r, pageSize);

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
        void unmap(PointerType &d, const RangeType &r) {

          munmap(d, (r.end - RangeType::align_to_page(r, pageSize)) * sizeof(T));
        }

    };

  } /* namespace functional */
} /* namespace bliss */
#endif /* FILE_LOADER_HPP_ */
