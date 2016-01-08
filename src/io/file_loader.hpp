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

#include "bliss-config.hpp"

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

#include <sys/sysinfo.h>  // for meminfo
#include <type_traits>

#include "partition/range.hpp"
#include "partition/partitioner.hpp"
#include "io/data_block.hpp"
#include "io/io_exception.hpp"
#include "utils/logging.h"
#include "common/sequence.hpp"
#include <mxx/comm.hpp> // for mxx::comm


namespace bliss
{
namespace io
{
  struct BaseFile {

  };

  // TODO: simplify  possibly by separating L1 and L2 partitioning.


  /**
   * @brief  Partitions input data into sequences, e.g. blocks or reads.
   * @details  Iterator template parameter refers to the iterator type that is currently being traversed.
   * 		 this applies to the sequence object, and the following methods:
   * 		 get_next_record, and find_first_record.
   * 		 
   * 		 init method can potentially use a different iterator type. this is useful when
   * 		   state information is extracted from perhaps a containing block (e.g. L1Block),
   * 		   then used with a sub block (e.g. L2Block).  
   * 		 	
   * @note   there are 2 (non-interchangeable) scenarios when we may need to change the iterator type.
   *    1. L1 at process level, L2 at thread level.  SeqParser state is calculated at L1.
   *    2. L1 at rank 0 of a sub communicator (processes), L2 is broadcast from rank0 of the rest of processes.  SeqParser state is calculated at L1.
   *
   *    to handle the first case: approaches are
   *    1a. to allow a parser to operate on different iterator types when invoking find_first_record, increment.
   *          additional template parameter for these functions.  SequenceType will need to change.
   *          API in SequenceIterator, BaseFileParser, FASTQParser, FASTAParser all need to change
   *    1b. convert the seq parser from L1 to L2, but copy the state of L1.
   *          can be shared between threads.  additional api for convert only
   *
   *    to handle the second case.
   *    2a. convert the seq parser from L1 to L2, communicate (and copy) the state of L1.
   *          additional api for convert and broadcast. minimal api change. requires L1 to be block (or cyclic with synchronization).
   *    2b. create seq parser for L2, and initialize using all processes.
   *          no api change, but works only when L1 and L2 are block partitioned (or cyclic with synchronization).
   *          potentially less communication, and Parser does not need to have a broadcast semantic.
   *
   *    1a is not general, potentially complicated, and a source of potential errors.
   *    1b and 2a share common apis.
   *    2b requires application to specifically call init at L2.  Also, already initialized L1 during data load.
   *
   *    choose 2b and 1b.  Sequence Iterator does not call init_for_iterator, mostly for flexibility and compatibility with case 1.
   *         
   * @tparam Iterator   The underlying iterator to be traversed to generate a Sequence
   *
   */
  template <typename Iterator = unsigned char* >
   class BaseFileParser {

      // let another BaseFileParser with some other iterator type be a friend so we can convert.
      template <typename Iterator2>
      friend class BaseFileParser;

     public:
//      using SequenceIdType = SeqIdType;
      using SequenceType = typename ::bliss::common::Sequence<Iterator>;
      using SequenceIdType = typename SequenceType::IdType;

      /// static constant for end of line.  note this is same for unicode as well
      static constexpr unsigned char eol = '\n';
      /// static constant for carriage return.  note this is same for unicode as well
      static constexpr unsigned char cr = '\r';

     protected:
       /**
        * @typedef RangeType
        * @brief   range object types
        */
      using RangeType = bliss::partition::range<size_t>;

       /**
        * @brief  search for first non-EOL character in a iterator, returns the stopping position as an iterator.  also records offset.
        * @details       iter can point to the previous EOL, or a nonEOL character.
        * @param[in/out] iter    iterator to advance
        * @param[in]     end     position to stop the traversal
        * @param[in/out] offset  the global offset of the sequence record within the file.  used as record id.
        * @return        iterator at the new position, where the Non EOL char is found, or end.
        */
       template <typename IT = Iterator,
           typename = typename ::std::enable_if<(::std::is_same<typename ::std::iterator_traits<IT>::value_type, char>::value ||
                                                ::std::is_same<typename ::std::iterator_traits<IT>::value_type, unsigned char>::value)
                                               >::type >
       inline IT findNonEOL(IT& iter, const IT& end, size_t &offset) const {
         while ((iter != end) && ((*iter == eol) || (*iter == cr))) {
           ++iter;
           ++offset;
         }
         return iter;
       }

       /**
        * @brief  search for first EOL character in a iterator, returns the stopping position as an iterator.  also records offset.
        * @details       iter can point to a nonEOL char, or an EOL char.
        * @param[in/out] iter    iterator to advance
        * @param[in]     end     position to stop the traversal
        * @param[in/out] offset  the global offset of the sequence record within the file.  used as record id.
        * @return        iterator at the new position, where the EOL char is found, or end.
        */
       template <typename IT = Iterator,
           typename = typename ::std::enable_if<(::std::is_same<typename ::std::iterator_traits<IT>::value_type, char>::value ||
                                                ::std::is_same<typename ::std::iterator_traits<IT>::value_type, unsigned char>::value)
                                               >::type >
       inline IT findEOL(IT& iter, const IT& end, size_t &offset) const {
         while ((iter != end) && ((*iter != eol) && (*iter != cr) ) ) {
           ++iter;
           ++offset;
         }
         return iter;
       }


       /**
        * @brief constructs an IOException object with the relevant debug data.
        * @param errType     string indicating the source of the error.  user choice
        * @param start       iterator pointing to beginning of the data in question
        * @param end         iterator pointing to end of the data in question
        * @param startOffset offset for the beginning of the data in question
        * @param endOffset   offset for the end of the data in question
        */
       void handleError(const std::string& errType, const Iterator &start, const Iterator &end, const size_t& startOffset, const size_t& endOffset) throw (bliss::io::IOException) {
         std::stringstream ss;
         ss << "ERROR: did not find "<< errType << " in " << startOffset << " to " << endOffset << std::endl;
         ss << "  offending string is \"" << std::endl;
         std::ostream_iterator<typename std::iterator_traits<Iterator>::value_type> oit(ss);
         std::copy(start, end, oit);
         ss << "\".";
         throw bliss::io::IOException(ss.str());
       }

       /**
        * @brief print out an Warning, for example when there is malformed partial data.
        * @param errType     string indicating the source of the error.  user choice
        * @param start       iterator pointing to beginning of the data in question
        * @param end         iterator pointing to end of the data in question
        * @param startOffset offset for the beginning of the data in question
        * @param endOffset   offset for the end of the data in question
        */
       void handleWarning(const std::string& errType, const Iterator &start, const Iterator &end, const size_t& startOffset, const size_t& endOffset) {
         BL_WARNING("WARNING: " << "did not find "<< errType << " in " << startOffset << " to " << endOffset);
       }


     public:

       /// default constructor.
       BaseFileParser() : offset(0) {};

       /// default destructor
       virtual ~BaseFileParser() {};


       /// converting constructor.
       template <typename Iterator2>
       BaseFileParser(BaseFileParser<Iterator2> const & other) : offset(other.offset) {}
       /// converting assignment operator that can transform the base iterator type.
       template <typename Iterator2>
       BaseFileParser<Iterator>& operator=(BaseFileParser<Iterator2> const & other) {
         offset = other.offset;
         return *this;
       }


       /**
        * @brief given a block/range, find the starting point of the first sequence object (here, just the actual start)
        * @note  not virtual for performance reason.
        *         used by getNextL1BlockRange to find exact range for an L1 Partition.
        *
        *         Also can be used to do other initializations.
        *
        * @param _data
        * @param parentRange
        * @param inMemRange
        * @param searchRange
        * @return
        */

       std::size_t find_first_record(const Iterator &_data, const RangeType &parentRange, const RangeType &inMemRange, const RangeType &searchRange)
       {
         //== range checking
         if(!parentRange.contains(inMemRange)) throw std::invalid_argument("ERROR: Parent Range does not contain inMemRange");

         // no op, other than to use intersect to make sure that the start is valid and within parent, and in memry ranges.
         return RangeType::intersect(searchRange, inMemRange).start;
       }


       // TODO: store offset for first record...

       virtual void reset() {}

       virtual std::size_t init_parser(const Iterator &_data, const RangeType &parentRange, const RangeType &inMemRange, const RangeType &searchRange) {
         return find_first_record(_data, parentRange, inMemRange, searchRange);
       }
#ifdef USE_MPI
       /// initializes the parser.  only useful for FASTA parser for now.  Assumes searchRange do NOT overlap between processes.
       virtual std::size_t init_parser(const Iterator &_data, const RangeType &parentRange, const RangeType &inMemRange, const RangeType &searchRange, const mxx::comm& comm)
       {
         return find_first_record(_data, parentRange, inMemRange, searchRange);
       };
#endif



       /**
        * @brief increment to get next sequence object
        * @note   not virtual because SequenceType is different for each subclass, so the signatures are different.  We rely on compiler to enforce the correct one.
        *         this is in general not an issue as we do not need polymorphic Sequence Parsers - they are specific to file types.
        *         this is probably better for performance anyways.
        *
        *         Used by FileLoader's getRecordSize to compute an average record size.
        *
        * @param iter           start of a sequence object.  this is the beginning of a record, not just DNA sequence
        * @param end            end of a sequence object.  this is the end of a record, not just DNA sequence
        * @param offset         offset in the file for the start of the record.
        * @param seq_offset     index in the local (shared) vector of sequence breaks, if there is one (e.g. FASTA).  used by FASTA to lookup the nearest complete sequence record info.
        * @param output         updated sequence object.
        * @return               next seq id offset.
        */
       SequenceType get_next_record(Iterator & iter, const Iterator & end, size_t & offset)
       {
         Iterator orig_iter = iter;
         size_t orig_offset = offset;

         size_t dist = std::distance(iter, end);

         offset += dist;
         iter = end;

         return SequenceType(SequenceIdType(orig_offset), dist, 0, orig_iter, end);
       }

   };
   /// template class' static variable definition (declared and initialized in class)
  template <typename Iterator>
   constexpr unsigned char bliss::io::BaseFileParser<Iterator >::eol;
  template <typename Iterator>
   constexpr unsigned char bliss::io::BaseFileParser<Iterator>::cr;

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
     *          Performance of mmap is good generally compared to traditional file open-read-close-process
     *          NOTE: reading past the mmapped length can cause SEGV
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
     * Alternatively, modify the range in BaseFileParser.  See fastq_loader.hpp for example.
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
     * @tparam  Overlap           Size of overlap between L1 or L2 blocks. in units of T.  this affects mmap.
     * @tparam  Parser            Sequence Parser for generating the the sequence objects.  NOTE: template template parameter, template being an Iterator with value type T.
     * @tparam  L1Buffering       bool indicating if L1 partition blocks should be buffered.  default to false
     * @tparam  L2Buffering       bool indicating if L2 partition blocks should be buffered.  default to true (to avoid contention between threads)
     * @tparam  L1PartitionerT    Type of the Level 1 Partitioner to generate the range of the file to load from disk
     * @tparam  L2PartitionerT    L2 partitioner, default to DemandDrivenPartitioner
     */
  // TODO: simplify and cleanup.

    template <typename T, size_t Overlap = 0, template <typename > class Parser = BaseFileParser,
          bool L2Buffering = true,
          bool L1Buffering = false,
          typename L2PartitionerT = bliss::partition::DemandDrivenPartitioner<bliss::partition::range<size_t> >,
          typename L1PartitionerT = bliss::partition::BlockPartitioner<bliss::partition::range<size_t> > >
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
                                          bliss::io::UnbufferedDataBlock<typename L1BlockType::iterator, RangeType> >::type           L2BlockType;

        /**
         * return overlap size
         */
        static constexpr size_t get_overlap_size() { return Overlap; }

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

        /**
         * @brief number of sequence characters in a sequence record
         */
        mutable size_t seqSizeInRecord;

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
         * @brief   The size of each L1 Blocks.  User tunable.  autocomputed when using BlockPartitioner.
         * @note    Constant for the lifetime of the File Loader.
         */
        mutable size_t L1BlockSize;

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


        Parser<typename L1BlockType::iterator> L1parser;

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
         *            Range is in units of T, overlap is INCLUDED in range.end
         *
         * @param _filename       input file name
         * @param _comm           MPI Communicator (defined only if CMake is configured with MPI).
         *                        default is MPI_COMM_WORLD, which means all nodes participate.
         * @param _nThreads       number of threads to consume L2 blocks
         * @param _L2BlockSize    size of each L2Block.
         * @param _L1BlockSize    size of each L1Block.  can default to 0 (PAGE_SIZE used instead)
         *
         * TODO: test on remotely mounted file system.
         *
         */
        explicit FileLoader(const std::string &_filename, const MPI_Comm& _comm, const size_t _nThreads = 1,
                            const size_t _L2BlockSize = 0, const size_t _L1BlockSize = 0) throw (bliss::io::IOException)
            : pageSize(sysconf(_SC_PAGE_SIZE)), recordSize(0), seqSizeInRecord(0),
              loaded(false), fileHandle(-1), fileRange(), comm(_comm),
              L1Partitioner(), mmapData(nullptr), L1Block(),
              nConsumingThreads(_nThreads), L2Partitioner(), L2Blocks(nullptr)
        {
          if (_filename.length() <= 0) throw std::invalid_argument("ERROR: Filename Length is less than 1");
          if (_nThreads <= 0) throw std::invalid_argument("ERROR: Number of threads is less than 1");

          L1BlockSize = (_L1BlockSize == 0) ? pageSize : _L1BlockSize;
          L2BlockSize = (_L2BlockSize == 0) ? L1BlockSize / _nThreads :
                        std::min(_L2BlockSize, (L1BlockSize / _nThreads));

          if (L1BlockSize == 0) throw std::invalid_argument("ERROR: L1Blocksize is 0");
          if (L2BlockSize == 0) throw std::invalid_argument("ERROR: L2Blocksize is 0");

          // open the file
          fileHandle = openFile(_filename);

          // get the file size.  all processes do this, since all have to open the file and read from it anyways.
          // compute the full range of the file
          fileRange = RangeType(0, getFileSize(fileHandle) / sizeof(T));   // range is in units of T

          // open the first part of the file, and estimate the record sizes.
          std::tie(recordSize, seqSizeInRecord) = getRecordSize(10);


          // get from the communicator the number of concurrent loaders, and the id of the current loader.
          int nConcurrentLoaders = 0;
          MPI_Comm_rank(comm, &loaderId);
          MPI_Comm_size(comm, &nConcurrentLoaders);


          // configure the L1 partitioner
          configL1Partitioner(L1Partitioner, fileRange, nConcurrentLoaders, L1BlockSize);

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
         *            Range is in units of T, overlap is INCLUDED in range.end
         *
         * @param _filename       input file name
         * @param _nConcurrentLoaders  number of concurrent file loaders operating on the file.
         * @param _loaderId       id of each file loader instance.
         * @param _nThreads       number of threads to consume L2 blocks
         * @param _L2BlockSize    size of each L2Block.
         * @param _L1BlockSize    size of each L1Block.  can default to 0 (PAGE_SIZE used instead)
         *
         * TODO: test on remotely mounted file system.
         *
         */
        FileLoader(const std::string &_filename, const size_t _nConcurrentLoaders, const size_t _loaderId, const size_t _nThreads = 1,
                   const size_t _L2BlockSize = 0, const size_t _L1BlockSize = 0 ) throw (bliss::io::IOException)
            : pageSize(sysconf(_SC_PAGE_SIZE)), recordSize(0), seqSizeInRecord(0),
              loaded(false), fileHandle(-1), fileRange(),
#ifdef USE_MPI
              comm(MPI_COMM_SELF),   // still needs to initialize this.
#endif
              loaderId(_loaderId), L1Partitioner(), mmapData(nullptr), L1Block(),
              nConsumingThreads(_nThreads), L2Partitioner(), L2Blocks(nullptr)
        {

          if (_filename.length() <= 0) throw std::invalid_argument("ERROR: Filename Length is less than 1");
          if (_nThreads <= 0) throw std::invalid_argument("ERROR: Number of threads is less than 1");
          //if (_loaderId < 0) throw std::invalid_argument("ERROR: Loader ID is less than 0");
          if (_nConcurrentLoaders <= _loaderId) throw std::invalid_argument("ERROR: Loader ID is greater than number of loaders");

          L1BlockSize = (_L1BlockSize == 0) ? pageSize : _L1BlockSize;
          L2BlockSize = (_L2BlockSize == 0) ? L1BlockSize / _nThreads :
                        std::min(_L2BlockSize, (L1BlockSize / _nThreads));
          if (L1BlockSize == 0) throw std::invalid_argument("ERROR: L1Blocksize is 0");
          if (L2BlockSize == 0) throw std::invalid_argument("ERROR: L2Blocksize is 0");


          // open the file
          fileHandle = openFile(_filename);

          // get the file size.  all processes do this, since all have to open the file and read from it anyways.
          // compute the full range of the file
          fileRange = RangeType(0, getFileSize(fileHandle) / sizeof(T));   // range is in units of T

          // open the first part of the file, and estimate the record sizes.
          std::tie(recordSize, seqSizeInRecord) = getRecordSize(10);


          // configure the L1 partitioner
          configL1Partitioner(L1Partitioner, fileRange, _nConcurrentLoaders, L1BlockSize);
        }



        /// Removed default copy constructor
        FileLoader(const FileLoader & other) = delete;
        /// Removed default copy assignment operator
        FileLoader& operator=(const FileLoader & other) = delete;


        /**
         * @brief move constructor.  here to ensure that the DataBlock instances are moved.
         * @param other FileLoader to move
         */
        FileLoader(FileLoader && other) :
          pageSize(other.pageSize), recordSize(other.recordSize), seqSizeInRecord(other.seqSizeInRecord), loaded(other.loaded),
          fileHandle(other.fileHandle), fileRange(other.fileRange),
          loaderId(other.loaderId), L1BlockSize(other.L1BlockSize), L1Partitioner(other.L1Partitioner), mmapData(other.mmapData),
          L1Block(std::move(other.L1Block)), nConsumingThreads(other.nConsumingThreads), L2BlockSize(other.L2BlockSize),
          L2Partitioner(other.L2Partitioner), L2Blocks(other.L2Blocks)
        {
          other.pageSize = 1;
          other.fileHandle = -1;
          other.fileRange = RangeType();

          other.mmapData = nullptr;

          other.loaded = false;
          other.nConsumingThreads = 1;
          other.L1BlockSize = 1;
          other.L2BlockSize = 1;
          other.loaderId = 0;
          other.L2Blocks = nullptr;
          other.L1Partitioner = L1PartitionerT();
          other.L2Partitioner = L2PartitionerT();
          other.recordSize = 0;
          other.seqSizeInRecord = 0;

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
        FileLoader& operator=(FileLoader && other) {
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
            L1BlockSize   = other.L1BlockSize;                         other.L1BlockSize = 1;
            L2BlockSize   = other.L2BlockSize;                         other.L2BlockSize = 1;
            loaderId        = other.loaderId;                              other.loaderId = 0;
  #if defined(USE_MPI)
            comm        = other.comm;                              other.comm = MPI_COMM_NULL;
  #endif
            L2Blocks  = other.L2Blocks;                        other.L2Blocks = nullptr;
            L1Partitioner = other.L1Partitioner;                       other.L1Partitioner = bliss::partition::BlockPartitioner<RangeType>();
            L2Partitioner = other.L2Partitioner;             other.L2Partitioner = L2PartitionerT();
            recordSize  = other.recordSize;                        other.recordSize = 0;
            seqSizeInRecord  = other.seqSizeInRecord;                        other.seqSizeInRecord = 0;
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

        Parser<typename L1BlockType::iterator> getSeqParser() {
          return L1parser;
        }

        /**
         * @brief     Get the L2 partition Block Size.
         * @return    L2 block size
         */
        const size_t& getL2BlockSize() const {
          return L2BlockSize;
        }

        /**
         * @brief     Get the L2 partition Block Size.
         * @return    L2 block size
         */
        const size_t& getL1BlockSize() const {
          return L1BlockSize;
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
          RangeType r = getNextL1BlockRange();

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

          RangeType r = getNextL2BlockRange(tid);

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
         * @note      If the default BlockPartitioner were used, then there is only 1 NextL1Block.  else there may be more than 1.
         * @details   This method is the CTRP implementation called by the base class' getNextL1BlockRange method.
         *            The method modifies the base range to align the boundary to record boundaries at both ends
         *
         *            This is accomplished by memmapping a Range and performing find_first_record to find the
         *              start, and doing the same on the next L1Block to find the end.
         *            Start and End may be at the same position.  since we use recordSize, there is some hope that we find a valid starting point.
         * @param pid   The L1 partition id.
         * @return      The range of the next L1 partition block.
         */
        RangeType getNextL1BlockRange(const size_t pid) {
        	// NOTE: that we CANNOT use mpi to communicate since we are not always using block partitioner.
        	// NOTE: when searching, need to extend by 1 record size - so that the search can succeed.  to be safe, since recordSize is an estimate, do 2

          // get basic range
          RangeType hint = this->L1Partitioner.getNext(pid);


          // output data structure
          RangeType output(hint);

          // TODO: SEARCH 1 and communicate, then record size does not matter.

          // find starting and ending positions.  do not search in overlap region.
          if (hint.size() > 0) {

            // extend by 2 record
            RangeType start_search_range(hint.start, hint.end + 2 * recordSize);
            start_search_range.intersect(this->fileRange);
            RangeType end_search_range(hint.end, hint.end + 2 * recordSize);
            end_search_range.intersect(this->fileRange);


            // get the combined ranges
            RangeType loadRange = start_search_range;

            // memmap the content
            auto block_start = RangeType::align_to_page(loadRange, this->pageSize);
            auto mappedData = this->map(loadRange);  // this is page aligned.
            auto searchData = mappedData + (loadRange.start - block_start);

            Parser<decltype(searchData)> parser;  // local parser, since iterator type may not be the same.

            // search for new start and end using find_first_record
            try {
              output.start = parser.init_parser(searchData, this->fileRange, loadRange, start_search_range);
              output.end = parser.find_first_record(searchData, this->fileRange, loadRange, end_search_range);   // TODO: can't call this again

              //std::cout << "start in " << hint << " end in " << end_search_range << " final " << output << std::endl;

            } catch (IOException& ex) {
              // either start or end are not found so return an empty range.

              // TODO: need to handle this scenario better - should keep search until end.
              BL_WARNINGF("%s\n", ex.what());

              BL_WARNINGF("curr range: partition hint %lu-%lu, next %lu-%lu, file_range %lu-%lu\n",
                     hint.start, hint.end, end_search_range.start, end_search_range.end, this->fileRange.start, this->fileRange.end);
              BL_WARNINGF("got an exception search for partition:  %s \n", ex.what());

              output.start = hint.end;
              output.end = hint.end;
            }

            // clean up and unmap
            this->unmap(mappedData, loadRange);
          }
          return output;
        }

        /**
         * @brief get the next available range for the L1 partition block, given the L1 partition id
         * @details   This method uses CRTP to allow calling Derived class' implementation, which can further refine the partitioned range.
         *            The loaderId is set at initialization, possibly as MPI comm rank (if the comm version of the constructor were called)
         * @note      If the default BlockPartitioner were used, then there is only 1 NextL1Block.  else there may be more than 1.
         *
         * @return      The range of the next L1 partition block.
         */
        RangeType getNextL1BlockRange() {
          // just FileLoader calls this.
          return getNextL1BlockRange(loaderId);
        }

        /**
         * @brief Loads the file by memory maps the specified range for reading,
         * @details  Optionally buffers the data into memory
         *           The input range SHOULD come from the call to getNextL1BlockRange
         *
         *
         * @note     OVERLAP AWARE.  generated L1 block contains the overlap at the end.
         * @param r   range that will be mapped and loaded into memory (with buffering or without)
         * @return    loaded L1 DataBlock for the specified Range
         */
        L1BlockType& loadL1DataForRange(const RangeType &range) throw (bliss::io::IOException)
        {
          // clean up any previous runs.
          unloadL1Data();

          // don't do anything if range is empty.
          if (range.size() == 0) {
            L1Block.assign(nullptr, nullptr, range);
            return L1Block;
          }

          // make sure the mmapRange is within the file range, and page align it.
          RangeType blockRange = RangeType::intersect(range, fileRange);

          // search for overlap. naive search since we are not MPI aware here.
          RangeType overlappedRange = blockRange;
          overlappedRange.intersect(fileRange);
          size_t mmap_start = RangeType::align_to_page(overlappedRange.start, pageSize);

          // TODO: Handle Overlap better.  overlap should belong to the current proc, not previous.
          if (Overlap > 0) {
            // search for overlap starts at end of the requested block
            RangeType olr(overlappedRange.end, fileRange.end);
            size_t mmap_olr_start = RangeType::align_to_page(olr.start, pageSize);

            // map
            mmapData = map(olr);

            // then advance Overlap characters, excluding eol characters.
            auto e = mmapData + (olr.start - mmap_olr_start);
            auto pos = olr.start;
            for (size_t i = 0; (pos < olr.end) && (i < Overlap); ++e, ++pos) {
              // TODO: probably should not skip eol characters.  this is searching the tail end.
              if ((*e != bliss::io::BaseFileParser<decltype(mmapData)>::cr) && (*e != bliss::io::BaseFileParser<decltype(mmapData)>::eol)) ++i;
            }

            // new end.
            overlappedRange.end = pos;
            overlappedRange.intersect(fileRange);

            // clean up mapping
            unmap(mmapData, olr);
          }

          // check to see if there is enough memory for preloading.
          // if not, then we can't use Buffering even if L1Block chooses Buffering, so throw exception
          if (L1Buffering) {
            // check if we can buffer.  use at most 1/20 of the available memory because of kmer data size expansion.
            struct sysinfo memInfo;
            sysinfo (&memInfo);

            // linux freeram is limited to 10 to 15 % of physical, hence we just divide by 2, and that gives us 1/20 of the available memory
            if ((overlappedRange.size() * sizeof(T)) > (memInfo.freeram * memInfo.mem_unit / 2)) {
              // not enough memory available
              std::stringstream ss;
              ss << "ERROR in file buffering: requires " << (overlappedRange.size() * sizeof(T)) <<
                  " bytes; but (1/20th) available: " << (memInfo.freeram * memInfo.mem_unit / 2) <<
                  " free: " << memInfo.freeram << " unit " << memInfo.mem_unit;
              throw IOException(ss.str());
            }
          }

          // map the region of file to memory
          mmapData = map(overlappedRange);

          // now assign the data to the L1Block object. (if buffering, will copy).
          // use "+" since mmapData is a pointer so it's a random access iterator.
          L1Block.assign(mmapData + (overlappedRange.start - mmap_start), mmapData + (overlappedRange.end - mmap_start), overlappedRange);

          // if we are buffering, then it's okay to unmap the data, as we already have a copy
          if (L1Buffering)
          {
            unmap(mmapData, overlappedRange);
          }

          // mark data as loaded.
          loaded = true;

          // now configure L2 partitioner to traverse this data.  work with overlappedblocks because the last L2 block still needs to read the whole overlap region.
          configL2Partitioner(L2Partitioner, L2Blocks, overlappedRange, nConsumingThreads, L2BlockSize);


//          // now configurre the parser.  work with exclusive blocks.  DOES NOT BELONG HERE BECAUSE L1Partitioner may not be block or cyclic.
//#ifdef USE_MPI
//          if (this->comm == MPI_COMM_SELF)
//            L1parser.init_for_iterator(L1Block.begin(), fileRange, RangeType(mmap_start, overlappedRange.end), blockRange);
//          else
//            L1parser.init_for_iterator(L1Block.begin(), fileRange, RangeType(mmap_start, overlappedRange.end), blockRange, this->comm);
//#else
//          L1parser.init_for_iterator(L1Block.begin(), fileRange, RangeType(mmap_start, overlappedRange.end), blockRange);
//#endif


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
         *
         *            The method modifies the base range to align the boundary to record boundaries at both ends
         *
         *            This is accomplished by performing find_first_record to find the
         *              start, and doing the same on the next L2Block to find the end.
         *
         *
         * @note  not providing a no-arg default impl because the threading is setup before these calls and after construction.,
         *        and we do not want to get the id in this method from the threading library.
         *        OVERLAP AGNOSTIC.  this generates exclusive ranges.
         * @param tid   Thread Id.  L2 Partition is between threads within a group
         * @return  The range of the next L2 Partition Block.
         */
        RangeType getNextL2BlockRange(const size_t tid) {

          // NOTE: this is using block partitioner
        	// NOTE: when searching, need to extend by 1 record size - so that the search can potentially succeed.

          // data has to be loaded
          if (!this->loaded) throw std::logic_error("ERROR: getting L2Block range before L1Block is loaded");


          // get the parent range
          RangeType parentRange = this->L1Block.getRange();

          // now get the next, unmodified range.  this is atomic
          RangeType hint = this->L2Partitioner.getNext(tid);

          RangeType output(hint);

          if (hint.size() > 0) {

            // extend by 2 record
            RangeType start_search_range(hint.start, hint.end + 2 * recordSize);
            start_search_range.intersect(parentRange);
            RangeType end_search_range(hint.end, hint.end + 2 * recordSize);
            end_search_range.intersect(parentRange);



            try {

              // search for start and end
              output.start = L1parser.init_parser(this->L1Block.begin(), parentRange, parentRange, start_search_range);
              output.end = L1parser.find_first_record(this->L1Block.begin(), parentRange, parentRange, end_search_range);
            } catch (IOException& ex) {
              // either start or end are not found, so return nothing.

              BL_WARNING(ex.what());

              BL_ERRORF("curr range: chunk %lu, hint %lu-%lu, next %lu-%lu, srcData range %lu-%lu, mmap_range %lu-%lu\n",
                     this->L2BlockSize, hint.start, hint.end, end_search_range.start, end_search_range.end,
                     parentRange.start, parentRange.end, this->L1Block.getRange().start, this->L1Block.getRange().end);
              BL_ERRORF("got an exception search:  %s \n", ex.what());

              output.start = hint.end;
              output.end = hint.end;
            }

          }

          return output;

        }




        /**
         * @brief   get the DataBlock from L2 Partitioner for the given thread Id.
         * @note    Depending on data type of L2Block, the data may be buffered (copied to) memory.
         *          OVERLAP AWARE.  generated L2 block contains the overlap.
         * @param tid         id of thread calling this function.
         * @param L2BlockRange   The range of the L2Block to retrieve. computed by getNextL2BlockRange.
         * @return    reference to loaded L2Block for the specified thread id and range.
         */
        L2BlockType& getL2DataForRange(const size_t tid, const RangeType &L2BlockRange) {
          if (!loaded) throw std::logic_error("ERROR: getting L2Block data before file is loaded");

          // make sure the range is within the mmaped region.
          RangeType r = L2BlockRange;
          r.intersect(L1Block.getRange());

          // now compute the start and end iterators of this range.  use std::advance in case the iterator is not a random access iterator.
          auto s = L1Block.begin();
          std::advance(s, (r.start - L1Block.getRange().start));
          auto e = s;
          std::advance(e, r.size());

          // now handle overlap, excluding eol etc.
          // TODO: handle Overlap better.  should not skip overlap.
          if (Overlap > 0) {
            auto oe = e;
            for (size_t i = 0; (oe != L1Block.end()) && (i < Overlap); ++oe) {
              // TODO: probably should not skip EOL characters.  this is searching the tail end.
              if ((*oe != bliss::io::BaseFileParser<typename L1BlockType::iterator>::cr) && (*oe != bliss::io::BaseFileParser<typename L1BlockType::iterator>::eol)) ++i;
            }
            r.end += std::distance(e, oe);
            e = oe;
          }

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
         * @note      AGNOSTIC OF OVERLAP.  Partitioner always generate exclusive partitions.  results not aligned to record boundaries.
         * @param[in/out] partitioner L1 Partitioner to configure
         * @param[in] fileRange   portion of file to be partitioned by partitioner (0 to file size by defautl).
         * @param[in] nPartitions  number of partitions to create.
         * @param[in] chunkSize   the size of each partitioned chunk (== L1Block size).  default to 0, which means computed.
         */
        void configL1Partitioner(L1PartitionerT& partitioner, const RangeType& fileRange, const size_t& nPartitions, size_t& chunkSize) {
          // Set up the partitioner to partition the full file range with nPartitions.  Overlap is set to 0.  FileLoader handles overlap itself.

          // update the chunkSize  to at least 2x record size, (passes back to called).
          chunkSize = std::max(chunkSize, 2 * recordSize);

          // then configure the partitinoer to use the right number of threads.  Overlap is set to 0.  FileLoader handles overlap itself.
          chunkSize = partitioner.configure(fileRange, nPartitions, chunkSize, 0);
        }




        /**
         * @brief configure the L2 Partitioning
         * @details: This function does the following
         *            get the approximate record size and update the chunkSize for the L2Partitioner
         *            configure L2 partitioner for the memmapped range, number of threads, and L2Block Size.
         *            initialize the L2Block cache, one per thread.
         * @note      AGNOSTIC OF OVERLAP.  Partitioner always generate exclusive partitions.  results not aligned to record boundaries.
         * @param[in/out] partitioner L2 Partitioner to configure
         * @param[in/out] cache   L2 Block cache, 1 per thread.
         * @param[in] mmapRange   memory mapped region's range, to be partitioned by partitioner.
         * @param[in] nThreads    number of threads to create.
         * @param[in/out] chunkSize   the size of each partitioned chunk (== L1Block size).
         */
        void configL2Partitioner(L2PartitionerT& partitioner, L2BlockType* &cache, const RangeType& mmapRange, const size_t& nThreads, size_t& chunkSize) {

          //=====  adjust the L2BlockSize adaptively.
          // update the chunkSize  to at least 2x record size, (passes back to called).
          chunkSize = std::max(chunkSize, 2 * recordSize);

          // then configure the partitinoer to use the right number of threads.  Overlap is set to 0.  FileLoader handles overlap itself.
          chunkSize = partitioner.configure(L1Block.getRange(), nThreads, chunkSize, 0);

          // also configure the cache.
          if (cache == nullptr)
            cache = new L2BlockType[nThreads];
          else
            for (size_t i = 0; i < nThreads; ++i) {
              cache[i].clear();
            }
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
         * @return  pair of values, first is estimated record size, second is estimated sequence length in the records
         */
        // TODO: try to get rid of this - RecordSize is not used outside of this class, and only works for some file types.
        std::pair<size_t, size_t> getRecordSize(const int iterations = 10) {
          std::size_t counts[2] = {0, 0};  // do average.

          int rank = 0;
#ifdef USE_MPI
          MPI_Comm_rank(this->comm, &rank);
          if (rank == 0) {
#endif

        	//  map file directly on proc 1, determine Record size, then broadcast.

            // try mapping the first 9 page size worth  (typical, 4KB each, so 36KB.  chosen so that we can see if sequence length is greater than 32K (16 bit))
            auto search_range = this->getFileRange();
            search_range.end = std::min(9 * pageSize, search_range.end);

            // map the region.  (no need to align.)
            auto search_data = map(search_range);  // no alignment issue because we're at beginning of file.
            auto sstart = search_data;
            auto send = search_data + search_range.size();

            RangeType inmem(search_range);

            // create a L1 Parser local instance
            Parser<decltype(sstart)> local_parser;

            // For FASTAParser, this computes the sequence demarcation vector.  For FASTQ, this does nothing.
            size_t offset = 0;
#ifdef USE_MPI
            offset = local_parser.init_parser(sstart, this->getFileRange(), inmem, search_range, this->comm);  // in the case of FASTAParser, this creates the sequence demarcation vector.
#else
            offset = local_parser.init_parser(sstart, this->getFileRange(), inmem, search_range);
#endif
            // then find start, and adjust search range.  NOT NEEDED HERE since we are at beginning of file.

            // now initialize the search.  no need to call the first increment separately since we know that we have to start from seq id 0.
            typename Parser<decltype(sstart)>::SequenceType seq;
            size_t seq_data_len = 0;

            int i = 0;
            while ((i < iterations) && (sstart != send)) {
              // repeat for iterations or until no more sequences (empty sequence is returned).

              // get the next sequence
              seq = local_parser.get_next_record(sstart, send, offset);

              // get the char count in a record
              seq_data_len = std::distance(seq.seq_begin, seq.seq_end);

              if (seq_data_len >= (1 << 16))  {// larger than 16 bit DNA length  so stop and just return this.
                i = 1;
                break;
              }
              // else sum up.
              counts[0] += seq.record_size;
              counts[1] += seq_data_len;
              ++i;

            }
            if (i == 0) throw std::logic_error("ERROR: no record found");

            counts[0] /= i;
            counts[1] /= i;

            if (counts[0] == 0) throw std::logic_error("ERROR: estimated sequence data size is 0");
            if (counts[1] == 0) throw std::logic_error("ERROR: estimated record size is 0");

            // release the data.
            unmap(search_data, search_range);

#ifdef USE_MPI
          }
          int p;
          MPI_Comm_size(this->comm, &p);

          if (p > 1) MPI_Bcast(counts, 2, MPI_UNSIGNED_LONG, 0, this->comm);
#endif

          return std::make_pair(counts[0], counts[1]);
        }


        /// get number of estimated kmers, given k.
        // TODO: try to do this without getRecordSize.
        size_t getKmerCountEstimate(const int k) {

          if (recordSize == 0 || seqSizeInRecord == 0) std::tie(recordSize, seqSizeInRecord) = this->getRecordSize(10);

         size_t numRecords = (this->getFileRange().size() + recordSize - 1) / recordSize;

         size_t kmersPerRecord = seqSizeInRecord == 0 ? seqSizeInRecord - k + 1 : 0;  // fastq has quality.


          BL_DEBUGF("file range [%lu, %lu], recordSize = %lu, kmersPerRecord = %lu, numRecords = %lu",
                 this->getFileRange().start, this->getFileRange().end, recordSize, kmersPerRecord, numRecords);

          return numRecords * kmersPerRecord;
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
   constexpr unsigned char bliss::io::BaseFileParser<Iterator>::cr;
         */
        size_t getFileSize(const int& fd) throw (bliss::io::IOException) {

          size_t file_size = 0;

          int rank = 0;
#ifdef USE_MPI
          MPI_Comm_rank(comm, &rank);
          if (rank == 0) {
#endif

          struct stat64 filestat;

            // get the file state
            int ret = fstat64(fd, &filestat);

            // handle any error
            if (ret < 0) {
              throw IOException("ERROR in file size detection");
            }

            file_size = static_cast<size_t>(filestat.st_size);

#ifdef USE_MPI
          }

          int p;
          MPI_Comm_size(this->comm, &p);

          if (p > 1) MPI_Bcast(&file_size, 1, MPI_UNSIGNED_LONG, 0, comm);
#endif

          // return file size.
            return file_size;
        }


      protected:

        /**
         * @brief     map the specified portion of the file to memory.
         * @note      AGNOSTIC of overlaps
         * @param r   range specifying the portion of the file to map.
         * @return    memory address (pointer) to where the data is mapped.
         *            TODO: should return std::pair<PointerType, page aligned offset>
         */
        PointerType map(const RangeType &r) throw (IOException) {

          /// memory map.  requires that the starting position is block aligned.
          size_t block_start = RangeType::align_to_page(r, pageSize);   // mixed use of range semantic and pageSize byte

          size_t block_size = (r.end - block_start + Overlap) * sizeof(T);

          // NOT using MAP_POPULATE.  it slows things done when testing on single node.  NOTE HUGETLB not supported for file mapping.
          PointerType result = (PointerType)mmap64(nullptr, block_size,
                                     PROT_READ,
                                     MAP_PRIVATE, fileHandle,
                                     block_start * sizeof(T));  // TODO: range means bytes or number of elements?

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


          int madv_result = madvise(result, block_size, MADV_SEQUENTIAL | MADV_WILLNEED);
          if ( madv_result == -1) {
            std::stringstream ss;
            int myerr = errno;
            ss << "ERROR in madvise: " << myerr << ": " << strerror(myerr);
            throw IOException(ss.str());
          }

          return result;
        }




        /**
         * @brief unmaps a file region from memory
         * @note      AGNOSTIC of overlaps
         * @param d   The pointer to the memory address
         * @param r   The range that was mapped.
         */
        void unmap(PointerType &d, const RangeType &r) {

          munmap(d, (r.end - RangeType::align_to_page(r, pageSize) + Overlap) * sizeof(T));
        }

    };

  } /* namespace functional */
} /* namespace bliss */
#endif /* FILE_LOADER_HPP_ */
