/**
 * @file    fastq_loader.hpp
 * @ingroup bliss::io
 * @author  tpan
 * @brief   contains a FileLoader subclass to perform distributed and concurrent FASTQ file access.
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */

#ifndef FASTQ_PARTITIONER_HPP_
#define FASTQ_PARTITIONER_HPP_

#include <sstream>
#include <iterator>  // for ostream_iterator
#include <algorithm> // for copy.

#include "config.hpp"

#include "common/base_types.hpp"
#include "io/file_loader.hpp"

namespace bliss
{
  namespace io
  {
    /// dummy class to indicate FASTQ format.
    struct FASTQ {};


    //==== subclass of FileLoader.  using static polymorphism via CRTP.


    /**
     * @class FASTQLoader
     * @brief FileLoader subclass specialized for the FASTQ file format
     * @details   FASTQLoader understands the FASTQ file format, and enforces that the
     *            L1 and L2 partition boundaries occur at FASTQ sequence record boundaries.
     *
     *            FASTQ files allow storage of per-base quality score and is typically used
     *              to store a collection of short sequences.
     *            Long sequences (such as complete genome) can be found in FASTA formatted files
     *              these CAN be read using standard FileLoader, provided an appropriate overlap
     *              is specified (e.g. for kmer reading, overlap should be k-1).
     *
     *          for reads, expected lengths are maybe up to 1K in length, and files contain >> 1 seq
     *
     *          Majority of the interface is inherited from FileLoader.
     *          This class modifies the behaviors of GetNextL1Block and getNextL2Block,
     *            as described in FileLoader's "Advanced Usage".  Specifically, the range is modified
     *            for each block by calling "findStart" to align the start and end to FASTA sequence
     *            record boundaries.
     *
     * @note    Processes 1 file at a time.  For paired end reads, a separate subclass is appropriate.
     *          Sequences are identified by the following set of attributes:
     *
     *             file id.          (via filename to id map)
     *             read id in file   (e.g. using the offset at whcih read occurs in the file as id)
     *             position in read  (useful for kmer)
     *
     *          using the offset as id allows us to avoid parallel id assignment via prefix sum.
     *
     * @tparam  T                 type of each element read from file
     * @tparam  L1Buffering       bool indicating if L1 partition blocks should be buffered.  default to false
     * @tparam  L2Buffering       bool indicating if L2 partition blocks should be buffered.  default to true (to avoid contention between threads)
     * @tparam  L1PartitionerT    Type of the Level 1 Partitioner to generate the range of the file to load from disk
     * @tparam  L2PartitionerT    L2 partitioner, default to DemandDrivenPartitioner
     */
    template<typename T,
        bool L2Buffering = true,
        bool L1Buffering = false,
        typename L2PartitionerT = bliss::partition::DemandDrivenPartitioner<bliss::partition::range<size_t> >,
        typename L1PartitionerT = bliss::partition::BlockPartitioner<bliss::partition::range<size_t> > >
    class FASTQLoader : public FileLoader<T, L2Buffering, L1Buffering, L2PartitionerT, L1PartitionerT,
                                          FASTQLoader<T, L2Buffering, L1Buffering, L2PartitionerT, L1PartitionerT> >
    {
      protected:
        /// base class type (FileLoader with the specific template parameters)
        typedef FileLoader<T, L2Buffering, L1Buffering, L2PartitionerT, L1PartitionerT,
            FASTQLoader<T, L2Buffering, L1Buffering, L2PartitionerT, L1PartitionerT> >    SuperType;

        friend class FileLoader<T, L2Buffering, L1Buffering, L2PartitionerT, L1PartitionerT,
                                FASTQLoader<T, L2Buffering, L1Buffering, L2PartitionerT, L1PartitionerT> >;

      public:
        //==== exposing types from base class
        /// Type of iterator for traversing the input data.
        typedef typename SuperType::InputIteratorType                   InputIteratorType;

        /// Type of range object.
        typedef typename SuperType::RangeType                            RangeType;

      protected:
        /**
         * @brief search for first occurrence of @ from an arbitrary offset in data, thus getting the position of the beginning of a FASTA sequence record.
         * @details
         *        A FASTA sequence record has the following general structure (expressed as regular expression.  not
         *          rigorously expressed)
         *
         *          [@][^\n]*\n[ATCGN]+\n[+].*\n[\^n]*\n
         *
         *        line 1 has "@" followed by the name of the sequence, which can contain any ascii character except newline
         *        line 2 is the biological sequence, with a constrained alphabet
         *        line 3 has "+" followed by either newline or by some other string, potentially the name of the sequence
         *        line 4 has the quality scores, encoded as ascii
         *
         *        because line 4 may start with "@" or "+", there can be ambiguity in identifying a line as the true start of a record.
         *        further, because lines 1, 3, and 4 may contain "@" and "+", there can be ambiguity when a L1 or L2 partitioning
         *          produces a DataBlock with first character in the middle of lines 1, 3, or 4.  This method implements a simple
         *          algorithm to disambiguate these cases by scanning multiple lines, to find the true start of a FASTQ record.
         *
         *        Specifically, we are looking for the pairing of a @ line and a + line separated by 1 line
         *          e.g. .*\n[@][^\n]*\n[ATCGN]+\n[+].*\n.*
         *
         *        this combination can occur at 2 places only: at beginning of read/sequence, or if read name contains @ and the L1 or L2 DataBlock happens to start there.
         *
         *        we can also look for + followed by @ in 2 lines.
         *          e.g. [ATCGN]*\n[+][^\n]*\n[^\n]+\n[@].*\n.*
         *
         *        this combination can also occur at 2 places only.
         *
         *    using just 1 of these patterns would require us to search through up to 7 lines.
         *      pattern 1 works better if partition starts on lines 3 or 4, and pattern 2 works better if partition starts on 1 or 2.
         *    using both patterns, we can just look at 4 lines.
         *
         *
         * Decoding logic:
         * standard pattern is @x*
         *                     l*
         *                     +x*
         *                     q*
         *                     @x*
         *                     l*
         *                     +x*
         *                     q*
         *                     ....
         *                  x could be @, +, or other char except "\n".
         *                  l is alphabet, so no @, +, or "\n"
         *                  q could be @, +, or other char within some ascii range., no "\n"
         *
         *
         *  Algorithm:
         *          ignore all characters before teh first "\n" in the DataBlock
         *          Read 4 lines starting from the first "\n" found and save the first characters.
         *            If the DataBlock is the first one in the parent DataBlock/Range,
         *            then treat the block as if it's prefixed with a "\n"
         *              e.g. first L1Block in the File's range,
         *              e.g. first L2Block in the parent L1Block (L1Block are record-aligned by the same process)
         *          inspect the first char, and look for either @..+ or +..@.
         *
         * @note    the range should be large enough so that we would not run into trouble with partition search not finding a starting position.
         *
         *
         * @tparam Iterator   type of iterator for data to traverse.  raw pointer if no L2Buffering, or vector's iterator if buffering.
         * @param _data       start of iterator.
         * @param parentRange   the "full" range to which the inMemRange belongs.  used to determine if the target is a "first" block and last block.  has to start and end with valid delimiters (@)
         * @param inMemRange the portion of the "full" range that's loaded in memory (e.g. in DataBlock).  does NOT have to start and end with @
         * @param searchRange the range in which to search for a record.  the start and end position in the file.  does NOT have to start and end with @. should be inside inMemRange
         * @return            position of the start of next read sequence (@).  if there is no complete sequence within the parentRange, return end of parentRange.
         * @throws            if no start is found, and search range does not cover the parent's end, then the search range does not include a complete record, throws IOException
         */
        template<typename Iterator>
        const std::size_t findStart(const Iterator &_data, const RangeType &parentRange, const RangeType &inMemRange, const RangeType &searchRange) const
        throw (bliss::io::IOException) {

          typedef typename std::iterator_traits<Iterator>::value_type  ValueType;

          //== range checking
          assert(parentRange.contains(inMemRange));
          RangeType t = RangeType::intersect(searchRange, inMemRange); // intersection to bound target to between parent's ends.
          if (t.start == t.end) return t.start;

          //== set up the iterator for later
          // initialize the counter
          std::size_t i = t.start;

          // set iterator for the data
          Iterator data(_data);

          // create flag indicating that the last character visited was a EOL char.
          bool wasEOL;

          //== if at beginning of parent partition, treat as if previous line was EOL
          if (t.start == parentRange.start) { // at beginning of parent range, treat specially, since there is no preceding \n for the @
            // all other partitions will lose the part before the first "\n@" (will be caught by the previous partition)
            // if this is called to generate L1Block, then parentRange is the whole file, with first part start with @ (implicit \n prior).
            // if this is called to generate L2Block, then parentRange is the loaded L1 Block, already aligned to @, so previous is \n.
            wasEOL = true;
          } else {
            std::advance(data, t.start - inMemRange.start);
            wasEOL = false;
          }
          //== at this point, data points to start of the search range.


          //=== remove leading newlines
          while ((i < t.end) && (*data == '\n')) {
            wasEOL = true;
            ++data;
            ++i;
          }

          //== if no more, end.
          if (i == t.end)  // this part only contained \n
            return t.end;





          //== now read 4 lines
          ValueType first[4] =
          { 0, 0, 0, 0 };
          std::size_t offsets[4] =
          { t.end, t.end, t.end, t.end };   // units:  offset from beginning of file.

          int currLineId = 0;        // current line id in the buffer

          // read the rest of the lines.
          ValueType c;
          while ((i < t.end) && currLineId < 4)
          {
            c = *data;

            // encountered a newline.  mark newline found.
            if (c == '\n')
            {
              wasEOL = true;  // toggle on
            }
            else
            {
              if (wasEOL) // previous char was newline, so c is first char in line
              {
                // save the first char and the offset.
                first[currLineId] = c;
                offsets[currLineId] = i;
                wasEOL = false;  // toggle off
                ++currLineId;
              }
            }
            //    else  // other characters in the line - don't care.

            ++data;
            ++i;
          }

//        DEBUGF("chars: %c %c %c %c, %lu - %lu, lines: %d\n", first[0], first[1], first[2], first[3], target.start, target.end, currLineId);

          //=== determine the position of a read by looking for @...+ or +...@
          // at this point, first[0] is pointing to first char after first newline, or in the case of first block, the first char.
          // everything in "first" array are first characters right after a newline.
          if (first[0] == '@' && first[2] == '+')  return offsets[0];
          if (first[1] == '@' && first[3] == '+')  return offsets[1];
          if (first[0] == '+' && first[2] == '@')  return offsets[2];
          if (first[1] == '+' && first[3] == '@')  return offsets[3];

          //=== If nothing was found, because we are at the end of the parentRange, then return the end.
          if (t.end == parentRange.end) {
            return t.end;
          }

          //=== Finally, if nothing was found, and we are not at the parentRange's end, then
          // the search range does not include a complete record.  this is an error and
          // exception is thrown.
          std::stringstream ss;
          ss << "WARNING in file processing: file segment \n" << "\t\t"
              << t << "\n\t\t(original " << searchRange << ")\n"
              << "\t\tdoes not contain valid FASTQ markers.\n String:";
          std::ostream_iterator<typename std::iterator_traits<Iterator>::value_type> oit(ss);
          Iterator s(_data);
          std::advance(s, (t.start - inMemRange.start));
          Iterator e(_data);
          std::advance(e, (t.end - inMemRange.start));
          std::copy(s, e, oit);

          throw bliss::io::IOException(ss.str());
        }



        /**
         * @brief get the next available range for the L1 partition block, given the L1 partition id
         * @details   This method is the CTRP implementation called by the base class' getNextL1BlockRange method.
         *            The method modifies the base range to align the boundary to record boundaries at both ends
         *
         *            This is accomplished by memmapping a Range and performing findStart to find the
         *              start, and doing the same on the next L1Block to find the end.
         *
         * @param pid   The L1 partition id.
         * @return      The range of the next L1 partition block.
         */
        RangeType getNextL1BlockRangeImpl(const size_t pid) {

          // get basic range
          RangeType hint = this->L1Partitioner.getNext(pid);

          // get the right shifted range
          size_t length = hint.size();

          // output data structure
          RangeType output(hint);

          if (length > 0) {

            RangeType next = RangeType::shiftRight(hint, length);
            next.intersect(this->fileRange);

            // get the combined ranges
            RangeType loadRange = RangeType::merge(hint, next);

            // memmap the content
            auto block_start = RangeType::align_to_page(loadRange, this->pageSize);
            auto mappedData = this->map(loadRange);  // this is page aligned.
            auto searchData = mappedData + (loadRange.start - block_start);


            // search for new start and end using findStart
            try {
              output.start = findStart(searchData, this->fileRange, loadRange, hint);
              output.end = findStart(searchData, this->fileRange, loadRange, next);

            } catch (IOException& ex) {
              // either start or end are not found so return an empty range.

              // TODO: need to handle this scenario better - should keep search until end.
              WARNINGF("%s\n", ex.what());

              WARNINGF("curr range: partition hint %lu-%lu, next %lu-%lu, file_range %lu-%lu\n",
                     hint.start, hint.end, next.start, next.end, this->fileRange.start, this->fileRange.end);
              WARNINGF("got an exception search for partition:  %s \n", ex.what());

              output.start = hint.end;
              output.end = hint.end;
            }

            // clean up and unmap
            this->unmap(mappedData, loadRange);
          }
          return output;
        }


        /**
         * @brief get the next L2 Partition Block given a thread id.
         * @details  This is the CTRP implementation method called by FileLoader's getNextL2BlockRange method
         *            This call needs to properly support concurrent range computation.
         *
         *            The method modifies the base range to align the boundary to record boundaries at both ends
         *
         *            This is accomplished by performing findStart to find the
         *              start, and doing the same on the next L2Block to find the end.
         *
         * @param tid   Thread Id.  L2 Partition is between threads within a group
         * @return  The range of the next L2 Partition Block.
         */
        RangeType getNextL2BlockRangeImpl(const size_t tid) {

          // data has to be loaded
          assert(this->loaded);

          // get the parent range
          RangeType parentRange = this->L1Block.getRange();

          // now get the next, unmodified range.  this is atomic
          RangeType hint = this->L2Partitioner.getNext(tid);

          size_t length = hint.size();
          RangeType output(hint);

          if (length > 0) {
            // get the RightShifted range for searching at the end.
            RangeType next = RangeType::shiftRight(hint, length);
            next.intersect(parentRange);


            try {
              // search for start and end
              output.start = findStart(this->L1Block.begin(), parentRange, parentRange, hint);
              output.end = findStart(this->L1Block.begin(), parentRange, parentRange, next);
            } catch (IOException& ex) {
              // either start or end are not found, so return nothing.

              WARNING(ex.what());

              printf("curr range: chunk %lu, hint %lu-%lu, next %lu-%lu, srcData range %lu-%lu, mmap_range %lu-%lu\n",
                     this->L2BlockSize, hint.start, hint.end, next.start, next.end, parentRange.start, parentRange.end, this->L1Block.getRange().start, this->L1Block.getRange().end);
              printf("got an exception search:  %s \n", ex.what());

              output.start = hint.end;
              output.end = hint.end;
            }

          }

          // DEBUG("range " << range << " chunk " << cs << " chunkPos " << chunkPos);

          return output;
        }

        /**
         * @brief   compute the approximate size of a data record in the file by reading a few records.
         * @details This is the CRTP implementation method that FileLoader baseclass's getRecordSize calls
         *
         *          This method provides an approximate size of a record, which is useful for caching
         *          and determining size of a DataBlock.
         *
         *          note that this is just the max of several records, found using the findStart method.
         *
         *          this relies on the data being loaded, so L1Block can be used for this search
         *
         * @note    at most as many threadss as specified in the c'tor will call this function.
         *          Also, each process calling this will be reading a different part of the file,
         *          so record size may vary from process to process.  but since L2 Partitioning
         *          affects each process independently, having different RecordSize, thus different
         *          L2BlockSize, should not be a big problem.
         *
         * @param count    number of records to read to compute the approximation.  default = 3.
         * @return  approximate size of a record.
         */
        size_t getRecordSizeImpl(int iterations = 3) {

          assert(this->loaded);

          std::size_t s, e;
          std::size_t ss = 1;
          RangeType loaded(this->L1Block.getRange());
          RangeType search(this->L1Block.getRange());

          s = findStart(this->L1Block.begin(), loaded, loaded, search);

          for (int i = 0; i < iterations; ++i) {
            search.start = s + ss;   // advance by 1, in order to search for next entry.
            search.intersect(loaded);
            e = findStart(this->L1Block.begin(), loaded, loaded, search);
            ss = std::max(ss, (e - s));
            s = e;
          }

          // DEBUG("Sequence size is " << ss);
          return ss;
        }

      public:
        //== public methods.

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
        explicit FASTQLoader(const MPI_Comm& _comm, const std::string &_filename, const size_t _nThreads = 1, const size_t _L2BlockSize = 1) throw (IOException)
                             : SuperType(_comm, _filename, _nThreads, _L2BlockSize)
        {}
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
        FASTQLoader(const size_t _nProcs, const int _rank, const std::string &_filename,
                    const size_t _nThreads = 1, const size_t _L2BlockSize = 1) throw (IOException)
                            : SuperType(_nProcs, _rank, _filename, _nThreads, _L2BlockSize)
        {};

        /**
         * @brief default destructor.  unloads the data and closes the file
         */
        virtual ~FASTQLoader() {};

        /// Removed default copy constructor
        FASTQLoader(const FASTQLoader<T, L2Buffering, L1Buffering, L2PartitionerT, L1PartitionerT>& other) = delete;
        /// Removed default copy assignment operator
        FASTQLoader& operator=(const FASTQLoader<T, L2Buffering, L1Buffering, L2PartitionerT, L1PartitionerT>& other) = delete;




        /**
         * @brief move constructor.  here to ensure that the DataBlock instances are moved.
         * @param other FASTQLoader to move
         */
        FASTQLoader(FASTQLoader<T, L2Buffering, L1Buffering, L2PartitionerT, L1PartitionerT>&& other) : SuperType(std::forward(other)) {};

        /**
         * @brief move assignment operator.  here to ensure that the DataBlock instances are moved.
         *
         * @param other   FASTQLoader to move
         * @return        updated object with moved data
         */
        FASTQLoader<T, L2Buffering, L1Buffering, L2PartitionerT, L1PartitionerT>& operator=(FASTQLoader<T, L2Buffering, L1Buffering, L2PartitionerT, L1PartitionerT>&& other) {
          if (this != &other) {
            this->SuperType::operator =(other);
          }
          return *this;
        }

    };

  } /* namespace io */
} /* namespace bliss */
#endif /* FASTQ_PARTITIONER_HPP_ */
