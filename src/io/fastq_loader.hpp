/**
 * @file    fastq_loader.hpp
 * @ingroup io
 * @author  Tony Pan <tpan7@gatech.edu>
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

#include "bliss-config.hpp"

#include "common/base_types.hpp"
#include "common/sequence.hpp"
#include "io/file_loader.hpp"


namespace bliss
{
  namespace io
  {
    /// dummy class to indicate FASTQ format.
    struct FASTQ {


      /**
       * @class     bliss::io::FASTQ::SequenceId
       * @brief     represents a fastq sequence's id, also used for id of the FASTQ file, and for position inside a FASTQ sequence..
       * @details    this is set up as a union to allow easy serialization
       *            and parsing of the content.
       *            this keeps a 40 bit sequence ID, broken up into a 32 bit id and 8 bit significant bit id
       *                       an 8 bit file id
       *                       a 16 bit position within the sequence.
       *
       *            A separate FASTA version will have a different fields but keeps the same 64 bit total length.
       *
       *            file size is at the moment limited to 1TB (40 bits) in number of bytes, and up to 256 files.
       */
      union SequenceId
      {
          /// MANDATORY FIELD.  the concatenation of the id components as a single unsigned 64 bit field.  should use only 40 bits.
          uint64_t file_pos;

          /// the id field components.  anonymous struct.  CUSTOMIZABLE
          struct
          {
              /// sequence's id, lower 32 of 40 bits (potentially as offset in the containing file)
              uint32_t seq_id;
              /// sequence's id, upper 8 of 40 bits (potentially as offset in the containing file)
              uint8_t seq_id_msb;
              /// id of fastq file
              uint8_t file_id;
              /// offset within the read.  Default 0 refers to the whole sequence
              uint16_t pos;
          };

          /**
           * @brief sets the file position from the global position.
           * @note  MANDATORY
           */
          void set_file_pos(uint64_t pos) { file_pos = pos; }
      };

    };



    /**
     * @class bliss::io::FASTQParser
     * @brief Functoid encapsulating increment functionality traversing data in an iterator record by record, by parsing and searching for record boundaries.
     * @details   The purpose of this class is to abstract the increment operations in
     *    an iterator that transforms a block of data into a collection of FASTQ records. Using this class allows us to reuse
     *    an sequence parsing iterator, but potentially for different sequence record formats.
     *
     *    This class assumes that the iterator at start is pointing to the beginning of a FASTQ record already.
     *    Increment is done by walking through the data and recording the positions of the start and end of each
     *    of the 4 lines in a FASTQ record.  Because we are always searching to the end of a record (by reading 4 lines always),
     *    successive calls to increment always are aligned correctly with record boundaries
     *
     * @note  NOT THREAD SAFE.
     *
     * @tparam Iterator   The underlying iterator to be traversed to generate a Sequence
     * @tparam Quality    data type for quality scores for each base.  defaults to void to mean No Quality Score.
     */
    template<bool Quality = false>
    class FASTQParser
    {

      protected:
    	/// static constant for end of line.  note this is same for unicode as well
        static constexpr unsigned char eol = '\n';
    	/// static constant for carriage return.  note this is same for unicode as well
        static constexpr unsigned char cr = '\r';


        /**
         * @typedef RangeType
         * @brief   range object types
         */
        using RangeType = bliss::partition::range<size_t>;

        // NOT defining BaseIteratorType because different methods can then have different baseIteratorTypes

        // NOT defining SequenceType because different methods can then populate different sequence types.


      public:

        static constexpr bool has_quality = Quality;
        using SequenceIdType = bliss::io::FASTQ::SequenceId;

        /// default constructor.
        FASTQParser() {};

        /// default destructor
        virtual ~FASTQParser() {};

      protected:
        /**
         * @brief  search for first non-EOL character in a iterator, returns the stopping position as an iterator.  also records offset.
         * @details       iter can point to the previous EOL, or a nonEOL character.
         * @param[in/out] iter    iterator to advance
         * @param[in]     end     position to stop the traversal
         * @param[in/out] offset  the global offset of the sequence record within the file.  used as record id.
         * @return        iterator at the new position, where the Non EOL char is found, or end.
         */
        template <typename Iterator>
        inline Iterator findNonEOL(Iterator& iter, const Iterator& end, size_t &offset) const {
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
        template <typename Iterator>
        inline Iterator findEOL(Iterator& iter, const Iterator& end, size_t &offset) const {
          while ((iter != end) && ((*iter != eol) && (*iter != cr) ) ) {
            ++iter;
            ++offset;
          }
          return iter;
        }


        /**
         * @brief function to populate the quality score. defined only when Quality type is not void
         * @param[in] start   beginning of the quality score character sequence.
         * @param[in] end     end of the quality score character sequence.
         * @param[out] output  The output sequence record type to modify
         */
        template <typename Iterator, typename SequenceType, bool Q = Quality >
        inline typename std::enable_if<Q && SequenceType::has_quality::value>::type
        setQualityIterators(const Iterator & start, const Iterator & end, SequenceType& output) {
          output.qualBegin = start;
          output.qualEnd = end;
        }

        /**
         * @brief function to populate the quality score. defined only when Quality type is void
         * @param[in] start   beginning of the quality score character sequence.
         * @param[in] end     end of the quality score character sequence.
         * @param[out] output  The output sequence record type to modify
         */
        template <typename Iterator, typename SequenceType, bool Q = Quality >
        inline typename std::enable_if<!Q>::type
        setQualityIterators(const Iterator & start, const Iterator & end, SequenceType& output) {}

        /**
         * @brief check if quality iterator is at end.
         * @param sequence
         * @return
         */
        template <typename SequenceType, bool Q = Quality>
        inline typename std::enable_if<Q && SequenceType::has_quality::value, bool>::type
        isQualityIteratorAtEnd(SequenceType& seq) {
          return seq.qualBegin == seq.qualEnd;
        }

        /**
         * @brief check if quality iterator is at end.
         * @param[in] seq
         */
        template <typename SequenceType, bool Q = Quality>
        inline typename std::enable_if<!Q, bool>::type
        isQualityIteratorAtEnd(SequenceType& seq) { return false; }  // return false, so that we don't change the seq iterator.


        /**
         * @brief function to populate the sequence iterators.
         * @param[in] start   beginning of the sequence.
         * @param[in] end     end of the sequence.
         * @param[out] output  The output sequence record type to modify
         */
        template <typename Iterator, typename SequenceType>
        inline void setSequenceIterators(const Iterator & start, const Iterator & end, SequenceType& output) {
          output.seqBegin = start;
          output.seqEnd = end;
        }

        /**
         * @brief constructs an IOException object with the relevant debug data.
         * @param errType     string indicating the source of the error.  user choice
         * @param start       iterator pointing to beginning of the data in question
         * @param end         iterator pointing to end of the data in question
         * @param startOffset offset for the beginning of the data in question
         * @param endOffset   offset for the end of the data in question
         */
        template <typename Iterator>
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
        template <typename Iterator>
        void handleWarning(const std::string& errType, const Iterator &start, const Iterator &end, const size_t& startOffset, const size_t& endOffset) throw (bliss::io::IOException) {
          WARNING("WARNING: " << "did not find "<< errType << " in " << startOffset << " to " << endOffset);
        }


      public:
        /**
         * @brief increments the iterator to the beginning of the next record, while saving the current record in memory and update the record id.
         * @details   Called by the containing iterator during operator++ or operator+=
         *            Internally, this function finds the begin and end of the text lines 4 times,
         *            and populate the output sequence type with start and end positions of lines 2 and 4.
         *
         *            Each call to findNonEOL and findEOL may not move the iterator, if iter is already at end, or if iter already points to
         *            a nonEOL or EOL character, respectively.  This means that the output may have Sequence or Quality Score start and end iterators
         *            point to the block's end.  We check these conditions and throw an exception so the incomplete Sequence record does not get used
         *            inadvertently.
         *
         *            Iterator can be a forward iterator or better.
         *
         * @param[in/out] iter          source iterator, pointing to data to traverse
         * @param[in]     end           end of the source iterator - not to go past.
         * @param[in/out] offset        starting position in units of source character types (e.g. char) from the beginning of the file.
         * @param[in/out] output        updated sequence type, values updated.
         * @throw IOException           if parse failed, throw exception
         */
        template <typename Iterator, typename SequenceType>
        void increment(Iterator &iter, const Iterator &end, size_t &offset, SequenceType& output) throw (bliss::io::IOException)
        {
          //== first make a copy of iter so we can later use the original in exception handling.
          Iterator orig_iter = iter;
          size_t orig_offset = offset;

          // local variables, temp storage.
          Iterator lstart = end, lend = end;

          //==== parse 4 lines (always, because FASTQ record has 4 lines) for start and end.
          //==== okay to call findNonEOL or findEOL 4 times - if iter at end, won't advance. and okay for output to have end.
          // each call to findNonEOL will find first char, and trim the leading \n

          //== find 1st line, also the starting point.
          // find start of line, \nX or end.
          findNonEOL(iter, end, offset);
          // offset pointing to first start, save it as the id.  - either \n or end
          output.id.set_file_pos(offset);
          // then find the end of that line  - either \n or end.
          findEOL(iter, end, offset);

          // == find 2nd line, and save it.
          lstart = findNonEOL(iter, end, offset);
          // then find the end of that line  - either \n or end.
          lend = findEOL(iter, end, offset);
          setSequenceIterators(lstart, lend, output);

          // == find 3rd line, and discard it.
          findNonEOL(iter, end, offset);
          // then find the end of that line  - either \n or end.
          findEOL(iter, end, offset);

          // == find 4th line, and save it.  this depends on conditional definition of setQualityIterators.
          lstart = findNonEOL(iter, end, offset);
          // then find the end of that line  - either \n or end.
          lend = findEOL(iter, end, offset);
          // don't even bother to call if Quality type is void - avoids copy operations.
          setQualityIterators(lstart, lend, output);

//          // lend at this point is pointing to the last \n or at end.
//          if (iter != end) {  // if at \n, advance it by 1 position
//            ++iter;
//            ++offset;
//          }

          //== now check for error conditions.
          // if either sequence or if required quality score were not found, then failed parsing
          if (output.seqBegin == output.seqEnd) {
            // if nothing is found, throw a warning.

            handleWarning("sequence", orig_iter, lend, orig_offset, offset);
          } else {
            // sequence parsing is fine.  now check quality parsing.  if quality is not parsed, then this is an error.
            if (isQualityIteratorAtEnd(output)) {  // if quality score is not specified, the code below is never called.
              output.seqBegin = output.seqEnd;  // force seq data not to be used.
              handleError("required quality score", orig_iter, lend, orig_offset, offset);
            }
          }
          // leave iter and offset advanced.


        }

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
         *          ignore all characters before the first "\n" in the DataBlock
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
        template <typename Iterator>
        std::size_t findStart(const Iterator &_data, const RangeType &parentRange, const RangeType &inMemRange, const RangeType &searchRange) const
        throw (bliss::io::IOException) {

          typedef typename std::iterator_traits<Iterator>::value_type  ValueType;

          //== range checking
          if(!parentRange.contains(inMemRange)) throw std::invalid_argument("ERROR: Parent Range does not contain inMemRange");
          RangeType t = RangeType::intersect(searchRange, inMemRange); // intersection to bound target to between parent's ends.

          //== set up the iterator for later
          // initialize the counter
          std::size_t i = t.start;

          if (i == t.end) return t.end;

          // set iterator for the data
          Iterator iter(_data);
          Iterator end(_data);
          std::advance(iter, t.start - inMemRange.start);
          std::advance(end, t.end - inMemRange.start);


          //== now read 4 lines
          ValueType first[4] =
          { 0, 0, 0, 0 };
          std::size_t offsets[4] =
          { t.end, t.end, t.end, t.end };   // units:  offset from beginning of file.


          //== if at beginning of parent partition, treat as if previous line was EOL
          if ((t.start > parentRange.start) && ((*iter != eol) && (*iter !=cr))) { // at beginning of parent range, treat specially, since there is no preceding \n for the @
            // all other partitions will lose the part before the first "\n@" (will be caught by the previous partition)
            // if this is called to generate L1Block, then parentRange is the whole file, with first part start with @ (implicit \n prior).
            // if this is called to generate L2Block, then parentRange is the loaded L1 Block, already aligned to @, so previous is \n.

            // not at beginning of parent, and not an eol character, so find the next eol
            iter = findEOL(iter, end, i);

            // else at beginning of parent == 1st,  or an eol character == 1st is next.
          }
          if (i == t.end) return t.end;

          // now find the first non eol
          iter = findNonEOL(iter, end, i);
          if (i == t.end) return t.end;
          // assign to first.
          first[0] = *iter;
          offsets[0] = i;

          // lines 2 through 4
          for (int j = 1; j < 4; ++j) {
            iter = findEOL(iter, end, i);
            if (i == t.end) return t.end;
            iter = findNonEOL(iter, end, i);
            if (i == t.end) return t.end;
            first[j] = *iter;
            offsets[j] = i;
          }

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


    };

    /// template class' static variable definition (declared and initialized in class)
    template<bool Quality>
    constexpr unsigned char bliss::io::FASTQParser<Quality>::eol;
    template<bool Quality>
    constexpr unsigned char bliss::io::FASTQParser<Quality>::cr;
    template<bool Quality>
    constexpr bool bliss::io::FASTQParser<Quality>::has_quality;


    /**
     * @class FASTQLoader
     * @brief FileLoader subclass specialized for the FASTQ file format, uses CRTP to enforce interface consistency.
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
     *          using the offset as id allows us to avoid d assignment via prefix sum.
     *
     * @tparam  T                 type of each element read from file
     * @tparam  L1Buffering       bool indicating if L1 partition blocks should be buffered.  default to false
     * @tparam  L2Buffering       bool indicating if L2 partition blocks should be buffered.  default to true (to avoid contention between threads)
     * @tparam  L1PartitionerT    Type of the Level 1 Partitioner to generate the range of the file to load from disk
     * @tparam  L2PartitionerT    L2 partitioner, default to DemandDrivenPartitioner
     */
    template<typename T, bool Quality = false,
        bool L2Buffering = true,
        bool L1Buffering = false,
        typename L2PartitionerT = bliss::partition::DemandDrivenPartitioner<bliss::partition::range<size_t> >,
        typename L1PartitionerT = bliss::partition::BlockPartitioner<bliss::partition::range<size_t> > >
    using FASTQLoader = FileLoader<T, FASTQParser<Quality>, L2Buffering, L1Buffering, L2PartitionerT, L1PartitionerT>;

  } /* namespace io */
} /* namespace bliss */
#endif /* FASTQ_PARTITIONER_HPP_ */
