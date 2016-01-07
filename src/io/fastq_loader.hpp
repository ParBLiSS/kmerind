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

    };

    template <typename Iterator>
    class FASTQSequence : public ::bliss::common::Sequence<Iterator> {

      public:
      /// Iterator type for traversing the sequence.
      typedef Iterator IteratorType;
      /// type for the id struct/union
      typedef typename ::bliss::common::Sequence<Iterator>::IdType IdType;

      /// begin iterator for the sequence
      Iterator qualBegin;
      /// end iterator for the sequence.
      Iterator qualEnd;



      FASTQSequence() = default;

      FASTQSequence(IdType const & _id, size_t const & _length, size_t const& _offset, size_t const& _local_offset, Iterator const & _begin, Iterator const & _end) :
        ::bliss::common::Sequence<Iterator>(_id, _length, _offset, _local_offset,  _begin, _end) {}

      FASTQSequence(FASTQSequence const & other) :
        ::bliss::common::Sequence<Iterator>(other.id, other.length, other.seq_begin_offset, other.local_offset, other.seqBegin, other.seqEnd), qualBegin(other.qualBegin), qualEnd(other.qualEnd) {}

      FASTQSequence& operator=(FASTQSequence const & other) {
        this->id = other.id;
        this->length = other.length;
        this->seq_begin_offset = other.seq_begin_offset;
        this->local_offset = other.local_offset;
        this->seqBegin = other.seqBegin;
        this->seqEnd = other.seqEnd;
        qualBegin = other.qualBegin;
        qualEnd = other.qualEnd;

        return *this;
      }

      void set_quality(Iterator const & _begin, Iterator const & _end) {
        qualBegin = _begin;
        qualEnd = _end;
      }

      static constexpr bool has_quality() { return true; }
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
     * @note  all sequences now have the quality score iterators.
     * @note  NOT THREAD SAFE.
     *
     * @tparam Iterator   The underlying iterator to be traversed to generate a Sequence
     */
    template <typename Iterator>
    class FASTQParser : public bliss::io::BaseFileParser<Iterator >
    {
        // let another BaseFileParser with some other iterator type be a friend so we can convert.
        template <typename Iterator2>
        friend class FASTAParser;

      public:
        using SequenceType = bliss::io::FASTQSequence<Iterator>;  // redefined SequenceType
        using SequenceIdType = typename SequenceType::IdType;  // redefined SequenceType

        /// default constructor.
        FASTQParser() : bliss::io::BaseFileParser<Iterator>() {};

        /// default destructor
        ~FASTQParser() {};

        /// converting constructor.
        template <typename Iterator2>
        FASTQParser(FASTQParser<Iterator2> const & other) {}
        /// converting assignment operator that can transform the base iterator type.
        template <typename Iterator2>
        FASTQParser<Iterator>& operator=(FASTQParser<Iterator2> const & other) { return *this; }




      protected:

        /**
         * @typedef RangeType
         * @brief   range object types
         */
        using RangeType = bliss::partition::range<size_t>;

        /**
         * @brief check if quality iterator is at end.
         * @param sequence
         * @return
         */
        bool isQualityIteratorAtEnd(SequenceType& seq) {
          return seq.qualBegin == seq.qualEnd;
        }

      public:

        /**
         * @brief given a block, find the starting point that best aligns with the first sequence object (here, the first @)
         * @details
         *     search for first occurrence of @ from an arbitrary offset in data, thus getting the position of the beginning of a FASTQ sequence record within the block.
         *
         *        A FASTQ sequence record has the following general structure (expressed as regular expression.  not
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
          if ((t.start > parentRange.start) && ((*iter != ::bliss::io::BaseFileParser<Iterator>::eol) && (*iter != ::bliss::io::BaseFileParser<Iterator>::cr))) { // at beginning of parent range, treat specially, since there is no preceding \n for the @
            // all other partitions will lose the part before the first "\n@" (will be caught by the previous partition)
            // if this is called to generate L1Block, then parentRange is the whole file, with first part start with @ (implicit \n prior).
            // if this is called to generate L2Block, then parentRange is the loaded L1 Block, already aligned to @, so previous is \n.

            // not at beginning of parent, and not an eol character, so find the next eol
            iter = this->findEOL(iter, end, i);

            // else at beginning of parent == 1st,  or an eol character == 1st is next.
          }
          if (i == t.end) return t.end;

          // now find the first non eol
          iter = this->findNonEOL(iter, end, i);
          if (i == t.end) return t.end;
          // assign to first.
          first[0] = *iter;
          offsets[0] = i;

          // lines 2 through 4
          for (int j = 1; j < 4; ++j) {
            iter = this->findEOL(iter, end, i);
            if (i == t.end) return t.end;
            iter = this->findNonEOL(iter, end, i);
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
        size_t increment(Iterator &iter, const Iterator &end, size_t &offset, size_t &seq_offset, SequenceType& output) throw (bliss::io::IOException)
        {
          //== first make a copy of iter so we can later use the original in exception handling.
          Iterator orig_iter = iter;
          size_t orig_offset = offset;

          // local variables, temp storage.
          Iterator sstart = end, send = end;
          Iterator lstart = end, lend = end;

          size_t record_start_offset, seq_start_offset;

          //==== parse 4 lines (always, because FASTQ record has 4 lines) for start and end.
          //==== okay to call findNonEOL or findEOL 4 times - if iter at end, won't advance. and okay for output to have end.
          // each call to findNonEOL will find first char, and trim the leading \n

          //== find 1st line, also the starting point.
          // find start of line, \nX or end.
          this->findNonEOL(iter, end, offset);
          // offset pointing to first start, save it as the id.  - either \n or end
          record_start_offset = offset;
          // then find the end of that line  - either \n or end.
          this->findEOL(iter, end, offset);

          // == find 2nd line, and save it.
          sstart = this->findNonEOL(iter, end, offset);
          seq_start_offset = offset;
          // then find the end of that line  - either \n or end.
          send = this->findEOL(iter, end, offset);


          // == find 3rd line, and discard it.
          this->findNonEOL(iter, end, offset);
          // then find the end of that line  - either \n or end.
          this->findEOL(iter, end, offset);

          // == find 4th line, and save it.  this depends on conditional definition of setQualityIterators.
          lstart = this->findNonEOL(iter, end, offset);
          // then find the end of that line  - either \n or end.
          lend = this->findEOL(iter, end, offset);


          // skip over any \n.
          this->findNonEOL(iter, end, offset);

          // store where the actual DNA sequence starts relative to the beginnig of the record.
          output = SequenceType(SequenceIdType(record_start_offset),
                                offset - record_start_offset,
                                seq_start_offset - record_start_offset,
                                seq_start_offset - record_start_offset,
                                sstart, send);
          // don't even bother to call if Quality type is void - avoids copy operations.
          output.set_quality(lstart, lend);


          //== now check for error conditions.
          // if either sequence or if required quality score were not found, then failed parsing
          if (output.seqBegin == output.seqEnd) {
            // if nothing is found, throw a warning.

            this->handleWarning("sequence", orig_iter, lend, orig_offset, offset);
          } else {
            // sequence parsing is fine.  now check quality parsing.  if quality is not parsed, then this is an error.
            if (isQualityIteratorAtEnd(output)) {  // if quality score is not specified, the code below is never called.
              output.seqBegin = output.seqEnd;  // force seq data not to be used.
              this->handleError("required quality score", orig_iter, lend, orig_offset, offset);
            }
          }

          // leave iter and offset advanced.
          ++seq_offset;
          return seq_offset;
        }


        /**
         * @brief initializes the Parser object for sequencesIterator.  primarily to search for index within global array (e.g. FASTA)
         * @note SequenceType is different than that in base type, so this is an overload, not override.
         *
         * @param iter
         * @param parentRange
         * @param inMemRange
         * @param searchRange
         * @return   index in process-local (shared) sequence vector.
         */
        size_t increment(Iterator &iter, const Iterator &end, size_t &offset, SequenceType& output) throw (bliss::io::IOException) {
          size_t result = 0;

          increment(iter, end, offset, result, output);

          return 0;
        };
    };


    /**
     * @class FASTQLoader
     * @brief FileLoader subclass specialized for the FASTQ file format, uses CRTP to enforce interface consistency.
     * @details   FASTQLoader understands the FASTQ file format, and enforces that the
     *            L1 and L2 partition boundaries occur at FASTQ sequence record boundaries.
     *
     *            FASTQ files allow storage of per-base quality score and is typically used
     *              to store a collection of short sequences.
     *            Long sequences (such as complete genome) can be found in FASTQ formatted files
     *              these CAN be read using standard FileLoader, provided an appropriate overlap
     *              is specified (e.g. for kmer reading, overlap should be k-1).
     *
     *          for reads, expected lengths are maybe up to 1K in length, and files contain >> 1 seq
     *
     *          Majority of the interface is inherited from FileLoader.
     *          This class modifies the behaviors of GetNextL1Block and getNextL2Block,
     *            as described in FileLoader's "Advanced Usage".  Specifically, the range is modified
     *            for each block by calling "findStart" to align the start and end to FASTQ sequence
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
    template<typename T,
        bool L2Buffering = false,
        bool L1Buffering = true,
        typename L2PartitionerT = bliss::partition::DemandDrivenPartitioner<bliss::partition::range<size_t> >,
        typename L1PartitionerT = bliss::partition::BlockPartitioner<bliss::partition::range<size_t> > >
    using FASTQLoader = FileLoader<T, 0, FASTQParser, L2Buffering, L1Buffering, L2PartitionerT, L1PartitionerT>;

  } /* namespace io */
} /* namespace bliss */
#endif /* FASTQ_PARTITIONER_HPP_ */
