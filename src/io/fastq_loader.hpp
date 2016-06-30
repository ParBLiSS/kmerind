/*
 * Copyright 2015 Georgia Institute of Technology
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/**
 * @file    fastq_loader.hpp
 * @ingroup io
 * @author  Tony Pan <tpan7@gatech.edu>
 * @brief   contains a FileLoader subclass to perform distributed and concurrent FASTQ file access.
 *
 */

#ifndef FASTQ_PARTITIONER_HPP_
#define FASTQ_PARTITIONER_HPP_

#include <cmath>

#include <sstream>
#include <iterator>  // for ostream_iterator
#include <algorithm> // for copy.

#include <mxx/comm.hpp>
#include <mxx/shift.hpp>

#include "bliss-config.hpp"

#include "common/base_types.hpp"
#include "common/sequence.hpp"
#include "io/file_loader.hpp"
#include <type_traits>

namespace bliss
{
  namespace io
  {
    /// dummy class to indicate FASTQ format.
    struct FASTQ {

    };


    template <typename Iterator, bool HasQuality = true>
    class FASTQSequence : public ::bliss::common::Sequence<Iterator> {
        // NOTE: uses the default SequenceId type instead of ShortSequenceKmerType, because we are dealing with file coordinate right now and not Kmer's position
      protected:
        using BaseType = ::bliss::common::Sequence<Iterator>;

      public:
      /// Iterator type for traversing the sequence.
      typedef Iterator IteratorType;
      /// type for the id struct/union
      typedef typename ::bliss::common::Sequence<Iterator>::IdType IdType;

      /// begin iterator for the sequence
      Iterator qual_begin;
      /// end iterator for the sequence.
      Iterator qual_end;

      /**
       * @brief << operator to write out SequenceId
       * @param[in/out] ost   output stream to which the content is directed.
       * @param[in]     seq   Sequence object to write out
       * @return              output stream object
       */
      friend std::ostream& operator<<(std::ostream& ost, const FASTQSequence & seq)
      {
        ost << " FASTQ Sequence: id=[" << seq.id << "] record_size=" << seq.record_size << " seq_offset=" << seq.seq_offset;
        return ost;
      }


      FASTQSequence() = default;

      template <bool Q = HasQuality>
      FASTQSequence(IdType const & _id, size_t const & _record_size, size_t const& _seq_offset, Iterator const & _begin, Iterator const & _end,
                    typename ::std::enable_if<!Q>::type* = 0) :
        BaseType(_id, _record_size, _seq_offset,  _begin, _end) {}

      FASTQSequence(IdType const & _id, size_t const & _record_size, size_t const& _seq_offset,
                    Iterator const & _seq_begin, Iterator const & _seq_end,
                    Iterator const & _qual_begin, Iterator const & _qual_end) :
        BaseType(_id, _record_size, _seq_offset,  _seq_begin, _seq_end), qual_begin(_qual_begin), qual_end(_qual_end) {}


      FASTQSequence(FASTQSequence const & other) :
        BaseType(other.id, other.record_size, other.seq_offset, other.seq_begin, other.seq_end), qual_begin(other.qual_begin), qual_end(other.qual_end) {}

      FASTQSequence& operator=(FASTQSequence const & other) {

        this->id = other.id;
        this->record_size = other.record_size;
        this->seq_offset = other.seq_offset;
        this->seq_begin = other.seq_begin;
        this->seq_end = other.seq_end;

        qual_begin = other.qual_begin;
        qual_end = other.qual_end;

        return *this;
      }

//      void set_quality(Iterator const & _begin, Iterator const & _end) {
//        qual_begin = _begin;
//        qual_end = _end;
//      }

      static constexpr bool has_quality() { return HasQuality; }
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
     * @note  DEFAULT DOES NOT SUPPORT BLOCK PARTITION.
     * @tparam Iterator   The underlying iterator to be traversed to generate a Sequence
     */
    template <typename Iterator>
    class SequentialFASTQParser : public bliss::io::BaseFileParser<Iterator >
    {
        // let another BaseFileParser with some other iterator type be a friend so we can convert.
        template <typename Iterator2>
        friend class SequentialFASTQParser;

      protected:

      public:
        using SequenceType = bliss::io::FASTQSequence<Iterator>;  // redefined SequenceType
        using SequenceIdType = typename SequenceType::IdType;  // redefined SequenceType

        /// default constructor.

        SequentialFASTQParser() : bliss::io::BaseFileParser<Iterator>() {};

        /// default destructor
        virtual ~SequentialFASTQParser() {};

        /// converting constructor.
        template <typename Iterator2>
        SequentialFASTQParser(SequentialFASTQParser<Iterator2> const & other) : SequentialFASTQParser() {}
        /// converting assignment operator that can transform the base iterator type.
        template <typename Iterator2>
        SequentialFASTQParser<Iterator>& operator=(SequentialFASTQParser<Iterator2> const & other) { return *this; }

        // inherited reset and init_parser and find_overlap_end
        using bliss::io::BaseFileParser<Iterator>::find_overlap_end;
        using bliss::io::BaseFileParser<Iterator>::reset;
        using bliss::io::BaseFileParser<Iterator>::init_parser;
        using bliss::io::BaseFileParser<Iterator>::should_parse;


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
          return seq.qual_begin == seq.qual_end;
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
        virtual std::size_t find_first_record(const Iterator &_data, const RangeType &parentRange,
                                              const RangeType &inMemRange, const RangeType &searchRange)
        {

          typedef typename std::iterator_traits<Iterator>::value_type  ValueType;

          //== range checking
          if(!parentRange.contains(inMemRange)) {
            ::std::cout << "parent: " << parentRange << " , in mem: " << inMemRange << ::std::endl;
            throw std::invalid_argument("ERROR: Parent Range does not contain inMemRange");
          }
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
          if ((t.start > parentRange.start) &&
              ((*iter != ::bliss::io::BaseFileParser<Iterator>::eol) &&
               (*iter != ::bliss::io::BaseFileParser<Iterator>::cr))) { // at beginning of parent range, treat specially, since there is no preceding \n for the @
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


          // lines 2 through 4.  reorganized so that O3 optimization by gcc 5.2.1 does not skip over the whole loop (compiler bug?)
          for (int j = 1; j < 4; ++j) {
            iter = this->findEOL(iter, end, i);
            iter = this->findNonEOL(iter, end, i);
            offsets[j] = i;
            if (i != t.end) first[j] = *iter;
          }
          if (i == t.end) return t.end;

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
          ss << "ERROR in file processing: file segment \n" << "\t\t"
              << t << "\n\t\t(original " << searchRange << ")\n"
              << "\t\t(actual search " << t << ")\n"
              << "\t\tdoes not contain valid FASTQ markers.\n String:"
              << "first chars are " << first[0] << "," << first[1] << "," << first[2] << "," << first[3];
          std::ostream_iterator<typename std::iterator_traits<Iterator>::value_type> oit(ss);
          Iterator s(_data);
          std::advance(s, (t.start - inMemRange.start));
          Iterator e(_data);
          std::advance(e, (t.end - inMemRange.start));
          std::copy(s, e, oit);

          throw ::std::logic_error(ss.str());
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
        SequenceType get_next_record(Iterator &iter, const Iterator &end, size_t &offset)
        {

          if ((iter != end) && (*iter != '@'))
            this->handleError("missing @ on first line. ", iter, end, offset, offset + ::std::distance(iter, end));

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

//#ifdef USE_MPI
//         int rank = 0;
//         MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//         std::cout << "rank: " << rank << " offset " << offset << std::endl;
//#endif

          //== find 1st line, also the starting point.
          // find start of line, \nX or end.
          this->findNonEOL(iter, end, offset);
          if ((iter != end) && (*iter != '@'))
            this->handleError("missing @ on first line. ", orig_iter, end, orig_offset, orig_offset + ::std::distance(orig_iter, end));
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
          if ((iter != end) && (*iter != '+'))
            this->handleError("missing + on third line. ", orig_iter, end, orig_offset, orig_offset + ::std::distance(orig_iter, end));
          // then find the end of that line  - either \n or end.
          this->findEOL(iter, end, offset);

          // == find 4th line, and save it.  this depends on conditional definition of setQualityIterators.
          lstart = this->findNonEOL(iter, end, offset);
          // then find the end of that line  - either \n or end.
          lend = this->findEOL(iter, end, offset);


          // skip over any remaining \n.
          this->findNonEOL(iter, end, offset);

          // store where the actual DNA sequence starts relative to the beginning of the record.

          //== now check for error conditions.
          // if either sequence or if required quality score were not found, then failed parsing
          // likely the record was truncated.
          if ((sstart == send) || (lstart == lend)) {
            this->handleWarning("truncated record? missing seq or quality", orig_iter, lend, orig_offset, offset);
            return SequenceType(SequenceIdType(record_start_offset),
                                offset - record_start_offset,
                                seq_start_offset - record_start_offset,
                                sstart, send, lstart, lend);
          } else if (::std::distance(sstart, send) != ::std::distance(lstart, lend)) {
            this->handleError("truncated record? seq and qual differ in length", orig_iter, lend, orig_offset, offset);
          }

          return SequenceType(SequenceIdType(record_start_offset),
                              offset - record_start_offset,
                              seq_start_offset - record_start_offset,
                              sstart, send, lstart, lend);
        }



        /**
         * @brief   get the average record size in the supplied range
         * @param _data.    points to start of in memory range.
         * @param parentRange     the complete range to which the inMemRange belongs
         * @param inMemRange    the part of parentRange that is in memory, represented by _data
         * @param searchRange   the part of inMemRange to search.
         * @return  return the records size and the internal data size
         */
        virtual ::std::pair<size_t, size_t> get_record_size(const Iterator &_data, const RangeType &parentRange, const RangeType &inMemRange, const RangeType &searchRange, size_t const count) {

          if (searchRange.size() == 0) throw std::logic_error("calling FASTQParser get_record_size with 0 sized searchRange for iterator.");
          if (inMemRange.size() == 0) throw std::logic_error("calling FASTQParser get_record_size with 0 sized inMemRange for iterator.");

          if (count == 0) throw std::invalid_argument("ERROR: called FASTQParser get_record_size with count == 0");
          if (!(parentRange.contains(inMemRange))) throw std::invalid_argument("ERROR: parentRange does not container inMemRange");
          if (!(inMemRange.contains(searchRange)))  throw std::invalid_argument("ERROR: inMemRange does not container searchRange");

          // find the first record that starts within search range
          size_t start = this->find_first_record(_data, parentRange, inMemRange, searchRange);

          Iterator local_data = _data + start - inMemRange.start;


          if (*local_data != '@') throw std::logic_error("calling FASTQParser get_record_size with _data not pointing to start of a record.  call init_parser first.");

          Iterator local_end = _data + searchRange.end - inMemRange.start;
          size_t local_offset = searchRange.start;

          // now initialize the search.  no need to call the first increment separately since we know that we have to start from seq id 0.
          size_t seq_data_len = 0;
          size_t record_size = 0;


          size_t i = 0;
          while ((i < count) && (local_data != local_end)) {
            // repeat for iterations or until no more sequences (empty sequence is returned).

            // get the next sequence
            auto seq = this->get_next_record(local_data, local_end, local_offset);

            if (seq.record_size == 0) {  // nothing found, time to stop.
              break;
            }

            record_size += seq.record_size;
            // get the char count in a record
            seq_data_len += std::distance(seq.seq_begin, seq.seq_end);

            ++i;
          }
          if (i == 0) throw std::logic_error("ERROR: no record found");

          record_size /= i;
          seq_data_len /= i;

          if (record_size == 0) throw std::logic_error("ERROR: estimated sequence data size is 0");
          if (seq_data_len == 0) throw std::logic_error("ERROR: estimated record size is 0");

          return std::make_pair(record_size, seq_data_len);
        }

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
    class FASTQParser : public bliss::io::SequentialFASTQParser<Iterator >
    {
        // let another BaseFileParser with some other iterator type be a friend so we can convert.
        template <typename Iterator2>
        friend class FASTQParser;

      protected:

      public:
        using SequenceType = bliss::io::FASTQSequence<Iterator>;  // redefined SequenceType
        using SequenceIdType = typename SequenceType::IdType;  // redefined SequenceType

        /// default constructor.

        FASTQParser() : bliss::io::SequentialFASTQParser<Iterator>() {};

        /// default destructor
        virtual ~FASTQParser() {};

        /// converting constructor.
        template <typename Iterator2>
        FASTQParser(FASTQParser<Iterator2> const & other) : FASTQParser() {}
        /// converting assignment operator that can transform the base iterator type.
        template <typename Iterator2>
        FASTQParser<Iterator>& operator=(FASTQParser<Iterator2> const & other) { return *this; }

        // inherited reset and init_parser and find_overlap_end
        using bliss::io::SequentialFASTQParser<Iterator>::find_overlap_end;
        using bliss::io::SequentialFASTQParser<Iterator>::reset;
        using bliss::io::SequentialFASTQParser<Iterator>::init_parser;
        using bliss::io::SequentialFASTQParser<Iterator>::find_first_record;
        using bliss::io::SequentialFASTQParser<Iterator>::get_next_record;
        using bliss::io::SequentialFASTQParser<Iterator>::get_record_size;
        using bliss::io::SequentialFASTQParser<Iterator>::isQualityIteratorAtEnd;
        using bliss::io::SequentialFASTQParser<Iterator>::should_parse;


      protected:

        /**
         * @typedef RangeType
         * @brief   range object types
         */
        using RangeType = bliss::partition::range<size_t>;


      public:

#ifdef USE_MPI
       /// initializes the parser.  only useful for FASTA parser for now.  Assumes searchRange do NOT overlap between processes.   _data points to first element in mem.
       virtual std::size_t init_parser(const Iterator &_data, const RangeType &parentRange, const RangeType &inMemRange, const RangeType &searchRange, const mxx::comm& _comm)
       {
         return this->find_first_record(_data, parentRange, inMemRange, searchRange, _comm);
       };
#endif

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
#ifdef USE_MPI
       virtual std::size_t find_first_record(const Iterator &_data, const RangeType &parentRange,
                                             const RangeType &inMemRange, const RangeType &searchRange,
                                             mxx::comm const & comm)
       {
         // require non-overlapping search ranges.  only works with block decomposition
         // NOTE: do not rely on user supplying non-overlapping search ranges.


         typedef typename std::iterator_traits<Iterator>::value_type  ValueType;

         //== range checking
         if(!parentRange.contains(inMemRange)) {
           ::std::cout << "parent: " << parentRange << " , in mem: " << inMemRange << ::std::endl;
           throw std::invalid_argument("ERROR: Parent Range does not contain inMemRange");
         }


         // ensure non-overlapping search ranges.
         // left shift to query for the previous character.
         size_t next_start = searchRange.size() > 0 ? searchRange.start : parentRange.end;
         next_start = ::mxx::exscan(next_start, [](size_t const & x, size_t const & y){
           return std::min(x, y);
         }, comm.reverse());
         if (comm.rank() == (comm.size() - 1)) // last rank
           next_start = searchRange.end;

         RangeType t;
         t.start = std::min(next_start, searchRange.start);
         t.end = std::min(next_start, searchRange.end);
         t = RangeType::intersect(t, inMemRange); // intersection to bound target to between parent's ends.


         std::cout << "rank " << comm.rank()<< "   parent " << parentRange << std::endl;
         std::cout << "rank " << comm.rank()<< "   inmem " << inMemRange << std::endl;
         std::cout << "rank " << comm.rank()<< "   search " << searchRange << std::endl;
         std::cout << "rank " << comm.rank()<< "   intersect " << t << std::endl;


         // set iterator for the data
         Iterator iter(_data);
         Iterator end(_data);
         std::advance(iter, t.start - inMemRange.start);
         std::advance(end, t.end - inMemRange.start);

//         if (inMemRange.size() > 0) std::cout << "rank " << comm.rank() << " orig " << *(end-1) << " |";


         // get the rank to send the previous char.  (from low to high rank, possibly skipping ranks if there is no data.
         int target_rank = t.size() > 0 ? comm.rank() : (comm.size() - 1);
         target_rank = ::mxx::exscan(target_rank, [](size_t const & x, size_t const & y){
           return std::min(x, y);
         }, comm.reverse());
         if (comm.rank() == (comm.size() - 1)) // last rank send to self, but not really.
           target_rank = (comm.size() - 1);


         // move the prev character.
         std::vector<size_t> send_counts(comm.size(), 0);

         std::vector<ValueType> tmps(1, 0);
         if (t.size() > 0)  // have something to send, so do it.
           tmps[0] = *(end - 1);
         if (comm.rank() < (comm.size() - 1)) send_counts[target_rank] = 1;
         std::vector<ValueType> recv = ::mxx::all2allv(tmps, send_counts, comm);

         ValueType tmp = (recv.size() > 0) ? recv.front() : 0;
         //== if rank zero, and t.start is after inMemRange, then get that character.  else assume it's \n
         if (comm.rank() == 0) {
           tmp = (t.start > inMemRange.start) ? *(iter - 1) : '\n';   // if there are characters before t.start, use it.  else assume EOL
         }


         //== set up the iterator for later
         std::size_t i = t.start;

         // storage.
         std::vector<ValueType> first;
         std::vector<size_t> offsets;

         if ((tmp != ::bliss::io::BaseFileParser<Iterator>::eol) && (tmp != ::bliss::io::BaseFileParser<Iterator>::cr)) {  // previous char is not newline,
           iter = this->findEOL(iter, end, i);
         }

         // lines 1 through 4.  reorganized so that O3 optimization by gcc 5.2.1 does not skip over the whole loop (compiler bug?)
         for (int j = 0; j < 4; ++j) {
           iter = this->findNonEOL(iter, end, i);
           if (i >= t.end) break;


           offsets.push_back(i);
           first.push_back(*iter);
           iter = this->findEOL(iter, end, i);
         }

         std::cout << " rank " << comm.rank() << " prev: " << tmp << " raw lines = " << offsets.size() << "|";
         for (uint8_t j = 0; j < offsets.size(); ++j) {
           std::cout << offsets[j] << ":" << first[j] << ", ";
         }
         std::cout << std::endl;


         // find @ if present.
         size_t at_cnt = std::count(first.begin(), first.end(), '@');

//#ifdef USE_MPI
         // set rank to send to.  then use exscan to propagate it to adjacent ranks without @
         target_rank = (at_cnt == 0) ? 0 : comm.rank();
         target_rank = ::mxx::exscan(target_rank, [](int const & x, int const & y){
           return ::std::max(x, y);
         }, comm);

         // send the lines to target rank. (to the left)
         send_counts.clear();
         send_counts.resize(comm.size(), 0);
         recv.clear();
         if (comm.rank() > 0) send_counts[target_rank] = first.size();
         recv = ::mxx::all2allv(first, send_counts, comm);

         // concatenate recv to first.
         first.insert(first.end(), recv.begin(), recv.end());
//#endif

         // now that what needs to be sent has been sent, can stop if no more work is needed.
         if (at_cnt == 0) return searchRange.end;  // this has to happen AFTER comm...



         std::cout << " rank " << comm.rank() << " lines = " << offsets.size() << "|";
         for (auto ii : first) {
           std::cout << ii << ",";
         }
         std::cout << std::endl;

         //=== determine the position of a read by looking for @...+ or +...@
         // at this point, first[0] is pointing to first char after first newline, or in the case of first block, the first char.
         // everything in "first" array are first characters right after a newline.
         if (first[0] == '@' && first[2] == '+')  return (offsets.size() > 0) ? offsets[0] : searchRange.end;
         if (first[1] == '@' && first[3] == '+')  return (offsets.size() > 1) ? offsets[1] : searchRange.end;
         if (first[0] == '+' && first[2] == '@')  return (offsets.size() > 2) ? offsets[2] : searchRange.end;
         if (first[1] == '+' && first[3] == '@')  return (offsets.size() > 3) ? offsets[3] : searchRange.end;

         //=== If nothing was found, because we are at the end of the parentRange, then return the end.
         if (t.end == parentRange.end) {
           return t.end;
         }

         //=== Finally, if nothing was found, and we are not at the parentRange's end, then
         // the search range does not include a complete record.  this is an error and
         // exception is thrown.
         std::stringstream ss;
         ss << "ERROR in MPI file processing: file segment \n" << "\t\t"
             << t << "\n\t\t(original " << searchRange << ")\n"
             << "\t\t(actual search " << t << ")\n"
             << "\t\tdoes not contain valid FASTQ markers.\n String:"
             << "first chars are " << first[0] << "," << first[1] << "," << first[2] << "," << first[3];
         std::ostream_iterator<typename std::iterator_traits<Iterator>::value_type> oit(ss);
         Iterator s(_data);
         std::advance(s, (t.start - inMemRange.start));
         Iterator e(_data);
         std::advance(e, (t.end - inMemRange.start));
         std::copy(s, e, oit);

         throw ::std::logic_error(ss.str());
         }


       virtual ::std::pair<size_t, size_t> get_record_size(const Iterator &_data, const RangeType &parentRange, const RangeType &inMemRange, const RangeType &validRange, mxx::comm const & comm, size_t const count) {

         if (count == 0) throw std::invalid_argument("ERROR: called FASTQParser get_record_size with count == 0");
         if (!(parentRange.contains(inMemRange))) throw std::invalid_argument("ERROR: parentRange does not container inMemRange");
         if (!(inMemRange.contains(validRange)))  throw std::invalid_argument("ERROR: inMemRange does not container validRange");

         // find the first record that starts within search range
         size_t local_offset = this->find_first_record(_data, parentRange, inMemRange, validRange, comm);
         Iterator local_data = _data + local_offset - inMemRange.start;
         Iterator local_end = _data + validRange.end - inMemRange.start;


         // now initialize the search.  no need to call the first increment separately since we know that we have to start from seq id 0.
         size_t seq_data_len = 0;
         size_t record_size = 0;

         size_t max = (count + comm.size() - 1) / comm.size();
         size_t i = 0;
         for (; (i < max) && (local_data != local_end); ++i) {
           auto seq = this->get_next_record(local_data, local_end, local_offset);

           if (seq.record_size == 0) {  // nothing found, time to stop.
             break;
           }

           record_size += seq.record_size;
           seq_data_len += std::distance(seq.seq_begin, seq.seq_end);
         }

         i = mxx::allreduce(i, comm);  // total count

         // aggregate.
         seq_data_len = mxx::allreduce(seq_data_len, comm) / i;
         record_size = mxx::allreduce(record_size, comm) / i;


         if (record_size == 0) throw std::logic_error("ERROR: estimated sequence data size is 0");
         if (seq_data_len == 0) throw std::logic_error("ERROR: estimated record size is 0");

         return std::make_pair(record_size, seq_data_len);

       }

       virtual bool should_parse(const RangeType & range, const mxx::comm & comm) {
         return mxx::any_of(range.size() > 0, comm);
       }

#endif

    };




//    /**
//     * @class FASTQLoader
//     * @brief FileLoader subclass specialized for the FASTQ file format, uses CRTP to enforce interface consistency.
//     * @details   FASTQLoader understands the FASTQ file format, and enforces that the
//     *            L1 and L2 partition boundaries occur at FASTQ sequence record boundaries.
//     *
//     *            FASTQ files allow storage of per-base quality score and is typically used
//     *              to store a collection of short sequences.
//     *            Long sequences (such as complete genome) can be found in FASTQ formatted files
//     *              these CAN be read using standard FileLoader, provided an appropriate overlap
//     *              is specified (e.g. for kmer reading, overlap should be k-1).
//     *
//     *          for reads, expected lengths are maybe up to 1K in length, and files contain >> 1 seq
//     *
//     *          Majority of the interface is inherited from FileLoader.
//     *          This class modifies the behaviors of GetNextL1Block and getNextL2Block,
//     *            as described in FileLoader's "Advanced Usage".  Specifically, the range is modified
//     *            for each block by calling "find_first_record" to align the start and end to FASTQ sequence
//     *            record boundaries.
//     *
//     * @note    Processes 1 file at a time.  For paired end reads, a separate subclass is appropriate.
//     *          Sequences are identified by the following set of attributes:
//     *
//     *             file id.          (via filename to id map)
//     *             read id in file   (e.g. using the offset at whcih read occurs in the file as id)
//     *             position in read  (useful for kmer)
//     *
//     *          using the offset as id allows us to avoid d assignment via prefix sum.
//     *
//     * @tparam  T                 type of each element read from file
//     * @tparam  L1Buffering       bool indicating if L1 partition blocks should be buffered.  default to false
//     * @tparam  L2Buffering       bool indicating if L2 partition blocks should be buffered.  default to true (to avoid contention between threads)
//     * @tparam  L1PartitionerT    Type of the Level 1 Partitioner to generate the range of the file to load from disk
//     * @tparam  L2PartitionerT    L2 partitioner, default to DemandDrivenPartitioner
//     */
//    template<typename T,
//        bool L2Buffering = false,
//        bool L1Buffering = true,
//        typename L2PartitionerT = bliss::partition::DemandDrivenPartitioner<bliss::partition::range<size_t> > >
//    using DemandDrivenFASTQLoader =
//        FileLoader<T, 0, SequentialFASTQParser, L2Buffering, L1Buffering, L2PartitionerT, bliss::partition::DemandDrivenPartitioner<bliss::partition::range<size_t> > >;
//
//    template<typename T,
//        bool L2Buffering = false,
//        bool L1Buffering = true,
//        typename L2PartitionerT = bliss::partition::DemandDrivenPartitioner<bliss::partition::range<size_t> > >
//    using RoundRobinFASTQLoader =
//        FileLoader<T, 0, SequentialFASTQParser, L2Buffering, L1Buffering, L2PartitionerT, bliss::partition::CyclicPartitioner<bliss::partition::range<size_t> > >;
//
//    template<typename T,
//        bool L2Buffering = false,
//        bool L1Buffering = true,
//        typename L2PartitionerT = bliss::partition::DemandDrivenPartitioner<bliss::partition::range<size_t> > >
//    using FASTQLoader =
//        FileLoader<T, 0, FASTQParser, L2Buffering, L1Buffering, L2PartitionerT, bliss::partition::BlockPartitioner<bliss::partition::range<size_t> > >;

  } /* namespace io */
} /* namespace bliss */
#endif /* FASTQ_PARTITIONER_HPP_ */
