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
 * @file    filtered_sequence_iterator.hpp
 * @ingroup io
 * @author  Tony Pan <tpan7@gatech.edu>
 * @date    Feb 5, 2014
 * @brief   Contains an iterator for traversing a sequence data, record by record, a FASTQ format specific parser for use by the iterator,
 *          Sequence Id datatype definition, and 2 sequence types (with and without quality score iterators)
 *
 *
 */


#ifndef Filtered_SequencesIterator_HPP_
#define Filtered_SequencesIterator_HPP_

// C includes
#include "io/sequence_iterator.hpp"
#include "iterators/filter_iterator.hpp"
#include "io/fastq_loader.hpp"
#include "io/file_loader.hpp"

#include <algorithm>

namespace bliss
{

  namespace io
  {


    /**
     * @class bliss::io::FilteredSequencesIterator
     * @brief Iterator for parsing and traversing a block of data to access individual sequence records (of some file format), filtered by a user supplied predicate.
     * @details  This is a convenience class that allows sequence records to be filtered on the fly based on a user supplied predicate function
     *            The predicate can be used to detect presence of "N" character, or to keep only high quality reads.
     *            Sequences that satisfies the predicate (returns true) are passed through.
     *
     *            Predicate should have operation with input parameter <sequenceType>
     *
     * @note this iterator is not a random access or bidirectional iterator, as traversing in reverse is not implemented.
     *
     * @tparam Iterator	  Base iterator type to be parsed into sequences
     * @tparam SeqParser     Functoid type to parse data pointed by Iterator into sequence objects..
     * @tparam SequenceType  Sequence Type.  replaceable as long as it supports Quality if SeqParser requires it.
     */
    template<typename Iterator, template <typename> class SeqParser, typename Predicate>
    class FilteredSequencesIterator : public ::bliss::iterator::filter_iterator<
      Predicate,
      ::bliss::io::SequencesIterator<Iterator, SeqParser> >
    {
        // TODO: for long sequence FASTA file, this likely will NOT behave correctly.
        static_assert(!std::is_same<SeqParser<Iterator>, ::bliss::io::BaseFileParser<Iterator>>::value, "SplitSequencesIterator currently works only with FASTQ and FASTA files where sequences are defined." );

      protected:
        using SeqIterType = ::bliss::io::SequencesIterator<Iterator, SeqParser>;
        using FilterIterType = ::bliss::iterator::filter_iterator<Predicate, SeqIterType>;

      public:
        using iterator_category = typename SeqIterType::iterator_category;
        using value_type        = typename SeqIterType::value_type;
        using difference_type   = typename SeqIterType::difference_type;
        using pointer           = typename SeqIterType::pointer;
        using reference         = typename SeqIterType::reference;

        /**
         * @brief constructor, initializes with a start and end.  this represents the start of the output iterator
         * @details   at construction, the internal _curr, _next iterators are set to the input start iterator, and _end is set to the input end.
         *            then immediately the first record is parsed, with output populated with the first entry and next pointer advanced.
         *
         * @param f       parser functoid for parsing the data.
         * @param start   beginning of the data to be parsed
         * @param end     end of the data to be parsed.
         * @param _range  the Range associated with the start and end of the source data.  coordinates relative to the file being processed.
         */
        explicit FilteredSequencesIterator(const SeqParser<Iterator> & f,
                                           Iterator start,
                                           Iterator end,
                                           const size_t &_offset,
                                           const Predicate & pred = Predicate())
            : FilterIterType(pred, SeqIterType(f, start, end, _offset), SeqIterType(end)) {};

        /**
         * @brief constructor, initializes with only the end.  this represents the end of the output iterator.
         * @details   at construction, the internal _curr, _next, and _end iterators are set to the input end iterator.
         *          this instance therefore does not traverse and only serves as marker for comparing to the traversing FilteredSequencesIterator.
         * @param end     end of the data to be parsed.
         */
        explicit FilteredSequencesIterator(Iterator end,
                                           const Predicate & pred = Predicate())
            : FilterIterType(pred, SeqIterType(end)) {}

        /**
         * @brief default copy constructor
         * @param Other   The FilteredSequencesIterator to copy from
         */
        FilteredSequencesIterator(const FilteredSequencesIterator& Other)
            : FilterIterType(Other) {}

        /**
         * @brief  default copy assignment operator
         * @param Other   The FilteredSequencesIterator to copy from
         * @return  modified copy of self.
         */
        FilteredSequencesIterator& operator=(const FilteredSequencesIterator& Other)
        {
          this->FilterIterType::operator=(Other);
          return *this;
        }

        /**
         * @brief default copy constructor
         * @param Other   The FilteredSequencesIterator to copy from
         */
        FilteredSequencesIterator(FilteredSequencesIterator && Other)
            : FilterIterType(::std::forward<FilteredSequencesIterator>(Other)) {}

        /**
         * @brief  default copy assignment operator
         * @param Other   The FilteredSequencesIterator to copy from
         * @return  modified copy of self.
         */
        FilteredSequencesIterator& operator=(FilteredSequencesIterator && Other)
        {
          this->FilterIterType::operator=(::std::forward<FilteredSequencesIterator>(Other));
          return *this;
        }

        /// default constructor deleted.
        FilteredSequencesIterator() = default;

    };


    /**
     * @class bliss::io::SequenceNPredicate
     * @brief   given a sequence, determines if it DOES NOT contain N.
     */
    struct NSequenceFilter {

        template <typename SEQ>
        bool operator()(SEQ const & x) {
          // scan through x to look for N.
          return ::std::find(x.seq_begin, x.seq_end, 'N') == x.seq_end;
        }
    };


    template <typename Iterator, template <typename> class SeqParser>
    using NFilterSequencesIterator = bliss::io::FilteredSequencesIterator<Iterator, SeqParser, bliss::io::NSequenceFilter>;

    /**
     * @class bliss::io::SplitSequencesIterator
     * @brief Iterator for parsing and traversing a block of data to access individual sequence records (of some file format), split when predicate fails.
     * @details  for each sequence, scan with predicate.  when predicate fails, split the sequence, and continue on.
     *            EFFECTIVELY BREAKS THE SEQUENCE INTO PARTS WHERE PREDICATE FAILS.
     *
     * @note this iterator is a forward iterator only (for now).
     *
     * @tparam Iterator   Base iterator type to be parsed into sequences
     * @tparam SeqParser     Functoid type to parse data pointed by Iterator into sequence objects..
     * @tparam SequenceType  Sequence Type.  replaceable as long as it supports Quality if SeqParser requires it.
     */
    template<typename Iterator, template <typename> class SeqParser, typename Predicate>
    class SplitSequencesIterator : public ::std::iterator<
      typename ::std::conditional<
            ::std::is_same<typename ::std::iterator_traits<Iterator>::iterator_category,
                           ::std::input_iterator_tag>::value,
            ::std::input_iterator_tag,
            ::std::forward_iterator_tag>::type,
      typename SeqParser<Iterator>::SequenceType,
      typename std::iterator_traits<Iterator>::difference_type
    >
    {
        static_assert(!std::is_same<SeqParser<Iterator>, ::bliss::io::BaseFileParser<Iterator>>::value, "SplitSequencesIterator currently works only with FASTQ and FASTA files where sequences are defined." );


      protected:
        /// predicate function
        Predicate pred;

        /// internal sequence iterator type
        using SeqIterType = ::bliss::io::SequencesIterator<Iterator, SeqParser>;

        /// current position
        SeqIterType _curr;

        /// max (end) position
        SeqIterType _end;

        /// cache of the current result to return
        typename SeqParser<Iterator>::SequenceType seq;

        /// cache of the next result to check
        typename SeqParser<Iterator>::SequenceType next;

        /// get next candidate sequence.
        void get_next() {
          // repeatedly get a valid seq.  seq should not be "empty", and _curr should not be at end.
          while ((_curr != _end) && (next.seq_begin == next.seq_end)) {
            // move to next position.
            ++_curr;

            // get the next sequence if curr is not at end.
            if (_curr != _end) next = *_curr;
          }
        }

        /// splits a sequence based on predicate on chars.
        void split_seq() {
          // first copy
          // leave the quality score as is.  traversal is confined by seq_begin/seq_end.
          seq = next;
//          std::cout << "RAW SEQ " << seq << " len " << std::distance(seq.seq_begin, seq.seq_end) << std::endl;
//          std::cout << "RAW NEXT " << next << " len " << std::distance(next.seq_begin, next.seq_end) << std::endl;

          // find the beginning of valid.
          seq.seq_begin = std::find_if(next.seq_begin, next.seq_end, pred);
//          std::cout << "first SEQ " << seq << " len " << std::distance(seq.seq_begin, seq.seq_end) << std::endl;

          // find the end of the valid range
          seq.seq_end = std::find_if_not(seq.seq_begin, next.seq_end, pred);
//          std::cout << "second SEQ " << seq << " len " << std::distance(seq.seq_begin, seq.seq_end) << std::endl;

          // now update the next seq object.
          next.seq_begin_offset += std::distance(next.seq_begin, seq.seq_end);
          next.seq_begin = seq.seq_end;

//          std::cout << "SEQ " << seq << " len " << std::distance(seq.seq_begin, seq.seq_end) << std::endl;
//          std::cout << "NEXT " << next << " len " << std::distance(next.seq_begin, next.seq_end) << std::endl;

        }


      public:
        using iterator_category = typename SeqIterType::iterator_category;
        using value_type        = typename SeqIterType::value_type;
        using difference_type   = typename SeqIterType::difference_type;
        using pointer           = typename SeqIterType::pointer;
        using reference         = typename SeqIterType::reference;

        /**
         * @brief constructor, initializes with a start and end.  this represents the start of the output iterator
         * @details   at construction, the internal _curr, _next iterators are set to the input start iterator, and _end is set to the input end.
         *            then immediately the first record is parsed, with output populated with the first entry and next pointer advanced.
         *
         * @param f       parser functoid for parsing the data.
         * @param start   beginning of the data to be parsed
         * @param end     end of the data to be parsed.
         * @param _range  the Range associated with the start and end of the source data.  coordinates relative to the file being processed.
         */
        explicit SplitSequencesIterator(const SeqParser<Iterator> & f,
                                        Iterator start,
                       Iterator end, const size_t &_offset,
                       const Predicate & _pred = Predicate())
            : pred(_pred), _curr(f, start, end, _offset), _end(end) {

          // first initialize one.
          if (_curr != _end) next = *_curr;

          // this would not actually increment unless next is empty, but it will split seq.
          this->operator++();
        };

        /**
         * @brief constructor, initializes with only the end.  this represents the end of the output iterator.
         * @details   at construction, the internal _curr, _next, and _end iterators are set to the input end iterator.
         *          this instance therefore does not traverse and only serves as marker for comparing to the traversing SplitSequencesIterator.
         * @param end     end of the data to be parsed.
         */
        SplitSequencesIterator(Iterator end,
                               const Predicate & _pred = Predicate())
            : pred(_pred), _curr(end), _end(end) {}

        /**
         * @brief default copy constructor
         * @param Other   The SplitSequencesIterator to copy from
         */
        explicit SplitSequencesIterator(const SplitSequencesIterator& Other)
            : pred(Other.pred), _curr(Other._curr), _end(Other._end),
              seq(Other.seq), next(Other.next) {}

        /**
         * @brief  default copy assignment operator
         * @param Other   The SplitSequencesIterator to copy from
         * @return  modified copy of self.
         */
        SplitSequencesIterator& operator=(const SplitSequencesIterator& Other)
        {
          pred  = Other.pred;
          _curr = Other._curr;
          _end  = Other._end;
          seq   = Other.seq;
          next  = Other.next;
          return *this;
        }

        /**
         * @brief default copy constructor
         * @param Other   The SplitSequencesIterator to copy from
         */
        explicit SplitSequencesIterator(SplitSequencesIterator && Other)
            : pred(::std::move(Other.pred)),
              _curr(::std::move(Other._curr)), _end(::std::move(Other._end)),
              seq(::std::move(Other.seq)), next(::std::move(Other.next))  {}

        /**
         * @brief  default copy assignment operator
         * @param Other   The SplitSequencesIterator to copy from
         * @return  modified copy of self.
         */
        SplitSequencesIterator& operator=(SplitSequencesIterator && Other)
        {
          pred  = std::move(Other.pred );
          _curr = std::move(Other._curr);
          _end  = std::move(Other._end );
          seq   = std::move(Other.seq  );
          next  = std::move(Other.next );
          return *this;
        }
        /// default constructor deleted.
        SplitSequencesIterator() = default;

        /**
         * @brief   iterator's pre increment operation
         * @details to increment, the class uses the parser functoid to parse and traverse the input iterator
         *          At the end of the traversal, a FASTQ sequence record has been made and is kept in the parser functoid
         *          as state variable, and the input iterator has advanced to where parsing would start for the next record.
         *          The corresponding starting of the current record is saved as _curr iterator.
         *
         *          conversely, to get the starting position of the next record, we need to parse the current record.
         *
         * @note    for end iterator, because _curr and _end are both pointing to the input's end, "increment" does not do anything.
         * @return  self, updated with internal position.
         */
        SplitSequencesIterator &operator++()
        {
          // TODO get the next, and apply the predicate.
          this->get_next();

          // split if not empty
          if ((_curr != _end) && (next.seq_begin != next.seq_end)) this->split_seq();

          //== states updated.  return self.
          return *this;
        }

        /**
         * @brief   iterator's post increment operation
         * @details to increment, the class uses the parser functoid to parse and traverse the input iterator
         *          At the end of the traversal, a FASTQ sequence record has been made and is kept in the parser functoid
         *          as state variable, and the input iterator has advanced to where parsing would start for the next record.
         *          The corresponding starting of the current record is saved as _curr iterator.
         *
         *          conversely, to get the starting position of the next record, we need to parse the current record.
         *
         * @note    for end iterator, because _curr and _end are both pointing to the input's end, "increment" does not do anything.
         * @param   unnamed
         * @return  copy of self incremented.
         */
        SplitSequencesIterator operator++(int)
        {
          SplitSequencesIterator output(*this);
          this->operator++();
          return output;
        }

        /**
         * @brief     comparison operator.
         * @details   compares the underlying base iterator positions.
         * @param rhs the SequenceIterator to compare to.
         * @return    true if this and rhs SequenceIterators match.
         */
        bool operator==(const SplitSequencesIterator& rhs) const
        {
          return (_curr == rhs._curr) && ((_curr == _end) ||
              (seq.seq_begin == rhs.seq.seq_begin));
        }

        /**
         * @brief     comparison operator.
         * @details   compares the underlying base iterator positions.
         * @param rhs the SequenceIterator to compare to.
         * @return    true if this and rhs SequenceIterators do not match.
         */
        bool operator!=(const SplitSequencesIterator& rhs) const
        {
          return !(this->operator==(rhs));
        }

        /**
         * @brief dereference operator
         * @return    a const reference to the cached sequence object
         */
        typename SeqParser<Iterator>::SequenceType & operator*()
        {
          return seq;
        }

        /**
         * @brief pointer dereference operator
         * @return    a const reference to the cached sequence object
         */
        typename SeqParser<Iterator>::SequenceType *operator->() {
          return &seq;
        }

    };

    /**
     * @class bliss::io::NCharFilter
     * @brief   given a character, return true if it's NOT N.
     */
    struct NCharFilter {

        template <typename T>
        bool operator()(T const & x) {
          // scan through x to look for NOT N
          return (x != 'N') && (x != 'n');
        }
    };


    template <typename Iterator, template <typename> class SeqParser>
    using NSplitSequencesIterator = bliss::io::SplitSequencesIterator<Iterator, SeqParser, bliss::io::NCharFilter>;


    // TODO: filter for other characters.
    // TODO: quality score filter for sequences.

  } // iterator
} // bliss
#endif /* Filtered_SequencesIterator_HPP_ */
