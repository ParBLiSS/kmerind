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
     * @tparam Parser     Functoid type to parse data pointed by Iterator into sequence objects..
     * @tparam SequenceType  Sequence Type.  replaceable as long as it supports Quality if Parser requires it.
     */
    template<typename Iterator, template <typename> class Parser, typename Predicate>
    class FilteredSequencesIterator : public ::bliss::iterator::filter_iterator<
      Predicate,
      ::bliss::io::SequencesIterator<Iterator, Parser> >
    {
      protected:
        using SeqIterType = ::bliss::io::SequencesIterator<Iterator, Parser>;
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
        explicit FilteredSequencesIterator(const Parser<Iterator> & f, const Iterator& start,
                       const Iterator& end, const size_t &_offset,
                       const Predicate & pred = Predicate())
            : FilterIterType(pred, SeqIterType(f, start, end, _offset), SeqIterType(end)) {};

        /**
         * @brief constructor, initializes with only the end.  this represents the end of the output iterator.
         * @details   at construction, the internal _curr, _next, and _end iterators are set to the input end iterator.
         *          this instance therefore does not traverse and only serves as marker for comparing to the traversing FilteredSequencesIterator.
         * @param end     end of the data to be parsed.
         */
        explicit FilteredSequencesIterator(const Iterator& end,
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
        FilteredSequencesIterator() = delete;

//        using operator++;
//        using operator==;
//        using operator*;
//        using operator->;
//        using operator();

    };

    /**
     * @class bliss::io::SequenceNPredicate
     * @brief   given a sequence, determines if it DOES NOT contain N.
     */
    struct SequenceNPredicate {

        template <typename SEQ>
        bool operator()(SEQ const & x) {
          // scan through x to look for N.
          return ::std::find(x.seq_begin, x.seq_end, 'N') == x.seq_end;
        }
    };


    template <typename Iterator, template <typename> class Parser>
    using NFilterSequencesIterator = bliss::io::FilteredSequencesIterator<Iterator, Parser, bliss::io::SequenceNPredicate>;

    // TODO: filter for other characters.
    // TODO: quality score filter for sequences.

  } // iterator
} // bliss
#endif /* Filtered_SequencesIterator_HPP_ */
