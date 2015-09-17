/**
 * @file    sequence_iterator.hpp
 * @ingroup io
 * @author  Tony Pan <tpan7@gatech.edu>
 * @date    Feb 5, 2014
 * @brief   Contains an iterator for traversing a sequence data, record by record, a FASTQ format specific parser for use by the iterator,
 *          Sequence Id datatype definition, and 2 sequence types (with and without quality score iterators)
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 *
 */


#ifndef SequencesIterator_HPP_
#define SequencesIterator_HPP_

// C includes
#include <cassert>
#include <cstdint>

// C++ STL includes
#include <iterator>
#include <algorithm>
#include <sstream>
#include <type_traits>

// own includes
#include "common/sequence.hpp"
#include "io/io_exception.hpp"
#include "utils/logging.h"
#include "utils/function_traits.hpp"
#include "partition/range.hpp"

namespace bliss
{

  namespace io
  {


    /**
     * @class bliss::io::SequencesIterator
     * @brief Iterator for parsing and traversing a block of data to access individual sequence records (of some file format)
     * @details The iterator is initialized with start and end iterators pointing to a block of character data.
     *
     *    The source data should be of primitive type such as unsigned char.  The input iterators should be aligned
     *    to record boundaries.
     *
     *    The iterator uses a Parser Functoid to walk through the data one sequence record (e.g. read) at a time.  Different Parser
     *    types may be used, e.g. FASTA, FASTQ, to parse different file formats.
     *
     *    The iterator has the following characteristics.
     *      transforming  (changes data type)
     *      delimiting    (trivial search for the end of the "record" / beginning of next record.)
     *      finite        (needs to know base iterator's end.)
     *
     *    Design considerations:
     *      1. Repeated calls to operator*() returns the same output, WITHOUT reparsing the underlying data
     *      2. Repeated calls to operator++() moves the current pointer forward each time
     *      3. Comparison of 2 SequencesIterators compares the underlying input iterators at the start positions of the current record
     *      4. Iterator states:  before first use, during use, after reaching the end.
     *          before:  operator* points to the first output, no need to call ++ first.  so this is same case as "during".
     *          during:  operator* gets the most current output that's a result of calling ++, ++ moves pointer forward
     *          after:   operator* returns empty result, ++ is same as no-op.
     *     Based on these, we need to
     *      1. cache the output
     *      2. keep a current iterator  (for comparison between SequencesIterators)
     *      3. keep a next iterator     (for parsing)
     *      4. have constructor parse the first entry, so that dereference functions properly
     *      5. check for end condition  (next == end) and return if at end.
     *
     *      The Parser functoid then only has 1 simple task:  to increment, as it does not need to store locally.
     *
     * @note this iterator is not a random access or bidirectional iterator, as traversing in reverse is not implemented.
     *
     * @tparam Iterator	  Base iterator type to be parsed into sequences
     * @tparam Parser     Functoid type to parse data pointed by Iterator into sequence objects..
     * @tparam SequenceType  Sequence Type.  replaceable as long as it supports Quality if Parser requires it.
     */
    template<typename Iterator, template<typename> class Parser>
    class SequencesIterator :
    	public ::std::iterator<
            typename ::std::conditional<
                  ::std::is_same<typename ::std::iterator_traits<Iterator>::iterator_category,
                   	   	   	   	 ::std::input_iterator_tag>::value,
                  ::std::input_iterator_tag,
                  ::std::forward_iterator_tag>::type,
            typename Parser<Iterator>::SequenceType,
            typename std::iterator_traits<Iterator>::difference_type
          >
    {
      public:

        /// type of data elements in the input iterator
        typedef typename std::iterator_traits<Iterator>::value_type BaseValueType;

      protected:

        /// type of SequencesIterator class
        typedef SequencesIterator<Iterator, Parser> type;

        /// cache of the current result.
        typename Parser<Iterator>::SequenceType seq;

        /// iterator at the current starting point in the input data from where the current FASTQ record was parsed.
        Iterator _curr;

        /// iterator at the next starting point in the input data from where the next FASTQ record should be parsed.
        Iterator _next;

        /// iterator pointing to the end of the input data, not to go beyond.
        Iterator _end;

        /// instance of the parser functoid, used during record parsing.
        Parser<Iterator> parser;

        /**
         * @brief  the offset from the start of the file.  This is used as sequenceId.
         * @details  for each file with file id (fid) between 0 and 255, the seq id starts at (fid << 40).
         *          if base iterator is pointer, then units are in bytes.
         */
        size_t file_offset;

        ///  sequence index in any shared data vector.
        size_t seq_offset;


      public:
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
        SequencesIterator(const Parser<Iterator> & f, const Iterator& start,
                       const Iterator& end, const size_t &_offset)
            : seq(), _curr(start), _next(start), _end(end), parser(f), file_offset(_offset), seq_offset(0)
        {
          // parse the first entry, if start != end.
          // effect: _next incremented to next parse start point.
          //         file_offset updated to next start point,
          //         output updated with either empty (if _next at end)  or the new output
          if (_next != _end)
            seq_offset = parser.increment(_next, _end, file_offset, seq);  // first increment.
        }

        /**
         * @brief constructor, initializes with only the end.  this represents the end of the output iterator.
         * @details   at construction, the internal _curr, _next, and _end iterators are set to the input end iterator.
         *          this instance therefore does not traverse and only serves as marker for comparing to the traversing SequencesIterator.
         * @param end     end of the data to be parsed.
         */
        explicit SequencesIterator(const Iterator& end)
            : seq(), _curr(end), _next(end), _end(end), parser(), file_offset(std::numeric_limits<size_t>::max()), seq_offset(std::numeric_limits<size_t>::max())
        {
        }

        /**
         * @brief default copy constructor
         * @param Other   The SequencesIterator to copy from
         */
        explicit SequencesIterator(const type& Other)
            : seq(Other.seq), _curr(Other._curr), _next(Other._next), _end(Other._end),
              parser(Other.parser),  file_offset(Other.file_offset), seq_offset(Other.seq_offset)
        {}

        /**
         * @brief  default copy assignment operator
         * @param Other   The SequencesIterator to copy from
         * @return  modified copy of self.
         */
        type& operator=(const type& Other)
        {
          seq = Other.seq;
          _curr = Other._curr;
          _next = Other._next;
          _end = Other._end;
          parser = Other.parser;
          file_offset = Other.file_offset;
          seq_offset = Other.seq_offset;
          return *this;
        }

        /**
         * @brief default copy constructor
         * @param Other   The SequencesIterator to copy from
         */
        explicit SequencesIterator(type && Other)
            : seq(std::move(Other.seq)), _curr(std::move(Other._curr)), _next(std::move(Other._next)), _end(std::move(Other._end)),
              parser(std::move(Other.parser)),  file_offset(std::move(Other.file_offset)), seq_offset(std::move(Other.seq_offset))
        {}

        /**
         * @brief  default copy assignment operator
         * @param Other   The SequencesIterator to copy from
         * @return  modified copy of self.
         */
        type& operator=(type && Other)
        {
          seq =   std::move(Other.seq);
          _curr = std::move(Other._curr);
          _next = std::move(Other._next);
          _end =  std::move(Other._end);
          parser = std::move(Other.parser);
          file_offset = std::move(Other.file_offset);
          seq_offset = std::move(Other.seq_offset);
          return *this;
        }

        /// default constructor deleted.
        SequencesIterator() = delete;


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
        type &operator++()
        {
          //== at end of the data block, can't parse any further, no movement.  (also true for end iterator)
          if (_next == _end)
          {
            //== check if reaching _end for the first time
            if (_curr != _end) {
              // advance to end
              _curr = _next;
              // and set output to empty.
              seq = typename Parser<Iterator>::SequenceType();
            }
          } else {
            //== else we can parse, so do it.

            //== move the current forward
            _curr = _next;

            //== then try parsing.  this advances _next and offset, ready for next ++ call, and updates output cache.
            parser.increment(_next, _end, file_offset, seq_offset, seq);
          }

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
        type operator++(int)
        {
          type output(*this);
          this->operator++();
          return output;
        }

        /**
         * @brief   accessor for the internal base iterator, pointing to the current position
         * @return  base iterator
         */
        Iterator& getBaseIterator()
        {
          return _curr;
        }
        /**
         * @brief   const accessor for the internal base iterator, pointing to the current position
         * @return  const base iterator
         */
        const Iterator& getBaseIterator() const
        {
          return _curr;
        }

        /**
         * @brief     comparison operator.
         * @details   compares the underlying base iterator positions.
         * @param rhs the SequenceIterator to compare to.
         * @return    true if this and rhs SequenceIterators match.
         */
        bool operator==(const type& rhs) const
        {
          return _curr == rhs._curr;
        }

        /**
         * @brief     comparison operator.
         * @details   compares the underlying base iterator positions.
         * @param rhs the SequenceIterator to compare to.
         * @return    true if this and rhs SequenceIterators do not match.
         */
        bool operator!=(const type& rhs) const
        {
          return _curr != rhs._curr;
        }

        /**
         * @brief dereference operator
         * @return    a const reference to the cached sequence object
         */
        const typename Parser<Iterator>::SequenceType &operator*() const
        {
          return seq;
        }

        /**
         * @brief pointer dereference operator
         * @return    a const reference to the cached sequence object
         */
        const typename Parser<Iterator>::SequenceType *operator->() const {
          return &seq;
        }


    };

  } // iterator
} // bliss
#endif /* SequencesIterator_HPP_ */
