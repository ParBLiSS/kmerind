/**
 * @file    concatenating_iterator.hpp
 * @ingroup iterators
 * @author  Tony Pan <tpan7@gatech.edu>
 * @brief   contains an iterator that concatenates multiple ranges represented as iterators.
 * @details
 *
 * Copyright (c) 2015 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef CONCATENATING_ITERATOR_HPP_
#define CONCATENATING_ITERATOR_HPP_

#include <vector>
#include <iterator>
#include <utility>
#include "utils/logging.h"
#include <type_traits>  // for testing const iterator

namespace bliss
{
  namespace iterator
  {

    /**
     * @class    bliss::iterator::ConcatenatingIterator
     * @brief    this class presents a single/sequential view of a series of underlying iterator ranges
     * @details  random access iterator is not supported. other iterator categories are okay.
     *
     */
    template<typename Iterator>
    class ConcatenatingIterator :
        public std::iterator<
        typename std::conditional<
        std::is_same<typename std::iterator_traits<Iterator>::iterator_category,
        std::random_access_iterator_tag>::value,
        std::bidirectional_iterator_tag,
        typename std::iterator_traits<Iterator>::iterator_category>::type,
        typename std::iterator_traits<Iterator>::value_type
        >
    {
      protected:
        /// Define iterator trait of the base iterator
        typedef std::iterator_traits<Iterator> base_traits;

        /// Define pair of iterators as start and end of a range
        typedef std::pair<Iterator, Iterator> RangeType;

        /// list of iterator ranges to be concatenated, in order
        std::vector<RangeType > ranges;

        /// the current position in the ranges list
        int64_t curr_iter_pos;

        /// the current iterator position in the range of interest.
        Iterator curr;

        typedef ConcatenatingIterator<Iterator> type;
      public:

        /// DEFINE iterator element's value type
        typedef typename base_traits::value_type value_type;

        /// constructor for start concatenating iterator using copy semantic
        ConcatenatingIterator(const std::vector<RangeType>& _ranges)
        : ranges(_ranges), curr_iter_pos(0)
        {
          ranges.erase(std::remove_if(ranges.begin(), ranges.end(), [](const RangeType & range){
                  return (range.first == range.second);
                }),
                ranges.end());

          curr = ranges.begin()->first;
        };

        /// constructor for start iterator, using move semantic.
        ConcatenatingIterator(std::vector<RangeType>&& _ranges)
        : ranges(std::forward<std::vector<RangeType> >(_ranges)),
          curr_iter_pos(0)
        {

          ranges.erase(std::remove_if(ranges.begin(), ranges.end(), [](const RangeType & range){
                  return (range.first == range.second);
                }),
                ranges.end());

          curr = ranges.begin()->first;
        };


        /// constructor for end iterator, copy semantic
        ConcatenatingIterator(const Iterator& end) : curr_iter_pos(-1), curr(end)
        {
        };

        /// copnstructor for end iterator, move semantic
        ConcatenatingIterator(Iterator&& end) : curr_iter_pos(-1), curr(std::move(end))
        {
        };

        /// default constructor if the base iterator is forward, bidirection, or random access.  useful for incremental concatenation.
        template<typename C = typename base_traits::iterator_category,
            typename std::enable_if<
            std::is_same<C, std::forward_iterator_tag>::value ||
            std::is_same<C, std::bidirectional_iterator_tag>::value ||
            std::is_same<C, std::random_access_iterator_tag>::value,
            int>::type = 0
            >
        ConcatenatingIterator() : curr_iter_pos(-1), curr()
        {
        };

        // note that explicit keyword cannot be on copy and move constructors else the constructors are not defined/found.

        /// copy constructor.  respects multi-pass requirement of forward iterator.
        ConcatenatingIterator(const type& other)
        : ranges(other.ranges), curr_iter_pos(other.curr_iter_pos)
        {
          // move curr forward the same distance as other.curr - other.starts[curr_iter_pos]
          curr = ranges[curr_iter_pos].first;
          for (auto it = other.ranges[curr_iter_pos].first;
              it != other.curr;
              ++it, ++curr);
        };

        /// copy constructor.  respects multi-pass requirement of forward iterator.
        type& operator=(const type& other)
        {
          ranges = other.ranges;
          curr_iter_pos(other.curr_iter_pos);
          // move curr forward the same distance as other.curr - other.starts[curr_iter_pos]
          curr = ranges[curr_iter_pos].first;
          for (auto it = other.ranges[curr_iter_pos].first;
              it != other.curr;
              ++it, ++curr);

          return *this;
        }

        /// move constructor.
        ConcatenatingIterator(type&& other)
        {
          curr_iter_pos = other.curr_iter_pos;

          // move curr forward the same distance as other.curr - other.starts[curr_iter_pos]
          int i = 0;
          for (auto it = other.ranges[curr_iter_pos].first;
              it != other.curr;
              ++it, ++i);

          ranges = std::move(other.ranges);
          curr = ranges[curr_iter_pos].first;
          for (int j = 0; j < i; ++i, ++curr);

        };

        /// move assignment operator
        type& operator=(type&& other)
        {
          curr_iter_pos = other.curr_iter_pos;

          // move curr forward the same distance as other.curr - other.starts[curr_iter_pos]
          int i = 0;
          for (auto it = other.ranges[curr_iter_pos].first;
              it != other.curr;
              ++it, ++i);

          ranges = std::move(other.ranges);
          curr = ranges[curr_iter_pos].first;
          for (int j = 0; j < i; ++i, ++curr);

          return *this;
        }


        /**
         * @brief add additional iterator range to the list, via copy semantic
         */
        bool addRange(const Iterator& start, const Iterator& end) {
          return this->addRange(std::move(std::make_pair(start, end)));
        }

        /**
         * @brief add additional iterator range to the list via move semantic
         */
        bool addRange(Iterator&& start, Iterator&& end) {
          return this->addRange(std::move(std::make_pair(std::forward<Iterator>(start), std::forward<Iterator>(end))));
        }



        /**
         * @brief add additional iterator range to the list via copy semantic
         */
        bool addRange(const std::pair<Iterator, Iterator>& range) {
          if (range.first == range.second) return false;  // nothing useful to add.

          // error when final position is greater than the size of range vector
          if (curr_iter_pos >= static_cast<int64_t>(ranges.size())) {
            ERRORF("position pointing to beyond available ranges: pos %ld, size %lu", curr_iter_pos, ranges.size());
            return false;
          }

          // check if at end
          bool atEnd = this->at_end();  // eval at_end before changing size of vector.
          ranges.push_back(range);      // but first add the range

          // if iterator was at the end of the concatenation, then move curr to the start of the new range.
          if (atEnd) { 
            ++curr_iter_pos;

            curr = ranges[curr_iter_pos].first;
            //WARNINGF("Converting empty or finished ConcatenatingIterator (end iterator) to non-empty (non-end).  pos %ld", curr_iter_pos);
            // end iterator, now converted to non-end
          }
          return true;

        }

        /**
         * @brief add additional iterator range to the list via move semantic
         */
        bool addRange(std::pair<Iterator, Iterator>&& range) {
          if (range.first == range.second) return false;  // not useful to add.

          // error when final position is greater than the size of range vector
          if (curr_iter_pos >= static_cast<int64_t>(ranges.size())) {
            ERRORF("position pointing to beyond available ranges: pos %ld, size %lu", curr_iter_pos, ranges.size());
            return false;
          }

          // check if at end
          bool atEnd = this->at_end();  // eval at_end before changing size of vector.
          ranges.push_back(std::move(std::forward<std::pair<Iterator,Iterator> >(range)));  // first add the range


          // if iterator was at the end of the concatenation, then move curr to the start of the new range.
          if (atEnd) { 
            ++curr_iter_pos;

            curr = ranges[curr_iter_pos].first;
            WARNINGF("Converting empty or finished ConcatenatingIterator (end iterator) to non-empty (non-end).  pos %ld", curr_iter_pos);
            // end iterator, now converted to non-end
          }
          return true;
        }



        /**
         * @brief increment:  move to the next position in the concatenating iterator, which may cross range boundaries.
         * @return
         */
        type& operator++()
            {
          //=== first set of conditions are when we are at the end.
          if (this->at_end()) return *this;

          // if not at end, 2 possibilities:
          // 1, was at very end, but a new range was added.
          // 2. in middle of a start-end range.
          // 3. moved onto end of a range. subcase, moved onto last end.

          // first case, curr would be at very end, so move it to next
          // if at an iterator boundary,  go to next start iterator (
          if (curr == ranges[curr_iter_pos].second) {
            ++curr_iter_pos;
            curr = ranges[curr_iter_pos].first;
            return *this;
          }

          // second case
          ++curr;

          // and third case.
          if (curr == ranges[curr_iter_pos].second) {
            if (this->at_end()) return *this;  // could be moved into last entry

            // else go to next range.
            ++curr_iter_pos;
            curr = ranges[curr_iter_pos].first;
          }
          return *this;
        }

        /**
         * post increment.  make a copy then increment that.
         */
        type operator++(int)
            {
          type output(*this);
          this->operator++();
          return output;
            }

        //=== input iterator specific

        /// comparison operator
        inline bool operator==(const type& rhs) const
            {
          if (this->at_end() && rhs.at_end()) return true;

          return curr == rhs.curr;
            }

        /// comparison operator
        inline bool operator!=(const type& rhs) const
            {
          return !(this->operator==(rhs));
            }

        /// get pointer
        inline typename base_traits::pointer operator->() const
        {
          return &(*curr);
        }

        /// dereference operator.  normal (not an output iterator) - returns a rvalue
        inline value_type operator*() const
        {
          return *curr;
        }

        /// dereference operator.  output, forward, bidirectional, and random access iterators supported.  return lvalue
        template<typename P = typename base_traits::pointer>
        inline typename std::enable_if<
          !std::is_const<typename std::remove_pointer<P>::type>::value,
          typename base_traits::reference
        >::type operator*()
        {
          return *curr;
        }


        //=== bidirectional iterator

        /**
         * @brief pre decrement operator:  semantics of -- does not have a bound on the start side.
         */
        template<typename C = typename base_traits::iterator_category>
        typename std::enable_if<
          std::is_same<C, std::bidirectional_iterator_tag>::value ||
          std::is_same<C, std::random_access_iterator_tag>::value,
          type
        >::type& operator--()
        {
          if (at_beginning()) return *this;

          // not at beginning.  a few cases:
          // 1. at beginning of a range.  go to previous
          if (curr == ranges[curr_iter_pos].first) {
            --curr_iter_pos;
            curr = ranges[curr_iter_pos].second;
            --curr;
            return *this;
          }

          // 2. not at beginning of a range.
          --curr;
          return *this;
        }

        /**
         * @brief post decrement operator:  semantics of -- does not have a bound on the start side.
         */
        template<typename C = typename base_traits::iterator_category>
        typename std::enable_if<
          std::is_same<C, std::bidirectional_iterator_tag>::value ||
          std::is_same<C, std::random_access_iterator_tag>::value,
          type
        >::type operator--(int)
        {
          type output(*this);
          this->operator--();
          return output;
        }

      protected:

        /**
         * @brief   check if current iterator is at the end of all ranges (if an end iterator, or end of start iterator, etc.)
         * @return
         */
        bool at_end() const
        {
          // end iterator
          if (curr_iter_pos < 0) return true;
          // start iterator ended
          if (curr_iter_pos >= static_cast<int64_t>(ranges.size())) {
            ERRORF("ERROR: concat iterator at one position past end of iterator vector");
            return true;
          }
          // last iterator and last entry in iterator, finished as well
          if (curr_iter_pos == (static_cast<int64_t>(ranges.size()) - 1) && curr == ranges[curr_iter_pos].second) return true;

          return false;
        }

        /**
         * @brief   check if current iterator is at the beginning of all ranges (if an end iterator, or end of start iterator, etc.)
         * @return
         */
        bool at_beginning() const
        {
          // end iterator
          if (curr_iter_pos < 0) return true;
          // last iterator and last entry in iterator, finished as well
          if (curr_iter_pos == 0 && curr == ranges[0].first) return true;

          return false;
        }

    };

  } /* namespace iterator */
} /* namespace bliss */

#endif /* CONCATENATING_ITERATOR_HPP_ */
