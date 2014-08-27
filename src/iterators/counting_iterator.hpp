/**
 * @file    counting_iterator.hpp
 * @ingroup bliss::iterators
 * @author  tpan
 * @brief   contains a counting iterator, which iterates over a range from start to end.
 * @details
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef COUNTING_ITERATOR_HPP_
#define COUNTING_ITERATOR_HPP_

#include <type_traits>
#include <iterator>

#include "partition/range.hpp"


namespace bliss
{
  namespace iterator
  {

    /**
     * @class    bliss::iterator::CountingIterator
     * @brief    iterator that iterates from start to end in specified steps size
     * @details  counting iterator counts from the start to the end in regular steps.
     *           this is useful to establish index values for some other array.
     * @tparam RangeType  the type of range to iterate over.
     */
    template<typename RangeType>
    class CountingIterator :  public std::iterator<std::random_access_iterator_tag,
                                                   typename RangeType::ValueType,
                                                   typename std::conditional<std::is_integral<typename RangeType::ValueType>::value,
                                                     int64_t, typename RangeType::ValueType>::type>
    {
      protected:
        /// output value type
        using T = typename std::iterator_traits<CountingIterator<RangeType> >::value_type;

        /// difference type
        using D = typename std::iterator_traits<CountingIterator<RangeType> >::difference_type;

        /// the range within which to iterate over
        const RangeType range;

        /// the stride for each iteration.
        const T stride;

        /// current position in range;
        const T val;

      public:

        /**
         * constructor
         * @param _range    range to iterate over
         * @param _stride     the distance traversed during each call to the increment/decrement method.
         */
        CountingIterator(const RangeType& _range, const T &_stride) : range(_range), stride(_stride), val(_range.start) {};

        /**
         * constructor, with stride defaults to 1.
         * @param _range    range to iterate over
         */
        CountingIterator(const RangeType& _range) : range(_range), stride(1), val(_range.start) {};

        /**
         * default constructor, sets range to empty, stride to 1
         */
        CountingIterator() : range(), stride(1), val(0) {};

        /**
         * default copy constructor
         * @param other  instance of CountingIterator to copy from
         */
        CountingIterator(const CountingIterator<RangeType> & other) = default;

        /**
         * default copy assignment operator
         * @param other  instance of CountingIterator to copy from
         * @return reference to self
         */
        CountingIterator<RangeType>& operator=(const CountingIterator<RangeType> & other) = default;

        /**
         * default move constructor
         * @param other  instance of CountingIterator to move from
         */
        CountingIterator(CountingIterator<RangeType> && other) = default;

        /**
         * default move assignment operator
         * @param other  instance of CountingIterator to move from
         * @return reference to self
         */
        CountingIterator<RangeType>& operator=(CountingIterator<RangeType> && other) = default;

        /**
         * default destructor
         */
        virtual ~CountingIterator();


        /**
         * @brief   pre-increment operator: ++iter
         * @return  reference to incremented iterator
         */
        CountingIterator<RangeType>& operator++() {
          val += stride;
          return *this;
        }

        /**
         * @brief   post-increment operator: iter++
         * @param   dummy for c++ to identify this as post increment.
         * @return  incremented copy of iterator
         */
        CountingIterator<RangeType> operator++(int) {
          CountingIterator<RangeType> out(*this);
          ++out;
          return out;
        }

        /**
         * @brief compare to other iterator for equality.  only inspect the val (position in range)
         * @param other   iterator to compare to.
         * @return  bool, true if equal, false otherwise.
         */
        bool operator==(const CountingIterator<RangeType>& other) {
          return val == other.val;
        }

        /**
         * @brief compare to other iterator for inequality
         * @param other   iterator to compare to.
         * @return  bool, true if not equal, false otherwise.
         */
        bool operator!=(const CountingIterator<RangeType>& other) {
          return !(this->operator==(other));
        }

        /**
         * @brief dereference function, *iter
         * @return  current value (position in range)
         */
        const T operator*() const {
          return val;
        }

        /**
         * @brief pointer access function, iter->val
         * @return  current value (position in range)
         */
        const T* operator->() const {
          return &val;
        }

        /**
         * @brief   pre-decrement operator: --iter
         * @return  reference to decremented iterator
         */
        CountingIterator<RangeType>& operator--() {
          val -= stride;
          return *this;
        }

        /**
         * @brief   post-decrement operator: iter--
         * @param   dummy for c++ to identify this as post decrement.
         * @return  decremented copy of iterator
         */
        CountingIterator<RangeType> operator--(int) {
          CountingIterator<RangeType> out(*this);
          --out;
          return out;
        }

        /**
         * @brief   arithmetic increment by multiple steps
         * @param diff    number of steps to increment by
         * @return  incremented copy of iterator
         */
        CountingIterator<RangeType> operator+(const D& diff) {
          CountingIterator<RangeType> out(*this);
          out += diff;
          return out;
        }


        /**
         * @brief   arithmetic decrement by multiple steps
         * @param diff    number of steps to decrement by
         * @return  decremented copy of iterator
         */
        CountingIterator<RangeType> operator-(const D& diff) {
          CountingIterator<RangeType> out(*this);
          out -= diff;
          return out;
        }



        /**
         * @brief compare to other iterator: >
         * @param other   iterator to compare to.
         * @return  bool, true if greater than, false otherwise.
         */
        bool operator>(const CountingIterator<RangeType>& other) {
          return val > other.val;
        }

        /**
         * @brief compare to other iterator: <
         * @param other   iterator to compare to.
         * @return  bool, true if less than, false otherwise.
         */
        bool operator<(const CountingIterator<RangeType>& other) {
          return val < other.val;
        }

        /**
         * @brief compare to other iterator: >=
         * @param other   iterator to compare to.
         * @return  bool, true if greater than or equal to, false otherwise.
         */
        bool operator>=(const CountingIterator<RangeType>& other) {
          return val >= other.val;
        }

        /**
         * @brief compare to other iterator: <=
         * @param other   iterator to compare to.
         * @return  bool, true if less than or equal to, false otherwise.
         */
        bool operator<=(const CountingIterator<RangeType>& other) {
          return val <= other.val;
        }

        /**
         * @brief   arithmetic increment by multiple steps
         * @param diff    number of steps to increment by
         * @return  reference to incremented iterator
         */
        CountingIterator<RangeType>& operator+=(const D& diff) {
          val += diff * stride;
          return *this;
        }

        /**
         * @brief   arithmetic decrement by multiple steps
         * @param diff    number of steps to decrement by
         * @return  reference to decremented iterator
         */
        CountingIterator<RangeType>& operator-=(const D& diff) {
          val -= diff * stride;
          return *this;
        }

        /**
         * @brief   offset dereference operator
         * @param i offset at which the value is retrieved.
         * @return  value (position in range) for ith offset
         */
        T operator[](const D& i) {
          return range.start + stride * i;
        }

    };

    /**
     * @brief increment operator, with first operand being a number and second being an iterator  n + iter;
     * @tparam RangeType    Range type for which to perform the increment operation
     * @param diff          number of steps to increment by
     * @param self          iterator to increment
     * @return              copy of incremented iterator
     */
    template<typename RangeType>
    CountingIterator<RangeType> operator+(const typename std::iterator_traits<CountingIterator<RangeType> >::difference_type& diff, CountingIterator<RangeType>& self) {
      return self + diff;
    }


  } /* namespace iterator */
} /* namespace bliss */

#endif /* COUNTING_ITERATOR_HPP_ */
