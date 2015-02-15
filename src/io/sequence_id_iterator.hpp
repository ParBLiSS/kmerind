/**
 * @file    sequence_id_iterator.hpp
 * @ingroup bliss::iterators
 * @author  Tony Pan <tpan7@gatech.edu>
 * @brief   contains a sequence_id iterator
 * @details
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef SEQUENCE_ID_ITERATOR_HPP_
#define SEQUENCE_ID_ITERATOR_HPP_

#include <cstddef>

#include <type_traits>
#include <iterator>
#include <utility>

namespace bliss
{
  namespace iterator
  {

    /**
     * @class    bliss::iterator::SequenceIdIterator
     * @brief    iterator that iterates from start to end in specified steps size
     * @details  sequence_id iterator counts from the start to the end in regular steps.
     *           this is useful to establish index values for some other array.
     * @tparam SequenceIdType the type of data to over.
     */
    template<typename SequenceIdType>
    class SequenceIdIterator :  public std::iterator<std::random_access_iterator_tag, SequenceIdType>
    {
      protected:

        /// difference type
        using D = std::ptrdiff_t;

        /// position type
        using T = decltype(std::declval<SequenceIdType>().pos);

        /// the stride for each iteration.
        mutable T stride;

        /// current value;
        SequenceIdType val;

      public:

        /**
         * constructor
         * @param _start      first value
         * @param _stride     the distance traversed during each call to the increment/decrement method.
         */
        SequenceIdIterator(const SequenceIdType& _start, const T &_stride) : stride(_stride), val(_start) {};

        /**
         * constructor, with stride defaults to 1.
         * @param _start    first value
         */
        SequenceIdIterator(const SequenceIdType& _start) : stride(1), val(_start) {};

        /**
         * default constructor, sets start to 0, stride to 1
         */
        SequenceIdIterator() : stride(1), val(0) {};

        /**
         * default copy constructor
         * @param other  instance of SequenceIdIterator to copy from
         */
        SequenceIdIterator(const SequenceIdIterator<SequenceIdType> & other) : stride(other.stride), val(other.val) {};

        /**
         * default copy assignment operator
         * @param other  instance of SequenceIdIterator to copy from
         * @return reference to self
         */
        SequenceIdIterator<SequenceIdType>& operator=(const SequenceIdIterator<SequenceIdType> & other) {
          stride = other.stride;
          val = other.val;
          return *this;
        }

        /**
         * default move constructor
         * @param other  instance of SequenceIdIterator to move from
         */
        SequenceIdIterator(SequenceIdIterator<SequenceIdType> && other) : stride(other.stride), val(std::move(other.val)) {
          other.stride = 1;
        };

        /**
         * default move assignment operator
         * @param other  instance of SequenceIdIterator to move from
         * @return reference to self
         */
        SequenceIdIterator<SequenceIdType>& operator=(SequenceIdIterator<SequenceIdType> && other) {
          stride = other.stride; other.stride = 1;
          val = std::move(other.val);
          return *this;
        };

        /**
         * default destructor
         */
        virtual ~SequenceIdIterator() {};


        /**
         * @brief   pre-increment operator: ++iter
         * @return  reference to incremented iterator
         */
        SequenceIdIterator<SequenceIdType>& operator++() {
          val.pos += stride;
          return *this;
        }

        /**
         * @brief   post-increment operator: iter++
         * @param   dummy for c++ to identify this as post increment.
         * @return  incremented copy of iterator
         */
        SequenceIdIterator<SequenceIdType> operator++(int) {
          SequenceIdIterator<SequenceIdType> out(*this);
          ++out;
          return out;
        }

        /**
         * @brief compare to other iterator for equality.  only inspect the val
         * @param other   iterator to compare to.
         * @return  bool, true if equal, false otherwise.
         */
        bool operator==(const SequenceIdIterator<SequenceIdType>& other) {
          return val == other.val;
        }

        /**
         * @brief compare to other iterator for inequality
         * @param other   iterator to compare to.
         * @return  bool, true if not equal, false otherwise.
         */
        bool operator!=(const SequenceIdIterator<SequenceIdType>& other) {
          return !(this->operator==(other));
        }

        /**
         * @brief dereference function, *iter
         * @return  current value
         */
        const SequenceIdType operator*() const {
          return val;
        }

        /**
         * @brief   pre-decrement operator: --iter
         * @return  reference to decremented iterator
         */
        SequenceIdIterator<SequenceIdType>& operator--() {
          val.pos -= stride;
          return *this;
        }

        /**
         * @brief   post-decrement operator: iter--
         * @param   dummy for c++ to identify this as post decrement.
         * @return  decremented copy of iterator
         */
        SequenceIdIterator<SequenceIdType> operator--(int) {
          SequenceIdIterator<SequenceIdType> out(*this);
          --out;
          return out;
        }

        /**
         * @brief   arithmetic increment by multiple steps
         * @param diff    number of steps to increment by
         * @return  incremented copy of iterator
         */
        SequenceIdIterator<SequenceIdType> operator+(const D& diff) {
          SequenceIdIterator<SequenceIdType> out(*this);
          out += diff;
          return out;
        }


        /**
         * @brief   arithmetic decrement by multiple steps
         * @param diff    number of steps to decrement by
         * @return  decremented copy of iterator
         */
        SequenceIdIterator<SequenceIdType> operator-(const D& diff) {
          SequenceIdIterator<SequenceIdType> out(*this);
          out -= diff;
          return out;
        }

        /**
         * @brief difference between 2 iterators;
         * @param other         the iterator to subtract by
         * @return              distance between the iterators
         */
        D operator-(const SequenceIdIterator<SequenceIdType>& other) {
          return (val.pos - other.val.pos) / stride;
        }


        /**
         * @brief compare to other iterator: >
         * @param other   iterator to compare to.
         * @return  bool, true if greater than, false otherwise.
         */
        bool operator>(const SequenceIdIterator<SequenceIdType>& other) {
          return val > other.val;
        }

        /**
         * @brief compare to other iterator: <
         * @param other   iterator to compare to.
         * @return  bool, true if less than, false otherwise.
         */
        bool operator<(const SequenceIdIterator<SequenceIdType>& other) {
          return val < other.val;
        }

        /**
         * @brief compare to other iterator: >=
         * @param other   iterator to compare to.
         * @return  bool, true if greater than or equal to, false otherwise.
         */
        bool operator>=(const SequenceIdIterator<SequenceIdType>& other) {
          return val >= other.val;
        }

        /**
         * @brief compare to other iterator: <=
         * @param other   iterator to compare to.
         * @return  bool, true if less than or equal to, false otherwise.
         */
        bool operator<=(const SequenceIdIterator<SequenceIdType>& other) {
          return val <= other.val;
        }

        /**
         * @brief   arithmetic increment by multiple steps
         * @param diff    number of steps to increment by
         * @return  reference to incremented iterator
         */
        SequenceIdIterator<SequenceIdType>& operator+=(const D& diff) {
          val.pos += diff * stride;
          return *this;
        }

        /**
         * @brief   arithmetic decrement by multiple steps
         * @param diff    number of steps to decrement by
         * @return  reference to decremented iterator
         */
        SequenceIdIterator<SequenceIdType>& operator-=(const D& diff) {
          val.pos -= diff * stride;
          return *this;
        }

        /**
         * @brief   offset dereference operator
         * @param i offset at which the value is retrieved.
         * @return  value for ith offset
         */
        SequenceIdType operator[](const D& i) {
          return val.pos + stride * i;
        }

    };

    /**
     * @brief increment operator, with first operand being a number and second being an iterator  n + iter;
     * @tparam SequenceIdType    value type for which to perform the increment operation
     * @param diff          number of steps to increment by
     * @param self          iterator to increment
     * @return              copy of incremented iterator
     */
    template<typename SequenceIdType>
    SequenceIdIterator<SequenceIdType> operator+(const typename std::iterator_traits<SequenceIdIterator<SequenceIdType> >::difference_type& diff, SequenceIdIterator<SequenceIdType>& self) {
      return self + diff;
    }




  } /* namespace iterator */
} /* namespace bliss */

#endif /* SEQUENCE_ID_ITERATOR_HPP_ */
