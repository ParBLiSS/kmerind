/**
 * @file    constant_iterator.hpp
 * @ingroup bliss::iterators
 * @author  tpan
 * @brief   contains a constant iterator
 * @details
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef CONSTANT_ITERATOR_HPP_
#define CONSTANT_ITERATOR_HPP_

#include <type_traits>
#include <iterator>


namespace bliss
{
  namespace iterator
  {

    /**
     * @class    bliss::iterator::ConstantIterator
     * @brief    iterator that returns a constant value for all positions.
     * @tparam T the type of data to over.
     */
    template<typename T>
    class ConstantIterator :  public std::iterator<std::random_access_iterator_tag,
                                                   T,
                                                   typename std::conditional<std::is_integral<T>::value,
                                                     ptrdiff_t, T>::type>
    {
      protected:

        /// difference type
        using D = typename std::iterator_traits<ConstantIterator<T> >::difference_type;

        /// current value;
        mutable T val;

      public:

        /**
         * constructor, with stride defaults to 1.
         * @param _val    the constant value that the iterator returns.
         */
        ConstantIterator(const T& _val) : val(_val) {};

        /**
         * default constructor, sets value to 0.
         */
        ConstantIterator() : val(0) {};

        /**
         * default copy constructor
         * @param other  instance of ConstantIterator to copy from
         */
        ConstantIterator(const ConstantIterator<T> & other) : val(other.val) {};

        /**
         * default copy assignment operator
         * @param other  instance of ConstantIterator to copy from
         * @return reference to self
         */
        ConstantIterator<T>& operator=(const ConstantIterator<T> & other) {
          val = other.val;
          return *this;
        }

        /**
         * default move constructor
         * @param other  instance of ConstantIterator to move from
         */
        ConstantIterator(ConstantIterator<T> && other) : val(other.val) {
           other.val = 0;
        };

        /**
         * default move assignment operator
         * @param other  instance of ConstantIterator to move from
         * @return reference to self
         */
        ConstantIterator<T>& operator=(ConstantIterator<T> && other) {
          val = other.val;         other.val = 0;
        };

        /**
         * default destructor
         */
        virtual ~ConstantIterator() {};


        /**
         * @brief   pre-increment operator: ++iter.  For constant iterator, no op.
         * @return  reference to incremented iterator
         */
        ConstantIterator<T>& operator++() {
          return *this;
        }

        /**
         * @brief   post-increment operator: iter++.  for constant iterator, no op.
         * @param   dummy for c++ to identify this as post increment.
         * @return  incremented copy of iterator
         */
        ConstantIterator<T> operator++(int) {
          return ConstantIterator<T>(*this);
        }

        /**
         * @brief compare to other iterator for equality.  only inspect the val
         * @param other   iterator to compare to.
         * @return  bool, true if equal, false otherwise.
         */
        bool operator==(const ConstantIterator<T>& other) {
          return val == other.val;
        }

        /**
         * @brief compare to other iterator for inequality
         * @param other   iterator to compare to.
         * @return  bool, true if not equal, false otherwise.
         */
        bool operator!=(const ConstantIterator<T>& other) {
          return !(this->operator==(other));
        }

        /**
         * @brief dereference function, *iter
         * @return  current value
         */
        const T operator*() const {
          return val;
        }

        /**
         * @brief   pre-decrement operator: --iter.  for constant iterator, no op.
         * @return  reference to decremented iterator
         */
        ConstantIterator<T>& operator--() {
          return *this;
        }

        /**
         * @brief   post-decrement operator: iter--.  for constant iterator, no op.
         * @param   dummy for c++ to identify this as post decrement.
         * @return  decremented copy of iterator
         */
        ConstantIterator<T> operator--(int) {
          return ConstantIterator<T>(*this);
        }

        /**
         * @brief   arithmetic increment by multiple steps.  for constant iterator, no op.
         * @param diff    number of steps to increment by
         * @return  incremented copy of iterator
         */
        ConstantIterator<T> operator+(const D& diff) {
          return ConstantIterator<T>(*this);
        }


        /**
         * @brief   arithmetic decrement by multiple steps.  for constant iterator, no op.
         * @param diff    number of steps to decrement by
         * @return  decremented copy of iterator
         */
        ConstantIterator<T> operator-(const D& diff) {
          return ConstantIterator<T>(*this);
        }

        /**
         * @brief difference between 2 iterators.  for constant iterator, depends on val.
         * @param other         the iterator to subtract by
         * @return              distance between the iterators.  for constant iterator, depends on val
         */
        D operator-(const ConstantIterator<T>& other) {
          return (val - other.val);
        }


        /**
         * @brief compare to other iterator: >.  for constant iterator, depends on val
         * @param other   iterator to compare to.
         * @return  bool, true if val is greater than, false otherwise.
         */
        bool operator>(const ConstantIterator<T>& other) {
          return val > other.val;
        }

        /**
         * @brief compare to other iterator: <.  for constant iterator, depends on val
         * @param other   iterator to compare to.
         * @return  bool, true if val is less than, false otherwise.
         */
        bool operator<(const ConstantIterator<T>& other) {
          return val < other.val;
        }

        /**
         * @brief compare to other iterator: >=.  for constant iterator, depends on val
         * @param other   iterator to compare to.
         * @return  bool, true if val is greater than or equal to, false otherwise.
         */
        bool operator>=(const ConstantIterator<T>& other) {
          return val >= other.val;
        }

        /**
         * @brief compare to other iterator: <=.  for constant iterator, depends on val
         * @param other   iterator to compare to.
         * @return  bool, true if val is less than or equal to, false otherwise.
         */
        bool operator<=(const ConstantIterator<T>& other) {
          return val <= other.val;
        }

        /**
         * @brief   arithmetic increment by multiple steps.  for constant iterator, no op
         * @param diff    number of steps to increment by
         * @return  reference to incremented iterator
         */
        ConstantIterator<T>& operator+=(const D& diff) {
          return *this;
        }

        /**
         * @brief   arithmetic decrement by multiple steps.  for constant iterator, no op
         * @param diff    number of steps to decrement by
         * @return  reference to decremented iterator
         */
        ConstantIterator<T>& operator-=(const D& diff) {
          return *this;
        }

        /**
         * @brief   offset dereference operator. for constant iterator, always the same.
         * @param i offset at which the value is retrieved.
         * @return  value for ith offset, always the same as val.
         */
        T operator[](const D& i) {
          return val;
        }

    };

    /**
     * @brief increment operator, with first operand being a number and second being an iterator  n + iter.  for constant iterator, no op
     * @tparam T    value type for which to perform the increment operation
     * @param diff          number of steps to increment by
     * @param self          iterator to increment
     * @return              copy of incremented iterator
     */
    template<typename T>
    ConstantIterator<T> operator+(const typename std::iterator_traits<ConstantIterator<T> >::difference_type& diff, ConstantIterator<T>& self) {
      return ConstantIterator<T>(self);
    }




  } /* namespace iterator */
} /* namespace bliss */

#endif /* CONSTANT_ITERATOR_HPP_ */
