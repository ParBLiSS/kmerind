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
 * @file    sequence_id_iterator.hpp
 * @ingroup bliss::iterator
 * @author  Tony Pan <tpan7@gatech.edu>
 * @brief   contains a sequence_id iterator
 * @details
 *
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
     * @brief    Iterator that takes an input sequenceId and increment the positions.
     * @details  sequence_id iterator counts from the start to the end in regular steps.
     *           this is useful to establish index values for some other array.
     *
     * @note     this is similar to a counting iterator, except that the data type is more complex.
     *            a sequence Id may have multiple parts to represent file/read/position inside read, etc.
     *
     *           for short sequences, id is the position in file, and the position field should be incremented from 0.
     *           for long sequences, the position field is the poisiton in file, so it should be incremented from that.
     *
     *           Note that we use standard arithmetic and logical operators on the SequenceIdType, so any types that have those operators would works as SequenceIdType.
     *           Also, this means that LongSequenceKmerId and ShortSequenceKmerId need to have these operators defined.
     *
     *           The constructor requires that the parameter be of type SequenceIdType.  If the input type is convertible to SequenceIdType, that should be okay as well.
     *
     * @tparam SequenceIdType the type of data to over.
     */
    template<typename SequenceIdType>
    class SequenceIdIterator :  public std::iterator<std::random_access_iterator_tag, SequenceIdType>
    {
      protected:

        /// difference type
        using D = std::ptrdiff_t;

        /// position type
        using T = size_t;

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
        SequenceIdIterator(const SequenceIdType& _start, const T &_stride = 1) : stride(_stride), val(_start) {};


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
            //printf("INCREMENT id: pos %lu, id %lu, file %u\n", val.get_pos(), val.get_id(), val.get_file_id());
          val += stride;
          return *this;
        }

        /**
         * @brief   post-increment operator: iter++
         * @param   dummy for c++ to identify this as post increment.
         * @return  incremented copy of iterator
         */
        SequenceIdIterator<SequenceIdType> operator++(int) {
          SequenceIdIterator<SequenceIdType> out(*this);
          this->operator++();
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
          val -= stride;
          return *this;
        }

        /**
         * @brief   post-decrement operator: iter--
         * @param   dummy for c++ to identify this as post decrement.
         * @return  decremented copy of iterator
         */
        SequenceIdIterator<SequenceIdType> operator--(int) {
          SequenceIdIterator<SequenceIdType> out(*this);
          this->operator--();
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
         * @brief difference between 2 iterators;  relies on regularity of the strides
         * @param other         the iterator to subtract by
         * @return              distance between the iterators
         */
        D operator-(const SequenceIdIterator<SequenceIdType>& other) {
          return (val - (*other)) / stride;
        }


        /**
         * @brief compare to other iterator: >.  relies on ordering of the sequence Ids correspond to ordering of iterators.
         * @param other   iterator to compare to.
         * @return  bool, true if greater than, false otherwise.
         */
        bool operator>(const SequenceIdIterator<SequenceIdType>& other) {
          return val > other.val;
        }

        /**
         * @brief compare to other iterator: <. relies on ordering of the sequence Ids correspond to ordering of iterators.
         * @param other   iterator to compare to.
         * @return  bool, true if less than, false otherwise.
         */
        bool operator<(const SequenceIdIterator<SequenceIdType>& other) {
          return val < other.val;
        }

        /**
         * @brief compare to other iterator: >=. relies on ordering of the sequence Ids correspond to ordering of iterators.
         * @param other   iterator to compare to.
         * @return  bool, true if greater than or equal to, false otherwise.
         */
        bool operator>=(const SequenceIdIterator<SequenceIdType>& other) {
          return ! (val < other.val);
        }

        /**
         * @brief compare to other iterator: <=. relies on ordering of the sequence Ids correspond to ordering of iterators.
         * @param other   iterator to compare to.
         * @return  bool, true if less than or equal to, false otherwise.
         */
        bool operator<=(const SequenceIdIterator<SequenceIdType>& other) {
          return ! (val > other.val);
        }

        /**
         * @brief   arithmetic increment by multiple steps
         * @param diff    number of steps to increment by
         * @return  reference to incremented iterator
         */
        SequenceIdIterator<SequenceIdType>& operator+=(const D& diff) {
          val += diff * stride;
          return *this;
        }

        /**
         * @brief   arithmetic decrement by multiple steps
         * @param diff    number of steps to decrement by
         * @return  reference to decremented iterator
         */
        SequenceIdIterator<SequenceIdType>& operator-=(const D& diff) {
          val -= diff * stride;
          return *this;
        }

        /**
         * @brief   offset dereference operator
         * @param i offset at which the value is retrieved.
         * @return  value for ith offset
         */
        SequenceIdType operator[](const D& i) {
          SequenceIdType v = val;
          v += stride * i;
          return v;
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
