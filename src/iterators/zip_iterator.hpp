/**
 * @file    zip_iterator.hpp
 * @ingroup bliss::iterators
 * @author  tpan
 * @brief   contains a zip iterator
 * @details
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef ZIP_ITERATOR_HPP_
#define ZIP_ITERATOR_HPP_

#include <type_traits>
#include <iterator>


namespace bliss
{
  namespace iterator
  {

    /**
     * @class    bliss::iterator::ZipIterator
     * @brief    iterator that returns std::pair of the elements from 2 underlying iterators
     * @details  zip iterator simultaneously traverses 2 iterators during increment,
     *           and returns a std::pair of the current elements from the 2 iterators.
     *
     * @tparam FirstIter    Type of first iterator to zip
     * @tparam SecondIter   Type of second iterator to zip
     *
     */
    template<typename FirstIter, typename SecondIter>
    class ZipIterator :  public std::iterator<std::input_iterator_tag,
                                              std::pair<std::iterator_traits<FirstIter>::value_type,
                                                        std::iterator_traits<SecondIter>::value_type> >
    {
      protected:

        /// value type
        using T = typename std::iterator_traits<ZipIterator<FirstIter, SecondIter> >::value_type;

        /// difference type
        using D = typename std::iterator_traits<ZipIterator<FirstIter, SecondIter> >::difference_type;

        /// first iterator to zip
        FirstIter iter1;

        /// second iterator to zip
        SecondIter iter2;

        /// current value;
        T val;

        /**
         * default constructor, sets start to 0, stride to 1, deleted
         */
        ZipIterator() = delete;

        /**
         * default move constructor, deleted
         * @param other  instance of ZipIterator to move from
         */
        ZipIterator(ZipIterator<FirstIter, SecondIter> && other) = delete;

        /**
         * default move assignment operator, deleted.
         * @param other  instance of ZipIterator to move from
         * @return reference to self
         */
        ZipIterator<FirstIter, SecondIter>& operator=(ZipIterator<FirstIter, SecondIter> && other) = delete;



      public:

        /**
         * constructor.  constructs from 2 iterators to zip
         * @param _first      first iterator, whose elements are used as first element of the pairs
         * @param _second     second iterator, whose elements are used as second element of the pairs
         */
        ZipIterator(const FirstIter& _first, const SecondIter &_second) : iter1(_first), iter2(_second), val() {};


        /**
         * default copy constructor
         * @param other  instance of ZipIterator to copy from
         */
        ZipIterator(const ZipIterator<FirstIter, SecondIter> & other) : iter1(other.iter1), iter2(other.iter2), val(other.val) {};

        /**
         * default copy assignment operator
         * @param other  instance of ZipIterator to copy from
         * @return reference to self
         */
        ZipIterator<FirstIter, SecondIter>& operator=(const ZipIterator<FirstIter, SecondIter> & other) {
          iter1 = other.iter1;
          iter2 = other.iter2;
          val = other.val;
          return *this;
        }

        /**
         * default destructor
         */
        virtual ~ZipIterator() {};


        /**
         * @brief   pre-increment operator: ++iter.  increments underlying iterators
         * @return  reference to incremented iterator
         */
        ZipIterator<FirstIter, SecondIter>& operator++() {
          ++iter1;
          ++iter2;
          return *this;
        }

        /**
         * @brief   post-increment operator: iter++.  calls ++iter
         * @param   dummy for c++ to identify this as post increment.
         * @return  incremented copy of iterator
         */
        ZipIterator<FirstIter, SecondIter> operator++(int) {
          ZipIterator<FirstIter, SecondIter> out(*this);
          ++out;
          return out;
        }

        /**
         * @brief compare to other iterator for equality.  inspect the underlying iterators
         * @param other   iterator to compare to.
         * @return  bool, true if equal, false otherwise.
         */
        bool operator==(const ZipIterator<FirstIter, SecondIter>& other) {
          return iter1 == other.iter1 && iter2 == other.iter2;
        }

        /**
         * @brief compare to other iterator for inequality
         * @param other   iterator to compare to.
         * @return  bool, true if not equal, false otherwise.
         */
        bool operator!=(const ZipIterator<FirstIter, SecondIter>& other) {
          return !(this->operator==(other));
        }

        /**
         * @brief dereference function, *iter
         * @return  current value
         */
        const T operator*() const {
          val.first = *iter1;
          val.second = *iter2;
          return val;
        }

        /**
         * @brief dereference function, iter->
         * @return  pointer to current value, so can access internal members
         */
        const T* operator->() const {
          this->operator*();
          return &val;
        }

    };


  } /* namespace iterator */
} /* namespace bliss */

#endif /* ZIP_ITERATOR_HPP_ */
