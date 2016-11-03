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
 * @file    zip_iterator.hpp
 * @ingroup iterators
 * @author  Tony Pan <tpan7@gatech.edu>
 * @brief   contains a zip iterator that allows co-iteration of 2 base iterators.
 * @details
 *
 */
#ifndef ZIP_ITERATOR_HPP_
#define ZIP_ITERATOR_HPP_

#include <cstddef>
#include <type_traits>
#include <iterator>

#include <iostream>  // cout

// TODO: tuple version that supports arbitrary number of base iterators,  via variadic templating.

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
     *           first iterator specified is the "primary", on which equality is tested, e.g. for end condition.
     *
     * @note     This is a forward iterator only.
     *
     * @tparam FirstIter    Type of first iterator to zip
     * @tparam SecondIter   Type of second iterator to zip
     *
     */
    template<typename FirstIter, typename SecondIter>
    class ZipIterator :  public std::iterator<std::input_iterator_tag,
                                              std::pair<typename std::iterator_traits<FirstIter>::value_type,
                                                        typename std::iterator_traits<SecondIter>::value_type > >
    {
      protected:

        /// value type
//        using T = typename std::iterator_traits<ZipIterator<FirstIter, SecondIter> >::value_type;
        using T = std::pair<typename std::iterator_traits<FirstIter>::value_type,
            typename std::iterator_traits<SecondIter>::value_type >;

        /// difference type
        using D = std::ptrdiff_t;

        /// first iterator to zip, "primary" iterator that's used for testing stopping criteria.
        FirstIter iter1;

        /// second iterator to zip
        SecondIter iter2;

        /// current value;
        mutable T val;


      public:
        /**
         * default constructor. does nothing.
         */
        ZipIterator() : val() {};


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
         * default move constructor
         * @param other  instance of ZipIterator to move from
         */
        ZipIterator(ZipIterator<FirstIter, SecondIter> && other) : iter1(std::move(other.iter1)), iter2(std::move(other.iter2)), val(std::move(other.val)) {};

        /**
         * default move assignment operator.
         * @param other  instance of ZipIterator to move from
         * @return reference to self
         */
        ZipIterator<FirstIter, SecondIter>& operator=(ZipIterator<FirstIter, SecondIter> && other) {
          iter1 = std::move(other.iter1);
          iter2 = std::move(other.iter2);
          val = std::move(other.val);
          return *this;
        };



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
          this->operator++();
          return out;
        }

        /**
         * @brief compare to other iterator for equality.  inspect the underlying iterators
         * @param other   iterator to compare to.
         * @return  bool, true if equal, false otherwise.
         */
        bool operator==(const ZipIterator<FirstIter, SecondIter>& other) const {
          return iter1 == other.iter1;
        }

        /**
         * @brief compare to other iterator for inequality
         * @param other   iterator to compare to.
         * @return  bool, true if not equal, false otherwise.
         */
        bool operator!=(const ZipIterator<FirstIter, SecondIter>& other) const {
          return !(this->operator==(other));
        }

        /**
         * @brief dereference function, *iter
         * @return  current value
         */

        inline T operator*() const {
          val.first = *iter1;
          val.second = *iter2;
          return val;
        }

        /**
         * @brief dereference function, iter->
         * @return  pointer to current value, so can access internal members
         */
        T* operator->() {
          this->operator*();
          return &val;
        }

        /// readonly accessor for first iterator
        FirstIter const & get_first_iterator() const {
          return iter1;
        }

        /// readonly accessor for second iterator
        SecondIter const & get_second_iterator() const {
          return iter2;
        }


    };


#if (0)
    /**
     * @class    bliss::iterator::ZipIterator
     * @brief    iterator that returns std::pair of the elements from 2 underlying iterators
     * @details  zip iterator simultaneously traverses 2 iterators during increment,
     *           and returns a std::pair of the current elements from the 2 iterators.
     *
     *           first iterator specified is the "primary", on which equality is tested, e.g. for end condition.
     *
     * @note     This is a forward iterator only.
     *
     * @tparam FirstIter    Type of first iterator to zip
     * @tparam SecondIter   Type of second iterator to zip
     *
     */
    template<typename... Iter, typename = typename std::enable_if<(sizeof...(Iter) > 2)>::type >
    class ZipIterator :  public std::iterator<std::input_iterator_tag,
                                              std::tuple<typename std::iterator_traits<Iter>::value_type... > >
    {
      protected:

        /// value type
//        using T = typename std::iterator_traits<ZipIterator<FirstIter, SecondIter> >::value_type;
        using T = std::tuple<typename std::iterator_traits<Iter>::value_type... >;

        /// difference type
        using D = std::ptrdiff_t;

        // index sequence so we can iterate over the iterators.
        static constexpr ::std::array<size_t, sizeof...(Iter)> index = ::bliss::utils::make_array<size_t, sizeof...(Iter)>();

        /// iterators
        using IT = ::std::tuple<Iter... >;
        IT iterators;

        /// current value;
        mutable T val;
//
//        // base cases for recursive increment
//        template <typename ... Iters, size_t N = sizeof...(Iters), typename = typename std::enable_if<(N==1)>::type>
//        void increment(std::tuple<Iters...> & iters) {
//          ++(std::get<0>(iters));
//        }
//        template <typename ... Iters, size_t N = sizeof...(Iters), typename = typename std::enable_if<(N==2)>::type >
//        void increment(std::tuple<Iters...> & iters) {
//          ++(std::get<0>(iters));
//          ++(std::get<1>(iters));
//        }
//        // recursive increment
//        template <typename ... Iters, size_t N = sizeof...(Iters), typename = typename std::enable_if<(N > 2)>::type >
//        void increment(std::tuple<Iters...> & iters) {
//          increment(iters1);
//          increment(iters2);
//        }




      public:
        /**
         * default constructor. does nothing.  here for copying.
         */
        ZipIterator() : val() {};


        /**
         * constructor.  constructs from 2 iterators to zip
         * @param _first      first iterator, whose elements are used as first element of the pairs
         * @param _second     second iterator, whose elements are used as second element of the pairs
         */
        ZipIterator(Iter const & ... iters) : iterators(iters...), val() {};


        /**
         * default copy constructor
         * @param other  instance of ZipIterator to copy from
         */
        ZipIterator(const ZipIterator<Iters... > & other) : iterators(other.iterators), val(other.val) {};

        /**
         * default copy assignment operator
         * @param other  instance of ZipIterator to copy from
         * @return reference to self
         */
        ZipIterator& operator=(const ZipIterator & other) {
          iterators = other.iterators;
          val = other.val;
          return *this;
        }

        /**
         * default move constructor
         * @param other  instance of ZipIterator to move from
         */
        ZipIterator(ZipIterator && other) : iterators(std::move(other.iterators)), val(std::move(other.val)) {};

        /**
         * default move assignment operator.
         * @param other  instance of ZipIterator to move from
         * @return reference to self
         */
        ZipIterator& operator=(ZipIterator && other) {
          iterators = std::move(other.iterators);
          val = std::move(other.val);
          return *this;
        };



        /**
         * default destructor
         */
        virtual ~ZipIterator() {};


        /**
         * @brief   pre-increment operator: ++iter.  increments underlying iterators
         * @return  reference to incremented iterator
         */
        ZipIterator& operator++() {
          // WIP, TODO.  increment<Iter...>();
          return *this;
        }

        /**
         * @brief   post-increment operator: iter++.  calls ++iter
         * @param   dummy for c++ to identify this as post increment.
         * @return  incremented copy of iterator
         */
        ZipIterator operator++(int) {
          ZipIterator out(*this);
          this->operator++();
          return out;
        }

        /**
         * @brief compare to other iterator for equality.  inspect the underlying iterators
         * @param other   iterator to compare to.
         * @return  bool, true if equal, false otherwise.
         */
        bool operator==(const ZipIterator& other) const {
          return iterators == other.iterators;
        }

        /**
         * @brief compare to other iterator for inequality
         * @param other   iterator to compare to.
         * @return  bool, true if not equal, false otherwise.
         */
        bool operator!=(const ZipIterator& other) const {
          return !(this->operator==(other));
        }

        /**
         * @brief dereference function, *iter
         * @return  current value
         */
//        inline T operator*() {
//          std::cout << "zip deref.  before " << val << std::endl;
//
//          val.first = *iter1;
//          val.second = *iter2;
//          std::cout << "zip deref.  after " << val << std::endl;
//          return val;
//        }

        inline T operator*() const {
          std::cout << "zip const deref.  before " << val << std::endl;

//          val = std::make_tuple()
//          for (int i = 0; i < iterators.size(); ++i)
//
//          val.first = *iter1;
//          val.second = *iter2;

          // WIP.  TODO.

          std::cout << "zip const deref.  after " << val << std::endl;
          return val;
        }

        /**
         * @brief dereference function, iter->
         * @return  pointer to current value, so can access internal members
         */
        T* operator->() {
          this->operator*();
          return &val;
        }

        /// readonly accessor for first iterator
        template <size_t I>
        std::tuple_element<I, IT> const & get_iterator() const {
          return iterators[I];
        }



    };

    template<typename... Iter>
    constexpr ::std::array<size_t, sizeof...(Iter)> ZipIterator<Iter...>::index;
#endif


  } /* namespace iterator */
} /* namespace bliss */

#endif /* ZIP_ITERATOR_HPP_ */
