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
 * @file    filter_iterator.hpp
 * @ingroup iterators
 * @author  Tony Pan <tpan7@gatech.edu>
 * @brief   an iterator that returns only elements that match a predicate
 * @details
 */

#ifndef FILTER_ITERATOR_HPP_
#define FILTER_ITERATOR_HPP_

#include <iterator>
#include "utils/function_traits.hpp"

namespace bliss
{

  namespace iterator
  {

    // careful with the use of enable_if.  evaluation should occur at function call time,
    //   i.e. class template params will be evaluated with no substitution.
    // instead, a function should declare a template parameter proxy for the class template parameter.
    //   then enable_if evaluates using the proxy.
    //  e.g. template< class c = C; typename std::enable_if<std::is_same<c, std::string>::value, int>::type x = 0>


    /**
     * @class   filter_iterator
     * @brief   filters the element in the list to return only specific ones that match the criteria.
     * @details each increment skips over elements that did not pass the predicate test
     *          dereference returns the base iterator's element that passes the predicate test.
     *  patterned after boost's filter iterator.  uses functor class for predicate.  allows variadic parameters
     *    support functor (possibly overloaded operator() ), function pointer.
     *
     *  Does NOT support member function pointer.  for those, recommend wrap in a functor.
     *  NOTE: function pointer and std::function usually performs far worse then functor.  probably due to compiler optimization
     *
     * it supports all iterator categories for base, but only implement bidirectional iterator operations.
     *
     * can't be a random access iterator itself (because it is unknown how many elements there are, so +/-/+=/-=/[] n operators all are non-sensical)
     *
     */
    template<typename Filter, typename Iterator>
    class filter_iterator : public std::iterator<
        typename std::conditional<
            std::is_same<
                typename std::iterator_traits<Iterator>::iterator_category,
                std::random_access_iterator_tag>::value,
            std::bidirectional_iterator_tag,
            typename std::iterator_traits<Iterator>::iterator_category>::type,
        typename std::iterator_traits<Iterator>::value_type>
    {
      protected:
        // base iterator traits
        typedef std::iterator_traits<Iterator> base_traits;

        // predicate fucntor traits/return type
        typedef bliss::functional::function_traits<Filter,
            typename std::iterator_traits<Iterator>::value_type> functor_traits;

        // curr position
        Iterator _curr;
        // start of base range
        Iterator _start;
        // end of base range
        Iterator _end;

        // predicate function
        Filter _f;
        // flag indicating that nothing has been called on this class yet.
        bool before_start;

        /// DEFINE filter tierator type
        typedef filter_iterator<Filter, Iterator> type;

        /// DEFINE filter_iterator's trait
        typedef std::iterator_traits<type> traits;

      public:

        using iter_cat = typename std::conditional<
            std::is_same<
                typename std::iterator_traits<Iterator>::iterator_category,
                std::random_access_iterator_tag>::value,
            std::bidirectional_iterator_tag,
            typename std::iterator_traits<Iterator>::iterator_category>::type;
//
//        /// DEFINE value type of iterator's elements.
//        using value_type = typename std::remove_reference<typename base_traits::reference>::type;
//
//        /// DEFINE base iterator's value type, should be same as filter_iterator's
//        typedef value_type base_value_type;

        // accessors
//        Iterator& getBaseIterator()
//        {
//          return _curr;
//        }

        Iterator const & getBaseIterator() const
        {
          return _curr;
        }

        /// constructor, for creating start iterator.
        filter_iterator(const Filter & f, Iterator curr,
                        Iterator end)
            : _curr(curr), _start(curr), _end(end), _f(f), before_start(false)
        {
          // find the first position that satisfies the predicate.
          while (_curr != _end && !_f(*_curr)) // need to check to make sure we are not pass the end of the base iterator.
            {
              ++_curr;
            }
        };
        filter_iterator(Iterator curr,
                        Iterator end)
            : _curr(curr), _start(curr), _end(end), before_start(false)
        {
          // find the first position that satisfies the predicate.
          while (_curr != _end && !_f(*_curr)) // need to check to make sure we are not pass the end of the base iterator.
            {
              ++_curr;
            }
        };


        /// constructor, for creating end iterator
        filter_iterator(const Filter & f, Iterator curr) // for end iterator.
            : _curr(curr), _start(curr), _end(curr), _f(f), before_start(false)
        {};
        explicit filter_iterator(Iterator curr) // for end iterator.
            : _curr(curr), _start(curr), _end(curr), before_start(false)
        {};

        /// copy constructor
        filter_iterator(const type& Other)
            : _curr(Other._curr), _start(Other._start), _end(Other._end),
              _f(Other._f), before_start(Other.before_start)
        {};


        /// copy assignment iterator
        type& operator=(const type& Other)
        {
          _curr = Other._curr;
          _start = Other._start;
          _end = Other._end;
          _f = Other._f;
          before_start = Other.before_start;
          return *this;
        }

        /// move constructor
        filter_iterator(type&& Other)
            : _curr(std::move(Other._curr)), _start(std::move(Other._start)), _end(std::move(Other._end)),
              _f(std::move(Other._f)), before_start(std::move(Other.before_start))
        {};


        /// move assignment iterator
        type& operator=(type&& Other)
        {
          _curr = std::move(Other._curr);
          _start = std::move(Other._start);
          _end = std::move(Other._end);
          _f = std::move(Other._f);
          before_start = std::move(Other.before_start);
          return *this;
        }


        /// increment to next matching element in base iterator
        type& operator++()
        {  // if _curr at end, subsequent calls should not move _curr.
           // on call, if not at end, need to move first then evaluate.
          if (_curr == _end)  // if at end, don't move it.
            return *this;

          ++_curr;            // else move forward 1, and check
          while (_curr != _end && !_f(*_curr)) // need to check to make sure we are not pass the end of the base iterator.
          {
            ++_curr;
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

        /// compare 2 filter iterators
        inline bool operator==(const type& rhs) const
        {
          return _curr == rhs._curr;
        }

        /// compare 2 filter iterators
        inline bool operator!=(const type& rhs) const
        {
          return _curr != rhs._curr;
        }

        /// dereference operator.  returned entry passes the predicate test.  guaranteed to be at a valid position
//        inline typename base_traits::value_type& operator*() {
//          return *_curr;
//        }
//
//        inline typename base_traits::value_type* operator->() {
//          return _curr.operator->();
//        }

        /// dereference operator.  returned entry passes the predicate test.  guaranteed to be at a valid position
        inline typename base_traits::value_type const & operator*() const {
          return *_curr;
        }

        inline const typename base_traits::value_type * operator->() const {
          return _curr.operator->();
        }


        /// referece operator.  returned entry passes the predicate test.  guaranteed to be at a valid position
        inline Iterator& operator()() const
        {
          return _curr;
        }

//        /// returns a pointer to tthe value held in iterator.
//        inline typename base_traits::pointer operator->()
//        {
//          return &(*_curr);
//        }

        // NO support for output iterator at this point, since modifying a value can prevent multipass requirement of forward iterator.

        // forward iterator required default constructor
        template<typename IterCat = iter_cat,
            typename = typename std::enable_if<
                std::is_same<IterCat, std::forward_iterator_tag>::value ||
                std::is_same<IterCat, std::bidirectional_iterator_tag>::value ||
                std::is_same<IterCat, std::random_access_iterator_tag>::value, int >::type >
        explicit filter_iterator()
            : _curr(), _start(), _end(), _f(), before_start(false)
        {
        }
        ;

        // bidirectional iterator
        /**
         * semantics of -- does not have a bound on the start side.
         */
        template <typename IterCat = iter_cat,
            typename = typename std::enable_if<
              std::is_same<IterCat, std::bidirectional_iterator_tag>::value ||
              std::is_same<IterCat, std::random_access_iterator_tag>::value, int>::type >
        type& operator--()
        {
          if (_curr == _start)  // at beginning.  don't move it.
            before_start = true;

          if (before_start)
            return *this;

          --_curr;            // else move back 1, and check
          while (_curr != _start && !_f(*_curr)) // need to check to make sure we are not pass the end of the base iterator.
          {
            --_curr;
          }
          return *this;
        }

        /**
         * make a copy, then move it back 1.
         */
        template <typename IterCat = iter_cat,
            typename = typename std::enable_if<
              std::is_same<IterCat, std::bidirectional_iterator_tag>::value ||
              std::is_same<IterCat, std::random_access_iterator_tag>::value, int>::type >
        type& operator--(int)
        {
          type output(*this);
          this->operator--();
          return output;
        }

        // random access iterator requirements.
        // DO NOT ALLOW RANDOM ACCESS ITERATOR operators for now.

    };

  } // iterator
} // bliss
#endif /* FILTER_ITERATOR_HPP_ */
