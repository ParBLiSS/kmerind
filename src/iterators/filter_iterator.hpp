/**
 * filter_iterator.hpp
 *
 *  Created on: Feb 5, 2014
 *      Author: Tony Pan <tpan7@gatech.edu>
 *
 *  patterned after boost's filter iterator.  uses functor class.  allows variadic parameters
 *  support functor (possibly overloaded operator() ), function pointer.
 *
 *  Does NOT support member function pointer.  for those, recommend wrap in a functor.
 *  NOTE: function pointer performs far worse then functor.  probably due to compiler optimization
 *
 *  Relies on compiler to cast the input parameter to the specific types needed by the functor, including constness and references.
 *
 *     - no std::function - general feeling is that it's slower.
 *
 *  TODO:
 *  commenting.
 */

#ifndef FILTER_ITERATOR_HPP_
#define FILTER_ITERATOR_HPP_

#include <iterator>
#include "iterators/function_traits.hpp"

namespace bliss
{

  namespace iterator
  {

    /**
     * filter iterator class
     *
     *
     * filters the element in the list to return only specific ones that match the criteria.
     * it supports all iterator categories for base, but only implement bidirectional iterator operations.
     *
     * can't be a random access iterator itself (because it is unknown how many elements there are, so +/-/+=/-=/[] n operators all are non-sensical)
     *
     *
     * inheriting from std::iterator ONLY to get iterator_traits support.
     *
     * careful with the use of enable_if.  false causes the "type" typedef inside not to be defined, creating an error.
     *   this may be a SFINAE usage or implementation issue.
     */
    template<typename Filter, typename Iterator>
    class filter_iterator : public std::iterator<
        typename std::conditional<
            std::is_same<
                typename std::iterator_traits<Iterator>::iterator_category,
                std::random_access_iterator_tag>::value,
            std::bidirectional_iterator_tag,
            typename std::iterator_traits<Iterator>::iterator_category>::type,
        typename std::iterator_traits<Iterator>::value_type, std::ptrdiff_t,
        typename std::iterator_traits<Iterator>::pointer,
        typename std::iterator_traits<Iterator>::reference>
    {
      protected:
        // define first, to avoid -Wreorder error (where the variables are initialized before filter_iterator::Filter, etc are defined.
        typedef std::iterator_traits<Iterator> base_traits;
        typedef bliss::functional::function_traits<Filter,
            typename std::iterator_traits<Iterator>::value_type> functor_traits;

        Iterator _curr;
        Iterator _start;
        Iterator _end;
        Filter _f;
        bool before_start;

      public:
        typedef filter_iterator<Filter, Iterator> type;
        typedef std::iterator_traits<type> traits;

        typedef typename std::conditional<
            std::is_same<
                typename std::iterator_traits<Iterator>::iterator_category,
                std::random_access_iterator_tag>::value,
            std::bidirectional_iterator_tag,
            typename std::iterator_traits<Iterator>::iterator_category>::type iterator_category;
        typedef typename base_traits::value_type value_type;
        typedef std::ptrdiff_t difference_type;
        typedef typename base_traits::reference reference;
        typedef typename base_traits::pointer pointer;

        typedef value_type base_value_type;

        // accessors
        Iterator& getBase()
        {
          return _curr;
        }

        // class specific constructor
        filter_iterator(const Filter & f, const Iterator& curr,
                        const Iterator& end)
            : _curr(curr), _start(curr), _end(end), _f(f), before_start(false)
        {
        }
        ;

        filter_iterator(const Filter & f, const Iterator& curr) // for end iterator.
            : _curr(curr), _start(curr), _end(curr), _f(f), before_start(false)
        {
        }
        ;

        // for all classes
        // no EXPLICIT so can copy assign.
        filter_iterator(const type& Other)
            : _curr(Other._curr), _start(Other._start), _end(Other._end),
              _f(Other._f), before_start(Other.before_start)
        {
        }
        ;

        type& operator=(const type& Other)
        {
          _curr = Other._curr;
          _start = Other._start;
          _end = Other._end;
          _f = Other._f;
          before_start = Other.before_start;
          return *this;
        }

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

          return ++output;
        }

        // input iterator specific
        inline bool operator==(const type& rhs)
        {
          return _curr == rhs._curr;
        }

        inline bool operator!=(const type& rhs)
        {
          return _curr != rhs._curr;
        }

        inline value_type operator*()
        {
          return *_curr;
        }
        inline Iterator& operator()()
        {
          return _curr;
        }

        // no -> operator.  -> returns a pointer to the value held in iterator.
        // however, in filter iterator, we are not holding data, so pointer does not make sense.

        // NO support for output iterator

        // forward iterator
        template<
            typename = typename std::enable_if<
                std::is_same<iterator_category, std::forward_iterator_tag>::value ||
                std::is_same<iterator_category, std::bidirectional_iterator_tag>::value ||
                std::is_same<iterator_category, std::random_access_iterator_tag>::value >::type >
        explicit filter_iterator()
            : _curr(), _start(), _end(), _f(), before_start(false)
        {
        }
        ;

        // bidirectional iterator
        /**
         * semantics of -- does not have a bound on the start side.
         */
        typename std::enable_if<
            std::is_same<iterator_category, std::bidirectional_iterator_tag>::value || std::is_same<
                iterator_category, std::random_access_iterator_tag>::value, type>::type&
        operator--()
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
        typename std::enable_if<
            std::is_same<iterator_category, std::bidirectional_iterator_tag>::value || std::is_same<
                iterator_category, std::random_access_iterator_tag>::value, type>::type operator--(
            int)
        {
          type output(*this);

          return --output;
        }

        // random access iterator requirements.
        // DO NOT ALLOW RANDOM ACCESS ITERATOR operators for now.
//       typename std::enable_if<std::is_same<iterator_category, std::random_access_iterator_tag>::value,
//                                        type>::type
//       operator+(difference_type n)
//       {
//         type output(*this);
//         output += n;
//         return output;
//       }
//
//       friend
//       typename std::enable_if<std::is_same<iterator_category, std::random_access_iterator_tag>::value,
//                                        type>::type
//       operator+(difference_type n, const type& right)
//       {  return right + n; }
//
//       typename std::enable_if<std::is_same<iterator_category, std::random_access_iterator_tag>::value,
//                                        type>::type
//       operator-(difference_type n)
//       {
//         type output(*this);
//         output += -n;
//         return output;
//       }
//
//       typename std::enable_if<std::is_same<iterator_category, std::random_access_iterator_tag>::value,
//                                        difference_type>::type
//       operator-(const type& other)
//       {
//         // count the number of entries between these two.
//         type lhs(*this);
//         type rhs(other);
//
//         difference_type dist = 0;
//         if (lhs > rhs)
//         {
//           while (lhs > rhs)
//           {
//             ++rhs;
//             ++dist;
//           }
//         }
//         else if (lhs < rhs)
//         {
//           while (lhs < rhs)
//           {
//             ++lhs;
//             --dist;
//           }
//         }
//         // else distance is 0.
//         return dist;
//       }
//
//       inline
//       typename std::enable_if<std::is_same<iterator_category, std::random_access_iterator_tag>::value,
//                                        bool>::type
//       operator<(const type& rhs)
//       { return _curr < rhs._curr; }
//
//       inline
//       typename std::enable_if<std::is_same<iterator_category, std::random_access_iterator_tag>::value,
//                                        bool>::type
//       operator>(const type& rhs)
//       { return _curr > rhs._curr; }
//
//       inline
//       typename std::enable_if<std::is_same<iterator_category, std::random_access_iterator_tag>::value,
//                                        bool>::type
//       operator<=(const type& rhs)
//       { return _curr <= rhs._curr; }
//
//       inline
//       typename std::enable_if<std::is_same<iterator_category, std::random_access_iterator_tag>::value,
//                                        bool>::type
//       operator>=(const type& rhs)
//       { return _curr >= rhs._curr; }
//
//
//       typename std::enable_if<std::is_same<iterator_category, std::random_access_iterator_tag>::value,
//                                        type>::type&
//       operator+=(difference_type n) {
//         if (n > 0)
//           for (difference_type i = 0; i < n; ++i)
//             ++(*this);
//         else if (n < 0)
//           for (difference_type i = 0; i < n; ++i)
//             if (_curr >= _start)
//               --(*this);
//         return *this;
//       }
//
//       typename std::enable_if<std::is_same<iterator_category, std::random_access_iterator_tag>::value,
//                                        type>::type&
//       operator-=(difference_type n)
//       {
//         *this += -n;
//         return *this;
//       }
//
//       typename std::enable_if<std::is_same<iterator_category, std::random_access_iterator_tag>::value,
//                                        value_type>::type
//       operator[](difference_type n) {
//         type output(*this);
//         output += n;
//         return *output;
//       }

    };

  } // iterator
} // bliss
#endif /* FILTER_ITERATOR_HPP_ */
