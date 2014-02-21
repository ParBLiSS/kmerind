/**
 * transform_iterator.hpp
 *
 *  Created on: Feb 5, 2014
 *      Author: tpan
 *
 *  patterned after boost's transform iterator.  no multiple inheritance.  uses functor class.  allows variadic parameters
 *  support functor (possibly multiple operator() ), function pointer.
 *
 *  Does NOT support member function pointer.  for those, recommend wrap in a functor.
 *  NOTE: function pointer performs far worse then functor.  probably due to compiler optimization
 *
 *
 *  Relies on compiler to cast the input parameter to the specific types needed by the functor, including constness and references.
 *
 *     - no std::function - general feeling is that it's slower.
 *
 *  TODO:
 *  commenting.
 */

#ifndef TRANSFORM_ITERATOR_HPP_
#define TRANSFORM_ITERATOR_HPP_

#include <iterator>
#include "iterators/function_traits.hpp"

namespace bliss
{

namespace iterator
{

  /**
   * transform iterator class
   *
   * transforms each element in the list. otherwise it's the same kind of iterator as the base.
   *
   * inheriting from std::iterator ONLY to get iterator_traits support.
   *
   * careful with the use of enable_if.  false causes the "type" typedef inside not to be defined, creating an error.
   *   this may be a SFINAE usage or implementation issue.
   */
  template<typename Transformer,
           typename Iterator
           >
  class transform_iterator
    : public std::iterator<typename std::iterator_traits<Iterator>::iterator_category,
                           typename std::remove_reference<typename bliss::functional::function_traits<Transformer, typename std::iterator_traits<Iterator>::value_type>::return_type>::type,
                           typename std::iterator_traits<Iterator>::difference_type,
                           typename std::add_pointer<typename std::remove_reference<typename bliss::functional::function_traits<Transformer, typename std::iterator_traits<Iterator>::value_type>::return_type>::type>::type,
                           typename std::add_rvalue_reference<typename std::remove_reference<typename bliss::functional::function_traits<Transformer, typename std::iterator_traits<Iterator>::value_type>::return_type>::type>::type
                           >
    {
      protected:
        // define first, to avoid -Wreorder error (where the variables are initialized before transform_iterator::Transformer, etc are defined.
        typedef std::iterator_traits<Iterator>                                                         base_traits;
        typedef bliss::functional::function_traits<Transformer, typename std::iterator_traits<Iterator>::value_type>                                    functor_traits;

        Iterator                                                                                       _base;
        Transformer                                                                                             _f;

      public:
        typedef transform_iterator< Transformer, Iterator >                                                type;
        typedef std::iterator_traits<type>                                                                  traits;

        typedef typename base_traits::iterator_category                                                     iterator_category;
        typedef typename std::remove_reference<typename functor_traits::return_type>::type                  value_type;
        typedef typename base_traits::difference_type                                                       difference_type;
        typedef typename std::add_rvalue_reference<value_type>::type                                        reference_type;
        typedef typename std::add_pointer<value_type>::type                                                 pointer_type;

        typedef typename base_traits::value_type                                                            base_value_type;


        // class specific constructor
        explicit
        transform_iterator(const Iterator& base_iter, const Transformer & f)
          : _base(base_iter), _f(f) {};


        // for all classes
        // no EXPLICIT so can copy assign.
        transform_iterator(const type& Other) : _base(Other._base), _f(Other._f) {};

        type& operator=(const type& Other) {
          _f = Other._f;
          _base = Other._base;
          return *this;
        }

        type& operator++() {
          ++_base;
          return *this;
        }

        type operator++(int) {
          return type(_base++, _f);
        }


        // accessor functions for internal state.
        Transformer& getTransformer() {
          return _f;
        }
        const Transformer& getTransformer() const {
          return _f;
        }
        Iterator& getBaseIterator() {
          return _base;
        }
        const Iterator& getBaseIterator() const {
          return _base;
        }

        // input iterator specific
        inline bool operator==(const type& rhs)
        { return _base == rhs._base; }

        inline bool operator!=(const type& rhs)
        { return _base != rhs._base; }

        inline value_type operator*() {
          return _f(*_base);
        }

        // no -> operator.  -> returns a pointer to the value held in iterator.
        // however, in transform iterator, we are not holding data, so pointer does not make sense.

        // NO support for output iterator


        // forward iterator
        template<typename = std::enable_if<std::is_same<iterator_category, std::forward_iterator_tag>::value ||
                                           std::is_same<iterator_category, std::bidirectional_iterator_tag>::value ||
                                           std::is_same<iterator_category, std::random_access_iterator_tag>::value > >
        explicit
        transform_iterator() : _f(), _base() {};





       // bidirectional iterator
        typename std::enable_if<std::is_same<iterator_category, std::bidirectional_iterator_tag>::value ||
                                std::is_same<iterator_category, std::random_access_iterator_tag>::value,
                                type
                               >::type&
       operator--() {
         --_base;
         return *this;
       }

       typename std::enable_if<std::is_same<iterator_category, std::bidirectional_iterator_tag>::value ||
                               std::is_same<iterator_category, std::random_access_iterator_tag>::value,
                               type
                              >::type
       operator--(int) {
         return type(_base--, _f);
       }


       // random access iterator requirements.
       typename std::enable_if<std::is_same<iterator_category, std::random_access_iterator_tag>::value,
                                        type>::type
       operator+(difference_type n)
       {
         type output(*this);
         std::advance(output._base, n);
         return output;
       }

       friend
       typename std::enable_if<std::is_same<iterator_category, std::random_access_iterator_tag>::value,
                                        type>::type
       operator+(difference_type n, const type& right)
       {  return right + n; }

       typename std::enable_if<std::is_same<iterator_category, std::random_access_iterator_tag>::value,
                                        type>::type
       operator-(difference_type n)
       {
         type output(*this);
         std::advance(output._base, -n);
         return output;
       }

       typename std::enable_if<std::is_same<iterator_category, std::random_access_iterator_tag>::value,
                                        difference_type>::type
       operator-(const type& other)
       {
         return _base - other._base;
       }

       inline
       typename std::enable_if<std::is_same<iterator_category, std::random_access_iterator_tag>::value,
                                        bool>::type
       operator<(const type& rhs)
       { return _base < rhs._base; }

       inline
              typename std::enable_if<std::is_same<iterator_category, std::random_access_iterator_tag>::value,
                                               bool>::type
              operator>(const type& rhs)
       { return _base > rhs._base; }

       inline
              typename std::enable_if<std::is_same<iterator_category, std::random_access_iterator_tag>::value,
                                               bool>::type
              operator<=(const type& rhs)
       { return _base <= rhs._base; }

       inline
              typename std::enable_if<std::is_same<iterator_category, std::random_access_iterator_tag>::value,
                                               bool>::type
              operator>=(const type& rhs)
       { return _base >= rhs._base; }


       typename std::enable_if<std::is_same<iterator_category, std::random_access_iterator_tag>::value,
                                        type>::type&
       operator+=(difference_type n) {
         std::advance(_base, n);
         return *this;
       }

       typename std::enable_if<std::is_same<iterator_category, std::random_access_iterator_tag>::value,
                                        type>::type&
       operator-=(difference_type n)
       {
          std::advance(_base, -n);
          return *this;
       }

       typename std::enable_if<std::is_same<iterator_category, std::random_access_iterator_tag>::value,
                                        value_type>::type
       operator[](difference_type n) {
         return _f(_base[n]);
       }



    };


} // iterator
} // bliss
#endif /* TRANSFORM_ITERATOR_HPP_ */
