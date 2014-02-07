/**
 * transform_iterator.hpp
 *
 *  Created on: Feb 5, 2014
 *      Author: tpan
 */

#ifndef TRANSFORM_ITERATOR_HPP_
#define TRANSFORM_ITERATOR_HPP_

#include <iterator>

namespace bliss
{

namespace iterator
{

  template<typename Functor,
           typename Base_Iterator,
           typename OT >
  class transform_iterator_base
    : public std::iterator<typename std::iterator_traits<Base_Iterator>::iterator_category,
                           OT,
                           typename std::iterator_traits<Base_Iterator>::difference_type,
                           OT*,
                           OT&
                           >
    {
      public:
        typedef transform_iterator_base<Functor, Base_Iterator, OT> iterator_type;

      protected:
        typedef std::iterator_traits<iterator_type>                     traits_type;
        typedef std::iterator_traits<Base_Iterator>                     base_traits_type;

      public:
        typedef typename base_traits_type::iterator_category iterator_category;
        typedef OT    value_type;
        typedef typename base_traits_type::difference_type   difference_type;
        typedef OT&   reference;
        typedef OT*     pointer;

         typedef typename base_traits_type::value_type input_value_type;

         typedef Base_Iterator base_iterator_type;
         typedef Functor functor_type;

      protected:
        functor_type _f;
        base_iterator_type _base_iter;

      public:
       functor_type& getFunctor() const {
         return _f;
       }
       base_iterator_type& getBaseIterator() {
         return _base_iter;
       }

       transform_iterator_base(const Base_Iterator& base_iter, const Functor & f)
         : _base_iter(base_iter), _f(f) {};

       transform_iterator_base(Base_Iterator& base_iter, Functor & f)
         : _base_iter(base_iter), _f(f) {};


       // forward iterator:
       OT operator*(const iterator_type& iter) {
         return iter._f(*_base_iter);
       }

       // no -> operator

       transform_iterator_base& operator++() {
         ++_base_iter;
         return *this;
       }

       transform_iterator_base operator++(int) {
         return transform_iterator_base(_base_iter++, _f);
       }

       // bidirectional iterator
       transform_iterator_base& operator--() {
         --_base_iter;
         return *this;
       }

       transform_iterator_base operator--(int) {
         return transform_iterator_base(_base_iter--, _f);
       }

       // random access iterator requirements.
       OT operator[](const difference_type& n) {
         return _f(_base_iter[n]);
       }

       transform_iterator_base& operator+=(const difference_type& n) {
         _base_iter += n;
         return *this;
       }

       transform_iterator_base operator+(const difference_type& n)
       {
         return transform_iterator_base(_base_iter + n, _f);
       }

       transform_iterator_base& operator-=(const difference_type& n)
       {
          _base_iter -= n;
          return *this;
       }

       transform_iterator_base operator-(const difference_type& n)
       {
         return transform_iterator_base(_base_iter - n, _f);
       }


       // Forward iterator requirements



       transform_iterator_base() {};

       explicit
       transform_iterator_base(const transform_iterator_base& Other) : _base_iter(Other.getBaseIterator()), _f(Other.getFunctor()) {};

       explicit
       transform_iterator_base(transform_iterator_base& Other) : _base_iter(Other.getBaseIterator()), _f(Other.getFunctor()) {};


       iterator_type& operator=(const iterator_type& Other) {
         _f = Other.getFunctor();
         _base_iter = Other.getBaseIterator();
         return *this;
       }

       iterator_type& operator=(iterator_type& Other) {
         _f = Other.getFunctor();
         _base_iter = Other.getBaseIterator();
         return *this;
       }



    };


  template<typename iterator_type>
  inline bool operator==(const iterator_type& lhs,
       const iterator_type& rhs)
  { return lhs.getBaseIterator() == rhs.getBaseIterator(); }

  template<typename iterator_type>
  inline bool operator!=(const iterator_type& lhs,
       const iterator_type& rhs)
  { return lhs.getBaseIterator() != rhs.getBaseIterator(); }

  // Random access iterator requirements
  template<typename iterator_type>
  inline bool operator<(const iterator_type& lhs,
      const iterator_type& rhs)
  { return lhs.getBaseIterator() < rhs.getBaseIterator(); }

  template<typename iterator_type>
  inline bool operator>(const iterator_type& lhs,
      const iterator_type& rhs)
  { return lhs.getBaseIterator() > rhs.getBaseIterator(); }

  template<typename iterator_type>
  inline bool operator<=(const iterator_type& lhs,
       const iterator_type& rhs)
  { return lhs.getBaseIterator() <= rhs.getBaseIterator(); }

  template<typename iterator_type>
  inline bool operator>=(const iterator_type& lhs,
       const iterator_type& rhs)
  { return lhs.getBaseIterator() >= rhs.getBaseIterator(); }

  template<typename Functor,
           typename Base_Iterator >
  class transform_iterator
    : public transform_iterator_base<Functor,
                                     Base_Iterator,
                                     typename std::result_of<Functor(typename Base_Iterator::value_type)>::type
                                     >
  {
  };



} // iterator
} // bliss
#endif /* TRANSFORM_ITERATOR_HPP_ */
