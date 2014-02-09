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
    : public Functor,
      public std::iterator<typename std::iterator_traits<Base_Iterator>::iterator_category,
                           OT,
                           typename std::iterator_traits<Base_Iterator>::difference_type,
                           OT*,
                           OT&
                           >
    {
      protected:
        Base_Iterator _base_iter;
        Functor _f;

      public:
        typedef transform_iterator_base<Functor, Base_Iterator, OT> iterator_type;
        typedef Base_Iterator base_iterator_type;
        typedef Functor functor_type;

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


       functor_type& getFunctor() {
         return _f;
       }
       const functor_type& getFunctor() const {
         return _f;
       }
       base_iterator_type& getBaseIterator() {
         return _base_iter;
       }

       const base_iterator_type& getBaseIterator() const {
         return _base_iter;
       }

       // forward iterator:

       // can't seem to get this way to work.
//       OT operator*(transform_iterator_base const & iter) {
//         return _f.operator()(*(iter.getBaseIterator()));
//       }
       // this way works, but not with "const" on the function.
       OT operator*() {
         return Functor::operator()(*_base_iter);
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
       OT operator[](difference_type n) {
         return _f.operator()(_base_iter[n]);
       }



       transform_iterator_base& operator+=(difference_type n) {
         std::advance(_base_iter, n);
         return *this;
       }

       transform_iterator_base operator+(difference_type n)
       {
         transform_iterator_base other(*this);
         std::advance(other._base_iter, n);
         return other;
       }

       friend transform_iterator_base operator+(difference_type n, const transform_iterator_base& right)
       {  return right + n; }

       transform_iterator_base& operator-=(difference_type n)
       {
          std::advance(_base_iter, -n);
          return *this;
       }

       transform_iterator_base operator-(difference_type n)
       {
         transform_iterator_base other(*this);
         std::advance(other._base_iter, -n);
         return other;
       }

       // Forward iterator requirements


       explicit
       transform_iterator_base(const Base_Iterator& base_iter, const Functor & f)
         : _base_iter(base_iter), _f(f) {};

       explicit
       transform_iterator_base() : _f(), _base_iter() {};

       explicit
       transform_iterator_base(const transform_iterator_base& Other) : _base_iter(Other.getBaseIterator()), _f(Other.getFunctor()) {};

       explicit
       transform_iterator_base(transform_iterator_base& Other) : _base_iter(Other.getBaseIterator()), _f(Other.getFunctor()) {};


       transform_iterator_base& operator=(const transform_iterator_base& Other) {
         this->_f = Other.getFunctor();
         this->_base_iter = Other.getBaseIterator();
         return *this;
       }





       inline friend bool operator==(const transform_iterator_base& lhs,
            const transform_iterator_base& rhs)
       { return lhs._base_iter == rhs._base_iter; }

       inline friend bool operator!=(const transform_iterator_base& lhs,
            const transform_iterator_base& rhs)
       { return lhs._base_iter != rhs._base_iter;  }

       // Random access iterator requirements
       inline friend bool operator<(const transform_iterator_base& lhs,
           const transform_iterator_base& rhs)
       { return lhs._base_iter < rhs._base_iter; }

       inline friend bool operator>(const transform_iterator_base& lhs,
           const transform_iterator_base& rhs)
       { return lhs._base_iter > rhs._base_iter; }

       inline friend bool operator<=(const transform_iterator_base& lhs,
            const transform_iterator_base& rhs)
       { return lhs._base_iter <= rhs._base_iter; }

       inline friend bool operator>=(const transform_iterator_base& lhs,
            const transform_iterator_base& rhs)
       { return lhs._base_iter >= rhs._base_iter; }


    };



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
