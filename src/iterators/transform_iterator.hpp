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
           typename Base_Iterator>
  class transform_iterator
    : public Functor,
      public std::iterator<typename std::iterator_traits<Base_Iterator>::iterator_category,
                           typename std::remove_reference<typename std::result_of<Functor(typename Base_Iterator::value_type)>::type>::type,
                           typename std::iterator_traits<Base_Iterator>::difference_type,
                           typename std::add_pointer<typename std::remove_reference<typename std::result_of<Functor(typename Base_Iterator::value_type)>::type>::type>::type,
                           typename std::add_rvalue_reference<typename std::remove_reference<typename std::result_of<Functor(typename Base_Iterator::value_type)>::type>::type>::type
                           >
    {
      protected:
        // define first, to avoid -Wreorder error (where the variables are initialized before transform_iterator::Functor, etc are defined.
        Base_Iterator _base_iter;
        Functor _f;
        typedef std::iterator_traits<Base_Iterator>                     base_traits_type;

      public:
        typedef transform_iterator<Functor, Base_Iterator>     iterator_type;
        typedef Base_Iterator                                           base_iterator_type;
        typedef Functor                                                 functor_type;
        typedef std::iterator_traits<iterator_type>                     traits_type;

        typedef typename base_traits_type::iterator_category iterator_category;
        typedef typename std::remove_reference<typename std::result_of<Functor(typename Base_Iterator::value_type)>::type>::type  value_type;
        typedef typename base_traits_type::difference_type   difference_type;
        typedef typename std::add_rvalue_reference<value_type>::type    reference;
        typedef typename std::add_pointer<value_type>::type    pointer;

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
//       OT operator*(transform_iterator const & iter) {
//         return _f.operator()(*(iter.getBaseIterator()));
//       }
       // this way works, but not with "const" on the function.
       value_type operator*() {
         return Functor::operator()(*_base_iter);
       }
       // no -> operator

       transform_iterator& operator++() {
         ++_base_iter;
         return *this;
       }

       transform_iterator operator++(int) {
         return transform_iterator(_base_iter++, _f);
       }

       // bidirectional iterator
       transform_iterator& operator--() {
         --_base_iter;
         return *this;
       }

       transform_iterator operator--(int) {
         return transform_iterator(_base_iter--, _f);
       }

       // random access iterator requirements.
       value_type operator[](difference_type n) {
         return Functor::operator()(_base_iter[n]);
       }



       transform_iterator& operator+=(difference_type n) {
         std::advance(_base_iter, n);
         return *this;
       }

       transform_iterator operator+(difference_type n)
       {
         transform_iterator other(*this);
         std::advance(other._base_iter, n);
         return other;
       }

       friend transform_iterator operator+(difference_type n, const transform_iterator& right)
       {  return right + n; }

       transform_iterator& operator-=(difference_type n)
       {
          std::advance(_base_iter, -n);
          return *this;
       }

       transform_iterator operator-(difference_type n)
       {
         transform_iterator other(*this);
         std::advance(other._base_iter, -n);
         return other;
       }

       // Forward iterator requirements


       explicit
       transform_iterator(const Base_Iterator& base_iter, const Functor & f)
         : _base_iter(base_iter), _f(f) {};

       explicit
       transform_iterator() : _f(), _base_iter() {};

       explicit
       transform_iterator(const transform_iterator& Other) : _base_iter(Other.getBaseIterator()), _f(Other.getFunctor()) {};

       explicit
       transform_iterator(transform_iterator& Other) : _base_iter(Other.getBaseIterator()), _f(Other.getFunctor()) {};


       transform_iterator& operator=(const transform_iterator& Other) {
         this->_f = Other.getFunctor();
         this->_base_iter = Other.getBaseIterator();
         return *this;
       }





       inline friend bool operator==(const transform_iterator& lhs,
            const transform_iterator& rhs)
       { return lhs._base_iter == rhs._base_iter; }

       inline friend bool operator!=(const transform_iterator& lhs,
            const transform_iterator& rhs)
       { return lhs._base_iter != rhs._base_iter;  }

       // Random access iterator requirements
       inline friend bool operator<(const transform_iterator& lhs,
           const transform_iterator& rhs)
       { return lhs._base_iter < rhs._base_iter; }

       inline friend bool operator>(const transform_iterator& lhs,
           const transform_iterator& rhs)
       { return lhs._base_iter > rhs._base_iter; }

       inline friend bool operator<=(const transform_iterator& lhs,
            const transform_iterator& rhs)
       { return lhs._base_iter <= rhs._base_iter; }

       inline friend bool operator>=(const transform_iterator& lhs,
            const transform_iterator& rhs)
       { return lhs._base_iter >= rhs._base_iter; }


    };




} // iterator
} // bliss
#endif /* TRANSFORM_ITERATOR_HPP_ */
