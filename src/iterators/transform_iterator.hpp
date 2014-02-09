/**
 * transform_iterator.hpp
 *
 *  Created on: Feb 5, 2014
 *      Author: tpan
 *
 *  patterned after boost's transform iterator.  had to use multiple inheritance. (else get unmatched operator*()) _f is not used...
 *  TODO:
 *  specialization for windowed iterator and block partitioned iterator.
 *    no need for cyclic partitioned iterator right now, and no need for skipping iterator (random access handles that)
 *  specialization for function pointers (can't use multiple inheritance)
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
           typename Value = typename std::remove_reference<typename std::result_of<Functor(typename Base_Iterator::value_type)>::type>::type
           >
  class transform_iterator
    : public Functor,
      public std::iterator<typename std::iterator_traits<Base_Iterator>::iterator_category,
                           Value,
                           typename std::iterator_traits<Base_Iterator>::difference_type,
                           typename std::add_pointer<Value>::type,
                           typename std::add_rvalue_reference<Value>::type
                           >
    {
      protected:
        // define first, to avoid -Wreorder error (where the variables are initialized before transform_iterator::Functor, etc are defined.
        Base_Iterator                                                   _base;
        Functor                                                         _f;
        typedef std::iterator_traits<Base_Iterator>                     base_traits;

      public:
        typedef transform_iterator<Functor, Base_Iterator>              type;
        typedef Base_Iterator                                           base_type;
        typedef Functor                                                 functor_type;
        typedef std::iterator_traits<type>                              traits;

        typedef typename base_traits::iterator_category                 iterator_category;
        typedef Value                                                   value_type;
        typedef typename base_traits::difference_type                   difference_type;
        typedef typename std::add_rvalue_reference<value_type>::type    reference_type;
        typedef typename std::add_pointer<value_type>::type             pointer_type;

        typedef typename base_traits::value_type                        base_value_type;


       functor_type& getFunctor() {
         return _f;
       }
       const functor_type& getFunctor() const {
         return _f;
       }
       base_type& getBaseIterator() {
         return _base;
       }

       const base_type& getBaseIterator() const {
         return _base;
       }

       // forward iterator:

       // can't seem to get this way to work.
//       OT operator*(transform_iterator const & iter) {
//         return _f.operator()(*(iter.getBaseIterator()));
//       }
       // this way works, but not with "const" on the function.
       value_type operator*() {
         return Functor::operator()(*_base);
       }
       // no -> operator

       transform_iterator& operator++() {
         ++_base;
         return *this;
       }

       transform_iterator operator++(int) {
         return transform_iterator(_base++, _f);
       }

       // bidirectional iterator
       transform_iterator& operator--() {
         --_base;
         return *this;
       }

       transform_iterator operator--(int) {
         return transform_iterator(_base--, _f);
       }

       // random access iterator requirements.
       value_type operator[](difference_type n) {
         return Functor::operator()(_base[n]);
       }



       transform_iterator& operator+=(difference_type n) {
         std::advance(_base, n);
         return *this;
       }

       transform_iterator operator+(difference_type n)
       {
         transform_iterator other(*this);
         std::advance(other._base, n);
         return other;
       }

       friend transform_iterator operator+(difference_type n, const transform_iterator& right)
       {  return right + n; }

       transform_iterator& operator-=(difference_type n)
       {
          std::advance(_base, -n);
          return *this;
       }

       transform_iterator operator-(difference_type n)
       {
         transform_iterator other(*this);
         std::advance(other._base, -n);
         return other;
       }

       // Forward iterator requirements


       explicit
       transform_iterator(const Base_Iterator& base_iter, const Functor & f)
         : _base(base_iter), _f(f) {};

       explicit
       transform_iterator() : _f(), _base() {};

       explicit
       transform_iterator(const transform_iterator& Other) : _base(Other.getBaseIterator()), _f(Other.getFunctor()) {};

       explicit
       transform_iterator(transform_iterator& Other) : _base(Other.getBaseIterator()), _f(Other.getFunctor()) {};


       transform_iterator& operator=(const transform_iterator& Other) {
         this->_f = Other.getFunctor();
         this->_base = Other.getBaseIterator();
         return *this;
       }





       inline friend bool operator==(const transform_iterator& lhs,
            const transform_iterator& rhs)
       { return lhs._base == rhs._base; }

       inline friend bool operator!=(const transform_iterator& lhs,
            const transform_iterator& rhs)
       { return lhs._base != rhs._base;  }

       // Random access iterator requirements
       inline friend bool operator<(const transform_iterator& lhs,
           const transform_iterator& rhs)
       { return lhs._base < rhs._base; }

       inline friend bool operator>(const transform_iterator& lhs,
           const transform_iterator& rhs)
       { return lhs._base > rhs._base; }

       inline friend bool operator<=(const transform_iterator& lhs,
            const transform_iterator& rhs)
       { return lhs._base <= rhs._base; }

       inline friend bool operator>=(const transform_iterator& lhs,
            const transform_iterator& rhs)
       { return lhs._base >= rhs._base; }


    };




} // iterator
} // bliss
#endif /* TRANSFORM_ITERATOR_HPP_ */
