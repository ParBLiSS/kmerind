/**
 * memiterator.hpp
 *
 *  Created on: Feb 17, 2014
 *      Author: tpan
 */

#ifndef MEMITERATOR_HPP_
#define MEMITERATOR_HPP_

namespace bliss
{
  namespace iterator
  {

    /**
     *
     */
    template<typename T>
    class mem_iterator
      : public std::iterator<std::random_access_iterator_tag,
                              T,
                              std::ptrdiff_t,
                              typename std::add_pointer<T>::type,
                              typename std::add_rvalue_reference<T>::type
                              >
    {
      protected:
        T* _orig;
        T* _base;
        size_t _size;

      public:
        typedef mem_iterator< T >                                       type;


        typedef std::random_access_iterator_tag                         iterator_category;
        typedef T                                                       value_type;
        typedef std::ptrdiff_t                                          difference_type;
        typedef typename std::add_rvalue_reference<value_type>::type    reference_type;
        typedef typename std::add_pointer<value_type>::type             pointer_type;


        mem_iterator(T* base, const size_t& size) : _orig(base), _base(base), _size(size) {};
        mem_iterator(T* orig, T* base, const size_t& size) : _orig(orig), _base(base), _size(size) {};
        virtual ~mem_iterator() {};

        // for all classes
        // no EXPLICIT so can copy assign.
        mem_iterator(const type& Other) : _orig(Other._orig), _base(Other._base), _size(Other._size) {};

        type& operator=(const type& Other) {
          _size = Other._size;
          _base = Other._base;
          _orig = Other._orig;
          return *this;
        }

        type& operator++() {
          ++_base;
          return *this;
        }

        type operator++(int) {
          return type(_orig, _base++, _size);
        }


        // accessor functions for internal state.
        size_t& getSize() {
          return _size;
        }
        const size_t& getSize() const {
          return _size;
        }
        T*& getBaseIterator() {
          return _base;
        }
        const T*& getBaseIterator() const {
          return _base;
        }

        // input iterator specific
        inline bool operator==(const type& rhs)
        { return _base == rhs._base && _size == rhs._size && _orig == orig; }

        inline bool operator!=(const type& rhs)
        { return _base != rhs._base || _size != rhs._size || _orig != orig; }

        inline value_type operator*() {
          return *_base;
        }

        // no -> operator.  -> returns a pointer to the value held in iterator.
        // however, in transform iterator, we are not holding data, so pointer does not make sense.

        // NO support for output iterator


        // forward iterator
        explicit
        mem_iterator() : _size(), _base(), _origin() {};


       // bidirectional iterator
        type& operator--() {
         --_base;
         return *this;
        }

       type
       operator--(int) {
         return type(_orig, _base--, _size);
       }


       // random access iterator requirements.
       type
       operator+(difference_type n)
       {
         type other(*this);
         std::advance(other._base, n);
         return other;
       }

       friend
       type
       operator+(difference_type n, const type& right)
       {  return right + n; }

       type
       operator-(difference_type n)
       {
         type other(*this);
         std::advance(other._base, -n);
         return other;
       }

       difference_type
       operator-(const type& other)
       {
         return _base - other._base;
       }

       inline
       bool
       operator<(const type& rhs)
       { return _base < rhs._base; }

       inline
       bool
              operator>(const type& rhs)
       { return _base > rhs._base; }

       inline
       bool
              operator<=(const type& rhs)
       { return _base <= rhs._base; }

       inline
       bool
              operator>=(const type& rhs)
       { return _base >= rhs._base; }


       type&
       operator+=(difference_type n) {
         std::advance(_base, n);
         return *this;
       }

       type&
       operator-=(difference_type n)
       {
          std::advance(_base, -n);
          return *this;
       }

       value_type
       operator[](difference_type n) {
         return _base[n];
       }




    };

  } /* namespace iterator */
} /* namespace bliss */
#endif /* MEMITERATOR_HPP_ */
