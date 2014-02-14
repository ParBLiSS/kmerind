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
 *
 *  Relies on compiler to cast the input parameter to the specific types needed by the functor, including constness and references.
 *
 *     - no std::function - general feeling is that it's slower.
 *
 *  TODO:
 *  commenting.
 *  specialization for windowed iterator and block partitioned iterator.
 *    no need for cyclic partitioned iterator right now, and no need for skipping iterator (random access handles that)
 */

#ifndef TRANSFORM_ITERATOR_HPP_
#define TRANSFORM_ITERATOR_HPP_

#include <iterator>

namespace bliss
{

namespace iterator
{
  ///////////// functor and function pointer traits support.
  //
  // Supports:
  //   overloaded operator() and functions : via parameter types.
  //     for functors, fortunately we can construct the list from the Base_Iterator (our use case).
  //     for function pointers, (member, nonmember) these are already explicitly specified by the user.
  //   static function and member function pointers.
  //   for overloaded function, need to construct the function pointer, which requires specifying the exact function input parameters for decltype.
  //        use static_cast<RET (*)(ARGS)>
  //        else create a typedef for the function pointer type, assign the function to it (compiler resolve it) and then do decltype on that.
  //   works with lambda functions (except that class type may be main function.  don't know how this affects evaluation.
  //   functor input const/ref:  user/iterator specified arg types, which may not match correctly.  when calling, need to supply
  //       references instead, then let the function call cast to const, ref, const ref, or just make copy.
  //       result_of appears to be smart enough to deal with const/ref during decltype step.
  //  function returning const type:  this is considered not useful and obselete with c++11.  not tested here.
  //
  // Possible Issues:
  //   if operator() overloading is based on input parameter const/ref, then decltype inside result_of will not be able to resolve.  does not happen with func ptrs.
  //
  // Use:  'dereference' function is implemented in the func_traits class as static functions.
  //
  //


  template<typename F, typename... Args>
  struct func_traits
  {
      typedef decltype(std::declval<F>()(std::declval<typename std::add_lvalue_reference<Args>::type>()...)) return_type;
  };


  /**
   * transform iterator class
   *
   * transforms each element in the list. otherwise it's a standard random access iterator.
   *
   * inheriting from std::iterator ONLY to get iterator_traits support.
   *
   * tried to use enable_if.  problem is that false causes the "type" typedef inside not to be defined, creating an error.
   *   this may be a SFINAE usage or implementation issue.
   */
  template<typename Functor,
           typename Base_Iterator
           >
  class transform_iterator
    : public std::iterator<typename std::iterator_traits<Base_Iterator>::iterator_category,
                           typename std::remove_reference<typename func_traits<Functor, typename Base_Iterator::value_type>::return_type>::type,
                           typename std::iterator_traits<Base_Iterator>::difference_type,
                           typename std::add_pointer<typename std::remove_reference<typename func_traits<Functor, typename Base_Iterator::value_type>::return_type>::type>::type,
                           typename std::add_rvalue_reference<typename std::remove_reference<typename func_traits<Functor, typename Base_Iterator::value_type>::return_type>::type>::type
                           >
    {
      protected:
        // define first, to avoid -Wreorder error (where the variables are initialized before transform_iterator::Functor, etc are defined.
        typedef std::iterator_traits<Base_Iterator>                                                         base_traits;
        typedef func_traits<Functor, typename Base_Iterator::value_type>                                    functor_traits;

        Base_Iterator                                                                                       _base;
        Functor                                                                                             _f;

      public:
        typedef transform_iterator< Functor, Base_Iterator >                                                type;
        typedef std::iterator_traits<type>                                                                  traits;

        typedef typename base_traits::iterator_category                                                     iterator_category;
        typedef typename std::remove_reference<typename functor_traits::return_type>::type                  value_type;
        typedef typename base_traits::difference_type                                                       difference_type;
        typedef typename std::add_rvalue_reference<value_type>::type                                        reference_type;
        typedef typename std::add_pointer<value_type>::type                                                 pointer_type;

        typedef typename base_traits::value_type                                                            base_value_type;


        // class specific constructor
        explicit
        transform_iterator(const Base_Iterator& base_iter, const Functor & f)
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
        Functor& getFunctor() {
          return _f;
        }
        const Functor& getFunctor() const {
          return _f;
        }
        Base_Iterator& getBaseIterator() {
          return _base;
        }
        const Base_Iterator& getBaseIterator() const {
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
         type other(*this);
         std::advance(other._base, n);
         return other;
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
         type other(*this);
         std::advance(other._base, -n);
         return other;
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
