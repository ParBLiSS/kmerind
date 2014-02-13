/**
 * transform_iterator.hpp
 *
 *  Created on: Feb 5, 2014
 *      Author: tpan
 *
 *  patterned after boost's transform iterator.  no multiple inheritance.  uses functor class.  allows variadic parameters
 *  support functor (with a single operator() declaration), function pointer, and member function pointer.
 *    in the last 2 cases, need to be specific about the argument types.  given a function ptr type, we deduce the argument types,
 *      but we do not check that the types from the base iterator is matching.
 *      suggest:  typedef Ret (Class::*fptr_t)(Args...);  fptr_t = &Class::func
 *    in the first case, functor should have a single operator(), with the matching parameters.  we use the iterator value type
 *      and assume the function is operator(), and get the return type.
 *
 *     - no std::function - general feeling is that it's slower.
 *
 *  TODO:
 *  more testing
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
  ///////////// functor and function traits support.  These rely on
  //   1. no default values: if Is_Class default value is supplied, then we end up with template parameters that are ambiguously mapped to multiple templating.
  //        default value also causes issues with variadic parameters, and partial specialization.
  //   2. specialization: an matching functor_traits specialization is then chosen base on T and Is_Class.  else fall back to default which is declared but not defined..
  //   3. type deduction: the specialization uses T, which is represented as a composition of the template parameters of the specialization.
  //        i.e. specialization of the generic template is itself templated with a different set of parameters.
  //        the specialization's own template parameters are completely deduced from T.
  //   4. use variadic template parameter for function's parameter list.  this reduces the number of specializations significantly.
  // the deduction then allows us to get the class and return types, even the types of function arguments.
  //
  //
  // Possible Issues:
  //   deal with static function?
  //   may not do well with const functions?
  //   may not work with lambda functions?
  //   does not work with overloaded operator().  return type is deduced incorrectly.  DO NOT DO OVERLOADED FUNCTION
  //
  // Use:  'dereference' function is implemented in the func_traits class as static functions.
  //
  // To deal with possible overloaded operators in functors, we need to specify parameters.
  //  for functors, fortunately we can construct the list from the Base_Iterator (our use case).
  //  for function pointers, (member, nonmember) these are already explicitly specified by the user.
  //
  // **BEST PRACTICE:  user should implement functors to have a single operator() (overloading operator() negates the simplicity of using a functor)
  //                   ELSE specify the exact function as a member function pointer.
  //
  // ***variadic templates is hard to get right.  see http://pic.dhe.ibm.com/infocenter/ratdevz/v8r0/index.jsp?topic=%2Fcom.ibm.xlcpp111.aix.doc%2Flanguage_ref%2Fvariadic_templates.html
  // especially when mixing 2 variadic template list (first one needs to be deducible.) variadic parameters normally should be last unless deducible.
  //

  // unspecialized functor traits class.
  // don't default init Is_Class.  else the generic template will look like the functor specialization version.
  template<typename F, bool Is_Class, typename... Args>
  struct func_traits;

  // functor object's specialization
  template<typename Functor, typename... Args>
  struct func_traits<Functor, true, Args...> {
      typedef std::true_type has_class;
      typedef Functor class_type;
      typedef typename std::result_of<Functor(Args...)>::type return_type;

      //   have class, even though not used.
      static inline return_type eval(Functor& f, class_type& c, Args... inputs) {
        return f(inputs...);
      }
  };

  // member function pointer's specialization.  DArgs are deduced.
  template<typename Class, typename Ret, typename... Args>
  struct func_traits<Ret (Class::*)(Args...), false, Args...>
  {
      typedef std::true_type has_class;
      typedef Class class_type;
      typedef Ret return_type;

      // get back the functor type
      typedef Ret (Class::*Functor)(Args...);
      static inline return_type eval(Functor& f, class_type& c, Args... inputs) {
        return (c.*f)(inputs...);
      }
  };

  // function pointer's specialization.
  template<typename Ret, typename... Args>
  struct func_traits<Ret (*)(Args...), false, Args...>
  {
      typedef std::false_type has_class;
      typedef std::nullptr_t class_type;
      typedef Ret return_type;

      // get back the functor type.  have class, even though not used.
      typedef Ret (*Functor)(Args...);
      static inline return_type eval(Functor& f, class_type& c, Args... inputs) {
        return (*f)(inputs...);
      }
  };


  /**
   * transform iterator class
   *
   * will deduce from Functor template param whether to specialize for Functor, Function Pointer, or Member Function Pointer.
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
                           typename std::remove_reference<typename func_traits<Functor, std::is_class<Functor>::value, typename Base_Iterator::value_type>::return_type>::type,
                           typename std::iterator_traits<Base_Iterator>::difference_type,
                           typename std::add_pointer<typename std::remove_reference<typename func_traits<Functor, std::is_class<Functor>::value, typename Base_Iterator::value_type>::return_type>::type>::type,
                           typename std::add_rvalue_reference<typename std::remove_reference<typename func_traits<Functor, std::is_class<Functor>::value, typename Base_Iterator::value_type>::return_type>::type>::type
                           >
    {
      protected:
        // define first, to avoid -Wreorder error (where the variables are initialized before transform_iterator::Functor, etc are defined.
        typedef std::iterator_traits<Base_Iterator>                                                         base_traits;
        typedef func_traits<Functor, std::is_class<Functor>::value, typename Base_Iterator::value_type>     functor_traits;

        Base_Iterator                                                                                       _base;
        Functor                                                                                             _f;
        typename functor_traits::class_type                                                                 _c;

      public:
        typedef transform_iterator< Functor, Base_Iterator >                                                type;
        typedef std::iterator_traits<type>                                                                  traits;

        typedef typename base_traits::iterator_category                                                     iterator_category;
        typedef typename std::remove_reference<typename functor_traits::return_type>::type                  value_type;
        typedef typename base_traits::difference_type                                                       difference_type;
        typedef typename std::add_rvalue_reference<value_type>::type                                        reference_type;
        typedef typename std::add_pointer<value_type>::type                                                 pointer_type;

        typedef typename base_traits::value_type                                                            base_value_type;


        explicit
        transform_iterator(const Base_Iterator& base_iter, const Functor & f)
          : _base(base_iter), _f(f), _c() {};

        explicit
        transform_iterator(const Base_Iterator& base_iter, const Functor & f,
                           typename std::add_lvalue_reference<
                             typename std::add_const<
                               typename std::remove_reference<
                                 typename func_traits<Functor,
                                                      std::is_class<Functor>::value,
                                                      typename Base_Iterator::value_type
                                                      >::class_type
                               >::type
                             >::type
                           >::type c)
          : _base(base_iter), _f(f), _c(c) {};

        explicit
        transform_iterator() : _f(), _base(), _c() {};

        // no EXPLICIT so can copy assign.
        transform_iterator(const type& Other) : _base(Other._base), _f(Other._f), _c(Other._c) {};

        type& operator=(const type& Other) {
          _f = Other._f;
          _base = Other._base;
          _c = Other._c;
          return *this;
        }

       Functor& getFunctor() {
         return _f;
       }
       const Functor& getFunctor() const {
         return _f;
       }

       typename std::add_rvalue_reference<typename functor_traits::class_type>::type& getClass() {
         return _c;
       }
       typename std::add_rvalue_reference<typename std::add_const<typename functor_traits::class_type>::type>::type& getClass() const {
         return _c;
       }

       Base_Iterator& getBaseIterator() {
         return _base;
       }
       const Base_Iterator& getBaseIterator() const {
         return _base;
       }

       // forward iterator:

       // no -> operator

       // this way works, but not with "const" on this function.
       inline value_type operator*() {
//         return functor_traits::eval(_f, _c, *_base);
         return functor_traits::eval(_f, _c, _base[0]);
       }

       type& operator++() {
         ++_base;
         return *this;
       }

       type operator++(int) {
         return type(_base++, _f, _c);
       }

       // bidirectional iterator
       type& operator--() {
         --_base;
         return *this;
       }

       type operator--(int) {
         return type(_base--, _f, _c);
       }

       // random access iterator requirements.
       inline value_type operator[](difference_type n) {
         return functor_traits::eval(_f, _c, _base[n]);
       }


       type& operator+=(difference_type n) {
         std::advance(_base, n);
         return *this;
       }

       type operator+(difference_type n)
       {
         type other(*this);
         std::advance(other._base, n);
         return other;
       }

       friend type operator+(difference_type n, const type& right)
       {  return right + n; }

       type& operator-=(difference_type n)
       {
          std::advance(_base, -n);
          return *this;
       }

       type operator-(difference_type n)
       {
         type other(*this);
         std::advance(other._base, -n);
         return other;
       }

       // Forward iterator requirements
       inline bool operator==(const type& rhs)
       { return _base == rhs._base; }

       inline bool operator!=(const type& rhs)
       { return _base != rhs._base; }

       // Random access iterator requirements
       inline bool operator<(const type& rhs)
       { return _base < rhs._base; }

       inline bool operator>(const type& rhs)
       { return _base > rhs._base; }

       inline bool operator<=(const type& rhs)
       { return _base <= rhs._base; }

       inline bool operator>=(const type& rhs)
       { return _base >= rhs._base; }

    };


} // iterator
} // bliss
#endif /* TRANSFORM_ITERATOR_HPP_ */
