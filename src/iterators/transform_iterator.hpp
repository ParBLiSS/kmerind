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
 *  specialization for function pointers (can't use multiple inheritance) - use std::function - general feeling is that it's slower.
 */

#ifndef TRANSFORM_ITERATOR_HPP_
#define TRANSFORM_ITERATOR_HPP_

#include <iterator>

namespace bliss
{

namespace iterator
{
  ///////////// functor and function traits support.  These rely on
  //   1. default values: functor_traits<T> causes Is_Class to be evaluated.
  //   2. specialization: an matching functor_traits specialization is then chosen base on T and Is_Class.  else fall back to default and generate assert.
  //   3. type deduction: the specialization uses T, which is represented as a composition of the template parameters of the specialization.
  //        i.e. specialization of the generic template is itself templated with a different set of parameters.
  //        the specialization's own template parameters are completely deduced from T.
  // the deduction then allows us to get the class and return types, even the types of function arguments.
  //
  // Possible Issues:
  //   does not deal with static function
  //   may not do well with const functions.
  //   may not work with lambda functions
  //   does not deal with overloaded operator()
  //
  // Use:  specialized implementation of 'dereference' based on function type (functor, function pointer, or member function pointer)

  // TODO: to deal with possible overloaded operators in functors, we need to specify parameters.
  //  for functors, fortunately we can construct the list from the Base_Iterator (our use case).
  //  for function pointers, (member, nonmember) these are already explicitly specified by the user.

  // **BEST PRACTICE:  user should implement functors to have a single operator() (overloading operator() negates the simplicity of using a functor)
  //                   ELSE specify the exact function as a member function pointer.




  // unspecialized functor traits class.  throws a static assert error.
  template<typename F, bool Is_Class>
  struct func_traits {
      static_assert(std::is_class<F>::value || std::is_function<F>::value || std::is_member_function_pointer<F>::value,
                    "functor_traits expects template parameter of type class/struct with operator(), function pointer, or member function pointer\n");
  };

  // functor object's specialization
  template<typename Functor>
  struct func_traits<Functor, true> {
      typedef std::true_type has_class;
      typedef Functor class_type;
      typedef typename func_traits<decltype(&Functor::operator()), false>::return_type return_type;
  };

  // member function pointer's specialization
  template<typename Class,typename Ret,typename... Args>
  struct func_traits<Ret (Class::*)(Args...), false>
  {
      typedef std::true_type has_class;
      typedef Class class_type;
      typedef Ret return_type;
  };

  // function pointer's specialization.
  template<typename Ret,typename... Args>
  struct func_traits<Ret (*)(Args...), false>
  {
      typedef std::false_type has_class;
      typedef Ret return_type;
  };


  /////////////////////
  // wrapper struct for types.  specialize based on whether F is a functor/member pointer (has Class), or a function pointer (no class).
  template<typename F, bool Has_Class = func_traits<F, std::is_class<F>::value>::has_class::value>
  struct functor_traits;

  template<typename F>
  struct functor_traits<F, true> {
      typedef func_traits<F, std::is_class<F>::value> trait_type;
      typedef typename trait_type::return_type return_type;
      typedef typename trait_type::class_type class_type;
  };

  template<typename F>
  struct functor_traits<F, false> {
      typedef func_traits<F, std::is_class<F>::value> trait_type;
      typedef typename trait_type::return_type return_type;
  };




  // no Value.  can't have default value for partial specialization, so no point in having it.
  // also causes instantiation problem - Value without default here causes 'provided' error, and value with default causes incomplete type error.
  template<typename Functor,
           typename Base_Iterator,
           bool Is_Class = std::is_class<Functor>::value,
           bool Is_Member_Function = std::is_member_function_pointer<Functor>::value
           >
  class transform_iterator;


  template<typename Functor,
           typename Base_Iterator
           >
  class transform_iterator<Functor, Base_Iterator, true, false>
    : public std::iterator<typename std::iterator_traits<Base_Iterator>::iterator_category,
                           typename std::remove_reference<typename std::result_of<Functor(typename Base_Iterator::value_type)>::type>::type,
                           typename std::iterator_traits<Base_Iterator>::difference_type,
                           typename std::add_pointer<typename std::remove_reference<typename std::result_of<Functor(typename Base_Iterator::value_type)>::type>::type>::type,
                           typename std::add_rvalue_reference<typename std::remove_reference<typename std::result_of<Functor(typename Base_Iterator::value_type)>::type>::type>::type
                           >
    {
      protected:
        // define first, to avoid -Wreorder error (where the variables are initialized before transform_iterator::Functor, etc are defined.
        Base_Iterator                                                   _base;
        Functor                                                         _f;
        typedef std::iterator_traits<Base_Iterator>                     base_traits;

      public:
        typedef transform_iterator<Functor,
            Base_Iterator,
            true, false
            >              type;
        typedef Base_Iterator                                           base_type;
        typedef Functor                                                 functor_type;
        typedef std::iterator_traits<type>                              traits;

        typedef typename base_traits::iterator_category                 iterator_category;
        typedef typename std::remove_reference<typename std::result_of<Functor(typename Base_Iterator::value_type)>::type>::type                                                   value_type;
        typedef typename base_traits::difference_type                   difference_type;
        typedef typename std::add_rvalue_reference<value_type>::type    reference_type;
        typedef typename std::add_pointer<value_type>::type             pointer_type;

        typedef typename base_traits::value_type                        base_value_type;



        explicit
        transform_iterator(const Base_Iterator& base_iter, const Functor & f)
          : _base(base_iter), _f(f) {};

        explicit
        transform_iterator() : _f(), _base() {};

        // no EXPLICIT so can copy assign.
        transform_iterator(const type& Other) : _base(Other.getBaseIterator()), _f(Other.getFunctor()) {};

        type& operator=(const type& Other) {
          this->_f = Other.getFunctor();
          this->_base = Other.getBaseIterator();
          return *this;
        }

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
         // Functor object
         //return Functor::operator()(*_base);
         return _f(*_base);  // this works now.  why?
       }
//       typename std::enable_if<std::is_function<Functor>::value, value_type>::type operator*() {
//         // normal function
//         return (*_f)(*_base);
//       }
       // no -> operator

       type& operator++() {
         ++_base;
         return *this;
       }

       type operator++(int) {
         return type(_base++, _f);
       }

       // bidirectional iterator
       type& operator--() {
         --_base;
         return *this;
       }

       type operator--(int) {
         return type(_base--, _f);
       }

       // random access iterator requirements.
       value_type operator[](difference_type n) {
         return _f(_base[n]);
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





       inline friend bool operator==(const type& lhs,
            const type& rhs)
       { return lhs._base == rhs._base; }

       inline friend bool operator!=(const type& lhs,
            const type& rhs)
       { return lhs._base != rhs._base;  }

       // Random access iterator requirements
       inline friend bool operator<(const type& lhs,
           const type& rhs)
       { return lhs._base < rhs._base; }

       inline friend bool operator>(const type& lhs,
           const type& rhs)
       { return lhs._base > rhs._base; }

       inline friend bool operator<=(const type& lhs,
            const type& rhs)
       { return lhs._base <= rhs._base; }

       inline friend bool operator>=(const type& lhs,
            const type& rhs)
       { return lhs._base >= rhs._base; }


    };


  // Functor type passed in has to be the reference of function.  also, is class and is function will both return false.  only is_constructible returns true.
  // so need the double negative checks.
  // need to use &func as function pointer to initialize.  need to use *func to call.
  template<typename Functor,
           typename Base_Iterator
           >
  class transform_iterator<Functor, Base_Iterator,
    false, false>
  : public std::iterator<typename std::iterator_traits<Base_Iterator>::iterator_category,
                         typename std::remove_reference<typename std::result_of< Functor(typename Base_Iterator::value_type)>::type>::type,
                         typename std::iterator_traits<Base_Iterator>::difference_type,
                         typename std::add_pointer<typename std::remove_reference<typename std::result_of<Functor(typename Base_Iterator::value_type)>::type>::type>::type,
                         typename std::add_rvalue_reference<typename std::remove_reference<typename std::result_of<Functor(typename Base_Iterator::value_type)>::type>::type>::type
                         >
      {
        protected:
          // define first, to avoid -Wreorder error (where the variables are initialized before transform_iterator::Functor, etc are defined.
          Base_Iterator                                                   _base;
          Functor                                                         _f;
          typedef std::iterator_traits<Base_Iterator>                     base_traits;

        public:
          typedef transform_iterator<Functor, Base_Iterator,
              false, false>                                                           type;
          typedef Base_Iterator                                           base_type;
          typedef Functor                                                 functor_type;
          typedef std::iterator_traits<type>                              traits;

          typedef typename base_traits::iterator_category                 iterator_category;
          typedef typename std::remove_reference<typename std::result_of<Functor(typename Base_Iterator::value_type)>::type>::type   value_type;
          typedef typename base_traits::difference_type                   difference_type;
          typedef typename std::add_rvalue_reference<value_type>::type    reference_type;
          typedef typename std::add_pointer<value_type>::type             pointer_type;

          typedef typename base_traits::value_type                        base_value_type;



          explicit
          transform_iterator(const Base_Iterator& base_iter, const Functor & f)
            : _base(base_iter) { _f = f; };

          explicit
          transform_iterator() : _f(), _base() {};

          // let COMPILER provide a copy constructor.
          transform_iterator(const type& Other) : _base(Other.getBaseIterator()) { _f = Other.getFunctor(); };

          type& operator=(const type& Other) {
            this->_f = Other.getFunctor();
            this->_base = Other.getBaseIterator();
            return *this;
          }




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

        // this way works, but not with "const" on the function.
         value_type operator*() {
           // Functor object
           //return Functor::operator()(*_base);
           return (*_f)(*_base);  // this works now.  why?
         }
  //       typename std::enable_if<std::is_function<Functor>::value, value_type>::type operator*() {
  //         // normal function
  //         return (*_f)(*_base);
  //       }

         // no -> operator

         type& operator++() {
           ++_base;
           return *this;
         }

         type operator++(int) {
           return type(_base++, _f);
         }

         // bidirectional iterator
         type& operator--() {
           --_base;
           return *this;
         }

         type operator--(int) {
           return type(_base--, _f);
         }

         // random access iterator requirements.
         value_type operator[](difference_type n) {
           return (*_f)(_base[n]);
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




         inline friend bool operator==(const type& lhs,
              const type& rhs)
         { return lhs._base == rhs._base; }

         inline friend bool operator!=(const type& lhs,
              const type& rhs)
         { return lhs._base != rhs._base;  }

         // Random access iterator requirements
         inline friend bool operator<(const type& lhs,
             const type& rhs)
         { return lhs._base < rhs._base; }

         inline friend bool operator>(const type& lhs,
             const type& rhs)
         { return lhs._base > rhs._base; }

         inline friend bool operator<=(const type& lhs,
              const type& rhs)
         { return lhs._base <= rhs._base; }

         inline friend bool operator>=(const type& lhs,
              const type& rhs)
         { return lhs._base >= rhs._base; }


      };



  // member function pointer passed in via reference, and also need an instance of the object. invoke by c.*f
  template<typename Functor,
           typename Base_Iterator
           >
  class transform_iterator<Functor, Base_Iterator,
    false, true>
  : public std::iterator<typename std::iterator_traits<Base_Iterator>::iterator_category,
                         typename std::remove_reference<typename functor_traits<Functor>::return_type>::type,
                         typename std::iterator_traits<Base_Iterator>::difference_type,
                         typename std::add_pointer<typename std::remove_reference<typename functor_traits<Functor>::return_type>::type>::type,
                         typename std::add_rvalue_reference<typename std::remove_reference<typename functor_traits<Functor>::return_type>::type>::type
                         >
      {
        protected:
          // define first, to avoid -Wreorder error (where the variables are initialized before transform_iterator::Functor, etc are defined.
          Base_Iterator                                                   _base;
          Functor                                                         _f;
          typedef std::iterator_traits<Base_Iterator>                     base_traits;
          typedef typename functor_traits<Functor>::class_type         class_type;
          class_type                                                      _c;

        public:
          typedef transform_iterator<Functor, Base_Iterator,
              false, true>                                                           type;
          typedef Base_Iterator                                           base_type;
          typedef Functor                                                 functor_type;
          typedef std::iterator_traits<type>                              traits;

          typedef typename base_traits::iterator_category                 iterator_category;
          typedef typename std::remove_reference<typename functor_traits<Functor>::return_type >::type   value_type;
          typedef typename base_traits::difference_type                   difference_type;
          typedef typename std::add_rvalue_reference<value_type>::type    reference_type;
          typedef typename std::add_pointer<value_type>::type             pointer_type;

          typedef typename base_traits::value_type                        base_value_type;

          explicit
          transform_iterator(const Base_Iterator& base_iter, const class_type & c, const Functor & f)
            : _base(base_iter) { _f = f; _c = c; };

          explicit
          transform_iterator() : _f(), _c(), _base() {};

          // NO EXPLICIT else can't copy, copy assign, etc.
          transform_iterator(const type& Other) : _base(Other.getBaseIterator()) { _f = Other.getFunctorPtr(); _c = Other.getFunctorObject(); };


          type& operator=(const type& Other) {
            this->_f = Other.getFunctorPtr();
            this->_base = Other.getBaseIterator();
            this->_c = Other.getFunctorObject();
            return *this;
          }

         functor_type& getFunctorPtr() {
           return _f;
         }
         const functor_type& getFunctorPtr() const {
           return _f;
         }
         class_type& getFunctorObject() {
           return _c;
         }
         const class_type& getFunctorObject() const {
           return _c;
         }
         base_type& getBaseIterator() {
           return _base;
         }

         const base_type& getBaseIterator() const {
           return _base;
         }

         // forward iterator:

         // can't seem to get this way to work.
  //       OT operator*(type const & iter) {
  //         return _f.operator()(*(iter.getBaseIterator()));
  //       }
         // this way works, but not with "const" on the function.
         value_type operator*() {
           // Functor object
           //return Functor::operator()(*_base);
           return (_c.*_f)(*_base);  // this works now.  why?
         }
  //       typename std::enable_if<std::is_function<Functor>::value, value_type>::type operator*() {
  //         // normal function
  //         return (*_f)(*_base);
  //       }
         // no -> operator

         type& operator++() {
           ++_base;
           return *this;
         }

         type operator++(int) {
           return type(_base++, _c, _f);
         }

         // bidirectional iterator
         type& operator--() {
           --_base;
           return *this;
         }

         type operator--(int) {
           return type(_base--, _c, _f);
         }

         // random access iterator requirements.
         value_type operator[](difference_type n) {
           return (_c.*_f)(_base[n]);
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

         inline friend bool operator==(const type& lhs,
              const type& rhs)
         { return lhs._base == rhs._base; }

         inline friend bool operator!=(const type& lhs,
              const type& rhs)
         { return lhs._base != rhs._base;  }

         // Random access iterator requirements
         inline friend bool operator<(const type& lhs,
             const type& rhs)
         { return lhs._base < rhs._base; }

         inline friend bool operator>(const type& lhs,
             const type& rhs)
         { return lhs._base > rhs._base; }

         inline friend bool operator<=(const type& lhs,
              const type& rhs)
         { return lhs._base <= rhs._base; }

         inline friend bool operator>=(const type& lhs,
              const type& rhs)
         { return lhs._base >= rhs._base; }


      };

} // iterator
} // bliss
#endif /* TRANSFORM_ITERATOR_HPP_ */
