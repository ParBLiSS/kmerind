/**
 * @file		container_traits.hpp
 * @ingroup
 * @author	tpan
 * @brief
 * @details modified from solution here: http://stackoverflow.com/questions/9407367/determine-if-a-type-is-an-stl-container-at-compile-time,
 *            which is in turn from here: http://stackoverflow.com/questions/4850473/pretty-print-c-stl-containers
 *
 *            updated to use C++ 11 stuff.
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef CONTAINER_TRAITS_HPP_
#define CONTAINER_TRAITS_HPP_


template<typename T>
struct has_const_iterator
{
private:
    typedef char                      yes;
    typedef struct { char array[2]; } no;

    template<typename C> static yes test(typename C::const_iterator*);
    template<typename C> static no  test(...);
public:
    static const bool value = sizeof(test<T>(0)) == sizeof(yes);
    typedef T type;
};

template <typename T>
struct has_begin_end
{
    template<typename C> static char (&f(typename std::enable_if<
      std::is_same<decltype(static_cast<typename C::const_iterator (C::*)() const>(&C::begin)),
      typename C::const_iterator(C::*)() const>::value, void>::type*))[1];

    template<typename C> static char (&f(...))[2];

    template<typename C> static char (&g(typename std::enable_if<
      std::is_same<decltype(static_cast<typename C::const_iterator (C::*)() const>(&C::end)),
      typename C::const_iterator(C::*)() const>::value, void>::type*))[1];

    template<typename C> static char (&g(...))[2];

    static bool const beg_value = sizeof(f<T>(0)) == 1;
    static bool const end_value = sizeof(g<T>(0)) == 1;
};


template<typename T>
struct is_container : std::integral_constant<bool, has_const_iterator<T>::value && has_begin_end<T>::beg_value && has_begin_end<T>::end_value>
{ };


/////// check if "assign" is available.
///// adopted from  http://stackoverflow.com/questions/9530928/checking-a-member-exists-possibly-in-a-base-class-c11-version

template <typename T, typename InputIterator>
constexpr auto has_assign_method(T*) -> decltype(std::declval<T>().assign(std::declval<InputIterator>(), std::declval<InputIterator>()), bool())
{
    return true;
}
constexpr bool has_assign_method(...) {return false;}

//
//template <typename, typename T>
//struct has_assign
//{
//    static_assert(std::integral_constant<T, false>::value,
//      "Second template parameter needs to be of function type.");
//};
//
//template<typename C, typename Ret, typename... Args>
//struct has_assign<C, Ret(Args...)> {
//  private:
//      template<typename T>
//      static constexpr auto check(T*)
//      -> typename
//          std::is_same<
//              decltype( std::declval<T>().assign( std::declval<Args>()... ) ),
//              Ret    // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
//          >::type;  // attempt to call it and see if the return type is correct
//
//      template<typename>
//      static constexpr std::false_type check(...);
//
//      typedef decltype(check<C>(0)) type;
//
//  public:
//      static constexpr bool value = type::value;
//};
//
//template<typename T, typename InputIterator>
//struct is_sequence_container : std::integral_constant<bool, has_assign<T, void(InputIterator, InputIterator) >::value>
//{};

template<typename T>
struct get_const_iterator_type {
    typedef typename std::enable_if<is_container<T>::value, typename T::const_iterator>::type type;
};

#endif /* CONTAINER_TRAITS_HPP_ */
