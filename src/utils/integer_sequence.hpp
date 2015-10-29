/**
 * @file    integer_sequence.hpp
 * @ingroup utils
 * @author  Tony Pan <tpan7@gatech.edu>
 * @brief   templates to support generating a integer/index sequence
 * @details adapted from http://stackoverflow.com/questions/13072359/c11-compile-time-array-with-logarithmic-evaluation-depth/13073076#13073076
 *          this is logarithmic depth recursive metatemplate instantiation algorithm for creating an integer sequence.
 *
 *          useful because Clang only supports 256 levels of recursion.
 *
 *          generalized to use arbitrary primitive type, and therefore also using variadic template types.
 *          NOT YET generalized to allow arbitrary starting offset.
 */
#ifndef INTEGER_SEQUENCE_HPP_
#define INTEGER_SEQUENCE_HPP_

#include <type_traits>    // for is_integral
#include <cstdint>

namespace bliss
{
  namespace utils
  {

    // using aliases for cleaner syntax
    template<class T> using Invoke = typename T::type;

    template<class S1, class S2> struct concat;                    // TCP: generic concat template

    template<class S1, class S2>
    using Concat = Invoke<concat<S1, S2> >;                        // TCP: alias for convenience - function like.


    template<typename T, T...> struct seq { using type = seq; };    // TCP: "variable" holding a sequence of values. values may not be used.


    template<typename T, T... I1, T... I2>                         // TCP: template specialization of concat that does some real work
    struct concat<seq<T, I1...>, seq<T, I2...> >
      : seq<T, I1..., (sizeof...(I1) + I2)...> {};                 // TCP: note that the values in I2 are offset by size of I1. limitation: I1 and I2 need to have same starting value.

    template<typename T, T N > struct gen_seq;                      // TCP: forward declare.  hard coded starting value for the sequence to be 0.
    template<typename T, T N > using GenSeq = Invoke<gen_seq<T, N > >;

    //========= type specializations for size_t
    template<size_t N> struct gen_seq<size_t, N > : Concat<GenSeq<size_t, N/2>, GenSeq<size_t, N - N/2> >{}; // TCP: actual recursive template instance.
    template<> struct gen_seq<size_t, 0 > : seq<size_t>{};                               // TCP: specializations to handle recursion termination.
    template<> struct gen_seq<size_t, 1 > : seq<size_t, 0>{};

    //========= type specializations for int64
    template<int64_t N> struct gen_seq<int64_t, N > : Concat<GenSeq<int64_t, N/2>, GenSeq<int64_t, N - N/2> >{}; // TCP: actual recursive template instance.
    template<> struct gen_seq<int64_t, 0 > : seq<int64_t>{};                               // TCP: specializations to handle recursion termination.
    template<> struct gen_seq<int64_t, 1 > : seq<int64_t, 0>{};

    //========= type specializations for uint8_t
    template<uint8_t N> struct gen_seq<uint8_t, N > : Concat<GenSeq<uint8_t, N/2>, GenSeq<uint8_t, N - N/2> >{}; // TCP: actual recursive template instance.
    template<> struct gen_seq<uint8_t, 0 > : seq<uint8_t>{};                               // TCP: specializations to handle recursion termination.
    template<> struct gen_seq<uint8_t, 1 > : seq<uint8_t, 0>{};

    //========= type specializations for int8_t
    template<int8_t N> struct gen_seq<int8_t, N > : Concat<GenSeq<int8_t, N/2>, GenSeq<int8_t, N - N/2> >{}; // TCP: actual recursive template instance.
    template<> struct gen_seq<int8_t, 0 > : seq<int8_t>{};                               // TCP: specializations to handle recursion termination.
    template<> struct gen_seq<int8_t, 1 > : seq<int8_t, 0>{};

    // IF offset is needed, compute it when using the integer sequence.


    /**
     * @class index_sequence
     * @brief define a template specialized type for index sequence, where the value data type is size_t.
     * @tparam I    variadic list of index values of type size_t
     */
    template<std::size_t... I>
      using index_sequence = seq<std::size_t, I...>;

    /**
     * @typedef make_integer_sequence
     * @brief   make_integer_sequence<T, N> is an alias for integer_sequence<T, 0,...N-1>
     * @tparam T  value type
     * @tparam N  number of entries
     * @tparam O  offset of the first value.
     */
    template<typename T, T N>
      using make_integer_sequence = gen_seq<T, N >; // typename detail::iota<T, N, N>::type;

    /**
     * @typedef make_index_sequence
     * @brief   specialization of integer sequence for index. type is set to size_t
     * @tparam N  number of entries
     */
    template<size_t N>
      using make_index_sequence = make_integer_sequence<std::size_t, N>;


    /**
     * @typedef index_sequence_for
     * @brief   create a integer index version of a sequence of types.
     * @details index_sequence_for<A, B, C> is an alias for index_sequence<0, 1, 2>
     * @tparam Args  variadic list of  type template arguments
     */
    template<typename... Args>
      using index_sequence_for = make_index_sequence<sizeof...(Args)>;


  } // namespace utils

} // namespace bliss



#endif /* INTEGER_SEQUENCE_HPP_ */
