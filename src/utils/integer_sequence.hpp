/*
 * Copyright 2015 Georgia Institute of Technology
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

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
 *          allows arbitrary offset via Offset template function.
 *
 *          this is instead of std::integer_sequence and std::make_integer_sequence in <utility>, which have linear depth instead of logarithmic.
 */
#ifndef INTEGER_SEQUENCE_HPP_
#define INTEGER_SEQUENCE_HPP_

#include <type_traits>    // for is_integral
#include <cstdint>
#include <tuple>

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


    /**
     * @typedef make_integer_sequence
     * @brief   make_integer_sequence<T, N> is an alias for integer_sequence<T, 0,...N-1>
     * @tparam T  value type
     * @tparam N  number of entries
     * @tparam O  offset of the first value.
     */
    template<typename T, T N>
      using make_integer_sequence = gen_seq<T, N >; // typename detail::iota<T, N, N>::type;


    //================= make offset integer sequence.
    /**
     * @typedef offset
     * @brief generic template declaration.  see specialization below.
     */
    template<class S1, class T, T O> struct offset;                    // TCP: generic concat template

    /**
     * @typedef Offset
     * @brief generic template declaration.  see specialization below.  make it function like
     */
    template<class S1, class T, T O> using Offset = Invoke<offset<S1, T, O> >;   // TCP: alias for convenience - function like.

    /**
     * @typedef offset
     * @brief   specialization to add an offset to all elements of an integer sequence type
     * @tparam T index value type
     * @tparam O offset to add to the integer sequence type
     * @tparam I values of original indices.
     */
    template<typename T, T OFFSET, T... I>      // TCP: template specialization of offset that does some real work
    struct offset<seq<T, I...>, T, OFFSET>
      : seq<T, (OFFSET + I)...> {};             // TCP: note that the values in I are now offset by OFFSET.


    // ================= index sequences

    /**
     * @typedef make_index_sequence
     * @brief   specialization of integer sequence for index. type is set to size_t
     * @tparam N  number of entries
     */
    template<size_t N>
      using make_index_sequence = make_integer_sequence<std::size_t, N>;

    /**
     * @class index_sequence
     * @brief define a template specialized type for index sequence, where the value data type is size_t.
     * @tparam I    variadic list of index values of type size_t
     */
    template<size_t... I>
      using index_sequence = seq<std::size_t, I...>;


    // ================== variadic param pack utils

    /**
     * @typedef index_sequence_for
     * @brief   create a integer index version of a sequence of types.
     * @details index_sequence_for<A, B, C> is an alias for index_sequence<0, 1, 2>
     * @tparam Args  variadic list of  type template arguments
     */
    template<typename... Args>
      using index_sequence_for = make_index_sequence<sizeof...(Args)>;


    /// generic select template
    template <typename Index, typename... Args> struct select;
    /// generic select template, emulating function
    template <typename Index, typename... Args> using Select = Invoke<select<Index, Args...> >;

    /**
     * @typedef select
     * @brief   select parameters from variadic param pack based on an index sequence.
     * @detail  might be possible to do this better, rather than using tuple over and over.
     * @tparam I     indices to be selected in a seq struct
     * @tparam Args  input types to be selected.
     */
     template <std::size_t... I, typename... Args>
     struct select<seq<std::size_t, I...>, Args...> :
       std::tuple<std::tuple_element<I, std::tuple<Args...> >... > {};






    // =============== for tuples.
    /// generic select template
    template <typename Index, typename Tuple> struct tuple_elements;


     /**
      * @typedef tuple_elements
      * @brief variadic extension to std::tuple_element.  select specific tuple elements and return as a tuple of selected element types.
      * @tparam I     indices to be selected in a seq struct
      * @tparam Args  input types to be selected.
      */
     template <std::size_t... I, typename... Args>
     struct tuple_elements< ::bliss::utils::seq<std::size_t, I...>, std::tuple<Args...> >
       : Select<::bliss::utils::seq<std::size_t, I...>, Args...> {};

     /**
      * @brief get element references from an instantiated tuple based on an index sequence
      * @tparam I  indices to select
      * @tparam T  tuple element types.
      */
     template < typename Seq, typename... Args>
     tuple_elements< Seq, std::tuple<Args...> > get(std::tuple<Args...> & vars) {
       //WIP, TODO

     }





  } // namespace utils

} // namespace bliss



#endif /* INTEGER_SEQUENCE_HPP_ */
