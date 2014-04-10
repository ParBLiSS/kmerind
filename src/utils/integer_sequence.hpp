/**
 * @file		integer_sequence.hpp
 * @ingroup
 * @author	tpan
 * @brief   templates to support generating a integer/index sequence
 * @details adapted from http://coliru.stacked-crooked.com/a/663dee82db2aef0a
 *
 *
  // Copyright Jonathan Wakely 2012-2013
  // Distributed under the Boost Software License, Version 1.0.
  // (See accompanying file LICENSE_1_0.txt or copy at
  // http://www.boost.org/LICENSE_1_0.txt)
  // https://gitorious.org/redistd/integer_seq/raw/fdce85d18df4ea99240524c626d5452a0f7c3faf:integer_seq.h

  // A C++11 implementation of std::integer_sequence from C++14
 *
 */
#ifndef INTEGER_SEQUENCE_HPP_
#define INTEGER_SEQUENCE_HPP_

#include <type_traits>

namespace bliss
{
  namespace utils
  {

    /// A type that represents a parameter pack of zero or more integers.
    template<typename T, T... I>
      struct integer_sequence
      {
        static_assert( std::is_integral<T>::value, "Integral type" );

        using type = T;

        static constexpr T size = sizeof...(I);

        /// Generate an integer_sequence with an additional element.
        template<T N>
          using append = integer_sequence<T, I..., N>;

        using next = append<size>;
      };

    template<typename T, T... I>
      constexpr T integer_sequence<T, I...>::size;

    template<std::size_t... I>
      using index_sequence = integer_sequence<std::size_t, I...>;

    namespace detail
    {
      // Metafunction that generates an integer_sequence of T containing [0, N)
      template<typename T, T Nt, std::size_t N>
        struct iota
        {
          static_assert( Nt >= 0, "N cannot be negative" );

          using type = typename iota<T, Nt-1, N-1>::type::next;
        };

      // Terminal case of the recursive metafunction.
      template<typename T, T Nt>
        struct iota<T, Nt, 0ul>
        {
          using type = integer_sequence<T>;
        };
    }


    // make_integer_sequence<T, N> is an alias for integer_sequence<T, 0,...N-1>
    template<typename T, T N>
      using make_integer_sequence = typename detail::iota<T, N, N>::type;

    template<int N>
      using make_index_sequence = make_integer_sequence<std::size_t, N>;


    // index_sequence_for<A, B, C> is an alias for index_sequence<0, 1, 2>
    template<typename... Args>
      using index_sequence_for = make_index_sequence<sizeof...(Args)>;


  } // namespace utils

} // namespace bliss



#endif /* INTEGER_SEQUENCE_HPP_ */
