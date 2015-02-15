/**
 * @file    integer_sequence.hpp
 * @ingroup bliss::utils
 * @author  Tony Pan <tpan7@gatech.edu>
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

#include <type_traits>    // for is_integral

namespace bliss
{
  namespace utils
  {

    /**
     * @class integer_sequence
     * @brief     A type that represents a parameter pack of zero or more integer values, each of type T.
     * @tparam T  the type of the values.  can only be integral type.
     * @tparam I  a variadic list of values.
     */
    template<typename T, T... I>
      struct integer_sequence
      {
        // T can only be integral
        static_assert( std::is_integral<T>::value, "Integral type" );

        /**
         * @typedef type
         * @brief type for the values in sequence.
         */
        using type = T;

        /**
         * @var   size
         * @brief number of values in the parameter pack.
         */
        static constexpr T size = sizeof...(I);

        /**
         *  @typedef    append
         *  @brief      Generate a new integer_sequence with an additional element of value N
         *  @details    Usage:  seq.append<4>
         *  @tparam N   a new value to append to current sequence, with type T.
         *
         */
        template<T N>
          using append = integer_sequence<T, I..., N>;

        /**
         * @typedef next
         * @brief   get the next sequence by appending the current size
         * @details repeated call creates a sequence with values 1, 2, 3, ...
         */
        using next = append<size>;
      };

    /*
     * define the size (declared and initialized in integer_sequence)
     */
    template<typename T, T... I>
      constexpr T integer_sequence<T, I...>::size;

    /**
     * @class index_sequence
     * @brief define a template specialized type for index sequence, where the value data type is size_t.
     * @tparam I    variadic list of index values of type size_t
     */
    template<std::size_t... I>
      using index_sequence = integer_sequence<std::size_t, I...>;

    /**
     * @namespace
     * @brief  namespace detail contains the logic for constructing a integer sequence at compile time.
     */
    namespace detail
    {
      /**
       *  @class iota
       *  @brief  Metafunction that generates an integer_sequence of T containing [0, N)
       *  @details recursively type substituted to construct the integer sequence.
       *  @tparam T   value's type
       *  @tparam Nt  values
       *  @tparam N   total number of elements (remaining)
       */
      template<typename T, T Nt, std::size_t N>
        struct iota
        {
          static_assert( Nt >= 0, "N cannot be negative" );

          using type = typename iota<T, Nt-1, N-1>::type::next;
        };

      /**
       *  @class iota
       *  @brief  Terminal Case.  Metafunction that generates an integer_sequence of T containing [0, N)
       *  @details Base class constructs an empty integer sequence (not even 0).
       *            from here, each recursions back up the stack adds one more value = size of current sequence.
       *  @tparam T   value's type
       *  @tparam Nt  total number of elements (remaining)
       */
      template<typename T, T Nt>
        struct iota<T, Nt, 0ul>
        {
          using type = integer_sequence<T>;
        };
    }


    /**
     * @typedef make_integer_sequence
     * @brief   make_integer_sequence<T, N> is an alias for integer_sequence<T, 0,...N-1>
     * @tparam T  value type
     * @tparam N  number of entries
     */
    template<typename T, T N>
      using make_integer_sequence = typename detail::iota<T, N, N>::type;

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
     * @tparams Args  variadic list of  type template arguments
     */
    template<typename... Args>
      using index_sequence_for = make_index_sequence<sizeof...(Args)>;


  } // namespace utils

} // namespace bliss



#endif /* INTEGER_SEQUENCE_HPP_ */
