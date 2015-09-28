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
 */
#ifndef INTEGER_SEQUENCE_HPP_
#define INTEGER_SEQUENCE_HPP_

#include <type_traits>    // for is_integral

namespace bliss
{
  namespace utils
  {

    // using aliases for cleaner syntax
    template<class T> using Invoke = typename T::type;

    template<typename T, T...> struct seq{ using type = seq; };    // TCP: "variable" holding a sequence of values. values may not be used.

    template<class S1, class S2> struct concat;                    // TCP: generic concat template

    template<typename T, T... I1, T... I2>                         // TCP: template specialization of concat that does some real work
    struct concat<seq<T, I1...>, seq<T, I2...> >
      : seq<T, I1..., (sizeof...(I1) + I2)...>{};                  // TCP: note that the values in I2 are offset by size of I1.

    template<class S1, class S2>
    using Concat = Invoke<concat<S1, S2> >;                        // TCP: alias for convenience - function like.

    template<typename T, T N> struct gen_seq;                      // TCP: forward declare
    template<typename T, T N> using GenSeq = Invoke<gen_seq<T, N> >;

    template<typename T, T N>                                      // TCP: actual recursive template instance.
    struct gen_seq : Concat<GenSeq<T, N/2>, GenSeq<T, N - N/2> >{};

    template<typename T>
    struct gen_seq<T, 0> : seq<T>{};                               // TCP: specializations to handle recursion termination.
    template<typename T>
    struct gen_seq<T, 1> : seq<T, 0>{};


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
        static_assert( std::is_integral<T>::value, "Integral type only" );

        /**
         * @typedef type
         * @brief type for the values in sequence.
         */
        using type = T;

        /**
         * @var   size
         * @brief number of values in the parameter pack.
         * @note can only support non-negative integral values.
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
     * @tparam Args  variadic list of  type template arguments
     */
    template<typename... Args>
      using index_sequence_for = make_index_sequence<sizeof...(Args)>;


  } // namespace utils

} // namespace bliss



#endif /* INTEGER_SEQUENCE_HPP_ */
