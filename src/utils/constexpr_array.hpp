/**
 * @file		constexpr_array.hpp
 * @ingroup utils
 * @author	Tony Pan <tpan7@gatech.edu>
 * @brief   utility functions to create an index array and transform it into a const value array, at compile time
 *
 * @details create std array of size N and transform the values via a constexpr function to create a const value array at compile time.
 *          adapted from http://stackoverflow.com/questions/19019252/c11-create-0-to-n-constexpr-array-in-c
 *
 *          essentially, recursive template parameter deduction and substitution via variadic templates
 *          to progressively build the final template.
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 *
 */
#ifndef CONSTEXPR_ARRAY_HPP_
#define CONSTEXPR_ARRAY_HPP_

#include <array>
#include "utils/integer_sequence.hpp"

namespace bliss {
  namespace utils {

    namespace detail {

      /**
       * @brief  Expands the index value sequence parameter pack (index_sequence) and apply function to transform them.
       * @tparam Function   function type for transforming from index value to output value
       * @tparam Indices    index values.
       * @param f           function to transform from index value to output value
       * @param             index value sequence parameter pack - for deducing the Indices template parameters only.
       * @return
       */
      template<class Functor, std::size_t ... Indices>
      constexpr auto make_array_helper(
          Functor f,
          bliss::utils::index_sequence<Indices...>)
          -> std::array<typename std::result_of<Functor(std::size_t)>::type, sizeof...(Indices)>
      {
          // double curl braces aggregate initializes the stdarray.
          // ... does a variadic expansion.
        return {{ f(Indices)... }};
      }


      /**
       * @brief  Expands the index value sequence parameter pack (index_sequence) and apply function to transform them.
       * @tparam Function   function type for transforming from index value to output value
       * @tparam Indices    index values.
       * @param f           function to transform from index value to output value
       * @param             index value sequence parameter pack - for deducing the Indices template parameters only.
       * @return
       */
      template<typename T, std::size_t ... Indices>
      constexpr auto make_array_helper(bliss::utils::index_sequence<Indices...>)
          -> std::array<T, sizeof...(Indices)>
      {
          // double curl braces aggregate initializes the stdarray.
          // ... does a variadic expansion.
        return {{ static_cast<T>(Indices)... }};
      }

    }  // namespace detail


    /**
     * @brief             construct and initialize a const_expr array.
     * @tparam N          number of values in array
     * @tparam Functor    constexpr functor type to transform an index (integer) value to a output value.
     * @param f           constexpr functor to transform an index (integer) value to a output value.
     * @return            std::array with the computed constexpr values
     */
    template<size_t N, class Functor>
    constexpr auto make_array(Functor f)
    -> std::array<typename std::result_of<Functor(std::size_t)>::type, N>
    {
        // expands N into a index sequence [0, .., N), then transform each via f.
      return detail::make_array_helper(f, bliss::utils::make_index_sequence<N>{});
    }

    /**
     * @brief             construct and initialize a const_expr array.
     * @tparam N          number of values in array
     * @tparam Functor    constexpr functor type to transform an index (integer) value to a output value.
     * @param f           constexpr functor to transform an index (integer) value to a output value.
     * @return            std::array with the computed constexpr values
     */
    template<typename T, size_t N>
    constexpr auto make_array()
    -> std::array<T, N>
    {
        // expands N into a index sequence [0, .., N), then transform each via f.
      return detail::make_array_helper<T>(bliss::utils::make_index_sequence<N>{});
    }

  } // namespace utils
} // namespace bliss



#endif /* CONSTEXPR_ARRAY_HPP_ */
