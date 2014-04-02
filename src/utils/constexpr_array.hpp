/**
 * @file		constexpr_array.hpp
 * @ingroup
 * @author	tpan
 * @brief   create std array of size N and initialized via a constexpr function based on the index value.
 * @details adapted from http://stackoverflow.com/questions/19019252/c11-create-0-to-n-constexpr-array-in-c
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

template<class Function, std::size_t ... Indices>
constexpr auto make_array_helper(
    Function f,
    bliss::utils::index_sequence<Indices...>)
    -> std::array<typename std::result_of<Function(std::size_t)>::type, sizeof...(Indices)>
{
  return {{ f(Indices)...}};
}

template<int N, class Functor>
constexpr auto make_array(Functor f)
-> std::array<typename std::result_of<Functor(std::size_t)>::type, N>
{
  return make_array_helper(f, bliss::utils::make_index_sequence<N>{});
}

#endif /* CONSTEXPR_ARRAY_HPP_ */
