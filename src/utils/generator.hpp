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
 * @file    generator.hpp
 * @ingroup utils
 * @author  Patrick Flick
 * @brief   Implements input string generators used for benchmarking and testing.
 *
 */
#ifndef BLISS_UTILS_GENERATOR_H
#define BLISS_UTILS_GENERATOR_H

// C++ STL includes[
#include <string>
#include <random>
#include <algorithm>
#include <functional>

namespace bliss
{
namespace utils
{

/**
 * @brief   Generates a random string of specified length using the given
 *          generator function.
 * @param   length    The length of the generated random string.
 * @param   generator The generator function of type `char (void)` which is
 *                    used to generate the random characters of the string.
 * @returns           A random string of length `length`.
 */
std::string random_string(const std::size_t length, const std::function<char(void)> generator)
{
  // create string of size `length`
  std::string result;
  result.resize(length);

  // fill with random characters
  std::generate(result.begin(), result.end(), generator);

  return result;
}

/**
 * @brief   Generates a random string of specified length containing nothing
 *          but the passed characters.
 *
 * @param length  The length of the generated random string.
 * @param chars   The characters to be used in the generated string.
 * @returns       A random string of length `length`.
 */
std::string random_string(const std::size_t length, const std::vector<char>& chars)
{
  // get random engine
  std::default_random_engine generator;
  std::uniform_int_distribution<std::size_t> distribution(0, chars.size() - 1);

  // call the random string generator with our distribution
  return random_string(length, [&](){return chars[distribution(generator)];});
}

/**
 * @brief   Generates a random string of dna characters.
 *
 * @param length  The length of the generated random DNA string.
 * @returns       The random DNA string of length `length`.
 */
std::string random_dna(const std::size_t length)
{
  std::vector<char> dna_chars = {'A','C','G','T'};
  return random_string(length, dna_chars);
}


}

} // namespace bliss

#endif // BLISS_UTILS_GENERATOR_H
