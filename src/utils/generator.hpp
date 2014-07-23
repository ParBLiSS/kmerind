
/**
 * @file    generator.hpp
 * @ingroup utils
 * @author  Patrick Flick
 * @brief   Implements input generators used for benchmarking and testing.
 *
 * Copyright (c) TODO
 *
 * TODO add Licence
 */
#ifndef BLISS_UTILS_GENERATOR_H
#define BLISS_UTILS_GENERATOR_H

// C++ STL includes[
#include <string>
#include <random>
#include <algorithm>

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



//
///**
// * generate the lengths
// */
//template<int K>
//void generateLengths(std::vector<ReadLengthType>& lengths)
//{
//
//  std::default_random_engine generator;
//  std::normal_distribution<double> distribution(L_MEAN, L_STDEV);
//
//  ReadLengthType l;
//
//  for (ReadCountType i = 0; i < N_R; ++i)
//  {
//    // generate the random lengths of the reads
//    l = static_cast<ReadLengthType>(round(distribution(generator)));
//
//    if (l < 80 || l > 120)
//    {
//      --i;
//      continue;
//    }
//
//    lengths[i] = l;
//  }
//
//}
//
///**
// * allocate input vector.
// * supports packed string with padding bits.
// */
//void generateStrings(std::vector<ReadLengthType> const & lengths,
//                     std::vector<SequenceT>& reads)
//{
//  std::default_random_engine generator;
//  std::uniform_int_distribution<WordType> distribution;
//
//  // use a random number generator to generate the "packed string"
//  for (ReadCountType i = 0; i < N_R; ++i)
//  {
//    // allocate a sequence of packed string blocks.
//    SequenceT read((lengths[i] + N_PACKED_CHARS - 1) / N_PACKED_CHARS, 0);
//
//    int idx = 0;
//    for (ReadLengthType j = 0; j < lengths[i]; j += N_PACKED_CHARS)
//    {
//      // packed bits, with padding.  low order bits only.
//      read[idx] = distribution(generator) & PADDING_MASK;
//      ++idx;
//    }
//    // the last one may not be completely filled.
//    if (lengths[i] % N_PACKED_CHARS > 0)
//    {
//      WordType keep_mask = std::numeric_limits<WordType>::max()
//          >> (WORD_BITS - (lengths[i] % N_PACKED_CHARS) * N_BITS);
//      read[idx - 1] = read[idx - 1] & keep_mask;
//    }
//    read.shrink_to_fit();
//    reads[i] = read;
//  }
//
//}
//


}

} // namespace bliss

#endif // BLISS_UTILS_GENERATOR_H
