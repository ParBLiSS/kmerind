// include google test
#include <gtest/gtest.h>
//#include <boost/concept_check.hpp>

// include classes to test
#include "common/kmer.hpp"
#include "common/alphabets.hpp"
#include "iterators/transform_iterator.hpp"
#include "common/kmer_iterators.hpp"
#include "utils/kmer_utils.hpp"
#include "utils/logging.h"

template<typename Alphabet, int K>
void compute_kmer_iter(std::string input) {

  using KmerType = bliss::common::Kmer<K, Alphabet>;

  using BaseIterator = std::string::const_iterator;

  using Decoder = bliss::common::ASCII2<Alphabet, typename BaseIterator::value_type>;
  using BaseCharIterator = bliss::iterator::transform_iterator<BaseIterator, Decoder>;

  BaseCharIterator charStart(input.cbegin(), Decoder());
  BaseCharIterator charEnd  (input.cend(),   Decoder());

  using KmerIterator = bliss::common::KmerGenerationIterator<BaseCharIterator, KmerType>;

  KmerIterator start(charStart, true);
  KmerIterator end(charEnd, false);

  int i = 0;
  KmerType kmer;
  for (; start != end; ++start, ++i) {
    kmer = *start;

    std::string gold = input.substr(i, K);

    int res = strncmp(gold.c_str(), bliss::utils::KmerUtils::toASCIIString(kmer).c_str(), K);

    if (res != 0) {
      INFOF("%d iterator input %s\n", i, gold.c_str());
      INFOF("kmer %s %s %s\n", kmer.toString().c_str(), kmer.toAlphabetString().c_str(), bliss::utils::KmerUtils::toASCIIString(kmer).c_str());
    }

    EXPECT_EQ(0, res);
  }
}


/**
 * Test k-mer generation with 2 bits for each character
 */
TEST(KmerIterator, TestKmerIterator)
{
  // test sequence: GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
  std::string input = "GATTTGGGGTTCAAAGCAGT"
                         "ATCGATCAAATAGTAAATCC"
                         "ATTTGTTCAACTCACAGTTT";

  compute_kmer_iter<bliss::common::DNA, 21>(input);
}

/**
 * Test k-mer generation with 2 bits for each character
 */
TEST(KmerIterator, TestKmerIteratorDNA5)
{
  // test sequence: GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
  std::string input = "GATTTGGGGTTCAAAGCAGT"
                         "ATCGATCAAATAGTAAATCC"
                         "ATTTGTTCAACTCACAGTTT";


  compute_kmer_iter<bliss::common::DNA5, 21>(input);

}

/**
 * Test k-mer generation with 2 bits for each character
 */
TEST(KmerIterator, TestKmerIteratorMultiWord)
{
  // test sequence: GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
  std::string input = "GATTTGGGGTTCAAAGCAGT"
                         "ATCGATCAAATAGTAAATCC"
                         "ATTTGTTCAACTCACAGTTT";

  compute_kmer_iter<bliss::common::DNA, 33>(input);

}

/**
 * Test k-mer generation with 2 bits for each character
 */
TEST(KmerIterator, TestKmerIteratorDNA5MultiWord)
{
  // test sequence: GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
  std::string input = "GATTTGGGGTTCAAAGCAGT"
                         "ATCGATCAAATAGTAAATCC"
                         "ATTTGTTCAACTCACAGTTT";
  compute_kmer_iter<bliss::common::DNA5, 33>(input);


}


