/**
 * @file    kmer_reverse.cpp
 * @ingroup
 * @author  tpan
 * @brief
 * @details
 *
 * Copyright (c) 2015 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */

#include <string>
#include <random>

#include "common/kmer.hpp"
#include "common/alphabets.hpp"
#include "common/alphabet_traits.hpp"
#include "utils/timer.hpp"
#include "iterators/transform_iterator.hpp"

template<unsigned int K, typename Alphabet, typename WORD_TYPE>
void make_kmer(std::string input) {

  using KmerType = bliss::common::Kmer<K, Alphabet, WORD_TYPE>;

  using Decoder = bliss::common::ASCII2<Alphabet, std::string::value_type>;
  using BaseCharIterator = bliss::iterator::transform_iterator<std::string::const_iterator, Decoder>;
  auto temp = BaseCharIterator(input.cbegin(), Decoder());

  KmerType kmer;
  kmer.fillFromChars(temp, false);
  unsigned int i = K;
  for (; i < input.length(); ++i) {
    std::string gold = input.substr(i-K, K);

    int res = strncmp(gold.c_str(), bliss::utils::KmerUtils::toASCIIString(kmer).c_str(), K);

    if (res != 0) {
      INFOF("%d iterator input %s\n", i, gold.c_str());
      INFOF("kmer %s %s %s\n", kmer.toString().c_str(), kmer.toAlphabetString().c_str(), bliss::utils::KmerUtils::toASCIIString(kmer).c_str());
    }

    kmer.nextFromChar(*temp);  ++temp;

  }
  std::string gold = input.substr(input.length()-K, K);

  int res = strncmp(gold.c_str(), bliss::utils::KmerUtils::toASCIIString(kmer).c_str(), K);

  if (res != 0) {
    INFOF("%d iterator input %s\n", i, gold.c_str());
    INFOF("kmer %s %s %s\n", kmer.toString().c_str(), kmer.toAlphabetString().c_str(), bliss::utils::KmerUtils::toASCIIString(kmer).c_str());
  }

}


template <unsigned int KMER_SIZE, typename ALPHABET, typename WORD_TYPE>
void test_reverse() {

  using KMER = bliss::common::Kmer<KMER_SIZE, ALPHABET, WORD_TYPE>;
  KMER kmer, rev, revcomp, rev2, revcomp2;
  int iterations = 100000;

  printf("testing kmer with k=%d, nbits = %d, word bits =%lu\n", KMER_SIZE, KMER::nBits, sizeof(WORD_TYPE) * 8);

  std::random_device r;  std::mt19937 gen(r());
  std::uniform_int_distribution<unsigned char> dist;
  for (int i = 0; i < KMER_SIZE; ++i) {
    kmer.nextFromChar(dist(gen) % ALPHABET::SIZE);
  }
//  printf("Input kmer %s\n", kmer.toAlphabetString().c_str());

  TIMER_INIT(test);

  TIMER_START(test);
  for (int i = 0; i  < iterations; ++i) {
    rev = kmer.reversed_kmer();
  }
  TIMER_END(test, "old reverse", iterations);

  TIMER_START(test);
  for (int i = 0; i  < iterations; ++i) {
    revcomp = kmer.reverse_complement();
  }
  TIMER_END(test, "old revcomp", iterations);

  TIMER_START(test);
  for (int i = 0; i  < iterations; ++i) {
    rev2 = kmer.reversed_kmer2();
  }
  TIMER_END(test, "new reverse", iterations);

  TIMER_START(test);
  for (int i = 0; i  < iterations; ++i) {
    revcomp2 = kmer.reverse_complement2();
  }
  TIMER_END(test, "new revcomp", iterations);

  if (rev != rev2) {
    printf("ERROR: old and new method of computing reverse complement produced different results:\n  rev %s, rev2 %s\n", rev.toAlphabetString().c_str(), rev2.toAlphabetString().c_str());
  }

  if (revcomp != revcomp2) {
    printf("ERROR: old and new method of computing reverse complement produced different results:\n  rev %s, rev2 %s\n", revcomp.toAlphabetString().c_str(), revcomp2.toAlphabetString().c_str());
  }


  TIMER_REPORT(test, 0);

}


int main(int argc, char** argv) {

  test_reverse< 31, bliss::common::DNA, uint64_t>();  // 1 word, not full
  test_reverse< 32, bliss::common::DNA, uint64_t>();  // 1 word, full
  test_reverse< 33, bliss::common::DNA, uint64_t>();  // 2 words, not full
  test_reverse< 64, bliss::common::DNA, uint64_t>();  // 2 words, full
  test_reverse< 80, bliss::common::DNA, uint64_t>();  // 3 words, not full
  test_reverse< 96, bliss::common::DNA, uint64_t>();  // 3 words, full

  test_reverse< 3,  bliss::common::DNA, uint8_t>();  // 1 word, not full
  test_reverse< 4,  bliss::common::DNA, uint8_t>();  // 1 word, full
  test_reverse< 5,  bliss::common::DNA, uint8_t>();  // 2 words, not full
  test_reverse< 8,  bliss::common::DNA, uint8_t>();  // 2 words, full
  test_reverse< 10, bliss::common::DNA, uint8_t>();  // 3 words, not full
  test_reverse< 12, bliss::common::DNA, uint8_t>();  // 3 words, full

  test_reverse< 21, bliss::common::DNA5, uint64_t>();  // 1 word, not full
  test_reverse< 22, bliss::common::DNA5, uint64_t>();  // 2 word, not full
  test_reverse< 42, bliss::common::DNA5, uint64_t>();  // 2 words, not full
  test_reverse< 43, bliss::common::DNA5, uint64_t>();  // 3 words, not full
  test_reverse< 64, bliss::common::DNA5, uint64_t>();  // 3 words, full

  test_reverse< 2,  bliss::common::DNA5, uint8_t>();  // 1 word, not full
  test_reverse< 3,  bliss::common::DNA5, uint8_t>();  // 2 word, not full
  test_reverse< 5,  bliss::common::DNA5, uint8_t>();  // 2 words, not full
  test_reverse< 6,  bliss::common::DNA5, uint8_t>();  // 3 words, not full
  test_reverse< 8,  bliss::common::DNA5, uint8_t>();  // 3 words, full

  test_reverse< 15, bliss::common::DNA16, uint64_t>();  // 1 word, not full
  test_reverse< 16, bliss::common::DNA16, uint64_t>();  // 1 word, full
  test_reverse< 17, bliss::common::DNA16, uint64_t>();  // 2 words, not full
  test_reverse< 32, bliss::common::DNA16, uint64_t>();  // 2 words, full
  test_reverse< 40, bliss::common::DNA16, uint64_t>();  // 3 words, not full
  test_reverse< 48, bliss::common::DNA16, uint64_t>();  // 3 words, full

  test_reverse< 1,  bliss::common::DNA16, uint8_t >();  // 1 word, not full
  test_reverse< 2,  bliss::common::DNA16, uint8_t >();  // 1 word, full
  test_reverse< 3,  bliss::common::DNA16, uint8_t >();  // 2 words, not full
  test_reverse< 4,  bliss::common::DNA16, uint8_t >();  // 2 words, full
  test_reverse< 5,  bliss::common::DNA16, uint8_t >();  // 3 words, not full
  test_reverse< 6,  bliss::common::DNA16, uint8_t >();  // 3 words, full

}
