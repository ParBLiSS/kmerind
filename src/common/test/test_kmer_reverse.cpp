/**
 * @file    test_kmer_reverse.cpp
 * @ingroup
 * @author  tpan
 * @brief
 * @details
 *
 * Copyright (c) 2015 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */



#include "utils/logging.h"

// include google test
#include <gtest/gtest.h>

#include <random>
#include <cstdint>

#include "common/kmer.hpp"
#include "common/alphabets.hpp"
#include "common/alphabet_traits.hpp"

// include files to test

//TESTS: Sequential, SWAR/BSWAP, SSSE3, AVX2 versions of kmer reverse.
//TESTS: for each, test kmer size, different word type, and Alphabet (bit group size)
//       different offsets, different bit group sizes, different word types, and different byte array lengths.

template <typename T>
class KmerReverseTest : public ::testing::Test {
  protected:

    T kmer;
    static const size_t iterations = 10001;

    virtual void SetUp()
    {
      srand(0);
      for (unsigned int i = 0; i < T::size; ++i) {

        kmer.nextFromChar(rand() % T::KmerAlphabet::SIZE);
      }

    }

};


// indicate this is a typed test
TYPED_TEST_CASE_P(KmerReverseTest);

TYPED_TEST_P(KmerReverseTest, reverse_seq)
{
  TypeParam km, gold;
  km = this->kmer;
  gold = this->kmer;

  bool same = true;
  bool local_same;

  for (size_t i = 0; i < this->iterations; ++i) {
    km = km.reverse_serial();
    km = km.reverse_complement_serial();

    km = km.reverse_serial();
    km = km.reverse_complement_serial();

    local_same = (km == gold);

    if (!local_same) {
      DEBUGF("ERROR: seq rev-revcomp-rev-revcomp diff at iter %lu:\n\tinput %s\n\toutput %s", i, gold.toAlphabetString().c_str(), km.toAlphabetString().c_str());
    }
    same &= local_same;

    gold.nextFromChar(rand() % TypeParam::KmerAlphabet::SIZE);
    km = gold;
  }
  EXPECT_TRUE(same);
}

TYPED_TEST_P(KmerReverseTest, reverse_bswap)
{
  if (TypeParam::bitsPerChar == 3) return;

  TypeParam km, rev, revcomp, rev_seq, revcomp_seq;
  km = this->kmer;

  bool rev_same = true;
  bool revcomp_same = true;
  bool local_rev_same = true;
  bool local_revcomp_same = true;

  for (size_t i = 0; i < this->iterations; ++i) {
    rev_seq = km.reverse_serial();
    revcomp_seq = km.reverse_complement_serial();

    rev = km.reverse_bswap();
    revcomp = km.reverse_complement_bswap();

    local_rev_same = (rev == rev_seq);
    local_revcomp_same = (revcomp == revcomp_seq);

    if (!local_rev_same) {
      DEBUGF("ERROR: rev diff at iter %lu:\n\tinput %s\n\toutput %s\n\tgold %s", i, km.toAlphabetString().c_str(), rev.toAlphabetString().c_str(), rev_seq.toAlphabetString().c_str());
    }
    if (!local_revcomp_same) {
      DEBUGF("ERROR: revcomp diff at iter %lu:\n\tinput %s\n\toutput %s\n\tgold %s", i, km.toAlphabetString().c_str(), revcomp.toAlphabetString().c_str(), revcomp_seq.toAlphabetString().c_str());
    }

    rev_same &= local_rev_same;
    revcomp_same &= local_revcomp_same;


    km.nextFromChar(rand() % TypeParam::KmerAlphabet::SIZE);
  }
  EXPECT_TRUE(rev_same);
  EXPECT_TRUE(revcomp_same);

}

TYPED_TEST_P(KmerReverseTest, reverse_swar)
{
  if (TypeParam::bitsPerChar == 3) return;

  TypeParam km, rev, revcomp, rev_seq, revcomp_seq;
  km = this->kmer;

  bool rev_same = true;
  bool revcomp_same = true;
  bool local_rev_same = true;
  bool local_revcomp_same = true;

  for (size_t i = 0; i < this->iterations; ++i) {
    rev_seq = km.reverse_serial();
    revcomp_seq = km.reverse_complement_serial();

    rev = km.reverse_swar();
    revcomp = km.reverse_complement_swar();

    local_rev_same = (rev == rev_seq);
    local_revcomp_same = (revcomp == revcomp_seq);

    if (!local_rev_same) {
      DEBUGF("ERROR: rev diff at iter %lu:\n\tinput %s\n\toutput %s\n\tgold %s", i, km.toAlphabetString().c_str(), rev.toAlphabetString().c_str(), rev_seq.toAlphabetString().c_str());
    }
    if (!local_revcomp_same) {
      DEBUGF("ERROR: revcomp diff at iter %lu:\n\tinput %s\n\toutput %s\n\tgold %s", i, km.toAlphabetString().c_str(), revcomp.toAlphabetString().c_str(), revcomp_seq.toAlphabetString().c_str());
    }

    rev_same &= local_rev_same;
    revcomp_same &= local_revcomp_same;

    km.nextFromChar(rand() % TypeParam::KmerAlphabet::SIZE);
  }
  EXPECT_TRUE(rev_same);
  EXPECT_TRUE(revcomp_same);

}

TYPED_TEST_P(KmerReverseTest, reverse_new)
{
  TypeParam km, rev, revcomp, rev_seq, revcomp_seq;
  km = this->kmer;

  bool rev_same = true;
  bool revcomp_same = true;
  bool local_rev_same = true;
  bool local_revcomp_same = true;

  for (size_t i = 0; i < this->iterations; ++i) {
    rev_seq = km.reverse_serial();
    revcomp_seq = km.reverse_complement_serial();

    rev = km.reverse_new();
    revcomp = km.reverse_complement_new();

    local_rev_same = (rev == rev_seq);
    local_revcomp_same = (revcomp == revcomp_seq);

    if (!local_rev_same) {
      DEBUGF("ERROR: rev diff at iter %lu:\n\tinput %s\n\toutput %s\n\tgold %s", i, km.toAlphabetString().c_str(), rev.toAlphabetString().c_str(), rev_seq.toAlphabetString().c_str());
    }
    if (!local_revcomp_same) {
      DEBUGF("ERROR: revcomp diff at iter %lu:\n\tinput %s\n\toutput %s\n\tgold %s", i, km.toAlphabetString().c_str(), revcomp.toAlphabetString().c_str(), revcomp_seq.toAlphabetString().c_str());
    }

    rev_same &= local_rev_same;
    revcomp_same &= local_revcomp_same;

    km.nextFromChar(rand() % TypeParam::KmerAlphabet::SIZE);
  }
  EXPECT_TRUE(rev_same);
  EXPECT_TRUE(revcomp_same);

}

#ifdef __SSSE3__
TYPED_TEST_P(KmerReverseTest, reverse_ssse3)
{
  if (TypeParam::bitsPerChar == 3) return;


  TypeParam km, rev, revcomp, rev_seq, revcomp_seq;
  km = this->kmer;

  bool rev_same = true;
  bool revcomp_same = true;
  bool local_rev_same = true;
  bool local_revcomp_same = true;

  for (size_t i = 0; i < this->iterations; ++i) {
    rev_seq = km.reverse_serial();
    revcomp_seq = km.reverse_complement_serial();

    rev = km.reverse_simd();
    revcomp = km.reverse_complement_simd();

    local_rev_same = (rev == rev_seq);
    local_revcomp_same = (revcomp == revcomp_seq);

    if (!local_rev_same) {
      DEBUGF("ERROR: rev diff at iter %lu:\n\tinput %s\n\toutput %s\n\tgold %s", i, km.toAlphabetString().c_str(), rev.toAlphabetString().c_str(), rev_seq.toAlphabetString().c_str());
    }
    if (!local_revcomp_same) {
      DEBUGF("ERROR: revcomp diff at iter %lu:\n\tinput %s\n\toutput %s\n\tgold %s", i, km.toAlphabetString().c_str(), revcomp.toAlphabetString().c_str(), revcomp_seq.toAlphabetString().c_str());
    }

    rev_same &= local_rev_same;
    revcomp_same &= local_revcomp_same;


    km.nextFromChar(rand() % TypeParam::KmerAlphabet::SIZE);
  }
  EXPECT_TRUE(rev_same);
  EXPECT_TRUE(revcomp_same);

}
#endif



#ifdef __SSSE3__
// now register the test cases
REGISTER_TYPED_TEST_CASE_P(KmerReverseTest, reverse_seq, reverse_bswap, reverse_swar, reverse_new, reverse_ssse3);
#else
// now register the test cases
REGISTER_TYPED_TEST_CASE_P(KmerReverseTest, reverse_seq, reverse_bswap, reverse_swar, reverse_new);
#endif

//////////////////// RUN the tests with different types.

// max of 50 cases
typedef ::testing::Types<
    ::bliss::common::Kmer< 31, bliss::common::DNA,   uint64_t>,  // 1 word, not full
    ::bliss::common::Kmer< 32, bliss::common::DNA,   uint64_t>,  // 1 word, full
    ::bliss::common::Kmer< 33, bliss::common::DNA,   uint64_t>,  // 2 words, not full
    ::bliss::common::Kmer< 64, bliss::common::DNA,   uint64_t>,  // 2 words, full
    ::bliss::common::Kmer< 80, bliss::common::DNA,   uint64_t>,  // 3 words, not full
    ::bliss::common::Kmer< 96, bliss::common::DNA,   uint64_t>,  // 3 words, full
    ::bliss::common::Kmer< 15, bliss::common::DNA,   uint32_t>,  // 1 word, not full
    ::bliss::common::Kmer< 16, bliss::common::DNA,   uint32_t>,  // 1 word, full
    ::bliss::common::Kmer< 17, bliss::common::DNA,   uint32_t>,  // 2 words, not full
    ::bliss::common::Kmer< 32, bliss::common::DNA,   uint32_t>,  // 2 words, full
    ::bliss::common::Kmer< 40, bliss::common::DNA,   uint32_t>,  // 3 words, not full
    ::bliss::common::Kmer< 48, bliss::common::DNA,   uint32_t>,  // 3 words, full
    ::bliss::common::Kmer<  7, bliss::common::DNA,   uint16_t>,  // 1 word, not full
    ::bliss::common::Kmer<  8, bliss::common::DNA,   uint16_t>,  // 1 word, full
    ::bliss::common::Kmer<  9, bliss::common::DNA,   uint16_t>,  // 2 words, not full
    ::bliss::common::Kmer< 16, bliss::common::DNA,   uint16_t>,  // 2 words, full
    ::bliss::common::Kmer<  3, bliss::common::DNA,    uint8_t>,  // 1 word, not full
    ::bliss::common::Kmer<  4, bliss::common::DNA,    uint8_t>,  // 1 word, full
    ::bliss::common::Kmer<  5, bliss::common::DNA,    uint8_t>,  // 2 words, not full
    ::bliss::common::Kmer<  8, bliss::common::DNA,    uint8_t>,  // 2 words, full
    ::bliss::common::Kmer< 21, bliss::common::DNA5,  uint64_t>,  // 1 word, not full
    ::bliss::common::Kmer< 22, bliss::common::DNA5,  uint64_t>,  // 2 word, not full
    ::bliss::common::Kmer< 42, bliss::common::DNA5,  uint64_t>,  // 2 words, not full
    ::bliss::common::Kmer< 43, bliss::common::DNA5,  uint64_t>,  // 3 words, not full
    ::bliss::common::Kmer< 64, bliss::common::DNA5,  uint64_t>,  // 3 words, full
    ::bliss::common::Kmer<  2, bliss::common::DNA5,   uint8_t>,  // 1 word, not full
    ::bliss::common::Kmer<  3, bliss::common::DNA5,   uint8_t>,  // 2 word, not full
    ::bliss::common::Kmer<  5, bliss::common::DNA5,   uint8_t>,  // 2 words, not full
    ::bliss::common::Kmer<  6, bliss::common::DNA5,   uint8_t>,  // 3 words, not full
    ::bliss::common::Kmer<  8, bliss::common::DNA5,   uint8_t>,  // 3 words, full
    ::bliss::common::Kmer< 15, bliss::common::DNA16, uint64_t>,  // 1 word, not full
    ::bliss::common::Kmer< 16, bliss::common::DNA16, uint64_t>,  // 1 word, full
    ::bliss::common::Kmer< 17, bliss::common::DNA16, uint64_t>,  // 2 words, not full
    ::bliss::common::Kmer< 32, bliss::common::DNA16, uint64_t>,  // 2 words, full
    ::bliss::common::Kmer< 40, bliss::common::DNA16, uint64_t>,  // 3 words, not full
    ::bliss::common::Kmer< 48, bliss::common::DNA16, uint64_t>,  // 3 words, full
    ::bliss::common::Kmer<  7, bliss::common::DNA16, uint32_t>,  // 1 word, not full
    ::bliss::common::Kmer<  8, bliss::common::DNA16, uint32_t>,  // 1 word, full
    ::bliss::common::Kmer<  9, bliss::common::DNA16, uint32_t>,  // 2 words, not full
    ::bliss::common::Kmer< 16, bliss::common::DNA16, uint32_t>,  // 2 words, full
    ::bliss::common::Kmer< 20, bliss::common::DNA16, uint32_t>,  // 3 words, not full
    ::bliss::common::Kmer< 24, bliss::common::DNA16, uint32_t>,  // 3 words, full
    ::bliss::common::Kmer<  3, bliss::common::DNA16, uint16_t>,  // 1 word, not full
    ::bliss::common::Kmer<  4, bliss::common::DNA16, uint16_t>,  // 1 word, full
    ::bliss::common::Kmer<  5, bliss::common::DNA16, uint16_t>,  // 2 words, not full
    ::bliss::common::Kmer<  8, bliss::common::DNA16, uint16_t>,  // 2 words, full
    ::bliss::common::Kmer<  1, bliss::common::DNA16,  uint8_t>,  // 1 word, not full
    ::bliss::common::Kmer<  2, bliss::common::DNA16,  uint8_t>,  // 1 word, full
    ::bliss::common::Kmer<  3, bliss::common::DNA16,  uint8_t>,  // 2 words, not full
    ::bliss::common::Kmer<  4, bliss::common::DNA16,  uint8_t>   // 2 words, full
> KmerReverseTestTypes;
INSTANTIATE_TYPED_TEST_CASE_P(Bliss, KmerReverseTest, KmerReverseTestTypes);

