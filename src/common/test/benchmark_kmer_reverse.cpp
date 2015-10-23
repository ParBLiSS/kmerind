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
#include "utils/bitgroup_ops.hpp"

#include <random>
#include <cstdint>

#include "common/kmer.hpp"
#include "common/test/kmer_reverse_helper.hpp"

#include "common/alphabets.hpp"
#include "common/alphabet_traits.hpp"



// include files to test

//TESTS: Sequential, SWAR/BSWAP, SSSE3, AVX2 versions of kmer reverse.
//TESTS: for each, test kmer size, different word type, and Alphabet (bit group size)
//       different offsets, different bit group sizes, different word types, and different byte array lengths.

template <typename T>
class KmerReverseBenchmark : public ::testing::Test {
  protected:

    T kmer;
    bliss::common::test::KmerReverseHelper<T> helper;

    static const size_t iterations = 1000001;

    virtual void SetUp()
    {
      srand(0);
      for (unsigned int i = 0; i < T::size; ++i) {

        kmer.nextFromChar(rand() % T::KmerAlphabet::SIZE);
      }

    }

};


// indicate this is a typed test
TYPED_TEST_CASE_P(KmerReverseBenchmark);

TYPED_TEST_P(KmerReverseBenchmark, rev_seq)
{
  TypeParam km, rev;
  km = this->kmer;

  for (size_t i = 0; i < this->iterations; ++i) {
    rev ^= this->helper.reverse_serial(km);


    km.nextFromChar(rand() % TypeParam::KmerAlphabet::SIZE);
  }
  EXPECT_FALSE(rev == km);
}

TYPED_TEST_P(KmerReverseBenchmark, rev_seq2)
{
  TypeParam km, rev, tmp;
  km = this->kmer;

  uint8_t* out = reinterpret_cast<uint8_t*>(tmp.getData());
  const uint8_t* in = reinterpret_cast<uint8_t const *>(km.getConstData());

  for (size_t i = 0; i < this->iterations; ++i) {

    memset(out, 0, TypeParam::nBytes);

    bliss::utils::bit_ops::reverse_seq<TypeParam::bitsPerChar>(out, in, TypeParam::nBytes);
    tmp.right_shift_bits(TypeParam::nBytes * 8 - TypeParam::nBits);  // shift by remainder/padding.

    rev ^= tmp;

    km.nextFromChar(rand() % TypeParam::KmerAlphabet::SIZE);
  }
  EXPECT_FALSE(rev == km);
}


TYPED_TEST_P(KmerReverseBenchmark, revcomp_seq)
{
  TypeParam km, rev;
  km = this->kmer;


  for (size_t i = 0; i < this->iterations; ++i) {
    rev ^= this->helper.reverse_complement_serial(km);

    km.nextFromChar(rand() % TypeParam::KmerAlphabet::SIZE);
  }
  EXPECT_FALSE(rev == km);
}


TYPED_TEST_P(KmerReverseBenchmark, rev_bswap)
{
  if (TypeParam::bitsPerChar == 3) return;


  TypeParam km, rev;
  km = this->kmer;


  for (size_t i = 0; i < this->iterations; ++i) {
    rev ^= this->helper.reverse_bswap(km);

    km.nextFromChar(rand() % TypeParam::KmerAlphabet::SIZE);
  }
  EXPECT_FALSE(rev == km);
}

TYPED_TEST_P(KmerReverseBenchmark, revcomp_bswap)
{
  if (TypeParam::bitsPerChar == 3) return;


  TypeParam km, rev;
  km = this->kmer;


  for (size_t i = 0; i < this->iterations; ++i) {
    rev ^= this->helper.reverse_complement_bswap(km);

    km.nextFromChar(rand() % TypeParam::KmerAlphabet::SIZE);
  }
  EXPECT_FALSE(rev == km);
}


TYPED_TEST_P(KmerReverseBenchmark, rev_swar)
{
  if (TypeParam::bitsPerChar == 3) return;


  TypeParam km, rev;
  km = this->kmer;


  for (size_t i = 0; i < this->iterations; ++i) {
    rev ^= this->helper.reverse_swar(km);

    km.nextFromChar(rand() % TypeParam::KmerAlphabet::SIZE);
  }
  EXPECT_FALSE(rev == km);
}

TYPED_TEST_P(KmerReverseBenchmark, revcomp_swar)
{
  if (TypeParam::bitsPerChar == 3) return;


  TypeParam km, rev;
  km = this->kmer;


  for (size_t i = 0; i < this->iterations; ++i) {
    rev ^= this->helper.reverse_complement_swar(km);

    km.nextFromChar(rand() % TypeParam::KmerAlphabet::SIZE);
  }
  EXPECT_FALSE(rev == km);
}

TYPED_TEST_P(KmerReverseBenchmark, rev)
{

  TypeParam km, rev;
  km = this->kmer;


  for (size_t i = 0; i < this->iterations; ++i) {
    rev ^= km.reverse();

    km.nextFromChar(rand() % TypeParam::KmerAlphabet::SIZE);
  }
  EXPECT_FALSE(rev == km);
}

TYPED_TEST_P(KmerReverseBenchmark, revcomp)
{

  TypeParam km, rev;
  km = this->kmer;


  for (size_t i = 0; i < this->iterations; ++i) {
    rev ^= km.reverse_complement();

    km.nextFromChar(rand() % TypeParam::KmerAlphabet::SIZE);
  }
  EXPECT_FALSE(rev == km);
}



TYPED_TEST_P(KmerReverseBenchmark, rev_ssse3)
{
#ifdef __SSSE3__
  if (TypeParam::bitsPerChar == 3) return;


  TypeParam km, rev;
  km = this->kmer;


  for (size_t i = 0; i < this->iterations; ++i) {
    rev ^= this->helper.reverse_simd(km);

    km.nextFromChar(rand() % TypeParam::KmerAlphabet::SIZE);
  }
  EXPECT_FALSE(rev == km);
}

TYPED_TEST_P(KmerReverseBenchmark, revcomp_ssse3)
{
  if (TypeParam::bitsPerChar == 3) return;


  TypeParam km, rev;
  km = this->kmer;


  for (size_t i = 0; i < this->iterations; ++i) {
    rev ^= this->helper.reverse_complement_simd(km);

    km.nextFromChar(rand() % TypeParam::KmerAlphabet::SIZE);
  }
  EXPECT_FALSE(rev == km);
#else
  WARNINGF("SSSE3 is not enabled or not available.");
#endif
}



REGISTER_TYPED_TEST_CASE_P(KmerReverseBenchmark, rev_seq, rev_seq2, revcomp_seq, rev_bswap, revcomp_bswap, rev_swar, revcomp_swar, rev, revcomp, rev_ssse3, revcomp_ssse3);

//////////////////// RUN the tests with different types.

// max of 50 cases
typedef ::testing::Types<
    ::bliss::common::Kmer< 31, bliss::common::DNA,   uint64_t>,  // 1 word, not full
    ::bliss::common::Kmer< 32, bliss::common::DNA,   uint64_t>,  // 1 word, full
    ::bliss::common::Kmer< 33, bliss::common::DNA,   uint64_t>,  // 2 words, not full
    ::bliss::common::Kmer< 80, bliss::common::DNA,   uint64_t>,  // 3 words, not full
    ::bliss::common::Kmer< 96, bliss::common::DNA,   uint64_t>,  // 3 words, full
    ::bliss::common::Kmer< 31, bliss::common::DNA,   uint32_t>,  // 1 word, not full
    ::bliss::common::Kmer< 32, bliss::common::DNA,   uint32_t>,  // 1 word, full
    ::bliss::common::Kmer< 33, bliss::common::DNA,   uint32_t>,  // 2 words, not full
    ::bliss::common::Kmer< 80, bliss::common::DNA,   uint32_t>,  // 3 words, not full
    ::bliss::common::Kmer< 96, bliss::common::DNA,   uint32_t>,  // 3 words, full
    ::bliss::common::Kmer< 31, bliss::common::DNA,   uint16_t>,  // 1 word, not full
    ::bliss::common::Kmer< 32, bliss::common::DNA,   uint16_t>,  // 1 word, full
    ::bliss::common::Kmer< 33, bliss::common::DNA,   uint16_t>,  // 2 words, not full
    ::bliss::common::Kmer< 80, bliss::common::DNA,   uint16_t>,  // 3 words, not full
    ::bliss::common::Kmer< 96, bliss::common::DNA,   uint16_t>,  // 3 words, full
    ::bliss::common::Kmer< 31, bliss::common::DNA,   uint8_t >,  // 1 word, not full
    ::bliss::common::Kmer< 32, bliss::common::DNA,   uint8_t >,  // 1 word, full
    ::bliss::common::Kmer< 33, bliss::common::DNA,   uint8_t >,  // 2 words, not full
    ::bliss::common::Kmer< 80, bliss::common::DNA,   uint8_t >,  // 3 words, not full
    ::bliss::common::Kmer< 96, bliss::common::DNA,   uint8_t >,  // 3 words, full
    ::bliss::common::Kmer< 21, bliss::common::DNA5,  uint64_t>,  // 1 word, not full
    ::bliss::common::Kmer< 22, bliss::common::DNA5,  uint64_t>,  // 2 words, not full
    ::bliss::common::Kmer< 43, bliss::common::DNA5,  uint64_t>,  // 3 words, not full
    ::bliss::common::Kmer< 64, bliss::common::DNA5,  uint64_t>,  // 3 words, full
    ::bliss::common::Kmer< 21, bliss::common::DNA5,  uint8_t >,  // 1 word, not full
    ::bliss::common::Kmer< 22, bliss::common::DNA5,  uint8_t >,  // 2 words, not full
    ::bliss::common::Kmer< 43, bliss::common::DNA5,  uint8_t >,  // 3 words, not full
    ::bliss::common::Kmer< 64, bliss::common::DNA5,  uint8_t >,  // 3 words, full
    ::bliss::common::Kmer< 15, bliss::common::DNA16, uint64_t>,  // 1 word, not full
    ::bliss::common::Kmer< 16, bliss::common::DNA16, uint64_t>,  // 1 word, full
    ::bliss::common::Kmer< 17, bliss::common::DNA16, uint64_t>,  // 2 words, not full
    ::bliss::common::Kmer< 40, bliss::common::DNA16, uint64_t>,  // 3 words, not full
    ::bliss::common::Kmer< 48, bliss::common::DNA16, uint64_t>,  // 3 words, full
    ::bliss::common::Kmer< 15, bliss::common::DNA16, uint32_t>,  // 1 word, not full
    ::bliss::common::Kmer< 16, bliss::common::DNA16, uint32_t>,  // 1 word, full
    ::bliss::common::Kmer< 17, bliss::common::DNA16, uint32_t>,  // 2 words, not full
    ::bliss::common::Kmer< 40, bliss::common::DNA16, uint32_t>,  // 3 words, not full
    ::bliss::common::Kmer< 48, bliss::common::DNA16, uint32_t>,  // 3 words, full
    ::bliss::common::Kmer< 15, bliss::common::DNA16, uint16_t>,  // 1 word, not full
    ::bliss::common::Kmer< 16, bliss::common::DNA16, uint16_t>,  // 1 word, full
    ::bliss::common::Kmer< 17, bliss::common::DNA16, uint16_t>,  // 2 words, not full
    ::bliss::common::Kmer< 40, bliss::common::DNA16, uint16_t>,  // 3 words, not full
    ::bliss::common::Kmer< 48, bliss::common::DNA16, uint16_t>,  // 3 words, full
    ::bliss::common::Kmer< 15, bliss::common::DNA16, uint8_t >,  // 1 word, not full
    ::bliss::common::Kmer< 16, bliss::common::DNA16, uint8_t >,  // 1 word, full
    ::bliss::common::Kmer< 17, bliss::common::DNA16, uint8_t >,  // 2 words, not full
    ::bliss::common::Kmer< 40, bliss::common::DNA16, uint8_t >,  // 3 words, not full
    ::bliss::common::Kmer< 48, bliss::common::DNA16, uint8_t >   // 3 words, full
> KmerReverseBenchmarkTypes;
INSTANTIATE_TYPED_TEST_CASE_P(Bliss, KmerReverseBenchmark, KmerReverseBenchmarkTypes);

