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
 * @file    test_kmer_reverse.cpp
 * @ingroup
 * @author  tpan
 * @brief
 * @details
 *
 */



#include "utils/logging.h"

// include google test
#include <gtest/gtest.h>
#include "utils/bitgroup_ops.hpp"

#include <random>
#include <cstdint>

#include <atomic>

#include "common/kmer.hpp"
#include "common/alphabets.hpp"
#include "common/alphabet_traits.hpp"

#include "common/test/kmer_reverse_helper.hpp"



// include files to test

//TESTS: Sequential, SWAR/BSWAP, SSSE3, AVX2 versions of kmer reverse.
//TESTS: for each, test kmer size, different word type, and Alphabet (bit group size)
//       different offsets, different bit group sizes, different word types, and different byte array lengths.

template <typename T>
class KmerReverseTest : public ::testing::Test {
  protected:

    T kmer;
    bliss::common::test::KmerReverseHelper<T> helper;

    static const size_t iterations = 100000;

    virtual void SetUp()
    {
      srand(0);
      for (unsigned int i = 0; i < T::size; ++i) {
//    	  printf(" setup next from char\n");
        kmer.nextFromChar(rand() % T::KmerAlphabet::SIZE);
      }

    }

	template <typename data_type, size_t data_size>
    void print(data_type (&w)[data_size]) {
      std::cout << "data type size " << std::dec << sizeof(data_type) << " len " << data_size << ": ";
      for (int k = data_size -1 ; k >= 0; --k) {
        std::cout << std::hex << std::setfill('0') << std::setw(sizeof(data_type) * 2) << static_cast<size_t>(w[k]) << " ";
      }
      std::cout << std::endl;
    }

};


// indicate this is a typed test
TYPED_TEST_CASE_P(KmerReverseTest);

TYPED_TEST_P(KmerReverseTest, reverse_seq_self)
{
  TypeParam km, gold;
  km = this->kmer;
  gold = this->kmer;

  bool same = true;
  bool local_same;

  for (size_t i = 0; i < this->iterations; ++i) {
    km = this->helper.reverse_serial(km);
    km = this->helper.reverse_complement_serial(km);

    km = this->helper.reverse_serial(km);
    km = this->helper.reverse_complement_serial(km);

    local_same = (km == gold);

    if (!local_same) {
      BL_DEBUGF("ERROR: seq rev-revcomp-rev-revcomp diff at iter %lu:\n\tinput %s\n\toutput %s", i, gold.toAlphabetString().c_str(), km.toAlphabetString().c_str());
      std::cout << "output: ";  this->print(km.getDataRef());
      std::cout << "gold: ";  this->print(gold.getDataRef());
    }
    ASSERT_TRUE(local_same);

    same &= local_same;

    gold.nextFromChar(rand() % TypeParam::KmerAlphabet::SIZE);
    km = gold;
  }
  EXPECT_TRUE(same);
}



TYPED_TEST_P(KmerReverseTest, reverse_seq)
{
  TypeParam km, rev, rev_seq, rev2;
  km = this->kmer;

  bool rev_same = true;
  bool local_rev_same = true;

  uint8_t* out = reinterpret_cast<uint8_t*>(rev.getDataRef());
  const uint8_t* in = reinterpret_cast<uint8_t const *>(km.getData());

  for (size_t i = 0; i < this->iterations; ++i) {
//	printf("rev_seq serial %lu\n", i);
    rev_seq = this->helper.reverse_serial(km);

//	printf("rev_seq BITSEQ reverse %lu\n", i);
    bliss::utils::bit_ops::reverse<TypeParam::bitsPerChar, bliss::utils::bit_ops::BIT_REV_SEQ>(out, in, TypeParam::nBytes);
//    bliss::utils::bit_ops::right_shift<bliss::utils::bit_ops::BITREV_AUTO_AGGRESSIVE<TypeParam::nWords * sizeof(typename TypeParam::KmerWordType)>,
//    		>(rev2.getDataRef(), rev.getDataRef());
//	printf("rev_seq shift %lu\n", i);

    ::std::atomic_thread_fence(::std::memory_order_seq_cst);
    rev.template right_shift_bits<(TypeParam::nBytes * 8 - TypeParam::nBits)>();  // shift by remainder bits..

    local_rev_same = (rev == rev_seq);

    if (!local_rev_same) {
    	printf("nwords %u word size %lu, shift by %lu\n",
    			TypeParam::nWords,
    			sizeof(typename TypeParam::KmerWordType),
    			(TypeParam::nWords * sizeof(typename TypeParam::KmerWordType) * 8 - TypeParam::nBits));
      BL_ERROR("ERROR: rev diff at iter " << i << ":");
//      printf("printing km\n");
      BL_ERROR("\tinput " << km.toAlphabetString().c_str());
//      printf("printing rev2\n");
      BL_ERROR("\toutput " << rev.toAlphabetString().c_str());
//      printf("printing rev\n");
      BL_ERROR("\tgold " << rev_seq.toAlphabetString().c_str());
      std::cout << "output: ";  this->print(rev.getDataRef());
      std::cout << "gold: ";  this->print(rev_seq.getDataRef());
    }
    ASSERT_TRUE(local_rev_same);

    rev_same &= local_rev_same;

//	printf("rev_seq BITSEQ next from char %lu\n", i);
    km.nextFromChar(rand() % TypeParam::KmerAlphabet::SIZE);
  }
  ASSERT_TRUE(rev_same);

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
//		printf("rev_bswap rev-seq %lu\n", i);
    rev_seq = this->helper.reverse_serial(km);
//	printf("rev_bswap revcomp-seq %lu\n", i);
    revcomp_seq = this->helper.reverse_complement_serial(km);

//	printf("rev_bswap rev-bswap %lu\n", i);
    rev = this->helper.reverse_bswap(km);
//	printf("rev_bswap revcomp-bswap %lu\n", i);
    revcomp = this->helper.reverse_complement_bswap(km);

    local_rev_same = (rev == rev_seq);
    local_revcomp_same = (revcomp == revcomp_seq);

    if (!local_rev_same) {
//    	printf("rev_bswap printing %lu\n", i);

      BL_ERRORF("ERROR: rev diff at iter %lu:\n\tinput %s\n\toutput %s\n\tgold %s", i,
    		  km.toAlphabetString().c_str(),
    		  rev.toAlphabetString().c_str(),
    		  rev_seq.toAlphabetString().c_str());
      std::cout << "output: ";  this->print(rev.getDataRef());
      std::cout << "gold: ";  this->print(rev_seq.getDataRef());
//      printf("rev_bswap printing done %lu\n", i);
    }
    if (!local_revcomp_same) {
//    	printf("rev_bswap printing %lu\n", i);
      BL_ERRORF("ERROR: revcomp diff at iter %lu:\n\tinput %s\n\toutput %s\n\tgold %s", i,
    		  km.toAlphabetString().c_str(),
    		  revcomp.toAlphabetString().c_str(),
    		  revcomp_seq.toAlphabetString().c_str());
      std::cout << "output: ";  this->print(revcomp.getDataRef());
      std::cout << "gold: ";  this->print(revcomp_seq.getDataRef());
//  	printf("rev_bswap printing done %lu\n", i);
    }
    ASSERT_TRUE(local_rev_same);
    ASSERT_TRUE(local_revcomp_same);

    rev_same &= local_rev_same;
    revcomp_same &= local_revcomp_same;

//	printf("rev_bswap next char %lu\n", i);
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
    rev_seq = this->helper.reverse_serial(km);
    revcomp_seq = this->helper.reverse_complement_serial(km);

    rev = this->helper.reverse_swar(km);
    revcomp = this->helper.reverse_complement_swar(km);

    local_rev_same = (rev == rev_seq);
    local_revcomp_same = (revcomp == revcomp_seq);

    if (!local_rev_same) {
      BL_ERRORF("ERROR: rev diff at iter %lu:\n\tinput %s\n\toutput %s\n\tgold %s", i, km.toAlphabetString().c_str(), rev.toAlphabetString().c_str(), rev_seq.toAlphabetString().c_str());
      std::cout << "output: ";  this->print(rev.getDataRef());
      std::cout << "gold: ";  this->print(rev_seq.getDataRef());
    }
    if (!local_revcomp_same) {
        BL_ERRORF("ERROR: revcomp diff at iter %lu:\n\tinput %s\n\toutput %s\n\tgold %s", i, km.toAlphabetString().c_str(), revcomp.toAlphabetString().c_str(), revcomp_seq.toAlphabetString().c_str());
        std::cout << "output: ";  this->print(revcomp.getDataRef());
        std::cout << "gold: ";  this->print(revcomp_seq.getDataRef());
    }
    ASSERT_TRUE(local_rev_same);
    ASSERT_TRUE(local_revcomp_same);

    rev_same &= local_rev_same;
    revcomp_same &= local_revcomp_same;

    km.nextFromChar(rand() % TypeParam::KmerAlphabet::SIZE);
  }
  EXPECT_TRUE(rev_same);
  EXPECT_TRUE(revcomp_same);

}


TYPED_TEST_P(KmerReverseTest, reverse_ssse3)
{
#ifdef __SSSE3__
  if (TypeParam::bitsPerChar == 3) return;


  TypeParam km, rev, revcomp, rev_seq, revcomp_seq;
  km = this->kmer;

  bool rev_same = true;
  bool revcomp_same = true;
  bool local_rev_same = true;
  bool local_revcomp_same = true;

  for (size_t i = 0; i < this->iterations; ++i) {
    rev_seq = this->helper.reverse_serial(km);
    revcomp_seq = this->helper.reverse_complement_serial(km);

    rev = this->helper.reverse_simd(km);
    revcomp = this->helper.reverse_complement_simd(km);

    local_rev_same = (rev == rev_seq);
    local_revcomp_same = (revcomp == revcomp_seq);

    if (!local_rev_same) {
      BL_ERRORF("ERROR: rev diff at iter %lu:\n\tinput %s\n\toutput %s\n\tgold %s", i, km.toAlphabetString().c_str(), rev.toAlphabetString().c_str(), rev_seq.toAlphabetString().c_str());
      std::cout << "output: ";  this->print(rev.getDataRef());
      std::cout << "gold: ";  this->print(rev_seq.getDataRef());
    }
    if (!local_revcomp_same) {
        BL_ERRORF("ERROR: revcomp diff at iter %lu:\n\tinput %s\n\toutput %s\n\tgold %s", i, km.toAlphabetString().c_str(), revcomp.toAlphabetString().c_str(), revcomp_seq.toAlphabetString().c_str());
        std::cout << "output: ";  this->print(revcomp.getDataRef());
        std::cout << "gold: ";  this->print(revcomp_seq.getDataRef());
    }
    ASSERT_TRUE(local_rev_same);
    ASSERT_TRUE(local_revcomp_same);

    rev_same &= local_rev_same;
    revcomp_same &= local_revcomp_same;


    km.nextFromChar(rand() % TypeParam::KmerAlphabet::SIZE);
  }
  EXPECT_TRUE(rev_same);
  EXPECT_TRUE(revcomp_same);
#else
  BL_WARNINGF("SSSE3 is not enabled or not available.");
#endif

}

TYPED_TEST_P(KmerReverseTest, reverse)
{
  TypeParam km, rev, revcomp, rev_seq, revcomp_seq;
  km = this->kmer;

  bool rev_same = true;
  bool revcomp_same = true;
  bool local_rev_same = true;
  bool local_revcomp_same = true;

  for (size_t i = 0; i < this->iterations; ++i) {
    rev_seq = this->helper.reverse_serial(km);
    revcomp_seq = this->helper.reverse_complement_serial(km);

    rev = km.reverse();
    revcomp = km.reverse_complement();

    local_rev_same = (rev == rev_seq);
    local_revcomp_same = (revcomp == revcomp_seq);

    if (!local_rev_same) {
      BL_ERRORF("ERROR: rev diff at iter %lu:\n\tinput %s\n\toutput %s\n\tgold %s", i, km.toAlphabetString().c_str(), rev.toAlphabetString().c_str(), rev_seq.toAlphabetString().c_str());
      std::cout << "output: ";  this->print(rev.getDataRef());
      std::cout << "gold: ";  this->print(rev_seq.getDataRef());
    }
    if (!local_revcomp_same) {
        BL_ERRORF("ERROR: revcomp diff at iter %lu:\n\tinput %s\n\toutput %s\n\tgold %s", i, km.toAlphabetString().c_str(), revcomp.toAlphabetString().c_str(), revcomp_seq.toAlphabetString().c_str());
        std::cout << "output: ";  this->print(revcomp.getDataRef());
        std::cout << "gold: ";  this->print(revcomp_seq.getDataRef());
    }
    rev_same &= local_rev_same;
    revcomp_same &= local_revcomp_same;
    ASSERT_TRUE(local_rev_same);
    ASSERT_TRUE(local_revcomp_same);

    km.nextFromChar(rand() % TypeParam::KmerAlphabet::SIZE);
  }
  EXPECT_TRUE(rev_same);
  EXPECT_TRUE(revcomp_same);

}

REGISTER_TYPED_TEST_CASE_P(KmerReverseTest, reverse_seq_self, reverse_seq, reverse_bswap, reverse_swar, reverse_ssse3, reverse);


