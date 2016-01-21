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
 * @file    test_kmer_hash.cpp
 * @ingroup
 * @author  tpan
 * @brief
 * @details
 *
 */



#include "utils/logging.h"

// include google test
#include <gtest/gtest.h>

#include <random>
#include <cstdint>

#include "common/kmer.hpp"
#include "common/alphabets.hpp"
#include "common/alphabet_traits.hpp"
#include "common/kmer_transform.hpp"



// include files to test

//TESTS: transform of kmer.  these are fairly simple, since the component functions are already tested - i.e. revcomplement, comparisons

template <typename T>
class KmerTransformTest : public ::testing::Test {
  protected:

    T kmer;

    static const size_t iterations = 100000;

    virtual void SetUp()
    {
      srand(0);
      for (unsigned int i = 0; i < T::size; ++i) {

        kmer.nextFromChar(rand() % T::KmerAlphabet::SIZE);
      }
    }


};


// indicate this is a typed test
TYPED_TEST_CASE_P(KmerTransformTest);

TYPED_TEST_P(KmerTransformTest, identity)
{

  auto km = this->kmer;
  auto gold = this->kmer;

  bliss::kmer::transform::identity<TypeParam> op;

  bool same = true;
  bool local_same;

  for (size_t i = 0; i < this->iterations; ++i) {
    km = op(km);

    local_same = (km == gold);

    if (!local_same) {
      BL_DEBUGF("ERROR: seq transform diff at iter %lu:\n\tinput %s\n\toutput %s", i, gold.toAlphabetString().c_str(), km.toAlphabetString().c_str());
    }
    same &= local_same;

    gold.nextFromChar(rand() % TypeParam::KmerAlphabet::SIZE);
    km = gold;
  }
  EXPECT_TRUE(same);
}


TYPED_TEST_P(KmerTransformTest, trans_xor)
{

  auto km = this->kmer;
  auto gold = this->kmer;
  auto rev = km.reverse_complement();

  bliss::kmer::transform::xor_rev_comp<TypeParam> op;

  bool same = true;
  bool local_same;


  for (size_t i = 0; i < this->iterations; ++i) {
    // check symmetry.
	local_same = (op(std::make_pair(km, rev)) == op(std::make_pair(rev, km)));
    if (!local_same) {
      BL_DEBUGF("ERROR: seq xor is not symmetric at iter %lu:\n\tkmer %s\n\trevcomp %s", i, km.toAlphabetString().c_str(), rev.toAlphabetString().c_str());
    }
    same &= local_same;

	// check reversable.
	local_same = (gold == op(std::make_pair(op(std::make_pair(km, rev)), rev)));
    if (!local_same) {
      BL_DEBUGF("ERROR: seq xor is not reversable with revcomp at iter %lu:\n\tkmer %s\n\trevcomp %s", i, km.toAlphabetString().c_str(), rev.toAlphabetString().c_str());
    }
    same &= local_same;

	local_same = (rev == op(std::make_pair(op(std::make_pair(rev, km)), km)));
    if (!local_same) {
      BL_DEBUGF("ERROR: seq xor is not reversable with km at iter %lu:\n\tkmer %s\n\trevcomp %s", i, km.toAlphabetString().c_str(), rev.toAlphabetString().c_str());
    }
    same &= local_same;

    // check single operand is same as double operand.
    local_same = (op(std::make_pair(km, rev)) == op(km));
    if (!local_same) {
      BL_DEBUGF("ERROR: seq xor single vs double operands are not the same at iter %lu:\n\tkmer %s\n\trevcomp %s", i, km.toAlphabetString().c_str(), rev.toAlphabetString().c_str());
    }
    same &= local_same;

    gold.nextFromChar(rand() % TypeParam::KmerAlphabet::SIZE);
    km = gold;
    rev = km.reverse_complement();
  }
  EXPECT_TRUE(same);
}


TYPED_TEST_P(KmerTransformTest, lex_less)
{

  auto km = this->kmer;
  auto rev = km;
  auto gold = km;


  bliss::kmer::transform::lex_less<TypeParam> op;

  auto tmp = km;
  bool same = true;
  bool local_same;

  for (size_t i = 0; i < this->iterations; ++i) {
	rev = km.reverse_complement();
	gold = km < rev ? km : rev;
	tmp = op(km);

    local_same = (tmp == gold);
    if (!local_same) {
      BL_DEBUGF("ERROR: seq lex_less diff at iter %lu:\n\tinput %s\n\trev %s\n\tresult %s\n\tgold %s", i, km.toAlphabetString().c_str(), rev.toAlphabetString().c_str(), tmp.toAlphabetString().c_str(), gold.toAlphabetString().c_str());
    }
    same &= local_same;

    km.nextFromChar(rand() % TypeParam::KmerAlphabet::SIZE);
  }
  EXPECT_TRUE(same);
}

TYPED_TEST_P(KmerTransformTest, lex_greater)
{

  auto km = this->kmer;
  auto rev = km;
  auto gold = km;


  bliss::kmer::transform::lex_greater<TypeParam> op;

  auto tmp = km;
  bool same = true;
  bool local_same;

  for (size_t i = 0; i < this->iterations; ++i) {
	rev = km.reverse_complement();
	gold = km > rev ? km : rev;
	tmp = op(km);

    local_same = (tmp == gold);
    if (!local_same) {
      BL_DEBUGF("ERROR: seq lex_greater diff at iter %lu:\n\tinput %s\n\trev %s\n\tresult %s\n\tgold %s", i, km.toAlphabetString().c_str(), rev.toAlphabetString().c_str(), tmp.toAlphabetString().c_str(), gold.toAlphabetString().c_str());
    }
    same &= local_same;

    km.nextFromChar(rand() % TypeParam::KmerAlphabet::SIZE);
  }
  EXPECT_TRUE(same);
}




REGISTER_TYPED_TEST_CASE_P(KmerTransformTest, identity, trans_xor, lex_less, lex_greater);

//////////////////// RUN the tests with different types.

// max of 50 cases
typedef ::testing::Types<
    ::bliss::common::Kmer< 31, bliss::common::DNA,   uint64_t>,  // 1 word, not full
    ::bliss::common::Kmer< 32, bliss::common::DNA,   uint64_t>,  // 1 word, full
    ::bliss::common::Kmer< 64, bliss::common::DNA,   uint64_t>,  // 2 words, full
    ::bliss::common::Kmer< 80, bliss::common::DNA,   uint64_t>,  // 3 words, not full
    ::bliss::common::Kmer< 96, bliss::common::DNA,   uint64_t>,  // 3 words, full
    ::bliss::common::Kmer< 15, bliss::common::DNA,   uint32_t>,  // 1 word, not full
    ::bliss::common::Kmer< 16, bliss::common::DNA,   uint32_t>,  // 1 word, full
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
    ::bliss::common::Kmer< 32, bliss::common::DNA16, uint64_t>,  // 2 words, full
    ::bliss::common::Kmer< 40, bliss::common::DNA16, uint64_t>,  // 3 words, not full
    ::bliss::common::Kmer< 48, bliss::common::DNA16, uint64_t>,  // 3 words, full
    ::bliss::common::Kmer<  7, bliss::common::DNA16, uint32_t>,  // 1 word, not full
    ::bliss::common::Kmer<  8, bliss::common::DNA16, uint32_t>,  // 1 word, full
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
    ::bliss::common::Kmer<  4, bliss::common::DNA16,  uint8_t>,   // 2 words, full
    ::bliss::common::Kmer< 21, bliss::common::ASCII,  uint64_t>,  // 3 words, not full
    ::bliss::common::Kmer< 21, bliss::common::ASCII,  uint8_t>  // 3 words, not full
> KmerTransformTestTypes;
INSTANTIATE_TYPED_TEST_CASE_P(Bliss, KmerTransformTest, KmerTransformTestTypes);

