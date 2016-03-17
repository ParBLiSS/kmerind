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
      BL_ERRORF("ERROR: seq transform diff at iter %lu:\n\tinput %s\n\toutput %s", i, gold.toAlphabetString().c_str(), km.toAlphabetString().c_str());
    }
    ASSERT_TRUE(local_same);
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
      BL_ERRORF("ERROR: seq xor is not symmetric at iter %lu:\n\tkmer %s\n\trevcomp %s", i, km.toAlphabetString().c_str(), rev.toAlphabetString().c_str());
    }
    ASSERT_TRUE(local_same);
    same &= local_same;

	// check reversable.
	local_same = (gold == op(std::make_pair(op(std::make_pair(km, rev)), rev)));
    if (!local_same) {
      BL_ERRORF("ERROR: seq xor is not reversable with revcomp at iter %lu:\n\tkmer %s\n\trevcomp %s", i, km.toAlphabetString().c_str(), rev.toAlphabetString().c_str());
    }
    ASSERT_TRUE(local_same);
    same &= local_same;

	local_same = (rev == op(std::make_pair(op(std::make_pair(rev, km)), km)));
    if (!local_same) {
      BL_ERRORF("ERROR: seq xor is not reversable with km at iter %lu:\n\tkmer %s\n\trevcomp %s", i, km.toAlphabetString().c_str(), rev.toAlphabetString().c_str());
    }
    ASSERT_TRUE(local_same);
    same &= local_same;

    // check single operand is same as double operand.
    local_same = (op(std::make_pair(km, rev)) == op(km));
    if (!local_same) {
      BL_ERRORF("ERROR: seq xor single vs double operands are not the same at iter %lu:\n\tkmer %s\n\trevcomp %s", i, km.toAlphabetString().c_str(), rev.toAlphabetString().c_str());
    }
    ASSERT_TRUE(local_same);
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
      BL_ERRORF("ERROR: seq lex_less diff at iter %lu:\n\tinput %s\n\trev %s\n\tresult %s\n\tgold %s", i, km.toAlphabetString().c_str(), rev.toAlphabetString().c_str(), tmp.toAlphabetString().c_str(), gold.toAlphabetString().c_str());
    }
    ASSERT_TRUE(local_same);
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
      BL_ERRORF("ERROR: seq lex_greater diff at iter %lu:\n\tinput %s\n\trev %s\n\tresult %s\n\tgold %s", i, km.toAlphabetString().c_str(), rev.toAlphabetString().c_str(), tmp.toAlphabetString().c_str(), gold.toAlphabetString().c_str());
    }
    ASSERT_TRUE(local_same);
    same &= local_same;

    km.nextFromChar(rand() % TypeParam::KmerAlphabet::SIZE);
  }
  EXPECT_TRUE(same);
}




REGISTER_TYPED_TEST_CASE_P(KmerTransformTest, identity, trans_xor, lex_less, lex_greater);

//////////////////// RUN the tests with different types.


