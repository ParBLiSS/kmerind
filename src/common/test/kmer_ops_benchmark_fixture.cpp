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
 */



#include "utils/logging.h"

// include google test
#include <gtest/gtest.h>
#include "utils/bitgroup_ops.hpp"

#include <random>
#include <cstdint>

#include "common/kmer.hpp"

#include "common/alphabets.hpp"
#include "common/alphabet_traits.hpp"
#include "utils/benchmark_utils.hpp"

// include files to test

//TESTS: Sequential, SWAR/BSWAP, SSSE3, AVX2 versions of kmer reverse.
//TESTS: for each, test kmer size, different word type, and Alphabet (bit group size)
//       different offsets, different bit group sizes, different word types, and different byte array lengths.

template <typename T>
class KmerOpsBenchmark : public ::testing::Test {
  protected:

    static constexpr size_t iterations = 1000000;

    static std::vector<T> kmers;
    static std::vector<T> kmers2;
    std::vector<T> outputs;

  public:
    static void SetUpTestCase()
    {
      kmers.resize(iterations);
      kmers2.resize(iterations);

      srand(23);
      for (size_t i = 0; i < iterations; ++i) {
        for (size_t j = 0; j < T::nWords; ++j) {
          kmers[i].getDataRef()[j] = static_cast<typename T::KmerWordType>(static_cast<long>(rand()) << 32) | static_cast<long>(rand());
          kmers2[i].getDataRef()[j] = static_cast<typename T::KmerWordType>(static_cast<long>(rand()) << 32) | static_cast<long>(rand());
        }
      }
    }

    static void TearDownTestCase() {
    	std::vector<T>().swap(kmers);
    	std::vector<T>().swap(kmers2);
    }

    virtual void SetUp() {
    	outputs.resize(iterations);
    }

    virtual void TearDown() {
    	std::vector<T>().swap(outputs);
    }


    inline bool old_equal(const T& lhs, const T& rhs) const {
		// slower.  for testing only
    	std::atomic_thread_fence(::std::memory_order_seq_cst);
		return (memcmp(lhs.getData(), rhs.getData(), T::nBytes) == 0);
    }

    inline bool old_less(const T& lhs, const T& rhs) const
    {
        auto first = lhs.getData() + T::nWords - 1;
        auto second = rhs.getData() + T::nWords - 1;

        for (size_t i = 0; i < T::nWords; ++i, --first, --second) {
          if (*first != *second) return (*first < *second);  // if equal, keep comparing. else decide.
        }
        return false;  // all equal
    }

    inline void old_xor(T & out, const T& lhs, const T& rhs)
    {
        for (uint32_t i = 0; i < T::nWords; ++i) out.getDataRef()[i] = lhs.getData()[i] ^ rhs.getData()[i];
    }
    inline void old_and(T & out, const T& lhs, const T& rhs)
    {
        for (uint32_t i = 0; i < T::nWords; ++i) out.getDataRef()[i] = lhs.getData()[i] & rhs.getData()[i];
    }
    inline void old_or(T & out, const T& lhs, const T& rhs)
    {
        for (uint32_t i = 0; i < T::nWords; ++i) out.getDataRef()[i] = lhs.getData()[i] | rhs.getData()[i];
    }


};

template <typename T>
constexpr size_t KmerOpsBenchmark<T>::iterations;
template <typename T>
std::vector<T> KmerOpsBenchmark<T>::kmers;
template <typename T>
std::vector<T> KmerOpsBenchmark<T>::kmers2;


// indicate this is a typed test
TYPED_TEST_CASE_P(KmerOpsBenchmark);

TYPED_TEST_P(KmerOpsBenchmark, left_shift)
{
  BL_TIMER_INIT(km);

  BL_TIMER_START(km);
  for (size_t i = 0; i < KmerOpsBenchmark<TypeParam>::iterations; ++i) {
	  this->outputs[i].left_shift_bits();
  }
	BL_TIMER_END(km, "bit << auto", KmerOpsBenchmark<TypeParam>::iterations);

	BL_TIMER_START(km);
  for (size_t i = 0; i < KmerOpsBenchmark<TypeParam>::iterations; ++i) {
	  this->outputs[i].left_shift_bits(TypeParam::bitsPerChar);
  }
  BL_TIMER_END(km, "word <<", KmerOpsBenchmark<TypeParam>::iterations);


  BL_TIMER_REPORT(km);
}

TYPED_TEST_P(KmerOpsBenchmark, right_shift)
{
  BL_TIMER_INIT(km);

    BL_TIMER_START(km);
    for (size_t i = 0; i < KmerOpsBenchmark<TypeParam>::iterations; ++i) {
    	  this->outputs[i].right_shift_bits();
    }
	BL_TIMER_END(km, "bit >> auto", KmerOpsBenchmark<TypeParam>::iterations);

	BL_TIMER_START(km);
	  for (size_t i = 0; i < KmerOpsBenchmark<TypeParam>::iterations; ++i) {
		  this->outputs[i].right_shift_bits(TypeParam::bitsPerChar);
	  }
	  BL_TIMER_END(km, "word >>", KmerOpsBenchmark<TypeParam>::iterations);


  BL_TIMER_REPORT(km);
}


TYPED_TEST_P(KmerOpsBenchmark, bit_and)
{
  BL_TIMER_INIT(km);

    BL_TIMER_START(km);
    for (size_t i = 0; i < KmerOpsBenchmark<TypeParam>::iterations; ++i) {
    	  this->outputs[i].bit_and(this->kmers[i], this->kmers2[i]);

    }
	BL_TIMER_END(km, "bit & auto", KmerOpsBenchmark<TypeParam>::iterations);

	  BL_TIMER_START(km);
	  for (size_t i = 0; i < KmerOpsBenchmark<TypeParam>::iterations; ++i) {
		  this->old_and(this->outputs[i], this->kmers[i], this->kmers2[i]);
	  }
	  BL_TIMER_END(km, "&", KmerOpsBenchmark<TypeParam>::iterations);

  BL_TIMER_REPORT(km);
}

TYPED_TEST_P(KmerOpsBenchmark, bit_or)
{
  BL_TIMER_INIT(km);

  BL_TIMER_START(km);
  for (size_t i = 0; i < KmerOpsBenchmark<TypeParam>::iterations; ++i) {
  	  this->outputs[i].bit_or(this->kmers[i], this->kmers2[i]);

  }
	BL_TIMER_END(km, "bit | auto", KmerOpsBenchmark<TypeParam>::iterations);

  BL_TIMER_START(km);
  for (size_t i = 0; i < KmerOpsBenchmark<TypeParam>::iterations; ++i) {
	  this->old_or(this->outputs[i], this->kmers[i], this->kmers2[i]);
  }
  BL_TIMER_END(km, "|", KmerOpsBenchmark<TypeParam>::iterations);

  BL_TIMER_REPORT(km);
}

TYPED_TEST_P(KmerOpsBenchmark, bit_xor)
{
  BL_TIMER_INIT(km);

  BL_TIMER_START(km);
  for (size_t i = 0; i < KmerOpsBenchmark<TypeParam>::iterations; ++i) {
	  this->outputs[i].bit_xor(this->kmers[i], this->kmers2[i]);
  }
	BL_TIMER_END(km, "bit^ auto", KmerOpsBenchmark<TypeParam>::iterations);

  BL_TIMER_START(km);
  for (size_t i = 0; i < KmerOpsBenchmark<TypeParam>::iterations; ++i) {
	  this->old_xor(this->outputs[i], this->kmers[i], this->kmers2[i]);
  }
  BL_TIMER_END(km, "^", KmerOpsBenchmark<TypeParam>::iterations);

  BL_TIMER_REPORT(km);
}

TYPED_TEST_P(KmerOpsBenchmark, equal)
{
  BL_TIMER_INIT(km);

  bool result = true;
  BL_TIMER_START(km);
  for (size_t i = 0; i < KmerOpsBenchmark<TypeParam>::iterations; ++i) {
	  result &= (this->kmers[i] == this->kmers2[i]);
  }
  BL_TIMER_END(km, "bit equal", KmerOpsBenchmark<TypeParam>::iterations);
   printf("equal? %s\n", (result ? "true" : "false"));

  bool result2 = true;
    BL_TIMER_START(km);
    for (size_t i = 0; i < KmerOpsBenchmark<TypeParam>::iterations; ++i) {
    	result2 &= (this->old_equal(this->kmers[i], this->kmers2[i]));
    }
	BL_TIMER_END(km, "equal", KmerOpsBenchmark<TypeParam>::iterations);
	  printf("equal? %s\n", (result2 ? "true" : "false"));

  BL_TIMER_REPORT(km);

  ASSERT_TRUE(result == result2);

}


TYPED_TEST_P(KmerOpsBenchmark, less)
{
  BL_TIMER_INIT(km);

  bool result = true;
  BL_TIMER_START(km);
  for (size_t i = 0; i < KmerOpsBenchmark<TypeParam>::iterations; ++i) {
	  result &= this->kmers[i] < this->kmers2[i];
  }
  BL_TIMER_END(km, "bit less", KmerOpsBenchmark<TypeParam>::iterations);
  // printf("less? %s\n", (result ? "true" : "false"));

  bool result2 = true;
    BL_TIMER_START(km);
    for (size_t i = 0; i < KmerOpsBenchmark<TypeParam>::iterations; ++i) {
    	result2 &= this->old_less(this->kmers[i], this->kmers2[i]);
    }
	BL_TIMER_END(km, "less", KmerOpsBenchmark<TypeParam>::iterations);
	//  printf("less? %s\n", (result2 ? "true" : "false"));

  BL_TIMER_REPORT(km);

  ASSERT_TRUE(result == result2);
}



//REGISTER_TYPED_TEST_CASE_P(KmerOpsBenchmark, rev_seq, rev_seq2, revcomp_seq, rev_bswap, revcomp_bswap, rev_swar, revcomp_swar, rev, revcomp, rev_ssse3, revcomp_ssse3);

REGISTER_TYPED_TEST_CASE_P(KmerOpsBenchmark,
		left_shift,
		right_shift,
		bit_and,
		bit_or,
		bit_xor,
		equal,
		less);

