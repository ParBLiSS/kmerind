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
#include "index/kmer_hash.hpp"

#include <random>
#include <cstdint>
#include <vector>
#include <unordered_set>
#include <set>
#include <cmath>

#include "common/kmer.hpp"
#include "common/alphabets.hpp"
#include "common/alphabet_traits.hpp"



// include files to test

//TESTS: Hash functions.  test boundary cases - the number of unique outputs should be relatively close to number of unique inputs.

template <typename T>
class KmerHashTest : public ::testing::Test {
  protected:

    std::vector<T> kmers;
    std::set<T> unique_kmers;
    T kmer;

    static const size_t iterations = 100000;

    virtual void SetUp()
    {
      srand(0);
      for (unsigned int i = 0; i < T::size; ++i) {

          kmer.nextFromChar(rand() % T::KmerAlphabet::SIZE);
      }

      kmers.resize(iterations);
      for (size_t i = 0; i < iterations; ++i) {
    	  kmers[i] = kmer;
    	  unique_kmers.emplace(kmer);
          kmer.nextFromChar(rand() % T::KmerAlphabet::SIZE);
      }
    }

    template <template <typename, bool> class H>
    void hash_vector(std::string name) {
    	std::unordered_set<size_t> hashes;

       H<T, false> op;

      bool same = true;

      for (size_t i = 0; i < this->iterations; ++i) {
        hashes.emplace(op(this->kmers[i]));
      }

      //double limit = this->iterations;
      //double max_uniq_kmers = std::pow(static_cast<double>(T::KmerAlphabet::SIZE), static_cast<double>(T::size));
      //limit = (limit < max_uniq_kmers) ? limit : max_uniq_kmers;

//      same = (limit - hashes.size()) < (limit * 0.2);
      same = (hashes.size() == this->unique_kmers.size());
      //printf(" hash size: %lu, unique_kmers %lu, maxUnique %f, iterations %lu\n", hashes.size(), this->unique_kmers.size(), max_uniq_kmers, this->iterations);
      if (!same)
    	  BL_DEBUGF("ERROR: hash %s size %lu is not same as unique kmers: %lu.  input size %lu", name.c_str(), hashes.size(), this->unique_kmers.size(), this->kmers.size());

      EXPECT_TRUE(same);

    }
};


// indicate this is a typed test
TYPED_TEST_CASE_P(KmerHashTest);

TYPED_TEST_P(KmerHashTest, hash)
{
	this->template hash_vector<bliss::kmer::hash::cpp_std >(std::string("cpp_std"));
	this->template hash_vector<bliss::kmer::hash::identity>(std::string("identity"));
	this->template hash_vector<bliss::kmer::hash::murmur  >(std::string("murmur"));
	this->template hash_vector<bliss::kmer::hash::farm    >(std::string("farm"));
}






REGISTER_TYPED_TEST_CASE_P(KmerHashTest, hash);

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
    ::bliss::common::Kmer< 21, bliss::common::ASCII,  uint8_t>,  // 3 words, not full
    ::bliss::common::Kmer< 5120, bliss::common::ASCII,  uint64_t>  // 80 words - cpp_std should see collisions

> KmerHashTestTypes;
INSTANTIATE_TYPED_TEST_CASE_P(Bliss, KmerHashTest, KmerHashTestTypes);

