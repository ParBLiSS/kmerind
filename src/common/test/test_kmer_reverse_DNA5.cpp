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



#include "common/test/kmer_reverse_test_fixture.cpp"

//////////////////// RUN the tests with different types.
// max of 50 cases
typedef ::testing::Types<
     ::bliss::common::Kmer< 21, bliss::common::DNA5,  uint64_t>,  // 1 word, not full
     ::bliss::common::Kmer< 22, bliss::common::DNA5,  uint64_t>,  // 2 word, not full
     ::bliss::common::Kmer< 42, bliss::common::DNA5,  uint64_t>,  // 2 words, not full
     ::bliss::common::Kmer< 43, bliss::common::DNA5,  uint64_t>,  // 3 words, not full
     ::bliss::common::Kmer< 64, bliss::common::DNA5,  uint64_t>,  // 3 words, full
     ::bliss::common::Kmer<  2, bliss::common::DNA5,   uint8_t>,  // 1 word, not full
     ::bliss::common::Kmer<  3, bliss::common::DNA5,   uint8_t>,  // 2 word, not full
     ::bliss::common::Kmer<  5, bliss::common::DNA5,   uint8_t>,  // 2 words, not full
     ::bliss::common::Kmer<  6, bliss::common::DNA5,   uint8_t>,  // 3 words, not full
     ::bliss::common::Kmer<  8, bliss::common::DNA5,   uint8_t>  // 3 words, full
> KmerReverseTestDNA5Types;
INSTANTIATE_TYPED_TEST_CASE_P(Bliss_DNA5, KmerReverseTest, KmerReverseTestDNA5Types);


