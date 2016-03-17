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

//    ::bliss::common::Kmer< 15, bliss::common::DNA_IUPAC,   uint64_t>,  // 1 word, not full.  not tested since it is not a bijection.
//    ::bliss::common::Kmer< 40, bliss::common::DNA_IUPAC,   uint64_t>,  // 3 word, not full
    ::bliss::common::Kmer< 21, bliss::common::ASCII,  uint64_t>,  // 3 words, not full
    ::bliss::common::Kmer< 21, bliss::common::ASCII,  uint8_t>  // 3 words, not full
> KmerReverseTestOtherTypes;
INSTANTIATE_TYPED_TEST_CASE_P(Bliss_other, KmerReverseTest, KmerReverseTestOtherTypes);


