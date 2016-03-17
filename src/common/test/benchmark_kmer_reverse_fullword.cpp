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


#include "common/test/template_benchmark_kmer_reverse.cpp"

//////////////////// RUN the tests with different types.

// max of 50 cases
typedef ::testing::Types<
    ::bliss::common::Kmer< 32, bliss::common::DNA,   uint64_t>,
    ::bliss::common::Kmer< 64, bliss::common::DNA,   uint64_t>,
    ::bliss::common::Kmer< 96, bliss::common::DNA,   uint64_t>,
    ::bliss::common::Kmer<128, bliss::common::DNA,   uint64_t>,
    ::bliss::common::Kmer<256, bliss::common::DNA,   uint64_t>,
    ::bliss::common::Kmer< 32, bliss::common::DNA5,  uint64_t>,
    ::bliss::common::Kmer< 64, bliss::common::DNA5,  uint64_t>,
    ::bliss::common::Kmer< 96, bliss::common::DNA5,  uint64_t>,
    ::bliss::common::Kmer<128, bliss::common::DNA5,  uint64_t>,
    ::bliss::common::Kmer<256, bliss::common::DNA5,  uint64_t>,
    ::bliss::common::Kmer< 32, bliss::common::DNA16, uint64_t>,
    ::bliss::common::Kmer< 64, bliss::common::DNA16, uint64_t>,
    ::bliss::common::Kmer< 96, bliss::common::DNA16, uint64_t>,
    ::bliss::common::Kmer<128, bliss::common::DNA16, uint64_t>,
    ::bliss::common::Kmer<256, bliss::common::DNA16, uint64_t>
> KmerReverseBenchmarkTypes_FullWord;
INSTANTIATE_TYPED_TEST_CASE_P(Bliss_Kmer_FullWord,
		KmerReverseBenchmark,
		KmerReverseBenchmarkTypes_FullWord);

