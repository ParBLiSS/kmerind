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
 * mpi_test_distribute.cpp
 *   No separate MPI test.  FileLoader has no logical differences between MPI and non-MPI versions.  Only difference is MPI_comm is stored for convenience.
 *   note that this may change in the future if MPI_IO is used.
 *
 *  Created on: Feb 18, 2014
 *      Author: Tony Pan <tpan7@gatech.edu>
 */

// include google test
#include <gtest/gtest.h>
#include "bliss-config.hpp"

#if defined(USE_MPI)
#include "mxx/env.hpp"
#include "mxx/algos.hpp"
#include "mxx/collective.hpp"
#include "mxx/samplesort.hpp"



#include <cstdint> // for uint64_t, etc.
#include <type_traits>  // for integral_constant

#include <io/incremental_mxx.hpp>

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <random>
#include <iterator>  // for ostream iterator.
#include <algorithm>  // for transform.
#include <random>
#include <cstdint>  // uint32_t
#include <utility>  // pair
#include <vector>




//=========================  DISTRIBUTE TESTS


struct SamplesortBenchmarkInfo {
    size_t input_size;
    bool stable;

    SamplesortBenchmarkInfo() = default;
    SamplesortBenchmarkInfo(size_t const & _input_size, bool const & _stable) :
      input_size(_input_size), stable(_stable) {};

    SamplesortBenchmarkInfo(SamplesortBenchmarkInfo const & other) = default;
    SamplesortBenchmarkInfo& operator=(SamplesortBenchmarkInfo const & other) = default;
    SamplesortBenchmarkInfo(SamplesortBenchmarkInfo && other) = default;
    SamplesortBenchmarkInfo& operator=(SamplesortBenchmarkInfo && other) = default;
};


/**
 * @brief test fixture for file parsers, which segments file into sequences.
 */
class SamplesortBenchmark : public ::testing::TestWithParam<SamplesortBenchmarkInfo>
{
  protected:
    SamplesortBenchmark() {};
    virtual ~SamplesortBenchmark() {};

    using T = std::pair<size_t, int>;

    virtual void SetUp()
    {
      p = GetParam();
      //std::cout << "p " << p.input_size << "," << p.bucket_count << "," << p.block_size << std::endl;

      sorted.clear();

      // allocate test data
      data.clear();

    }

    void init(mxx::comm const & comm) {

      // populate data with random numbers
      srand((comm.rank() + 1) * (comm.rank() + 1) - 1);

      size_t count = (p.input_size + comm.size() - 1) / comm.size();

      data.reserve(count);
      sorted.reserve(count);


      size_t val;
      for (size_t i = 0; i < count; ++i) {
        val = rand();
        val <<= 32;
        val |= rand();
        data.emplace_back(val, comm.rank());
      }

      _comm = comm.copy();
    }


    ::mxx::comm _comm;


    SamplesortBenchmarkInfo p;

    // input data
    std::vector<T> data;

    // sorted sorted
    std::vector<T> sorted;

};



TEST_P(SamplesortBenchmark, samplesort)
{

  ::mxx::comm comm;

  this->init(comm);

//  for (int i = 0; i < 1; ++i) 
  if (this->p.stable)
	  auto splitters = imxx::samplesort<true>(this->data, this->sorted, [](const T& x, const T& y){ return x.first < y.first; }, comm);
  else
	  auto splitters = imxx::samplesort<false>(this->data, this->sorted, [](const T& x, const T& y){ return x.first < y.first; }, comm);

}



TEST_P(SamplesortBenchmark, samplesort_buf)
{

  ::mxx::comm comm;

  this->init(comm);

//  for (int i = 0; i < 1; ++i)
  if (this->p.stable)
	auto splitters = imxx::samplesort_buf<true>(this->data, this->sorted, [](const T& x, const T& y){ return x.first < y.first; }, comm, 0, 0);
  else
	  auto splitters = imxx::samplesort_buf<false>(this->data, this->sorted, [](const T& x, const T& y){ return x.first < y.first; }, comm, 0, 0);
}
TEST_P(SamplesortBenchmark, samplesort_buf_sort_full)
{

  ::mxx::comm comm;

  this->init(comm);

//  for (int i = 0; i < 1; ++i)
  if (this->p.stable)
	auto splitters = imxx::samplesort_buf<true>(this->data, this->sorted, [](const T& x, const T& y){ return x.first < y.first; }, comm, 1,1);
  else
	  auto splitters = imxx::samplesort_buf<false>(this->data, this->sorted, [](const T& x, const T& y){ return x.first < y.first; }, comm, 1,1);
}
TEST_P(SamplesortBenchmark, samplesort_buf_sort_part)
{

  ::mxx::comm comm;

  this->init(comm);

//  for (int i = 0; i < 1; ++i)
  if (this->p.stable)
	auto splitters = imxx::samplesort_buf<true>(this->data, this->sorted, [](const T& x, const T& y){ return x.first < y.first; }, comm, 1,-1);
  else
	  auto splitters = imxx::samplesort_buf<false>(this->data, this->sorted, [](const T& x, const T& y){ return x.first < y.first; }, comm, 1,-1);
}
TEST_P(SamplesortBenchmark, samplesort_buf_merge_full)
{

  ::mxx::comm comm;

  this->init(comm);

//  for (int i = 0; i < 1; ++i)
  if (this->p.stable)
	auto splitters = imxx::samplesort_buf<true>(this->data, this->sorted, [](const T& x, const T& y){ return x.first < y.first; }, comm, -1,1);
  else
	  auto splitters = imxx::samplesort_buf<false>(this->data, this->sorted, [](const T& x, const T& y){ return x.first < y.first; }, comm, -1,1);
}
TEST_P(SamplesortBenchmark, samplesort_buf_merge_part)
{

  ::mxx::comm comm;

  this->init(comm);

//  for (int i = 0; i < 1; ++i)
  if (this->p.stable)
	auto splitters = imxx::samplesort_buf<true>(this->data, this->sorted, [](const T& x, const T& y){ return x.first < y.first; }, comm, -1,-1);
  else
	  auto splitters = imxx::samplesort_buf<false>(this->data, this->sorted, [](const T& x, const T& y){ return x.first < y.first; }, comm, -1,-1);
}




TEST_P(SamplesortBenchmark, mxx_samplesort)
{

  ::mxx::comm comm;

  this->init(comm);

  auto comp = [](const T& x, const T& y){ return x.first < y.first; };
  bool empty = mxx::all_of(this->data.size() == 0, comm);

  if (!empty) {

  this->sorted.resize(this->data.size());
  std::copy(this->data.begin(), this->data.end(), this->sorted.begin());

  if (this->p.stable)
	  ::mxx::impl::samplesort<decltype(this->sorted.begin()), decltype(comp), true>(this->sorted.begin(), this->sorted.end(), comp, comm);
  else
	  ::mxx::impl::samplesort<decltype(this->sorted.begin()), decltype(comp), false>(this->sorted.begin(), this->sorted.end(), comp, comm);
  }
}




INSTANTIATE_TEST_CASE_P(Bliss, SamplesortBenchmark, ::testing::Values(
    // base cases
	    SamplesortBenchmarkInfo((1UL << 24), false),  // 1
	    SamplesortBenchmarkInfo((1UL << 24), true)  // 1

));










#endif

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

#if defined(USE_MPI)
  ::mxx::env e(argc, argv);
  ::mxx::comm comm;

  result = RUN_ALL_TESTS();

  comm.barrier();
#endif

  return result;
}


