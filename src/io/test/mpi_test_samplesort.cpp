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


struct SamplesortTestInfo {
    size_t input_size;
    bool stable;

    SamplesortTestInfo() = default;
    SamplesortTestInfo(size_t const & _input_size, bool const & _stable) :
      input_size(_input_size), stable(_stable) {};

    SamplesortTestInfo(SamplesortTestInfo const & other) = default;
    SamplesortTestInfo& operator=(SamplesortTestInfo const & other) = default;
    SamplesortTestInfo(SamplesortTestInfo && other) = default;
    SamplesortTestInfo& operator=(SamplesortTestInfo && other) = default;
};


/**
 * @brief test fixture for file parsers, which segments file into sequences.
 */
class SamplesortTest : public ::testing::TestWithParam<SamplesortTestInfo>
{
  protected:
    SamplesortTest() {};
    virtual ~SamplesortTest() {};

    using T = std::pair<size_t, int>;

    virtual void SetUp()
    {
      p = GetParam();
      //std::cout << "p " << p.input_size << "," << p.bucket_count << "," << p.block_size << std::endl;

      sorted.clear();
      sorted.reserve(p.input_size);

      // allocate test data
      data.clear();
      data.reserve(p.input_size);

      gold.clear();

    }

    void init(mxx::comm const & comm) {

      // populate data with random numbers
      srand((comm.rank() + 1) * (comm.rank() + 1) - 1);
      size_t val;
      for (size_t i = 0; i < p.input_size; ++i) {
        val = rand();
        val <<= 32;
        val |= rand();
        data.emplace_back(val, comm.rank());
      }

      // make a copy of data.
      mxx::gatherv(data, 0, comm).swap(gold);

      if (comm.rank() == 0) {
    	  if (p.stable)
    		  std::stable_sort(gold.begin(), gold.end(), [](const T& x, const T& y){ return x.first < y.first; });
    	  else
    		  std::sort(gold.begin(), gold.end(), [](const T& x, const T& y){ return x.first < y.first; });
      }
      _comm = comm.copy();
    }

    virtual void TearDown() {

      std::cout << "comm " << _comm.rank() <<  " data size = " << data.size() << " sorted size = " << sorted.size() << std::endl;

      std::vector<T> temp;
      mxx::gatherv(sorted, 0, _comm).swap(temp);

	bool same = true;
      if (_comm.rank() == 0) {
	std::cout << "comm " << _comm.rank() << " gold size = " << gold.size() << " temp size = " << temp.size() << std::endl;

	// check sorted same as gold
//      EXPECT_EQ(temp.size(), gold.size());
      same = (temp.size() == gold.size());
      for (size_t i = 0; i < temp.size(); ++i) {
          same &= (temp[i] == gold[i]);

          if ((temp[i] != gold[i]) && ((i < 10) || ( i > temp.size() - 10))) {
            std::cout << "ERROR comm " << _comm.rank() <<  " pos " << i << " sorted " << temp[i].first << "," << temp[i].second << "; gold " << gold[i].first << "," << gold[i].second << std::endl;
//            EXPECT_EQ(temp[i], gold[i]);
          }

        }

    }

	same = ::mxx::all_of(same, _comm);
        EXPECT_TRUE(same);
	}
    ::mxx::comm _comm;


    SamplesortTestInfo p;

    // input data
    std::vector<T> data;

    // sorted sorted
    std::vector<T> sorted;

    // gold standard results
    std::vector<T> gold;
};



TEST_P(SamplesortTest, samplesort)
{

  ::mxx::comm comm;

  this->init(comm);

//  if (this->p.use_sort_override == 0 && this->p.full_buffer_override == 0) {

  if (this->p.stable)
	  imxx::samplesort<true>(this->data, this->sorted, [](const T& x, const T& y){ return x.first < y.first; }, comm);
  else
	  imxx::samplesort<false>(this->data, this->sorted, [](const T& x, const T& y){ return x.first < y.first; }, comm);
//  }
}


TEST_P(SamplesortTest, samplesort_buf)
{

  ::mxx::comm comm;

  this->init(comm);

  if (this->p.stable)
	  imxx::samplesort_buf<true>(this->data, this->sorted, [](const T& x, const T& y){ return x.first < y.first; }, comm, 0, 0);
  else
	  imxx::samplesort_buf<false>(this->data, this->sorted, [](const T& x, const T& y){ return x.first < y.first; }, comm, 0, 0);

}
TEST_P(SamplesortTest, samplesort_buf_sort_full)
{

  ::mxx::comm comm;

  this->init(comm);

  if (this->p.stable)
	  imxx::samplesort_buf<true>(this->data, this->sorted, [](const T& x, const T& y){ return x.first < y.first; }, comm, 1, 1);
  else
	  imxx::samplesort_buf<false>(this->data, this->sorted, [](const T& x, const T& y){ return x.first < y.first; }, comm, 1, 1);

}
TEST_P(SamplesortTest, samplesort_buf_sort_part)
{

  ::mxx::comm comm;

  this->init(comm);

  if (this->p.stable)
	  imxx::samplesort_buf<true>(this->data, this->sorted, [](const T& x, const T& y){ return x.first < y.first; }, comm, 1, -1);
  else
	  imxx::samplesort_buf<false>(this->data, this->sorted, [](const T& x, const T& y){ return x.first < y.first; }, comm, 1, -1);

}
TEST_P(SamplesortTest, samplesort_buf_merge_full)
{

  ::mxx::comm comm;

  this->init(comm);

  if (this->p.stable)
	  imxx::samplesort_buf<true>(this->data, this->sorted, [](const T& x, const T& y){ return x.first < y.first; }, comm, -1, 1);
  else
	  imxx::samplesort_buf<false>(this->data, this->sorted, [](const T& x, const T& y){ return x.first < y.first; }, comm, -1, 1);

}
TEST_P(SamplesortTest, samplesort_buf_merge_part)
{

  ::mxx::comm comm;

  this->init(comm);

  if (this->p.stable)
	  imxx::samplesort_buf<true>(this->data, this->sorted, [](const T& x, const T& y){ return x.first < y.first; }, comm, -1, -1);
  else
	  imxx::samplesort_buf<false>(this->data, this->sorted, [](const T& x, const T& y){ return x.first < y.first; }, comm, -1, -1);

}





TEST_P(SamplesortTest, mxx_samplesort)
{

  ::mxx::comm comm;

  this->init(comm);

  this->sorted.resize(this->data.size());
  std::copy(this->data.begin(), this->data.end(), this->sorted.begin());

  auto comp = [](const T& x, const T& y){ return x.first < y.first; };

   bool empty = mxx::all_of(this->data.size() == 0, comm);
  if (!empty) {
	  if (this->p.stable)
		  ::mxx::impl::samplesort<decltype(this->sorted.begin()), decltype(comp), true>(this->sorted.begin(), this->sorted.end(), comp, comm);
	  else
		  ::mxx::impl::samplesort<decltype(this->sorted.begin()), decltype(comp), false>(this->sorted.begin(), this->sorted.end(), comp, comm);
  }
}



INSTANTIATE_TEST_CASE_P(Bliss, SamplesortTest, ::testing::Values(
    // base cases
    SamplesortTestInfo(0UL, false),   //  0, boundary case
    SamplesortTestInfo(0UL, true),   //  0, boundary case

    SamplesortTestInfo((1UL <<  8), false),  // 1
    SamplesortTestInfo((1UL <<  8), true)  // 1
//    SamplesortTestInfo((1UL << 11), false),  // 1
//    SamplesortTestInfo((1UL << 11), true),  // 1
//    SamplesortTestInfo((1UL << 12), false),  // 1
//    SamplesortTestInfo((1UL << 12), true)  // 1

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


