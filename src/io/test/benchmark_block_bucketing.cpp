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


// include google test
#include <gtest/gtest.h>
#include "io/incremental_mxx.hpp"

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

struct BlockBucketBenchmarkInfo {
    size_t input_size;
    size_t bucket_count;
    size_t block_size;

    BlockBucketBenchmarkInfo() = default;
    BlockBucketBenchmarkInfo(size_t const & _input_size, size_t const& _bucket_count, size_t const & _block_size) :
      input_size(_input_size), bucket_count(_bucket_count), block_size(_block_size) {};

    BlockBucketBenchmarkInfo(BlockBucketBenchmarkInfo const & other) = default;
    BlockBucketBenchmarkInfo& operator=(BlockBucketBenchmarkInfo const & other) = default;
    BlockBucketBenchmarkInfo(BlockBucketBenchmarkInfo && other) = default;
    BlockBucketBenchmarkInfo& operator=(BlockBucketBenchmarkInfo && other) = default;
};


/**
 * @brief test fixture for file parsers, which segments file into sequences.
 */
class BlockBucketBenchmark : public ::testing::TestWithParam<BlockBucketBenchmarkInfo>
{
  protected:
    BlockBucketBenchmark() {};
    virtual ~BlockBucketBenchmark() {};

    virtual void SetUp()
    {
      p = GetParam();
      //std::cout << "p " << p.input_size << "," << p.bucket_count << "," << p.block_size << std::endl;

      mapping.clear();
      bcounts.clear();
      bucketed.clear();
      unbucketed.clear();

      // allocate test data
      data.clear();
      data.reserve(p.input_size);

      // populate with random numbers
      srand(17);
      size_t val;
      for (size_t i = 0; i < p.input_size; ++i) {
        val = rand();
        val <<= 32;
        val |= rand();
        data.emplace_back(val, i);
      }

    }

    virtual void TearDown() {
    }

    BlockBucketBenchmarkInfo p;
    using T = std::pair<size_t, size_t>;

    // input data
    std::vector<T> data;

    // bucket counts
    std::vector<size_t> bcounts;

    // one of 2 possible outputs.
    std::vector<size_t> mapping;
    std::vector<T> bucketed;
    std::vector<T> unbucketed;

};

TEST_P(BlockBucketBenchmark, assign)
{
  this->unbucketed.clear();
  this->bucketed.clear();

  // allocate.
  this->mapping.reserve(this->p.input_size);

  BlockBucketBenchmarkInfo pp = this->p;

  imxx::local::assign_to_buckets(this->data, [&pp](std::pair<size_t, size_t> const & x){ return x.first % pp.bucket_count; },
                              this->p.bucket_count,
                                 this->bcounts, this->mapping);

}

TEST_P(BlockBucketBenchmark, make_permutation)
{
	this->unbucketed.clear();
  this->bucketed.clear();

  // allocate.
  this->mapping.reserve(this->p.input_size);

  BlockBucketBenchmarkInfo pp = this->p;

  imxx::local::assign_to_buckets(this->data, [&pp](std::pair<size_t, size_t> const & x){ return x.first % pp.bucket_count; },
		  	  	  	  	  	  	  this->p.bucket_count,
                                 this->bcounts, this->mapping);

  size_t min_bucket = *(std::min_element(this->bcounts.begin(), this->bcounts.end()));
  size_t nblocks = p.block_size > 0 ? (min_bucket / p.block_size) : 0;
  imxx::local::bucket_to_block_permutation(this->p.block_size, nblocks, this->bcounts, this->mapping);

}


TEST_P(BlockBucketBenchmark, permute)
{
  this->mapping.clear();
  this->bucketed.clear();
  this->unbucketed.clear();

  BlockBucketBenchmarkInfo pp = this->p;

  // allocate.
  this->mapping.reserve(this->p.input_size);


  imxx::local::assign_to_buckets(this->data, [&pp](std::pair<size_t, size_t> const & x){ return x.first % pp.bucket_count; },
		  	  	  	  	  	  	  this->p.bucket_count,
                                 this->bcounts, this->mapping);

  size_t min_bucket = *(std::min_element(this->bcounts.begin(), this->bcounts.end()));
  size_t nblocks = p.block_size > 0 ? (min_bucket / p.block_size) : 0;
  imxx::local::bucket_to_block_permutation(this->p.block_size, nblocks, this->bcounts, this->mapping);


  if (pp.bucket_count > 0) {

	  // allocate.
	  this->bucketed.resize(this->p.input_size);

	  imxx::local::permute(this->data,
						   this->mapping,
							  this->bucketed);

  }
}

TEST_P(BlockBucketBenchmark, inplace_permute)
{
	this->unbucketed.clear();
	  this->mapping.clear();
	  this->bucketed.clear();
	  // allocate.
	  this->mapping.reserve(this->p.input_size);

	  BlockBucketBenchmarkInfo pp = this->p;

	  imxx::local::assign_to_buckets(this->data, [&pp](std::pair<size_t, size_t> const & x){ return x.first % pp.bucket_count; },
			  this->p.bucket_count,
	                                 this->bcounts, this->mapping);

    size_t min_bucket = *(std::min_element(this->bcounts.begin(), this->bcounts.end()));
    size_t nblocks = p.block_size > 0 ? (min_bucket / p.block_size) : 0;
    imxx::local::bucket_to_block_permutation(this->p.block_size, nblocks, this->bcounts, this->mapping);

	  if (pp.bucket_count > 0) {


		  // allocate.
		  this->bucketed.resize(this->p.input_size);

		  // copy
		  std::copy(this->data.begin(), this->data.end(), this->bucketed.begin());

		  imxx::local::permute_inplace(this->bucketed,
							   this->mapping);
	  }

}

TEST_P(BlockBucketBenchmark, unpermute)
{
	  this->mapping.clear();
	  this->bucketed.clear();
	  this->unbucketed.clear();

	  // allocate.
	  this->mapping.reserve(this->p.input_size);

	  BlockBucketBenchmarkInfo pp = this->p;

	  imxx::local::assign_to_buckets(this->data, [&pp](std::pair<size_t, size_t> const & x){ return x.first % pp.bucket_count; },
			  this->p.bucket_count,
	                                 this->bcounts, this->mapping);

    size_t min_bucket = *(std::min_element(this->bcounts.begin(), this->bcounts.end()));
    size_t nblocks = p.block_size > 0 ? (min_bucket / p.block_size) : 0;
    imxx::local::bucket_to_block_permutation(this->p.block_size, nblocks, this->bcounts, this->mapping);

	  if (pp.bucket_count > 0) {


		  // allocate.
		  this->bucketed.resize(this->p.input_size);

		  imxx::local::permute(this->data,
							   this->mapping, this->bucketed);



		  // allocate.
		  this->unbucketed.resize(this->p.input_size);

		  // copy

		  imxx::local::unpermute(this->bucketed,
							   this->mapping, this->unbucketed);
	  }
}

TEST_P(BlockBucketBenchmark, inplace_unpermute)
{
	  this->mapping.clear();
	  this->bucketed.clear();
	  this->unbucketed.clear();

	  // allocate.
	  this->mapping.reserve(this->p.input_size);

	  BlockBucketBenchmarkInfo pp = this->p;

	  imxx::local::assign_to_buckets(this->data, [&pp](std::pair<size_t, size_t> const & x){ return x.first % pp.bucket_count; },
	                                 this->p.bucket_count,
	                                 this->bcounts, this->mapping);


	  size_t min_bucket = *(std::min_element(this->bcounts.begin(), this->bcounts.end()));
	  size_t nblocks = p.block_size > 0 ? (min_bucket / p.block_size) : 0;
	  imxx::local::bucket_to_block_permutation(this->p.block_size, nblocks, this->bcounts, this->mapping);

	  if (pp.bucket_count > 0) {

		  // allocate.
		  this->bucketed.resize(this->p.input_size);

		  imxx::local::permute(this->data,
							   this->mapping, this->bucketed);

		  // allocate.
		  this->unbucketed.resize(this->p.input_size);

		  // copy
		  std::copy(this->bucketed.begin(), this->bucketed.end(), this->unbucketed.begin());


		  imxx::local::unpermute_inplace(this->unbucketed,
							   this->mapping);
	  }
}



INSTANTIATE_TEST_CASE_P(Bliss, BlockBucketBenchmark, ::testing::Values(
    BlockBucketBenchmarkInfo((1UL << 26), 1UL << 16, (1UL <<  8))  // 1, full
));










