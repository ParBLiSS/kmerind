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

#include "mxx/algos.hpp"

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

struct BucketBenchmarkInfo {
    size_t input_size;
    size_t bucket_count;
    size_t first;
    size_t last;

    BucketBenchmarkInfo() = default;
    BucketBenchmarkInfo(size_t const & _input_size, size_t const& _bucket_count, size_t const & _first, size_t const & _last) :
      input_size(_input_size), bucket_count(_bucket_count), first(_first), last(_last) {};

    BucketBenchmarkInfo(BucketBenchmarkInfo const & other) = default;
    BucketBenchmarkInfo& operator=(BucketBenchmarkInfo const & other) = default;
    BucketBenchmarkInfo(BucketBenchmarkInfo && other) = default;
    BucketBenchmarkInfo& operator=(BucketBenchmarkInfo && other) = default;
};


/**
 * @brief test fixture for file parsers, which segments file into sequences.
 */
class BucketBenchmark : public ::testing::TestWithParam<BucketBenchmarkInfo>
{
  protected:
    BucketBenchmark() {};
    virtual ~BucketBenchmark() {};

    virtual void SetUp()
    {
      p = GetParam();

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

      mapping.clear();
      bcounts.clear();
    }

    virtual void TearDown() {
    }

    BucketBenchmarkInfo p;
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


TEST_P(BucketBenchmark, assign)
{
	this->unbucketed.clear();
  this->bucketed.clear();

  // allocate.
  this->mapping.reserve(this->p.input_size);

  BucketBenchmarkInfo pp = this->p;

  imxx::local::assign_to_buckets(this->data, [&pp](std::pair<size_t, size_t> const & x){ return x.first % pp.bucket_count; },
		  this->p.bucket_count,
                                 this->bcounts, this->mapping, this->p.first, this->p.last);


}

TEST_P(BucketBenchmark, make_permutation)
{
  this->unbucketed.clear();
  this->bucketed.clear();

  // allocate.
  this->mapping.reserve(this->p.input_size);

  BucketBenchmarkInfo pp = this->p;

  imxx::local::assign_to_buckets(this->data, [&pp](std::pair<size_t, size_t> const & x){ return x.first % pp.bucket_count; },
      this->p.bucket_count,
                                 this->bcounts, this->mapping, this->p.first, this->p.last);

  imxx::local::bucket_to_permutation(this->bcounts, this->mapping, this->p.first, this->p.last);

}

TEST_P(BucketBenchmark, destructive_bucket)
{
	this->unbucketed.clear();
  this->mapping.clear();

  // allocate.
  this->bucketed.resize(this->p.input_size);

  // copy
  std::copy(this->data.begin(), this->data.end(), this->bucketed.begin());


  BucketBenchmarkInfo pp = this->p;

  imxx::local::bucketing_impl(this->bucketed, [&pp](std::pair<size_t, size_t> const & x){ return x.first % pp.bucket_count; },
		  this->p.bucket_count,
                           this->bcounts, this->p.first, this->p.last);

}

TEST_P(BucketBenchmark, bucket)
{
	this->bcounts.clear();
	this->unbucketed.clear();
  this->mapping.clear();

  // allocate.
  this->bucketed.resize(this->p.input_size);

  BucketBenchmarkInfo pp = this->p;

  imxx::local::bucketing_impl(this->data, [&pp](std::pair<size_t, size_t> const & x){ return x.first % pp.bucket_count; },
		  this->p.bucket_count, this->bcounts,   this->bucketed,
                                   this->p.first, this->p.last);
}

TEST_P(BucketBenchmark, mxx_inplace_bucket)
{
  this->unbucketed.clear();
  this->mapping.clear();

  // allocate.
  this->bucketed.resize(this->p.input_size);

  // copy
  std::copy(this->data.begin(), this->data.end(), this->bucketed.begin());

  BucketBenchmarkInfo pp = this->p;

  ::mxx::bucketing_inplace(this->bucketed, [&pp](std::pair<size_t, size_t> const & x){ return x.first % pp.bucket_count; },
		  this->p.bucket_count).swap(this->bcounts);

}

TEST_P(BucketBenchmark, mxx_bucket)
{
	this->bcounts.clear();
	this->unbucketed.clear();
  this->mapping.clear();

  // allocate.
  this->bucketed.resize(this->p.input_size);
  // copy
  std::copy(this->data.begin(), this->data.end(), this->bucketed.begin());

  BucketBenchmarkInfo pp = this->p;

  ::mxx::bucketing(this->bucketed, [&pp](std::pair<size_t, size_t> const & x){ return x.first % pp.bucket_count; },
		  this->p.bucket_count).swap(this->bcounts);
}


TEST_P(BucketBenchmark, permute)
{
	this->unbucketed.clear();
  this->mapping.clear();
  this->bucketed.clear();

  BucketBenchmarkInfo pp = this->p;

  // allocate.
  this->mapping.reserve(this->p.input_size);

  imxx::local::assign_to_buckets(this->data, [&pp](std::pair<size_t, size_t> const & x){ return x.first % pp.bucket_count; },
		  this->p.bucket_count,
                                 this->bcounts, this->mapping, this->p.first, this->p.last);

  imxx::local::bucket_to_permutation(this->bcounts, this->mapping, this->p.first, this->p.last);


  if (pp.bucket_count > 0) {

	  // allocate.
	  this->bucketed.resize(this->p.input_size);

	  imxx::local::permute(this->data.begin() + this->p.first, this->data.begin() + this->p.last,
						   this->mapping.begin() + this->p.first,
							  this->bucketed.begin() + this->p.first,
						   this->p.first);

  }
}

TEST_P(BucketBenchmark, inplace_permute)
{
	this->unbucketed.clear();
	  this->mapping.clear();
	  this->bucketed.clear();
	  // allocate.
	  this->mapping.reserve(this->p.input_size);

	  BucketBenchmarkInfo pp = this->p;

	  imxx::local::assign_to_buckets(this->data, [&pp](std::pair<size_t, size_t> const & x){ return x.first % pp.bucket_count; },
			  this->p.bucket_count,
	                                 this->bcounts, this->mapping, this->p.first, this->p.last);

	  imxx::local::bucket_to_permutation(this->bcounts, this->mapping, this->p.first, this->p.last);

	  if (pp.bucket_count > 0) {


		  // allocate.
		  this->bucketed.resize(this->p.input_size);

		  // copy
		  std::copy(this->data.begin(), this->data.end(), this->bucketed.begin());

		  imxx::local::permute_inplace(this->bucketed,
							   this->mapping, this->p.first, this->p.last);
	  }

}

TEST_P(BucketBenchmark, unpermute)
{
	  this->mapping.clear();
	  this->bucketed.clear();
	  this->unbucketed.clear();

	  // allocate.
	  this->mapping.reserve(this->p.input_size);

	  BucketBenchmarkInfo pp = this->p;

	  imxx::local::assign_to_buckets(this->data, [&pp](std::pair<size_t, size_t> const & x){ return x.first % pp.bucket_count; },
			  this->p.bucket_count,
	                                 this->bcounts, this->mapping, this->p.first, this->p.last);

	  imxx::local::bucket_to_permutation(this->bcounts, this->mapping, this->p.first, this->p.last);



	  if (pp.bucket_count > 0) {


		  // allocate.
		  this->bucketed.resize(this->p.input_size);

		  imxx::local::permute(this->data.begin() + this->p.first, this->data.begin() + this->p.last,
							   this->mapping.begin() + this->p.first,
								  this->bucketed.begin() + this->p.first,
							   this->p.first);



		  // allocate.
		  this->unbucketed.resize(this->p.input_size);

		  // copy
		  imxx::local::unpermute(this->bucketed.begin() + this->p.first, this->bucketed.begin() + this->p.last,
							   this->mapping.begin() + this->p.first,
							   this->unbucketed.begin() + this->p.first,
							   this->p.first);

	  }
}

TEST_P(BucketBenchmark, inplace_unpermute)
{
	  this->mapping.clear();
	  this->bucketed.clear();
	  this->unbucketed.clear();

	  // allocate.
	  this->mapping.reserve(this->p.input_size);

	  BucketBenchmarkInfo pp = this->p;

	  imxx::local::assign_to_buckets(this->data, [&pp](std::pair<size_t, size_t> const & x){ return x.first % pp.bucket_count; },
			  this->p.bucket_count,
	                                 this->bcounts, this->mapping, this->p.first, this->p.last);

	  imxx::local::bucket_to_permutation(this->bcounts, this->mapping, this->p.first, this->p.last);

	  if (pp.bucket_count > 0) {

		  // allocate.
		  this->bucketed.resize(this->p.input_size);


		  imxx::local::permute(this->data.begin() + this->p.first, this->data.begin() + this->p.last,
							   this->mapping.begin() + this->p.first,
								  this->bucketed.begin() + this->p.first,
							   this->p.first);


		  // allocate.
		  this->unbucketed.resize(this->p.input_size);

		  // copy
		  std::copy(this->bucketed.begin(), this->bucketed.end(), this->unbucketed.begin());


		  imxx::local::unpermute_inplace(this->unbucketed,
							   this->mapping, this->p.first, this->p.last);
	  }
}



INSTANTIATE_TEST_CASE_P(Bliss, BucketBenchmark, ::testing::Values(
    BucketBenchmarkInfo((1UL << 22), 1UL << 16, 0, (1UL << 22))  // 1, full

));










