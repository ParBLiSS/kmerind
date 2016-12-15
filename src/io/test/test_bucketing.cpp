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

struct BucketTestInfo {
    size_t input_size;
    size_t bucket_count;
    size_t first;
    size_t last;

    BucketTestInfo() = default;
    BucketTestInfo(size_t const & _input_size, size_t const& _bucket_count, size_t const & _first, size_t const & _last) :
      input_size(_input_size), bucket_count(_bucket_count), first(_first), last(_last) {};

    BucketTestInfo(BucketTestInfo const & other) = default;
    BucketTestInfo& operator=(BucketTestInfo const & other) = default;
    BucketTestInfo(BucketTestInfo && other) = default;
    BucketTestInfo& operator=(BucketTestInfo && other) = default;
};


/**
 * @brief test fixture for file parsers, which segments file into sequences.
 */
class BucketTest : public ::testing::TestWithParam<BucketTestInfo>
{
  protected:
    BucketTest() {};
    virtual ~BucketTest() {};

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

      // compute the gold - use identity function
      gold.clear();
      gold.resize(p.bucket_count);   // create bucket_count number of empty vectors.

      if ((p.bucket_count > 0) && (p.input_size > 0)) {
        size_t pid;
        for (size_t i = p.first; i < p.last; ++i) {
          pid = data[i].first % p.bucket_count;

          gold[pid].emplace_back(data[i]);
        }
      }

      mapping.clear();
      bcounts.clear();
    }

    virtual void TearDown() {

      // check that the bucket counts are consistent with gold
      bool same = true;
      for (size_t i = 0; i < p.bucket_count; ++i) {
        same &= (bcounts[i] == gold[i].size());

      }
      EXPECT_TRUE(same);

//      if (!same) {
//		  for (size_t i = 0; i < p.bucket_count; ++i) {
//				  std::cout << "i " << i << " test " << bcounts[i] << " gold " << gold[i].size() << std::endl;
//		  }
//      }

      same = true;
      size_t sid = p.first;
      if (mapping.size() > 0) {

        for (size_t i = 0; i < p.bucket_count; ++i) {

          for (size_t j = 0; j < gold[i].size(); ++j, ++sid) {
            //get the stored src id (in second field),  lookup its mapping, which should equal to the current position in output.
            same &= (mapping[gold[i][j].second] == sid) &&
                ((gold[i][j].first % p.bucket_count) == i);
          }
        }

      }
      EXPECT_TRUE(same);

      same = true;
      sid = p.first;
      if (bucketed.size() > 0) {

        for (size_t i = 0; i < p.bucket_count; ++i) {

          for (size_t j = 0; j < gold[i].size(); ++j, ++sid) {
            if (gold[i][j].first != bucketed[sid].first) {
              std::cout << "NO MATCH gold bucket " << i << " pos " << j << " " << gold[i][j].first << " from " << gold[i][j].second << "; bucketed id " << sid << " " << bucketed[sid].first << " from " << bucketed[sid].second << std::endl;
            }

            same &= (gold[i][j] == bucketed[sid]);
          }
        }

      }
      EXPECT_TRUE(same);

      same = true;
      if (unbucketed.size() > 0) {

          for (size_t i = p.first; i < p.last; ++i) {
			  if (unbucketed[i].first != data[i].first) {
				std::cout << "NO MATCH data " << i << " " << data[i].first << " from " << data[i].second << "; unbucketed " << unbucketed[i].first << " from " << bucketed[i].second << std::endl;
			  }
        	  same &= (unbucketed[i] == data[i]);
          }
      }
      EXPECT_TRUE(same);

    }

    BucketTestInfo p;
    using T = std::pair<size_t, size_t>;

    // input data
    std::vector<T> data;

    // bucket counts
    std::vector<size_t> bcounts;

    // one of 2 possible outputs.
    std::vector<size_t> mapping;
    std::vector<T> bucketed;
    std::vector<T> unbucketed;

    // gold standard results
    std::vector<std::vector<T>> gold;
};


TEST_P(BucketTest, assign)
{
	this->unbucketed.clear();
  this->bucketed.clear();

  // allocate.
  this->mapping.reserve(this->p.input_size);

  BucketTestInfo pp = this->p;

  imxx::local::assign_to_buckets(this->data, [&pp](std::pair<size_t, size_t> const & x){ return x.first % pp.bucket_count; },
		  this->p.bucket_count,
                                 this->bcounts, this->mapping, this->p.first, this->p.last);

  imxx::local::bucket_to_permutation(this->bcounts, this->mapping, this->p.first, this->p.last);

}

TEST_P(BucketTest, inplace_bucket)
{
	this->unbucketed.clear();
  this->mapping.clear();

  // allocate.
  this->bucketed.resize(this->p.input_size);

  // copy
  std::copy(this->data.begin(), this->data.end(), this->bucketed.begin());


  BucketTestInfo pp = this->p;

  imxx::local::bucketing_impl(this->bucketed, [&pp](std::pair<size_t, size_t> const & x){ return x.first % pp.bucket_count; },
		  this->p.bucket_count,
                           this->bcounts, this->p.first, this->p.last);

}

TEST_P(BucketTest, bucket)
{
	this->bcounts.clear();
	this->unbucketed.clear();
  this->mapping.clear();

  // allocate.
  this->bucketed.resize(this->p.input_size);

  BucketTestInfo pp = this->p;

  imxx::local::bucketing_impl(this->data, [&pp](std::pair<size_t, size_t> const & x){ return x.first % pp.bucket_count; },
		  this->p.bucket_count, this->bcounts,   this->bucketed,
                                   this->p.first, this->p.last);
}

TEST_P(BucketTest, mxx_bucket)
{
	this->bcounts.clear();
	this->unbucketed.clear();
  this->mapping.clear();

  // allocate.
  this->bucketed.resize(this->p.input_size);

  BucketTestInfo pp = this->p;

  if ((this->p.first == 0) && (this->p.last == this->p.input_size)) {
	  std::copy(this->data.begin(), this->data.end(), this->bucketed.begin());
	  this->bcounts = mxx::bucketing(this->bucketed, [&pp](std::pair<size_t, size_t> const & x){ return x.first % pp.bucket_count; },
			  this->p.bucket_count);
  } else {
	  std::cout << "not doing whole range, using imxx bucketing." << std::endl;
			  imxx::local::bucketing_impl(this->data, [&pp](std::pair<size_t, size_t> const & x){ return x.first % pp.bucket_count; },
					  this->p.bucket_count, this->bcounts,   this->bucketed,
			                                   this->p.first, this->p.last);
  }
}

TEST_P(BucketTest, permute )
{
	this->unbucketed.clear();
  this->mapping.clear();
  this->bucketed.clear();

  BucketTestInfo pp = this->p;

  // allocate.
  this->mapping.reserve(this->p.input_size);


  imxx::local::assign_to_buckets(this->data, [&pp](std::pair<size_t, size_t> const & x){ return x.first % pp.bucket_count; },
		  this->p.bucket_count,
                                 this->bcounts, this->mapping, this->p.first, this->p.last);

  imxx::local::bucket_to_permutation(this->bcounts, this->mapping, this->p.first, this->p.last);


  if ((pp.bucket_count > 0) && (this->p.last <= this->p.input_size) && (this->p.first <= this->p.last) ) {

	  // allocate.
	  this->bucketed.resize(this->p.input_size);


	  imxx::local::permute(this->data.begin() + this->p.first, this->data.begin() + this->p.last,
						   this->mapping.begin() + this->p.first,
						   this->bucketed.begin() + this->p.first,
						   this->p.first);
  }
}

TEST_P(BucketTest, inplace_permute)
{
	this->unbucketed.clear();
	  this->mapping.clear();
	  this->bucketed.clear();
	  // allocate.
	  this->mapping.reserve(this->p.input_size);

	  BucketTestInfo pp = this->p;

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

TEST_P(BucketTest, unpermute)
{
	  this->mapping.clear();
	  this->bucketed.clear();
	  this->unbucketed.clear();

	  // allocate.
	  this->mapping.reserve(this->p.input_size);

	  BucketTestInfo pp = this->p;

	  imxx::local::assign_to_buckets(this->data, [&pp](std::pair<size_t, size_t> const & x){ return x.first % pp.bucket_count; },
			  this->p.bucket_count,
	                                 this->bcounts, this->mapping, this->p.first, this->p.last);

	  imxx::local::bucket_to_permutation(this->bcounts, this->mapping, this->p.first, this->p.last);



	  if ((pp.bucket_count > 0) && (this->p.last <= this->p.input_size) && (this->p.first <= this->p.last) ) {


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

TEST_P(BucketTest, inplace_unpermute)
{
	  this->mapping.clear();
	  this->bucketed.clear();
	  this->unbucketed.clear();

	  // allocate.
	  this->mapping.reserve(this->p.input_size);

	  BucketTestInfo pp = this->p;

	  imxx::local::assign_to_buckets(this->data, [&pp](std::pair<size_t, size_t> const & x){ return x.first % pp.bucket_count; },
			  this->p.bucket_count,
	                                 this->bcounts, this->mapping, this->p.first, this->p.last);

	  imxx::local::bucket_to_permutation(this->bcounts, this->mapping, this->p.first, this->p.last);

	  if ((pp.bucket_count > 0) && (this->p.last <= this->p.input_size) && (this->p.first <= this->p.last) ) {

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



INSTANTIATE_TEST_CASE_P(Bliss, BucketTest, ::testing::Values(
    // base cases
//    BucketTestInfo(0UL, 0UL,  0UL, 0UL),   //  0, boundary case
//    BucketTestInfo(0UL, 0UL,  0UL, 1UL),   //  1, boundary case
//    BucketTestInfo(0UL, 0UL,  1UL, 0UL),   //  2, boundary case
//    BucketTestInfo(0UL, 0UL,  1UL, 1UL),   //  3, boundary case
    BucketTestInfo(0UL, 1UL,  0UL, 0UL),   //  4, boundary case
    BucketTestInfo(0UL, 1UL,  0UL, 1UL),   //  5, boundary case
    BucketTestInfo(0UL, 1UL,  1UL, 0UL),   //  6, boundary case
    BucketTestInfo(0UL, 1UL,  1UL, 1UL),   //  7, boundary case
//    BucketTestInfo(8UL, 0UL,  0UL, 0UL),   //  8, boundary case
//    BucketTestInfo(8UL, 0UL,  0UL, 8UL),   //  9, boundary case
//    BucketTestInfo(8UL, 0UL,  8UL, 0UL),   // 10, boundary case
//    BucketTestInfo(8UL, 0UL,  8UL, 8UL),   // 11, boundary case
    BucketTestInfo(8UL, 1UL,  0UL, 0UL),   // 12, boundary case
    BucketTestInfo(8UL, 1UL,  0UL, 8UL),   // 13, boundary case
//    //BucketTestInfo(8UL, 1UL,  8UL, 0UL), // ??, boundary case
    BucketTestInfo(8UL, 1UL,  8UL, 8UL),   // 14, boundary case

    // data size = 255, full and partial
    BucketTestInfo((1UL <<  4),       1UL,         0UL, (1UL <<  4)),  // 15, full
    BucketTestInfo((1UL <<  4),       1UL,         0UL, (1UL <<  2)),  // 16, front partial
    BucketTestInfo((1UL <<  4),       1UL, (1UL <<  3), (1UL <<  4)),  // 17, back partial
    BucketTestInfo((1UL <<  4),       1UL, (1UL <<  2), (1UL <<  3)),  // 18, middle partial
    BucketTestInfo((1UL <<  4), 1UL <<  2,         0UL, (1UL <<  4)),  // 19, full
    BucketTestInfo((1UL <<  4), 1UL <<  2,         0UL, (1UL <<  2)),  // 20, front partial
    BucketTestInfo((1UL <<  4), 1UL <<  2, (1UL <<  3), (1UL <<  4)),  // 21, back partial
    BucketTestInfo((1UL <<  4), 1UL <<  2, (1UL <<  2), (1UL <<  3)),  // 22, middle partial

    BucketTestInfo((1UL <<  4), 1UL <<  4,         0UL, (1UL <<  4)),  // 23, full
    BucketTestInfo((1UL <<  4), 1UL <<  4, (1UL <<  2), (1UL <<  3)),  // 24, middle partial
    BucketTestInfo((1UL <<  4), 1UL <<  8,         0UL, (1UL <<  4)),  // 25, full
    BucketTestInfo((1UL <<  4), 1UL <<  8, (1UL <<  2), (1UL <<  3)),  // 26, middle partial

    BucketTestInfo((1UL <<  8), 1UL <<  4,         0UL, (1UL <<  8)),  // 27, full
    BucketTestInfo((1UL <<  8), 1UL <<  4, (1UL <<  2), (1UL <<  3)),  // 28, middle partial
    BucketTestInfo((1UL <<  8), 1UL <<  8,         0UL, (1UL <<  8)),  // 29, full
    BucketTestInfo((1UL <<  8), 1UL <<  8, (1UL <<  2), (1UL <<  3)),  // 30, middle partial
    BucketTestInfo((1UL <<  8), 1UL << 16,         0UL, (1UL <<  8)),  // 31, full
    BucketTestInfo((1UL <<  8), 1UL << 16, (1UL <<  2), (1UL <<  3)),  // 32, middle partial

    BucketTestInfo((1UL << 16), 1UL <<  8,         0UL, (1UL << 16)),  // 33, full
    BucketTestInfo((1UL << 16), 1UL <<  8, (1UL <<  2), (1UL <<  3)),  // 34, middle partial
    BucketTestInfo((1UL << 16), 1UL << 16,         0UL, (1UL << 16)),  // 35, full
    BucketTestInfo((1UL << 16), 1UL << 16, (1UL <<  2), (1UL <<  3)),  // 36, middle partial
    BucketTestInfo((1UL << 16), 1UL << 17,         0UL, (1UL << 16)),  // 37, full
    BucketTestInfo((1UL << 16), 1UL << 17, (1UL <<  2), (1UL <<  3)),   // 38, middle partial

	BucketTestInfo((1UL << 16), 960,         0UL, (1UL << 16))  // 33, full

// //    BucketTestInfo((1UL << 16), 1UL << 32,         0UL, (1UL << 16)),  // full
// //    BucketTestInfo((1UL << 16), 1UL << 32, (1UL <<  2), (1UL <<  3)),  // middle partial
// //    BucketTestInfo((1UL << 32), 1UL << 16,         0UL, (1UL << 32)),  // full
// //    BucketTestInfo((1UL << 32), 1UL << 16, (1UL <<  2), (1UL <<  3)),  // middle partial
// //    BucketTestInfo((1UL << 32), 1UL << 32,         0UL, (1UL << 32)),  // full
// //    BucketTestInfo((1UL << 32), 1UL << 32, (1UL <<  2), (1UL <<  3))   // middle partial
));










