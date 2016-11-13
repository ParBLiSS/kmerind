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

struct BlockBucketTestInfo {
    size_t input_size;
    size_t bucket_count;
    size_t block_size;

    BlockBucketTestInfo() = default;
    BlockBucketTestInfo(size_t const & _input_size, size_t const& _bucket_count, size_t const & _block_size) :
      input_size(_input_size), bucket_count(_bucket_count), block_size(_block_size) {};

    BlockBucketTestInfo(BlockBucketTestInfo const & other) = default;
    BlockBucketTestInfo& operator=(BlockBucketTestInfo const & other) = default;
    BlockBucketTestInfo(BlockBucketTestInfo && other) = default;
    BlockBucketTestInfo& operator=(BlockBucketTestInfo && other) = default;
};


/**
 * @brief test fixture for file parsers, which segments file into sequences.
 */
class BlockBucketTest : public ::testing::TestWithParam<BlockBucketTestInfo>
{
  protected:
    BlockBucketTest() {};
    virtual ~BlockBucketTest() {};

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

      gold_bcounts.clear();
      gold_bcounts.resize(p.bucket_count, 0);

      // populate with random numbers
      srand(17);
      size_t val;
      for (size_t i = 0; i < p.input_size; ++i) {
        val = rand();
        val <<= 32;
        val |= rand();
        data.emplace_back(val, i);

        if (p.bucket_count > 0) ++gold_bcounts[val % p.bucket_count];
      }

      if (p.bucket_count == 0) return;

      // get the number of blocks and the size of each bucket in the first blocks.
      size_t min_bucket = p.bucket_count > 0 ? *(std::min_element(gold_bcounts.begin(), gold_bcounts.end())) : 0;
      //size_t min_block = std::min(min_bucket, p.block_size);
      size_t nblocks = p.block_size > 0 ? (min_bucket / p.block_size) : 0;

      // compute the gold - use identity function
      gold.clear();
      gold.resize((nblocks + 1) * p.bucket_count);   // create bucket_count number of empty vectors.
      bcounts.clear();
      bcounts.resize(p.bucket_count);

      //std::cout << "min bucket " << min_bucket << " min block " << p.block_size << " nblocks " << nblocks << std::endl;

      if ((p.bucket_count > 0) && (p.input_size > 0)) {
        size_t pid;  // process (bucket) id
        size_t bid;  // block id.
        // first blocks
        for (size_t i = 0; i < p.input_size; ++i) {
        	//which bucket?
          pid = data[i].first % p.bucket_count;

          // which block does for the current bucket?
          bid = bcounts[pid]++;
          bid = p.block_size > 0 ? (bid / p.block_size) : 0;
          bid = (bid > nblocks) ? nblocks : bid;  // cap at nblocks

          gold[pid + bid * p.bucket_count].emplace_back(data[i]);
//          std::cout << " insert " << i << " src " << data[i].second << " bucket " << pid << " bid " << bid << " into " <<
//        		  (pid + bid * p.bucket_count) << " pos " << (gold[pid + bid * p.bucket_count].size() - 1) <<
//				  " bcounts " << (bcounts[pid] - 1) << std::endl;
        }

      }

      // check that the first nblocks all have same count in each bucket.
      bool same = true;
      for (size_t i = 0; i < nblocks * p.bucket_count; ++i) {
        same &= gold[i].size() == p.block_size;
      }
      EXPECT_TRUE(same);

      bcounts.clear();

    }

    virtual void TearDown() {

        // get the number of blocks and the size of each bucket in the first blocks.
        size_t min_bucket = p.bucket_count > 0 ? *(std::min_element(gold_bcounts.begin(), gold_bcounts.end())) : 0;
//        size_t min_block = std::min(min_bucket, p.block_size);
//        size_t nblocks = min_block > 0 ? (min_bucket / min_block) : 0;
        size_t nblocks = p.block_size > 0 ? (min_bucket / p.block_size) : 0;

      // check that the bucket counts (for the remainder block) are consistent with gold
      bool same = true;
      for (size_t i = 0; i < p.bucket_count; ++i) {
        same &= (bcounts[i] == (gold_bcounts[i] - nblocks * p.block_size));
      }
      EXPECT_TRUE(same);

      // compare mapping.
      same = true;
      size_t sid = 0;
      size_t bid = 0;   // bucket id
      if (mapping.size() > 0) {

        // the regular part.
        size_t n = 0;
        for (; n <= nblocks; ++n) {

          for (size_t i = 0; i < p.bucket_count; ++i) {
            bid = n * p.bucket_count + i;

            for (size_t j = 0; j < gold[bid].size(); ++j, ++sid) {
              //get the stored src id (in second field),  lookup its mapping, which should equal to the current position in output.
				  same &= (mapping[gold[bid][j].second] == sid) &&
					  ((gold[bid][j].first % p.bucket_count) == i);
            	if ((mapping[gold[bid][j].second] != sid) ||
                  ((gold[bid][j].first % p.bucket_count) != i)) {
            		std::cout << "block,bucket(linear),pos " << n << "," << i << "(" << bid << ")," << j <<
            				": val,src " << gold[bid][j].first << "," << gold[bid][j].second << " map_to " << mapping[gold[bid][j].second] <<
							" sid " << sid << " bucket_count " << p.bucket_count << " bucket_size " << p.block_size << std::endl << std::flush;
            		ASSERT_EQ(mapping[gold[bid][j].second], sid);
            		ASSERT_EQ((gold[bid][j].first % p.bucket_count), i);
            	}
            }
          }

        }

      }
      EXPECT_TRUE(same);

      same = true;
      sid = 0;
      if (bucketed.size() > 0) {

        for (size_t i = 0; i < gold.size(); ++i) {

          for (size_t j = 0; j < gold[i].size(); ++j, ++sid) {
//            if (gold[i][j].first != bucketed[sid].first) {
//              std::cout << "NO MATCH gold bucket " << i << " pos " << j << " " << gold[i][j].first << " from " << gold[i][j].second << "; bucketed id " << sid << " " << bucketed[sid].first << " from " << bucketed[sid].second << std::endl;
//            }

            same &= (gold[i][j] == bucketed[sid]);
          }
        }

      }
      EXPECT_TRUE(same);

      same = true;
      if (unbucketed.size() > 0) {

          for (size_t i = 0; i < p.input_size; ++i) {
            if (unbucketed[i].first != data[i].first) {
              std::cout << "NO MATCH data " << i << " " << data[i].first << " from " << data[i].second << "; unbucketed " << unbucketed[i].first << " from " << bucketed[i].second << std::endl;
            }
        	  same &= (unbucketed[i] == data[i]);
          }
      }
      EXPECT_TRUE(same);

    }

    BlockBucketTestInfo p;
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
    std::vector<size_t> gold_bcounts;
};



TEST_P(BlockBucketTest, assign)
{
	this->unbucketed.clear();
  this->bucketed.clear();

  // allocate.
  this->mapping.reserve(this->p.input_size);

  BlockBucketTestInfo pp = this->p;

  imxx::local::assign_to_buckets(this->data, [&pp](std::pair<size_t, size_t> const & x){ return x.first % pp.bucket_count; },
		  	  	  	  	  	  	  this->p.bucket_count,
                                 this->bcounts, this->mapping);

  imxx::local::bucket_to_block_permutation(this->p.block_size, this->bcounts, this->mapping);

}


TEST_P(BlockBucketTest, permute)
{
  this->mapping.clear();
  this->bucketed.clear();
  this->unbucketed.clear();

  BlockBucketTestInfo pp = this->p;

  // allocate.
  this->mapping.reserve(this->p.input_size);


  imxx::local::assign_to_buckets(this->data, [&pp](std::pair<size_t, size_t> const & x){ return x.first % pp.bucket_count; },
		  	  	  	  	  	  	  this->p.bucket_count,
                                 this->bcounts, this->mapping);

  imxx::local::bucket_to_block_permutation(this->p.block_size, this->bcounts, this->mapping);


  if (pp.bucket_count > 0) {

	  // allocate.
	  this->bucketed.resize(this->p.input_size);

	  imxx::local::permute(this->data,
						   this->mapping,
							  this->bucketed);

  }
}

TEST_P(BlockBucketTest, inplace_permute)
{
	this->unbucketed.clear();
	  this->mapping.clear();
	  this->bucketed.clear();
	  // allocate.
	  this->mapping.reserve(this->p.input_size);

	  BlockBucketTestInfo pp = this->p;

	  imxx::local::assign_to_buckets(this->data, [&pp](std::pair<size_t, size_t> const & x){ return x.first % pp.bucket_count; },
			  this->p.bucket_count,
	                                 this->bcounts, this->mapping);

	  imxx::local::bucket_to_block_permutation(this->p.block_size, this->bcounts, this->mapping);

	  if (pp.bucket_count > 0) {


		  // allocate.
		  this->bucketed.resize(this->p.input_size);

		  // copy
		  std::copy(this->data.begin(), this->data.end(), this->bucketed.begin());

		  imxx::local::permute_inplace(this->bucketed,
							   this->mapping);
	  }

}

TEST_P(BlockBucketTest, unpermute)
{
	  this->mapping.clear();
	  this->bucketed.clear();
	  this->unbucketed.clear();

	  // allocate.
	  this->mapping.reserve(this->p.input_size);

	  BlockBucketTestInfo pp = this->p;

	  imxx::local::assign_to_buckets(this->data, [&pp](std::pair<size_t, size_t> const & x){ return x.first % pp.bucket_count; },
			  this->p.bucket_count,
	                                 this->bcounts, this->mapping);

	  imxx::local::bucket_to_block_permutation(this->p.block_size, this->bcounts, this->mapping);

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

TEST_P(BlockBucketTest, inplace_unpermute)
{
	  this->mapping.clear();
	  this->bucketed.clear();
	  this->unbucketed.clear();

	  // allocate.
	  this->mapping.reserve(this->p.input_size);

	  BlockBucketTestInfo pp = this->p;

	  imxx::local::assign_to_buckets(this->data, [&pp](std::pair<size_t, size_t> const & x){ return x.first % pp.bucket_count; },
			  this->p.bucket_count,
	                                 this->bcounts, this->mapping);

	  imxx::local::bucket_to_block_permutation(this->p.block_size, this->bcounts, this->mapping);

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



INSTANTIATE_TEST_CASE_P(Bliss, BlockBucketTest, ::testing::Values(
    // base cases
//    BlockBucketTestInfo(0UL, 0UL,  0UL),   //  0, boundary case
//    BlockBucketTestInfo(0UL, 0UL,  1UL),   //  1, boundary case
//    BlockBucketTestInfo(0UL, 1UL,  0UL),   //  2, boundary case
    BlockBucketTestInfo(0UL, 1UL,  1UL),   //  3, boundary case
//    BlockBucketTestInfo(8UL, 0UL,  0UL),   //  4, boundary case
//    BlockBucketTestInfo(8UL, 0UL,  1UL),   //  5, boundary case
//    BlockBucketTestInfo(8UL, 0UL,  8UL),   //  6, boundary case
//    BlockBucketTestInfo(8UL, 1UL,  0UL),   //  7, boundary case
    BlockBucketTestInfo(8UL, 1UL,  1UL),   //  8, boundary case
    BlockBucketTestInfo(8UL, 1UL,  8UL),   //  9, boundary case

    // data size = 255, full and partial
    BlockBucketTestInfo((1UL <<  4),       1UL, (1UL <<  2)),  // 10, full
    BlockBucketTestInfo((1UL <<  4), 1UL <<  2, (1UL <<  2)),  // 11, full
    BlockBucketTestInfo((1UL <<  4), 1UL <<  4, (1UL <<  2)),  // 12, full
    BlockBucketTestInfo((1UL <<  4), 1UL <<  8, (1UL <<  2)),  // 13, full

    BlockBucketTestInfo((1UL <<  8), 1UL <<  4, (1UL <<  4)),  // 14, full
    BlockBucketTestInfo((1UL <<  8), 1UL <<  8, (1UL <<  4)),  // 15, full
    BlockBucketTestInfo((1UL <<  8), 1UL << 16, (1UL <<  4)),  // 16, full

    BlockBucketTestInfo((1UL << 16), 1UL <<  8, (1UL <<  8)),  // 17, full
    BlockBucketTestInfo((1UL << 16), 1UL << 16, (1UL <<  8)),  // 18, full
    BlockBucketTestInfo((1UL << 16), 1UL << 17, (1UL <<  8)),  // 19, full
    BlockBucketTestInfo((1UL << 16), 1UL << 14, (1UL <<  8)),  // 20, front partial
    BlockBucketTestInfo((1UL << 16), 1UL << 14, (1UL <<  8)),  // 21, back partial
    BlockBucketTestInfo((1UL << 16), 1UL << 14, (1UL << 12)),  // 22, middle partial
    BlockBucketTestInfo((1UL << 16), 1UL << 14, (1UL << 14)),  // 23, middle partial
    BlockBucketTestInfo((1UL << 16), 1UL << 14, (1UL << 15))   // 24, middle partial

// //    BlockBucketTestInfo((1UL << 16), 1UL << 32,         0UL, (1UL << 16)),  // full
// //    BlockBucketTestInfo((1UL << 16), 1UL << 32, (1UL <<  2), (1UL <<  3)),  // middle partial
// //    BlockBucketTestInfo((1UL << 32), 1UL << 16,         0UL, (1UL << 32)),  // full
// //    BlockBucketTestInfo((1UL << 32), 1UL << 16, (1UL <<  2), (1UL <<  3)),  // middle partial
// //    BlockBucketTestInfo((1UL << 32), 1UL << 32,         0UL, (1UL << 32)),  // full
// //    BlockBucketTestInfo((1UL << 32), 1UL << 32, (1UL <<  2), (1UL <<  3))   // middle partial
));










