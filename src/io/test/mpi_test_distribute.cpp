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


//===============  BLOCK All2All tests

struct A2ADistributeTestInfo {
    size_t input_size;
    size_t input_offset;
    size_t output_size;
    size_t output_offset;
    size_t block_size;

    A2ADistributeTestInfo() = default;
    A2ADistributeTestInfo(size_t const & _input_size, size_t const& _input_offset, size_t const & _output_size, size_t const & _output_offset, size_t const & _block_size) :
      input_size(_input_size), input_offset(_input_offset), output_size(_output_size), output_offset(_output_offset), block_size(_block_size) {};

    A2ADistributeTestInfo(A2ADistributeTestInfo const & other) = default;
    A2ADistributeTestInfo& operator=(A2ADistributeTestInfo const & other) = default;
    A2ADistributeTestInfo(A2ADistributeTestInfo && other) = default;
    A2ADistributeTestInfo& operator=(A2ADistributeTestInfo && other) = default;
};


/**
 * @brief test fixture for file parsers, which segments file into sequences.
 */
class A2ADistributeTest : public ::testing::TestWithParam<A2ADistributeTestInfo>
{
  protected:
    A2ADistributeTest() {};
    virtual ~A2ADistributeTest() {};

    using T = std::pair<size_t, int>;

    virtual void SetUp()
    {
      p = GetParam();
      //std::cout << "p " << p.input_size << "," << p.bucket_count << "," << p.block_size << std::endl;

      distributed.clear();
      distributed.resize(p.output_size);
      roundtripped.clear();
      roundtripped.resize(p.input_size);


      // allocate test data
      data.clear();
      data.reserve(p.input_size);

      gold.clear();
      gold.resize(p.output_size);

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

      // move the data around using standard mpi calls.
      if ((p.input_size > 0) && (p.output_size > 0)) {
        ::mxx::datatype dt = ::mxx::get_datatype<T>();
        MPI_Alltoall(static_cast<T*>(&(data[p.input_offset])), p.block_size, dt.type(), static_cast<T*>(&(gold[p.output_offset])), p.block_size, dt.type(), comm);

      // check that the first nblocks all have same count in each bucket.
        bool same = true;
        size_t offset = p.output_offset;
        for (int i = 0; i < comm.size(); ++i) {
          for (size_t j = 0; j < p.block_size; ++j, ++offset) {
            same &= gold[offset].second == i;
          }
        }
        ASSERT_TRUE(same);
      }
    }

    virtual void TearDown() {

      // check distributed same as gold
      EXPECT_EQ(distributed.size(), gold.size());
      if (distributed.size() > 0) {
        bool same = true;
        for (size_t i = p.output_offset; i < p.output_offset + p.block_size; ++i) {
          same &= (distributed[i] == gold[i]);

          if (distributed[i] != gold[i]) {
            std::cout << " pos " << i << " distributed " << distributed[i].first << "," << distributed[i].second << "; gold " << gold[i].first << "," << gold[i].second << std::endl;
          }

        }
        EXPECT_TRUE(same);
      }

      // check roundtripped same as data
      if ((roundtripped.size() > 0) && (distributed.size() > 0)) {
        EXPECT_EQ(roundtripped.size(), data.size());
        bool same = true;
        for (size_t i = p.input_offset; i < p.input_offset + p.block_size; ++i) {
          same &= (roundtripped[i] == data[i]);
          if (roundtripped[i] != data[i]) {
            std::cout << " pos " << i << " roundtripped " << roundtripped[i].first << "," << roundtripped[i].second << "; data " << data[i].first << "," << data[i].second << std::endl;
          }

        }
        EXPECT_TRUE(same);
      }
    }

    A2ADistributeTestInfo p;

    // input data
    std::vector<T> data;

    // one of 2 possible outputs.
    std::vector<T> distributed;
    std::vector<T> roundtripped;

    // gold standard results
    std::vector<T> gold;
};



TEST_P(A2ADistributeTest, block_a2a)
{

  ::mxx::comm comm;

  this->init(comm);
  this->roundtripped.clear();

  // allocate.
  A2ADistributeTestInfo pp = this->p;

  this->distributed.clear();
  this->distributed.resize(pp.output_size);

  imxx::block_all2all(this->data, pp.block_size, this->distributed,
                      pp.input_offset, pp.output_offset, comm);

}


TEST_P(A2ADistributeTest, block_a2a_inplace)
{
  ::mxx::comm comm;

  this->init(comm);
  this->roundtripped.clear();

  // allocate.
  A2ADistributeTestInfo pp = this->p;

  this->distributed.clear();
  this->distributed.resize(pp.output_size);

  if ((pp.input_size > 0) && (pp.output_size > 0))
    std::copy(this->data.begin() + pp.input_offset, this->data.begin() + pp.input_offset + pp.block_size * comm.size(), this->distributed.begin() + pp.output_offset);

  imxx::block_all2all_inplace(this->distributed, pp.block_size,
                      pp.output_offset, comm);
}



TEST_P(A2ADistributeTest, block_roundtrip)
{
  ::mxx::comm comm;

  this->init(comm);

  // allocate.
  A2ADistributeTestInfo pp = this->p;

  this->distributed.clear();
  this->distributed.resize(pp.output_size);

  imxx::block_all2all(this->data, pp.block_size, this->distributed,
                      pp.input_offset, pp.output_offset, comm);


  this->roundtripped.clear();
  this->roundtripped.resize(pp.input_size);

  imxx::block_all2all(this->distributed, pp.block_size, this->roundtripped,
                      pp.output_offset, pp.input_offset, comm);

}


TEST_P(A2ADistributeTest, block_roundtrip_inplace)
{
  ::mxx::comm comm;

  this->init(comm);

  // allocate.
  A2ADistributeTestInfo pp = this->p;

  this->distributed.clear();
  this->distributed.resize(pp.output_size);

  if ((pp.input_size > 0) && (pp.output_size > 0))
    std::copy(this->data.begin() + pp.input_offset, this->data.begin() + pp.input_offset + pp.block_size * comm.size(), this->distributed.begin() + pp.output_offset);

  imxx::block_all2all_inplace(this->distributed, pp.block_size,
                      pp.output_offset, comm);


  this->roundtripped.clear();
  this->roundtripped.resize(pp.input_size);

  if ((pp.input_size > 0) && (pp.output_size > 0))
    std::copy(this->distributed.begin() + pp.output_offset, this->distributed.begin() + pp.output_offset + pp.block_size * comm.size(),
            this->roundtripped.begin() + pp.input_offset);

  imxx::block_all2all_inplace(this->roundtripped, pp.block_size,
                      pp.input_offset, comm);

}



INSTANTIATE_TEST_CASE_P(Bliss, A2ADistributeTest, ::testing::Values(
    // base cases
    A2ADistributeTestInfo(0UL, 0UL,  0UL, 0UL, 0UL),   //  0, boundary case
    A2ADistributeTestInfo(1UL, 0UL,  0UL, 0UL, 0UL),   //  1, boundary case
    A2ADistributeTestInfo(0UL, 0UL,  1UL, 0UL, 0UL),   //  2, boundary case

    A2ADistributeTestInfo((1UL << 16), 0UL,   (1UL <<  16), 0UL,   (1UL <<  8)),  // 3
    A2ADistributeTestInfo((1UL << 16), 256UL, (1UL <<  16), 0UL,   (1UL <<  8)),  // 4
    A2ADistributeTestInfo((1UL << 16), 0UL,   (1UL <<  16), 256UL, (1UL <<  8)),  // 5
    A2ADistributeTestInfo((1UL << 16), 255UL, (1UL <<  16), 255UL, (1UL <<  8)),  // 6

    A2ADistributeTestInfo((1UL << 15), 255UL, (1UL <<  16), 255UL, (1UL <<  8)),  // 7
    A2ADistributeTestInfo((1UL << 16), 255UL, (1UL <<  15), 255UL, (1UL <<  8)),  // 8


    A2ADistributeTestInfo(0UL, 0UL,  0UL, 0UL, 1UL)     //  9, boundary case
//    A2ADistributeTestInfo(1UL, 0UL,  0UL, 0UL, 1UL),  // 10, boundary case.  assert will fail
//    A2ADistributeTestInfo(0UL, 0UL,  1UL, 0UL, 1UL)   // 11, boundary case.  assert will fail

));


//=========================  DISTRIBUTE TESTS


template <typename IIT, typename OIT>
struct copy {
    void operator()(IIT begin, IIT end, OIT out) const {
      for (; begin != end; ++begin, ++out) {
        *out = *begin;
      }
    }
};


struct DistributeTestInfo {
    size_t input_size;

    DistributeTestInfo() = default;
    DistributeTestInfo(size_t const & _input_size) :
      input_size(_input_size) {};

    DistributeTestInfo(DistributeTestInfo const & other) = default;
    DistributeTestInfo& operator=(DistributeTestInfo const & other) = default;
    DistributeTestInfo(DistributeTestInfo && other) = default;
    DistributeTestInfo& operator=(DistributeTestInfo && other) = default;
};


/**
 * @brief test fixture for file parsers, which segments file into sequences.
 */
class DistributeTest : public ::testing::TestWithParam<DistributeTestInfo>
{
  protected:
    DistributeTest() {};
    virtual ~DistributeTest() {};

    using T = std::pair<size_t, int>;

    virtual void SetUp()
    {
      p = GetParam();
      //std::cout << "p " << p.input_size << "," << p.bucket_count << "," << p.block_size << std::endl;

      distributed.clear();
      distributed.reserve(p.input_size);
      roundtripped.clear();
      roundtripped.resize(p.input_size);

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
      std::vector<T> temp(data.begin(), data.end());

      // bucket it.
      int p = comm.size();
      std::vector<size_t> send_counts = ::mxx::bucketing(temp, [&p](T const & x ){ return x.first % p; }, p);

      // mxx all2all
      ::mxx::all2allv(temp, send_counts, comm).swap(gold);

    }

    virtual void TearDown() {

      // check distributed same as gold
      if (distributed.size() > 0) {
        EXPECT_EQ(distributed.size(), gold.size());
        bool same = true;
        for (size_t i = 0; i < distributed.size(); ++i) {
          same &= (distributed[i] == gold[i]);

          if (distributed[i] != gold[i]) {
            std::cout << " pos " << i << " distributed " << distributed[i].first << "," << distributed[i].second << "; gold " << gold[i].first << "," << gold[i].second << std::endl;
            ASSERT_TRUE(false);
          }

        }
        EXPECT_TRUE(same);
      }

      // check roundtripped same as data
      if (roundtripped.size() > 0) {
        EXPECT_EQ(roundtripped.size(), data.size());
        bool same = true;
        for (size_t i = 0; i < roundtripped.size(); ++i) {
          same &= (roundtripped[i] == data[i]);
          if (roundtripped[i] != data[i]) {
            std::cout << " pos " << i << " roundtripped " << roundtripped[i].first << "," << roundtripped[i].second << "; data " << data[i].first << "," << data[i].second << std::endl;
            ASSERT_TRUE(false);
          }

        }
        EXPECT_TRUE(same);
      }
    }

    DistributeTestInfo p;

    // input data
    std::vector<T> data;

    // one of 2 possible outputs.
    std::vector<T> distributed;
    std::vector<T> roundtripped;

    // gold standard results
    std::vector<T> gold;
};



TEST_P(DistributeTest, distribute_preserve_input)
{

  ::mxx::comm comm;

  this->init(comm);


  // copy data into roundtripped.
  this->roundtripped.resize(this->data.size());
  std::copy(this->data.begin(), this->data.end(), this->roundtripped.begin());

  // distribute
  int p = comm.size();
  std::vector<size_t> recv_counts;
  std::vector<size_t> mapping;

  imxx::distribute(this->roundtripped, [&p](T const & x ){ return x.first % p; },
                   recv_counts, mapping, this->distributed, comm, true);

}

TEST_P(DistributeTest, distribute)
{

  ::mxx::comm comm;

  this->init(comm);


  // copy data into roundtripped.
  this->roundtripped.resize(this->data.size());
  std::copy(this->data.begin(), this->data.end(), this->roundtripped.begin());

  // distribute
  int p = comm.size();
  std::vector<size_t> recv_counts;
  std::vector<size_t> mapping;

  imxx::distribute(this->roundtripped, [&p](T const & x ){ return x.first % p; },
                   recv_counts, mapping, this->distributed, comm, false);

  this->roundtripped.clear();
}



TEST_P(DistributeTest, distribute_preserve_input_rt)
{

  ::mxx::comm comm;

  this->init(comm);


  // copy data into roundtripped.
  this->roundtripped.resize(this->data.size());
  std::copy(this->data.begin(), this->data.end(), this->roundtripped.begin());

  // distribute
  int p = comm.size();
  std::vector<size_t> recv_counts;
  std::vector<size_t> mapping;

  imxx::distribute(this->roundtripped, [&p](T const & x ){ return x.first % p; },
                   recv_counts, mapping, this->distributed, comm, true);

  imxx::undistribute(distributed, recv_counts, mapping, this->roundtripped, comm, true);

}

TEST_P(DistributeTest, distribute_rt)
{

  ::mxx::comm comm;

  this->init(comm);


  // copy data into roundtripped.
  this->roundtripped.resize(this->data.size());
  std::copy(this->data.begin(), this->data.end(), this->roundtripped.begin());

  // distribute
  int p = comm.size();
  std::vector<size_t> recv_counts;
  std::vector<size_t> mapping;

  imxx::distribute(this->roundtripped, [&p](T const & x ){ return x.first % p; },
                   recv_counts, mapping, this->distributed, comm, false);

  imxx::undistribute(distributed, recv_counts, mapping, this->roundtripped, comm, true);
}

TEST_P(DistributeTest, scatter_compute_gather)
{

  ::mxx::comm comm;

  this->init(comm);


  // copy data into roundtripped.
  this->roundtripped.resize(this->data.size());

  this->distributed.resize(this->data.size());
  std::copy(this->data.begin(), this->data.end(), this->distributed.begin());


  // distribute
  int p = comm.size();
  std::vector<size_t> mapping;

  std::vector<T> inbuf;
  std::vector<T> outbuf;

  imxx::scatter_compute_gather(this->distributed, [&p](T const & x ){ return x.first % p; },
                               copy<typename std::vector<T>::const_iterator,
                                    typename std::vector<T>::iterator>(),
                   mapping, this->roundtripped, inbuf, outbuf, comm, false);

  this->distributed.clear();

  imxx::local::unpermute_inplace(this->roundtripped, mapping);

}




INSTANTIATE_TEST_CASE_P(Bliss, DistributeTest, ::testing::Values(
    // base cases
    DistributeTestInfo(0UL),   //  0, boundary case

    DistributeTestInfo((1UL <<  8)),  // 1
    DistributeTestInfo((1UL << 15)),  // 2
    DistributeTestInfo((1UL << 16))   // 3

));





//=========================  2 Part DISTRIBUTE TESTS

struct Distribute2PartTestInfo {
    size_t input_size;

    Distribute2PartTestInfo() = default;
    Distribute2PartTestInfo(size_t const & _input_size) :
      input_size(_input_size) {};

    Distribute2PartTestInfo(Distribute2PartTestInfo const & other) = default;
    Distribute2PartTestInfo& operator=(Distribute2PartTestInfo const & other) = default;
    Distribute2PartTestInfo(Distribute2PartTestInfo && other) = default;
    Distribute2PartTestInfo& operator=(Distribute2PartTestInfo && other) = default;
};


/**
 * @brief test fixture for file parsers, which segments file into sequences.
 */
class Distribute2PartTest : public ::testing::TestWithParam<Distribute2PartTestInfo>
{
  protected:
    Distribute2PartTest() {};
    virtual ~Distribute2PartTest() {};

    using T = std::pair<size_t, int>;

    virtual void SetUp()
    {
      p = GetParam();
      //std::cout << "p " << p.input_size << "," << p.bucket_count << "," << p.block_size << std::endl;

      distributed.clear();
      distributed.reserve(p.input_size);
      roundtripped.clear();
      roundtripped.resize(p.input_size);

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
      std::vector<T> temp(data.begin(), data.end());

      // bucket it.
      int p = comm.size();
      std::vector<size_t> send_counts = ::mxx::bucketing(temp, [&p](T const & x ){ return x.first % p; }, p);

      // get the minimum
      size_t min_bucket_size = *(std::min_element(send_counts.begin(), send_counts.end()));
      min_bucket_size = ::mxx::allreduce(min_bucket_size, mxx::min<size_t>(), comm);
      size_t first_part = min_bucket_size * comm.size();

      std::vector<size_t> offsets(send_counts.size() + 1, 0);
      for (size_t i = 1; i <= send_counts.size(); ++i) {
        offsets[i] = offsets[i-1] + send_counts[i-1];
        send_counts[i-1] -= min_bucket_size;
      }
      EXPECT_EQ(offsets.back(), temp.size());
      size_t second_part = std::accumulate(send_counts.begin(), send_counts.end(), static_cast<size_t>(0));
      EXPECT_EQ(second_part + first_part, temp.size());

      std::vector<size_t> recv_counts(send_counts.size(), 0);
      ::mxx::all2all(send_counts.data(), 1, recv_counts.data(), comm);
      second_part = std::accumulate(recv_counts.begin(), recv_counts.end(), static_cast<size_t>(0));

      gold.resize(first_part + second_part);

      // shuffle the buckets
      size_t offset1 = 0;
      size_t offset2 = 0;
      std::vector<T> temp2(temp.size());
      for (size_t i = 0; i < send_counts.size(); ++i) {
        std::copy(temp.begin() + offsets[i], temp.begin() + offsets[i] + min_bucket_size, temp2.begin() + offset1);
        offset1 += min_bucket_size;

        std::copy(temp.begin() + offsets[i] + min_bucket_size, temp.begin() + offsets[i+1], temp2.begin() + first_part + offset2);
        offset2 += send_counts[i];
      }

      // mxx all2all
      ::mxx::all2all(temp2.data(), min_bucket_size, gold.data(), comm);
      // mxx all2allv
      ::mxx::all2allv(temp2.data() + first_part, send_counts, gold.data() + first_part, recv_counts, comm);

    }

    virtual void TearDown() {

      // check distributed same as gold
      if (distributed.size() > 0) {
        EXPECT_EQ(distributed.size(), gold.size());
        bool same = true;
        for (size_t i = 0; i < distributed.size(); ++i) {
          same &= (distributed[i] == gold[i]);

          if (distributed[i] != gold[i]) {
            std::cout << " pos " << i << " distributed " << distributed[i].first << "," << distributed[i].second << "; gold " << gold[i].first << "," << gold[i].second << std::endl;
//            ASSERT_TRUE(false);
          }

        }
        EXPECT_TRUE(same);
      }

      // check roundtripped same as data
      if (roundtripped.size() > 0) {
        EXPECT_EQ(roundtripped.size(), data.size());
        bool same = true;
        for (size_t i = 0; i < roundtripped.size(); ++i) {
          same &= (roundtripped[i] == data[i]);
          if (roundtripped[i] != data[i]) {
            std::cout << " pos " << i << " roundtripped " << roundtripped[i].first << "," << roundtripped[i].second << "; data " << data[i].first << "," << data[i].second << std::endl;
            ASSERT_TRUE(false);
          }

        }
        EXPECT_TRUE(same);
      }
    }

    Distribute2PartTestInfo p;

    // input data
    std::vector<T> data;

    // one of 2 possible outputs.
    std::vector<T> distributed;
    std::vector<T> roundtripped;

    // gold standard results
    std::vector<T> gold;
};



TEST_P(Distribute2PartTest, distribute_preserve_input_2part)
{

  ::mxx::comm comm;

  this->init(comm);


  // copy data into roundtripped.
  this->roundtripped.resize(this->data.size());
  std::copy(this->data.begin(), this->data.end(), this->roundtripped.begin());

  // distribute
  int p = comm.size();
  std::vector<size_t> recv_counts;
  std::vector<size_t> mapping;

  imxx::distribute_2part(this->roundtripped, [&p](T const & x ){ return x.first % p; },
                   recv_counts, mapping, this->distributed, comm, true);

}

TEST_P(Distribute2PartTest, distribute_2part)
{

  ::mxx::comm comm;

  this->init(comm);


  // copy data into roundtripped.
  this->roundtripped.resize(this->data.size());
  std::copy(this->data.begin(), this->data.end(), this->roundtripped.begin());

  // distribute
  int p = comm.size();
  std::vector<size_t> recv_counts;
  std::vector<size_t> mapping;

  imxx::distribute_2part(this->roundtripped, [&p](T const & x ){ return x.first % p; },
                   recv_counts, mapping, this->distributed, comm, false);

  this->roundtripped.clear();
}



TEST_P(Distribute2PartTest, distribute_preserve_input_2part_rt)
{

  ::mxx::comm comm;

  this->init(comm);


  // copy data into roundtripped.
  this->roundtripped.resize(this->data.size());
  std::copy(this->data.begin(), this->data.end(), this->roundtripped.begin());

  // distribute
  int p = comm.size();
  std::vector<size_t> recv_counts;
  std::vector<size_t> mapping;

  imxx::distribute_2part(this->roundtripped, [&p](T const & x ){ return x.first % p; },
                   recv_counts, mapping, this->distributed, comm, true);

  imxx::undistribute_2part(this->distributed, recv_counts, mapping, this->roundtripped, comm, true);

}

TEST_P(Distribute2PartTest, distribute_2part_rt)
{

  ::mxx::comm comm;

  this->init(comm);


  // copy data into roundtripped.
  this->roundtripped.resize(this->data.size());
  std::copy(this->data.begin(), this->data.end(), this->roundtripped.begin());

  // distribute
  int p = comm.size();
  std::vector<size_t> recv_counts;
  std::vector<size_t> mapping;

  imxx::distribute_2part(this->roundtripped, [&p](T const & x ){ return x.first % p; },
                   recv_counts, mapping, this->distributed, comm, false);

  imxx::undistribute_2part(this->distributed, recv_counts, mapping, this->roundtripped, comm, true);
}


TEST_P(Distribute2PartTest, scatter_compute_gather_2part)
{

  ::mxx::comm comm;

  this->init(comm);


  // copy data into roundtripped.
  this->roundtripped.resize(this->data.size());

  this->distributed.resize(this->data.size());
  std::copy(this->data.begin(), this->data.end(), this->distributed.begin());


  // distribute
  int p = comm.size();
  std::vector<size_t> mapping;

  std::vector<T> inbuf;
  std::vector<T> outbuf;

  imxx::scatter_compute_gather_2part(this->distributed, [&p](T const & x ){ return x.first % p; },
                                     copy<typename std::vector<T>::const_iterator,
                                          typename std::vector<T>::iterator>(),
                   mapping, this->roundtripped, inbuf, outbuf, comm, false);

  this->distributed.clear();
  imxx::local::unpermute_inplace(this->roundtripped, mapping);

}

TEST_P(Distribute2PartTest, scatter_compute_gather_lowmem)
{

  ::mxx::comm comm;

  this->init(comm);


  // copy data into roundtripped.
  this->roundtripped.resize(this->data.size());

  this->distributed.resize(this->data.size());
  std::copy(this->data.begin(), this->data.end(), this->distributed.begin());


  // distribute
  int p = comm.size();
  std::vector<size_t> mapping;

  std::vector<T> inbuf;
  std::vector<T> outbuf;

  imxx::scatter_compute_gather_lowmem(this->distributed, [&p](T const & x ){ return x.first % p; },
                                     copy<typename std::vector<T>::const_iterator,
                                          typename std::vector<T>::iterator>(),
                   mapping, this->roundtripped, inbuf, outbuf, comm, false);

  this->distributed.clear();
  imxx::local::unpermute_inplace(this->roundtripped, mapping);

}



INSTANTIATE_TEST_CASE_P(Bliss, Distribute2PartTest, ::testing::Values(
    // base cases
    Distribute2PartTestInfo(0UL),   //  0, boundary case
    Distribute2PartTestInfo((1UL <<  2)),  // 1

    Distribute2PartTestInfo((1UL <<  8)),  // 1
    Distribute2PartTestInfo((1UL << 15)),  // 2
    Distribute2PartTestInfo((1UL << 16))   // 3

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


