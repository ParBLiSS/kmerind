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

struct A2ADistributeBenchmarkInfo {
    size_t input_size;
    size_t input_offset;
    size_t output_size;
    size_t output_offset;
    size_t block_size;

    A2ADistributeBenchmarkInfo() = default;
    A2ADistributeBenchmarkInfo(size_t const & _input_size, size_t const& _input_offset, size_t const & _output_size, size_t const & _output_offset, size_t const & _block_size) :
      input_size(_input_size), input_offset(_input_offset), output_size(_output_size), output_offset(_output_offset), block_size(_block_size) {};

    A2ADistributeBenchmarkInfo(A2ADistributeBenchmarkInfo const & other) = default;
    A2ADistributeBenchmarkInfo& operator=(A2ADistributeBenchmarkInfo const & other) = default;
    A2ADistributeBenchmarkInfo(A2ADistributeBenchmarkInfo && other) = default;
    A2ADistributeBenchmarkInfo& operator=(A2ADistributeBenchmarkInfo && other) = default;
};


/**
 * @brief test fixture for file parsers, which segments file into sequences.
 */
class A2ADistributeBenchmark : public ::testing::TestWithParam<A2ADistributeBenchmarkInfo>
{
  protected:
    A2ADistributeBenchmark() {};
    virtual ~A2ADistributeBenchmark() {};

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

    }

    virtual void TearDown() {

    }

    A2ADistributeBenchmarkInfo p;

    // input data
    std::vector<T> data;

    // one of 2 possible outputs.
    std::vector<T> distributed;
    std::vector<T> roundtripped;

};



TEST_P(A2ADistributeBenchmark, block_a2a)
{

  ::mxx::comm comm;

  this->init(comm);
  this->roundtripped.clear();

  // allocate.
  A2ADistributeBenchmarkInfo pp = this->p;

  this->distributed.clear();
  this->distributed.resize(pp.output_size);

  imxx::block_all2all(this->data, pp.block_size, this->distributed,
                      pp.input_offset, pp.output_offset, comm);

}


TEST_P(A2ADistributeBenchmark, block_a2a_inplace)
{
  ::mxx::comm comm;

  this->init(comm);
  this->roundtripped.clear();

  // allocate.
  A2ADistributeBenchmarkInfo pp = this->p;

  this->distributed.clear();
  this->distributed.resize(pp.output_size);

  if ((pp.input_size > 0) && (pp.output_size > 0))
    std::copy(this->data.begin() + pp.input_offset, this->data.begin() + pp.input_offset + pp.block_size * comm.size(), this->distributed.begin() + pp.output_offset);

  imxx::block_all2all_inplace(this->distributed, pp.block_size,
                      pp.output_offset, comm);
}



TEST_P(A2ADistributeBenchmark, block_roundtrip)
{
  ::mxx::comm comm;

  this->init(comm);

  // allocate.
  A2ADistributeBenchmarkInfo pp = this->p;

  this->distributed.clear();
  this->distributed.resize(pp.output_size);

  imxx::block_all2all(this->data, pp.block_size, this->distributed,
                      pp.input_offset, pp.output_offset, comm);


  this->roundtripped.clear();
  this->roundtripped.resize(pp.input_size);

  imxx::block_all2all(this->distributed, pp.block_size, this->roundtripped,
                      pp.output_offset, pp.input_offset, comm);

}


TEST_P(A2ADistributeBenchmark, block_roundtrip_inplace)
{
  ::mxx::comm comm;

  this->init(comm);

  // allocate.
  A2ADistributeBenchmarkInfo pp = this->p;

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



INSTANTIATE_TEST_CASE_P(Bliss, A2ADistributeBenchmark, ::testing::Values(
    // base cases
    A2ADistributeBenchmarkInfo((1UL << 24), 0UL, (1UL <<  24), 0UL, (1UL <<  8))  // 8


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


struct DistributeBenchmarkInfo {
    size_t input_size;

    DistributeBenchmarkInfo() = default;
    DistributeBenchmarkInfo(size_t const & _input_size) :
      input_size(_input_size) {};

    DistributeBenchmarkInfo(DistributeBenchmarkInfo const & other) = default;
    DistributeBenchmarkInfo& operator=(DistributeBenchmarkInfo const & other) = default;
    DistributeBenchmarkInfo(DistributeBenchmarkInfo && other) = default;
    DistributeBenchmarkInfo& operator=(DistributeBenchmarkInfo && other) = default;
};


/**
 * @brief test fixture for file parsers, which segments file into sequences.
 */
class DistributeBenchmark : public ::testing::TestWithParam<DistributeBenchmarkInfo>
{
  protected:
    DistributeBenchmark() {};
    virtual ~DistributeBenchmark() {};

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

    }

    virtual void TearDown() {
    }

    DistributeBenchmarkInfo p;

    // input data
    std::vector<T> data;

    // one of 2 possible outputs.
    std::vector<T> distributed;
    std::vector<T> roundtripped;

};



TEST_P(DistributeBenchmark, distribute_preserve_input)
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

TEST_P(DistributeBenchmark, distribute)
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


TEST_P(DistributeBenchmark, distribute_preserve_input_rt)
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

TEST_P(DistributeBenchmark, distribute_rt)
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

TEST_P(DistributeBenchmark, scatter_compute_gather)
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

}

INSTANTIATE_TEST_CASE_P(Bliss, DistributeBenchmark, ::testing::Values(

    DistributeBenchmarkInfo((1UL << 24))   // 3

));





//=========================  2 Part DISTRIBUTE TESTS

struct Distribute2PartBenchmarkInfo {
    size_t input_size;

    Distribute2PartBenchmarkInfo() = default;
    Distribute2PartBenchmarkInfo(size_t const & _input_size) :
      input_size(_input_size) {};

    Distribute2PartBenchmarkInfo(Distribute2PartBenchmarkInfo const & other) = default;
    Distribute2PartBenchmarkInfo& operator=(Distribute2PartBenchmarkInfo const & other) = default;
    Distribute2PartBenchmarkInfo(Distribute2PartBenchmarkInfo && other) = default;
    Distribute2PartBenchmarkInfo& operator=(Distribute2PartBenchmarkInfo && other) = default;
};


/**
 * @brief test fixture for file parsers, which segments file into sequences.
 */
class Distribute2PartBenchmark : public ::testing::TestWithParam<Distribute2PartBenchmarkInfo>
{
  protected:
    Distribute2PartBenchmark() {};
    virtual ~Distribute2PartBenchmark() {};

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


    }

    virtual void TearDown() {

    }

    Distribute2PartBenchmarkInfo p;

    // input data
    std::vector<T> data;

    // one of 2 possible outputs.
    std::vector<T> distributed;
    std::vector<T> roundtripped;

};



TEST_P(Distribute2PartBenchmark, distribute_preserve_input_2part)
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

TEST_P(Distribute2PartBenchmark, distribute_2part)
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


TEST_P(Distribute2PartBenchmark, distribute_preserve_input_2part_rt)
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

TEST_P(Distribute2PartBenchmark, distribute_2part_rt)
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


TEST_P(Distribute2PartBenchmark, scatter_compute_gather_2part)
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

}


TEST_P(Distribute2PartBenchmark, scatter_compute_gather_lowmem)
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

}


INSTANTIATE_TEST_CASE_P(Bliss, Distribute2PartBenchmark, ::testing::Values(

    Distribute2PartBenchmarkInfo((1UL << 24))   // 3

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


