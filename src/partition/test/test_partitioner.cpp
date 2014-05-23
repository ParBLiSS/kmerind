/**
 * range_test.cpp
 * Test range class
 *  Created on: Feb 18, 2014
 *      Author: tpan
 */

#include "partition/range.hpp"

#include <gtest/gtest.h>
#include <cstdint>  // for uint64_t, etc.
#include <vector>
#include <limits>

using namespace bliss::partition;

/*
 * test class holding some information.  Also, needed for the typed tests
 */
template<typename T>
class PartitionTest : public ::testing::Test
{
  protected:
    size_t page_size;

    virtual void SetUp()
    {
      page_size = sysconf(_SC_PAGE_SIZE);
    }
};

// indicate this is a typed test
TYPED_TEST_CASE_P(PartitionTest);




// test the block partitioning operation.  Only testing the base function that all other overloaded functions calls
TYPED_TEST_P(PartitionTest, partition){
  typedef bliss::partition::range<TypeParam> RangeType;
  typedef bliss::partition::BlockPartitioner<range<TypeParam> > PartitionerType;

  RangeType src, r;
  TypeParam e;


  std::vector<TypeParam> starts =
  { std::numeric_limits<TypeParam>::min(), std::numeric_limits<TypeParam>::lowest(),
    0, 1, 2,
    std::numeric_limits<TypeParam>::max()-2, (std::numeric_limits<TypeParam>::max() >> 1) + 1};

  std::vector<TypeParam> lens =
  { 0, 1, 2};

  std::vector<int> partitionCount =
  { 1, 2, std::numeric_limits<int>::max()};

  PartitionerType part;
  for (auto start : starts)
  {
    for (auto len : lens)
    {
      src = RangeType(start, start + len);

      for (auto p : partitionCount)
      {
//        if (len < i)
//          continue;
        part.configure(src, p, 1);



        auto rem = len % p;
        auto div = len / p;

        //      printf("%ld, %d\n", static_cast<size_t>(len), i);
        int block = 0;

        // first block
        r = part.getNext(block);
        EXPECT_EQ(start, r.start);
        e = (rem == 0 ? (div) : (div + 1)) + start;
        EXPECT_EQ(e, r.end);

        // middle block
        block = (p-1)/2;
        r = part.getNext(block);
        e = (rem == 0 ? block * div :
            ( block >= rem ? block * div + rem : block * (div + 1))) + start;
        EXPECT_EQ(e, r.start);
        e = (rem == 0 ? (block + 1) * div :
            ( (block + 1) >= rem ? (block + 1) * div + rem : (block + 1) * (div + 1))) + start;
        EXPECT_EQ(e, r.end);

        // last block
        block = p-1;
        r = part.getNext(block);
        e = (rem == 0 ? block * div :
            ( block >= rem ? block * div + rem : block * (div + 1))) + start;
        EXPECT_EQ(e, r.start);
        EXPECT_EQ(start+len, r.end);
      }
    }
  }
}


// now register the test cases
REGISTER_TYPED_TEST_CASE_P(PartitionTest, partition);

////////////////////////
//  DEATH TESTS - test class named BlahDeathTest so gtest will not run these in a threaded context. and will run first

// typedef PartitionTest PartitionDeathTest
template<typename T>
class PartitionDeathTest : public ::testing::Test
{
  protected:
    size_t page_size;

    virtual void SetUp()
    {
      page_size = sysconf(_SC_PAGE_SIZE);
    }
};

// annotate that PartitionDeathTest is typed
TYPED_TEST_CASE_P(PartitionDeathTest);


// failed partitions due to asserts.
TYPED_TEST_P(PartitionDeathTest, badNumPartitions){
  typedef bliss::partition::range<TypeParam> RangeType;
  typedef bliss::partition::BlockPartitioner<range<TypeParam> > PartitionerType;


  RangeType src, r;

  std::string err_regex = ".*partitioner.hpp.*Partitioner.* Assertion .* failed.*";

  // end is before start
  std::vector<TypeParam> starts =
  { std::numeric_limits<TypeParam>::min()+1, std::numeric_limits<TypeParam>::lowest()+1, 1, 2, std::numeric_limits<TypeParam>::max(), (std::numeric_limits<TypeParam>::max() >> 1) + 1};

  std::vector<int> partitionCount =
  { 0, std::numeric_limits<int>::lowest(), std::numeric_limits<int>::min()};

  // proc id is too big
  for (auto start : starts)
  {

    src = RangeType(start-1, start);

    for (auto i : partitionCount)
    {
      ;
      //printf("%ld, %ld\n", static_cast<int64_t>(start), i);
      EXPECT_EXIT(PartitionType(src, i, 1), ::testing::KilledBySignal(SIGABRT), err_regex);
    }
  }
}

// failed partitions due to asserts.
TYPED_TEST_P(PartitionDeathTest, badChunkSize){
  typedef bliss::partition::range<TypeParam> RangeType;
  typedef bliss::partition::BlockPartitioner<range<TypeParam> > PartitionerType;


  RangeType src, r;

  std::string err_regex = ".*partitioner.hpp.*Partitioner.* Assertion .* failed.*";

  // end is before start
  std::vector<TypeParam> starts =
  { std::numeric_limits<TypeParam>::min()+1, std::numeric_limits<TypeParam>::lowest()+1, 1, 2, std::numeric_limits<TypeParam>::max(), (std::numeric_limits<TypeParam>::max() >> 1) + 1};

  std::vector<size_t> chunkSizes =
  { 0, std::numeric_limits<size_t>::lowest(), std::numeric_limits<size_t>::min()};

  // proc id is too big
  for (auto start : starts)
  {

    src = RangeType(start-1, start);

    for (auto i : chunkSizes)
    {
      ;
      //printf("%ld, %ld\n", static_cast<int64_t>(start), i);
      EXPECT_EXIT(PartitionType part(src, 2, i), ::testing::KilledBySignal(SIGABRT), err_regex);
    }
  }

}


// failed partitions due to asserts.
TYPED_TEST_P(PartitionDeathTest, badPartitionId){
  typedef bliss::partition::range<TypeParam> RangeType;
  typedef bliss::partition::BlockPartitioner<range<TypeParam> > PartitionerType;


  RangeType src, r;

  std::string err_regex = ".*partitioner.hpp.*Partitioner.* Assertion .* failed.*";

  // end is before start
  std::vector<TypeParam> starts =
  { std::numeric_limits<TypeParam>::min()+1, std::numeric_limits<TypeParam>::lowest()+1, 1, 2, std::numeric_limits<TypeParam>::max(), (std::numeric_limits<TypeParam>::max() >> 1) + 1};


  // negative or zero parition sizes
  std::vector<int> partitionIds =
  { 0, std::numeric_limits<int>::lowest(), std::numeric_limits<int>::min(), std::numeric_limits<int>::max() };

  for (auto start : starts)
  {

    src = RangeType(start-1, start);

    for (auto i : partitionIds)
    {
      //printf("%ld, %ld\n", static_cast<int64_t>(start), i);
      PartitionType part(src, 4, 1);

      EXPECT_EXIT(part.getNext(i), ::testing::KilledBySignal(SIGABRT), err_regex);
    }
  }
}


// register the death test cases
REGISTER_TYPED_TEST_CASE_P(PartitionDeathTest, badNumPartitions, badPartitionId, badChunkSize);

//////////////////// RUN the tests with different types.

typedef ::testing::Types<int8_t, uint8_t, int16_t, uint16_t, int32_t, uint32_t,
    int64_t, uint64_t, size_t> PartitionTestTypes;
INSTANTIATE_TYPED_TEST_CASE_P(Bliss, PartitionTest, PartitionTestTypes);
INSTANTIATE_TYPED_TEST_CASE_P(Bliss, PartitionDeathTest, PartitionTestTypes);
