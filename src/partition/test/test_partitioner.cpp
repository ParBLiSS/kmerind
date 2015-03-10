/**
 * range_test.cpp
 * Test range class
 *  Created on: Feb 18, 2014
 *      Author: Tony Pan <tpan7@gatech.edu>
 */

#include <gtest/gtest.h>
#include <cstdint>  // for uint64_t, etc.
#include <vector>
#include <limits>

#include "partition/range.hpp"
#include "partition/partitioner.hpp"



using namespace bliss::partition;

/*
 * test class holding some information.  Also, needed for the typed tests
 */
template<typename T>
class PartitionTest : public ::testing::Test
{
  protected:
    long int page_size;

    virtual void SetUp()
    {
      page_size = sysconf(_SC_PAGE_SIZE);
    }
};

// indicate this is a typed test
TYPED_TEST_CASE_P(PartitionTest);




// test the block partitioning operation.  Only testing the base function that all other overloaded functions calls
TYPED_TEST_P(PartitionTest, blockPartition){
  typedef bliss::partition::range<TypeParam> RangeType;
  typedef bliss::partition::BlockPartitioner<range<TypeParam> > PartitionerType;

  typedef decltype(std::declval<RangeType>().size()) SizeType;

  RangeType src, r;
  TypeParam e;


  std::vector<TypeParam> starts =
  { std::numeric_limits<TypeParam>::min(), std::numeric_limits<TypeParam>::lowest(),
    0, 1, 2,
    std::numeric_limits<TypeParam>::max()-2, (std::numeric_limits<TypeParam>::max() / 2) + 1};

  std::vector<SizeType> lens =
  { 0, 1, 2};

  std::vector<size_t> partitionCount =
  { 1, 2, std::numeric_limits<size_t>::max()};

  PartitionerType part;
  for (auto start : starts)
  {
    for (auto len : lens)
    {
      src = RangeType(start, (std::numeric_limits<TypeParam>::max() - start >= len) ? start + len : std::numeric_limits<TypeParam>::max());

      for (auto p : partitionCount)
      {
//        if (len < i)
//          continue;
        part.configure(src, p);


        auto div = len / static_cast<SizeType>(p);
        auto rem = len - div * static_cast<SizeType>(p);

        //      INFOF("%ld, %d\n", static_cast<size_t>(len), i);
        size_t block = 0;

        // first block
        r = part.getNext(block);
//        INFOF("here 1 p = %lu/%lu ", block, p);
//        INFO( "src: " << src << " part: " << r );
        EXPECT_EQ(start, r.start);
        e = (rem == 0 ? (div) : (div + 1)) + start;
        EXPECT_EQ(e, r.end);

        part.reset();

        // middle block
        block = (p-1)/2;
        r = part.getNext(block);
//        INFOF("here 2 p = %lu/%lu ", block, p);
//        INFO( "src: " << src << " part: " << r );
        e = (rem == 0 ? block * div :
            ( block >= rem ? block * div + rem : block * (div + 1))) + start;
        EXPECT_EQ(e, r.start);
        e = (rem == 0 ? (block + 1) * div :
            ( (block + 1) >= rem ? (block + 1) * div + rem : (block + 1) * (div + 1))) + start;
        EXPECT_EQ(e, r.end);

        part.reset();

        // last block
        block = p-1;
        r = part.getNext(block);
//        INFOF("here 3 p = %lu/%lu ", block, p);
//        INFO( "src: " << src << " part: " << r );
        e = (rem == 0 ? block * div :
            ( block >= rem ? block * div + rem : block * (div + 1))) + start;
        EXPECT_EQ(e, r.start);
        EXPECT_EQ(src.end, r.end);

      }
    }
  }
}


TYPED_TEST_P(PartitionTest, cyclicPartition){
  typedef bliss::partition::range<TypeParam> RangeType;
  typedef bliss::partition::CyclicPartitioner<range<TypeParam> > PartitionerType;
  typedef decltype(std::declval<RangeType>().size()) SizeType;

  RangeType src, r;
  TypeParam e;


  std::vector<TypeParam> starts =
  { std::numeric_limits<TypeParam>::min(), std::numeric_limits<TypeParam>::lowest(),
    0, 1, 2,
    std::numeric_limits<TypeParam>::max()-2, (std::numeric_limits<TypeParam>::max() / 2) + 1};

  std::vector<SizeType> lens =
  { 0, 1, 2, 8 };

  std::vector<size_t> partitionCount =
  { 1, 2, 4, std::numeric_limits<size_t>::max()};

  PartitionerType part;
  for (auto start : starts)
  {
//    INFO( "start: " << start );

    for (auto len : lens)
    {
//      INFO( "len: " << len );

      src = RangeType(start, (std::numeric_limits<TypeParam>::max() - start >= len) ? start + len : std::numeric_limits<TypeParam>::max());

      for (auto p : partitionCount)
      {
//        INFO( "parts: " << p );


        //        if (len < i)
//          continue;
        part.configure(src, p, 1);



        SizeType div = 1;

        //      INFOF("%ld, %d\n", static_cast<size_t>(len), i);
        size_t nChunks = static_cast<size_t>(std::ceil(src.size() / div));

        // middle block
        size_t block = (p-1)/2;
        r = part.getNext(block);
//        INFOF("here 1 p = %lu/%lu ", block, p);
//        INFO( "src: " << src << " part: " << r );
        e = std::min(static_cast<SizeType>(src.end),
                     static_cast<SizeType>(std::min(nChunks, block) * div + start));
        EXPECT_EQ(e, r.start);
        e = std::min(static_cast<SizeType>(src.end),
                     static_cast<SizeType>(std::min(nChunks, block + 1) * div + start));
        EXPECT_EQ(e, r.end);



        // middle block
        r = part.getNext(block);
//        INFOF("here 2 p = %lu/%lu ", block, p);
//        INFO( "src: " << src << " part: " << r );
        e = std::min(static_cast<SizeType>(src.end),
                     static_cast<SizeType>(std::min(nChunks, block + p) * div + start));
        EXPECT_EQ(e, r.start);
//        INFOF("block %lu, partition %lu, start %lu, len %lu, nchunks %lu, div %lu \n", block, p, start, len, nChunks, div);
        e = std::min(static_cast<SizeType>(src.end),
                     static_cast<SizeType>(std::min(nChunks, block + p + 1) * div + start));
        EXPECT_EQ(e, r.end);

        part.reset();

      }
    }
  }
}


TYPED_TEST_P(PartitionTest, demandPartition){
  typedef bliss::partition::range<TypeParam> RangeType;
  typedef bliss::partition::DemandDrivenPartitioner<range<TypeParam> > PartitionerType;
  typedef decltype(std::declval<RangeType>().size()) SizeType;

  RangeType src, r;
  TypeParam e;


  std::vector<TypeParam> starts =
  { std::numeric_limits<TypeParam>::min(), std::numeric_limits<TypeParam>::lowest(),
    0, 1, 2,
    std::numeric_limits<TypeParam>::max()-2, (std::numeric_limits<TypeParam>::max() / 2) + 1};

  std::vector<SizeType> lens =
  { 0, 1, 2, 8};

  std::vector<size_t> partitionCount =
  { 1, 2, 4, std::numeric_limits<size_t>::max()};

  PartitionerType part;
  for (auto start : starts)
  {
    for (auto len : lens)
    {
      src = RangeType(start, (std::numeric_limits<TypeParam>::max() - start >= static_cast<TypeParam>(len)) ? start + static_cast<TypeParam>(len) : std::numeric_limits<TypeParam>::max());

      for (auto p : partitionCount)
      {
//        if (len < i)
//          continue;
        part.configure(src, p, 1);


        SizeType div = 1;

        //      INFOF("%ld, %d\n", static_cast<size_t>(len), i);
        size_t block = 0;

        // first block
        r = part.getNext(block);
//        INFOF("here 1 block/p = %lu/%lu, len %lu, div %lu, start %d\n", block, p, len, div, start);
//        INFO( "src: " << src << " part: " << r );
        EXPECT_EQ(start, r.start);
        e = (start >= src.end) ? src.end :
            std::min(static_cast<SizeType>(src.end),
                     static_cast<SizeType>(div + start));
        EXPECT_EQ(e, r.end);

        // middle block
        block = (p-1)/2;
        r = part.getNext(block);
//        INFOF("here 1 block/p = %lu/%lu, len %lu, div %lu, start %d\n", block, p, len, div, start);
//        INFO( "src: " << src << " part: " << r );
        EXPECT_EQ(e, r.start);
        e = (e >= src.end) ? src.end :
            std::min(static_cast<SizeType>(src.end),
                     static_cast<SizeType>(div + e));
        EXPECT_EQ(e, r.end);

        // last block
        block = p-1;
        r = part.getNext(block);
//        INFOF("here 1 block/p = %lu/%lu, len %lu, div %lu, start %d\n", block, p, len, div, start);
//        INFO( "src: " << src << " part: " << r );
        EXPECT_EQ(e, r.start);
        e = (e >= src.end) ? src.end :
            std::min(static_cast<SizeType>(src.end),
                     static_cast<SizeType>(div + e));
        EXPECT_EQ(e, r.end);

      }
    }
  }
}


// failed partitions due to asserts.
TYPED_TEST_P(PartitionTest, badPartitionId){
  typedef bliss::partition::range<TypeParam> RangeType;
  typedef bliss::partition::BlockPartitioner<range<TypeParam> > PartitionerType;


  RangeType src, r;

  // end is before start
  std::vector<TypeParam> starts =
  { std::numeric_limits<TypeParam>::min()+1, std::numeric_limits<TypeParam>::lowest()+1, 1, 2, std::numeric_limits<TypeParam>::max(), (std::numeric_limits<TypeParam>::max() / 2) + 1};
  size_t size = 1;

  // negative or zero parition sizes
  std::vector<size_t> partitionIds =
  { std::numeric_limits<size_t>::max() };

  for (auto start : starts)
  {

    src = RangeType(start-size, start);

    for (auto i : partitionIds)
    {
      //INFOF("%ld, %ld\n", static_cast<int64_t>(start), i);
      PartitionerType part;
      part.configure(src, 4);

      try {
        part.getNext(i);
        ADD_FAILURE();

      } catch (const std::invalid_argument &e) {
        // okay, correct exception thrown.
      }
    }
  }
}




// now register the test cases
REGISTER_TYPED_TEST_CASE_P(PartitionTest, badPartitionId, blockPartition, cyclicPartition, demandPartition);




//////////////////// RUN the tests with different types.
typedef ::testing::Types<int8_t, uint8_t, int16_t, uint16_t, int32_t, uint32_t,
    int64_t, size_t, float, double> PartitionTestTypes;
//typedef ::testing::Types<size_t> PartitionTestTypes;
INSTANTIATE_TYPED_TEST_CASE_P(Bliss, PartitionTest, PartitionTestTypes);
