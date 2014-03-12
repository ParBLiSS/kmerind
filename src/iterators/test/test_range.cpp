/**
 * range_test.cpp
 *
 *  Created on: Feb 18, 2014
 *      Author: tpan
 */

#include "iterators/range.hpp"

#include <dirent.h>

#include <gtest/gtest.h>
#include <cstdint>
#include <unistd.h>
#include <vector>
#include <limits>

#include "config.hpp"

using namespace bliss::iterator;

template<typename T>
class RangeTest : public ::testing::Test
{
  protected:
    size_t page_size;

    virtual void SetUp()
    {
      page_size = sysconf(_SC_PAGE_SIZE);
    }
};

TYPED_TEST_CASE_P(RangeTest);

TYPED_TEST_P(RangeTest, equal){
range<TypeParam> r(0, 100, 3);
range<TypeParam> r2(0, 100, 0);
EXPECT_TRUE(r == r2);

range<TypeParam> r3(10, 100, 3);
range<TypeParam> r4(10, 100, 0);
EXPECT_TRUE(r3 == r4);

if (std::is_signed<TypeParam>::value)
{
  range<TypeParam> r(-10, 100, 3);
  range<TypeParam> r2(-10, 100, 0);
  EXPECT_TRUE(r == r2);

  range<TypeParam> r3(-101, -100, 3);
  range<TypeParam> r4(-101, -100, 0);
  EXPECT_TRUE(r3 == r4);
}
}

TYPED_TEST_P(RangeTest, assignment){
range<TypeParam> r;
range<TypeParam> r2(10, 100, 0);
r = r2;
EXPECT_TRUE(r == r2);

if (std::is_signed<TypeParam>::value)
{
  range<TypeParam> r3(-10, 100, 3);
  r = r3;
  EXPECT_TRUE(r == r3);
}

}

TYPED_TEST_P(RangeTest, copyConstruct){
range<TypeParam> r2(10, 100, 0);
range<TypeParam> r(r2);
EXPECT_TRUE(r == r2);

if (std::is_signed<TypeParam>::value)
{
  range<TypeParam> r3(-10, 100, 3);
  range<TypeParam> r4(r3);
  EXPECT_TRUE(r3 == r4);
}

}

TYPED_TEST_P(RangeTest, partition){
range<TypeParam> r;
TypeParam e;

std::vector<TypeParam> starts =
{ std::numeric_limits<TypeParam>::min(), std::numeric_limits<TypeParam>::lowest(), 0, 1, 2, std::numeric_limits<TypeParam>::max()-2, (std::numeric_limits<TypeParam>::max() >> 1) + 1};

std::vector<TypeParam> lens =
{ 0, 1, 2};
std::vector<size_t> partitionCount =
{ 1, 2, 3, 4, 5, 6, 7, 8, std::numeric_limits<size_t>::max()};

for (auto start : starts)
{
  for (auto len : lens)
  {
    for (size_t i : partitionCount)
    {
      if (len < i)
      continue;

      auto rem = len % i;
      auto div = len / i;

  //      printf("%ld, %d\n", static_cast<size_t>(len), i);
      int block = 0;
      r = range<TypeParam>::block_partition(i, block, start, start+len);
      EXPECT_EQ(start, r.start);
      e = (rem == 0 ? (div) : (div + 1)) + start;
      EXPECT_EQ(e, r.end);

      block = (i-1)/2;
      r = range<TypeParam>::block_partition(i, block, start, start+len);
      e = (rem == 0 ? block * div :
          ( block >= rem ? block * div + rem : block * (div + 1))) + start;
      EXPECT_EQ(e, r.start);
      e = (rem == 0 ? (block + 1) * div :
          ( (block + 1) >= rem ? (block + 1) * div + rem : (block + 1) * (div + 1))) + start;
      EXPECT_EQ(e, r.end);

      block = i-1;
      r = range<TypeParam>::block_partition(i, block, start, start+len);
      e = (rem == 0 ? block * div :
          ( block >= rem ? block * div + rem : block * (div + 1))) + start;
      EXPECT_EQ(e, r.start);
      EXPECT_EQ(start+len, r.end);
    }
  }
}
}

TYPED_TEST_P(RangeTest, align){
range<TypeParam> r;

std::vector<TypeParam> starts =
{ 0, 1, (std::numeric_limits<TypeParam>::max() >> 1) + 1, std::numeric_limits<TypeParam>::max()-1};

if (std::is_signed<TypeParam>::value)
{
  starts.push_back(-1);
}

std::vector<size_t> pageSizes =
{ 1, 64};

for (auto s : starts)
{
  for (auto p : pageSizes)
  {

    if (s - std::numeric_limits<TypeParam>::lowest() < p)
    continue;
    //printf("1 %ld %ld\n", static_cast<int64_t>(s), static_cast<int64_t>(p));

    r = range<TypeParam>(s, s+1);
    r = r.align_to_page(p);
    EXPECT_TRUE(r.is_page_aligned(p));
  }
}

}

REGISTER_TYPED_TEST_CASE_P(RangeTest, equal, assignment, copyConstruct,
                           partition, align);

//////////////////////// DEATH TESTS - test class named BlahDeathTest so gtest will not run these in a threaded context.

// typedef RangeTest RangeDeathTest
template<typename T>
class RangeDeathTest : public ::testing::Test
{
  protected:
    size_t page_size;

    virtual void SetUp()
    {
      page_size = sysconf(_SC_PAGE_SIZE);
    }
};

TYPED_TEST_CASE_P(RangeDeathTest);

TYPED_TEST_P(RangeDeathTest, constructFails){

std::string err_regex = ".*range.hpp.* Assertion .* failed.*";

EXPECT_EXIT(range<TypeParam>(std::numeric_limits<TypeParam>::max(), 100, 0), ::testing::KilledBySignal(SIGABRT), err_regex);

if (std::is_signed<TypeParam>::value)
{
  EXPECT_EXIT(range<TypeParam>(0, -100, 0), ::testing::KilledBySignal(SIGABRT), err_regex);
  EXPECT_EXIT(range<TypeParam>(0, std::numeric_limits<TypeParam>::max() + 1, 0), ::testing::KilledBySignal(SIGABRT), err_regex);
}
}

TYPED_TEST_P(RangeDeathTest, partitionFails){
range<TypeParam> r;

std::string err_regex = ".*range.hpp.*block_partition.* Assertion .* failed.*";

// end is before start
std::vector<TypeParam> starts =
{ std::numeric_limits<TypeParam>::min()+1, std::numeric_limits<TypeParam>::lowest()+1, 1, 2, std::numeric_limits<TypeParam>::max(), (std::numeric_limits<TypeParam>::max() >> 1) + 1};

std::vector<size_t> partitionCount =
{ 1, 2, 3, 4, 5, 6, 7, 8, std::numeric_limits<size_t>::max()};
for (auto start : starts)
{
  for (size_t i : partitionCount)
  {
    //printf("%ld, %ld\n", static_cast<int64_t>(start), i);
    EXPECT_EXIT(range<TypeParam>::block_partition(i, 0,         start, start - 1), ::testing::KilledBySignal(SIGABRT), err_regex);
    EXPECT_EXIT(range<TypeParam>::block_partition(i, (i-1) / 2, start, start - 1), ::testing::KilledBySignal(SIGABRT), err_regex);
    EXPECT_EXIT(range<TypeParam>::block_partition(i, i-1,       start, start - 1), ::testing::KilledBySignal(SIGABRT), err_regex);
  }
}

// negative or zero parition sizes
starts =
{ std::numeric_limits<TypeParam>::min(), std::numeric_limits<TypeParam>::lowest(), 0, 1, 2, std::numeric_limits<TypeParam>::max()-1, (std::numeric_limits<TypeParam>::max() >> 1) };
partitionCount =
{ 0, std::numeric_limits<size_t>::lowest(), std::numeric_limits<size_t>::min()};
for (auto start : starts)
{
  for (size_t i : partitionCount)
  {
    //printf("%ld, %ld\n", static_cast<int64_t>(start), i);

    EXPECT_EXIT(range<TypeParam>::block_partition(i, 0,         start, start + 1), ::testing::KilledBySignal(SIGABRT), err_regex);
    EXPECT_EXIT(range<TypeParam>::block_partition(i, (i-1) / 2, start, start + 1), ::testing::KilledBySignal(SIGABRT), err_regex);
    EXPECT_EXIT(range<TypeParam>::block_partition(i, i-1,       start, start + 1), ::testing::KilledBySignal(SIGABRT), err_regex);
  }
}

// proc id is not right
starts =
{ std::numeric_limits<TypeParam>::min(), std::numeric_limits<TypeParam>::lowest(), 0, 1, 2, std::numeric_limits<TypeParam>::max()-1, (std::numeric_limits<TypeParam>::max() >> 1) };
partitionCount =
{ 1, 2, 3, 4, 5, 6, 7, 8, std::numeric_limits<size_t>::max()};
for (auto start : starts)
{
  for (size_t i : partitionCount)
  {
    //printf("%ld, %ld\n", static_cast<int64_t>(start), i);
    EXPECT_EXIT(range<TypeParam>::block_partition(i,i, start, start + 1), ::testing::KilledBySignal(SIGABRT), err_regex);
  }
}

}

TYPED_TEST_P(RangeDeathTest, alignFails){
range<TypeParam> r;

std::vector<TypeParam> starts =
{ 0, 1, (std::numeric_limits<TypeParam>::max() >> 1) + 1, std::numeric_limits<TypeParam>::max()-1};

if (std::is_signed<TypeParam>::value)
{
  starts.push_back(std::numeric_limits<TypeParam>::min());
  starts.push_back(std::numeric_limits<TypeParam>::lowest());
}
if (std::is_floating_point<TypeParam>::value)
{
  starts.push_back(std::numeric_limits<TypeParam>::lowest());
}

std::vector<size_t> pageSizes =
{ 0, std::numeric_limits<size_t>::lowest(), std::numeric_limits<size_t>::min()};

std::string err_regex = ".*range.hpp.*align_to_page.* Assertion .* failed.*";

for (auto s : starts)
{
  for (auto p : pageSizes)
  {
    //printf("align fail processing start %ld page size %lud\n", static_cast<int64_t>(s), static_cast<uint64_t>(p));

    r = range<TypeParam>(s, s+1);
    EXPECT_EXIT(r.align_to_page(p), ::testing::KilledBySignal(SIGABRT), err_regex);

  }

  // if s is negative, also check to make sure we fail with large p.
  if (s < 0) {
    r = range<TypeParam>(s, s+1);
    EXPECT_EXIT(r.align_to_page(std::numeric_limits<size_t>::max()), ::testing::KilledBySignal(SIGABRT), err_regex);
  }
}

}

REGISTER_TYPED_TEST_CASE_P(RangeDeathTest, constructFails, partitionFails,
                           alignFails);

//////////////////// RUN the tests with different types.

typedef ::testing::Types<char, uint8_t, int16_t, uint16_t, int, uint32_t,
    int64_t, uint64_t, size_t> RangeTestTypes;
INSTANTIATE_TYPED_TEST_CASE_P(Bliss, RangeTest, RangeTestTypes);
INSTANTIATE_TYPED_TEST_CASE_P(Bliss, RangeDeathTest, RangeTestTypes);
