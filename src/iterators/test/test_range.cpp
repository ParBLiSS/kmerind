/**
 * range_test.cpp
 *
 *  Created on: Feb 18, 2014
 *      Author: tpan
 */

#include "iterators/range.hpp"

#include <gtest/gtest.h>
#include <cstdint>
#include <unistd.h>
#include <vector>
#include <limits>

#include "config.hpp"

using namespace bliss::iterator;

template <typename T>
class RangeTest : public ::testing::Test {
  protected:
    size_t page_size;

    virtual void SetUp() {
      page_size = sysconf(_SC_PAGE_SIZE);
    }
};


TYPED_TEST_CASE_P(RangeTest);

TYPED_TEST_P(RangeTest, equal) {
  range<TypeParam> r(0, 100, 3, 1);
  range<TypeParam> r2(0, 100, 0, 1);
  ASSERT_TRUE(r == r2);
}


TYPED_TEST_P(RangeTest, assignment) {
  range<TypeParam> r;
  range<TypeParam> r2(0, 100, 0, 1);
  r = r2;
  ASSERT_TRUE(r == r2);
}

TYPED_TEST_P(RangeTest, copyConstruct) {
  range<TypeParam> r2(0, 100, 0, 1);
  range<TypeParam> r(r2);
  ASSERT_TRUE(r == r2);
}

TYPED_TEST_P(RangeTest, partition) {
  range<TypeParam> r;

  std::vector<TypeParam> lens = {1, std::numeric_limits<TypeParam>::lowest(), std::numeric_limits<TypeParam>::max(), 127};
  std::vector<int> partitionCount = {1, 2, 3, 4, 5, 7, 9, 11, 13};
  for (auto len : lens) {
    for (int i : partitionCount) {
      if (len < i)
        continue;

      r = range<TypeParam>::block_partition(len, i, 0);
      ASSERT_EQ(0, r.start);
      ASSERT_EQ((len/i + 1), r.end);

      r = range<TypeParam>::block_partition(len, i, i-1);
      ASSERT_EQ(len - len/i, r.start);
      ASSERT_EQ(len, r.end);
    }
  }
}

TYPED_TEST_P(RangeTest, align) {
  range<TypeParam> r;

  std::vector<TypeParam> starts = {0, 1, std::numeric_limits<TypeParam>::lowest(), std::numeric_limits<TypeParam>::max()/11, std::numeric_limits<TypeParam>::max()-1};
  std::vector<size_t> pageSizes = {1, 256, 4096};

  for (auto s : starts) {
    for (auto p : pageSizes) {
      r = range<TypeParam>(s, s+1);
      r = r.align_to_page(this->page_size);
      ASSERT_TRUE(r.is_page_aligned(this->page_size));

      r = range<TypeParam>::block_partition(len, i, i-1);
      r = r.align_to_page(this->page_size);
      ASSERT_TRUE(r.is_page_aligned(this->page_size));
    }
  }

}


TYPED_TEST_P(RangeTest, partitionFails) {
  range<TypeParam> r;

  std::vector<TypeParam> lens = {1, std::numeric_limits<TypeParam>::lowest(), std::numeric_limits<TypeParam>::max(), 127};
  std::vector<int> partitionCount = {1, 2, 3, 4, 5, 7, 9, 11, 13};
  for (auto len : lens) {
    for (int i : partitionCount) {
      if (len < i)
        continue;


      r = range<TypeParam>::block_partition(len, i, 0);
      ASSERT_EQ(0, r.start);
      ASSERT_EQ((len/i + 1), r.end);

      r = range<TypeParam>::block_partition(len, i, i-1);
      ASSERT_EQ(len - len/i, r.start);
      ASSERT_EQ(len, r.end);
    }
  }
}

TYPED_TEST_P(RangeTest, alignFails) {
  range<TypeParam> r;

  std::vector<TypeParam> lens = {1, std::numeric_limits<TypeParam>::lowest(), std::numeric_limits<TypeParam>::max(), 127};
  std::vector<int> partitionCount = {1, 2, 3, 4, 5, 7, 9, 11, 13};

  for (auto len : lens) {
    for (int i : partitionCount) {
      if (len < i)
        continue;

      r = range<TypeParam>::block_partition(len, i, 0);
      r = r.align_to_page(this->page_size);
      ASSERT_TRUE(r.is_page_aligned(this->page_size));

      r = range<TypeParam>::block_partition(len, i, i-1);
      r = r.align_to_page(this->page_size);
      ASSERT_TRUE(r.is_page_aligned(this->page_size));
    }
  }
}


REGISTER_TYPED_TEST_CASE_P(RangeTest, equal, assignment, copyConstruct, partition, align);

typedef ::testing::Types<char, uint8_t, int16_t, uint16_t, int, uint32_t, int64_t, uint64_t, size_t> RangeTestTypes;
INSTANTIATE_TYPED_TEST_CASE_P(Bliss, RangeTest, RangeTestTypes);
