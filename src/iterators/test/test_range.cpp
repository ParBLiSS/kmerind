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
  EXPECT_TRUE(r == r2);
}


TYPED_TEST_P(RangeTest, assignment) {
  range<TypeParam> r;
  range<TypeParam> r2(0, 100, 0, 1);
  r = r2;
  EXPECT_TRUE(r == r2);
}

TYPED_TEST_P(RangeTest, copyConstruct) {
  range<TypeParam> r2(0, 100, 0, 1);
  range<TypeParam> r(r2);
  EXPECT_TRUE(r == r2);
}

TYPED_TEST_P(RangeTest, partition) {
  range<TypeParam> r;
  TypeParam e;

  std::vector<TypeParam> lens = {1, 2, 4, 8, 16, 32, 64, 127, std::numeric_limits<TypeParam>::max()};
  std::vector<int> partitionCount = {1, 2, 4, 8, 16, 32, 64, 127, std::numeric_limits<int>::max()};
  for (auto len : lens) {
    for (int i : partitionCount) {
      if (len < i)
        continue;

      //printf("%ld, %d\n", static_cast<size_t>(len), i);
      r = range<TypeParam>::block_partition(len, i, 0);
      EXPECT_EQ(0, r.start);
      e = ((len % i) == 0 ? (len / i) : (len / i + 1));
      EXPECT_EQ(e, r.end);

      r = range<TypeParam>::block_partition(len, i, (i-1) / 2);
      e = ((len % i) == 0 ? (i-1) / 2 * len / i :
              ( (i-1)/2 > (len % i) ? (i -1) / 2 * len / i + len % i : (i - 1) / 2 * (len / i + 1)));
      EXPECT_EQ(e, r.start);
      e = ((len % i) == 0 ? (i+1) / 2 * len / i :
              ( (i-1)/2 > (len % i) ? (i + 1) / 2 * len / i + len % i : (i + 1) / 2 * (len / i + 1)));
      EXPECT_EQ(e, r.end);

      r = range<TypeParam>::block_partition(len, i, i-1);
      e = (len % i == 0 ? (i-1) * len / i :
          ( i-1 > (len % i) ? (i - 1) * len / i + len % i : (i - 1) * (len / i + 1)));
      EXPECT_EQ(e, r.start);
      EXPECT_EQ(len, r.end);
    }
  }
}

TYPED_TEST_P(RangeTest, align) {
  range<TypeParam> r;

  std::vector<TypeParam> starts = {0, 1, std::numeric_limits<TypeParam>::lowest(), std::numeric_limits<TypeParam>::min(),
                                   std::numeric_limits<TypeParam>::max()/3, std::numeric_limits<TypeParam>::max()/2, std::numeric_limits<TypeParam>::max()-1};
  std::vector<TypeParam> pageSizes = {1, std::numeric_limits<TypeParam>::max() / 2, std::numeric_limits<TypeParam>::max()};

  for (auto s : starts) {
    for (auto p : pageSizes) {

      r = range<TypeParam>(s, s+1);
      r = r.align_to_page(p);
      EXPECT_TRUE(r.is_page_aligned(p));
    }
  }

}


TYPED_TEST_P(RangeTest, partitionFails) {
  range<TypeParam> r;

  std::vector<TypeParam> lens = {std::numeric_limits<TypeParam>::lowest(), std::numeric_limits<TypeParam>::min(), 0};
  std::vector<int> partitionCount = {1, 2, 4, 8, 16, 32, 64, 127, std::numeric_limits<int>::max()};
  for (auto len : lens) {
    for (int i : partitionCount) {
      EXPECT_EXIT(range<TypeParam>::block_partition(len, i, 0), ::testing::KilledBySignal(SIGABRT), "error on line .* of block_partition()");
      EXPECT_EXIT(range<TypeParam>::block_partition(len, i, (i-1) / 2), ::testing::KilledBySignal(SIGABRT), "error on line .* of block_partition()");
      EXPECT_EXIT(range<TypeParam>::block_partition(len, i, i-1), ::testing::KilledBySignal(SIGABRT), "error on line .* of block_partition()");
    }
  }

  lens = {1, 2, 4, 8, 16, 32, 64, 127, std::numeric_limits<TypeParam>::max()};
  partitionCount = {0, std::numeric_limits<int>::lowest(), std::numeric_limits<int>::max(), std::numeric_limits<int>::min()};
  for (auto len : lens) {
    for (int i : partitionCount) {
      EXPECT_EXIT(range<TypeParam>::block_partition(len, i, 0), ::testing::KilledBySignal(SIGABRT), "error on line .* of block_partition()");
      EXPECT_EXIT(range<TypeParam>::block_partition(len, i, (i-1) / 2), ::testing::KilledBySignal(SIGABRT), "error on line .* of block_partition()");
      EXPECT_EXIT(range<TypeParam>::block_partition(len, i, i-1), ::testing::KilledBySignal(SIGABRT), "error on line .* of block_partition()");
    }
  }

  lens = {std::numeric_limits<TypeParam>::lowest(), std::numeric_limits<TypeParam>::min(), 0};
  partitionCount = {0, std::numeric_limits<int>::lowest(), std::numeric_limits<int>::max(), std::numeric_limits<int>::min()};
  for (auto len : lens) {
    for (int i : partitionCount) {
      EXPECT_EXIT(range<TypeParam>::block_partition(len, i, 0), ::testing::KilledBySignal(SIGABRT), "error on line .* of block_partition()");
      EXPECT_EXIT(range<TypeParam>::block_partition(len, i, (i-1) / 2), ::testing::KilledBySignal(SIGABRT), "error on line .* of block_partition()");
      EXPECT_EXIT(range<TypeParam>::block_partition(len, i, i-1), ::testing::KilledBySignal(SIGABRT), "error on line .* of block_partition()");
    }
  }

}

TYPED_TEST_P(RangeTest, alignFails) {
  range<TypeParam> r;

  std::vector<TypeParam> starts = {0, 1, std::numeric_limits<TypeParam>::lowest(), std::numeric_limits<TypeParam>::min(),
                                   std::numeric_limits<TypeParam>::max()/3, std::numeric_limits<TypeParam>::max()/2, std::numeric_limits<TypeParam>::max()-1};
  std::vector<TypeParam> pageSizes = {0, std::numeric_limits<TypeParam>::lowest(), std::numeric_limits<TypeParam>::min()};

  for (auto s : starts) {
    for (auto p : pageSizes) {
      r = range<TypeParam>(s, s+1);
      EXPECT_EXIT(r.align_to_page(p), ::testing::KilledBySignal(SIGABRT), ".*test-bliss-iterators: /home/tpan/src/bliss/src/iterators/range.hpp:125: bliss::iterator::range<T> bliss::iterator::range<T>::align_to_page(const T&) const [with T = char]: Assertion `page_size > 0' failed.*");
//      EXPECT_DEATH(r.align_to_page(p), "test-bliss-iterators: /home/tpan/src/bliss/src/iterators/range.hpp:125: bliss::iterator::range<T> bliss::iterator::range<T>::align_to_page(const T&) const [with T = char]: Assertion `page_size > 0' failed.\n");
    }
  }

  starts = {std::numeric_limits<TypeParam>::max()};
  pageSizes = {1, std::numeric_limits<TypeParam>::max() / 2, std::numeric_limits<TypeParam>::max()};

  for (auto s : starts) {
    for (auto p : pageSizes) {
      r = range<TypeParam>(s, s+1);
//      EXPECT_EXIT(r.align_to_page(p), ::testing::KilledBySignal(SIGABRT), "error on line .* of align_to_page()");
      EXPECT_DEATH(r.align_to_page(p), "error on line .* of align_to_page()");
    }
  }

  starts = {std::numeric_limits<TypeParam>::max()};
  pageSizes = {0, std::numeric_limits<TypeParam>::lowest(), std::numeric_limits<TypeParam>::min()};

  for (auto s : starts) {
    for (auto p : pageSizes) {
      r = range<TypeParam>(s, s+1);
//      EXPECT_EXIT(r.align_to_page(p), ::testing::KilledBySignal(SIGABRT), "error on line .* of align_to_page()");
      EXPECT_DEATH(r.align_to_page(p), "error on line .* of align_to_page()");
    }
  }
}


REGISTER_TYPED_TEST_CASE_P(RangeTest, equal, assignment, copyConstruct, partition, align, partitionFails, alignFails);

//typedef ::testing::Types<char, uint8_t, int16_t, uint16_t, int, uint32_t, int64_t, uint64_t, size_t> RangeTestTypes;
typedef ::testing::Types<char> RangeTestTypes;
INSTANTIATE_TYPED_TEST_CASE_P(Bliss, RangeTest, RangeTestTypes);
