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
 * range_test.cpp
 * Test range class
 *  Created on: Feb 18, 2014
 *      Author: Tony Pan <tpan7@gatech.edu>
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
class RangeTest : public ::testing::Test
{
  protected:
    size_t page_size;

    virtual void SetUp()
    {
      page_size = sysconf(_SC_PAGE_SIZE);
    }
};

// indicate this is a typed test
TYPED_TEST_CASE_P(RangeTest);


// testing the equal function
TYPED_TEST_P(RangeTest, equal){
  range<TypeParam> r(0, 100);
  range<TypeParam> r2(0, 100);
  EXPECT_TRUE(r == r2);

  range<TypeParam> r3(10, 100);
  range<TypeParam> r4(10, 100);
  EXPECT_TRUE(r3 == r4);

  EXPECT_FALSE(r == r4);

  // if signed type, then range can have negative start and end
  if (std::is_signed<TypeParam>::value)
  {
    range<TypeParam> r(-10, 100);
    range<TypeParam> r2(-10, 100);
    EXPECT_TRUE(r == r2);

    range<TypeParam> r3(-101, -100);
    range<TypeParam> r4(-101, -100);
    EXPECT_TRUE(r3 == r4);

    EXPECT_FALSE(r == r4);
  }
}

// testing the assignment operator
TYPED_TEST_P(RangeTest, assignment){
  range<TypeParam> r;
  range<TypeParam> r2(10, 100);
  r = r2;
  EXPECT_TRUE(r == r2);

  if (std::is_signed<TypeParam>::value)
  {
    range<TypeParam> r3(-10, 100);
    r = r3;
    EXPECT_TRUE(r == r3);
  }

}


// testing the copy constructor
TYPED_TEST_P(RangeTest, copyConstruct){
  range<TypeParam> r2(10, 100);
  range<TypeParam> r(r2);
  EXPECT_TRUE(r == r2);

  if (std::is_signed<TypeParam>::value)
  {
    range<TypeParam> r3(-10, 100);
    range<TypeParam> r4(r3);
    EXPECT_TRUE(r3 == r4);
  }

}



// test page alignment
TYPED_TEST_P(RangeTest, align){
  range<TypeParam> r;

  std::vector<TypeParam> starts =
  { 0, 1, (std::numeric_limits<TypeParam>::max() >> 1) + 1, std::numeric_limits<TypeParam>::max()-1};

  if (std::is_signed<TypeParam>::value)
  {
    starts.push_back(-1);
  }

  size_t size = 1;
  std::vector<size_t> pageSizes =
  { 1, 64};

  for (auto s : starts)
  {
    for (auto p : pageSizes)
    {
      // don't test cases where alignment would push us lower than the lowest possible value for the type.
      if (static_cast<size_t>(s - std::numeric_limits<TypeParam>::lowest()) < p)
        continue;

      //BL_INFOF("1 %ld %ld\n", static_cast<int64_t>(s), static_cast<int64_t>(p));
      r = range<TypeParam>(s, s+size);
      TypeParam val = range<TypeParam>::align_to_page(r, p);
      EXPECT_TRUE(range<TypeParam>::is_page_aligned(val, p));
    }
  }

}


// failed construction due to asserts
TYPED_TEST_P(RangeTest, constructFails){

  // basically, if start is larger than end.
  try {
    range<TypeParam>(std::numeric_limits<TypeParam>::max(), std::numeric_limits<TypeParam>::min());
    ADD_FAILURE();
  } catch (const std::invalid_argument &e) {
    // did throw exception, so okay.
  }

  try {
    range<TypeParam>(std::numeric_limits<TypeParam>::max(), std::numeric_limits<TypeParam>::lowest());
    ADD_FAILURE();
  } catch (const std::invalid_argument &e) {
    // did throw exception, so okay.
  }



}




// failed alignment
TYPED_TEST_P(RangeTest, alignFails){
  range<TypeParam> r;

  std::vector<TypeParam> starts =
  { 0, 1, (std::numeric_limits<TypeParam>::max() >> 1) + 1, std::numeric_limits<TypeParam>::max() - 1};

  if (std::is_signed<TypeParam>::value)
  {
    starts.push_back(std::numeric_limits<TypeParam>::min());
    starts.push_back(std::numeric_limits<TypeParam>::lowest());
  }
  size_t size= 1;
  std::vector<size_t> pageSizes =
  { 0, std::numeric_limits<size_t>::lowest(), std::numeric_limits<size_t>::min()};

  std::string err_regex = ".*range.hpp.*align_to_page.* Assertion .* failed.*";

  for (auto s : starts)
  {
    for (auto p : pageSizes)
    {
      //BL_INFOF("align fail processing start %ld page size %lud\n", static_cast<int64_t>(s), static_cast<uint64_t>(p));

      // align fails because of bad page sizes (0 or negative
      r = range<TypeParam>(s, s+size);
      try {
        range<TypeParam>::align_to_page(r, p);
        ADD_FAILURE();
      } catch (const std::range_error &e) {
        // okay.  threw the right error
      } catch (const std::invalid_argument &e) {
        // okay.  threw the right error
      }

    }

    // if s is negative, also check to make sure we fail with large p.
    if (s < 0) {
      r = range<TypeParam>(s, s+size);

      try {
        range<TypeParam>::align_to_page(r, std::numeric_limits<size_t>::max());
        ADD_FAILURE();
      } catch (const std::range_error &e) {
        // okay.  threw the right error
      } catch (const std::invalid_argument &e) {
        // okay.  threw the right error
      }
    }
  }

}



// now register the test cases
REGISTER_TYPED_TEST_CASE_P(RangeTest, equal, assignment, copyConstruct, align, constructFails, alignFails);


//////////////////// RUN the tests with different types.

typedef ::testing::Types<int8_t, uint8_t, int16_t, uint16_t, int32_t, uint32_t,
    int64_t, uint64_t, size_t> RangeTestTypes;
INSTANTIATE_TYPED_TEST_CASE_P(Bliss, RangeTest, RangeTestTypes);
