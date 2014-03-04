/**
 * range_test.cpp
 *
 *  Created on: Feb 18, 2014
 *      Author: tpan
 */

#include "iterators/range.hpp"

#include <gtest/gtest.h>
#include <cstdint>

#include "config.hpp"

template <typename T>
class RangeTest : public ::testing::Test {

};


TYPED_TEST_CASE_P(RangeTest);


TYPED_TEST_P(RangeTest, test1) {

}

TYPED_TEST_P(RangeTest, test1) {

}


REGISTER_TYPED_TEST_CASE_P(RangeTest, test1, test2);

typedef ::testing::Types<char, uint8_t, int16_t, uint16_t, int, uint32_t, int64_t, uint64_t, size_t> RangeTestTypes;
INSTANTIATE_TYPED_TEST_CASE_P(Bliss, RangeTest, RangeTestTypes);
