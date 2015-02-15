/**
 * @file		test_constexpr_array.cpp
 * @ingroup
 * @author	Tony Pan <tpan7@gatech.edu>
 * @brief   tests const expre array using google test
 * @details
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */


#include "utils/constexpr_array.hpp"
#include <gtest/gtest.h>
#include <iostream>
#include <iterator>
#include <cmath>
#include <limits>

#include "config.hpp"

//// adapted from from http://stackoverflow.com/questions/10060110/how-does-gtest-compare-the-values-in-two-arrays
//template<typename T, size_t size>
//::testing::AssertionResult ArraysMatch(const std::array<T, size> & expected,
//                                       const std::array<T, size> & actual){
//    for (size_t i(0); i < size; ++i){
//        if (expected[i] != actual[i]){
//            return ::testing::AssertionFailure() << "array[" << i
//                << "] (" << actual[i] << ") != expected[" << i
//                << "] (" << expected[i] << ")";
//        }
//    }
//
//    return ::testing::AssertionSuccess();
//}


// set of functions
struct square_fn {
    constexpr double operator()(size_t x) {
      return x * x;
    }
};
struct sqrt_fn {
    constexpr double operator()(size_t x) {
      return sqrt(x);
    }
};

struct inverse_fn {
    constexpr double operator()(size_t x) {
      return 1.0/static_cast<double>(x + 1);
    }
};
struct log_fn {
    constexpr double operator()(size_t x) {
      return x == 0 ? std::numeric_limits<double>::lowest() : log(x);
    }
};



template<typename Func>
class ConstexprArrayTest : public ::testing::Test
{
  protected:
    virtual void SetUp() {
    };
};


// indicate this is a typed test
TYPED_TEST_CASE_P(ConstexprArrayTest);


TYPED_TEST_P(ConstexprArrayTest, compute) {
  constexpr auto N=10;
  constexpr auto a = bliss::utils::make_array<N>(TypeParam());  // OKAY to init one here, but not ahead of time and pass in object.

  constexpr TypeParam f;
  for (int i = 0; i < N; ++i) {
    EXPECT_EQ(f(i), a[i]);
  }
}

// now register the test cases
REGISTER_TYPED_TEST_CASE_P(ConstexprArrayTest, compute);

//////////////////// RUN the tests with different types.

typedef ::testing::Types<square_fn, inverse_fn, log_fn, sqrt_fn> ConstexprArrayTestTypes;
INSTANTIATE_TYPED_TEST_CASE_P(Bliss, ConstexprArrayTest, ConstexprArrayTestTypes);
