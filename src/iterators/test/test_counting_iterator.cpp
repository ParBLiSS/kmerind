/**
 * test_counting_iterator.cpp
 * Test CountingIterator class
 *  Created on: Aug 27, 2014
 *      Author: tpan
 */

#include "iterators/counting_iterator.hpp"

#include <gtest/gtest.h>
#include <cstdint>  // for uint64_t, etc.
#include <limits>

using namespace bliss::iterator;

/*
 * test class holding some information.  Also, needed for the typed tests
 */
template<typename T>
class CountingIteratorTest : public ::testing::Test
{
  protected:

    virtual void SetUp()
    {
    }
};

// indicate this is a typed test
TYPED_TEST_CASE_P(CountingIteratorTest);


// testing the equal function
TYPED_TEST_P(CountingIteratorTest, equal){
}

// testing the assignment operator
TYPED_TEST_P(CountingIteratorTest, assignment){
}


// testing the copy constructor
TYPED_TEST_P(CountingIteratorTest, copyConstruct){
}



// test page alignment
TYPED_TEST_P(CountingIteratorTest, align){
}


// failed construction due to asserts
TYPED_TEST_P(CountingIteratorTest, constructFails){
}




// failed alignment
TYPED_TEST_P(CountingIteratorTest, alignFails){
}



// now register the test cases
REGISTER_TYPED_TEST_CASE_P(CountingIteratorTest, equal, assignment, copyConstruct, align, constructFails, alignFails);


//////////////////// RUN the tests with different types.

typedef ::testing::Types<int8_t, uint8_t, int16_t, uint16_t, int32_t, uint32_t,
    int64_t, uint64_t, size_t> CountingIteratorTestTypes;
INSTANTIATE_TYPED_TEST_CASE_P(Bliss, CountingIteratorTest, CountingIteratorTestTypes);
