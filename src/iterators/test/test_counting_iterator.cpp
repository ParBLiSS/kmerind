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
TYPED_TEST_P(CountingIteratorTest, construct){
  // default construct
  {
    CountingIterator<TypeParam> iter;

    ASSERT_EQ(0, *iter);
  }
  // construct
  {
    CountingIterator<TypeParam> iter(static_cast<TypeParam>(2));
    ASSERT_EQ(2, *iter);
    ++iter;
    ASSERT_EQ(3, *iter);
  }

  // construct
  {
    CountingIterator<TypeParam> iter(static_cast<TypeParam>(3), 2);
    ASSERT_EQ(3, *iter);
    ++iter;
    ASSERT_EQ(5, *iter);
  }

}




// testing the assignment operator
TYPED_TEST_P(CountingIteratorTest, copy){
  // copy construct
  {
    CountingIterator<TypeParam> iter(3, 3);
    CountingIterator<TypeParam> iter2(iter);

    ASSERT_EQ(3, *iter2);
  }

  // assignment
  {
    CountingIterator<TypeParam> iter(2, 3);
    CountingIterator<TypeParam> iter2;
    iter2 = iter;

    ASSERT_EQ(2, *iter2);
  }

  // multipass
  {
    CountingIterator<TypeParam> iter(3, 3);
    CountingIterator<TypeParam> iter2;
    iter2 = iter;

    ASSERT_EQ(3, *iter);
    ASSERT_EQ(6, *(++iter));   // dereference pre increment
    ASSERT_EQ(9, *iter++);     // dereference post increment (creates a copy)
    ASSERT_EQ(6, *iter);       // check original hasn't been incremented after the post increment dereference
    ASSERT_EQ(3, *iter2);      // check copy is still same.
  }

}


// testing the copy constructor
TYPED_TEST_P(CountingIteratorTest, increment){
  CountingIterator<TypeParam> iter(3, 3);

  // increment
  ++iter;
  ASSERT_EQ(6, *iter);  // pre increment
  iter++;
  ASSERT_EQ(6, *iter);  // post increment but did not save
  ASSERT_EQ(9, *iter++);  // post increment and use data
  ASSERT_EQ(6, *iter);  // check original still same before post increment.


  // decrement
  ++iter;               // get it to 9.
  ASSERT_EQ(9, *iter);  // pre increment

  --iter;
  ASSERT_EQ(6, *iter);

  iter--;
  ASSERT_EQ(6, *iter);  // post increment but did not save
  ASSERT_EQ(3, *iter--);  // post increment and use data
  ASSERT_EQ(6, *iter);  // check original still same before post increment.

  // arithmetic
  CountingIterator<TypeParam> iter2;
  iter2 = iter;                // *iter = 6

  ASSERT_EQ(15, *(iter + 3));  // 6 + 3 *3
  ASSERT_EQ(0, *(iter - 2));   // 6 - 3 * 2
  ASSERT_EQ(15, *(3 + iter));  // 3*3 + 6
  ASSERT_EQ(0, iter - iter2);  // 6 - 6

  // compound assignment
  iter -= 2;                   // 6 - 3 *2
  ASSERT_EQ(0, *iter);

  iter += 3;                    // 0 + 3 * 3
  ASSERT_EQ(9, *iter);

  ASSERT_EQ(1, iter - iter2);   // (9 - 6) / 3
}

// failed construction due to asserts
TYPED_TEST_P(CountingIteratorTest, dereference){
  CountingIterator<TypeParam> iter(4, 2);
  // *a
  ASSERT_EQ(4, *iter);
  ++iter;

  // a[]
  ASSERT_EQ(30, iter[12]);
}


// test page alignment
TYPED_TEST_P(CountingIteratorTest, compare){
  // equal
  CountingIterator<TypeParam> iter(4, 2);
  CountingIterator<TypeParam> iter2(4, 3);

  ASSERT_TRUE(iter == iter2);

  iter += 3;
  iter2 += 2;

  ASSERT_TRUE(iter == iter2);

  // not equal
  ++iter2;
  ASSERT_TRUE(iter != iter2);

  // gt, ls, ge, le
  ASSERT_TRUE(iter < iter2);
  ASSERT_TRUE(iter2 > iter);
  ASSERT_TRUE(iter <= iter2);
  ASSERT_TRUE(iter2 >= iter);

  --iter2;
  ASSERT_TRUE(iter <= iter2);
  ASSERT_TRUE(iter2 >= iter);


}







// now register the test cases
REGISTER_TYPED_TEST_CASE_P(CountingIteratorTest, construct, copy, increment, compare, dereference );


//////////////////// RUN the tests with different types.

typedef ::testing::Types<int8_t, uint8_t, int16_t, uint16_t, int32_t, uint32_t,
    int64_t, uint64_t, size_t> CountingIteratorTestTypes;
INSTANTIATE_TYPED_TEST_CASE_P(Bliss, CountingIteratorTest, CountingIteratorTestTypes);
