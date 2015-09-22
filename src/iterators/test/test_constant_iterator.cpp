/**
 * test_counting_iterator.cpp
 * Test ConstantIterator class
 *  Created on: Aug 27, 2014
 *      Author: Tony Pan <tpan7@gatech.edu>
 */

#include "iterators/constant_iterator.hpp"

#include <gtest/gtest.h>
#include <cstdint>  // for uint64_t, etc.
#include <limits>

using namespace bliss::iterator;

/*
 * test class holding some information.  Also, needed for the typed tests
 */
template<typename T>
class ConstantIteratorTest : public ::testing::Test
{
  protected:

    virtual void SetUp()
    {
    }
};

// indicate this is a typed test
TYPED_TEST_CASE_P(ConstantIteratorTest);


// testing the equal function
TYPED_TEST_P(ConstantIteratorTest, construct){
  // default construct
  {
    ConstantIterator<TypeParam> iter;

    ASSERT_EQ(static_cast<TypeParam>(0), *iter);
  }
  // construct
  {
    ConstantIterator<TypeParam> iter(static_cast<TypeParam>(2));
    ASSERT_EQ(static_cast<TypeParam>(2), *iter);
  }

}




// testing the assignment operator
TYPED_TEST_P(ConstantIteratorTest, copy){
  // copy construct
  {
    ConstantIterator<TypeParam> iter(3);
    ConstantIterator<TypeParam> iter2(iter);

    ASSERT_EQ(static_cast<TypeParam>(3), *iter2);
  }

  // assignment
  {
    ConstantIterator<TypeParam> iter(2);
    ConstantIterator<TypeParam> iter2;
    iter2 = iter;

    ASSERT_EQ(static_cast<TypeParam>(2), *iter2);
  }

  // multipass
  {
    ConstantIterator<TypeParam> iter(3);
    ConstantIterator<TypeParam> iter2;
    iter2 = iter;

    ASSERT_EQ(static_cast<TypeParam>(3), *iter);
    ASSERT_EQ(static_cast<TypeParam>(3), *(++iter));   // dereference pre increment
    ASSERT_EQ(static_cast<TypeParam>(3), *iter++);     // dereference post increment (creates a copy)
    ASSERT_EQ(static_cast<TypeParam>(3), *iter);       // check original hasn't been incremented after the post increment dereference
    ASSERT_EQ(static_cast<TypeParam>(3), *iter2);      // check copy is still same.
  }

}


// testing the copy constructor
TYPED_TEST_P(ConstantIteratorTest, increment){
  ConstantIterator<TypeParam> iter(3);

  // increment
  ++iter;
  ASSERT_EQ(static_cast<TypeParam>(3), *iter);  // pre increment
  iter++;
  ASSERT_EQ(static_cast<TypeParam>(3), *iter);  // post increment but did not save
  ASSERT_EQ(static_cast<TypeParam>(3), *iter++);  // post increment and use data
  ASSERT_EQ(static_cast<TypeParam>(3), *iter);  // check original still same before post increment.


  // decrement
  ++iter;
  ASSERT_EQ(static_cast<TypeParam>(3), *iter);  // pre increment

  --iter;
  ASSERT_EQ(static_cast<TypeParam>(3), *iter);

  iter--;
  ASSERT_EQ(static_cast<TypeParam>(3), *iter);  // post increment but did not save
  ASSERT_EQ(static_cast<TypeParam>(3), *iter--);  // post increment and use data
  ASSERT_EQ(static_cast<TypeParam>(3), *iter);  // check original still same before post increment.

  // arithmetic
  ConstantIterator<TypeParam> iter2;
  iter2 = iter;                // *iter = 3

  ASSERT_EQ(static_cast<TypeParam>(3), *(iter + 3));  // 6 + 3 *3
  ASSERT_EQ(static_cast<TypeParam>(3), *(iter - 2));   // 6 - 3 * 2
  ASSERT_EQ(static_cast<TypeParam>(3), *(3 + iter));  // 3*3 + 6
  ASSERT_EQ(static_cast<ptrdiff_t>(0), iter - iter2);  // 6 - 6

  ConstantIterator<TypeParam> iter3(5);
  ASSERT_EQ(static_cast<ptrdiff_t>(2), iter3 - iter2);  // 6 - 6



  // compound assignment
  iter -= 2;                   // 6 - 3 *2
  ASSERT_EQ(static_cast<TypeParam>(3), *iter);

  iter += 3;                    // 0 + 3 * 3
  ASSERT_EQ(static_cast<TypeParam>(3), *iter);

  ASSERT_EQ(static_cast<ptrdiff_t>(0), iter - iter2);   // (9 - 6) / 3
}

// failed construction due to asserts
TYPED_TEST_P(ConstantIteratorTest, dereference){
  ConstantIterator<TypeParam> iter(4);
  // *a
  ASSERT_EQ(static_cast<TypeParam>(4), *iter);
  ++iter;

  // a[]
  ASSERT_EQ(static_cast<TypeParam>(4), iter[12]);
}


// test page alignment
TYPED_TEST_P(ConstantIteratorTest, compare){
  // equal
  ConstantIterator<TypeParam> iter(4);
  ConstantIterator<TypeParam> iter2(4);

  ASSERT_TRUE(iter == iter2);

  iter += 3;
  iter2 += 2;

  ASSERT_TRUE(iter == iter2);

  // not equal
  ++iter2;
  ASSERT_TRUE(iter == iter2);

  // gt, ls, ge, le
  ASSERT_FALSE(iter < iter2);
  ASSERT_FALSE(iter2 > iter);
  ASSERT_TRUE(iter <= iter2);
  ASSERT_TRUE(iter2 >= iter);

  --iter2;
  ASSERT_TRUE(iter <= iter2);
  ASSERT_TRUE(iter2 >= iter);


  ConstantIterator<TypeParam> iter3(5);

  ASSERT_FALSE(iter3 == iter2);
  ASSERT_TRUE(iter3 != iter2);

  // gt, ls, ge, le
  ASSERT_FALSE(iter3 < iter2);
  ASSERT_FALSE(iter2 > iter3);
  ASSERT_FALSE(iter3 <= iter2);
  ASSERT_FALSE(iter2 >= iter3);

}

// now register the test cases
REGISTER_TYPED_TEST_CASE_P(ConstantIteratorTest, construct, copy, increment, compare, dereference );


//////////////////// RUN the tests with different types.

typedef ::testing::Types<int8_t, uint8_t, int16_t, uint16_t, int32_t, uint32_t,
    int64_t, uint64_t, size_t> ConstantIteratorTestTypes;
INSTANTIATE_TYPED_TEST_CASE_P(Bliss, ConstantIteratorTest, ConstantIteratorTestTypes);
