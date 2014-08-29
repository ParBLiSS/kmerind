/**
 * test_zip_iterator.cpp
 * Test ZipIterator class
 *  Created on: Aug 27, 2014
 *      Author: tpan
 */

#include "iterators/zip_iterator.hpp"
#include "iterators/constant_iterator.hpp"
#include "iterators/counting_iterator.hpp"

#include <gtest/gtest.h>
#include <cstdint>  // for uint64_t, etc.
#include <limits>

using namespace bliss::iterator;

/*
 * test class holding some information.  Also, needed for the typed tests
 */
template<typename T>
class ZipIteratorTest : public ::testing::Test
{
  protected:

    bliss::iterator::ZipIterator<CountingIterator, ConstantIterator> iter1;
    bliss::iterator::ZipIterator<CountingIterator, CountingIterator> iter2;

    virtual void SetUp()
    {
      iter1 = ZipIterator<CountingIterator, ConstantIterator>(CountingIterator<T>(4, 2), ConstantIterator<T>(3));
      iter2 = ZipIterator<CountingIterator, CountingIterator>(CountingIterator<T>(4, 2), CountingIterator<T>(3, 3));
    }
};

// indicate this is a typed test
TYPED_TEST_CASE_P(ZipIteratorTest);


// testing the assignment operator
TYPED_TEST_P(ZipIteratorTest, copy){
  // copy construct
  {
    typename decltype(this->iter1) iter(this->iter1);

    ASSERT_EQ(4, iter->first);
    ASSERT_EQ(3, iter->second);
  }

  // assignment
  {
    typename decltype(this->iter1) iter;
    iter = this->iter1;

    ASSERT_EQ(4, iter->first);
    ASSERT_EQ(3, iter->second);
  }

  // multipass
  {
    typename decltype(this->iter1) iter;
    iter = this->iter1;

    ASSERT_EQ(4, iter->first);
    ASSERT_EQ(6, (++iter)->first);   // dereference pre increment
    ASSERT_EQ(8, (iter++)->first);     // dereference post increment (creates a copy)
    ASSERT_EQ(4, this->iter1->first);       // check original hasn't been incremented after the post increment dereference
    ASSERT_EQ(8, iter->first);      // check copy is still same.
  }

}


// testing the copy constructor
TYPED_TEST_P(ZipIteratorTest, increment){

  // increment
  ++this->iter2;
  ASSERT_EQ(6, *this->iter2.first);  // pre increment
  this->iter2++;
  ASSERT_EQ(6, *this->iter2.first);  // post increment but did not save
  ASSERT_EQ(9, *this->iter2++.first);  // post increment and use data
  ASSERT_EQ(6, *this->iter2.first);  // check original still same before post increment.


  // decrement
  ++this->iter2;               // get it to 9.
  ASSERT_EQ(9, *this->iter2.first);  // pre increment

  --this->iter2;
  ASSERT_EQ(6, *this->iter2.first);

  this->iter2--;
  ASSERT_EQ(6, *this->iter2.first);  // post increment but did not save
  ASSERT_EQ(3, *this->iter2--.first);  // post increment and use data
  ASSERT_EQ(6, *this->iter2.first);  // check original still same before post increment.

  // arithmetic
  ASSERT_EQ(15, *(iter + 3));  // 6 + 3 *3
  ASSERT_EQ(0, *(iter - 2));   // 6 - 3 * 2
  ASSERT_EQ(15, *(3 + iter));  // 3*3 + 6
  ASSERT_EQ(0, iter - iter2);  // 6 - 6

  // compound assignment
  iter -= 2;                   // 6 - 3 *2
  ASSERT_EQ(0, *iter);

  iter += 3;                    // 0 + 3 * 3
  ASSERT_EQ(9, *iter);

  ASSERT_EQ(3, iter - iter2);   // 9 - 6
}

// failed construction due to asserts
TYPED_TEST_P(ZipIteratorTest, dereference){
  ZipIterator<TypeParam> iter(4, 2);
  // *a
  ASSERT_EQ(4, *iter);
  ++iter;

  // a[]
  ASSERT_EQ(30, iter[12]);
}


// test page alignment
TYPED_TEST_P(ZipIteratorTest, compare){
  // equal
  ZipIterator<TypeParam> iter(4, 2);
  ZipIterator<TypeParam> iter2(4, 3);

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
REGISTER_TYPED_TEST_CASE_P(ZipIteratorTest, copy, increment, compare, dereference );


//////////////////// RUN the tests with different types.

typedef ::testing::Types<int8_t, uint8_t, int16_t, uint16_t, int32_t, uint32_t,
    int64_t, uint64_t, size_t> ZipIteratorTestTypes;
INSTANTIATE_TYPED_TEST_CASE_P(Bliss, ZipIteratorTest, ZipIteratorTestTypes);
