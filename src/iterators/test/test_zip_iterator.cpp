/**
 * test_zip_iterator.cpp
 * Test ZipIterator class
 *  Created on: Aug 27, 2014
 *      Author: Tony Pan <tpan7@gatech.edu>
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

    bliss::iterator::ZipIterator<CountingIterator<T>, ConstantIterator<T> > iter1;
    bliss::iterator::ZipIterator<CountingIterator<T>, CountingIterator<T> > iter2;

    virtual void SetUp()
    {
      iter1 = ZipIterator<CountingIterator<T>, ConstantIterator<T> >(CountingIterator<T>(4, 2), ConstantIterator<T>(3));
      iter2 = ZipIterator<CountingIterator<T>, CountingIterator<T> >(CountingIterator<T>(4, 2), CountingIterator<T>(3, 3));
    }
};

// indicate this is a typed test
TYPED_TEST_CASE_P(ZipIteratorTest);


// testing the assignment operator
TYPED_TEST_P(ZipIteratorTest, copy){
  // copy construct
  {
    decltype(this->iter1) iter(this->iter1);

    ASSERT_EQ(static_cast<TypeParam>(4), iter->first);
    ASSERT_EQ(static_cast<TypeParam>(3), iter->second);
  }

  // assignment
  {
    decltype(this->iter1) iter;
    iter = this->iter1;

    ASSERT_EQ(static_cast<TypeParam>(4), iter->first);
    ASSERT_EQ(static_cast<TypeParam>(3), iter->second);
  }

  // multipass
  {
    decltype(this->iter1) iter;
    iter = this->iter1;
    decltype(this->iter1) iter3;
    iter3 = iter;

    ASSERT_EQ(static_cast<TypeParam>(6), (++iter3)->first);   // dereference pre increment
    ASSERT_EQ(static_cast<TypeParam>(4), iter->first);
    ASSERT_EQ(static_cast<TypeParam>(6), (++iter)->first);   // dereference pre increment
    ASSERT_EQ(static_cast<TypeParam>(6), (iter++)->first);     // dereference post increment (creates a copy)
    ASSERT_EQ(static_cast<TypeParam>(6), iter3->first);       // check original hasn't been incremented after the post increment dereference
    ASSERT_EQ(static_cast<TypeParam>(8), iter->first);      // check copy is still same.
  }

}


// testing the copy constructor
TYPED_TEST_P(ZipIteratorTest, increment){

  decltype(this->iter2) iter(this->iter2);
  // increment
  ++iter;
  ASSERT_EQ(static_cast<TypeParam>(6), (*iter).first);  // pre increment
  iter++;
  ASSERT_EQ(static_cast<TypeParam>(8), (*iter).first);  // post increment but did not save
  ASSERT_EQ(static_cast<TypeParam>(8), (*(iter++)).first);  // post increment and use data
  ASSERT_EQ(static_cast<TypeParam>(10), (*iter).first);  // check original still same before post increment.
}

// failed construction due to asserts
TYPED_TEST_P(ZipIteratorTest, dereference){
  decltype(this->iter2) iter(this->iter2);

  // *a
  ASSERT_EQ(static_cast<TypeParam>(4), (*iter).first);
  ++iter;

  // a->
  ASSERT_EQ(static_cast<TypeParam>(6), iter->first);

}


// test page alignment
TYPED_TEST_P(ZipIteratorTest, compare){
  // equal
  decltype(this->iter1) iter(this->iter1);
  decltype(this->iter1) iter3(this->iter1);

  ASSERT_TRUE(iter == iter3);

  // not equal
  std::ptrdiff_t x;
  x = 3;
  std::advance<decltype(this->iter1), std::ptrdiff_t>(iter, x);
  x = 2;
  std::advance<decltype(this->iter1), std::ptrdiff_t>(iter3, x);

  ASSERT_TRUE(iter != iter3);

  ++iter3;
  ASSERT_TRUE(iter == iter3);


}


TYPED_TEST_P(ZipIteratorTest, traits) {

  bool same_type = std::is_same<typename std::iterator_traits<ZipIterator<CountingIterator<TypeParam>, ConstantIterator<TypeParam> > >::value_type,
      std::pair<TypeParam, TypeParam> >::value;

  ASSERT_TRUE(same_type);


}




// now register the test cases
REGISTER_TYPED_TEST_CASE_P(ZipIteratorTest, copy, increment, dereference, compare, traits);


//////////////////// RUN the tests with different types.

typedef ::testing::Types<int8_t, uint8_t, int16_t, uint16_t, int32_t, uint32_t,
    int64_t, uint64_t, size_t> ZipIteratorTestTypes;
INSTANTIATE_TYPED_TEST_CASE_P(Bliss, ZipIteratorTest, ZipIteratorTestTypes);
