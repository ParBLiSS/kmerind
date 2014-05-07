/**
 * range_test.cpp
 * Test range class
 *  Created on: Feb 18, 2014
 *      Author: tpan
 */

#include "io/data_block.hpp"
#include "iterators/range.hpp"

#include <gtest/gtest.h>
#include <cstdint>  // for uint64_t, etc.
#include <vector>
#include <limits>
#include <algorithm>

using namespace bliss::io;

/*
 * test class holding some information.  Also, needed for the typed tests
 */
template<typename T>
class DataBlockTest : public ::testing::Test
{
  protected:
    std::vector<T> src;
    std::vector<T> dest;
    size_t slen;
    size_t dlen;

    virtual void SetUp()
    {
      slen = (std::numeric_limits<T>::max() < 1000000) ? std::numeric_limits<T>::max() : 1000000;

      src.reserve(slen);
      for (T i = 0; i < slen; ++i) {
        src.emplace_back(i);
      }

      dlen = (std::numeric_limits<T>::max() < 1000) ? std::numeric_limits<T>::max() : 1000;
      dest.reserve(dlen);
    }
};

// indicate this is a typed test
TYPED_TEST_CASE_P(DataBlockTest);


// testing the equal function
TYPED_TEST_P(DataBlockTest, NoBuffer){
  typedef bliss::io::DataBlock<typename std::vector<TypeParam>::iterator, bliss::iterator::range<size_t>> NoBuffer_DB_Type;

  bliss::iterator::range<size_t> r(0, this->slen);
  NoBuffer_DB_Type db(this->src.begin(), this->src.end(), r);
  EXPECT_EQ(this->src.begin(), db.begin());
  EXPECT_EQ(this->src.end(),   db.end()  );

  int i = r.start;
  for (auto it = db.begin(); it != db.end(); ++it, ++i) {
    EXPECT_EQ(this->src[i], *it) << "Vectors x and y differ at index " << i;
  }

  bliss::iterator::range<size_t> r2(10, this->slen / 2);
  NoBuffer_DB_Type db2(this->src.begin() + r2.start, this->src.begin() + r2.end, r2);
  EXPECT_EQ(this->src.begin()+ r2.start, db2.begin());
  EXPECT_EQ(this->src.begin()+ r2.end  , db2.end()  );

  i = r2.start;
  for (auto it = db2.begin(); it != db2.end(); ++it, ++i) {
    EXPECT_EQ(this->src[i], *it) << "Vectors x and y differ at index " << i;
  }

}

// testing the equal function
TYPED_TEST_P(DataBlockTest, Buffer){
  typedef bliss::io::DataBlock<typename std::vector<TypeParam>::iterator, bliss::iterator::range<size_t>> Buffer_DB_Type;

  bliss::iterator::range<size_t> r2(10, 10 + this->dlen);
  Buffer_DB_Type db2(this->src.begin() + r2.start, this->src.begin() + r2.end, r2, this->dest.begin(), this->dlen);
  EXPECT_NE(this->src.begin()+ r2.start, db2.begin());
  EXPECT_NE(this->src.begin()+ r2.end  , db2.end()  );

  int i = r2.start;
  for (auto it = db2.begin(); it != db2.end(); ++it, ++i) {
    EXPECT_EQ(this->src[i], *it) << "Vectors x and y differ at index " << i;
  }
}

// now register the test cases
REGISTER_TYPED_TEST_CASE_P(DataBlockTest, NoBuffer, Buffer);


/*
 * test class holding some information.  Also, needed for the typed tests
 */
template<typename T>
class DataBlockPtrTest : public ::testing::Test
{
  protected:
    T* src;
    T* dest;
    size_t slen;
    size_t dlen;

    virtual void SetUp()
    {
      slen = (std::numeric_limits<T>::max() < 1000000) ? std::numeric_limits<T>::max() : 1000000;

      src = new T[slen];
      for (T i = 0; i < slen; ++i) {
        src[i] = (i);
      }

      dlen = (std::numeric_limits<T>::max() < 1000) ? std::numeric_limits<T>::max() : 1000;
      dest = new T[dlen];
    }

    virtual void TearDown()
    {
      delete [] src;
      delete [] dest;
    }
};

// indicate this is a typed test
TYPED_TEST_CASE_P(DataBlockPtrTest);


// testing the equal function
TYPED_TEST_P(DataBlockPtrTest, NoBuffer){
  typedef bliss::io::DataBlock<TypeParam*, bliss::iterator::range<size_t>> NoBuffer_DB_Type;

  bliss::iterator::range<size_t> r(0, this->slen);
  NoBuffer_DB_Type db(this->src, this->src + this->slen, r);
  EXPECT_EQ(this->src,                db.begin());
  EXPECT_EQ(this->src + this->slen,   db.end()  );

  int i = r.start;
  for (auto it = db.begin(); it != db.end(); ++it, ++i) {
    EXPECT_EQ(this->src[i], *it) << "Vectors x and y differ at index " << i;
  }

  bliss::iterator::range<size_t> r2(10, this->slen / 2);
  NoBuffer_DB_Type db2(this->src + r2.start, this->src + r2.end, r2);
  EXPECT_EQ(this->src + r2.start, db2.begin());
  EXPECT_EQ(this->src + r2.end  , db2.end()  );

  i = r2.start;
  for (auto it = db2.begin(); it != db2.end(); ++it, ++i) {
    EXPECT_EQ(this->src[i], *it) << "Vectors x and y differ at index " << i;
  }

}

// testing the equal function
TYPED_TEST_P(DataBlockPtrTest, Buffer){
  typedef bliss::io::DataBlock<TypeParam*, bliss::iterator::range<size_t>> Buffer_DB_Type;

  bliss::iterator::range<size_t> r2(10, 10 + this->dlen);
  Buffer_DB_Type db2(this->src + r2.start, this->src + r2.end, r2, this->dest, this->dlen);
  EXPECT_NE(this->src + r2.start, db2.begin());
  EXPECT_NE(this->src + r2.end  , db2.end()  );

  int i = r2.start;
  for (auto it = db2.begin(); it != db2.end(); ++it, ++i) {
    EXPECT_EQ(this->src[i], *it) << "Vectors x and y differ at index " << i;
  }
}

// now register the test cases
REGISTER_TYPED_TEST_CASE_P(DataBlockPtrTest, NoBuffer, Buffer);



////////////////////////
//  DEATH TESTS - test class named BlahDeathTest so gtest will not run these in a threaded context. and will run first

// typedef RangeTest DataBlockDeathTest
template<typename T>
class DataBlockDeathTest : public ::testing::Test
{
  protected:
    std::vector<T> src;
    std::vector<T> dest;
    size_t slen;
    size_t dlen;

    virtual void SetUp()
    {
      slen = (std::numeric_limits<T>::max() < 1000000) ? std::numeric_limits<T>::max() : 1000000;

      src.reserve(slen);
      for (T i = 0; i < slen; ++i) {
        src.emplace_back(i);
      }

      dlen = (std::numeric_limits<T>::max() < 1000) ? std::numeric_limits<T>::max() : 1000;
      dlen /= 2;
      dest.reserve(dlen);
    }
};

// annotate that DataBlockDeathTest is typed
TYPED_TEST_CASE_P(DataBlockDeathTest);

// failed construction due to asserts
TYPED_TEST_P(DataBlockDeathTest, badRange){

  typedef bliss::io::DataBlock<typename std::vector<TypeParam>::iterator, bliss::iterator::range<size_t>> Buffer_DB_Type;

  bliss::iterator::range<size_t> r2(0, this->slen);
  std::string err_regex = ".data_block.hpp.* Assertion .* failed.*";

  // basically, if start is larger than end.
  EXPECT_EXIT(Buffer_DB_Type(this->src.begin(), this->src.begin() + this->dlen, r2, this->dest.begin(), this->dlen), ::testing::KilledBySignal(SIGABRT), err_regex);
}



// register the death test cases
REGISTER_TYPED_TEST_CASE_P(DataBlockDeathTest, badRange);


//////////////////// RUN the tests with different types.

typedef ::testing::Types<int8_t, uint8_t, int16_t, uint16_t, int32_t, uint32_t,
    int64_t, uint64_t, size_t> DataBlockTestTypes;
INSTANTIATE_TYPED_TEST_CASE_P(Bliss, DataBlockTest, DataBlockTestTypes);
INSTANTIATE_TYPED_TEST_CASE_P(Bliss, DataBlockPtrTest, DataBlockTestTypes);
INSTANTIATE_TYPED_TEST_CASE_P(Bliss, DataBlockDeathTest, DataBlockTestTypes);
