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
 * test_filter_iterator.cpp
 * Test filter_iterator class
 *  Created on: Aug 27, 2014
 *      Author: Tony Pan <tpan7@gatech.edu>
 */

#include "iterators/filter_iterator.hpp"
#include "iterators/counting_iterator.hpp"

#include <gtest/gtest.h>
#include <cstdint>  // for uint64_t, etc.
#include <limits>

using namespace bliss::iterator;

/*
 * test class holding some information.  Also, needed for the typed tests
 */
template<typename T>
class FilterIteratorTest : public ::testing::Test
{
    static_assert(std::is_integral<T>::value, "data type has to be integral.");

  protected:

    struct Even {
      bool operator()(T const & v) {
        return (v % 2) == 0;
      }
    };

    struct Odd {
      bool operator()(T const & v) {
        return (v % 2) == 1;
      }
    };

    CountingIterator<T> base_iter_begin;
    CountingIterator<T> base_iter_end;
    filter_iterator<Even, CountingIterator<T> > even_filter_begin;
    filter_iterator<Even, CountingIterator<T> > even_filter_end;
    filter_iterator<Odd, CountingIterator<T> > odd_filter_begin;
    filter_iterator<Odd, CountingIterator<T> > odd_filter_end;

    virtual void SetUp()
    {
      base_iter_begin = CountingIterator<T>(0, 1);
      base_iter_end = CountingIterator<T>(1000, 1);
      even_filter_begin = filter_iterator<Even, CountingIterator<T> >(Even(), base_iter_begin, base_iter_end);
      even_filter_end = filter_iterator<Even, CountingIterator<T> >(Even(), base_iter_end);
      odd_filter_begin = filter_iterator<Odd, CountingIterator<T> >(Odd(), base_iter_begin, base_iter_end);
      odd_filter_end = filter_iterator<Odd, CountingIterator<T> >(Odd(), base_iter_end);
    }

};

// indicate this is a typed test
TYPED_TEST_CASE_P(FilterIteratorTest);



// testing the copy constructor
TYPED_TEST_P(FilterIteratorTest, increment){

	std::cout << "base: [" << *(this->base_iter_begin) << ", " << *(this->base_iter_end) << "]" << std::endl;


  TypeParam i = 0;
  auto it = this->even_filter_begin;
  for (; it != this->even_filter_end; ++it, i += 2) {
    ASSERT_TRUE(*it % 2 == 0);
    ASSERT_EQ(*it, i);
  }

  // decrement
  --it;
  i -= 2;
  for (; it != this->even_filter_begin; --it, i -= 2) {
    ASSERT_TRUE(*it % 2 == 0);
    ASSERT_EQ(*it, i);
  }


  i = 1;
  auto it2 = this->odd_filter_begin;
  for (; it2 != this->odd_filter_end; ++it2, i += 2) {
    ASSERT_TRUE(*it2 % 2 == 1);
    ASSERT_EQ(*it2, i);
  }

  --it2;
  i -= 2;
  for (; it2 != this->odd_filter_begin; --it2, i -= 2) {
    ASSERT_TRUE(*it2 % 2 == 1);
    ASSERT_EQ(*it2, i);
  }


}






// now register the test cases
REGISTER_TYPED_TEST_CASE_P(FilterIteratorTest, increment);


//////////////////// RUN the tests with different types.

typedef ::testing::Types<int32_t, uint32_t> FilterIteratorTestTypes;
INSTANTIATE_TYPED_TEST_CASE_P(Bliss, FilterIteratorTest, FilterIteratorTestTypes);
