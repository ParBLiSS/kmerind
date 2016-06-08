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


// include google test
#include <gtest/gtest.h>
#include "containers/densehash_map.hpp"

#include <string>
#include <unordered_map>
#include <random>
#include <iterator>  // for ostream iterator.
#include <algorithm>  // for transform.

// include files to test
#include "utils/logging.h"

  template <typename T>
  struct IsLower {

      static_assert(std::is_integral<T>::value, "only supporting integral types in tests right now.");

      static constexpr T high_bit_mask = static_cast<T>(0x80);   // only allowing 256 values

      bool operator()(T const & x) const {
        return (x & high_bit_mask) == 0; // lower is when high bit is 0 - so works with signed and unsigned.
      }
      bool operator()(::std::pair<T,T> const & x) const {
        return (x.first & high_bit_mask) == 0; // lower is when high bit is 0 - so works with signed and unsigned.
      }
    };


/*
 * test class holding some information.  Also, needed for the typed tests
 */
template<typename T>
class DenseHashMapTest : public ::testing::Test
{
    static_assert(std::is_integral<T>::value, "only supporting integral types in tests right now.");
  protected:


    ::std::unordered_map<T, T> gold;
    ::std::vector<std::pair<T, T>> temp;


    size_t iters = 10000;
    int max_val = 254;

    virtual void SetUp()
    { // generate some inputs


      std::default_random_engine generator;
      std::uniform_int_distribution<T> distribution(2, max_val);

      for (size_t i=0; i< iters; ++i) {
        T key = distribution(generator);
        T val = distribution(generator);
        gold.emplace(key, val);
        temp.emplace_back(::std::move(key), ::std::move(val));
      }

    }
};

// indicate this is a typed test
TYPED_TEST_CASE_P(DenseHashMapTest);

TYPED_TEST_P(DenseHashMapTest, insert_partial)
{
  bool same = false;

  ::fsc::densehash_map<TypeParam, TypeParam, ::fsc::TruePredicate> test(this->temp.begin(), this->temp.end(), 255, 254);

      ::std::vector<::std::pair<TypeParam, TypeParam> > test_vals = test.to_vector();
      ::std::vector<::std::pair<TypeParam, TypeParam> > gold_vals(this->gold.begin(), this->gold.end());


      ::std::sort(test_vals.begin(), test_vals.end(), [](::std::pair<TypeParam, TypeParam> const & x, ::std::pair<TypeParam, TypeParam> const &y) {
        return (x.first == y.first) ? (x.second < y.second) : (x.first < y.first);
      } );
      ::std::sort(gold_vals.begin(), gold_vals.end(), [](::std::pair<TypeParam, TypeParam> const & x, ::std::pair<TypeParam, TypeParam> const &y) {
        return (x.first == y.first) ? (x.second < y.second) : (x.first < y.first);
      } );

      same = ::std::equal(test_vals.begin(), test_vals.end(), gold_vals.begin());
//
//      if (!same) {
//        for (int i = 0; i < this->iters; ++i) {
//          printf("%ld->%ld\t%ld->%ld\n", test_vals[i].first, test_vals[i].second, gold_vals[i].first, gold_vals[i].second);
//        }
//      }

      EXPECT_TRUE(same);
}


TYPED_TEST_P(DenseHashMapTest, equal_range_partial)
{
  ::fsc::densehash_map<TypeParam, TypeParam, ::fsc::TruePredicate> test(this->temp.begin(), this->temp.end(), 255, 254);

  bool same = false;
	  for (int i = 0; i < this->max_val; ++i) {
	    auto test_range = test.equal_range(i);
	    auto gold_range = this->gold.equal_range(i);


	    ::std::vector<TypeParam> test_vals;
	    ::std::vector<TypeParam> gold_vals;

	    int jmax = 0;
	    for (auto it = test_range.first; it != test_range.second; ++it) {
	      test_vals.push_back((*it).second);
	      ++jmax;
	    }

	    int kmax = 0;
      for (auto it = gold_range.first; it != gold_range.second; ++it) {
        gold_vals.push_back(it->second);
        ++kmax;
      }

      EXPECT_EQ(jmax, kmax);

      ::std::sort(test_vals.begin(), test_vals.end());
      ::std::sort(gold_vals.begin(), gold_vals.end());

	    same = ::std::equal(test_vals.begin(), test_vals.end(), gold_vals.begin());

//	    if (!same) {
//	      printf("test %d\n", i);
//	      for (int j = 0; j < jmax; ++j) {
//	        printf("%ld\t%ld\n", test_vals[j], gold_vals[i]);
//	      }
//	    }

		  EXPECT_TRUE(same);
	  }
}


TYPED_TEST_P(DenseHashMapTest, count_partial)
{
  ::fsc::densehash_map<TypeParam, TypeParam, ::fsc::TruePredicate> test(this->temp.begin(), this->temp.end(), 255, 254);

    for (int i = 0; i < this->max_val; ++i) {
      EXPECT_EQ(this->gold.count(i), test.count(i));
    }
}


TYPED_TEST_P(DenseHashMapTest, insert_full)
{
  bool same = false;

  ::fsc::densehash_map<TypeParam, TypeParam, IsLower<TypeParam>> test(this->temp.begin(), this->temp.end(), 255, 254, 0, 1);

      ::std::vector<::std::pair<TypeParam, TypeParam> > test_vals = test.to_vector();
      ::std::vector<::std::pair<TypeParam, TypeParam> > gold_vals(this->gold.begin(), this->gold.end());

      EXPECT_EQ(test_vals.size(), gold_vals.size());

      ::std::sort(test_vals.begin(), test_vals.end(), [](::std::pair<TypeParam, TypeParam> const & x, ::std::pair<TypeParam, TypeParam> const &y) {
        return (x.first == y.first) ? (x.second < y.second) : (x.first < y.first);
      } );
      ::std::sort(gold_vals.begin(), gold_vals.end(), [](::std::pair<TypeParam, TypeParam> const & x, ::std::pair<TypeParam, TypeParam> const &y) {
        return (x.first == y.first) ? (x.second < y.second) : (x.first < y.first);
      } );

      same = ::std::equal(test_vals.begin(), test_vals.end(), gold_vals.begin());

//      if (!same) {
//        for (size_t i = 0; i < test_vals.size(); ++i) {
//          printf("%ld->%ld\t%ld->%ld\n", test_vals[i].first, test_vals[i].second, gold_vals[i].first, gold_vals[i].second);
//        }
//      }

      EXPECT_TRUE(same);
}


TYPED_TEST_P(DenseHashMapTest, equal_range_full)
{
  ::fsc::densehash_map<TypeParam, TypeParam, IsLower<TypeParam> > test(this->temp.begin(), this->temp.end(), 255, 254, 0, 1);

  bool same = false;
    for (int i = 0; i < this->max_val; ++i) {
      auto test_range = test.equal_range(i);
      auto gold_range = this->gold.equal_range(i);


      ::std::vector<TypeParam> test_vals;
      ::std::vector<TypeParam> gold_vals;

      int jmax = 0;
      for (auto it = test_range.first; it != test_range.second; ++it) {
        test_vals.push_back((*it).second);
        ++jmax;
      }

      int kmax = 0;
      for (auto it = gold_range.first; it != gold_range.second; ++it) {
        gold_vals.push_back(it->second);
        ++kmax;
      }

      EXPECT_EQ(jmax, kmax);

      ::std::sort(test_vals.begin(), test_vals.end());
      ::std::sort(gold_vals.begin(), gold_vals.end());

      same = ::std::equal(test_vals.begin(), test_vals.end(), gold_vals.begin());

//      if (!same) {
//        printf("test %d\n", i);
//        for (int j = 0; j < jmax; ++j) {
//          printf("%ld\t%ld\n", test_vals[j], gold_vals[i]);
//        }
//      }

      ASSERT_TRUE(same);
    }
}


TYPED_TEST_P(DenseHashMapTest, count_full)
{
  ::fsc::densehash_map<TypeParam, TypeParam, IsLower<TypeParam>> test(this->temp.begin(), this->temp.end(), 255, 254, 0, 1);

    for (int i = 0; i < this->max_val; ++i) {
      EXPECT_EQ(this->gold.count(i), test.count(i));
    }
}




// now register the test cases
REGISTER_TYPED_TEST_CASE_P(DenseHashMapTest, insert_partial, equal_range_partial, count_partial, insert_full, equal_range_full, count_full);


//////////////////// RUN the tests with different types.

typedef ::testing::Types<uint8_t, int16_t, int32_t,
    int64_t, uint64_t> DenseHashMapTestTypes;
INSTANTIATE_TYPED_TEST_CASE_P(Bliss, DenseHashMapTest, DenseHashMapTestTypes);









/*
 * test class holding some information.  Also, needed for the typed tests
 */
template<typename T>
class DenseHashMultimapTest : public ::testing::Test
{
    static_assert(std::is_integral<T>::value, "only supporting integral types in tests right now.");
  protected:


    ::std::unordered_multimap<T, T> gold;
    ::std::vector<std::pair<T, T>> temp;


    size_t iters = 10000;
    int max_val = 254;

    virtual void SetUp()
    { // generate some inputs


      std::default_random_engine generator;
      std::uniform_int_distribution<T> distribution(0, max_val);

      for (size_t i=0; i< iters; ++i) {
        T key = distribution(generator);
        T val = distribution(generator);
        gold.emplace(key, val);
        temp.emplace_back(::std::move(key), ::std::move(val));
      }

    }
};

// indicate this is a typed test
TYPED_TEST_CASE_P(DenseHashMultimapTest);

TYPED_TEST_P(DenseHashMultimapTest, insert_partial)
{
  bool same = false;

  ::fsc::densehash_multimap<TypeParam, TypeParam, ::fsc::TruePredicate> test(this->temp.begin(), this->temp.end(), 255, 254);

      ::std::vector<::std::pair<TypeParam, TypeParam> > test_vals = test.to_vector();
      ::std::vector<::std::pair<TypeParam, TypeParam> > gold_vals(this->gold.begin(), this->gold.end());


      EXPECT_EQ(gold_vals.size(), test_vals.size());


      ::std::sort(test_vals.begin(), test_vals.end(), [](::std::pair<TypeParam, TypeParam> const & x, ::std::pair<TypeParam, TypeParam> const &y) {
        return (x.first == y.first) ? (x.second < y.second) : (x.first < y.first);
      } );
      ::std::sort(gold_vals.begin(), gold_vals.end(), [](::std::pair<TypeParam, TypeParam> const & x, ::std::pair<TypeParam, TypeParam> const &y) {
        return (x.first == y.first) ? (x.second < y.second) : (x.first < y.first);
      } );

      same = ::std::equal(test_vals.begin(), test_vals.end(), gold_vals.begin());

//      if (!same) {
//        for (size_t i = 0; i < 100; ++i) {
//          printf("%ld->%ld\t%ld->%ld\n", test_vals[i].first, test_vals[i].second, gold_vals[i].first, gold_vals[i].second);
//        }
//        printf("\n...\n\n");
//        for (size_t i = test_vals.size() - std::min(100UL, test_vals.size()); i < test_vals.size(); ++i) {
//          printf("%ld->%ld\t%ld->%ld\n", test_vals[i].first, test_vals[i].second, gold_vals[i].first, gold_vals[i].second);
//        }
//      }
//


      EXPECT_TRUE(same);
}


TYPED_TEST_P(DenseHashMultimapTest, equal_range_partial)
{
  ::fsc::densehash_multimap<TypeParam, TypeParam, ::fsc::TruePredicate> test(this->temp.begin(), this->temp.end(), 255, 254);

  bool same = false;
    for (int i = 0; i < this->max_val; ++i) {
      auto test_range = test.equal_range(i);
      auto gold_range = this->gold.equal_range(i);


      ::std::vector<TypeParam> test_vals;
      ::std::vector<TypeParam> gold_vals;

      int jmax = 0;
      for (auto it = test_range.first; it != test_range.second; ++it) {
        test_vals.push_back((*it).second);
        ++jmax;
      }

      int kmax = 0;
      for (auto it = gold_range.first; it != gold_range.second; ++it) {
        gold_vals.push_back(it->second);
        ++kmax;
      }

      EXPECT_EQ(jmax, kmax);

      ::std::sort(test_vals.begin(), test_vals.end());
      ::std::sort(gold_vals.begin(), gold_vals.end());

      same = ::std::equal(test_vals.begin(), test_vals.end(), gold_vals.begin());

//      if (!same) {
//        printf("test %d\n", i);
//        for (int j = 0; j < jmax; ++j) {
//          printf("%ld\t%ld\n", test_vals[j], gold_vals[i]);
//        }
//      }

      EXPECT_TRUE(same);
    }
}


TYPED_TEST_P(DenseHashMultimapTest, count_partial)
{
  ::fsc::densehash_multimap<TypeParam, TypeParam, ::fsc::TruePredicate> test(this->temp.begin(), this->temp.end(), 255, 254);

    for (int i = 0; i < this->max_val; ++i) {
      EXPECT_EQ(this->gold.count(i), test.count(i));
    }
}


TYPED_TEST_P(DenseHashMultimapTest, insert_full)
{
  bool same = false;

  ::fsc::densehash_multimap<TypeParam, TypeParam, IsLower<TypeParam>> test(this->temp.begin(), this->temp.end(), 255, 254, 0, 1);

      ::std::vector<::std::pair<TypeParam, TypeParam> > test_vals = test.to_vector();
      ::std::vector<::std::pair<TypeParam, TypeParam> > gold_vals(this->gold.begin(), this->gold.end());


      ::std::sort(test_vals.begin(), test_vals.end(), [](::std::pair<TypeParam, TypeParam> const & x, ::std::pair<TypeParam, TypeParam> const &y) {
        return (x.first == y.first) ? (x.second < y.second) : (x.first < y.first);
      } );
      ::std::sort(gold_vals.begin(), gold_vals.end(), [](::std::pair<TypeParam, TypeParam> const & x, ::std::pair<TypeParam, TypeParam> const &y) {
        return (x.first == y.first) ? (x.second < y.second) : (x.first < y.first);
      } );

      same = ::std::equal(test_vals.begin(), test_vals.end(), gold_vals.begin());
//
//      if (!same) {
//        for (int i = 0; i < this->iters; ++i) {
//          printf("%ld->%ld\t%ld->%ld\n", test_vals[i].first, test_vals[i].second, gold_vals[i].first, gold_vals[i].second);
//        }
//      }

      EXPECT_TRUE(same);
}


TYPED_TEST_P(DenseHashMultimapTest, equal_range_full)
{
  ::fsc::densehash_multimap<TypeParam, TypeParam, IsLower<TypeParam>> test(this->temp.begin(), this->temp.end(), 255, 254, 0, 1);

  bool same = false;
    for (int i = 0; i < this->max_val; ++i) {
      auto test_range = test.equal_range(i);
      auto gold_range = this->gold.equal_range(i);


      ::std::vector<TypeParam> test_vals;
      ::std::vector<TypeParam> gold_vals;

      int jmax = 0;
      for (auto it = test_range.first; it != test_range.second; ++it) {
        test_vals.push_back((*it).second);
        ++jmax;
      }

      int kmax = 0;
      for (auto it = gold_range.first; it != gold_range.second; ++it) {
        gold_vals.push_back(it->second);
        ++kmax;
      }

      EXPECT_EQ(jmax, kmax);

      ::std::sort(test_vals.begin(), test_vals.end());
      ::std::sort(gold_vals.begin(), gold_vals.end());

      same = ::std::equal(test_vals.begin(), test_vals.end(), gold_vals.begin());

//      if (!same) {
//        printf("test %d\n", i);
//        for (int j = 0; j < jmax; ++j) {
//          printf("%ld\t%ld\n", test_vals[j], gold_vals[i]);
//        }
//      }

      EXPECT_TRUE(same);
    }
}


TYPED_TEST_P(DenseHashMultimapTest, count_full)
{
  ::fsc::densehash_multimap<TypeParam, TypeParam, IsLower<TypeParam>> test(this->temp.begin(), this->temp.end(), 255, 254, 0, 1);

    for (int i = 0; i < this->max_val; ++i) {
      EXPECT_EQ(this->gold.count(i), test.count(i));
    }
}




// now register the test cases
REGISTER_TYPED_TEST_CASE_P(DenseHashMultimapTest, insert_partial, equal_range_partial, count_partial, insert_full, equal_range_full, count_full);


//////////////////// RUN the tests with different types.

typedef ::testing::Types<uint8_t, int16_t, int32_t,
    int64_t, uint64_t> DenseHashMultimapTestTypes;
INSTANTIATE_TYPED_TEST_CASE_P(Bliss, DenseHashMultimapTest, DenseHashMultimapTestTypes);








