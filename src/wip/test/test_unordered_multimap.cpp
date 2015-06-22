
// include google test
#include <gtest/gtest.h>

#include <string>
#include <unordered_map>
#include <random>
#include <iterator>  // for ostream iterator.
#include <algorithm>  // for transform.

// include files to test
#include "wip/unordered_multimap.hpp"
#include "utils/logging.h"

/*
 * test class holding some information.  Also, needed for the typed tests
 */
template<typename T>
class UnorderedMultimapTest : public ::testing::Test
{
  protected:
    ::std::unordered_multimap<T, T> gold;
    ::fsc::unordered_multimap<T, T> test;

    virtual void SetUp()
    { // generate some inputs
      const int nrolls = 10000; // number of experiments

      std::default_random_engine generator;
      std::uniform_int_distribution<T> distribution(0,99);


      for (int i=0; i<nrolls; ++i) {
        T key = distribution(generator);
        T val = distribution(generator);
        test.emplace(::std::move(key), ::std::move(val));
        gold.emplace(key, val);
      }

    }
};

// indicate this is a typed test
TYPED_TEST_CASE_P(UnorderedMultimapTest);



TYPED_TEST_P(UnorderedMultimapTest, equal_range)
{
  bool same = false;
	  for (int i = 0; i < 99; ++i) {
	    auto test_range = this->test.equal_range(i);
	    auto gold_range = this->gold.equal_range(i);


	    ::std::vector<TypeParam> test_vals;
	    ::std::vector<TypeParam> gold_vals;

	    for (auto it = test_range.first; it != test_range.second; ++it) {
	      test_vals.push_back(it->second);
	    }
      for (auto it = gold_range.first; it != gold_range.second; ++it) {
        gold_vals.push_back(it->second);
      }

      ::std::sort(test_vals.begin(), test_vals.end());
      ::std::sort(gold_vals.begin(), gold_vals.end());

	    same = ::std::equal(test_vals.begin(), test_vals.end(), gold_vals.begin());

//	    if (!same) {
//	      printf("test %d\n", i);
//	      for (auto it = test_range.first; it != test_range.second; ++it) {
//	        printf("%ld->%ld\n", it->first, it->second);
//	      }
//        printf("gold %d\n", i);
//        for (auto it = gold_range.first; it != gold_range.second; ++it) {
//          printf("%ld->%ld\n", it->first, it->second);
//        }
//	    }

		  EXPECT_TRUE(same);
	  }
}


TYPED_TEST_P(UnorderedMultimapTest, count)
{
    for (int i = 0; i < 99; ++i) {
      EXPECT_EQ(this->gold.count(i), this->test.count(i));
    }
}

TYPED_TEST_P(UnorderedMultimapTest, iterator)
{
  using valType = ::std::pair<TypeParam, TypeParam>;


  ::std::vector<valType> test_vals;
  ::std::vector<valType> gold_vals;

  for (auto it = this->test.begin(), max = this->test.end(); it != max; ++it) {
    test_vals.push_back(*it);
  }
  for (auto it = this->gold.begin(), max = this->gold.end(); it != max; ++it) {
    gold_vals.push_back(*it);
  }

  ::std::sort(test_vals.begin(), test_vals.end(), [](valType const &x, valType const &y) {
    if (x.first < y.first) return true;
    else if (x.first == y.first) return (x.second < y.second);
    else return false;
  });
  ::std::sort(gold_vals.begin(), gold_vals.end(), [](valType const &x, valType const &y) {
    if (x.first < y.first) return true;
    else if (x.first == y.first) return (x.second < y.second);
    else return false;
  });

  bool same = ::std::equal(test_vals.begin(), test_vals.end(), gold_vals.begin(), [](valType const &x, valType const &y) {
    return (x.first == y.first) && (x.second == y.second);
  });

//      if (!same) {
//        printf("test %d\n", i);
//        for (auto it = test_range.first; it != test_range.second; ++it) {
//          printf("%ld->%ld\n", it->first, it->second);
//        }
//        printf("gold %d\n", i);
//        for (auto it = gold_range.first; it != gold_range.second; ++it) {
//          printf("%ld->%ld\n", it->first, it->second);
//        }
//      }

  EXPECT_TRUE(same);
}


TYPED_TEST_P(UnorderedMultimapTest, copy)
{
  using valType = ::std::pair<TypeParam, TypeParam>;


  ::std::vector<valType> test_vals(this->test.size());
  ::std::vector<valType> gold_vals(this->gold.size());

  ::std::copy(this->test.begin(), this->test.end(), test_vals.begin());
  ::std::copy(this->gold.begin(), this->gold.end(), gold_vals.begin());

  ::std::sort(test_vals.begin(), test_vals.end(), [](valType const &x, valType const &y) {
    if (x.first < y.first) return true;
    else if (x.first == y.first) return (x.second < y.second);
    else return false;
  });
  ::std::sort(gold_vals.begin(), gold_vals.end(), [](valType const &x, valType const &y) {
    if (x.first < y.first) return true;
    else if (x.first == y.first) return (x.second < y.second);
    else return false;
  });

  bool same = ::std::equal(test_vals.begin(), test_vals.end(), gold_vals.begin(), [](valType const &x, valType const &y) {
    return (x.first == y.first) && (x.second == y.second);
  });

//      if (!same) {
//        printf("test %d\n", i);
//        for (auto it = test_range.first; it != test_range.second; ++it) {
//          printf("%ld->%ld\n", it->first, it->second);
//        }
//        printf("gold %d\n", i);
//        for (auto it = gold_range.first; it != gold_range.second; ++it) {
//          printf("%ld->%ld\n", it->first, it->second);
//        }
//      }

  EXPECT_TRUE(same);
}




// now register the test cases
REGISTER_TYPED_TEST_CASE_P(UnorderedMultimapTest, equal_range, count, iterator, copy);


//////////////////// RUN the tests with different types.

typedef ::testing::Types<int8_t, int16_t, int32_t,
    int64_t, uint64_t> UnorderedMultimapTestTypes;
INSTANTIATE_TYPED_TEST_CASE_P(Bliss, UnorderedMultimapTest, UnorderedMultimapTestTypes);
