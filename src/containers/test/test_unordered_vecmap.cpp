
// include google test
#include <gtest/gtest.h>
#include "containers/unordered_vecmap.hpp"

#include <string>
#include <unordered_map>
#include <random>
#include <iterator>  // for ostream iterator.
#include <algorithm>  // for transform.

// include files to test
#include "utils/logging.h"

/*
 * test class holding some information.  Also, needed for the typed tests
 */
template<typename T>
class UnorderedCompactVecMapTest : public ::testing::Test
{
  protected:
    ::std::unordered_multimap<T, T> gold;
    ::fsc::unordered_compact_vecmap<T, T> test;

    size_t iters = 10000;

    virtual void SetUp()
    { // generate some inputs


      std::default_random_engine generator;
      std::uniform_int_distribution<T> distribution(0,99);


      for (size_t i=0; i< iters; ++i) {
        T key = distribution(generator);
        T val = distribution(generator);
        test.emplace(::std::move(key), ::std::move(val));
        gold.emplace(key, val);
      }

    }
};

// indicate this is a typed test
TYPED_TEST_CASE_P(UnorderedCompactVecMapTest);

TYPED_TEST_P(UnorderedCompactVecMapTest, insert)
{
  bool same = false;


  ::fsc::unordered_compact_vecmap<TypeParam, TypeParam> test2(this->gold.begin(), this->gold.end());

      ::std::vector<::std::pair<TypeParam, TypeParam> > test_vals(test2.begin(), test2.end());
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



TYPED_TEST_P(UnorderedCompactVecMapTest, equal_range_value_only)
{
  bool same = false;
    for (int i = 0; i < 99; ++i) {
      auto test_range = this->test.equal_range_value_only(i);
      auto gold_range = this->gold.equal_range(i);


      ::std::vector<TypeParam> test_vals;
      ::std::vector<TypeParam> gold_vals;

      for (auto it = test_range.first; it != test_range.second; ++it) {
        test_vals.push_back((*it));
      }
      for (auto it = gold_range.first; it != gold_range.second; ++it) {
        gold_vals.push_back(it->second);
      }

      ::std::sort(test_vals.begin(), test_vals.end());
      ::std::sort(gold_vals.begin(), gold_vals.end());

      same = ::std::equal(test_vals.begin(), test_vals.end(), gold_vals.begin());

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
}

TYPED_TEST_P(UnorderedCompactVecMapTest, equal_range)
{
  bool same = false;
	  for (int i = 0; i < 99; ++i) {
	    auto test_range = this->test.equal_range(i);
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


TYPED_TEST_P(UnorderedCompactVecMapTest, count)
{
    for (int i = 0; i < 99; ++i) {
      EXPECT_EQ(this->gold.count(i), this->test.count(i));
    }
}

TYPED_TEST_P(UnorderedCompactVecMapTest, iterator)
{
  using valType = ::std::pair<TypeParam, TypeParam>;


  ::std::vector<valType> test_vals;
  ::std::vector<valType> gold_vals;

  auto mx = this->test.end();
  for (auto it = this->test.begin(); it != mx; ++it) {
    test_vals.push_back(*it);
  }
  auto mx2 = this->gold.end();
  for (auto it = this->gold.begin(); it != mx2; ++it) {
    gold_vals.push_back(*it);
  }

  EXPECT_EQ(test_vals.size(), this->iters);
  EXPECT_EQ(this->test.size(), this->iters);

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


TYPED_TEST_P(UnorderedCompactVecMapTest, rand_iterator)
{
  auto b = this->test.begin();
  auto e = this->test.end();

  auto dist = e - b;
  EXPECT_EQ( dist, static_cast<ptrdiff_t>(this->iters));

  auto dist2 = ::std::distance(b, e);
  EXPECT_EQ( dist2, static_cast<ptrdiff_t>(this->iters));

  b += dist;
  dist = e - b;
  EXPECT_EQ( dist, static_cast<ptrdiff_t>(0));

  b = this->test.begin();
  for (size_t i = 0; i < this->iters / 2; ++i) {
    ++b;
  }

  dist = b - this->test.begin();
  EXPECT_EQ( dist, static_cast<ptrdiff_t>(this->iters / 2));

  auto c = this->test.begin();
  auto v = c[this->iters / 2];
  EXPECT_EQ(v, *b);

  auto d = this->test.begin() + this->iters/2;
  EXPECT_EQ(*d, *b);


}



TYPED_TEST_P(UnorderedCompactVecMapTest, copy)
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
REGISTER_TYPED_TEST_CASE_P(UnorderedCompactVecMapTest, insert, equal_range_value_only, equal_range, count, iterator, rand_iterator, copy);


//////////////////// RUN the tests with different types.

typedef ::testing::Types<int8_t, int16_t, int32_t,
    int64_t, uint64_t> UnorderedCompactVecMapTestTypes;
INSTANTIATE_TYPED_TEST_CASE_P(Bliss, UnorderedCompactVecMapTest, UnorderedCompactVecMapTestTypes);






/*
 * test class holding some information.  Also, needed for the typed tests
 */
template<typename T>
class UnorderedVecMapTest : public ::testing::Test
{
  protected:
    ::std::unordered_multimap<T, T> gold;
    ::fsc::unordered_vecmap<T, T> test;

    size_t iters = 10000;

    virtual void SetUp()
    { // generate some inputs


      std::default_random_engine generator;
      std::uniform_int_distribution<T> distribution(0,99);


      for (size_t i=0; i< iters; ++i) {
        T key = distribution(generator);
        T val = distribution(generator);
        test.emplace(::std::move(key), ::std::move(val));
        gold.emplace(key, val);
      }

    }
};

// indicate this is a typed test
TYPED_TEST_CASE_P(UnorderedVecMapTest);

TYPED_TEST_P(UnorderedVecMapTest, insert)
{
  bool same = false;


  ::fsc::unordered_vecmap<TypeParam, TypeParam> test2(this->gold.begin(), this->gold.end());

      ::std::vector<::std::pair<TypeParam, TypeParam> > test_vals(test2.begin(), test2.end());
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
//        for (size_t i = 0; i < this->iters; ++i) {
//          printf("%ld->%ld\t%ld->%ld\n", test_vals[i].first, test_vals[i].second, gold_vals[i].first, gold_vals[i].second);
//        }
//      }

      EXPECT_TRUE(same);
}




TYPED_TEST_P(UnorderedVecMapTest, equal_range)
{
  bool same = false;
    for (int i = 0; i < 99; ++i) {
      auto test_range = this->test.equal_range(i);
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


TYPED_TEST_P(UnorderedVecMapTest, count)
{
    for (int i = 0; i < 99; ++i) {
      EXPECT_EQ(this->gold.count(i), this->test.count(i));
    }
}

TYPED_TEST_P(UnorderedVecMapTest, iterator)
{
  using valType = ::std::pair<TypeParam, TypeParam>;


  ::std::vector<valType> test_vals;
  ::std::vector<valType> gold_vals;

  auto mx = this->test.end();
  for (auto it = this->test.begin(); it != mx; ++it) {
    test_vals.push_back(*it);
  }
  auto mx2 = this->gold.end();
  for (auto it = this->gold.begin(); it != mx2; ++it) {
    gold_vals.push_back(*it);
  }

  EXPECT_EQ(test_vals.size(), this->iters);
  EXPECT_EQ(this->test.size(), this->iters);

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


TYPED_TEST_P(UnorderedVecMapTest, rand_iterator)
{
  auto b = this->test.begin();
  auto e = this->test.end();

  auto dist = e - b;
  EXPECT_EQ( dist, static_cast<ptrdiff_t>(this->iters));

  auto dist2 = ::std::distance(b, e);
  EXPECT_EQ( dist2, static_cast<ptrdiff_t>(this->iters));

  b += dist;
  dist = e - b;
  EXPECT_EQ( dist, static_cast<ptrdiff_t>(0));

  b = this->test.begin();
  for (size_t i = 0; i < this->iters / 2; ++i) {
    ++b;
  }

  dist = b - this->test.begin();
  EXPECT_EQ( dist, static_cast<ptrdiff_t>(this->iters / 2));

  auto c = this->test.begin();
  auto v = c[this->iters / 2];
  EXPECT_EQ(v, *b);

  auto d = this->test.begin() + this->iters/2;
  EXPECT_EQ(*d, *b);


}



TYPED_TEST_P(UnorderedVecMapTest, copy)
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
REGISTER_TYPED_TEST_CASE_P(UnorderedVecMapTest, insert, equal_range, count, iterator, rand_iterator, copy);


//////////////////// RUN the tests with different types.

typedef ::testing::Types<int8_t, int16_t, int32_t,
    int64_t, uint64_t> UnorderedVecMapTestTypes;
INSTANTIATE_TYPED_TEST_CASE_P(Bliss, UnorderedVecMapTest, UnorderedVecMapTestTypes);

