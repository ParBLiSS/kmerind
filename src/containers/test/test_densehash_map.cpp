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
#include <unordered_set>
#include <random>
#include <iterator>  // for ostream iterator.
#include <algorithm>  // for transform.
#include <cstdint>  // uint32_t
#include <utility>  // pair
#include <vector>

// include files to test
#include "utils/logging.h"
#include "utils/transform_utils.hpp"
#include "utils/filter_utils.hpp"


#include "common/kmer.hpp"
#include "common/alphabets.hpp"
#include "common/kmer_transform.hpp"
#include "index/kmer_hash.hpp"
#include "iterators/transform_iterator.hpp"
#include "containers/fsc_container_utils.hpp"

  template <typename Key>
  struct full_special_keys {
	static_assert(::std::is_integral<Key>::value && !::std::is_signed<Key>::value, "example imple only supports unsigned int");

	inline Key generate(uint8_t id = 0) {
		return ::std::numeric_limits<Key>::max() - id;
	}

	inline Key invert(Key const &x) {
		return static_cast<Key>(~x);
	}

	inline Key get_splitter() {
		return static_cast<Key>(~(::std::numeric_limits<Key>::max() >> 2));
	}

	static constexpr bool need_to_split = true;
  };

/*
 * test class holding some information.  Also, needed for the typed tests
 */
template<typename T>
class DenseHashMapPartialTest : public ::testing::Test
{
    static_assert(std::is_integral<T>::value, "only supporting integral types in tests right now.");
  protected:


    ::std::unordered_map<T, T> gold;
    ::std::vector<std::pair<T, T>> temp;


    size_t iters = 100000;
    T min_val = 2;
    T max_val = ::std::numeric_limits<T>::max() - 2;

    virtual void SetUp()
    { // generate some inputs


      std::default_random_engine generator;
      std::uniform_int_distribution<T> distribution(min_val, max_val);

      for (size_t i=0; i< iters; ++i) {
        T key = distribution(generator);
        T val = distribution(generator);
        gold.emplace(key, val);
        temp.emplace_back(::std::move(key), ::std::move(val));
      }

    }
};

// indicate this is a typed test
TYPED_TEST_CASE_P(DenseHashMapPartialTest);

TYPED_TEST_P(DenseHashMapPartialTest, insert_partial)
{
  bool same = false;

  using MAP = ::fsc::densehash_map<TypeParam, TypeParam>;


   MAP test(this->temp.begin(), this->temp.end());

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


TYPED_TEST_P(DenseHashMapPartialTest, equal_range_partial)
{
	  using MAP = ::fsc::densehash_map<TypeParam, TypeParam>;


	   MAP test(this->temp.begin(), this->temp.end());


	   ::std::vector<::std::pair<TypeParam, TypeParam> > unique(this->temp.begin(), this->temp.end());
	   std::sort(unique.begin(), unique.end(), [](::std::pair<TypeParam, TypeParam> const & x, ::std::pair<TypeParam, TypeParam> const & y){
		   return x.first < y.first;
	   });
	   auto newend = std::unique(unique.begin(), unique.end(), [](::std::pair<TypeParam, TypeParam> const & x, ::std::pair<TypeParam, TypeParam> const & y){
		   return x.first == y.first;
	   });
	   unique.erase(newend, unique.end());


  bool same = false;
	  for (auto i : unique) {
	    auto test_range = test.equal_range(i.first);
	    auto gold_range = this->gold.equal_range(i.first);


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


TYPED_TEST_P(DenseHashMapPartialTest, count_partial)
{
	  using MAP = ::fsc::densehash_map<TypeParam, TypeParam>;


	   MAP test(this->temp.begin(), this->temp.end());

	   ::std::vector<::std::pair<TypeParam, TypeParam> > unique(this->temp.begin(), this->temp.end());
	   std::sort(unique.begin(), unique.end(), [](::std::pair<TypeParam, TypeParam> const & x, ::std::pair<TypeParam, TypeParam> const & y){
		   return x.first < y.first;
	   });
	   auto newend = std::unique(unique.begin(), unique.end(), [](::std::pair<TypeParam, TypeParam> const & x, ::std::pair<TypeParam, TypeParam> const & y){
		   return x.first == y.first;
	   });
	   unique.erase(newend, unique.end());



		  for (auto i : unique) {
      EXPECT_EQ(this->gold.count(i.first), test.count(i.first));
    }
}

// now register the test cases
REGISTER_TYPED_TEST_CASE_P(DenseHashMapPartialTest, insert_partial, equal_range_partial, count_partial);


//////////////////// RUN the tests with different types.

typedef ::testing::Types<uint8_t, uint16_t,
    uint32_t, uint64_t> DenseHashMapPartialTestTypes;
INSTANTIATE_TYPED_TEST_CASE_P(Bliss, DenseHashMapPartialTest, DenseHashMapPartialTestTypes);



/*
 * test class holding some information.  Also, needed for the typed tests
 */
template<typename T>
class DenseHashMapFullTest : public ::testing::Test
{
    static_assert(std::is_integral<T>::value, "only supporting integral types in tests right now.");
  protected:


    ::std::unordered_map<T, T> gold;
    ::std::vector<std::pair<T, T>> temp;


    size_t iters = 100000;
    T min_val = 0;
    T max_val = 255;

    virtual void SetUp()
    { // generate some inputs


      std::default_random_engine generator;
      std::uniform_int_distribution<T> distribution(min_val, max_val);

      for (size_t i=0; i< iters; ++i) {
        T key = distribution(generator);
        T val = distribution(generator);
        gold.emplace(key, val);
        temp.emplace_back(::std::move(key), ::std::move(val));
      }

      gold.emplace(0, 0);
      gold.emplace(1, 1);
      gold.emplace(max_val - 1, 0);
      gold.emplace(max_val, 0);

      temp.emplace_back(0, 0);
      temp.emplace_back(1, 1);
      temp.emplace_back(max_val - 1, 0);
      temp.emplace_back(max_val, 0);


    }
};

// indicate this is a typed test
TYPED_TEST_CASE_P(DenseHashMapFullTest);


TYPED_TEST_P(DenseHashMapFullTest, insert_full)
{
  bool same = false;

  using MAP = ::fsc::densehash_map<TypeParam, TypeParam, full_special_keys<TypeParam> >;


   MAP test(this->temp.begin(), this->temp.end());


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


TYPED_TEST_P(DenseHashMapFullTest, equal_range_full)
{
	  using MAP = ::fsc::densehash_map<TypeParam, TypeParam, full_special_keys<TypeParam> >;


	   MAP test(this->temp.begin(), this->temp.end());

	   ::std::vector<::std::pair<TypeParam, TypeParam> > unique(this->temp.begin(), this->temp.end());
	   std::sort(unique.begin(), unique.end(), [](::std::pair<TypeParam, TypeParam> const & x, ::std::pair<TypeParam, TypeParam> const & y){
		   return x.first < y.first;
	   });
	   auto newend = std::unique(unique.begin(), unique.end(), [](::std::pair<TypeParam, TypeParam> const & x, ::std::pair<TypeParam, TypeParam> const & y){
		   return x.first == y.first;
	   });
	   unique.erase(newend, unique.end());


  bool same = false;
  for (auto i : unique) {
      auto test_range = test.equal_range(i.first);
      auto gold_range = this->gold.equal_range(i.first);


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


TYPED_TEST_P(DenseHashMapFullTest, count_full)
{
	  using MAP = ::fsc::densehash_map<TypeParam, TypeParam, full_special_keys<TypeParam> >;


	   MAP test(this->temp.begin(), this->temp.end());

	   ::std::vector<::std::pair<TypeParam, TypeParam> > unique(this->temp.begin(), this->temp.end());
	   std::sort(unique.begin(), unique.end(), [](::std::pair<TypeParam, TypeParam> const & x, ::std::pair<TypeParam, TypeParam> const & y){
		   return x.first < y.first;
	   });
	   auto newend = std::unique(unique.begin(), unique.end(), [](::std::pair<TypeParam, TypeParam> const & x, ::std::pair<TypeParam, TypeParam> const & y){
		   return x.first == y.first;
	   });
	   unique.erase(newend, unique.end());



	   for (auto i : unique) {
      EXPECT_EQ(this->gold.count(i.first), test.count(i.first));
    }
}




// now register the test cases
REGISTER_TYPED_TEST_CASE_P(DenseHashMapFullTest, insert_full, equal_range_full, count_full);


//////////////////// RUN the tests with different types.

typedef ::testing::Types<uint8_t, uint16_t,
	    uint32_t, uint64_t> DenseHashMapFullTestTypes;
INSTANTIATE_TYPED_TEST_CASE_P(Bliss, DenseHashMapFullTest, DenseHashMapFullTestTypes);









/*
 * test class holding some information.  Also, needed for the typed tests
 */
template<typename T>
class DenseHashMultimapPartialTest : public ::testing::Test
{
    static_assert(std::is_integral<T>::value, "only supporting integral types in tests right now.");
  protected:


    ::std::unordered_multimap<T, T> gold;
    ::std::vector<std::pair<T, T>> temp;


    size_t iters = 100000;
    T min_val = 2;
    T max_val = ::std::numeric_limits<T>::max() - 2;

    virtual void SetUp()
    { // generate some inputs


      std::default_random_engine generator;
      std::uniform_int_distribution<T> distribution(min_val, max_val);

      for (size_t i=0; i< iters; ++i) {
        T key = distribution(generator);
        T val = distribution(generator);
        gold.emplace(key, val);
        temp.emplace_back(::std::move(key), ::std::move(val));
      }

    }
};

// indicate this is a typed test
TYPED_TEST_CASE_P(DenseHashMultimapPartialTest);

TYPED_TEST_P(DenseHashMultimapPartialTest, insert_partial)
{
  bool same = false;

  using MAP = ::fsc::densehash_multimap<TypeParam, TypeParam>;


   MAP test(this->temp.begin(), this->temp.end());


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


TYPED_TEST_P(DenseHashMultimapPartialTest, equal_range_partial)
{
	  using MAP = ::fsc::densehash_multimap<TypeParam, TypeParam>;


	   MAP test(this->temp.begin(), this->temp.end());

	   ::std::vector<::std::pair<TypeParam, TypeParam> > unique(this->temp.begin(), this->temp.end());
	   std::sort(unique.begin(), unique.end(), [](::std::pair<TypeParam, TypeParam> const & x, ::std::pair<TypeParam, TypeParam> const & y){
		   return x.first < y.first;
	   });
	   auto newend = std::unique(unique.begin(), unique.end(), [](::std::pair<TypeParam, TypeParam> const & x, ::std::pair<TypeParam, TypeParam> const & y){
		   return x.first == y.first;
	   });
	   unique.erase(newend, unique.end());


  bool same = false;
  for (auto i : unique) {
	  auto test_range = test.equal_range(i.first);
      auto gold_range = this->gold.equal_range(i.first);


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


TYPED_TEST_P(DenseHashMultimapPartialTest, count_partial)
{
	  using MAP = ::fsc::densehash_multimap<TypeParam, TypeParam>;


	   MAP test(this->temp.begin(), this->temp.end());

	   ::std::vector<::std::pair<TypeParam, TypeParam> > unique(this->temp.begin(), this->temp.end());
	   std::sort(unique.begin(), unique.end(), [](::std::pair<TypeParam, TypeParam> const & x, ::std::pair<TypeParam, TypeParam> const & y){
		   return x.first < y.first;
	   });
	   auto newend = std::unique(unique.begin(), unique.end(), [](::std::pair<TypeParam, TypeParam> const & x, ::std::pair<TypeParam, TypeParam> const & y){
		   return x.first == y.first;
	   });
	   unique.erase(newend, unique.end());


	   for (auto i : unique) {
      EXPECT_EQ(this->gold.count(i.first), test.count(i.first));
    }
}




// now register the test cases
REGISTER_TYPED_TEST_CASE_P(DenseHashMultimapPartialTest, insert_partial, equal_range_partial, count_partial);


//////////////////// RUN the tests with different types.

typedef ::testing::Types<uint8_t, uint16_t,
	    uint32_t,  uint64_t> DenseHashMultimapPartialTestTypes;
INSTANTIATE_TYPED_TEST_CASE_P(Bliss, DenseHashMultimapPartialTest, DenseHashMultimapPartialTestTypes);






/*
 * test class holding some information.  Also, needed for the typed tests
 */
template<typename T>
class DenseHashMultimapFullTest : public ::testing::Test
{
    static_assert(std::is_integral<T>::value, "only supporting integral types in tests right now.");
  protected:


    ::std::unordered_multimap<T, T> gold;
    ::std::vector<std::pair<T, T>> temp;


    size_t iters = 100000;
    T min_val = 0;
    T max_val = ::std::numeric_limits<T>::max();

    virtual void SetUp()
    { // generate some inputs


      std::default_random_engine generator;
      std::uniform_int_distribution<T> distribution(min_val, max_val);

      for (size_t i=0; i< iters; ++i) {
        T key = distribution(generator);
        T val = distribution(generator);
        gold.emplace(key, val);
        temp.emplace_back(::std::move(key), ::std::move(val));
      }

// workorder      607874632


      gold.emplace(0, 0);
      gold.emplace(1, 1);
//      gold.emplace(max_val - 1, 0);
      gold.emplace(max_val, 0);

      temp.emplace_back(0, 0);
      temp.emplace_back(1, 1);
//      temp.emplace_back(max_val - 1, 0);
      temp.emplace_back(max_val, 0);

    }
};

// indicate this is a typed test
TYPED_TEST_CASE_P(DenseHashMultimapFullTest);


TYPED_TEST_P(DenseHashMultimapFullTest, insert_full)
{
  bool same = false;

  using MAP = ::fsc::densehash_multimap<TypeParam, TypeParam, full_special_keys<TypeParam> >;


   MAP test(this->temp.begin(), this->temp.end());


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


TYPED_TEST_P(DenseHashMultimapFullTest, equal_range_full)
{
	  using MAP = ::fsc::densehash_multimap<TypeParam, TypeParam, full_special_keys<TypeParam> >;
	   MAP test(this->temp.begin(), this->temp.end());

	   ::std::vector<::std::pair<TypeParam, TypeParam> > unique(this->temp.begin(), this->temp.end());
	   std::sort(unique.begin(), unique.end(), [](::std::pair<TypeParam, TypeParam> const & x, ::std::pair<TypeParam, TypeParam> const & y){
		   return x.first < y.first;
	   });
	   auto newend = std::unique(unique.begin(), unique.end(), [](::std::pair<TypeParam, TypeParam> const & x, ::std::pair<TypeParam, TypeParam> const & y){
		   return x.first == y.first;
	   });
	   unique.erase(newend, unique.end());


  bool same = false;
  for (auto i : unique) {
      auto test_range = test.equal_range(i.first);
      auto gold_range = this->gold.equal_range(i.first);


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


TYPED_TEST_P(DenseHashMultimapFullTest, count_full)
{
	  using MAP = ::fsc::densehash_multimap<TypeParam, TypeParam, full_special_keys<TypeParam> >;


	   MAP test(this->temp.begin(), this->temp.end());

	   ::std::vector<::std::pair<TypeParam, TypeParam> > unique(this->temp.begin(), this->temp.end());
	   std::sort(unique.begin(), unique.end(), [](::std::pair<TypeParam, TypeParam> const & x, ::std::pair<TypeParam, TypeParam> const & y){
		   return x.first < y.first;
	   });
	   auto newend = std::unique(unique.begin(), unique.end(), [](::std::pair<TypeParam, TypeParam> const & x, ::std::pair<TypeParam, TypeParam> const & y){
		   return x.first == y.first;
	   });
	   unique.erase(newend, unique.end());


   for (auto i : unique) {
      EXPECT_EQ(this->gold.count(i.first), test.count(i.first));
    }
}




// now register the test cases
REGISTER_TYPED_TEST_CASE_P(DenseHashMultimapFullTest, insert_full, equal_range_full, count_full);


//////////////////// RUN the tests with different types.

typedef ::testing::Types<uint8_t, uint16_t,
	    uint32_t,  uint64_t> DenseHashMultimapFullTestTypes;
INSTANTIATE_TYPED_TEST_CASE_P(Bliss, DenseHashMultimapFullTest, DenseHashMultimapFullTestTypes);



//TODO change data above to test split case.


/*
 * test class holding some information.  Also, needed for the typed tests
 */
template<typename T>
class DenseHashKmerMultimapTest : public ::testing::Test
{
  protected:


    ::std::vector<std::pair<T, uint32_t> > temp;

    using ALPHA = typename T::KmerAlphabet;
    using CANONICAL_ITER = ::bliss::iterator::transform_iterator<
    		typename ::std::vector<std::pair<T, uint32_t> >::iterator,
    		::bliss::kmer::transform::lex_less<T> >;

    size_t iters = 100000;
    int min_val = 2;
    int max_val = 253;

    virtual void SetUp()
    { // generate some inputs


      std::default_random_engine generator1, generator2;
      std::uniform_int_distribution<int> distribution1(min_val, max_val);
      std::uniform_int_distribution<int> distribution2(0, 3);

      for (size_t i=0; i< iters; ++i) {
        T key;
        for (size_t j = 0; j < T::size; ++j) {
        	key.nextFromChar(ALPHA::FROM_INDEX[distribution2(generator2)]);
        }

        uint32_t val = distribution1(generator1);
        temp.emplace_back(::std::move(key), ::std::move(val));
      }

    }

    template <bool canonical, typename Less, typename MAP, typename Kmer, typename Hash, typename Equal>
    void map_insert(MAP & test,
                    ::std::unordered_map<Kmer, uint32_t, Hash, Equal> & gold,
                     ::std::vector<std::pair<Kmer, uint32_t> > & entries) {

      test.clear();
      gold.clear();
      entries.clear();

      if (canonical) {
        entries.insert(entries.end(),
                     CANONICAL_ITER(this->temp.begin(), ::bliss::kmer::transform::lex_less<Kmer>()),
                     CANONICAL_ITER(this->temp.end(), ::bliss::kmer::transform::lex_less<Kmer>()));
        test.insert(entries.begin(), entries.end());
        gold.insert(entries.begin(), entries.end());
        printf("canonical insert.  sizes input %lu, test %lu, gold %lu\n", entries.size(), test.size(), gold.size());
      } else {
        entries.insert(entries.end(), this->temp.begin(), this->temp.end());
        test.insert(entries.begin(), entries.end());
        gold.insert(entries.begin(), entries.end());
        printf("raw insert.  sizes input %lu, test %lu, gold %lu\n", entries.size(), test.size(), gold.size());
      }

      // check unique items in list.
      std::stable_sort(entries.begin(), entries.end(), Less());
      auto new_end = std::unique(entries.begin(), entries.end(), Equal());
      entries.erase(new_end, entries.end());

      ASSERT_EQ(gold.size(), entries.size());

      ASSERT_EQ(test.size(), gold.size());

    }


    template <bool canonical, typename Less, typename MAP, typename Kmer, typename Hash, typename Equal>
    void multimap_insert(MAP & test,
                    ::std::unordered_multimap<Kmer, uint32_t, Hash, Equal> & gold,
                     ::std::vector<std::pair<Kmer, uint32_t> > & entries) {
        test.clear();
        gold.clear();
        entries.clear();

      if (canonical) {
        entries.insert(entries.end(),
                     CANONICAL_ITER(this->temp.begin(), ::bliss::kmer::transform::lex_less<Kmer>()),
                     CANONICAL_ITER(this->temp.end(), ::bliss::kmer::transform::lex_less<Kmer>()));
        test.insert(entries.begin(), entries.end());
        gold.insert(entries.begin(), entries.end());
        printf("canonical insert.  sizes input %lu, test %lu, gold %lu\n", entries.size(), test.size(), gold.size());
      } else {
        entries.insert(entries.end(), this->temp.begin(), this->temp.end());
        test.insert(entries.begin(), entries.end());
        gold.insert(entries.begin(), entries.end());
        printf("raw insert.  sizes input %lu, test %lu, gold %lu\n", entries.size(), test.size(), gold.size());
      }

      ASSERT_EQ(gold.size(), entries.size());
      ASSERT_EQ(test.size(), entries.size());


      // check unique items in list.
      std::stable_sort(entries.begin(), entries.end(), Less());
      auto new_end = std::unique(entries.begin(), entries.end(), Equal());
      entries.erase(new_end, entries.end());

      printf("number of unique entries is %lu\n", entries.size());


      ASSERT_EQ(test.size(), gold.size());

    }

    template <typename Kmer = T, bool canonical,
    		template <typename> class Transform = ::bliss::transform::identity,
			typename Hash, typename Equal>
    ::fsc::densehash_map<Kmer, uint32_t, ::bliss::kmer::hash::sparsehash::special_keys<Kmer, canonical>, Transform, Hash, Equal>
    make_kmer_map() {
//    	::bliss::kmer::hash::sparsehash::special_keys<Kmer, canonical> specials;
//    	if (canonical) std::cout << "CANONICAL ";
//    	Kmer t;
//    	std::cout << "keys:\t" << std::hex;
//    	for (int i = Kmer::nWords - 1; i >= 0; --i) {
//    		t = specials.generate(0);
//    		std::cout << t.getData()[i] << " ";
//    	}
//    	std::cout << std::endl << "\t" << std::hex;
//    	for (int i = Kmer::nWords - 1; i >= 0; --i) {
//    		t = specials.generate(1);
//    		std::cout << t.getData()[i] << " ";
//    	}
//    	std::cout << std::endl;
//    	if (!canonical) {
//    		std::cout << "\t" << std::hex;
//			for (int i = Kmer::nWords - 1; i >= 0; --i) {
//				t = specials.invert(specials.generate(0));
//				std::cout << t.getData()[i] << " ";
//			}
//			std::cout << std::endl << "\t" << std::hex;
//			for (int i = Kmer::nWords - 1; i >= 0; --i) {
//				t = specials.invert(specials.generate(1));
//				std::cout << t.getData()[i] << " ";
//			}
//			std::cout << std::endl;
//    	}

		return ::fsc::densehash_map<Kmer, uint32_t, ::bliss::kmer::hash::sparsehash::special_keys<Kmer, canonical>, Transform, Hash, Equal >();
    }


    template <typename Kmer = T, bool canonical,
    		template <typename> class Transform = ::bliss::transform::identity,
			typename Hash, typename Equal>
    ::fsc::densehash_multimap<Kmer, uint32_t, ::bliss::kmer::hash::sparsehash::special_keys<Kmer, canonical>, Transform, Hash, Equal>
    make_kmer_multimap() {
//    	::bliss::kmer::hash::sparsehash::special_keys<Kmer, canonical> specials;
//    	if (canonical) std::cout << "CANONICAL ";
//    	Kmer t;
//    	std::cout << "keys:\t" << std::hex;
//    	for (int i = Kmer::nWords - 1; i >= 0; --i) {
//    		t = specials.generate(0);
//    		std::cout << t.getData()[i] << " ";
//    	}
//    	std::cout << std::endl << "\t" << std::hex;
//    	for (int i = Kmer::nWords - 1; i >= 0; --i) {
//    		t = specials.generate(1);
//    		std::cout << t.getData()[i] << " ";
//    	}
//    	std::cout << std::endl;
//    	if (!canonical) {
//    		std::cout << "\t" << std::hex;
//			for (int i = Kmer::nWords - 1; i >= 0; --i) {
//				t = specials.invert(specials.generate(0));
//				std::cout << t.getData()[i] << " ";
//			}
//			std::cout << std::endl << "\t" << std::hex;
//			for (int i = Kmer::nWords - 1; i >= 0; --i) {
//				t = specials.invert(specials.generate(1));
//				std::cout << t.getData()[i] << " ";
//			}
//			std::cout << std::endl;
//    	}

    	return ::fsc::densehash_multimap<Kmer, uint32_t, ::bliss::kmer::hash::sparsehash::special_keys<Kmer, canonical>,  Transform, Hash, Equal >();
    }


    template <typename Kmer = T, bool canonical = false,
    		template <typename> class Transform = ::bliss::transform::identity,
			template <typename> class Hash = std::hash,
			template <typename> class Equal = std::equal_to,
			template <typename> class Less = std::less
			>
    void test_map_insert() {

    	using THash = ::fsc::TransformedHash<Kmer, Hash, Transform>;
    	using TLess = ::fsc::TransformedComparator<Kmer, Less, Transform>;
    	using Equal1 = ::fsc::TransformedComparator<Kmer, Equal, Transform>;
    	using Equal2 = ::fsc::sparsehash::compare<Kmer, Equal, Transform>;


    	auto test = make_kmer_map<Kmer, canonical, Transform, THash, Equal2>();
		::std::unordered_map<Kmer, uint32_t, THash, Equal1> gold;
		::std::vector<std::pair<Kmer, uint32_t> > entries;

		this->map_insert<canonical, TLess>(test, gold, entries);


		::std::vector<::std::pair<Kmer, uint32_t> > test_vals = test.to_vector();
		::std::vector<::std::pair<Kmer, uint32_t> > gold_vals(gold.begin(), gold.end());


		ASSERT_EQ(gold_vals.size(), test_vals.size());


		::std::sort(test_vals.begin(), test_vals.end(), [](::std::pair<Kmer, uint32_t> const & x, ::std::pair<Kmer, uint32_t> const &y) {
			return (x.first == y.first) ? (x.second < y.second) : (x.first < y.first);
		} );
		::std::sort(gold_vals.begin(), gold_vals.end(), [](::std::pair<Kmer, uint32_t> const & x, ::std::pair<Kmer, uint32_t> const &y) {
			return (x.first == y.first) ? (x.second < y.second) : (x.first < y.first);
		} );

		bool same = ::std::equal(test_vals.begin(), test_vals.end(), gold_vals.begin());

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

		ASSERT_TRUE(same);
    }


    template <typename Kmer = T, bool canonical = false,
    		template <typename> class Transform = ::bliss::transform::identity,
			template <typename> class Hash = std::hash,
			template <typename> class Equal = std::equal_to,
			template <typename> class Less = std::less
			>
    void test_map_equal_range() {

    	using THash = ::fsc::TransformedHash<Kmer, Hash, Transform>;
    	using TLess = ::fsc::TransformedComparator<Kmer, Less, Transform>;
    	using Equal1 = ::fsc::TransformedComparator<Kmer, Equal, Transform>;
    	using Equal2 = ::fsc::sparsehash::compare<Kmer, Equal, Transform>;


    	auto test = make_kmer_map<Kmer, canonical, Transform, THash, Equal2>();
		::std::unordered_map<Kmer, uint32_t, THash, Equal1> gold;
		::std::vector<std::pair<Kmer, uint32_t> > entries;

		this->map_insert<canonical, TLess>(test, gold, entries);


		// get list of unique k-mers
		std::vector<Kmer> keys = test.keys();

		// assert that there are no duplicates
		std::unordered_set<Kmer, THash, Equal1> testset(keys.begin(), keys.end());
		ASSERT_EQ(testset.size(), keys.size());

		// assert that all entries are unique
		ASSERT_EQ(gold.size(), keys.size());
		ASSERT_EQ(entries.size(), keys.size());

    	// check one by one
		bool same = false;
    	for (auto i : entries) {
			auto test_range = test.equal_range(i.first);
			auto gold_range = gold.equal_range(i.first);

			::std::vector<uint32_t> test_vals;
			::std::vector<uint32_t> gold_vals;

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

			ASSERT_EQ(jmax, kmax);

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


    template <typename Kmer = T, bool canonical = false,
    		template <typename> class Transform = ::bliss::transform::identity,
			template <typename> class Hash = std::hash,
			template <typename> class Equal = std::equal_to,
			template <typename> class Less = std::less
			>
    void test_map_count() {


    	using THash = ::fsc::TransformedHash<Kmer, Hash, Transform>;
    	using TLess = ::fsc::TransformedComparator<Kmer, Less, Transform>;
    	using Equal1 = ::fsc::TransformedComparator<Kmer, Equal, Transform>;
    	using Equal2 = ::fsc::sparsehash::compare<Kmer, Equal, Transform>;

    	auto test = make_kmer_map<Kmer, canonical, Transform, THash, Equal2>();
		::std::unordered_map<Kmer, uint32_t, THash, Equal1> gold;
		::std::vector<std::pair<Kmer, uint32_t> > entries;

		this->map_insert<canonical, TLess>(test, gold, entries);

		// get list of unique k-mers
		std::vector<Kmer> keys = test.keys();

		// assert that there are no duplicates
		std::unordered_set<Kmer, THash, Equal1> testset(keys.begin(), keys.end());
		ASSERT_EQ(testset.size(), keys.size());

		// assert that all entries are unique
		ASSERT_EQ(gold.size(), keys.size());
		ASSERT_EQ(entries.size(), keys.size());

    	// check one by one
    	for (auto i : entries) {
    		ASSERT_EQ(test.count(i.first), gold.count(i.first));
    	}
    }

    template <typename Kmer = T, bool canonical = false,
    		template <typename> class Transform = ::bliss::transform::identity,
			template <typename> class Hash = std::hash,
			template <typename> class Equal = std::equal_to,
			template <typename> class Less = std::less
			>
    void test_multimap_insert() {


    	using THash = ::fsc::TransformedHash<Kmer, Hash, Transform>;
    	using TLess = ::fsc::TransformedComparator<Kmer, Less, Transform>;
    	using Equal1 = ::fsc::TransformedComparator<Kmer, Equal, Transform>;
    	using Equal2 = ::fsc::sparsehash::compare<Kmer, Equal, Transform>;


    	auto test = make_kmer_multimap<Kmer, canonical, Transform, THash, Equal2>();
		::std::unordered_multimap<Kmer, uint32_t, THash, Equal1> gold;
		::std::vector<std::pair<Kmer, uint32_t> > entries;

		this->multimap_insert<canonical, TLess>(test, gold, entries);


		::std::vector<::std::pair<Kmer, uint32_t> > test_vals = test.to_vector();
		::std::vector<::std::pair<Kmer, uint32_t> > gold_vals(gold.begin(), gold.end());


		ASSERT_EQ(gold_vals.size(), test_vals.size());


		::std::sort(test_vals.begin(), test_vals.end(), [](::std::pair<Kmer, uint32_t> const & x, ::std::pair<Kmer, uint32_t> const &y) {
			return (x.first == y.first) ? (x.second < y.second) : (x.first < y.first);
		} );
		::std::sort(gold_vals.begin(), gold_vals.end(), [](::std::pair<Kmer, uint32_t> const & x, ::std::pair<Kmer, uint32_t> const &y) {
			return (x.first == y.first) ? (x.second < y.second) : (x.first < y.first);
		} );

		bool same = ::std::equal(test_vals.begin(), test_vals.end(), gold_vals.begin());

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

		ASSERT_TRUE(same);
    }


    template <typename Kmer = T, bool canonical = false,
    		template <typename> class Transform = ::bliss::transform::identity,
			template <typename> class Hash = std::hash,
			template <typename> class Equal = std::equal_to,
			template <typename> class Less = std::less
			>
    void test_multimap_equal_range() {


    	using THash = ::fsc::TransformedHash<Kmer, Hash, Transform>;
    	using TLess = ::fsc::TransformedComparator<Kmer, Less, Transform>;
    	using Equal1 = ::fsc::TransformedComparator<Kmer, Equal, Transform>;
    	using Equal2 = ::fsc::sparsehash::compare<Kmer, Equal, Transform>;


    	auto test = make_kmer_multimap<Kmer, canonical, Transform, THash, Equal2>();
		::std::unordered_multimap<Kmer, uint32_t, THash, Equal1> gold;
		::std::vector<std::pair<Kmer, uint32_t> > entries;

		this->multimap_insert<canonical, TLess>(test, gold, entries);


		// get list of unique k-mers
		std::vector<Kmer> keys = test.keys();

		// assert that there are no duplicates
		std::unordered_set<Kmer, THash, Equal1> unique_keys(keys.begin(), keys.end());
		ASSERT_EQ(unique_keys.size(), keys.size());

		std::unordered_map<Kmer, uint32_t, THash, Equal1> unique_entries(gold.begin(), gold.end());
		ASSERT_EQ(unique_entries.size(), unique_keys.size());
		ASSERT_EQ(entries.size(), unique_keys.size());


    	// check one by one
		bool same = false;
    	for (auto i : entries) {
			auto test_range = test.equal_range(i.first);
			auto gold_range = gold.equal_range(i.first);

			::std::vector<uint32_t> test_vals;
			::std::vector<uint32_t> gold_vals;

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

			ASSERT_EQ(jmax, kmax);

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


    template <typename Kmer = T, bool canonical = false,
    		template <typename> class Transform = ::bliss::transform::identity,
			template <typename> class Hash = std::hash,
			template <typename> class Equal = std::equal_to,
			template <typename> class Less = std::less
			>
    void test_multimap_count() {


    	using THash = ::fsc::TransformedHash<Kmer, Hash, Transform>;
    	using TLess = ::fsc::TransformedComparator<Kmer, Less, Transform>;
    	using Equal1 = ::fsc::TransformedComparator<Kmer, Equal, Transform>;
    	using Equal2 = ::fsc::sparsehash::compare<Kmer, Equal, Transform>;


    	auto test = make_kmer_multimap<Kmer, canonical, Transform, THash, Equal2>();
		::std::unordered_multimap<Kmer, uint32_t, THash, Equal1> gold;
		::std::vector<std::pair<Kmer, uint32_t> > entries;

		this->multimap_insert<canonical, TLess>(test, gold, entries);

		// get list of unique k-mers
		std::vector<Kmer> keys = test.keys();

		// assert that there are no duplicates
		std::unordered_set<Kmer, THash, Equal1> unique_keys(keys.begin(), keys.end());
		ASSERT_EQ(unique_keys.size(), keys.size());

		std::unordered_map<Kmer, uint32_t, THash, Equal1> unique_entries(gold.begin(), gold.end());
		ASSERT_EQ(unique_entries.size(), unique_keys.size());

		ASSERT_EQ(entries.size(), unique_keys.size());


    	// check one by one
    	for (auto i : entries) {
    		ASSERT_EQ(test.count(i.first), gold.count(i.first));
    	}
    }

};

// indicate this is a typed test
TYPED_TEST_CASE_P(DenseHashKmerMultimapTest);


template<typename K>
using HASH_K = ::bliss::kmer::hash::farm<K, false>;


TYPED_TEST_P(DenseHashKmerMultimapTest, single_map_insert)
{
	this->template test_map_insert<TypeParam, false,
									  ::bliss::transform::identity,
									  HASH_K, std::equal_to, std::less>();
}


TYPED_TEST_P(DenseHashKmerMultimapTest, single_map_equal_range)
{

  this->template test_map_equal_range<TypeParam, false,
								  ::bliss::transform::identity,
								  HASH_K, std::equal_to, std::less>();

}

TYPED_TEST_P(DenseHashKmerMultimapTest, single_map_count)
{
  this->template test_map_count<TypeParam, false,
								  ::bliss::transform::identity,
								  HASH_K, std::equal_to, std::less>();
}

TYPED_TEST_P(DenseHashKmerMultimapTest, canonical_map_insert)
{
	this->template test_map_insert<TypeParam, true,
									  ::bliss::transform::identity,
									  HASH_K, std::equal_to, std::less>();
}



TYPED_TEST_P(DenseHashKmerMultimapTest, canonical_map_equal_range)
{
	  this->template test_map_equal_range<TypeParam, true,
									  ::bliss::transform::identity,
									  HASH_K, std::equal_to, std::less>();
}

TYPED_TEST_P(DenseHashKmerMultimapTest, canonical_map_count)
{
	  this->template test_map_count<TypeParam, true,
									  ::bliss::transform::identity,
									  HASH_K, std::equal_to, std::less>();

}


TYPED_TEST_P(DenseHashKmerMultimapTest, bimolecule_map_insert)
{
	this->template test_map_insert<TypeParam, false,
									  ::bliss::kmer::transform::lex_less,
									  HASH_K, std::equal_to, std::less>();

	//  using SPLITTER = typename std::conditional<((TypeParam::nWords * sizeof(typename TypeParam::KmerWordType) * 8 - TypeParam::nBits) > 1),
//		  ::bliss::filter::TruePredicate,
//		   ::fsc::TransformedPredicate<TypeParam, ::bliss::utils::KmerInLowerSpace, ::bliss::kmer::transform::lex_less> >::type;
//
//
//  using THASH = ::fsc::TransformedHash<TypeParam, HASH_K, ::bliss::kmer::transform::lex_less >;
//
//  using EQUAL = ::bliss::kmer::hash::sparsehash::TransformedComparator<TypeParam, ::std::equal_to, ::bliss::kmer::transform::lex_less >;
//  using LESS = ::fsc::TransformedComparator<TypeParam, ::std::less, ::bliss::kmer::transform::lex_less >;
//
//  this->template test_map_insert<TypeParam, false, ((TypeParam::nWords * sizeof(typename TypeParam::KmerWordType) * 8 - TypeParam::nBits) <= 1), THASH, EQUAL, LESS>();
}


TYPED_TEST_P(DenseHashKmerMultimapTest, bimolecule_map_equal_range)
{
	 this->template test_map_equal_range<TypeParam, false,
	 ::bliss::kmer::transform::lex_less,
	 									  HASH_K, std::equal_to, std::less>();

//  using SPLITTER = typename std::conditional<((TypeParam::nWords * sizeof(typename TypeParam::KmerWordType) * 8 - TypeParam::nBits) > 1),
//      ::bliss::filter::TruePredicate,
//       ::fsc::TransformedPredicate<TypeParam, ::bliss::utils::KmerInLowerSpace, ::bliss::kmer::transform::lex_less> >::type;
//
//
//  using THASH = ::fsc::TransformedHash<TypeParam, HASH_K, ::bliss::kmer::transform::lex_less >;
//
//  using EQUAL = ::bliss::kmer::hash::sparsehash::TransformedComparator<TypeParam, ::std::equal_to, ::bliss::kmer::transform::lex_less >;
//  using LESS = ::fsc::TransformedComparator<TypeParam, ::std::less, ::bliss::kmer::transform::lex_less >;
//
//  this->template test_map_equal_range<TypeParam, false, ((TypeParam::nWords * sizeof(typename TypeParam::KmerWordType) * 8 - TypeParam::nBits) <= 1), THASH, EQUAL, LESS>();
}


TYPED_TEST_P(DenseHashKmerMultimapTest, bimolecule_map_count)
{
	  this->template test_map_count<TypeParam, false,
	  ::bliss::kmer::transform::lex_less,
	  									  HASH_K, std::equal_to, std::less>();

//  using SPLITTER = typename std::conditional<((TypeParam::nWords * sizeof(typename TypeParam::KmerWordType) * 8 - TypeParam::nBits) > 1),
//      ::bliss::filter::TruePredicate,
//       ::fsc::TransformedPredicate<TypeParam, ::bliss::utils::KmerInLowerSpace, ::bliss::kmer::transform::lex_less> >::type;
//
//
//  using THASH = ::fsc::TransformedHash<TypeParam, HASH_K, ::bliss::kmer::transform::lex_less >;
//
//  using EQUAL = ::bliss::kmer::hash::sparsehash::TransformedComparator<TypeParam, ::std::equal_to, ::bliss::kmer::transform::lex_less >;
//  using LESS = ::fsc::TransformedComparator<TypeParam, ::std::less, ::bliss::kmer::transform::lex_less >;
//
//  this->template test_map_count<TypeParam, false, ((TypeParam::nWords * sizeof(typename TypeParam::KmerWordType) * 8 - TypeParam::nBits) <= 1), THASH, EQUAL, LESS>();
}

TYPED_TEST_P(DenseHashKmerMultimapTest, single_multimap_insert)
{
//  using SPLITTER = typename std::conditional<((TypeParam::nWords * sizeof(typename TypeParam::KmerWordType) * 8 - TypeParam::nBits) > 1),
//      ::bliss::filter::TruePredicate,  ::bliss::utils::KmerInLowerSpace<TypeParam> >::type;
//
//  using HASH = ::bliss::kmer::hash::farm<TypeParam, false>;
//
//  using EQUAL = ::fsc::TransformedComparator<TypeParam, ::std::equal_to, ::bliss::transform::identity >;
//  using LESS = ::fsc::TransformedComparator<TypeParam, ::std::less, ::bliss::transform::identity >;
//
//  this->template test_multimap_insert<    TypeParam, false, SPLITTER, HASH, EQUAL, LESS>();

	this->template test_multimap_insert<TypeParam, false,
									  ::bliss::transform::identity,
									  HASH_K, std::equal_to, std::less>();

}
TYPED_TEST_P(DenseHashKmerMultimapTest, single_multimap_equal_range)
{
	  this->template test_multimap_equal_range<TypeParam, false,
									  ::bliss::transform::identity,
									  HASH_K, std::equal_to, std::less>();

}
TYPED_TEST_P(DenseHashKmerMultimapTest, single_multimap_count)
{
	  this->template test_multimap_count<TypeParam, false,
									  ::bliss::transform::identity,
									  HASH_K, std::equal_to, std::less>();

}


TYPED_TEST_P(DenseHashKmerMultimapTest, canonical_multimap_insert)
{
	this->template test_multimap_insert<TypeParam, true,
									  ::bliss::transform::identity,
									  HASH_K, std::equal_to, std::less>();

}


TYPED_TEST_P(DenseHashKmerMultimapTest, canonical_multimap_equal_range)
{
	  this->template test_multimap_equal_range<TypeParam, true,
									  ::bliss::transform::identity,
									  HASH_K, std::equal_to, std::less>();

}

TYPED_TEST_P(DenseHashKmerMultimapTest, canonical_multimap_count)
{
	  this->template test_multimap_count<TypeParam, true,
									  ::bliss::transform::identity,
									  HASH_K, std::equal_to, std::less>();

}

TYPED_TEST_P(DenseHashKmerMultimapTest, bimolecule_multimap_insert)
{
	this->template test_multimap_insert<TypeParam, false,
									  ::bliss::kmer::transform::lex_less,
									  HASH_K, std::equal_to, std::less>();


	//  using SPLITTER = typename std::conditional<((TypeParam::nWords * sizeof(typename TypeParam::KmerWordType) * 8 - TypeParam::nBits) > 1),
//      ::bliss::filter::TruePredicate,
//       ::fsc::TransformedPredicate<TypeParam, ::bliss::utils::KmerInLowerSpace, ::bliss::kmer::transform::lex_less> >::type;
//
//
//  using THASH = ::fsc::TransformedHash<TypeParam, HASH_K, ::bliss::kmer::transform::lex_less >;
//
//  using EQUAL = ::bliss::kmer::hash::sparsehash::TransformedComparator<TypeParam, ::std::equal_to, ::bliss::kmer::transform::lex_less >;
//  using LESS = ::fsc::TransformedComparator<TypeParam, ::std::less, ::bliss::kmer::transform::lex_less >;
//
//  this->template test_multimap_insert<TypeParam, false, ((TypeParam::nWords * sizeof(typename TypeParam::KmerWordType) * 8 - TypeParam::nBits) <= 1), THASH, EQUAL, LESS>();
}



TYPED_TEST_P(DenseHashKmerMultimapTest, bimolecule_multimap_equal_range)
{
	this->template test_multimap_equal_range<TypeParam, false,
									  ::bliss::kmer::transform::lex_less,
									  HASH_K, std::equal_to, std::less>();
}


TYPED_TEST_P(DenseHashKmerMultimapTest, bimolecule_multimap_count)
{
	this->template test_multimap_count<TypeParam, false,
									  ::bliss::kmer::transform::lex_less,
									  HASH_K, std::equal_to, std::less>();
}

// now register the test cases
REGISTER_TYPED_TEST_CASE_P(DenseHashKmerMultimapTest,
                           single_multimap_insert, single_multimap_equal_range, single_multimap_count,
                           canonical_multimap_insert, canonical_multimap_equal_range, canonical_multimap_count,
                           bimolecule_multimap_insert, bimolecule_multimap_equal_range, bimolecule_multimap_count,
                           single_map_insert, single_map_equal_range, single_map_count,
                           canonical_map_insert, canonical_map_equal_range, canonical_map_count,
                           bimolecule_map_insert, bimolecule_map_equal_range, bimolecule_map_count
                           );


//////////////////// RUN the tests with different types.

typedef ::testing::Types<
		::bliss::common::Kmer<7, ::bliss::common::DNA, uint16_t>,
		 ::bliss::common::Kmer<8, ::bliss::common::DNA, uint16_t>,
			::bliss::common::Kmer<5, ::bliss::common::DNA6, uint16_t>,
				::bliss::common::Kmer<3, ::bliss::common::DNA16, uint16_t>,
				 ::bliss::common::Kmer<4, ::bliss::common::DNA16, uint16_t>
		> DenseHashKmerMultimapTestTypes;
INSTANTIATE_TYPED_TEST_CASE_P(Bliss, DenseHashKmerMultimapTest, DenseHashKmerMultimapTestTypes);



