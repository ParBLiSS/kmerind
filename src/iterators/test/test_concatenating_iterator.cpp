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

#include <string>
#include <iterator>
#include <iostream>
#include <algorithm>
#include <vector>

// include files to test
#include "iterators/concatenating_iterator.hpp"
#include "iterators/container_concatenating_iterator.hpp"

TEST(ConcatenatingIteratorTests, TestConcatenatingString)
{
  std::string gold = "123443210000bam!";
  std::string input1 = "1234";
  std::string input2 = "4321";
  std::string input3 = "0000";
  std::string input4 = "bam!";

  bliss::iterator::ConcatenatingIterator<std::string::const_iterator> start;
  bliss::iterator::ConcatenatingIterator<std::string::const_iterator> end;

  start.addRange(input1.cbegin(), input1.cend());
  start.addRange(input2.cbegin(), input2.cend());
  start.addRange(input3.cbegin(), input3.cend());
  start.addRange(input4.cbegin(), input4.cend());

  int i = 0;
  for (; start != end; ++start, ++i) {
    EXPECT_EQ(*start, gold[i]);
  }  
  EXPECT_EQ(i, 16);
}

TEST(ConcatenatingIteratorTests, TestRestartString)
{
  std::string gold1 = "12344321";
  std::string gold2 = "0000bam!";
  std::string input1 = "1234";
  std::string input2 = "4321";
  std::string input3 = "0000";
  std::string input4 = "bam!";

  bliss::iterator::ConcatenatingIterator<std::string::const_iterator> start;
  bliss::iterator::ConcatenatingIterator<std::string::const_iterator> end;

  start.addRange(input1.cbegin(), input1.cend());
  start.addRange(input2.cbegin(), input2.cend());


  int i = 0;
  for (; start != end; ++start, ++i) {
    EXPECT_EQ(gold1[i], *start);
  }  
  EXPECT_EQ(i, 8);

  start.addRange(input3.cbegin(), input3.cend());
  start.addRange(input4.cbegin(), input4.cend());

  EXPECT_TRUE(start != end);

  i = 0;
  for (; start != end; ++start, ++i) {
    EXPECT_EQ(gold2[i], *start);
  }  
  EXPECT_EQ(i, 8);
}


TEST(ConcatenatingIteratorTests, TestConcatenatingInts)
{
  std::vector<int> gold {1, 2, 3, 4, 4, 3, 2, 1, 0, 0, 0, 0 };
  std::vector<int> input1 {1, 2, 3, 4 };
  std::vector<int> input2 {4, 3, 2, 1 };
  std::vector<int> input3 {0, 0, 0, 0 };

  bliss::iterator::ConcatenatingIterator<std::vector<int>::iterator> start;
  bliss::iterator::ConcatenatingIterator<std::vector<int>::iterator> end;

  start.addRange(input1.begin(), input1.end());
  start.addRange(input2.begin(), input2.end());
  start.addRange(input3.begin(), input3.end());

  int i = 0;
  for (; start != end; ++start, ++i) {
    EXPECT_EQ(*start, gold[i]);
  }  
  EXPECT_EQ(i, 12);
}

TEST(ConcatenatingIteratorTests, TestConcatenatingIntConstIter)
{
  std::vector<int> gold {1, 2, 3, 4, 4, 3, 2, 1, 0, 0, 0, 0 };
  std::vector<int> input1 {1, 2, 3, 4 };
  std::vector<int> input2 {4, 3, 2, 1 };
  std::vector<int> input3 {0, 0, 0, 0 };

  bliss::iterator::ConcatenatingIterator<std::vector<int>::const_iterator> start;
  bliss::iterator::ConcatenatingIterator<std::vector<int>::const_iterator> end;

  start.addRange(input1.cbegin(), input1.cend());
  start.addRange(input2.cbegin(), input2.cend());
  start.addRange(input3.cbegin(), input3.cend());

  int i = 0;
  for (; start != end; ++start, ++i) {
    EXPECT_EQ(*start, gold[i]);
  }  
  EXPECT_EQ(i, 12);
}
//====  cannot test vector<const int>, since the valuetype has to be has to be assignable.




TEST(ConcatenatingIteratorTests, TestRestartInts)
{
  std::vector<int> gold1 {1, 2, 3, 4, 4, 3, 2, 1};
  std::vector<int> gold2 { 0, 0, 0, 0 };
  std::vector<int> input1 {1, 2, 3, 4 };
  std::vector<int> input2 {4, 3, 2, 1 };
  std::vector<int> input3 {0, 0, 0, 0 };

  bliss::iterator::ConcatenatingIterator<std::vector<int>::iterator> start;
  bliss::iterator::ConcatenatingIterator<std::vector<int>::iterator> end;

  start.addRange(input1.begin(), input1.end());
  start.addRange(input2.begin(), input2.end());


  int i = 0;
  for (; start != end; ++start, ++i) {
    EXPECT_EQ(gold1[i], *start);
  }  
  EXPECT_EQ(i, 8);

  start.addRange(input3.begin(), input3.end());

  EXPECT_TRUE(start != end);

  i = 0;
  for (; start != end; ++start, ++i) {
    EXPECT_EQ(gold2[i], *start);
  }  
  EXPECT_EQ(i, 4);
}


TEST(ConcatenatingIteratorTests, TestConcatenatingReverse)
{
  std::vector<int> gold {1, 2, 3, 4, 4, 3, 2, 1, 0, 0, 0, 0 };
  std::vector<int> input1 {1, 2, 3, 4 };
  std::vector<int> input2 {4, 3, 2, 1 };
  std::vector<int> input3 {0, 0, 0, 0 };

  bliss::iterator::ConcatenatingIterator<std::vector<int>::iterator> start;
  bliss::iterator::ConcatenatingIterator<std::vector<int>::iterator> end;

  start.addRange(input1.begin(), input1.end());
  start.addRange(input2.begin(), input2.end());
  start.addRange(input3.begin(), input3.end());

  bliss::iterator::ConcatenatingIterator<std::vector<int>::iterator> it(start);
  
  for (; it != end; ++it);
  
  --it;
  int i = 11;
  for (; it != start; --it, --i) 
  {
    EXPECT_EQ(*it, gold[i]);
  }  
  EXPECT_EQ(*it, gold[i]);  // first entry.
  EXPECT_EQ(i, 0);
}

TEST(ConcatenatingIteratorTests, TestConcatenatingModify)
{
  std::vector<int> gold {1, 2, 3, 4, 4, 3, 2, 1, 0, 0, 0, 0 };
  std::vector<int> gold2 {2, 4, 6, 8, 8, 6, 4, 2, 0, 0, 0, 0 };
  std::vector<int> input1 {1, 2, 3, 4 };
  std::vector<int> input2 {4, 3, 2, 1 };
  std::vector<int> input3 {0, 0, 0, 0 };

  bliss::iterator::ConcatenatingIterator<std::vector<int>::iterator> start;
  bliss::iterator::ConcatenatingIterator<std::vector<int>::iterator> end;

  start.addRange(input1.begin(), input1.end());
  start.addRange(input2.begin(), input2.end());
  start.addRange(input3.begin(), input3.end());

  bliss::iterator::ConcatenatingIterator<std::vector<int>::iterator> it(start);
  
  for (; it != end; ++it);
  
  --it;
  int i = 11;
  for (; it != start; --it, --i) 
  {
    *it *= 2;
  }
  *it *= 2;

  i = 0;
  for (; it != end; ++it, ++i) {
    EXPECT_EQ(*it, gold2[i]);
  }
  EXPECT_EQ(i, 12);
}


TEST(ConcatenatingIteratorTests, TestSameIterator)
{
  std::string input1 = "1234";

  bliss::iterator::ConcatenatingIterator<std::string::const_iterator> start1;
  bliss::iterator::ConcatenatingIterator<std::string::const_iterator> start2;
  bliss::iterator::ConcatenatingIterator<std::string::const_iterator> end;

  start1.addRange(input1.cbegin(), input1.cend());
  start2.addRange(input1.cbegin(), input1.cend());

  int i = 0;
  for (; start1 != end; ++start1, ++i, ++start2);
  EXPECT_EQ(i, 4);
  EXPECT_TRUE(start1 == start2);
}

TEST(ConcatenatingIteratorTests, TestEmptyIterator)
{

  bliss::iterator::ConcatenatingIterator<std::string::const_iterator> start;
  bliss::iterator::ConcatenatingIterator<std::string::const_iterator> end;


  int i = 0;
  for (; start != end; ++start, ++i);
  EXPECT_EQ(i, 0);
}


class ContainerConcatenatingIteratorTests : public ::testing::Test
{
protected:
	using ContainerAdapter = bliss::iterator::ConcatenatingIteratorContainerAdapter<int>;  // int is dummy

	using InnerContainer = std::vector<size_t>;
	using OuterContainer = std::vector<InnerContainer>;

	using Iter = bliss::iterator::ContainerConcatenatingIterator<typename OuterContainer::const_iterator, ContainerAdapter>;

	OuterContainer outer;
	ContainerAdapter adapter;

	size_t inner_container_count = 100;

	std::vector<size_t> counts;
	std::vector<size_t> cumulative;

	std::vector<size_t> gold;

    virtual ~ContainerConcatenatingIteratorTests() {};

    virtual void SetUp()
    {
    	// get counts.
    	counts.reserve(inner_container_count);
    	for (size_t i = 0; i < inner_container_count; ++i) {
    		counts.emplace_back(i);
    		std::cout << i << " ";
    	}
    	std::random_shuffle(counts.begin(), counts.end());
    	std::cout << std::endl;

    	// exclusive prefix sum of counts
    	cumulative.resize(inner_container_count + 1);
    	cumulative[0] = 0;
    	for (size_t i = 1; i <= inner_container_count; ++i) {
    		std::cout << counts[i-1] << " ";
    		cumulative[i] = cumulative[i-1] + counts[i - 1];
    	}
    	std::cout << std::endl;


    	gold.reserve(cumulative[inner_container_count]);
    	auto bi = ::std::back_inserter(gold);

    	// populate the outer container.
    	outer.reserve(inner_container_count);
    	for (size_t i = 0; i < inner_container_count; ++i) {
    		std::vector<size_t> inner;
    		if (counts[i] > 0) {
    			inner.resize(counts[i]);
    			std::fill(inner.begin(), inner.end(), i);
    			std::fill_n(bi, counts[i], i);
    		}
    		outer.emplace_back(std::move(inner));
    	}

    }

    template <typename Iter>
    bool compare_to_gold(Iter start, Iter end, size_t start_pos, size_t end_pos) {
    	size_t dist = std::distance(start, end);
    	EXPECT_EQ(end_pos - start_pos, dist);

    	bool same = true;
    	auto it = start;
    	size_t pos = start_pos;
    	for (; it != end && pos != end_pos; ++it, ++pos) {
    		if (*it != gold[pos]) std::cout << " pos " << pos << " test: [" << *it << "] gold [" << gold[pos] << "]" << std::endl;
    		same &= (*it == gold[pos]);
    	}
    	for (; it != end; ++it) {
    		std::cout << " test only: [" << *it << "]"  << std::endl;
    	}
    	for (; pos != end_pos; ++pos) {
    		std::cout << " pos " << pos << " gold only [" << gold[pos] << "]" << std::endl;
    	}

    	EXPECT_TRUE(same);
    	EXPECT_EQ(it, end);
    	EXPECT_EQ(pos, end_pos);

    	return same;
    }
};


TEST_F(ContainerConcatenatingIteratorTests, TestDefaultConstructed)
{
	using ContainerAdapter = bliss::iterator::ConcatenatingIteratorContainerAdapter<int>;  // int is dummy

	using InnerContainer = std::vector<size_t>;
	using OuterContainer = std::vector<InnerContainer>;

  bliss::iterator::ContainerConcatenatingIterator<typename OuterContainer::const_iterator, ContainerAdapter> start;
  bliss::iterator::ContainerConcatenatingIterator<typename OuterContainer::const_iterator, ContainerAdapter> end;

  this->compare_to_gold(start, end, 0, 0);

}


TEST_F(ContainerConcatenatingIteratorTests, TestEmptyVectorVector)
{
	using ContainerAdapter = bliss::iterator::ConcatenatingIteratorContainerAdapter<size_t>;  // int is dummy

	ContainerAdapter adapter;

	using InnerContainer = std::vector<size_t>;
	using OuterContainer = std::vector<InnerContainer>;

	OuterContainer outer;

	using Iter = bliss::iterator::ContainerConcatenatingIterator<typename OuterContainer::const_iterator, ContainerAdapter>;

	Iter end(adapter, outer.cend());

	{
	  Iter start(adapter, outer.cbegin(), outer.cend());

	  this->compare_to_gold(start, end, 0, 0);
	}
}

TEST_F(ContainerConcatenatingIteratorTests, ConstructorTest1)
{

	// test different sets of constructor and destructors
	{
		auto oend = this->outer.cend();

		Iter end(this->adapter, oend);

		{
		  Iter start(this->adapter, this->outer.cbegin(), oend);

		  this->compare_to_gold(start, end, 0, this->cumulative[this->inner_container_count]);
		}

		{
		  Iter start(this->adapter, oend, oend);

		  this->compare_to_gold(start, end, this->cumulative[this->inner_container_count], this->cumulative[this->inner_container_count]);
		}

		{
		  Iter start(this->adapter, this->outer.cbegin(), this->outer[0].cbegin(), oend);

		  this->compare_to_gold(start, end, 0, this->cumulative[this->inner_container_count]);
		}

		{
		  Iter start(this->adapter, this->outer.cbegin() + 4, oend);

		  this->compare_to_gold(start, end, this->cumulative[4], this->cumulative[this->inner_container_count]);
		}

		{
		  Iter start(this->adapter, this->outer.cbegin() + 4, this->outer[4].cbegin(), oend);

		  this->compare_to_gold(start, end, this->cumulative[4], this->cumulative[this->inner_container_count]);
		}

		{
		  Iter start(this->adapter, this->outer.cbegin() + 3, this->outer[3].cend(), oend);

		  this->compare_to_gold(start, end, this->cumulative[4], this->cumulative[this->inner_container_count]);
		}

		{
		  Iter start(this->adapter, this->outer.cbegin() + 4, this->outer[4].cbegin() + this->counts[4] / 2, oend);

		  this->compare_to_gold(start, end, this->cumulative[4] + this->counts[4] / 2, this->cumulative[this->inner_container_count]);
		}
	}


}


TEST_F(ContainerConcatenatingIteratorTests, ConstructorTest2)
{

	{
		auto oend = this->outer.cbegin() + 10;
		Iter end(this->adapter, oend);

		{
		  Iter start(this->adapter, this->outer.cbegin(), oend);

		  this->compare_to_gold(start, end, 0, this->cumulative[10]);
		}

		{
		  Iter start(this->adapter, oend, oend);

		  this->compare_to_gold(start, end, this->cumulative[10], this->cumulative[10]);
		}

		{
		  Iter start(this->adapter, oend, this->adapter.begin(*oend), oend);

		  this->compare_to_gold(start, end, this->cumulative[10], this->cumulative[10]);
		}


		{
		  Iter start(this->adapter, this->outer.cbegin(), this->outer[0].cbegin(), oend);

		  this->compare_to_gold(start, end, 0, this->cumulative[10]);
		}

		{
		  Iter start(this->adapter, this->outer.cbegin() + 4, oend);

		  this->compare_to_gold(start, end, this->cumulative[4], this->cumulative[10]);
		}

		{
		  Iter start(this->adapter, this->outer.cbegin() + 4, this->outer[4].cbegin(), oend);

		  this->compare_to_gold(start, end, this->cumulative[4], this->cumulative[10]);
		}

		{
		  Iter start(this->adapter, this->outer.cbegin() + 3, this->outer[3].cend(), oend);

		  this->compare_to_gold(start, end, this->cumulative[4], this->cumulative[10]);
		}

		{
		  Iter start(this->adapter, this->outer.cbegin() + 4, this->outer[4].cbegin() + this->counts[4] / 2, oend);

		  this->compare_to_gold(start, end, this->cumulative[4] + this->counts[4] / 2, this->cumulative[10]);
		}

	}

}


TEST_F(ContainerConcatenatingIteratorTests, ConstructorTest3)
{

	{
		auto oend = this->outer.cbegin() + 10;
		auto iend = this->outer[10].cbegin();

		Iter end(this->adapter, oend, iend);

		{
			std::cout << "1" << std::endl;
		  Iter start(this->adapter, this->outer.cbegin(), oend, iend);

		  this->compare_to_gold(start, end, 0, this->cumulative[10]);
		}

		{
			std::cout << "1" << std::endl;
		  Iter start(this->adapter, oend, oend, iend);

		  this->compare_to_gold(start, end, this->cumulative[10], this->cumulative[10]);
		}

		{
			std::cout << "1" << std::endl;
		  Iter start(this->adapter, oend, iend, oend, iend);

		  this->compare_to_gold(start, end, this->cumulative[10], this->cumulative[10]);
		}

		{
			std::cout << "2" << std::endl;
		  Iter start(this->adapter, this->outer.cbegin(), this->outer[0].cbegin(), oend, iend);

		  this->compare_to_gold(start, end, 0, this->cumulative[10]);
		}

		{
			std::cout << "3" << std::endl;

		  Iter start(this->adapter, this->outer.cbegin() + 4, oend, iend);

		  std::cout << "3. prepped." << std::endl;
		  this->compare_to_gold(start, end, this->cumulative[4], this->cumulative[10]);
		}

		{
			std::cout << "4" << std::endl;

		  Iter start(this->adapter, this->outer.cbegin() + 4, this->outer[4].cbegin(), oend, iend);

		  this->compare_to_gold(start, end, this->cumulative[4], this->cumulative[10]);
		}

		{
			std::cout << "5" << std::endl;

		  Iter start(this->adapter, this->outer.cbegin() + 3, this->outer[3].cend(), oend, iend);

		  this->compare_to_gold(start, end, this->cumulative[4], this->cumulative[10]);
		}

		{
			std::cout << "6" << std::endl;

		  Iter start(this->adapter, this->outer.cbegin() + 4, this->outer[4].cbegin() + this->counts[4] / 2, oend, iend);

		  this->compare_to_gold(start, end, this->cumulative[4] + this->counts[4] / 2, this->cumulative[10]);
		}

	}

}

TEST_F(ContainerConcatenatingIteratorTests, ConstructorTest4)
{

	{
		auto oend = this->outer.cbegin() + 9;
		auto iend = this->outer[9].cend();

		Iter end(this->adapter, oend, iend);

		{
			std::cout << "1" << std::endl;
		  Iter start(this->adapter, this->outer.cbegin(), oend, iend);

		  this->compare_to_gold(start, end, 0, this->cumulative[10]);
		}

		{
			std::cout << "1" << std::endl;
		  Iter start(this->adapter, oend, iend, oend, iend);

		  this->compare_to_gold(start, end, this->cumulative[10], this->cumulative[10]);
		}

		{
			std::cout << "1" << std::endl;
		  Iter start(this->adapter, oend, oend, iend);

		  this->compare_to_gold(start, end, this->cumulative[9], this->cumulative[10]);
		}

		{
			std::cout << "2" << std::endl;
		  Iter start(this->adapter, this->outer.cbegin(), this->outer[0].cbegin(), oend, iend);

		  this->compare_to_gold(start, end, 0, this->cumulative[10]);
		}

		{
			std::cout << "3" << std::endl;

		  Iter start(this->adapter, this->outer.cbegin() + 4, oend, iend);

		  this->compare_to_gold(start, end, this->cumulative[4], this->cumulative[10]);
		}

		{
			std::cout << "4" << std::endl;
		  Iter start(this->adapter, this->outer.cbegin() + 4, this->outer[4].cbegin(), oend, iend);

		  this->compare_to_gold(start, end, this->cumulative[4], this->cumulative[10]);
		}

		{
			std::cout << "5" << std::endl;
		  Iter start(this->adapter, this->outer.cbegin() + 3, this->outer[3].cend(), oend, iend);

		  this->compare_to_gold(start, end, this->cumulative[4], this->cumulative[10]);
		}

		{
			std::cout << "6" << std::endl;
		  Iter start(this->adapter, this->outer.cbegin() + 4, this->outer[4].cbegin() + this->counts[4] / 2, oend, iend);

		  this->compare_to_gold(start, end, this->cumulative[4] + this->counts[4] / 2, this->cumulative[10]);
		}

	}
}



TEST_F(ContainerConcatenatingIteratorTests, ConstructorTest5)
{

	{
		auto oend = this->outer.cbegin() + 10;
		auto iend = this->outer[10].cbegin() + this->counts[10] / 2;

		Iter end(this->adapter, oend, iend);

		{
			std::cout << "1" << std::endl;

		  Iter start(this->adapter, this->outer.cbegin(), oend, iend);

		  this->compare_to_gold(start, end, 0, this->cumulative[10] + this->counts[10] / 2);
		}

		{
			std::cout << "1" << std::endl;

		  Iter start(this->adapter, oend, iend, oend, iend);

		  this->compare_to_gold(start, end, this->cumulative[10] + this->counts[10] / 2, this->cumulative[10] + this->counts[10] / 2);
		}

		{
			std::cout << "2" << std::endl;
		  Iter start(this->adapter, this->outer.cbegin(), this->outer[0].cbegin(), oend, iend);

		  this->compare_to_gold(start, end, this->cumulative[0], this->cumulative[10] + this->counts[10] / 2);
		}

		{
			std::cout << "1" << std::endl;

		  Iter start(this->adapter, oend, oend, iend);

		  this->compare_to_gold(start, end, this->cumulative[10], this->cumulative[10] + this->counts[10] / 2);
		}

		{
			std::cout << "3" << std::endl;
		  Iter start(this->adapter, this->outer.cbegin() + 4, oend, iend);

		  this->compare_to_gold(start, end, this->cumulative[4], this->cumulative[10] + this->counts[10] / 2);
		}

		{
			std::cout << "4" << std::endl;

		  Iter start(this->adapter, this->outer.cbegin() + 4, this->outer[4].cbegin(), oend, iend);

		  this->compare_to_gold(start, end, this->cumulative[4], this->cumulative[10] + this->counts[10] / 2);
		}

		{
			std::cout << "5" << std::endl;

		  Iter start(this->adapter, this->outer.cbegin() + 3, this->outer[3].cend(), oend, iend);

		  this->compare_to_gold(start, end, this->cumulative[4], this->cumulative[10] + this->counts[10] / 2);
		}

		{
			std::cout << "6" << std::endl;

		  Iter start(this->adapter, this->outer.cbegin() + 4, this->outer[4].cbegin() + this->counts[4] / 2, oend, iend);

		  this->compare_to_gold(start, end, this->cumulative[4] + this->counts[4] / 2, this->cumulative[10] + this->counts[10] / 2);
		}

	}

}





