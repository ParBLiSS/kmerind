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

// include files to test
#include "iterators/concatenating_iterator.hpp"


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
