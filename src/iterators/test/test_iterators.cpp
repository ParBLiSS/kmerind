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

// include files to test
#include "iterators/many2one_iterator.hpp"
#include "utils/logging.h"


TEST(IteratorTests, TestCompressingIterator)
{
  std::string input = "2394" "1292" "1293" "4129" "4";
  typedef std::string::iterator it_t;
  struct MyFunctor
  {
    const int m = 4;
    int operator()(it_t& begin, it_t& end)
    {
      int sum = 0;
      int i = 0;
      while (begin != end && i++ < m)
      {
        sum += static_cast<int>(*begin - '0');
        ++begin;
      }
      return sum;
    }
  } f;
  typedef bliss::iterator::many2one_iterator<it_t, MyFunctor> comp_it_t;
  comp_it_t comp_it(input.begin(), input.end(), f, 4);
  comp_it_t comp_it_end(input.end(), input.end(), f, 4);

  while (comp_it != comp_it_end)
  {
    comp_it--;
    BL_INFO( *(comp_it+2) << ", " );
    comp_it++;
    comp_it += 1;

  }
}
