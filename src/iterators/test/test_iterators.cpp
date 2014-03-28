
// include google test
#include <gtest/gtest.h>

#include <string>

// include files to test
#include <iterators/compression_iterator.hpp>


TEST(IteratorTests, TestCompressingIterator)
{
  std::string input = "23941292129341294214";
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
      }
      return sum;
    }
  } f;
  bliss::iterator::compressing_iterator<MyFunctor, it_t> comp_it(input.begin(), input.end(), f, 4);
}
