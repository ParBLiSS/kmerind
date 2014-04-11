
// include google test
#include <gtest/gtest.h>

#include <string>

// include files to test
#include <iterators/many2one_iterator.hpp>


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
    std::cout << *(comp_it+2) << ", " << std::endl;
    comp_it++;
    comp_it += 1;

  }
}
