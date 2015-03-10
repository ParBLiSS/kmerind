/**
 * @file    iterator_test_utils.hpp
 * @ingroup utils
 * @author  Tony Pan <tpan7@gatech.edu>
 * @brief   convenience functions to compare sequences of items
 * @details
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef TEST_UTILS_HPP_
#define TEST_UTILS_HPP_

#include <vector>
#include <iterator>
#include <algorithm>  // for sort
#include <sstream>

#include "utils/logging.h"

/// assign a sequence of numeric values to an iterator.
template<typename Iterator>
void assignSequence(Iterator data, size_t count,
                    typename std::iterator_traits<Iterator>::value_type offset = 0,
                    typename std::iterator_traits<Iterator>::value_type step = 1) {
  for (size_t i = 0; i < count; ++i) {
    data[i] = (offset + i * step);
  }
}

/**
 * @brief  compare 2 iterators up to count number of entries, and print out the number of items that are different.
 */
template<typename Iterator1, typename Iterator2>
bool compareSequences(Iterator1 data1, Iterator2 data2, size_t count) {
  static_assert(std::is_same<typename std::iterator_traits<Iterator1>::value_type,
                typename std::iterator_traits<Iterator2>::value_type >::value, "iterators have to have same value type.");

  // t1 - t2
  std::vector<typename std::iterator_traits<Iterator1>::value_type> t1only;
  std::set_difference(data1, data1 + count, data2, data2 + count, std::inserter(t1only, t1only.begin()));

  // t2 - t1
  std::vector<typename std::iterator_traits<Iterator2>::value_type> t2only;
  std::set_difference(data2, data2 + count, data1, data1 + count, std::inserter(t2only, t2only.begin()));

  bool same = t1only.size() == 0 && t2only.size() == 0;

  if (!same) {
    std::ostringstream ss;
    std::ostream_iterator<int> oit(ss, ",");
    if (t1only.size() > 0) {
      ss.str("");
      ss.clear();
      std::copy(t1only.begin(), t1only.end(), oit);
      INFOF("\t first buffer exclusive: size=%lu, entries=[%s]\n", t1only.size(), ss.str().c_str());
    }
    if (t2only.size() > 0) {
      ss.str("");
      ss.clear();
      std::copy(t2only.begin(), t2only.end(), oit);
      INFOF("\t second buffer exclusive: size=%lu, entries=[%s]\n", t2only.size(), ss.str().c_str());
    }
  }
  return same;
}


/// compare a input iterator with a numeric sequence
template<typename Iterator>
bool checkSequence(Iterator data, size_t count,
                   typename std::iterator_traits<Iterator>::value_type offset = 0,
                   typename std::iterator_traits<Iterator>::value_type step = 1) {
  std::vector<typename std::iterator_traits<Iterator>::value_type> seq(count);

  for (size_t i = 0; i < count; ++i) {
    seq[i] = (offset + i * step);
  }
  return compareSequences(data, seq.begin(), count);
}

/// compare an iterator's content to a numeric sequence, first sorting the content.
template<typename Iterator>
bool checkUnorderedSequence(Iterator data, size_t count,
                            typename std::iterator_traits<Iterator>::value_type offset = 0,
                            typename std::iterator_traits<Iterator>::value_type step = 1) {

  std::vector<typename std::iterator_traits<Iterator>::value_type> temp(data, data + count);
  std::sort(temp.begin(), temp.end());

  return checkSequence(temp.begin(), count, offset, step);
}


/// compare 2 iterators contents to each other, first sorting the content.
template<typename Iterator1, typename Iterator2>
bool compareUnorderedSequences(Iterator1 data, Iterator2 data2, size_t count) {
  static_assert(std::is_same<typename std::iterator_traits<Iterator1>::value_type,
                typename std::iterator_traits<Iterator2>::value_type >::value, "iterators have to have same value type.");

  // copy
  std::vector<typename std::iterator_traits<Iterator1>::value_type> t1(data, data + count);
  std::vector<typename std::iterator_traits<Iterator2>::value_type> t2(data2, data2 + count);

  // sort
  std::sort(t1.begin(), t1.end());
  std::sort(t2.begin(), t2.end());

  return compareSequences(t1.begin(), t2.begin(), count);
}


#endif /* TEST_UTILS_HPP_ */
