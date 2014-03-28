/**
 * @file    iterator_tools.hpp
 * @ingroup iterators
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   Implements common, shared utility functions for iterators.
 *
 * Copyright (c) TODO
 *
 * TODO add Licence
 */

#ifndef BLISS_ITERATORS_ITERATOR_TOOLS_HPP
#define BLISS_ITERATORS_ITERATOR_TOOLS_HPP

#include <iterator>

namespace bliss
{
namespace iterator
{

// implementation of iter_tools for random access itertors
template<typename Iterator>
struct _iter_tools_ra
{
  /// for convenience, keep the iterator_traits in here
  typedef typename std::iterator_traits<Iterator> traits;

  /// The difference type of this iterator
  typedef typename traits::difference_type difference_type;

  static void _advance_at_most(Iterator& it, const Iterator& end,
                              const difference_type m)
  {
    // distance and advance is both possible in constant time for
    // random access iterators:
    std::advance(it, std::min(std::distance(it, end), m));
  }
};


// implementation of iter_tools for NON random access iterators
template<typename Iterator>
struct _iter_tools_dir
{
  /// The difference type of this iterator
  typedef typename std::iterator_traits<Iterator>::difference_type difference_type;

  static void _advance_at_most(Iterator& it, const Iterator& end,
                              const difference_type m)
  {
    // no random access, need to iterate forward manually
    difference_type remaining = m;
    while (it != end && (remaining--) > 0)
      ++it;
  }
};

template<typename Iterator>
struct iter_tools
{
  typedef typename std::conditional<std::is_same<typename std::iterator_traits<Iterator>::iterator_category,
            std::random_access_iterator_tag>::value,
            _iter_tools_ra<Iterator>,
            _iter_tools_dir<Iterator> >::type base_class_type;

  /// The difference type of this iterator
  typedef typename std::iterator_traits<Iterator>::difference_type difference_type;

  static void advance_at_most(Iterator& it, const Iterator& end,
                              const difference_type m)
  {
    base_class_type::_advance_at_most(it, end, m);
  }
};


} // namespace iterator
} // namespace bliss

#endif // BLISS_ITERATORS_ITERATOR_TOOLS_HPP
