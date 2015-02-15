/**
 * @file    iterator_utils.hpp
 * @ingroup interators
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   Implements helper functions for implementation of iterators.
 *
 * Copyright (c) TODO
 *
 * TODO add Licence
 */


#ifndef BLISS_ITERATORS_ITERATOR_UTILS_HPP
#define BLISS_ITERATORS_ITERATOR_UTILS_HPP

#include <iterator>

namespace bliss
{
namespace iterator
{

/*********************************************************************
 *         Iterator tools (functions: advance_at_most, ...)
 *********************************************************************/


// implementation of iter_tools for random access itertors
template<typename Iterator>
struct _iter_tools_ra
{
  /// for convenience, keep the iterator_traits in here
  typedef typename std::iterator_traits<Iterator> traits;

  /// The difference type of this iterator
  typedef typename traits::difference_type difference_type;

  /// random access iterator jump ahead by m position, or to end iterator
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

  /// non-random access iterator jump ahead by m position, or to end iterator
  static void _advance_at_most(Iterator& it, const Iterator& end,
                              const difference_type m)
  {
    // no random access, need to iterate forward manually
    difference_type remaining = m;
    while (it != end && (remaining--) > 0)
      ++it;
  }
};

/**
 * @struct  Iterator tools.
 * @brief   Implements utility functions depending on the given iterator type.
 *
 * @tparam Iterator  The iterator type used in the functions.
 */
template<typename Iterator>
struct iter_tools
{
private:
  // Get the base class implementing the functionality for the
  // given iterator type
  typedef typename std::conditional<
              std::is_same<typename std::iterator_traits<Iterator>::iterator_category,
              std::random_access_iterator_tag>::value,
          _iter_tools_ra<Iterator>,
          _iter_tools_dir<Iterator> >::type base_class_type;
public:
  /// The difference type of this iterator
  typedef typename std::iterator_traits<Iterator>::difference_type difference_type;

  /**
   * @brief Advances the given iterator by `m` positions, but at most till the
   *        end iterator given by `end`.
   *
   * @param it[in|out]  The iterator to be advanced.
   * @param end         The maximum iterator position to iterate to.
   * @param m           The distance to iterate by.
   */
  static void advance_at_most(Iterator& it, const Iterator& end,
                              const difference_type m)
  {
    base_class_type::_advance_at_most(it, end, m);
  }
};



/*********************************************************************
 *      Iterator tags helper (inherit at most, tag types, etc)       *
 *********************************************************************/

/**
 * @brief A utility structure to ease handleing of iterator_tags.
 *
 * This structure supplies a variety of templated typedefs and constexpressions
 * to make handleing of iterator_tags easier. This includes tag `inheritance`.
 * The members `inherit_bidir` and `inherit_forward` will return an iterator
 * tag which is the same as the iterator tag from the given Iterator, but at
 * most bidirectional or forward (respectively).
 *
 * @tparam Iterator     The iterator whose iterator_tag is used in this helper
 *                      struct.
 */
template<typename Iterator>
struct iterator_tags
{
  /// The iterator traits of the given Iterator
  typedef typename std::iterator_traits<Iterator> traits;

  /// Whether the iterator is a random access iterator
  static constexpr bool is_random_access =
    std::is_same<typename traits::iterator_category,
                 std::random_access_iterator_tag>::value;

  /// Whether the iterator is a bidirectional iterator
  static constexpr bool is_bidirectional =
    std::is_same<typename traits::iterator_category,
                 std::bidirectional_iterator_tag>::value;

  /// Whether the iterator is a forward iterator
  static constexpr bool is_forward =
    std::is_same<typename traits::iterator_category,
                 std::forward_iterator_tag>::value;

  /// Whehter the iterator is at least a bidirectional iterator
  static constexpr bool is_min_bidir = is_random_access || is_bidirectional;

  /// Inherit the base iterator type but at most till forward iterator
  typedef typename std::conditional<is_min_bidir,
                // then: become a forward iterator
                std::forward_iterator_tag,
                // else: `steal` the underlying iterator type
                typename traits::iterator_category>::type inherit_forward;

  /// Inherit the base iterator type but at most till bidirectional iterator
  typedef typename std::conditional<is_random_access,
                // then: become a forward iterator
                std::bidirectional_iterator_tag,
                // else: `steal` the underlying iterator type
                typename traits::iterator_category>::type inherit_bidir;
};

} // iterator
} // bliss
#endif /* BLISS_ITERATORS_ITERATOR_UTILS_HPP */
