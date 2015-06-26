/**
 * @file    vector_utils.hpp
 * @ingroup
 * @author  tpan
 * @brief
 * @details
 *
 * Copyright (c) 2015 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef SRC_WIP_CONTAINER_UTILS_HPP_
#define SRC_WIP_CONTAINER_UTILS_HPP_

#include <iterator>  // iterator_traits

namespace fsc {

  /// append to container via emplace.
  /// modified based on http://stackoverflow.com/questions/18724999/why-no-emplacement-iterators-in-c11-or-c14
  template<class Container>
  class back_emplace_iterator : public std::iterator< std::output_iterator_tag,
                                                      void, void, void, void >
  {
  protected:
      Container* container;
  public:
      typedef Container container_type;

      explicit back_emplace_iterator(Container& x) : container(&x) {}

      template<class... Args>
      back_emplace_iterator<Container>& operator=(Args&&... args)
      {
          container->emplace_back(::std::forward<Args>(args)...);
          return *this;
      }

      back_emplace_iterator& operator*() { return *this; }
      back_emplace_iterator& operator++() { return *this; }
      back_emplace_iterator& operator++(int) { return *this; }
  };


  /// linear or logarithmic search for lowerbound.  assumes input is sorted.
  template <bool linear, class Iterator,
    class Less = ::std::less<typename ::std::iterator_traits<Iterator>::value_type>,
    class V = typename ::std::iterator_traits<Iterator>::value_type>
  inline Iterator lower_bound(Iterator b, Iterator e, V& v, Less lt = Less()) {
      // compiler choose one.
      if (linear) while ((b != e) && lt(*b, v)) ++b;
      else b = ::std::lower_bound(b, e, v, lt);
      return b;
  }

  /// linear or logarithmic search for upperbound.  assumes input is sorted.
  template <bool linear, class Iterator,
    class Less = ::std::less<typename ::std::iterator_traits<Iterator>::value_type>,
    class V = typename ::std::iterator_traits<Iterator>::value_type>
  inline Iterator upper_bound(Iterator b, Iterator e,  V& v, Less lt = Less()) {
      // compiler choose one.
      if (linear) while ((b != e) && !lt(v, *b)) ++b;
      else b = ::std::upper_bound(b, e, v, lt);
      return b;
  }

  ///  keep the unique keys in the input. primarily for reducing comm volume.   output is sorted.  equal operator forces comparison to Key
  template <typename V, typename Hs, typename Eq>
  void hash_unique(::std::vector<V> & input, bool sorted_input = false, const Hs & hash = Hs(), const Eq & equal = Eq()) {
    if (input.size() == 0) return;
    if (sorted_input) {  // already sorted, then just get the unique stuff and remove rest.
      auto end = ::std::unique(input.begin(), input.end(), equal);
      input.erase(end, input.end());
    } else {  // not sorted, so use an unordered_set to keep the first occurence.

      // sorting is SLOW and not scalable.  use unordered set instead.
      ::std::unordered_set<V, Hs, Eq > temp(input.begin(), input.end(), input.size());
      input.clear();
      input.assign(temp.begin(), temp.end());
    }
  }

  ///  keep the unique keys in the input. primarily for reducing comm volume.   output is sorted.  equal operator forces comparison to Key
  template <typename V, typename Comp, typename Eq>
  void sort_unique(::std::vector<V> & input, bool sorted_input = false, const Comp & comp = Comp(), const Eq & equal = Eq()) {
    if (input.size() == 0) return;
    if (!sorted_input) ::std::sort(input.begin(), input.end(), comp);
    // then just get the unique stuff and remove rest.
    auto end = ::std::unique(input.begin(), input.end(), equal);
    input.erase(end, input.end());
  }

}


#endif /* SRC_WIP_CONTAINER_UTILS_HPP_ */
