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

/**
 * @file    fsc_container_utils.hpp
 * @ingroup
 * @author  tpan
 * @brief
 * @details
 *
 */
#ifndef SRC_CONTAINERS_FSC_CONTAINER_UTILS_HPP_
#define SRC_CONTAINERS_FSC_CONTAINER_UTILS_HPP_

#include <iterator>  // iterator_traits
#include <unordered_set>
#include <algorithm>  // upper bound, unique, sort, etc.

#include "utils/benchmark_utils.hpp"

namespace fsc {

  struct TruePredicate {
      template <typename T>
      bool operator()(T const & x) const { return true; }
      template <typename Iter>
      bool operator()(Iter b, Iter e) const { return true; }
  };


  // identity transform
  template <typename Key>
  struct identity {
      inline Key operator()(Key const & x) const {
        return x;
      }
      template <typename VAL>
      inline ::std::pair<Key, VAL> operator()(std::pair<Key, VAL> & x) const {
        return x;
      }
  };


//template <typename T>
//struct IdentityTransform {
//    inline T operator()(T const& v) const {return v;};
//};


  template <typename Key, template <typename> class Hash, template <typename> class Transform>
  struct TransformedHash {
      Hash<Key> h;
      Transform<Key> trans;

      TransformedHash(Hash<Key> const & _hash = Hash<Key>(),
    		  Transform<Key> const &_trans = Transform<Key>()) : h(_hash), trans(_trans) {};

      inline uint64_t operator()(Key const& k) const {
        return h(trans(k));
      }
      template<typename V>
      inline uint64_t operator()(::std::pair<Key, V> const& x) const {
        return this->operator()(x.first);
      }
      template<typename V>
      inline uint64_t operator()(::std::pair<const Key, V> const& x) const {
        return this->operator()(x.first);
      }
  };


  template <typename Key, template <typename> class Predicate, template <typename> class Transform>
  struct TransformedPredicate {
      Predicate<Key> p;
      Transform<Key> trans;

      TransformedPredicate(Predicate<Key> const & _pred = Predicate<Key>(),
    		  Transform<Key> const &_trans = Transform<Key>()) : p(_pred), trans(_trans) {};

      inline bool operator()(Key const& k) const {
        return p(trans(k));
      }
      template<typename V>
      inline bool operator()(::std::pair<Key, V> const& x) const {
        return this->operator()(x.first);
      }
      template<typename V>
      inline bool operator()(::std::pair<const Key, V> const& x) const {
        return this->operator()(x.first);
      }
  };


  template <typename Key, template <typename> class Comparator, template <typename> class Transform>
  struct TransformedComparator {
      Comparator<Key> comp;
      Transform<Key> trans;

      TransformedComparator(Comparator<Key> const & _cmp = Comparator<Key>(),
    		  Transform<Key> const &_trans = Transform<Key>()) : comp(_cmp), trans(_trans) {};

      inline bool operator()(Key const & x, Key const & y) const {
        return comp(trans(x), trans(y));
      }
      template<typename V>
      inline bool operator()(::std::pair<Key, V> const & x, Key const & y) const {
        return this->operator()(x.first, y);
      }
      template<typename V>
      inline bool operator()(::std::pair<const Key, V> const & x, Key const & y) const {
        return this->operator()(x.first, y);
      }
      template<typename V>
      inline bool operator()(Key const & x, ::std::pair<Key, V> const & y) const {
        return this->operator()(x, y.first);
      }
      template<typename V>
      inline bool operator()(Key const & x, ::std::pair<const Key, V> const & y) const {
        return this->operator()(x, y.first);
      }
      template<typename V>
      inline bool operator()(::std::pair<Key, V> const & x, ::std::pair<Key, V> const & y) const {
        return this->operator()(x.first, y.first);
      }
      template<typename V>
      inline bool operator()(::std::pair<const Key, V> const & x, ::std::pair<const Key, V> const & y) const {
        return this->operator()(x.first, y.first);
      }
  };



  /// append to container via emplace.
  /// modified based on http://stackoverflow.com/questions/18724999/why-no-emplacement-iterators-in-c11-or-c14
  template<class Container>
  class back_emplace_iterator : public std::iterator< std::output_iterator_tag,
                                                      typename::std::iterator_traits<decltype(std::declval<Container>().begin())>::value_type,
                                                      void, void, void >
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

//  ///  keep the unique keys in the input. primarily for reducing comm volume.
//  ///  when sorted, ordering is unchanged.  when unsorted, ordering is not preserved.
//  ///  equal operator forces comparison to Key only (not pairs or tuples)
//  template <typename V, typename Hash, typename Equal>
//  struct unique_op {
//	  Hash h;
//	  Equal eq;
//
//	  unique_op(const Hash & hash = Hash(), const Equal & equal = Equal()) :
//	  	  h(hash), eq(equal) {};
//
//	  inline void process(::std::vector<V> & input, bool & sorted_input) {
//		  auto end = process(input.begin(), input.end(), sorted_input);
//		  input.erase(end, input.end());
//	  }
//
//	  template <typename Iter>
//	  inline Iter process(Iter first, Iter last, bool & sorted_input) {
//
//		if (first == last) return;
//		if (sorted_input) {  // already sorted, then just get the unique stuff and remove rest.
//		  return ::std::unique(first, last, eq);
//		} else {  // not sorted, so use an unordered_set to keep the first occurence.
//		  // sorting is SLOW and not scalable.  use unordered set instead.
//		  ::std::unordered_set<V, Hash, Equal> temp(first, last, std::distance(first, last), h, eq);
//		  return std::copy(temp.begin(), temp.end(), first);
//		}
//	  }
//  };
//
//  template <typename V, typename Hash, typename Equal>
//  void unique(::std::vector<V> & input, bool & sorted_input,
//                   const Hash & h = Hash(), const Equal & eq = Equal) {
//	  unique_op<V, Hash, Equal>(h, eq).process(input, sorted_input);
//  }
//
//
//  ///  keep the unique keys in the input. primarily for reducing comm volume.
//  /// output is SORTED.  when input is sorted, ordering is unchanged.
//  /// equal operator forces comparison to Key only (not pairs or tuples)
//  template <typename V, typename Less>
//  struct sort_op {
//	  Less ls;
//
//	  sort_op(const Less & less = Less()) :
//	  	  ls(less) {};
//
//	  inline void process(::std::vector<V> & input, bool & sorted_input) {
//		  process(input.begin(), input.end(), sorted_input);
//	  }
//
//	  template <typename Iter>
//	  inline Iter process(Iter first, Iter last, bool & sorted_input) {
//		  if (!sorted_input) {
//			  ::std::sort(first, last, ls);
//			  sorted_input = true;
//		  }
//		  return last;
//	  }
//  };
//
//  template <typename V, typename Less>
//  void sort(::std::vector<V> & input, bool & sorted_input,
//                   const Less & less = Less()) {
//	  sort_op<V, Less>(less).process(input, sorted_input);
//  }
//
//
//  ///  keep the unique keys in the input. primarily for reducing comm volume.
//  /// output is SORTED.  when input is sorted, ordering is unchanged.
//  /// equal operator forces comparison to Key only (not pairs or tuples)
//  template <typename V, typename Less, typename Equal>
//  struct sorted_unique_op {
//	  Less ls;
//	  Equal eq;
//
//	  sorted_unique_op(const Less & less = Less(), const Equal & equal = Equal()) :
//	  	  ls(less), eq(equal) {};
//
//	  inline void process(::std::vector<V> & input, bool & sorted_input) {
//		  auto end = process(input.begin(), input.end(), sorted_input);
//		  input.erase(end, input.end());
//	  }
//
//	  template <typename Iter>
//	  inline Iter process(Iter first, Iter last, bool & sorted_input) {
//
//		if (first == last) return last;
//
//	    if (!sorted_input) {
//	    	::std::sort(first, last, ls);
//		    sorted_input = true;
//	    }
//	    // then just get the unique stuff and remove rest.
//
//	    return ::std::unique(first, last, eq);
//	  }
//  };
//
//  template <typename V, typename Hash, typename Equal>
//  void sorted_unique(::std::vector<V> & input, bool & sorted_input,
//                   const Hash & h = Hash(), const Equal & eq = Equal) {
//	  sorted_unique_op<V, Hash, Equal>(h, eq).process(input, sorted_input);
//  }

  ///  keep the unique keys in the input. primarily for reducing comm volume.
  ///  when sorted, ordering is unchanged.  when unsorted, ordering is not preserved.
  ///  equal operator forces comparison to Key only (not pairs or tuples)
  template <typename V, typename Hash, typename Eq>
  void unique(::std::vector<V> & input, bool & sorted_input,
                   const Hash & hash = Hash(), const Eq & equal = Eq()) {
    if (input.size() == 0) return;
    if (sorted_input) {  // already sorted, then just get the unique stuff and remove rest.
      auto end = ::std::unique(input.begin(), input.end(), equal);
      input.erase(end, input.end());
    } else {  // not sorted, so use an unordered_set to keep the first occurence.

      // sorting is SLOW and not scalable.  use unordered set instead.
      ::std::unordered_set<V, Hash, Eq> temp(input.begin(), input.end(), input.size(), hash, equal);
      input.clear();
      input.assign(temp.begin(), temp.end());
    }
  }


  ///  keep the unique keys in the input. primarily for reducing comm volume.
  /// output is SORTED.  when input is sorted, ordering is unchanged.
  /// equal operator forces comparison to Key only (not pairs or tuples)
  template <typename V, typename Less>
  void sort(::std::vector<V> & input, bool & sorted_input,
                   const Less & less = Less()) {
    if (!sorted_input) ::std::sort(input.begin(), input.end(), less);

    sorted_input = true;
  }


  ///  keep the unique keys in the input. primarily for reducing comm volume.
  /// output is SORTED.  when input is sorted, ordering is unchanged.
  /// equal operator forces comparison to Key only (not pairs or tuples)
  template <typename V, typename Less, typename Eq>
  void sorted_unique(::std::vector<V> & input, bool & sorted_input,
                   const Less & less = Less(), const Eq & equal = Eq()) {
    if (input.size() == 0) return;
    if (!sorted_input) ::std::sort(input.begin(), input.end(), less);
    // then just get the unique stuff and remove rest.
    auto end = ::std::unique(input.begin(), input.end(), equal);
    input.erase(end, input.end());

    sorted_input = true;
  }

  /// keep the unique entries within each bucket.  complexity is b * O(N/b), where b is the bucket size, and O(N/b) is complexity of inserting into and copying from set.
  /// when used within bucket, scales with O(N/b), not with b.  this is as good as it gets wrt complexity.
  /// sortedness is MAINTAINED within buckets
  template<typename T, typename count_t, typename Hash, typename Eq>
  void bucket_unique(std::vector<T>& input, std::vector<count_t> &send_counts, bool & sorted_input,
                          const Hash & hash = Hash(), const Eq & equal = Eq()) {

    auto newstart = input.begin();
    auto newend = newstart;
    auto start = input.begin();
    auto end = start;


    if (sorted_input) {

      for (size_t i = 0; i < send_counts.size(); ++i) {
        end = start + send_counts[i];

        if (i == 0)
          newend = ::std::unique(start, end, equal);
        else
          newend = ::std::unique_copy(start, end, newstart, equal);

        send_counts[i] = ::std::distance(newstart, newend);

        start = end;
        newstart = newend;
      }

    } else {

      count_t max = *(::std::max_element(send_counts.begin(), send_counts.end()));
      ::std::unordered_set<T, Hash, Eq> set(max, hash, equal);

      for (size_t i = 0; i < send_counts.size(); ++i) {
        end = start + send_counts[i];

        // sorting is SLOW and not scalable.  use unordered set instead.
        // unordered_set for large data is memory intensive.  depending on use, bucket per processor first.
        set.clear();
        set.insert(start, end);
        newend = ::std::copy(set.begin(), set.end(), newstart);

        send_counts[i] = ::std::distance(newstart, newend);

        start = end;
        newstart = newend;
      }

    }

    // compact.
    input.erase(newend, input.end());
  }


  /// keep the unique entries within each bucket.  complexity is b * O(N/b), where b is the bucket size, and O(N/b) is complexity of inserting into and copying from set.
  /// when used within bucket, scales with O(N/b), not with b.  this is as good as it gets wrt complexity.
  /// sortedness is MAINTAINED within buckets
  template<typename T, typename count_t, typename Less>
  void bucket_sort(std::vector<T>& input, std::vector<count_t> &send_counts, bool & sorted_input,
                          const Less & less = Less()) {

    if (!sorted_input) {
      auto newstart = input.begin();
      auto newend = newstart;
      auto start = input.begin();
      auto end = start;

      for (size_t i = 0; i < send_counts.size(); ++i) {
        end = start + send_counts[i];

        ::std::sort(start, end, less);

        start = end;
      }
      sorted_input = true;
    }

  }

  /// keep the unique entries within each bucket.  complexity is b * O(N/b), where b is the bucket size, and O(N/b) is complexity of inserting into and copying from set.
  /// when used within bucket, scales with O(N/b), not with b.  this is as good as it gets wrt complexity.
  /// sortedness is MAINTAINED within buckets
  template<typename T, typename count_t, typename Less, typename Eq>
  void bucket_sorted_unique(std::vector<T>& input, std::vector<count_t> &send_counts, bool & sorted_input,
                          const Less & less = Less(), const Eq & equal = Eq()) {

    auto newstart = input.begin();
    auto newend = newstart;
    auto start = input.begin();
    auto end = start;

    for (size_t i = 0; i < send_counts.size(); ++i) {
      end = start + send_counts[i];

      if (!sorted_input) {
        ::std::sort(start, end, less);
      }

      if (i == 0)
        newend = ::std::unique(start, end, equal);
      else
        newend = ::std::unique_copy(start, end, newstart, equal);

      send_counts[i] = ::std::distance(newstart, newend);

      start = end;
      newstart = newend;
    }


    // compact.
    input.erase(newend, input.end());

    sorted_input = true;
  }


  /// keep the unique entries within each bucket.  complexity is b * O(N/b), where b is the bucket size, and O(N/b) is complexity of inserting into and copying from set.
  /// when used within bucket, scales with O(N/b), not with b.  this is as good as it gets wrt complexity.
  /// sortedness is MAINTAINED within buckets
  template<typename T, typename count_t, typename Reduc>
  void bucket_reduce(std::vector<T>& input, std::vector<count_t> &send_counts, bool & sorted_input,
                          const Reduc & reducer = Reduc()) {

	    auto newstart = input.begin();
	    auto newend = newstart;
	    auto start = input.begin();
	    auto end = start;

	    bool sorted = sorted_input;

	    for (size_t i = 0; i < send_counts.size(); ++i) {

	      end = start + send_counts[i];

	      sorted = sorted_input;
	      if (i == 0)
	        newend = reducer(start, end, start, sorted);
	      else
	        newend = reducer(start, end, newstart, sorted);

	      send_counts[i] = ::std::distance(newstart, newend);

	      start = end;
	      newstart = newend;
	    }


	    // compact.
	    input.erase(newend, input.end());

	    sorted_input = true;
  }


  // finds splitter positions for the given splitters (`map`)
  // in the sorted buffer, sorted according to the `comp` comparator, and maps to communicator rank order.
  template<typename T, typename Key, typename Comparator>
  std::vector<size_t> sorted_bucketing(std::vector<T>& buffer, ::std::vector<::std::pair<Key, int> > map, Comparator comp, int comm_size) {
    //== track the send count
    std::vector<size_t> send_counts(comm_size, 0);

    if (buffer.size() == 0) return send_counts;

    // check the map a little.
    {
      assert(map.size() < static_cast<size_t>(comm_size));
      std::vector<int> ranks(map.size(), 0);
      std::transform(map.begin(), map.end(), ranks.begin(), [](::std::pair<Key, int> const & x){
        return x.second;
      });
      // check min max values are within comm rank values.
      assert(*(std::min_element(ranks.begin(), ranks.end())) >= 0);
      assert(*(std::max_element(ranks.begin(), ranks.end())) < comm_size);
      // check ranks are sorted (since it was generated that way).
      assert(std::is_sorted(ranks.begin(), ranks.end()));
      // check no duplicates.
      assert(std::adjacent_find(ranks.begin(), ranks.end()) == ranks.end());
    }

    //== since both sorted, search map entry in buffer - mlog(n).  the otherway around is nlog(m).
    // note that map contains splitters.  range defined as [map[i], map[i+1]), so when searching in buffer using entry in map, search via lower_bound
    auto b = buffer.begin();
    auto e = b;

    // approach below allows skipping processes, which in turn minimizes shifting.
    // however, the ranks should still be organized in monotonically increasing order.
    // since the buffer is organized by increasing comm rank.

    // if mlog(n) is larger than n, then linear would be faster
    bool linear = (map.size() * std::log(buffer.size()) ) > buffer.size();

    if (linear)
      for (size_t i = 0; i < map.size(); ++i) {
        e = ::fsc::lower_bound<true>(b, buffer.end(), map[i].first, comp);
        send_counts[map[i].second] = ::std::distance(b, e);
        b = e;
      }
    else
      for (size_t i = 0; i < map.size(); ++i) {
        e = ::fsc::lower_bound<false>(b, buffer.end(), map[i].first, comp);
        send_counts[map[i].second] = ::std::distance(b, e);
        b = e;
      }

    // last splitter's rank + 1 gets all the remainder.
    send_counts[map.back().second + 1] = ::std::distance(b, buffer.end());

    return send_counts;
  }


}  // namespace fsc


#endif /* SRC_CONTAINERS_FSC_CONTAINER_UTILS_HPP_ */
