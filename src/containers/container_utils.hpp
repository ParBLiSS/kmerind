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
 * @file    vector_utils.hpp
 * @ingroup
 * @author  tpan
 * @brief
 * @details
 *
 */
#ifndef SRC_CONTAINERS_CONTAINER_UTILS_HPP_
#define SRC_CONTAINERS_CONTAINER_UTILS_HPP_

#include <iterator>  // iterator_traits
#include <unordered_set>
#include <algorithm>  // upper bound, unique, sort, etc.

#include "utils/benchmark_utils.hpp"

namespace fsc {

  template <typename Key, typename Hash, template <typename> class Transform>
  struct TransformedHash {
      Hash h;
      Transform<Key> trans;

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



  template <typename T>
  struct IdentityTransform {
      T operator()(T const& v) {return v;};
  };

  template <typename Key, typename Comparator, template <typename> class Transform>
  struct TransformedComparator {
      Comparator comp;
      Transform<Key> trans;

      TransformedComparator(Comparator const & _cmp = Comparator(), Transform<Key> const _trans = Transform<Key>()) : comp(_cmp), trans(_trans) {};

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

  template <typename Key, typename Comparator>
  struct TransformedComparator<Key, Comparator, IdentityTransform> {
      Comparator comp;

      TransformedComparator(Comparator const & cmp = Comparator()) : comp(cmp) {};

      inline bool operator()(Key const & x, Key const & y) const {
        return comp(x, y);
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


  template <typename Key, typename Less = ::std::less<Key> >
  struct Greater {
      Less lt;
      inline bool operator()(Key const &x, Key const &y) const {
        return lt(y, x);
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


  // finds splitter positions for the given splitters (`map`)
  // in the sorted buffer, sorted according to the `comp` comparator.
  template<typename T, typename Key, typename Comp>
  std::vector<size_t> get_bucket_sizes(std::vector<T>& buffer, ::std::vector<::std::pair<Key, int> > map, Comp comp) {
    //== track the send count
    std::vector<size_t> send_counts(map.size() + 1, 0);

    if (buffer.size() == 0) return send_counts;

    //== since both sorted, search map entry in buffer - mlog(n).  the otherway around is nlog(m).
    // note that map contains splitters.  range defined as [map[i], map[i+1]), so when searching in buffer using entry in map, search via lower_bound
    auto b = buffer.begin();
    auto e = b;

    // if mlog(n) is larger than n, then linear would be faster
    bool linear = (map.size() * std::log(buffer.size()) ) > buffer.size();

    if (linear)
      for (size_t i = 0; i < map.size(); ++i) {
        e = ::fsc::lower_bound<true>(b, buffer.end(), map[i].first, comp);
        send_counts[i] = ::std::distance(b, e);
        b = e;
      }
    else
      for (size_t i = 0; i < map.size(); ++i) {
        e = ::fsc::lower_bound<false>(b, buffer.end(), map[i].first, comp);
        send_counts[i] = ::std::distance(b, e);
        b = e;
      }

    // last 1
    send_counts[map.size()] = ::std::distance(b, buffer.end());

    return send_counts;
  }


}  // namespace fsc


namespace dsc {


  // =============== convenience functions for distribution of vector via all to all and a rank mapping function

  /**
   * @brief       distribute.  speed version.  no guarantee of output ordering, but actually is same.
   * @param[IN/OUT] vals    vals to distribute.  sortedness is NOT kept because of inplace bucketing.
   * @param[IN/OUT] sorted_input    indicates if input is sorted.  and whether each bucket is sorted.
   *
   * @return received counts.
   */
  template <typename V, typename ToRank>
  ::std::vector<size_t> distribute(::std::vector<V>& vals, ToRank const & to_rank, bool & sorted_input,
                                   mxx::comm const &_comm) {
    BL_BENCH_INIT(distribute);

      // speed over mem use.  mxx all2allv already has to double memory usage. same as stable distribute.

      BL_BENCH_START(distribute);
      // distribute (communication part)
      std::vector<size_t> send_counts = mxx::bucketing(vals, to_rank, _comm.size());
      BL_BENCH_END(distribute, "bucket", vals.size());

      // distribute (communication part)
      BL_BENCH_COLLECTIVE_START(distribute, "a2a", _comm);
      vals = mxx::all2allv(vals, send_counts, _comm);
      BL_BENCH_END(distribute, "a2a", vals.size());

      BL_BENCH_START(distribute);
      std::vector<size_t> recv_counts= mxx::all2all(send_counts, _comm);
      BL_BENCH_END(distribute, "a2a_counts", vals.size());

      BL_BENCH_REPORT_MPI_NAMED(distribute, "map_base:distribute", _comm);

      return recv_counts;
  }


    /**
     * @brief       distribute and ensure uniqueness within each bucket.  speed version.  no guarantee of output ordering.
     * @param[IN/OUT] keys    keys to distribute. sortedness is NOT kept because of inplace bucketing.,
     * @param[IN/OUT] sorted_input    indicates if input is sorted.  and whether each bucket is sorted.
     * @return received counts
     */
    template <typename V, typename ToRank, typename Hash, typename Eq>
    ::std::vector<size_t> distribute_unique(::std::vector<V>& vals, ToRank const & to_rank, bool & sorted_input,
                                            mxx::comm const &_comm, const Hash & hash = Hash(), const Eq & equal = Eq()) {

      BL_BENCH_INIT(distribute);
        // go for speed.   mxx::bucketing uses extra copy.  okay, since all2all also does.

        BL_BENCH_START(distribute);
        // distribute (communication part)
        std::vector<size_t> send_counts = mxx::bucketing(vals, to_rank, _comm.size());
        BL_BENCH_END(distribute, "bucket", vals.size());

        // using set is okay.
        BL_BENCH_START(distribute);
        // distribute (communication part)
        ::fsc::bucket_unique(vals, send_counts, sorted_input, hash, equal);
        BL_BENCH_END(distribute, "unique", vals.size());

        // distribute (communication part)
        BL_BENCH_COLLECTIVE_START(distribute, "a2a", _comm);
        vals = mxx::all2allv(vals, send_counts, _comm);
        BL_BENCH_END(distribute, "a2a", vals.size());

        BL_BENCH_START(distribute);
        std::vector<size_t> recv_counts= mxx::all2all(send_counts, _comm);
        BL_BENCH_END(distribute, "a2a_counts", vals.size());

        BL_BENCH_REPORT_MPI_NAMED(distribute, "map_base:distribute", _comm);


        return recv_counts;
    }


    /**
     * @brief       distribute.  order preserving version.  guarantees output ordering.
     * @param[IN/OUT] vals    vals to distribute.  sortedness is KEPT within each bucket
     * @param[IN/OUT] sorted_input    indicates if input is sorted.  and whether each bucket is sorted.
     *
     * @return received counts.
     */
    template <typename V, typename ToRank, typename Less>
    ::std::vector<size_t> distribute_sorted(::std::vector<V>& vals, ToRank const & to_rank, bool & sorted_input,
                                            mxx::comm const &_comm, const Less & less = Less()) {
      BL_BENCH_INIT(distribute);

      // ordering over speed/mem use.

        BL_BENCH_START(distribute);
        // distribute (communication part)
        std::vector<size_t> send_counts = mxx::bucketing(vals, to_rank, _comm.size());
        BL_BENCH_END(distribute, "bucket", vals.size());

        // using set is okay.
        BL_BENCH_START(distribute);
        // distribute (communication part)
        ::fsc::bucket_sort(vals, send_counts, sorted_input, less);
        BL_BENCH_END(distribute, "sort", vals.size());

        // distribute (communication part)
        BL_BENCH_COLLECTIVE_START(distribute, "a2a", _comm);
        vals = mxx::all2allv(vals, send_counts, _comm);
        BL_BENCH_END(distribute, "a2a", vals.size());

        BL_BENCH_START(distribute);
        std::vector<size_t> recv_counts= mxx::all2all(send_counts, _comm);
        BL_BENCH_END(distribute, "a2a_counts", vals.size());

        BL_BENCH_REPORT_MPI_NAMED(distribute, "map_base:distribute", _comm);


        return recv_counts;
    }

    /**
     * @brief       distribute.  order preserving version.  guarantees output ordering.
     *
     * @param[IN/OUT] keys    keys to distribute. sortedness is KEPT.,
     * @param[IN/OUT] sorted_input    indicates if input is sorted.  and whether each bucket is sorted.
     * @return received counts
     */
    template <typename V, typename ToRank, typename Less, typename Eq>
    ::std::vector<size_t> distribute_sorted_unique(::std::vector<V>& vals, ToRank const & to_rank, bool & sorted_input,
                                                   mxx::comm const &_comm, const Less & less = Less(), const Eq & equal = Eq()) {
      BL_BENCH_INIT(distribute);
      // ordering over speed/mem use.

        BL_BENCH_START(distribute);
        // distribute (communication part)
        std::vector<size_t> send_counts = mxx::bucketing(vals, to_rank, _comm.size());
        BL_BENCH_END(distribute, "bucket", vals.size());

        BL_BENCH_START(distribute);
        // distribute (communication part)
        ::fsc::bucket_sorted_unique(vals, send_counts, sorted_input, less, equal);
        BL_BENCH_END(distribute, "sort_unique", vals.size());


        // distribute (communication part)
        BL_BENCH_COLLECTIVE_START(distribute, "a2a", _comm);
        vals = mxx::all2allv(vals, send_counts, _comm);
        BL_BENCH_END(distribute, "a2a", vals.size());

        BL_BENCH_START(distribute);
        std::vector<size_t> recv_counts= mxx::all2all(send_counts, _comm);
        BL_BENCH_END(distribute, "a2a_counts", vals.size());

        BL_BENCH_REPORT_MPI_NAMED(distribute, "map_base:distribute", _comm);


        return recv_counts;
	}






}  // namespace dsc


#endif /* SRC_CONTAINERS_CONTAINER_UTILS_HPP_ */
