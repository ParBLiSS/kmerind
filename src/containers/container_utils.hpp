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

//template <typename T>
//struct IdentityTransform {
//    inline T operator()(T const& v) const {return v;};
//};


  template <typename Key, template <typename> class Hash, template <typename> class Transform>
  struct TransformedHash {
      Hash<Key> h;
      Transform<Key> trans;

      TransformedHash(Hash<Key> const & _hash = Hash<Key>(),
    		  Transform<Key> const _trans = Transform<Key>()) : h(_hash), trans(_trans) {};

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


//  template <typename Key, template <typename> class Hash>
//  struct TransformedHash<Key, Hash, IdentityTransform> {
//      Hash<Key> h;
//
//      TransformedHash(Hash<Key> const & _hash = Hash<Key>(),
//    		  IdentityTransform<Key> const &_trans = IdentityTransform<Key>()) : h(_hash) {};
//
//      inline uint64_t operator()(Key const& k) const {
//        return h(k);
//      }
//      template<typename V>
//      inline uint64_t operator()(::std::pair<Key, V> const& x) const {
//        return this->operator()(x.first);
//      }
//      template<typename V>
//      inline uint64_t operator()(::std::pair<const Key, V> const& x) const {
//        return this->operator()(x.first);
//      }
//  };
//

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

//  template <typename Key, template <typename> class Comparator>
//  struct TransformedComparator<Key, Comparator, IdentityTransform> {
//      Comparator<Key> comp;
//
//      TransformedComparator(Comparator<Key> const & cmp = Comparator<Key>(),
//    		  IdentityTransform<Key> const &_trans = IdentityTransform<Key>()) : comp(cmp) {};
//
//      inline bool operator()(Key const & x, Key const & y) const {
//        return comp(x, y);
//      }
//      template<typename V>
//      inline bool operator()(::std::pair<Key, V> const & x, Key const & y) const {
//        return this->operator()(x.first, y);
//      }
//      template<typename V>
//      inline bool operator()(::std::pair<const Key, V> const & x, Key const & y) const {
//        return this->operator()(x.first, y);
//      }
//      template<typename V>
//      inline bool operator()(Key const & x, ::std::pair<Key, V> const & y) const {
//        return this->operator()(x, y.first);
//      }
//      template<typename V>
//      inline bool operator()(Key const & x, ::std::pair<const Key, V> const & y) const {
//        return this->operator()(x, y.first);
//      }
//      template<typename V>
//      inline bool operator()(::std::pair<Key, V> const & x, ::std::pair<Key, V> const & y) const {
//        return this->operator()(x.first, y.first);
//      }
//      template<typename V>
//      inline bool operator()(::std::pair<const Key, V> const & x, ::std::pair<const Key, V> const & y) const {
//        return this->operator()(x.first, y.first);
//      }
//  };


//  template <typename Key, typename Less = ::std::less<Key> >
//  struct Greater {
//      Less lt;
//      inline bool operator()(Key const &x, Key const &y) const {
//        return lt(y, x);
//      }
//  };



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


namespace dsc {

	template <typename CONTAINER>
	bool empty(CONTAINER const & c, mxx::comm const & comm) {
		bool local_empty = (c.size() == 0);
	  if (comm.size() == 1) {
			if (local_empty) printf("rank 0/1 input is EMPTY.\n");
		return local_empty;
	} else { // all reduce
		local_empty =  mxx::all_of(local_empty, comm);
		if (local_empty && comm.rank() == 0) printf("rank ALL/%d inputs are all EMPTY.\n", comm.size());
		return local_empty;
	  }
	}


  // =============== convenience functions for distribution of vector via all to all and a rank mapping function
	// TODO: make this cleaner...


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
     *				to_rank uses distribution hash and transform.  unique should use the storage hash and transform
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

        BL_BENCH_REPORT_MPI_NAMED(distribute, "map_base:distribute_unique", _comm);


        return recv_counts;
    }



    /**
     * @brief       distribute.  order preserving version.  guarantees output ordering.
     *				to_rank uses distribution hash and transform.  unique should use the storage hash and transform
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

        BL_BENCH_REPORT_MPI_NAMED(distribute, "map_base:distribute_sorted", _comm);


        return recv_counts;
    }

    /**
     * @brief       distribute.  order preserving version.  guarantees output ordering.
     *				to_rank uses distribution hash and transform.  unique should use the storage hash and transform
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

        BL_BENCH_REPORT_MPI_NAMED(distribute, "map_base:distribute_sorted_unique", _comm);


        return recv_counts;
	}


    /**
     * @brief       distribute and ensure uniqueness within each bucket.  speed version.  no guarantee of output ordering.
     *				to_rank uses distribution hash and transform.  unique should use the storage hash and transform
     * @param[IN/OUT] keys    keys to distribute. sortedness is NOT kept because of inplace bucketing.,
     * @param[IN/OUT] sorted_input    indicates if input is sorted.  and whether each bucket is sorted.
     * @return received counts
     */
    template <typename V, typename ToRank, typename Reduc>
    ::std::vector<size_t> distribute_reduce(::std::vector<V>& vals, ToRank const & to_rank, bool & sorted_input,
                                            mxx::comm const &_comm, const Reduc & reducer = Reduc()) {

      BL_BENCH_INIT(distribute);
        // go for speed.   mxx::bucketing uses extra copy.  okay, since all2all also does.

  	if (_comm.rank() == 0) printf("start\n");  fflush(stdout);


        BL_BENCH_START(distribute);
        // distribute (communication part)
        std::vector<size_t> send_counts = mxx::bucketing_inplace(vals, to_rank, _comm.size());
        BL_BENCH_END(distribute, "bucket", vals.size());

    	if (_comm.rank() == 0) printf("bucket\n");  fflush(stdout);

    	double mean, stdev;
    	for (auto it = send_counts.begin(), max = send_counts.end(); it != max; ++it) {
    		mean += *it;
    		stdev += (*it) * (*it);
    	}
    	mean /= _comm.size();
    	stdev -= sqrt((stdev / _comm.size()) - (mean * mean));
    	if (_comm.rank() == 0) printf("mean: %f, stdev %f\n", mean, stdev);


        // using set is okay.
        BL_BENCH_START(distribute);
        // distribute (communication part)
        ::fsc::bucket_reduce(vals, send_counts, sorted_input, reducer);
        BL_BENCH_END(distribute, "reduce", vals.size());

    	if (_comm.rank() == 0) printf("reduce\n");  fflush(stdout);


        // distribute (communication part)
        BL_BENCH_COLLECTIVE_START(distribute, "a2a", _comm);
        vals = mxx::all2allv(vals, send_counts, _comm);
        BL_BENCH_END(distribute, "a2a", vals.size());


    	if (_comm.rank() == 0) printf("a2a\n");  fflush(stdout);

        BL_BENCH_START(distribute);
        std::vector<size_t> recv_counts= mxx::all2all(send_counts, _comm);
        BL_BENCH_END(distribute, "a2a_counts", vals.size());

    	if (_comm.rank() == 0) printf("a2acounts\n");  fflush(stdout);

        BL_BENCH_REPORT_MPI_NAMED(distribute, "map_base:distribute_reduce", _comm);


        return recv_counts;
    }

    /**
     * @brief samples a vector to be evenly split, sorting the samples to get splitters.
     * @note  goal is to make the sample sorting process parallel.
     * 		NOTE that duplicates are removed.
     * 		NOTE sampling uses random number generator.
     *
     *      random sample p.
	 *		local sort the p samples.
	 *		use all2allv with overlapping range to send 2 to each proc.
	 *			proc 1:  1,2
	 *			proc 2:  2,3
	 *			...
	 *			proc p:  p
	 *		local sort the 2p samples and reduce.
	 *		pick the middle, first p-1 procs.
	 *
	 *		allgatherv to form splitters.
	 *
	 *		sort splitters, and reduce
	 *
	 *		use splitters to split (get send_counts
	 *
	 *		all2allv
	 *
	 *		sort and reduce.
	 *
	 *		block decomp.
	 *		get 1st of each as final splitter.
	 *		allgather
	 *  NOT WELL TESTED.
	 * @param vals   vector (distributed on multiple procs. that needs to have its splitters generated.
     */
    template <typename V, typename Less, typename Eq>
    ::std::vector<V> sample_for_splitters(::std::vector<V>& vals, bool & sorted_input,
                    mxx::comm const &_comm, const Less & less = Less(), const Eq & equal = Eq()) {

    	// random number generator setup.
		std::default_random_engine generator;
		std::uniform_int_distribution<size_t> distribution(0, (vals.size() - 1));
		auto long_rand = std::bind ( distribution, generator );

		// random sample p.
		std::vector< V > samples;
		for (int i = 0; i < _comm.size(); ++i) {
			samples.emplace_back(vals[long_rand()]);
		}

//		// local sort p
//		std::sort(samples.begin(), samples.end(), less);
//
//		// shuffle, 2 items each
//		std::vector<size_t> send_counts(_comm.size(), 2);
//		send_counts[_comm.size() - 1] = 1;
//		std::vector<size_t> send_displs(_comm.size(), 0);
//		std::iota(send_displs.begin(), send_displs.end(), 0);
//
//		std::vector<size_t> recv_counts = ::mxx::all2all(send_counts, _comm);
//		std::vector<size_t> recv_displs = ::mxx::impl::get_displacements(send_counts);
//		size_t recv_total = recv_dspls.back() + recv_counts.back();
//		::std::vector<Key> moved_samples(recv_total);
//		mxx::all2allv(&(samples[0]), send_counts, send_displs,
//				&(moved_samples[0]), recv_counts, recv_displs,
//				, _comm);
		::std::vector<V> moved_samples = mxx::all2all(samples, _comm);

		// local sort again, with reduction this time.
		std::sort(moved_samples.begin(), moved_samples.end(), less);
		auto mend = std::unique(moved_samples.begin(), moved_samples.end(), equal);
		moved_samples.erase(mend, moved_samples.end());

//		// at this point, we expect overlaps, with worst case
//		// being p(p-1) for the split.
		std::vector< V > splitters;
//		if (_comm.rank() < (_comm.size() - 1)) splitters.emplace_back(moved_samples.back());
		if (_comm.rank() > 0) splitters.emplace_back(moved_samples.front());
		splitters = ::mxx::allgatherv(splitters, _comm);

		// allgather splitters and sort it
		std::sort(splitters.begin(), splitters.end(), less);
		auto s_end = std::unique(splitters.begin(), splitters.end(), equal);
		splitters.erase(s_end, splitters.end());

		// compute the send counts for samples - no identical splitters, so won't see things crossing boundaries.
		std::vector<size_t> send_counts =
				::mxx::impl::split(samples.begin(), samples.end(), less, splitters, _comm);

		// and all2all the samples (so there are no repeated values between nodes)
		samples = mxx::all2allv(samples, send_counts, _comm);

		// sort and reduce the split samples
		std::sort(samples.begin(), samples.end(), less);
		mend = std::unique(samples.begin(), samples.end(), equal);
		samples.erase(mend, samples.end());

		// block decompose again
		moved_samples = mxx::stable_distribute(samples, _comm);

		splitters.clear();
		if (_comm.rank() > 0) splitters.emplace_back(moved_samples.front());
		splitters = mxx::allgatherv(splitters, _comm);

		return splitters;

    }

}  // namespace dsc


#endif /* SRC_CONTAINERS_CONTAINER_UTILS_HPP_ */
