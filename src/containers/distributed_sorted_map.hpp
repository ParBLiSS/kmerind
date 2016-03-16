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
 * @file    distributed_sorted_map.hpp
 * @ingroup index
 * @author  Tony Pan <tpan7@gatech.edu>
 * @brief   Implements the distributed_multimap, distributed map, and distributed_reduction_map
 *          data structures.
 *
 *          implementation is sort-based (load balanced).  assumption is that map is built once.
 *          distribution is via sort.  local storage is via hash.
 *
 *          for now, input and output via local vectors.
 *          (later, allow remote vectors,
 *            which can have remote ranges  (all to all to "sort" remote ranges to src proc,
 *            then src proc compute target procs for each elements in the remote ranges,
 *            communicate remote ranges to target procs.  target proc can then materialize the data.
 *            may not be efficient though if we don't have local spatial coherence..
 *          )
 *
 *          most create-find-delete operations support remote filtering via predicates.
 *          most create-find-delete oeprations support remote transformation.
 *
 *          signature of predicate is bool pred(T&).  if predicate needs to access the local map, it should be done via its constructor.
 *
 */

#ifndef BLISS_DISTRIBUTED_SORTED_MAP_HPP
#define BLISS_DISTRIBUTED_SORTED_MAP_HPP


#include <vector>  // local storage hash table
#include <utility> 			  // for std::pair

//#include <sparsehash/dense_hash_map>  // local storage hash table, google dense hash map.
#include <functional> 		// for std::function and std::hash
#include <algorithm> 		// for sort, stable_sort, unique, is_sorted
#include <iterator>  // advance, distance

#include <cstdint>  // for uint8, etc.


#include <mxx/collective.hpp>
#include <mxx/reduction.hpp>
#include <mxx/sort.hpp>
#include "utils/benchmark_utils.hpp"  // for timing.
#include "utils/logging.h"
#include "containers/distributed_map_base.hpp"


namespace dsc  // distributed std container
{

  /**
   * @brief  distributed unordered map following std unordered map's interface.
   * @details   This class is modeled after the std::unordered_map.
   *         it has as much of the same methods of std::unordered_map as possible.  however, all methods consider the fact
   *         that the data are in distributed memory space, so to access the data, "communication" is needed.  Also since we
   *         are working with 'distributed' data, batched operations are preferred.
   *
   *         Note that "communication" is a weak concept here meaning that we are accessing a different local container.
   *         as such, communicator may be defined for MPI, UPC, OpenMP, etc.
   *
   *         This allows the possibility of using distributed unordered map as local storage for coarser grain distributed container.
   *
   *         Note that communicator requires a mapping strategy between a key and the target processor/thread/partition.  The mapping
   *         may be done using a hash, similar to the local distributed unordered map, or it may be done via sorting/lookup or other mapping
   *         mechanisms.  The choice may be constrained by the communication approach, e.g. global sorting  does not work well with
   *         incremental async communication
   *
   *  this class and its subclasses rely on 2 hash function for data distribution and 1 equal and 1 less comparators.  these are specified
   *  as template parameters.  the less comparator is for preprocessing queries and inserts.  keys (e.g.) kmers are transformed before they
   *  are hashed/compared.
   *    an alternative approach is to hold only canonical keys in the map.
   *    yet another alternative approach is to perform 2 queries for every key.- 2x computation but communication is spread out.
   *
   *  note: KeyTransform is applied before Hash, and Equal operators.  These operators should have NO KNOWLEDGE of any transform applied, including kmolecule to kmer mapping.
   *
   *  any operation that uses sort is not going to scale well.  this includes "retain_unique_key, retain_unique_tuple, local_reduction"...
   *  to reduce this, we can try by using a hash set instead.  http://www.vldb.org/pvldb/2/vldb09-257.pdf, http://www.vldb.org/pvldb/vol7/p85-balkesen.pdf
   *
   *  conditional version of insert/erase/find/count supports predicate that operate on INTERMEDIATE RESULTS.  input (key,key-value pair) can be pre-filtered.
   *    output (query result, e.g.) can be post filtered (and optionally reduce comm volume).
   *    intermediate results (such as counting in multimap only if T has certain value) can only be filtered at during local_operation.
   *
   * key to proc assignment can be done as hash or splitters in sorted range.
   * tuples can be sotred in hash table or as sorted array.
   *   hash-hash combination works
   *  sort-sort combination works as well
   *  hash-sort combination can work.  advantage is in range query.
   *  sort-hash combination would be expensive for updating splitters, and may not make sense.
   *
   * This version is the sort-sort.
   *
   * there exists a few problems with this approach.  1. we NEED to accumulate then sort.  so incremental update to sorted list is costly.
   * 	2. perfect load balance is not possible - need to send to multiple targets but we don't have enough info in the key_to_rank map.
   * 	3. imperfect load balance suffers the same issue as hash.
   *
   *
   * note: 1. incremental insert searching unique for each block of entry is less efficient than insert all then unique once at the end.  (if merge used, m^2n merge complexity)
   * 				insert : direct insert into local container
   * 				rehash : heavy lifting
   * 				find/cnt   : unique input, distribute, search, return results.
   * 				erase  : unique input, distribute, delete.  user can then call rehash with sorted set to true: (multimap: balance, sample and redistribute, keep pivots) (map: balance, get pivots)
   *       2. rehash - redistribute carefully.  map: should balance, sort(unique first, sample and redistribute, k-way merge, unique, balance (entries are unique so don't span multiple procs)), get pivots).
   *       										multimap:  balance, mxx::sort(sort, sample and redistribute, k-way merge), sample and redistribute (due to nonunique entries), get pivots.
   *       3. sampling:  should it be based on nyquist, i.e. 2p-1?
   *
   * predicate : operates on container content.  query can be filtered before hand, and results can be filtered after (on source proc, or on dest proc).  container content should be filtered during.
   *
   * @tparam Key
   * @tparam T
   * @tparam Container  default to vector, requiring 2 template params.
   * @tparam Comm   default to mpi_collective_communicator       communicator for global communication. may hash or sort.
   * @tparam KeyTransform   transform function for the key.  can supply identity.  requires a single template argument (Key).  useful for mapping kmolecule to kmer.
   * @tparam Hash   hash function for local and distribution.  requires a template arugment (Key), and a bool (prefix, chooses the MSBs of hash instead of LSBs)
   * @tparam Equal   default to ::std::equal_to<Key>   equal function for the local storage.
   * @tparam Alloc  default to ::std::allocator< ::std::pair<Key, T> >    allocator for local storage.
   */
  template<typename Key, typename T,
  template <typename, typename> class Container,
  class Comm,
  template <typename> class KeyTransform,
  class Less = ::std::less<Key>,
  class Equal = ::std::equal_to<Key>,
  class Alloc = ::std::allocator< ::std::pair<Key, T> >
  >  class sorted_map_base : public ::dsc::map_base<Key, T, Comm, KeyTransform, Less, Equal, Alloc> {

      // TODO: this should not be public but clang complains.
    protected:
      using Base = ::dsc::map_base<Key, T, Comm, KeyTransform, Less, Equal, Alloc>;

      // splitters.  there should be p-1 entries.

      // uses sorted lookup table to map to target ranks.  assume we store in a vector a pair:  <first kmer on proc, proc id>.  then we can use lower_bound to search.
      // note that this makes in-between values go with the larger proc id.
      struct KeyToRank {
          ::std::vector<::std::pair<Key, int> > map;  // the splitters need to support [map[i], map[i+1]), so they need to be constructed from the first element of the next range, and map to curr range = next-1.
          int p;
          KeyToRank(int _comm_size) : p(_comm_size) {};

          /// return id of selected element based on lookup table.  ranges are [map[i], map[i+1])
          inline int operator()(Key const & x) const {
            auto pos = ::std::upper_bound(map.begin(), map.end(), x, Base::less);  // if equal, goes to next range.
            return (pos == map.end()) ? (p-1) : pos->second;
          }
          template<typename V>
          inline int operator()(::std::pair<Key, V> const & x) const {
            return this->operator()(x.first);
          }
          template<typename V>
          inline int operator()(::std::pair<const Key, V> const & x) const {
            return this->operator()(x.first);
          }
      } key_to_rank;


      /**
       * @brief count elements with the specified keys in the distributed sorted_multimap.
       * @note  input cannot have duplicate elements.
       *
       * @param first
       * @param last
       */
      template <bool unique = false>
      struct Intersect {


          // get the overlapping range in container.  container must be sorted.
          template<typename QueryIter, typename DBIter>
          static ::std::pair<DBIter, DBIter> intersect(DBIter db_begin, DBIter db_end,
                                                       QueryIter query_begin, QueryIter query_end,
                                                       bool sorted_query = false) {
            if (db_begin == db_end) return ::std::make_pair(db_begin, db_end);  // no target input
            if (query_begin == query_end) return ::std::make_pair(db_end, db_end);  // no src input

            auto dist_query = ::std::distance(query_begin, query_end);

            // make sure both are sorted.
            if (!sorted_query) Base::sort_ascending(query_begin, query_end);
            //            if (!sorted_dest) Base::sort_ascending(db_begin, db_end);  can't sort the container because its begin()/end() keeps returning const iterators if "this" is const.

            // find the bounds.
            auto range_begin = ::std::lower_bound(db_begin, db_end, *query_begin, Base::less);
            auto query_last = query_begin; ::std::advance(query_last, dist_query - 1);  // last real entry.
            auto range_end = ::std::upper_bound(range_begin, db_end, *query_last, Base::less);

            return ::std::make_pair(range_begin, range_end);
          }

          // assumes that container is sorted. and exact overlap region is provided.  do not filter output here since it's an output iterator.
          template <class DBIter, class QueryIter, class OutputIter, class Operator, class Predicate = Identity>
          static size_t process(DBIter range_begin, DBIter range_end,
                                QueryIter query_begin, QueryIter query_end,
                                OutputIter &output, Operator const & op,
                                bool sorted_query = false, Predicate const &pred = Predicate()) {

              // no matches in container.
              if (range_begin == range_end) return 0;
              if (query_begin == query_end) return 0;  // no input

              //auto output_start = output;

              auto dist_range = ::std::distance(range_begin, range_end);
              auto dist_query = ::std::distance(query_begin, query_end);

              //if (!sorted_target) Base::sort_ascending(range_begin, range_end);  range_begin and range_end often are const iterators.
              if (!sorted_query) Base::sort_ascending(query_begin, query_end);

              auto el_end = range_begin;
              size_t count = 0;
              bool linear = (static_cast<double>(dist_query) * ::std::log2(dist_range)) > dist_range;
              typename ::std::iterator_traits<QueryIter>::value_type v;

              if (linear) {  // based on number of input and search source, choose a method to search.

                // iterate through the input and search in output -
                for (auto it = query_begin; it != query_end;) {
                  v = *it;
                  if (!::std::is_same<Predicate, Identity>::value)
                	  count += op.template operator()<true>(range_begin, el_end, range_end, v, output, pred);
                  else
                	  count += op.template operator()<true>(range_begin, el_end, range_end, v, output);

                  // compiler optimize out the conditional.
                  if (unique) it = upper_bound<true>(it, query_end, v);
                  else ++it;
                }
              } else {
                // use logarithmic search

                // iterate through the input and search in output -
                for (auto it = query_begin; it != query_end;) {
                  v = *it;
                  if (!::std::is_same<Predicate, Identity>::value)
                	  count += op.template operator()<false>(range_begin, el_end, range_end, v, output, pred);
                  else
                	  count += op.template operator()<false>(range_begin, el_end, range_end, v, output);

                  // compiler optmizes out the conditional
                  if (unique) it = upper_bound<true>(it, query_end, v);
                  else ++it;
                }
              }

              return count;
          }

      };



    public:
      using local_container_type = Container<::std::pair<Key, T>, Alloc>;

      // std::sorted_multimap public members.
      using key_type              = Key;
      using mapped_type           = T;
      using value_type            = ::std::pair<Key, T>;
      using iterator              = typename local_container_type::iterator;
      using const_iterator        = typename local_container_type::const_iterator;
      using size_type             = typename local_container_type::size_type;
      using difference_type       = typename local_container_type::difference_type;

    protected:
      local_container_type c;
      bool sorted;
      bool balanced;
      bool globally_sorted;


      /// reserve space.  n is the local container size.  this allows different processes to individually adjust its own size.
      void local_reserve( size_t n) {
        // vector's reserve will only do something if n > capacity.
        c.reserve(n);
      }


      template <bool linear, class DBIter, class Query = typename ::std::iterator_traits<DBIter>::value_type>
      static inline DBIter lower_bound(DBIter b, DBIter e, Query& v) {
          // compiler choose one.
          if (linear) while ((b != e) && Base::less(*b, v)) ++b;
          else b = ::std::lower_bound(b, e, v, Base::less);
          return b;
      }


      template <bool linear, class DBIter, class Query = typename ::std::iterator_traits<DBIter>::value_type>
      static inline DBIter upper_bound(DBIter b, DBIter e,  Query& v) {
          // compiler choose one.
          if (linear) while ((b != e) && !Base::less(v, *b)) ++b;
          else b = ::std::upper_bound(b, e, v, Base::less);
          return b;
      }


      /// rehash the local container.  n is the local container size.  this allows different processes to individually adjust its own size.
      void local_sort() {
        if (sorted) return;
        Base::sort_ascending(c.begin(), c.end());
        sorted = true;
      }

      void assert_sorted_locally() const { if (!sorted) throw ::std::logic_error("local_sort needed to be called to sort local vector"); }
      bool is_sorted_locally() const { return sorted; }

      struct LocalCount {
          // unfiltered.
          template<bool linear, class DBIter, typename Query, class OutputIter>
          size_t operator()(DBIter &range_begin, DBIter &el_end, DBIter const &range_end, Query const &v, OutputIter &output) const {
              // TODO: LINEAR SEARCH, O(n).  vs BINARY SEARCH, O(mlogn)
              range_begin = lower_bound<linear>(el_end, range_end, v);  // range_begin at equal or greater than v.
              el_end = upper_bound<linear>(range_begin, range_end, v);  // el_end at greater than v.
              // difference between the 2 iterators is the part that's equal.

              // add the output entry.
              size_t count = ::std::distance(range_begin, el_end);

              *output = ::std::move(::std::make_pair(v, count));
              ++output;
              return 1;
          }
          // filtered element-wise.
          template<bool linear, class DBIter, typename Query, class OutputIter, class Predicate = Identity>
          size_t operator()(DBIter &range_begin, DBIter &el_end, DBIter const &range_end, Query const &v, OutputIter &output,
                            Predicate const& pred) const {
              // TODO: LINEAR SEARCH, O(n).  vs BINARY SEARCH, O(mlogn)
              range_begin = lower_bound<linear>(el_end, range_end, v);  // range_begin at equal or greater than v.
              el_end = upper_bound<linear>(range_begin, range_end, v);  // el_end at greater than v.
              // difference between the 2 iterators is the part that's equal.

              // add the output entry.
              size_t count = 0;
              if (pred(range_begin, el_end))  // operator to decide if range matches.
            	  count = ::std::count_if(range_begin, el_end, pred);  // operator for each element in range.

              *output = ::std::move(::std::make_pair(v, count));
              ++output;
              return 1;
          }
          // no filter by range AND elemenet for now.
      } count_element;


      struct LocalErase {
          /// Return how much was KEPT.
          template<bool linear, class DBIter, typename Query>
          size_t operator()(DBIter &range_begin, DBIter &el_end, DBIter const &range_end, Query const &v, DBIter &output) {
              // find start of segment to delete == end of prev segment to keep
              range_begin = lower_bound<linear>(el_end, range_end, v);

              // if the keep range is larger than 0, then move data and update insert pos.
              auto dist = ::std::distance(el_end, range_begin);
              if (dist > 0) {
                // condense a portion.
                ::std::move(el_end, range_begin, output);
                ::std::advance(output, dist);
              }
              // find of end of the segment to delete == start of next segment to keep
              el_end = upper_bound<linear>(range_begin, range_end, v);
              return dist;
          }
          /// Return how much was KEPT.
          template<bool linear, class DBIter, typename Query, class Predicate = Identity>
          size_t operator()(DBIter &range_begin, DBIter &el_end, DBIter const &range_end, Query const &v, DBIter &output,
                            Predicate const & pred) {
              DBIter output_orig = output;

              // find start of segment to delete == end of prev segment to keep
              range_begin = lower_bound<linear>(el_end, range_end, v);

              if (::std::distance(el_end, range_begin) > 0) {
                // condense a portion.
                output = ::std::move(el_end, range_begin, output);
              }
              // find of end of the segment to delete == start of next segment to keep
              el_end = upper_bound<linear>(range_begin, range_end, v);


              if (!pred(range_begin, range_end)) {  // if NOTHING in range matches, then no erase.
            	  output = ::std::move(range_begin, el_end, output);
            	  return 0;
              } else {  // some elements match, need to erase those
            	  size_t count = 0;

                //output = ::std::copy_if(range_begin, el_end, output, ::std::unary_negate<Predicate>(pred));
            	  for (auto it = range_begin; it < el_end; ++it) {
            	    if (!pred(*it)) {
            	      *output = *it;
            	      ++output;
            	      ++count;
            	    }
            	  }
            	  return count;
              }

              // all erased.
              return ::std::distance(range_begin, el_end);
          }
          // no filter by range AND elemenet for now.
      } erase_element;

      /**
       * @brief erase elements with the specified keys in the distributed sorted_multimap.  return how much was erased.
       * @note  this method is here because need to move the end segment.
       * @param first
       * @param last
       */
      template <typename Predicate = Identity>
      size_t local_erase(::std::vector<Key>& keys, bool sorted_input = false, Predicate const & pred = Predicate() ) {
        if (keys.size() == 0) return 0;
        if (c.size() == 0) return 0;

        this->balanced = false;

        this->local_sort();

        //== now call local remove.
        // first get the range of intersection
        auto overlap = Intersect<true>::intersect(this->c.begin(), this->c.end(), keys.begin(), keys.end(), sorted_input);

        if (::std::distance(overlap.first, overlap.second) == 0) return 0;
        size_t before = c.size();

        // do the work.  unique = true
        size_t kept = Intersect<true>::process(overlap.first, overlap.second, keys.begin(), keys.end(), overlap.first, erase_element, true, pred);

        // move the last part.
        auto end = overlap.first;
        ::std::advance(end, kept);
        if (end == overlap.second) return 0;  // nothing was removed.

        end = ::std::move(overlap.second, this->c.end(), end);
        this->c.erase(end, this->c.end());

        return before - c.size();
      }



      /// clears the sorted_map
      virtual void local_clear() noexcept {
        c.clear();
        this->sorted = true; this->balanced = true; this->globally_sorted = true;
      }

      ///  keep the unique keys in the input.   output is sorted.  equal operator forces comparison to Key
      template <typename V>
      void retain_unique(::std::vector< V >& input, bool sorted_input = false) const {
        if (input.size() == 0) return;
        if (!sorted_input) Base::sort_ascending(input.begin(), input.end());
        auto end = ::std::unique(input.begin(), input.end(), this->equal);
        input.erase(end, input.end());
      }


      /**
       * @brief find elements with the specified keys in the distributed sorted_multimap.
       * @param first
       * @param last
       */
      template <class LocalFind, class Predicate = Identity >
      ::std::vector<::std::pair<Key, T> > find_a2a(LocalFind const & local_find, ::std::vector<Key>& keys, bool sorted_input = false,
    		  Predicate const& pred = Predicate() ) const {
          BL_BENCH_INIT(find);

//          for (int j = 0; j < keys.size(); ++j) {
//            printf("fa2a rank %d originally has key %s\n", this->comm.rank(), keys[j].toAlphabetString().c_str());
//          }

          BL_BENCH_START(find);
          ::std::vector<::std::pair<Key, T> > results;
          // even if count is 0, still need to participate in mpi calls.  if (keys.size() == 0) return results;
          ::fsc::back_emplace_iterator<::std::vector<::std::pair<Key, T> > > emplace_iter(results);
          BL_BENCH_END(find, "begin", keys.size());

          this->assert_sorted_locally();

          // keep unique keys
          BL_BENCH_START(find);
          this->retain_unique(keys, sorted_input);
          BL_BENCH_END(find, "uniq1", keys.size());

          if (this->comm_size > 1) {

            BL_BENCH_START(find);
            // first remove duplicates.  sort, then get unique, finally remove the rest.  may not be needed
            ::std::vector<size_t> send_counts = fsc::get_bucket_sizes(keys, this->key_to_rank.map, Base::less);  // keys are sorted
            BL_BENCH_END(find, "bucket", keys.size());

            BL_BENCH_COLLECTIVE_START(find, "a2a1", this->comm);
            std::vector<size_t> recv_counts = mxx::all2all(send_counts, this->comm);
            keys = mxx::all2allv(keys, send_counts, this->comm);
            BL_BENCH_END(find, "a2a1", keys.size());

//            for (int j = 0; j < keys.size(); ++j) {
//              printf("fa2a rank %d after a2a has key %s\n", this->comm.rank(), keys[j].toAlphabetString().c_str());
//            }

            // local find. memory utilization a potential problem.
            // do for each src proc one at a time.

            BL_BENCH_START(find);
            results.reserve(keys.size() );
            BL_BENCH_END(find, "reserve", keys.size() );

            BL_BENCH_START(find);
            auto start = keys.begin();
            auto end = start;
            for (int i = 0; i < this->comm_size; ++i) {
              ::std::advance(end, recv_counts[i]);

              // work on query from process i.
              auto overlap = Intersect<false>::intersect(this->c.begin(), this->c.end(), start, end, true);

              // within start-end, values are unique, so don't need to set unique to true.
              send_counts[i] = Intersect<false>::process(overlap.first, overlap.second, start, end, emplace_iter, local_find, true, pred);

//              if (this->comm.rank() == 0) BL_DEBUGF("R %d added %d results for %d queries for process %d\n", this->comm.rank(), send_counts[i], recv_counts[i], i);

              start = end;
            }
            BL_BENCH_END(find, "local_find", results.size());

//            for (int j = 0; j < results.size(); ++j) {
//              printf("rank %d found %s\n", this->comm.rank(), results[j].first.toAlphabetString().c_str());
//            }

            // send back using the constructed recv count
            BL_BENCH_COLLECTIVE_START(find, "a2a2", this->comm);
            results = mxx::all2allv(results, send_counts, this->comm);
            BL_BENCH_END(find, "a2a2", results.size());

//            for (int j = 0; j < results.size(); ++j) {
//              printf("rank %d moved results %s\n", this->comm.rank(), results[j].first.toAlphabetString().c_str());
//            }

          } else {

            BL_BENCH_START(find);
            results.reserve(keys.size());  // 1 result per key.
            BL_BENCH_END(find, "reserve", keys.size() );


            BL_BENCH_START(find);
            auto overlap = Intersect<false>::intersect(this->c.begin(), this->c.end(), keys.begin(), keys.end(), true);

            // within start-end, values are unique, so don't need to set unique to true.
            Intersect<false>::process(overlap.first, overlap.second, keys.begin(), keys.end(), emplace_iter, local_find, true, pred);

            BL_BENCH_END(find, "local_find", results.size());

          }

          BL_BENCH_REPORT_MPI(find, this->comm.rank(), this->comm);

          return results;
      }


      /**
       * @brief find elements with the specified keys in the distributed sorted_multimap.
       * @param first
       * @param last
       */
      template <class LocalFind, class Predicate = Identity >
      ::std::vector<::std::pair<Key, T> > find(LocalFind const & local_find, ::std::vector<Key>& keys, bool sorted_input = false,
          Predicate const& pred = Predicate() ) const {
          BL_BENCH_INIT(find);
//
//          for (int j = 0; j < keys.size(); ++j) {
//            printf("rank %d originally has key %s\n", this->comm.rank(), keys[j].toAlphabetString().c_str());
//          }

          BL_BENCH_START(find);
          ::std::vector<::std::pair<Key, T> > results;
          // even if count is 0, still need to participate in mpi calls.  if (keys.size() == 0) return results;
          ::fsc::back_emplace_iterator<::std::vector<::std::pair<Key, T> > > emplace_iter(results);

          ::std::vector<::std::pair<Key, T> > local_results;
          ::fsc::back_emplace_iterator<::std::vector<::std::pair<Key, T> > > local_emplace_iter(local_results);
          BL_BENCH_END(find, "begin", keys.size());

          this->assert_sorted_locally();

          // keep unique keys
          BL_BENCH_START(find);
          this->retain_unique(keys, sorted_input);
          BL_BENCH_END(find, "uniq1", keys.size());

          if (this->comm_size > 1) {

            BL_BENCH_START(find);
            // first remove duplicates.  sort, then get unique, finally remove the rest.  may not be needed
            ::std::vector<size_t> send_counts = fsc::get_bucket_sizes(keys, this->key_to_rank.map, Base::less);
            BL_BENCH_END(find, "bucket", keys.size());

            BL_BENCH_COLLECTIVE_START(find, "a2a1", this->comm);
            std::vector<size_t> recv_counts = mxx::all2all(send_counts, this->comm);
            keys = mxx::all2allv(keys, send_counts, this->comm);
            BL_BENCH_END(find, "a2a1", keys.size());


            // local count to determine amount of memory to allocate at destination.
            BL_BENCH_START(find);

            ::std::vector<::std::pair<Key, size_t> > count_results;
            size_t max_key_count = *(::std::max_element(recv_counts.begin(), recv_counts.end()));
            count_results.reserve(max_key_count);
            ::fsc::back_emplace_iterator<::std::vector<::std::pair<Key, size_t> > > count_emplace_iter(count_results);

            auto start = keys.begin();
            auto end = start;
            size_t total = 0;
            for (int i = 0; i < this->comm_size; ++i) {
              ::std::advance(end, recv_counts[i]);

              // count results for process i
              count_results.clear();
              auto overlap = Intersect<false>::intersect(this->c.begin(), this->c.end(), start, end, true);
              Intersect<false>::process(overlap.first, overlap.second, start, end, count_emplace_iter, count_element, true, pred);
              send_counts[i] = ::std::accumulate(count_results.begin(), count_results.end(), static_cast<size_t>(0),
                                                 [](size_t v, ::std::pair<Key, size_t> const & x) {
                             return v + x.second;
                           });
//              for (auto it = count_results.begin(), max = count_results.end(); it != max; ++it) {
//                send_counts[i] += it->second;
//              }
              total += send_counts[i];

              start = end;
              //printf("Rank %d local count for src rank %d:  recv %d send %d\n", this->comm.rank(), i, recv_counts[i], send_counts[i]);
            }
            ::std::vector<::std::pair<Key, size_t> >().swap(count_results);
            BL_BENCH_END(find, "local_count", total);


            BL_BENCH_COLLECTIVE_START(find, "a2a_count", this->comm);
            std::vector<size_t> resp_counts = mxx::all2all(send_counts, this->comm);  // compute counts of response to receive
            BL_BENCH_END(find, "a2a_count", keys.size());



            BL_BENCH_START(find);
            auto resp_displs = mxx::impl::get_displacements(resp_counts);  // compute response displacements.

            auto resp_total = resp_displs[this->comm_size - 1] + resp_counts[this->comm_size - 1];
            auto max_send_count = *(::std::max_element(send_counts.begin(), send_counts.end()));
            results.resize(resp_total);   // allocate, not just reserve
            local_results.reserve(max_send_count);
            BL_BENCH_END(find, "reserve", resp_total);

            BL_BENCH_START(find);
            auto recv_displs = mxx::impl::get_displacements(recv_counts);  // compute response displacements.
            int recv_from, send_to;
            size_t found;
            total = 0;
            std::vector<MPI_Request> reqs(2 * this->comm_size);

            mxx::datatype dt = mxx::get_datatype<::std::pair<Key, T>>();
            for (int i = 0; i < this->comm_size; ++i) {
              recv_from = (this->comm.rank() + (this->comm_size - i)) % this->comm_size; // rank to recv data from

              // set up receive.
              MPI_Irecv(&results[resp_displs[recv_from]], resp_counts[recv_from], dt.type(),
                        recv_from, i, this->comm, &reqs[2 * i]);


              send_to = (this->comm.rank() + i) % this->comm_size;    // rank to send data to

              //== get data for the dest rank
              start = keys.begin();                                   // keys for the query for the dest rank
              ::std::advance(start, recv_displs[send_to]);
              end = start;
              ::std::advance(end, recv_counts[send_to]);

              local_results.clear();
              // work on query from process i.
              auto overlap = Intersect<false>::intersect(this->c.begin(), this->c.end(), start, end, true);
              found = Intersect<false>::process(overlap.first, overlap.second, start, end, local_emplace_iter, local_find, true, pred);
             // if (this->comm.rank() == 0) BL_DEBUGF("R %d added %d results for %d queries for process %d\n", this->comm.rank(), send_counts[i], recv_counts[i], i);
              total += found;
              //== now send the results immediately - minimizing data usage so we need to wait for both send and recv to complete right now.

//              for (int j = 0; j < local_results.size(); ++j) {
//                printf("rank %d -> %d sent %s\n", this->comm.rank(), send_to, local_results[j].first.toAlphabetString().c_str());
//              }

              MPI_Isend(&(local_results[0]), found, dt.type(), send_to,
                        i, this->comm, &reqs[2 * i + 1]);

              // wait for both requests to complete.
              MPI_Waitall(2, &reqs[2 * i], MPI_STATUSES_IGNORE);

//              for (int j = 0; j < resp_counts[recv_from]; ++j) {
//                printf("rank %d -> %d recv %s\n", recv_from, this->comm.rank(), results[resp_displs[recv_from] + j].first.toAlphabetString().c_str());
//              }

              // within start-end, values are unique, so don't need to set unique to true.

            }

            BL_BENCH_END(find, "local_find", results.size());
//
//            // send back using the constructed recv count
//            BL_BENCH_COLLECTIVE_START(find, "a2a2", this->comm);
//            mxx2::all2all(results, send_counts, this->comm);
//            BL_BENCH_END(find, "a2a2", results.size());
//

          } else {

            BL_BENCH_START(find);

            ::std::vector<::std::pair<Key, size_t> > count_results;
            count_results.reserve(keys.size());
            ::fsc::back_emplace_iterator<::std::vector<::std::pair<Key, size_t> > > count_emplace_iter(count_results);

            // count now.
            auto overlap = Intersect<false>::intersect(this->c.begin(), this->c.end(), keys.begin(), keys.end(), sorted_input);
            Intersect<false>::process(overlap.first, overlap.second, keys.begin(), keys.end(), count_emplace_iter, count_element, true, pred);
            size_t count = ::std::accumulate(count_results.begin(), count_results.end(), static_cast<size_t>(0),
                                          [](size_t v, ::std::pair<Key, size_t> const & x) {
                      return v + x.second;
                    });
//            for (auto it = count_results.begin(), max = count_results.end(); it != max; ++it) {
//              count += it->second;
//            }
            BL_BENCH_END(find, "local_count", count);

            BL_BENCH_START(find);
            results.reserve(count);  // 1 result per key.
            BL_BENCH_END(find, "reserve", count);


            BL_BENCH_START(find);

            // within start-end, values are unique, so don't need to set unique to true.
            Intersect<false>::process(overlap.first, overlap.second, keys.begin(), keys.end(), emplace_iter, local_find, true, pred);

            BL_BENCH_END(find, "local_find", results.size());

          }

          BL_BENCH_REPORT_MPI(find, this->comm.rank(), this->comm);

          return results;
      }


      template <class LocalFind, class Predicate = Identity>
      ::std::vector<::std::pair<Key, T> > find(LocalFind const & local_find,
    		  Predicate const & pred = Predicate()) const {
          ::std::vector<::std::pair<Key, T> > results;
          ::fsc::back_emplace_iterator<::std::vector<::std::pair<Key, T> > > emplace_iter(results);

          auto keys = this->keys();

          ::std::vector<::std::pair<Key, size_t> > count_results;
          count_results.reserve(keys.size());
          ::fsc::back_emplace_iterator<::std::vector<::std::pair<Key, size_t> > > count_emplace_iter(count_results);

          // count now.
          Intersect<false>::process(this->c.begin(), this->c.end(), keys.begin(), keys.end(), count_emplace_iter, count_element, true, pred);
          size_t count = ::std::accumulate(count_results.begin(), count_results.end(), static_cast<size_t>(0),
                                           [](size_t v, ::std::pair<Key, size_t> const & x) {
                       return v + x.second;
                     });
//          for (auto it = count_results.begin(), max = count_results.end(); it != max; ++it) {
//            count += it->second;
//          }

          results.reserve(count);  // 1 result per key.

          // within start-end, values are unique, so don't need to set unique to true.
          Intersect<false>::process(this->c.begin(), this->c.end(), keys.begin(), keys.end(), emplace_iter, local_find, true, pred);

          if (this->comm_size > 1) MPI_Barrier(this->comm);
          return results;
      }


      virtual void rehash() = 0;



      sorted_map_base(const mxx::comm& _comm) : Base(_comm),
          key_to_rank(_comm.size()), sorted(false), balanced(false), globally_sorted(false) {}

    public:

      virtual ~sorted_map_base() {};

      /// returns the local storage.  please use sparingly.
      local_container_type& get_local_container() { return c; }

      const_iterator cbegin() const
      {
        return c.cbegin();
      }

      const_iterator cend() const {
        return c.cend();
      }

      /// update the multiplicity.  only multimap needs to do this.
      virtual size_t update_multiplicity() {
        this->rehash();
        return this->key_multiplicity;
      }

      /// convert the map to a vector.
      virtual std::vector<std::pair<Key, T> > to_vector() const {
        std::vector<std::pair<Key, T> > result(c.begin(), c.end());
        return result;
      }

      /// convert the map to a vector
      virtual void to_vector(std::vector<std::pair<Key, T> > & result) const {
        result.clear();
        if (c.empty()) return;

        result.assign(c.begin(), c.end());
      }

      /// extract the keys of a map.
      virtual std::vector<Key> keys() const {
        std::vector<Key> result;
        this->keys(result);
        return result;
      }

      /// extract the unique keys of a map.
      virtual void keys(std::vector<Key> & result) const {
        result.clear();
        if (c.empty()) return;

        // copy the keys
        auto end = c.end();
        for (auto it = c.begin(); it != end; ++it) {
          result.emplace_back(it->first);
        }

        // and then find unique.
        retain_unique(result, this->sorted);
      }

      // note that for each method, there is a local version of the operartion.
      // this is for use by the asynchronous version of communicator as callback for any messages received.
      /// check if empty.
      virtual bool local_empty() const noexcept {
        return c.empty();
      }

      /// get size of local container
      virtual size_t local_size() const noexcept {
        return c.size();
      }


      /// reserve space.  n is the local container size.  this allows different processes to individually adjust its own size.
      void reserve( size_t n) {
        // direct reserve + barrier
        this->local_reserve(n);
        if (this->comm_size > 1) MPI_Barrier(this->comm);
      }





      /**
       * @brief count elements with the specified keys in the distributed sorted_multimap.
       * @param first
       * @param last
       */
      template <typename Predicate = Identity>
      ::std::vector<::std::pair<Key, size_type> > count(::std::vector<Key>& keys, bool sorted_input = false,
    		  Predicate const & pred = Predicate()) const {


        BL_BENCH_INIT(count);

        BL_BENCH_START(count);

                  // keep unique keys
        ::std::vector<::std::pair<Key, size_type> > results;
        // even if count is 0, still need to participate in mpi calls.  if (keys.size() == 0) return results;
        ::fsc::back_emplace_iterator<::std::vector<::std::pair<Key, size_type> > > emplace_iter(results);

        this->assert_sorted_locally();
        BL_BENCH_END(count, "begin", keys.size());

        if (this->comm_size > 1) {
          BL_BENCH_START(count);
          // keep unique keys
          retain_unique(keys, sorted_input);
          BL_BENCH_END(count, "uniq1", keys.size());

          BL_BENCH_START(count);
          // first remove duplicates.  sort, then get unique, finally remove the rest.  may not be needed
          ::std::vector<size_t> send_counts = fsc::get_bucket_sizes(keys, this->key_to_rank.map, Base::less);
          BL_BENCH_END(count, "bucket", keys.size());

          BL_BENCH_COLLECTIVE_START(count, "a2a1", this->comm);
          std::vector<size_t> recv_counts = mxx::all2all(send_counts, this->comm);
          keys = mxx::all2allv(keys, send_counts, this->comm);
          BL_BENCH_END(count, "a2a1", keys.size());


          BL_BENCH_START(count);
          // local find. memory utilization a potential problem.
          // do for each src proc one at a time.

          results.reserve(keys.size());                   // TODO:  should estimate coverage.

          BL_BENCH_END(count, "reserve", keys.size());

          BL_BENCH_START(count);

          auto start = keys.begin();
          auto end = start;
          for (int i = 0; i < this->comm_size; ++i) {
            ::std::advance(end, recv_counts[i]);

            // work on query from process i.
            auto overlap = Intersect<false>::intersect(this->c.begin(), this->c.end(), start, end, true);

            // within start-end, values are unique, so don't need to set unique to true.
            Intersect<false>::process(overlap.first, overlap.second, start, end, emplace_iter, count_element, true, pred);

            if (this->comm.rank() == 0) BL_DEBUGF("R %d added %lu results for %lu queries for process %d\n", this->comm.rank(), send_counts[i], recv_counts[i], i);


            start = end;
          }
          BL_BENCH_END(count, "local_count", results.size());

          BL_BENCH_COLLECTIVE_START(count, "a2a2", this->comm);

          // send back using the constructed recv count
          results = mxx::all2allv(results, recv_counts, this->comm);
          BL_BENCH_END(count, "a2a2", results.size());


        } else {

          BL_BENCH_START(count);

          results.reserve(keys.size());
          BL_BENCH_END(count, "reserve", keys.size());

          BL_BENCH_START(count);
          // work on query from process i.
          auto overlap = Intersect<true>::intersect(this->c.begin(), this->c.end(), keys.begin(), keys.end(), sorted_input);

          // within key, values may not be unique,
          Intersect<true>::process(overlap.first, overlap.second, keys.begin(), keys.end(), emplace_iter, count_element, true, pred);
          BL_BENCH_END(count, "local_count", results.size());

        }
        BL_BENCH_REPORT_MPI(count, this->comm.rank(), this->comm);

        return results;

      }


      template <typename Predicate = Identity>
      ::std::vector<::std::pair<Key, size_type> > count(Predicate const & pred = Predicate()) const {
        ::std::vector<::std::pair<Key, size_type> > results;
        ::fsc::back_emplace_iterator<::std::vector<::std::pair<Key, size_type> > > emplace_iter(results);

        auto keys = this->keys();
        results.reserve(keys.size());

        // keys already unique
        Intersect<false>::process(c.begin(), c.end(), keys.begin(), keys.end(), emplace_iter, count_element, true, pred);

        if (this->comm_size > 1) MPI_Barrier(this->comm);
        return results;
      }


//      /**
//       * @brief insert new elements in the distributed sorted_multimap.  example use: stop inserting if more than x entries.
//       * @param src_begin
//       * @param src_end
//       */
//      template <class InputIter, class Predicate = Identity>
//      size_t insert(InputIter src_begin, InputIter src_end, bool sorted_input = false, Predicate const &pred = Predicate()) {
//          if (src_begin == src_end) return 0;
//          BL_BENCH_INIT(insert);
//
//          BL_BENCH_START(insert);
//
//          this->sorted = false; this->balanced = false; this->globally_sorted = false;
//          ::fsc::back_emplace_iterator<local_container_type> emplace_iter(c);
//
//          if (::std::is_same<Predicate, Identity>::value) ::std::copy(src_begin, src_end, emplace_iter);
//          else ::std::copy_if(src_begin, src_end, emplace_iter, pred);
//
//          size_t count = ::std::distance(src_begin, src_end);
//          BL_BENCH_END(insert, "insert", count);
//
//
//          BL_BENCH_REPORT_MPI(insert, this->comm.rank(), this->comm);
//
//          return count;
//      }

      /**
       * @brief insert new elements in the distributed sorted_multimap.  example use: stop inserting if more than x entries.
       * @param first
       * @param last
       */
      template <class Predicate = Identity>
      size_t insert(::std::vector<::std::pair<Key, T> > &input, bool sorted_input = false, Predicate const &pred = Predicate()) {
          if (input.size() == 0) return 0;
          BL_BENCH_INIT(insert);
          this->sorted = false; this->balanced = false; this->globally_sorted = false;

          size_t before = c.size();
          BL_BENCH_START(insert);

          ::fsc::back_emplace_iterator<local_container_type> emplace_iter(c);
          if (::std::is_same<Predicate, Identity>::value) {
            if (c.size() == 0)   // container is empty, so swap it in.
              c.swap(input);
            else {
            	this->local_reserve(before + input.size());
                ::std::move(input.begin(), input.end(), emplace_iter);    // else move it in.
            }
          }
          else {
        	  this->local_reserve(before + input.size());
        	  ::std::copy_if(::std::make_move_iterator(input.begin()),
                      ::std::make_move_iterator(input.end()), emplace_iter, pred);  // predicate needed.  move it though.
          }

          size_t count = c.size() - before;
          BL_BENCH_END(insert, "insert", count);


          BL_BENCH_REPORT_MPI(insert, this->comm.rank(), this->comm);

          return count;
      }


      /**
       * @brief erase elements with the specified keys in the distributed sorted_multimap.  return how much was erased.
       * @param first
       * @param last
       */
      template <typename Predicate = Identity>
      size_t erase(::std::vector<Key>& keys, bool sorted_input = false, Predicate const & pred = Predicate() ) {
        // even if count is 0, still need to participate in mpi calls.  if (keys.size() == 0) return;

        bool si = sorted_input;
        if (this->comm_size > 1) {
          // remove duplicates
          retain_unique(keys, si);

          // first remove duplicates.  sort, then get unique, finally remove the rest.  may not be needed
          ::std::vector<size_t> send_counts = fsc::get_bucket_sizes(keys, this->key_to_rank.map, Base::less);

          keys = mxx::all2allv(keys, send_counts, this->comm);

          si = false;
        }

        return this->local_erase(keys, si, pred);

      }

      template <typename Predicate = Identity>
      size_t erase(Predicate const & pred = Predicate()) {
        size_t before = c.size();
        if (!::std::is_same<Predicate, Identity>::value) {

          local_sort();

          auto end = ::std::partition(c.begin(), c.end(), [pred](value_type const &x) { return !pred(x.first); });
          c.erase(end, c.end());

          this->balanced = false;
        } else {
          // identity predicate - all are erased.
          this->local_clear();
        }

        if (this->comm_size > 1) MPI_Barrier(this->comm);

        return before - c.size();
      }

      // update done via erase/insert.


  };



  /**
   * @brief  distributed unordered map following std unordered map's interface.
   * @details   This class is modeled after the std::unordered_map.
   *         it has as much of the same methods of std::unordered_map as possible.  however, all methods consider the fact
   *         that the data are in distributed memory space, so to access the data, "communication" is needed.
   *
   *         Note that "communication" is a weak concept here meaning that we are accessing a different local container.
   *         as such, communicator may be defined for MPI, UPC, OpenMP, etc.
   *
   *         This allows the possibility of using distributed unordered map as local storage for coarser grain distributed container.
   *
   *         Note that communicator requires a mapping strategy between a key and the target processor/thread/partition.  The mapping
   *         may be done using a hash, similar to the local distributed unordered map, or it may be done via sorting/lookup or other mapping
   *         mechanisms.  The choice may be constrained by the communication approach, e.g. global sorting  does not work well with
   *         incremental async communication
   *
   * @tparam Key
   * @tparam T
   * @tparam Comm   default to mpi_collective_communicator       communicator for global communication. may hash or sort.
   * @tparam KeyTransform   transform function for the key.  can supply identity.  requires a single template argument (Key).  useful for mapping kmolecule to kmer.
   * @tparam Hash   hash function for local and distribution.  requires a template arugment (Key), and a bool (prefix, chooses the MSBs of hash instead of LSBs)
   * @tparam Equal   default to ::std::equal_to<Key>   equal function for the local storage.
   * @tparam Alloc  default to ::std::allocator< ::std::pair<Key, T> >    allocator for local storage.
   */
  template<typename Key, typename T,
  class Comm,
  template <typename> class KeyTransform,
  class Less = ::std::less<Key>,
  class Equal = ::std::equal_to<Key>,
  class Alloc = ::std::allocator< ::std::pair<Key, T> >
  >
  class sorted_map : public sorted_map_base<Key, T, ::std::vector, Comm, KeyTransform, Less, Equal, Alloc> {
    protected:
      using Base = sorted_map_base<Key, T, ::std::vector, Comm, KeyTransform, Less, Equal, Alloc>;


    public:
      using local_container_type = typename Base::local_container_type;

      // std::unordered_multimap public members.
      using key_type              = Key;
      using mapped_type           = T;
      using value_type            = ::std::pair<Key, T>;
      using iterator              = typename local_container_type::iterator;
      using const_iterator        = typename local_container_type::const_iterator;
      using size_type             = typename local_container_type::size_type;
      using difference_type       = typename local_container_type::difference_type;


    protected:

      // defined Communicator as a friend
      friend Comm;

      struct LocalFind {
          template<bool linear, class DBIter, typename Query, class OutputIter>
          size_t operator()(DBIter &range_begin, DBIter &el_end, DBIter const &range_end, Query const &v, OutputIter &output) const {
              // TODO: LINEAR SEARCH, O(n).  vs BINARY SEARCH, O(mlogn)

              // map, so only 1 entry.
              range_begin = Base::template lower_bound<linear>(el_end, range_end, v);
              el_end = range_begin;

              // add the output entry, if found.
              if (range_begin == range_end) return 0;
              if (Base::equal(*range_begin, v)) {
                *output = *range_begin;
                ++output;
                ++el_end;
                return 1;
              } else
            	  return 0;

          }
          template<bool linear, class DBIter, typename Query, class OutputIter, class Predicate = Identity>
          size_t operator()(DBIter &range_begin, DBIter &el_end, DBIter const &range_end,
        		  Query const &v, OutputIter &output, Predicate const &pred) const {
              // TODO: LINEAR SEARCH, O(n).  vs BINARY SEARCH, O(mlogn)

              // map, so only 1 entry.
              range_begin = Base::template lower_bound<linear>(el_end, range_end, v);
              el_end = range_begin;

              // add the output entry, if found.
              if (range_begin == range_end) return 0;
              if (Base::equal(*range_begin, v) && pred(v)) {
                *output = *range_begin;
                ++output;
                ++el_end;
                return 1;
              } else
            	  return 0;
          }
      } find_element;


      //      /**
      //       * @brief insert new elements in the distributed sorted_map.  example use: stop inserting if more than x entries.
      //       * @param first
      //       * @param last
      //       */
      //      template <class InputIterator, class Predicate>
      //      void local_insert_if(InputIterator first, InputIterator last, Predicate const &pred, bool sorted_input = false) {
      //          if (first == last) return;
      //
      //          // sort both
      //          this->local_rehash();
      //          if (!sorted_input) Base::Base::sort_ascending(first, last);
      //
      //          // walk through.  if exiting entry, remove from input.
      //          // now walk through the input to count stuff.
      //          auto t_start = ::std::lower_bound(this->c.begin(), this->c.end(), *first, Base::Base::less);
      //          auto last2 = first; ::std::advance(last2, ::std::distance(first, last) - 1);
      //          auto t_end = ::std::upper_bound(t_start, this->c.end(), *last2, Base::Base::less);
      //
      //          // go through the input and search the container.
      //          // iterate through the input and search in output,  moving items up as we go.
      //          auto it = first;
      //          auto end = first;
      //            for (it = first; it != last;) {
      //              auto v = *it;
      //
      //              if (pred(v)) {
      //                // find entry in container.
      //                t_start = this->template lower_bound<true>(t_start, t_end, v);
      //
      //                if ((t_start == t_end) || (!this->equal(v, *t_start))) {
      //                  // not matched.  copy the value, move pos up by 1.
      //                  if (end != it) *end = v;
      //                  ++end;
      //                } // else matched.  so need to skip it.  move it, but don't move pos
      //              }
      //
      //              it = this->upper_bound<true>(it, last, v);
      //
      //            }
      //          // insert at the end.
      //          if (first != end) {
      //            this->c.insert(this->c.end(), first, end);
      //            this->sorted = false;
      //          }
      //      }

      virtual void local_reduction(::std::vector<::std::pair<Key, T> > &input, bool sorted_input = false) {
        this->Base::retain_unique(input, sorted_input);
      }

    public:

      sorted_map(const mxx::comm& _comm) : Base(_comm) {}

      virtual ~sorted_map() {};

      template <class Predicate = Identity>
      ::std::vector<::std::pair<Key, T> > find(::std::vector<Key>& keys, bool sorted_input = false,
    		  Predicate const& pred = Predicate()) const {
          return Base::find(find_element, keys, sorted_input, pred);
      }

      template <class Predicate = Identity>
      ::std::vector<::std::pair<Key, T> > find(Predicate const& pred = Predicate()) const {
          return Base::find(find_element, pred);
      }


      // default is for map
      virtual void rehash() {
        BL_BENCH_INIT(rehash);

        //printf("c size before: %lu\n", this->c.size());
        BL_BENCH_START(rehash);
        BL_BENCH_END(rehash, "begin", this->c.size());

        if (this->comm_size > 1) {
          // first balance


          if (!this->balanced) {
            BL_BENCH_START(rehash);
            this->c = ::mxx::stable_distribute(this->c, this->comm);
            BL_BENCH_END(rehash, "block1", this->c.size());
          }


          // sort if needed
//          if (!this->globally_sorted) {
//          BL_BENCH_START(rehash);
//          // kway merge / sort
//          Base::Base::sort_ascending(this->c.begin(), this->c.end());
//          BL_BENCH_END(rehash, "sort1", this->c.size());
//
          // sort if needed
          BL_BENCH_START(rehash);
          if (!this->globally_sorted) ::mxx::sort(this->c.begin(), this->c.end(), Base::Base::less, this->comm);
          BL_BENCH_END(rehash, "mxxsort", this->c.size());


            BL_BENCH_START(rehash);
            // local unique
            this->local_reduction(this->c, true);
            BL_BENCH_END(rehash, "reduc1", this->c.size());

            // check to see if any value cross the boundary:
            // get the values from both sides of the boundaries, gather to rank 0
            // merge duplicate boundary values.
            // use a vector to indicate which to delete and which to keep.  (keep entry at start of a partition)
            // apply to remote.

            BL_BENCH_START(rehash);
            // allgather the left and right elements.  then each proc makes own decision.  relies on uniqueness within a proc.
            ::std::vector<value_type > boundary_values;
            ::std::vector<int > boundary_ids;
            if (this->c.size() > 1) {  // insert left if we have more than 1 entry.  (if 0, no insert.  if 1, will be inserted as right)
              boundary_values.emplace_back(this->c.front());
              boundary_ids.emplace_back(this->comm.rank());
            }
            if (this->c.size() > 0) {  // insert right if we have at least 1 entry
              boundary_values.emplace_back(this->c.back());
              boundary_ids.emplace_back(this->comm.rank());
            }
            // allgather the boundary entries and ranks.
            boundary_values = ::mxx::allgatherv(boundary_values, this->comm);
            boundary_ids = ::mxx::allgatherv(boundary_ids, this->comm);

            if (this->c.size() > 0) {
              // now reduce the range matching the left element and assign to last rank in that range.
              auto range = ::std::equal_range(boundary_values.begin(), boundary_values.end(), this->c.front(),
                                              Base::Base::less);
              size_t dist = ::std::distance(range.first, range.second);
              if (dist > 1) {  // if > 1, then value crossed boundary and need to reduce.

                size_t offset = ::std::distance(boundary_values.begin(), range.second);
                if (boundary_ids[offset - 1] == this->comm.rank()) { // this processor owns this value
                  // copy and reduce this range.
                  ::std::vector<value_type > reduced(range.first, range.second);
                  this->local_reduction(reduced, true);  //
                  this->c.front() = reduced[0];    // replace the left element
                } // else processor does not own the value so no update
              } // else there is only 1 entry, done. (guaranteed at least 1)

              // now remove the other entries from the local container
              range = ::std::equal_range(boundary_values.begin(), boundary_values.end(), this->c.back(),
                                                          Base::Base::less);
              dist = ::std::distance(range.first, range.second);
              if (dist > 1) {  // if > 1, the value crossed boundary and we have extra entries to remove.

                size_t offset = ::std::distance(boundary_values.begin(), range.second);
                if (boundary_ids[offset - 1] < this->comm.rank()) {  // not an owner, so need to remove this element
                  // now remove this element.  there is only 1, since it's unique.  also c is at least 1 in size before.
                  this->c.pop_back();

                }
              }  // back entry is unique.  keep.
            }
            BL_BENCH_END(rehash, "reduced boundaries", this->c.size());

            // then reblock.



//            // goal of the next block of code is to rebalance of the blocks, while keeping all of same elements in the same processor
//            // first thought is to resample using the sample sort logic.  however this is flawed and will introduce imbalance.
//            // for sorted input using the sample_arbit_decomp would result in imbalance that increases with decreasing p:
//            // first p-1 get (p-1) partitions, last p get p+1 partitions.  we can increase sampling rate,
//            // but that still does not remove the imbalance.  in addition, sampling is even within the blocks, which is not good for sorted blocks.
//
//            // alternative is to reblock first.  for there to be a run > 2 of same element, a processor would have had a single element.
//            // this means that the number of elements for the processor with the run post rebalance is not going to be very different,
//            // so we can just let the some procs idle (after their size 1 partitions are moved entirely.)
//
//            // rebalance block first
//            BL_BENCH_START(rehash);
//            this->c = ::mxx::stable_block_decompose(this->c, this->comm);
//            BL_BENCH_END(rehash, "block1", this->c.size());
//
// //            BL_BENCH_START(rehash);
// //            // sample
// //            mxx::datatype<::std::pair<Key, T> > dt;
// //            MPI_Datatype mpi_dt = dt.type();
// //            this->key_to_rank.map = ::mxx::impl::sample_arbit_decomp(this->c.begin(), this->c.end(), Base::Base::less, this->comm_size - 1, this->comm, mpi_dt);
// //            for (int i = 0; i < this->key_to_rank.map.size(); ++i) {
// //              // modify the splitters destinations
// //              //printf("R %d splitters %s -> %d\n", this->comm.rank(), this->key_to_rank.map[i].first.toAlphabetString().c_str(), this->key_to_rank.map[i].second);
// //              this->key_to_rank.map[i].second = i;
// //            }
// //            BL_BENCH_END(rehash, "splitter1", this->key_to_rank.map.size());
//
//            // get pivots using the last elements.
//            BL_BENCH_START(rehash);
//            this->key_to_rank.map.clear();
//            if ((this->comm.rank() > 0) && (this->c.size() > 0)) {  // splitters need to be the first entry of the next partition.
//              // only send for the first p-1 proc, and only if they have a kmer to split with.
//              this->key_to_rank.map.emplace_back(this->c.front().first, this->comm.rank() - 1);
//            }
//            this->key_to_rank.map = ::mxx::allgatherv(this->key_to_rank.map, this->comm);
//            BL_BENCH_END(rehash, "splitter1", this->c.size());
//
//            // rebucket by pivots, so no value crosses boundaries
//            BL_BENCH_START(rehash);
//            ::std::vector<size_t> send_counts = mxx2::bucketing<size_t>(this->c, this->key_to_rank.map, Base::Base::less);
//            BL_BENCH_END(rehash, "bucket", this->c.size());
//
//            BL_BENCH_COLLECTIVE_START(rehash, "a2a", this->comm);
//            mxx2::all2all(this->c, send_counts, this->comm);
//            BL_BENCH_END(rehash, "a2a", this->c.size());
//
//
////          BL_BENCH_START(rehash);
////          // kway merge / sort
////          Base::Base::sort_ascending(this->c.begin(), this->c.end());
////          BL_BENCH_END(rehash, "sort2", this->c.size());
////
//
//            // final reduction - nothing crosses boundaries now.
//            BL_BENCH_START(rehash);
//            // local unique
//            this->local_reduction(this->c, true);
//            BL_BENCH_END(rehash, "reduc2", this->c.size());

            // and final rebalance
            BL_BENCH_START(rehash);
            // rebalance
            this->c = ::mxx::stable_distribute(this->c, this->comm);
            BL_BENCH_END(rehash, "block2", this->c.size());

//          }
//
          BL_BENCH_START(rehash);
          // get new pivots
          // next compute the splitters.
          this->key_to_rank.map.clear();
          if ((this->comm.rank() > 0) && (this->c.size() > 0)) {
            // only send for the first p-1 proc, and only if they have a kmer to split with.
            this->key_to_rank.map.emplace_back(this->c.front().first, this->comm.rank() - 1);
          }
          this->key_to_rank.map = ::mxx::allgatherv(this->key_to_rank.map, this->comm);
          // note that key_to_rank.map needs to be unique.
          auto map_end = std::unique(this->key_to_rank.map.begin(), this->key_to_rank.map.end(), Base::Base::equal);
          this->key_to_rank.map.erase(map_end, this->key_to_rank.map.end());


//          for (int i = 0; i < this->key_to_rank.map.size(); ++i) {
//            printf("R %d key to rank %s -> %d\n", this->comm.rank(), this->key_to_rank.map[i].first.toAlphabetString().c_str(), this->key_to_rank.map[i].second);
//          }
          BL_BENCH_END(rehash, "splitter2", this->c.size());
          // no need to redistribute - each entry is unique so nothing is going to span processor boundaries.

        } else {
          BL_BENCH_START(rehash);
          // local unique
          this->local_reduction(this->c, this->sorted);

          BL_BENCH_END(rehash, "reduc", this->c.size());
        }
        this->sorted = true; this->balanced = true; this->globally_sorted = true;
        //printf("c size after: %lu\n", this->c.size());

        BL_BENCH_REPORT_MPI(rehash, this->comm.rank(), this->comm);

      }
  };



  /**
   * @brief  distributed unordered multimap following std unordered multimap's interface.
   * @details   This class is modeled after the std::unordered_multimap.
   *         it does not have all the methods of std::unordered_multimap.  Whatever methods that are present considers the fact
   *         that the data are in distributed memory space, so to access the data, "communication" is needed.
   *
   *         Iterators are assumed to be local rather than distributed, so methods that returns iterators are not provided.
   *         as an alternative, vectors are returned.
   *         methods that accept iterators as input assume that the input data is local.
   *
   *         Note that "communication" is a weak concept here meaning that we are accessing a different local container.
   *         as such, communicator may be defined for MPI, UPC, OpenMP, etc.
   *
   *         This allows the possibility of using distributed unordered map as local storage for coarser grain distributed container.
   *
   *         Note that communicator requires a mapping strategy between a key and the target processor/thread/partition.  The mapping
   *         may be done using a hash, similar to the local distributed unordered map, or it may be done via sorting/lookup or other mapping
   *         mechanisms.  The choice may be constrained by the communication approach, e.g. global sorting  does not work well with
   *         incremental async communication
   *
   * @tparam Key
   * @tparam T
   * @tparam Comm   default to mpi_collective_communicator       communicator for global communication. may hash or sort.
   * @tparam KeyTransform   transform function for the key.  can supply identity.  requires a single template argument (Key).  useful for mapping kmolecule to kmer.
   * @tparam Hash   hash function for local and distribution.  requires a template arugment (Key), and a bool (prefix, chooses the MSBs of hash instead of LSBs)
   * @tparam Equal   default to ::std::equal_to<Key>   equal function for the local storage.
   * @tparam Alloc  default to ::std::allocator< ::std::pair<Key, T> >    allocator for local storage.
   */
  template<typename Key, typename T,
  class Comm,
  template <typename> class KeyTransform,
  class Less = ::std::less<Key>,
  class Equal = ::std::equal_to<Key>,
  class Alloc = ::std::allocator< ::std::pair<Key, T> >
  >
  class sorted_multimap : public sorted_map_base<Key, T, ::std::vector, Comm, KeyTransform, Less, Equal, Alloc> {

    protected:
      using Base = sorted_map_base<Key, T, ::std::vector, Comm, KeyTransform, Less, Equal, Alloc>;


    public:
      using local_container_type = typename Base::local_container_type;

      // std::unordered_multimap public members.
      using key_type              = Key;
      using mapped_type           = T;
      using value_type            = ::std::pair<Key, T>;
      using iterator              = typename local_container_type::iterator;
      using const_iterator        = typename local_container_type::const_iterator;
      using size_type             = typename local_container_type::size_type;
      using difference_type       = typename local_container_type::difference_type;

    protected:

      // defined Communicator as a friend
      friend Comm;

      struct LocalFind {
          template<bool linear, class DBIter, typename Query, class OutputIter>
          size_t operator()(DBIter &range_begin, DBIter &el_end, DBIter const &range_end, Query const &v, OutputIter &output) const {
              // TODO: LINEAR SEARCH, O(n).  vs BINARY SEARCH, O(mlogn)


              range_begin = Base::template lower_bound<linear>(el_end, range_end, v);

              el_end = Base::template upper_bound<linear>(range_begin, range_end, v);

              // difference between the 2 iterators is the part that's equal.
              if (range_begin == range_end) return 0;

              // add the output entry, if found.
              output = ::std::copy(range_begin, el_end, output);

              return ::std::distance(range_begin, el_end);
          }
          template<bool linear, class DBIter, typename Query, class OutputIter, class Predicate = Identity>
          size_t operator()(DBIter &range_begin, DBIter &el_end, DBIter const &range_end,
        		  Query const &v, OutputIter &output, Predicate const &pred) const {
              // TODO: LINEAR SEARCH, O(n).  vs BINARY SEARCH, O(mlogn)
        	  //OutputIter output_orig = output;

              range_begin = Base::template lower_bound<linear>(el_end, range_end, v);

              el_end = Base::template upper_bound<linear>(range_begin, range_end, v);

              // difference between the 2 iterators is the part that's equal.
              if (range_begin == range_end) return 0;

              // add the output entry, if found.
              size_t count = 0;
              for (auto it = range_begin; it != el_end; ++it) {
                if (pred(*it)) {
                  *output = *it;
                  ++output;
                  ++count;
                }
              }

              return count;
              //output = ::std::copy_if(range_begin, el_end, output, pred);
              // return ::std::distance(output_orig, output);
          }
      } find_element;


    public:


      sorted_multimap(const mxx::comm& _comm) : Base(_comm) {
        this->key_multiplicity = 50;
      }

      virtual ~sorted_multimap() {}


      // default is for multimap
      virtual void rehash() {
        BL_BENCH_INIT(rehash);


        if (this->comm_size > 1) {
          // first balance

          // TODO: stable_block_decompose uses all2all internally.  is it better to move the deltas ourselves?
          if (!this->balanced) {
            BL_BENCH_START(rehash);
            ::mxx::stable_distribute(this->c, this->comm).swap(this->c);

            BL_BENCH_END(rehash, "block1", this->c.size());
          }

          // sort if needed
          if (!this->globally_sorted) {
            BL_BENCH_START(rehash);
            ::mxx::sort(this->c.begin(), this->c.end(), Base::Base::less, this->comm);
            BL_BENCH_END(rehash, "mxxsort", this->c.size());
          }

          // get new pivots
          BL_BENCH_START(rehash);
          this->key_to_rank.map.clear();
          if ((this->comm.rank() > 0) && (this->c.size() > 0)) {
            // only send for the first p-1 proc, and only if they have a kmer to split with.
            this->key_to_rank.map.emplace_back(this->c.front().first, this->comm.rank() - 1);
          }

          // using gatherv because some processes may not be sending anything.
          this->key_to_rank.map = ::mxx::allgatherv(this->key_to_rank.map, this->comm);

          assert(this->key_to_rank.map.size() > 0);

          // note that key_to_rank.map needs to be unique.
          auto map_end = std::unique(this->key_to_rank.map.begin(), this->key_to_rank.map.end(),
                                     Base::Base::equal);
          this->key_to_rank.map.erase(map_end, this->key_to_rank.map.end());

          if (this->comm.rank() == 0)
            for (size_t i = 0; i < this->key_to_rank.map.size(); ++i) {
              BL_DEBUGF("R %d unique key_to_rank.map[%lu] = (%s->%d)", this->comm.rank(), i, bliss::utils::KmerUtils::toASCIIString(this->key_to_rank.map[i].first).c_str() , this->key_to_rank.map[i].second );
            }
          assert(this->key_to_rank.map.size() > 0);


//          for (int i = 0; i < this->key_to_rank.map.size(); ++i) {
//            printf("R %d key to rank %s -> %d\n", this->comm.rank(), this->key_to_rank.map[i].first.toAlphabetString().c_str(), this->key_to_rank.map[i].second);
//          }
          BL_BENCH_END(rehash, "splitter1", this->key_to_rank.map.size());

          //			  mxx::datatype<::std::pair<Key, T> > dt;
          //			  MPI_Datatype mpi_dt = dt.type();
          //			  this->key_to_rank.map = ::mxx::sample_block_decomp(d.begin(), d.end(), Base::Base::less, this->comm_size - 1, this->comm, mpi_dt);

          // redistribute.  in trouble if we have a value that takes up a large part of the partition.
          BL_BENCH_START(rehash);
          ::std::vector<size_t> send_counts = fsc::get_bucket_sizes(this->c, this->key_to_rank.map, Base::Base::less);
          BL_BENCH_END(rehash, "bucket", this->c.size());

          // this should always be true.
          assert(send_counts.size() == static_cast<size_t>(this->comm_size));

          for (size_t i = 0; i < send_counts.size(); ++i) {
            BL_DEBUGF("R %d send_counts[%lu] = %lu", this->comm.rank(), i, send_counts[i]);
          }

          BL_BENCH_COLLECTIVE_START(rehash, "a2a", this->comm);
          // TODO: readjust boundaries using all2all.  is it better to move the deltas ourselves?
          this->c = mxx::all2allv(this->c, send_counts, this->comm);
          BL_BENCH_END(rehash, "a2a", this->c.size());

        } else {
          BL_BENCH_START(rehash);
          this->local_sort();
          BL_BENCH_END(rehash, "stdsort", this->c.size());
        }
        this->sorted = true; this->balanced = true; this->globally_sorted = true;

        BL_BENCH_REPORT_MPI(rehash, this->comm.rank(), this->comm);

      }


      /// update the multiplicity.  only multimap needs to do this.
      virtual size_t update_multiplicity() {

        // one approach is to add up the number of repeats for the key of each entry, then divide by total count.
        //  sum(count per key) / c.size.
        // problem with this approach is that for unordered map, to get the count for a key is essentially O(count), so we get quadratic time.
        // The approach is VERY SLOW for large repeat count.  - (0.0078125 human: 52 sec, synth: FOREVER.)

        // a second approach is to count the number of unique key then divide the map size by that.
        //  c.size / #unique.  requires unique set
        // To find unique set, we take each bucket, copy to vector, sort it, and then count unique.
        // This is precise, and is faster than the approach above.  (0.0078125 human: 54 sec.  synth: 57sec.)
        // but the n log(n) sort still grows with the duplicate count

        this->rehash();
/*
        size_t uniq_count = 0;
        ::std::pair<Key, T> v;
        for (auto it = this->c.begin(), max = this->c.end(); it != max;) {
          v = *it;
          ++uniq_count;
          it = Base::template upper_bound<true>(it, max, v);
        }
        if (uniq_count == 0)
          this->key_multiplicity = 1;
        else
          this->key_multiplicity = (this->c.size() + uniq_count - 1) / uniq_count + 1;
        //printf("%lu elements, %lu unique, key multiplicity = %lu\n", this->c.size(), uniq_count, this->key_multiplicity);
*/

        //        // third approach is to assume each bucket contains only 1 kmer/kmolecule.
        //        // This is not generally true for all hash functions, so this is an over estimation of the repeat count.
        //        // we equate bucket size to the number of repeats for that key.
        //        // we can use mean, max, or mean+stdev.
        //        // max overestimates significantly with potentially value > 1000, so don't use max.  (0.0078125 human: 50 sec. synth  32 sec)
        //        // mean may be underestimating for well behaving hash function.   (0.0078125 human: 50 sec. synth  32 sec)
        //        // mean + 2 stdev gets 95% of all entries.  1 stdev covers 67% of all entries, which for high coverage genome is probably better.
        //        //    (1 stdev:  0.0078125 human: 49 sec. synth  32 sec;  2stdev: 0.0078125 human 49s synth: 33 sec)

        // finally, hard coding.  (0.0078125 human:  50 sec.  synth:  32 s)
        // this->key_multiplicity = 50;

        return this->key_multiplicity;
      }

      template <class Predicate = Identity>
      ::std::vector<::std::pair<Key, T> > find(::std::vector<Key>& keys, bool sorted_input = false,
    		  Predicate const& pred = Predicate()) const {
          return Base::find(find_element, keys, sorted_input, pred);

/*         // DEBUG

          auto keys2 = keys;
          auto result = Base::find(find_element, keys, sorted_input, pred);
          auto result_a2a =  Base::find_a2a(find_element, keys2, sorted_input, pred);


          // DEBUG
          if (result.size() != result_a2a.size()) {
              throw ::std::logic_error("ERROR: not same size");

          }

          ::std::stable_sort(result.begin(), result.end(), this->less);
          ::std::stable_sort(result_a2a.begin(), result_a2a.end(), this->less);


          for (int i = 0; i < result_a2a.size(); ++i) {
            if (!this->equal(result[i], result_a2a[i])) {
              printf("rank %d failing at %d:  result: %s, result_a2a: %s\n", this->comm.rank(), i, result[i].first.toAlphabetString().c_str(), result_a2a[i].first.toAlphabetString().c_str());
              if ( i > 0)
                printf("rank %d   before   %d:  result: %s, result_a2a: %s\n", this->comm.rank(), i-1, result[i-1].first.toAlphabetString().c_str(), result_a2a[i-1].first.toAlphabetString().c_str());
              if (i < result_a2a.size() - 1)
                printf("rank %d   after    %d:  result: %s, result_a2a: %s\n", this->comm.rank(), i+1, result[i+1].first.toAlphabetString().c_str(), result_a2a[i+1].first.toAlphabetString().c_str());
              throw ::std::logic_error("ERROR: not same.");
            }
          }
//          bool same = ::std::equal(result.begin(), result.end(), result_a2a.begin(), typename Base::Base::TransformedEqual());
//          if (!same) throw ::std::logic_error("ERROR: not same.");

          return result;
*/
      }

      template <class Predicate = Identity>
      ::std::vector<::std::pair<Key, T> > find(Predicate const& pred = Predicate()) const {
          return Base::find(find_element, pred);
      }


  };



  /**
   * @brief  distributed unordered reduction map following std unordered map's interface.  Insertion applies the binary reduction operator between the existing and inserted element (in that order).
   * @details   This class is modeled after the std::unordered_map, but allows a binary reduction operator to be used during insertion.
   *
   *         the reduction operator is not assumed to be associative.  The operator is called with parameters existing element, then new element to insert.
   *
   *         it has as much of the same methods of std::unordered_map as possible.  however, all methods consider the fact
   *         that the data are in distributed memory space, so to access the data, "communication" is needed.
   *
   *         Note that "communication" is a weak concept here meaning that we are accessing a different local container.
   *         as such, communicator may be defined for MPI, UPC, OpenMP, etc.
   *
   *         This allows the possibility of using distributed unordered map as local storage for coarser grain distributed container.
   *
   *         Note that communicator requires a mapping strategy between a key and the target processor/thread/partition.  The mapping
   *         may be done using a hash, similar to the local distributed unordered map, or it may be done via sorting/lookup or other mapping
   *         mechanisms.  The choice may be constrained by the communication approach, e.g. global sorting  does not work well with
   *         incremental async communication
   *
   * @tparam Key
   * @tparam T
   * @tparam Comm   default to mpi_collective_communicator       communicator for global communication. may hash or sort.
   * @tparam KeyTransform   transform function for the key.  can supply identity.  requires a single template argument (Key).  useful for mapping kmolecule to kmer.
   * @tparam Hash   hash function for local and distribution.  requires a template arugment (Key), and a bool (prefix, chooses the MSBs of hash instead of LSBs)
   * @tparam Reduc  default to ::std::plus<key>    reduction operator
   * @tparam Equal   default to ::std::equal_to<Key>   equal function for the local storage.
   * @tparam Alloc  default to ::std::allocator< ::std::pair<Key, T> >    allocator for local storage.
   */
  template<typename Key, typename T,
  class Comm,
  template <typename> class KeyTransform,
  class Less = ::std::less<Key>,
  typename Reduc = ::std::plus<T>,
  class Equal = ::std::equal_to<Key>,
  class Alloc = ::std::allocator< ::std::pair<Key, T> >
  >
  class reduction_sorted_map : public sorted_map<Key, T, Comm, KeyTransform, Less, Equal, Alloc> {
      static_assert(::std::is_arithmetic<T>::value, "mapped type has to be arithmetic");

    protected:
      using Base = sorted_map<Key, T, Comm, KeyTransform, Less, Equal, Alloc>;


    public:
      using local_container_type = typename Base::local_container_type;

      // std::unordered_multimap public members.
      using key_type              = Key;
      using mapped_type           = T;
      using value_type            = ::std::pair<Key, T>;
      using iterator              = typename local_container_type::iterator;
      using const_iterator        = typename local_container_type::const_iterator;
      using size_type             = typename local_container_type::size_type;
      using difference_type       = typename local_container_type::difference_type;

    protected:
      Reduc r;

      // defined Communicator as a friend
      friend Comm;


      //=== below does local insert while reducing. doing reduction after, and letting insert be simple appending vector has lower complexity overall.
      //      /**
      //       * @brief insert new elements in the distributed sorted_map.  example use: stop inserting if more than x entries.
      //       * @param first
      //       * @param last
      //       */
      //      template <class InputIterator, class Predicate>
      //      void local_insert_if(InputIterator first, InputIterator last, Predicate const &pred, bool sorted_input = false) {
      //          if (first == last) return;
      //
      //          this->local_rehash();
      //
      //          // has to be locally reduced first.
      //          auto newlast = local_reduction(first, last, sorted_input);
      //
      //          // walk through.  if exiting entry, remove from input.
      //          // now walk through the input to count stuff.
      //          auto t_start = ::std::lower_bound(this->c.begin(), this->c.end(), *first, Base::Base::less);
      //          auto t_end = ::std::upper_bound(t_start, this->c.end(), *(newlast - 1), Base::Base::less);
      //
      //          // go through the input and search the container.
      //          // iterate through the input and search in output,  moving items up as we go.
      //          auto end = first;
      //          for (auto it = first; it != newlast; ++it) {  // walk though all input entries
      //              auto v = *it;
      //
      //              if (pred(v)) {
      //                // find entry in container.
      //                t_start = this->template lower_bound<true>(t_start, t_end, v);
      //
      //                if ((t_start == t_end) || (!this->equal(v, *t_start))) {
      //                  // not matched.  copy the value, move pos up by 1.
      //                  if (end != it) *end = v;
      //                  ++end;
      //                } else {
      //                  // matched.  so need to reduce
      //                  t_start->second = r(t_start->second, v.second);
      //                }
      //              }
      //            }
      //          // insert at the end.
      //          if (first != end) {
      //            this->c.insert(this->c.end(), first, end);
      //            this->sorted = false;
      //          }
      //
      //      }

      virtual void local_reduction(::std::vector<::std::pair<Key, T> >& input, bool sorted_input = false) {
        if (input.size() == 0) return;

        auto b = input.begin();
        auto e = input.end();
        if (!sorted_input) Base::Base::sort_ascending(b, e);

        // then do reduction
        auto reduc_target = b;
        auto curr = b;  ++curr;
        while (curr != e) {
          if (this->equal(*reduc_target, *curr))  // if same, do reduction
            reduc_target->second = r(reduc_target->second, curr->second);
          else {  // else reset b.
            ++reduc_target;
            *reduc_target = *curr;
          }
          ++curr;  // increment second.
        }
        ++reduc_target;
        input.erase(reduc_target, input.end());
      }


    public:


      reduction_sorted_map(const mxx::comm& _comm) : Base(_comm) {}

      virtual ~reduction_sorted_map() {};

  };


  /**
   * @brief  distributed unordered counting map following std unordered map's interface.  Insertion applies the binary reduction operator between the existing and inserted element (in that order).
   * @details   This class is modeled after the std::unordered_map, but allows a binary reduction operator to be used during insertion.
   *
   *         the reduction operator is not assumed to be associative.  The operator is called with parameters existing element, then new element to insert.
   *
   *         it has as much of the same methods of std::unordered_map as possible.  however, all methods consider the fact
   *         that the data are in distributed memory space, so to access the data, "communication" is needed.
   *
   *         Note that "communication" is a weak concept here meaning that we are accessing a different local container.
   *         as such, communicator may be defined for MPI, UPC, OpenMP, etc.
   *
   *         This allows the possibility of using distributed unordered map as local storage for coarser grain distributed container.
   *
   *         Note that communicator requires a mapping strategy between a key and the target processor/thread/partition.  The mapping
   *         may be done using a hash, similar to the local distributed unordered map, or it may be done via sorting/lookup or other mapping
   *         mechanisms.  The choice may be constrained by the communication approach, e.g. global sorting  does not work well with
   *         incremental async communication
   *
   * @tparam Key
   * @tparam T
   * @tparam Comm   default to mpi_collective_communicator       communicator for global communication. may hash or sort.
   * @tparam KeyTransform   transform function for the key.  can supply identity.  requires a single template argument (Key).  useful for mapping kmolecule to kmer.
   * @tparam Hash   hash function for local and distribution.  requires a template arugment (Key), and a bool (prefix, chooses the MSBs of hash instead of LSBs)
   * @tparam Equal   default to ::std::equal_to<Key>   equal function for the local storage.
   * @tparam Alloc  default to ::std::allocator< ::std::pair<Key, T> >    allocator for local storage.
   */
  template<typename Key, typename T,
  class Comm,
  template <typename> class KeyTransform,
  class Less = ::std::less<Key>,
  class Equal = ::std::equal_to<Key>,
  class Alloc = ::std::allocator< ::std::pair<Key, T> >
  >
  class counting_sorted_map : public reduction_sorted_map<Key, T, Comm, KeyTransform, Less, ::std::plus<T>, Equal,Alloc> {
      static_assert(::std::is_integral<T>::value, "count type has to be integral");

    protected:
      using Base = reduction_sorted_map<Key, T, Comm, KeyTransform, Less, ::std::plus<T>, Equal, Alloc>;

    public:
      using local_container_type = typename Base::local_container_type;

      // std::unordered_multimap public members.
      using key_type              = Key;
      using mapped_type           = T;
      using value_type            = ::std::pair<Key, T>;
      using iterator              = typename local_container_type::iterator;
      using const_iterator        = typename local_container_type::const_iterator;
      using size_type             = typename local_container_type::size_type;
      using difference_type       = typename local_container_type::difference_type;


    protected:

      // defined Communicator as a friend
      friend Comm;


    public:
      counting_sorted_map(const mxx::comm& _comm) : Base(_comm) {}

      virtual ~counting_sorted_map() {};


      /**
       * @brief insert new elements in the distributed sorted_multimap.
       * @param first
       * @param last
       */
      template <class Predicate = Identity>
      size_t insert(::std::vector<Key> &input, bool sorted_input = false, Predicate const &pred = Predicate()) {

        // even if count is 0, still need to participate in mpi calls.  if (input.size() == 0) return;
        BL_BENCH_INIT(count_insert);

        BL_BENCH_START(count_insert);
        ::std::vector<::std::pair<Key, T> > temp;
        temp.reserve(input.size());
        ::fsc::back_emplace_iterator<::std::vector<::std::pair<Key, T> > > emplace_iter(temp);
        ::std::transform(input.begin(), input.end(), emplace_iter, [](Key const & x) { return ::std::make_pair(x, T(1)); });
        BL_BENCH_END(count_insert, "convert", input.size());

        // distribute
        BL_BENCH_START(count_insert);
        // local compute part.  called by the communicator.
        size_t count = this->Base::insert(temp, sorted_input, pred);
        ::std::vector<::std::pair<Key, T> >().swap(temp);  // clear the temp.

        BL_BENCH_END(count_insert, "insert", this->c.size());

        // distribute
        BL_BENCH_REPORT_MPI(count_insert, this->comm.rank(), this->comm);
        return count;
      }

      /**
       * @brief insert new elements in the distributed sorted_multimap.
       * @param first
       * @param last
       */
      template <class Predicate = Identity>
      size_t insert(::std::vector<std::pair<Key, T> > &input, bool sorted_input = false, Predicate const &pred = Predicate()) {

        // even if count is 0, still need to participate in mpi calls.  if (input.size() == 0) return;
        BL_BENCH_INIT(count_insert);

        // distribute
        BL_BENCH_START(count_insert);
        // local compute part.  called by the communicator.
        size_t count = this->Base::insert(input, sorted_input, pred);

        BL_BENCH_END(count_insert, "insert", this->c.size());

        // distribute
        BL_BENCH_REPORT_MPI(count_insert, this->comm.rank(), this->comm);
        return count;
      }

  };



} /* namespace dsc */


#endif // BLISS_DISTRIBUTED_SORTED_MAP_HPP
