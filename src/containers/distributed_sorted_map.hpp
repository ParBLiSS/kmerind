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
  >
  class sorted_map_base : public ::dsc::map_base<Key, T, Comm, KeyTransform, Less, Equal, Alloc> {

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
       * @tparam unique 		indicates that only first of duplicated queries should be processed.
       */
      template <bool skip_duplicate_query = false>
      struct QueryProcessor {


          // get the overlapping range in container.  container must be sorted.
          template<typename QueryIter, typename DBIter>
          static ::std::pair<DBIter, DBIter> intersect(DBIter db_begin, DBIter db_end,
                                                       QueryIter query_begin, QueryIter query_end,
                                                       bool & sorted_query) {
            if (db_begin == db_end) return ::std::make_pair(db_begin, db_end);  // no target input
            if (query_begin == query_end) return ::std::make_pair(db_end, db_end);  // no src input

            auto dist_query = ::std::distance(query_begin, query_end);

            // make sure query is sorted sorted.
            if (!sorted_query) {
            	::std::sort(query_begin, query_end, Base::less);
            	sorted_query = true;
            }

            // find the bounds.
            auto range_begin = ::std::lower_bound(db_begin, db_end, *query_begin, Base::less);
            auto query_last = query_begin; ::std::advance(query_last, dist_query - 1);  // last real entry because query end may not be derefernceable.
            auto range_end = ::std::upper_bound(range_begin, db_end, *query_last, Base::less);

            return ::std::make_pair(range_begin, range_end);
          }

          // assumes that container is sorted. and exact overlap region is provided.  do not filter output here since it's an output iterator.
          template <class DBIter, class QueryIter, class OutputIter, class Operator, class Predicate = TruePredicate>
          static size_t process(DBIter range_begin, DBIter range_end,
                                QueryIter query_begin, QueryIter query_end,
                                OutputIter &output, Operator & op,
                                bool sorted_query = false, Predicate const &pred = Predicate()) {

              // no matches in container.
              if (range_begin == range_end) return 0;
              if (query_begin == query_end) return 0;  // no input

              //auto output_start = output;

              auto dist_range = ::std::distance(range_begin, range_end);
              auto dist_query = ::std::distance(query_begin, query_end);

              //if (!sorted_target) Base::sort_ascending(range_begin, range_end);  range_begin and range_end often are const iterators.
              if (!sorted_query)
            	  ::std::sort(query_begin, query_end, Base::less);

              auto el_end = range_begin;
              size_t count = 0;
              bool linear = (static_cast<double>(dist_query) * ::std::log2(dist_range)) > dist_range;
              typename ::std::iterator_traits<QueryIter>::value_type v;

              if (linear) {  // based on number of input and search source, choose a method to search.

                // iterate through the input and search in output -
                for (auto it = query_begin; it != query_end;) {
                  v = *it;
                  if (!::std::is_same<Predicate, TruePredicate>::value)
                	  count += op.template operator()<true>(range_begin, el_end, range_end, v, output, pred);
                  else
                	  count += op.template operator()<true>(range_begin, el_end, range_end, v, output);

                  // compiler optimize out the conditional.
                  if (skip_duplicate_query) it = ::fsc::upper_bound<true>(it, query_end, v, Base::less);
                  else ++it;
                }
              } else {
                // use logarithmic search

                // iterate through the input and search in output -
                for (auto it = query_begin; it != query_end;) {
                  v = *it;
                  if (!::std::is_same<Predicate, TruePredicate>::value)
                	  count += op.template operator()<false>(range_begin, el_end, range_end, v, output, pred);
                  else
                	  count += op.template operator()<false>(range_begin, el_end, range_end, v, output);

                  // compiler optmizes out the conditional
                  if (skip_duplicate_query) it = ::fsc::upper_bound<true>(it, query_end, v, Base::less);
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

    private:
      // ========== collective variables
      /**
       * @brief flag indicating if globally balanced.  changed by insert(), delete(), and redistribute().
       * @note  independent of sorted or globally sorted flags.
       */
      mutable bool balanced;   // this should be considered a global variable

      /**
       * @brief  flag indicating if globally sorted.  changed by insert(), delete(), and redistribute()
       * @note   globally sorted implies sorted.  converse not true
       */
      mutable bool globally_sorted;   // this should be considered a global variable.


    protected:
      local_container_type c;

      // ========== individual variables
      /**
       * @brief  flag indicating if local container is sorted.  changed by insert(), and any call to local_sort() (in all query functions).
       * @note    not sorted implies globally not sorted.  converse not true.
       */
      bool sorted;   // this is a local variable.


      // =========== accessors to change the local state of the container
      void set_balanced(bool v) const {
        balanced = v;
      }
      // =========== accessors to change the local state of the container
      void set_globally_sorted(bool v) const {
        globally_sorted = v;
      }

      // =========== collective operations to get distribution state of the container
      /**
       * @brief is the sorted vector balanced?
       */
      bool is_balanced() const {
        if (this->comm.size() == 1)
          return this->balanced;
        else // all reduce
          return mxx::all_of(this->balanced, this->comm);
      }

      bool is_globally_sorted() const {
        if (this->comm.size() == 1)
          return globally_sorted;
        else // all reduce
          globally_sorted = mxx::all_of(globally_sorted, this->comm);
          return globally_sorted;
      }


      // ============ shared functors.

      struct LocalCount {
          // unfiltered.
          template<bool linear, class DBIter, typename Query, class OutputIter>
          size_t operator()(DBIter &range_begin, DBIter &el_end, DBIter const &range_end, Query const &v, OutputIter &output) const {
              // TODO: LINEAR SEARCH, O(n).  vs BINARY SEARCH, O(mlogn)
              range_begin = ::fsc::lower_bound<linear>(el_end, range_end, v, Base::less);  // range_begin at equal or greater than v.
              el_end = ::fsc::upper_bound<linear>(range_begin, range_end, v, Base::less);  // el_end at greater than v.
              // difference between the 2 iterators is the part that's equal.

              // add the output entry.
              size_t count = ::std::distance(range_begin, el_end);

              *output = ::std::move(::std::make_pair(v, count));
              ++output;
              return 1;
          }
          // filtered element-wise.
          template<bool linear, class DBIter, typename Query, class OutputIter, class Predicate = TruePredicate>
          size_t operator()(DBIter &range_begin, DBIter &el_end, DBIter const &range_end, Query const &v, OutputIter &output,
                            Predicate const& pred) const {
              // TODO: LINEAR SEARCH, O(n).  vs BINARY SEARCH, O(mlogn)
              range_begin = ::fsc::lower_bound<linear>(el_end, range_end, v, Base::less);  // range_begin at equal or greater than v.
              el_end = ::fsc::upper_bound<linear>(range_begin, range_end, v, Base::less);  // el_end at greater than v.
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
              range_begin = ::fsc::lower_bound<linear>(el_end, range_end, v, Base::less);

              // if the keep range is larger than 0, then move data and update insert pos.
              auto dist = ::std::distance(el_end, range_begin);
              if (dist > 0) {
                // condense a portion.
                ::std::move(el_end, range_begin, output);
                ::std::advance(output, dist);
              }
              // find of end of the segment to delete == start of next segment to keep
              el_end = ::fsc::upper_bound<linear>(range_begin, range_end, v, Base::less);
              return dist;
          }
          /// Return how much was KEPT.
          template<bool linear, class DBIter, typename Query, class Predicate = TruePredicate>
          size_t operator()(DBIter &range_begin, DBIter &el_end, DBIter const &range_end, Query const &v, DBIter &output,
                            Predicate const & pred) {
              // find start of segment to delete == end of prev segment to keep
              range_begin = ::fsc::lower_bound<linear>(el_end, range_end, v, Base::less);

              if (::std::distance(el_end, range_begin) > 0) {
                // condense a portion.
                output = ::std::move(el_end, range_begin, output);
              }
              // find of end of the segment to delete == start of next segment to keep
              el_end = ::fsc::upper_bound<linear>(range_begin, range_end, v, Base::less);


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


      // ======================== global map functions that requires subclass defined functors

      /**
       * @brief find elements with the specified keys in the distributed sorted_multimap.
       * @param first
       * @param last
       */
      template <class LocalFind, class Predicate = TruePredicate >
      ::std::vector<::std::pair<Key, T> > find_a2a(LocalFind & lf,
    		  ::std::vector<Key>& keys, bool & sorted_input,
    		  Predicate const& pred = Predicate() ) const {

          BL_BENCH_INIT(find);
          ::std::vector<::std::pair<Key, T> > results;

          if (this->empty()) {
              BL_BENCH_REPORT_MPI_NAMED(find, "base_sorted_map:find_a2a", this->comm);
              return results;
          }


          BL_BENCH_START(find);
          // even if count is 0, still need to participate in mpi calls.  if (keys.size() == 0) return results;
          ::fsc::back_emplace_iterator<::std::vector<::std::pair<Key, T> > > emplace_iter(results);
          BL_BENCH_END(find, "begin", keys.size());

          if (this->comm.size() > 1) {
              // ensure that the container splitters are setup properly, and load balanced.  this also ensures that locally we are sorted.
              BL_BENCH_COLLECTIVE_START(find, "global_sort", this->comm);
              this->redistribute();
              BL_BENCH_END(find, "global_sort", this->local_size());

              BL_BENCH_COLLECTIVE_START(find, "dist_query", this->comm);
              // distribute (communication part)
              std::vector<size_t> recv_counts(::dsc::distribute_sorted_unique(keys,
                                                                              this->key_to_rank, sorted_input, this->comm,
                                                                              typename Base::TransformedLess(), typename Base::TransformedEqual()));
              BL_BENCH_END(find, "dist_query", keys.size());


            // local find. memory utilization a potential problem.
            // do for each src proc one at a time.
            BL_BENCH_START(find);
            results.reserve(keys.size() * this->get_multiplicity());                   // TODO:  should estimate coverage.
            BL_BENCH_END(find, "reserve", results.capacity());


            BL_BENCH_START(find);
            std::vector<size_t> send_counts(this->comm.size(), 0);

            auto start = keys.begin();
            auto end = start;
            for (int i = 0; i < this->comm.size(); ++i) {
              ::std::advance(end, recv_counts[i]);  // recv queries are unique in each bucket.

              // work on query from process i.  specify no skip_duplicate.
              auto overlap = QueryProcessor<false>::intersect(this->c.begin(), this->c.end(), start, end, sorted_input);

              // within start-end, values are unique, so don't need to set unique to true.
              send_counts[i] = QueryProcessor<false>::process(overlap.first, overlap.second,
            		  start, end, emplace_iter, lf, sorted_input, pred);

              start = end;
            }
            BL_BENCH_END(find, "local_find", results.size());

            // send back using the constructed recv count
            BL_BENCH_COLLECTIVE_START(find, "a2a2", this->comm);
            results = mxx::all2allv(results, send_counts, this->comm);
            BL_BENCH_END(find, "a2a2", results.size());

          } else {
        	  // ensure data is sorted locally.
              // ensure that the container splitters are setup properly, and load balanced.
              BL_BENCH_COLLECTIVE_START(find, "local_sort", this->comm);
              this->local_sort();
              BL_BENCH_END(find, "local_sort", this->local_size());


              // keep unique keys
              BL_BENCH_START(find);
              ::fsc::sorted_unique(keys, sorted_input, typename Base::TransformedLess(), typename Base::TransformedEqual());
              BL_BENCH_END(find, "uniq1", keys.size());

        	  BL_BENCH_START(find);
              results.reserve(keys.size() * this->get_multiplicity());                   // TODO:  should estimate coverage.
              BL_BENCH_END(find, "reserve", results.capacity());


            BL_BENCH_START(find);
            auto overlap = QueryProcessor<false>::intersect(this->c.begin(), this->c.end(),
            		keys.begin(), keys.end(), sorted_input);

            // within start-end, values are unique, so don't need to set unique to true.
            QueryProcessor<false>::process(overlap.first, overlap.second, keys.begin(), keys.end(),
            		emplace_iter, lf, sorted_input, pred);

            BL_BENCH_END(find, "local_find", results.size());

          }

          BL_BENCH_REPORT_MPI_NAMED(find, "base_sorted_map:find_a2a", this->comm);

          return results;
      }

      /**
       * @brief find elements with the specified keys in the distributed sorted_multimap.
       *
       * process request from a single processor at a time.
       *
       * @param first
       * @param last
       */
      template <class LocalFind, class Predicate = TruePredicate >
      ::std::vector<::std::pair<Key, T> > find(LocalFind & lf, ::std::vector<Key>& keys, bool sorted_input = false,
          Predicate const& pred = Predicate() ) const {
          BL_BENCH_INIT(find);
          ::std::vector<::std::pair<Key, T> > results;

          if (this->empty()) {
              BL_BENCH_REPORT_MPI_NAMED(find, "base_sorted_map:find", this->comm);
              return results;
          }

          BL_BENCH_START(find);
          // even if count is 0, still need to participate in mpi calls.  if (keys.size() == 0) return results;
          ::fsc::back_emplace_iterator<::std::vector<::std::pair<Key, T> > > emplace_iter(results);

          ::std::vector<::std::pair<Key, T> > local_results;
          ::fsc::back_emplace_iterator<::std::vector<::std::pair<Key, T> > > local_emplace_iter(local_results);
          BL_BENCH_END(find, "begin", keys.size());


          if (this->comm.size() > 1) {

              // ensure that the container splitters are setup properly, and load balanced.
              BL_BENCH_COLLECTIVE_START(find, "global_sort", this->comm);
              this->redistribute();
              BL_BENCH_END(find, "global_sort", this->local_size());

              BL_BENCH_COLLECTIVE_START(find, "dist_query", this->comm);
              // distribute (communication part)
              std::vector<size_t> recv_counts(::dsc::distribute_sorted_unique(keys,
                                                                              this->key_to_rank, sorted_input, this->comm,
                                                                              typename Base::TransformedLess(), typename Base::TransformedEqual()));
              BL_BENCH_END(find, "dist_query", keys.size());


            //====== local count to determine amount of memory to allocate at destination.
            BL_BENCH_START(find);

            ::std::vector<::std::pair<Key, size_t> > count_results;
            size_t max_key_count = *(::std::max_element(recv_counts.begin(), recv_counts.end()));
            count_results.reserve(max_key_count);
            ::fsc::back_emplace_iterator<::std::vector<::std::pair<Key, size_t> > > count_emplace_iter(count_results);
            std::vector<size_t> send_counts(this->comm.size(), 0);

            auto start = keys.begin();
            auto end = start;
            size_t total = 0;
            for (int i = 0; i < this->comm.size(); ++i) {
              ::std::advance(end, recv_counts[i]);

              // count results for process i
              count_results.clear();
              auto overlap = QueryProcessor<false>::intersect(this->c.begin(), this->c.end(), start, end,
            		  sorted_input);
              QueryProcessor<false>::process(overlap.first, overlap.second, start, end,
            		  count_emplace_iter, count_element, sorted_input, pred);
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


            //==== reserve

            BL_BENCH_START(find);
            auto resp_displs = mxx::impl::get_displacements(resp_counts);  // compute response displacements.

            auto resp_total = resp_displs[this->comm.size() - 1] + resp_counts[this->comm.size() - 1];
            auto max_send_count = *(::std::max_element(send_counts.begin(), send_counts.end()));
            results.resize(resp_total);   // allocate, not just reserve
            local_results.reserve(max_send_count);
            BL_BENCH_END(find, "reserve", resp_total);

            //=== process queries and send results.  O(p) iterations
            BL_BENCH_START(find);
            auto recv_displs = mxx::impl::get_displacements(recv_counts);  // compute response displacements.
            int recv_from, send_to;
            size_t found;
            total = 0;
            std::vector<MPI_Request> reqs(2 * this->comm.size());

            mxx::datatype dt = mxx::get_datatype<::std::pair<Key, T>>();
            for (int i = 0; i < this->comm.size(); ++i) {
              recv_from = (this->comm.rank() + (this->comm.size() - i)) % this->comm.size(); // rank to recv data from

              // set up receive.
              MPI_Irecv(&results[resp_displs[recv_from]], resp_counts[recv_from], dt.type(),
                        recv_from, i, this->comm, &reqs[2 * i]);


              send_to = (this->comm.rank() + i) % this->comm.size();    // rank to send data to

              //== get data for the dest rank
              start = keys.begin();                                   // keys for the query for the dest rank
              ::std::advance(start, recv_displs[send_to]);
              end = start;
              ::std::advance(end, recv_counts[send_to]);

              local_results.clear();
              // work on query from process i.
              auto overlap = QueryProcessor<false>::intersect(this->c.begin(), this->c.end(),
            		  start, end, sorted_input);
              found = QueryProcessor<false>::process(overlap.first, overlap.second,
            		  start, end, local_emplace_iter, lf, sorted_input, pred);
              total += found;
              //== now send the results immediately - minimizing data usage so we need to wait for both send and recv to complete right now.

              MPI_Isend(&(local_results[0]), found, dt.type(), send_to,
                        i, this->comm, &reqs[2 * i + 1]);

              // wait for both requests to complete.
              MPI_Waitall(2, &reqs[2 * i], MPI_STATUSES_IGNORE);

              // within start-end, values are unique, so don't need to set unique to true.

            }

            BL_BENCH_END(find, "find_send", results.size());

          } else {
              // ensure that the container splitters are setup properly, and load balanced.
              BL_BENCH_COLLECTIVE_START(find, "local_sort", this->comm);
              this->local_sort();  // ensure data is locally sorted
              BL_BENCH_END(find, "local_sort", this->local_size());

              // keep unique keys
              BL_BENCH_START(find);
              ::fsc::sorted_unique(keys, sorted_input, typename Base::TransformedLess(), typename Base::TransformedEqual());
              BL_BENCH_END(find, "uniq1", keys.size());



        	  // memory is constrained.  find EXACT count.
            BL_BENCH_START(find);

            ::std::vector<::std::pair<Key, size_t> > count_results;
            count_results.reserve(keys.size());
            ::fsc::back_emplace_iterator<::std::vector<::std::pair<Key, size_t> > > count_emplace_iter(count_results);

            // count now.
            auto overlap = QueryProcessor<false>::intersect(this->c.begin(), this->c.end(),
            		keys.begin(), keys.end(), sorted_input);
            QueryProcessor<false>::process(overlap.first, overlap.second,
            		keys.begin(), keys.end(), count_emplace_iter, count_element, sorted_input, pred);
            size_t count = ::std::accumulate(count_results.begin(), count_results.end(), static_cast<size_t>(0),
                                          [](size_t v, ::std::pair<Key, size_t> const & x) {
                      return v + x.second;
                    });
            BL_BENCH_END(find, "local_count", count);

            BL_BENCH_START(find);
            results.reserve(count);  // 1 result per key.
            BL_BENCH_END(find, "reserve", count);


            BL_BENCH_START(find);
            // within start-end, values are unique, so don't need to set unique to true.
            QueryProcessor<false>::process(overlap.first, overlap.second,
            		keys.begin(), keys.end(), emplace_iter, lf, sorted_input, pred);
            BL_BENCH_END(find, "local_find", results.size());

          }

          BL_BENCH_REPORT_MPI_NAMED(find, "base_sorted_map:find_isend", this->comm);

          return results;
      }

      /// version using predicate, applies to entire container.
      template <class LocalFind, class Predicate = TruePredicate>
      ::std::vector<::std::pair<Key, T> > find(LocalFind & lf,
    		  Predicate const & pred = Predicate()) const {
          ::std::vector<::std::pair<Key, T> > results;

          // TODO: directly iterate over the vector?

          if (this->local_empty()) return results;

          this->local_sort();

          ::fsc::back_emplace_iterator<::std::vector<::std::pair<Key, T> > > emplace_iter(results);

          auto keys = this->keys();

          ::std::vector<::std::pair<Key, size_t> > count_results;
          count_results.reserve(keys.size());
          ::fsc::back_emplace_iterator<::std::vector<::std::pair<Key, size_t> > > count_emplace_iter(count_results);

          // count now.
          QueryProcessor<false>::process(this->c.begin(), this->c.end(),
        		  keys.begin(), keys.end(), count_emplace_iter, count_element, true, pred);
          size_t count = ::std::accumulate(count_results.begin(), count_results.end(), static_cast<size_t>(0),
                                           [](size_t v, ::std::pair<Key, size_t> const & x) {
                       return v + x.second;
                     });
          results.reserve(count);  // 1 result per key.

          // within start-end, values are unique, so don't need to set unique to true.
          QueryProcessor<false>::process(this->c.begin(), this->c.end(),
        		  keys.begin(), keys.end(), emplace_iter, lf, true, pred);

          if (this->comm.size() > 1) this->comm.barrier();
          return results;
      }


      // ======================== local operations
      /**
       * @brief insert new elements in the distributed sorted_multimap.  example use: stop inserting if more than x entries.
       * TODO: split this into insert (into empty), and an append (into existing)
       * @param first
       * @param last
       */
      template <class Predicate = TruePredicate>
      size_t local_insert(::std::vector<::std::pair<Key, T> > &input, bool sorted_input = false, Predicate const &pred = Predicate()) {
          if (input.size() == 0) return 0;  // OKAY HERE ONLY BECAUSE NO COMMUNICATION IS HERE.
          BL_BENCH_INIT(insert);

          size_t before = c.size();
          BL_BENCH_START(insert);

          ::fsc::back_emplace_iterator<local_container_type> emplace_iter(c);
          if (::std::is_same<Predicate, TruePredicate>::value) {
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


          BL_BENCH_REPORT_MPI_NAMED(insert, "base_sorted_map:insert", this->comm);

          this->sorted = false;

          return count;
      }




      // ==================== constructor

      /// constructor
      sorted_map_base(const mxx::comm& _comm) : Base(_comm),
          key_to_rank(_comm.size()), balanced(false), globally_sorted(false), sorted(false) {}

      // ===================  sorted map specific virtual functions
      /// ensures container is globally sorted/organized and balanced, and splitters are capatured.  also ensures local sortedness.
      virtual void redistribute() = 0;
      virtual void redistribute() const {
        const_cast<typename std::remove_cv<typename std::remove_reference<decltype(*this)>::type>::type *>(this)->redistribute();
      }

      // ===================== overrides

      /// clears the sorted_map
      virtual void local_clear() {
        c.clear();

        this->sorted = true;
        this->set_balanced(false);
      }

      /// reserve space.  n is the local container size.  this allows different processes to individually adjust its own size.
      virtual void local_reserve( size_t n) {
        // vector's reserve will only do something if n > capacity.
        c.reserve(n);
      }


      // ==================== sorted vector specific functions.

      /// rehash the local container.  n is the local container size.  this allows different processes to individually adjust its own size.
      void local_sort() {
        ::fsc::sort(c, sorted, Base::less);
      }

      /// const version that sorts the local container.
      void local_sort() const {
        const_cast<typename std::remove_cv<typename std::remove_reference<decltype(*this)>::type>::type *>(this)->local_sort();
      }



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

      /// extract the unique keys of a map.
      virtual void keys(std::vector<Key> & result) const {
        result.clear();
        if (this->local_empty()) return;

        // copy the keys
        auto end = c.end();
        for (auto it = c.begin(); it != end; ++it) {
          result.emplace_back(it->first);
        }

        // and then find unique.
        bool temp = this->sorted;
        ::fsc::sorted_unique(result, temp, typename Base::TransformedLess(), typename Base::TransformedEqual());

      }



      /**
       * @brief count elements with the specified keys in the distributed sorted_multimap.
       * @param first
       * @param last
       */
      template <typename Predicate = TruePredicate>
      ::std::vector<::std::pair<Key, size_type> > count(::std::vector<Key>& keys, bool sorted_input = false,
    		  Predicate const & pred = Predicate()) const {
        ::std::vector<::std::pair<Key, size_type> > results;


        BL_BENCH_INIT(count);

        if (this->empty()) {
            BL_BENCH_REPORT_MPI_NAMED(count, "base_sorted_map:count", this->comm);
            return results;
        }


        BL_BENCH_START(count);
        // keep unique keys
        // even if count is 0, still need to participate in mpi calls.  if (keys.size() == 0) return results;
        ::fsc::back_emplace_iterator<::std::vector<::std::pair<Key, size_type> > > emplace_iter(results);
        BL_BENCH_END(count, "begin", keys.size());


        if (this->comm.size() > 1) {
            // ensure that the container splitters are setup properly, and load balanced.
            BL_BENCH_COLLECTIVE_START(count, "global_sort", this->comm);
            this->redistribute();
            BL_BENCH_END(count, "global_sort", this->local_size());


            BL_BENCH_START(count);
            std::vector<size_t> recv_counts(::dsc::distribute_sorted_unique(keys,
                                                                            this->key_to_rank, sorted_input, this->comm,
                                                                            typename Base::TransformedLess(), typename Base::TransformedEqual()));
            BL_BENCH_END(count, "dist_query", keys.size());


          BL_BENCH_START(count);
          // local find. memory utilization a potential problem.
          // do for each src proc one at a time.
          results.reserve(keys.size());                   // TODO:  should estimate coverage.
          BL_BENCH_END(count, "reserve", results.capacity());

          BL_BENCH_START(count);
          auto start = keys.begin();
          auto end = start;
          for (int i = 0; i < this->comm.size(); ++i) {
            ::std::advance(end, recv_counts[i]);

            // work on query from process i.
            auto overlap = QueryProcessor<false>::intersect(this->c.begin(), this->c.end(), start, end, sorted_input);

            // within start-end, values are unique, so don't need to set unique to true.
            QueryProcessor<false>::process(overlap.first, overlap.second,
            		start, end, emplace_iter, count_element, sorted_input, pred);

            start = end;
          }
          BL_BENCH_END(count, "local_count", results.size());

          BL_BENCH_COLLECTIVE_START(count, "a2a2", this->comm);
          // send back using the constructed recv count
          results = mxx::all2allv(results, recv_counts, this->comm);
          BL_BENCH_END(count, "a2a2", results.size());


        } else {
            // ensure that the container splitters are setup properly, and load balanced.
            BL_BENCH_COLLECTIVE_START(count, "local_sort", this->comm);
            this->local_sort();
            BL_BENCH_END(count, "local_sort", this->local_size());

            BL_BENCH_START(count);
            // keep unique keys
            ::fsc::sorted_unique(keys, sorted_input, typename Base::TransformedLess(), typename Base::TransformedEqual());
            BL_BENCH_END(count, "uniq1", keys.size());

          BL_BENCH_START(count);
          results.reserve(keys.size());
          BL_BENCH_END(count, "reserve", results.capacity());

          BL_BENCH_START(count);
          // work on query from process i.
          auto overlap = QueryProcessor<true>::intersect(this->c.begin(), this->c.end(),
        		  keys.begin(), keys.end(), sorted_input);

          // within key, values may not be unique,
          QueryProcessor<true>::process(overlap.first, overlap.second,
        		  keys.begin(), keys.end(), emplace_iter, count_element, sorted_input, pred);
          BL_BENCH_END(count, "local_count", results.size());

        }
        BL_BENCH_REPORT_MPI_NAMED(count, "base_sorted_map:count", this->comm);

        return results;

      }


      template <typename Predicate = TruePredicate>
      ::std::vector<::std::pair<Key, size_type> > count(Predicate const & pred = Predicate()) const {
        ::std::vector<::std::pair<Key, size_type> > results;

        if (this->local_empty()) return results;

        // ensure that the container splitters are setup properly, and load balanced.
        this->local_sort();


        ::fsc::back_emplace_iterator<::std::vector<::std::pair<Key, size_type> > > emplace_iter(results);

        auto keys = this->keys();
        results.reserve(keys.size());

        // keys already unique
        QueryProcessor<false>::process(c.begin(), c.end(), keys.begin(), keys.end(),
        		emplace_iter, count_element, true, pred);

        if (this->comm.size() > 1) this->comm.barrier();
        return results;
      }


//      /**
//       * @brief insert new elements in the distributed sorted_multimap.  example use: stop inserting if more than x entries.
//       * @param src_begin
//       * @param src_end
//       */
//      template <class InputIter, class Predicate = TruePredicate>
//      size_t insert(InputIter src_begin, InputIter src_end, bool sorted_input = false, Predicate const &pred = Predicate()) {
//          if (src_begin == src_end) return 0;
//          BL_BENCH_INIT(insert);
//
//          BL_BENCH_START(insert);
//
//          this->sorted = false; this->balanced = false; this->globally_sorted = false;
//          ::fsc::back_emplace_iterator<local_container_type> emplace_iter(c);
//
//          if (::std::is_same<Predicate, TruePredicate>::value) ::std::copy(src_begin, src_end, emplace_iter);
//          else ::std::copy_if(src_begin, src_end, emplace_iter, pred);
//
//          size_t count = ::std::distance(src_begin, src_end);
//          BL_BENCH_END(insert, "insert", count);
//
//
//          BL_BENCH_REPORT_MPI_NAMED(insert, "map:", this->comm);
//
//          return count;
//      }

      /**
       * @brief insert new elements in the distributed sorted_multimap.  example use: stop inserting if more than x entries.
       * TODO: split this into insert (into empty), and an append (into existing)
       * @param first
       * @param last
       */
      template <class Predicate = TruePredicate>
      size_t insert(::std::vector<::std::pair<Key, T> > &input, bool sorted_input = false, Predicate const &pred = Predicate()) {

          this->set_balanced(false);
          this->set_globally_sorted(false);

          return this->local_insert(input, sorted_input, pred);
      }


      /**
       * @brief erase elements with the specified keys in the distributed sorted_multimap.  return how much was erased.
       * @note  does not change global or local sort orders.
       * @param first
       * @param last
       *
       * TODO: return size should be global value
       */
      template <typename Predicate = TruePredicate>
      size_t erase(::std::vector<Key>& keys, bool sorted_input = false, Predicate const & pred = Predicate() ) {
        // even if count is 0, still need to participate in mpi calls.  if (keys.size() == 0) return;
          BL_BENCH_INIT(erase);

          if (this->empty()) {
              BL_BENCH_REPORT_MPI_NAMED(erase, "base_sorted_map:erase", this->comm);
        	  return 0;
          }

        size_t before = c.size();

        if (this->comm.size() > 1) {
          // ensure that the container splitters are setup properly, and load balanced.
          BL_BENCH_COLLECTIVE_START(erase, "global_sort", this->comm);
          this->redistribute();
          BL_BENCH_END(erase, "global_sort", this->local_size());


          // remove duplicates
          BL_BENCH_START(erase);
          auto recv_counts(::dsc::distribute(keys, this->key_to_rank, sorted_input, this->comm));
          BLISS_UNUSED(recv_counts);
          BL_BENCH_END(erase, "dist_query", keys.size());

          sorted_input = false;  // keys not sorted across buckets.


        } else {
          BL_BENCH_START(erase);
          this->local_sort();  // sort container before deleting.
          BL_BENCH_END(erase, "local_sort", keys.size());

        }

        BL_BENCH_START(erase);
        // keep unique keys
        ::fsc::sorted_unique(keys, sorted_input, typename Base::TransformedLess(), typename Base::TransformedEqual());
        BL_BENCH_END(erase, "unique_keys", keys.size());


        BL_BENCH_START(erase);
        //== now call local remove.
        // first get the range of intersection
        auto overlap = QueryProcessor<false>::intersect(this->c.begin(), this->c.end(),
            keys.begin(), keys.end(), sorted_input);

        // do the work.  skip duplicates = true
        size_t kept = QueryProcessor<false>::process(overlap.first, overlap.second,
            keys.begin(), keys.end(), overlap.first, erase_element, sorted_input, pred);

        // move the last part.
        auto end = overlap.first;
        ::std::advance(end, kept);   // last valid position within overlap range.
        end = ::std::move(overlap.second, this->c.end(), end);   // move everything after the overlap range
        this->c.erase(end, this->c.end());  // and clear the rest.
        BL_BENCH_END(erase, "erase", keys.size());


        this->set_balanced(false);


        BL_BENCH_REPORT_MPI_NAMED(erase, "base_sorted_map:erase", this->comm);

        return before - c.size();
      }

      template <typename Predicate = TruePredicate>
      size_t erase(Predicate const & pred = Predicate()) {

        if (this->local_empty()) return 0;

        size_t before = c.size();
        if (!::std::is_same<Predicate, TruePredicate>::value) {

          this->local_sort();

          auto end = ::std::partition(c.begin(), c.end(), [pred](value_type const &x) { return !pred(x.first); });
          c.erase(end, c.end());

        } else {
          // identity predicate - all are erased.
          this->local_clear();
        }

        this->set_balanced(false);  // TODO: except when all are erasing with same identity predicate.

        if (this->comm.size() > 1) this->comm.barrier();

        return before - c.size();
      }


      // =============================  overrides.

      // note that for each method, there is a local version of the operartion.
      // this is for use by the asynchronous version of communicator as callback for any messages received.
      /// check if empty.
      virtual bool local_empty() const {
        return c.empty();
      }

      /// get size of local container
      virtual size_t local_size() const {
        return c.size();
      }

      /// get the size of unique keys.
      virtual size_t local_unique_size() const {
        return this->local_size();
      }




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


      struct LocalFind {
          template<bool linear, class DBIter, typename Query, class OutputIter>
          size_t operator()(DBIter &range_begin, DBIter &el_end, DBIter const &range_end, Query const &v, OutputIter &output) const {
              // TODO: LINEAR SEARCH, O(n).  vs BINARY SEARCH, O(mlogn)

              // map, so only 1 entry.
              range_begin = ::fsc::lower_bound<linear>(el_end, range_end, v, Base::less);
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
          template<bool linear, class DBIter, typename Query, class OutputIter, class Predicate = TruePredicate>
          size_t operator()(DBIter &range_begin, DBIter &el_end, DBIter const &range_end,
        		  Query const &v, OutputIter &output, Predicate const &pred) const {
              // TODO: LINEAR SEARCH, O(n).  vs BINARY SEARCH, O(mlogn)

              // map, so only 1 entry.
              range_begin = ::fsc::lower_bound<linear>(el_end, range_end, v, Base::less);
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

      // ============= local reduction override.
      virtual void local_reduction(::std::vector<::std::pair<Key, T> > &input, bool sorted_input = false) {

        ::fsc::sorted_unique(input, sorted_input, typename Base::Base::TransformedLess(), typename Base::Base::TransformedEqual());
      }

      using Base::redistribute;

      // default is for map
      virtual void redistribute() {

        BL_BENCH_INIT(rehash);

        bool balanced = this->is_balanced();
        bool gsorted = this->is_globally_sorted();

        if (balanced && gsorted) // already balanced and globally_sorted.
        {
          // so the splitters are correct as well.

          // ensure locally sorted.
          BL_BENCH_START(rehash);
          this->local_reduction(this->c, this->sorted);
          BL_BENCH_END(rehash, "local_sort", this->c.size());

          // then return.
          BL_BENCH_REPORT_MPI_NAMED(rehash, "sorted_map:rehash", this->comm);
          return;
        }

        //printf("c size before: %lu\n", this->c.size());
        BL_BENCH_START(rehash);
        BL_BENCH_END(rehash, "begin", this->c.size());

        // stop if there are no data to rehash on any of the nodes.
        if (this->empty()) {
          this->key_to_rank.map.clear();
          BL_BENCH_REPORT_MPI_NAMED(rehash, "sorted_map:rehash", this->comm);

          this->sorted = true;
          this->set_balanced(true);
          this->set_globally_sorted(true);

          return;
        }


        if (this->comm.size() > 1) {
          // first balance

          if (!balanced) {
            BL_BENCH_START(rehash);
            this->c = ::mxx::stable_distribute(this->c, this->comm);
            BL_BENCH_END(rehash, "block1", this->c.size());
          }

          // global sort if needed
          BL_BENCH_START(rehash);
          if (!gsorted)
            ::mxx::sort(this->c.begin(), this->c.end(), Base::Base::less, this->comm);
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
                if (boundary_ids[offset - 1] == this->comm.rank()) { // this processor owns the value
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

            // at this point, each kmer is globally unique.  we can just run stable_distribute.

            // and final rebalance
            BL_BENCH_START(rehash);
            // rebalance
            this->c = ::mxx::stable_distribute(this->c, this->comm);
            BL_BENCH_END(rehash, "block2", this->c.size());

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


          BL_BENCH_END(rehash, "splitter2", this->c.size());
          // no need to redistribute - each entry is unique so nothing is going to span processor boundaries.

        } else {
          BL_BENCH_START(rehash);
          // local reduction
          this->local_reduction(this->c, this->sorted);
          BL_BENCH_END(rehash, "reduc", this->c.size());
        }

        this->sorted = true;
        this->set_balanced(true);
        this->set_globally_sorted(true);
        //printf("c size after: %lu\n", this->c.size());


        BL_BENCH_REPORT_MPI_NAMED(rehash, "sorted_map:rehash", this->comm);

      }


    public:

      sorted_map(const mxx::comm& _comm) : Base(_comm) {}

      virtual ~sorted_map() {};

      // specialized here so that the local_find functor can be used.
      template <class Predicate = TruePredicate>
      ::std::vector<::std::pair<Key, T> > find(::std::vector<Key>& keys, bool sorted_input = false,
    		  Predicate const& pred = Predicate()) const {
          return Base::find(find_element, keys, sorted_input, pred);
      }
      template <class Predicate = TruePredicate>
      ::std::vector<::std::pair<Key, T> > find_collective(::std::vector<Key>& keys, bool sorted_input = false,
    		  Predicate const& pred = Predicate()) const {
          return Base::find_a2a(find_element, keys, sorted_input, pred);
      }
      template <class Predicate = TruePredicate>
      ::std::vector<::std::pair<Key, T> > find(Predicate const& pred = Predicate()) const {
          return Base::find(find_element, pred);
      }

      // explicitly get the base class version of insert.
      using Base::insert;
      using Base::erase;
      using Base::count;

      /// update the multiplicity.  only multimap needs to do this.
      virtual float get_multiplicity() const {
        this->redistribute();

        return 1.0f;
      }

      virtual size_t unique_size() const {

        this->redistribute();

        return Base::unique_size();
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

      mutable size_t local_unique_count;

      struct LocalFind {
          template<bool linear, class DBIter, typename Query, class OutputIter>
          size_t operator()(DBIter &range_begin, DBIter &el_end, DBIter const &range_end, Query const &v, OutputIter &output) const {
              // TODO: LINEAR SEARCH, O(n).  vs BINARY SEARCH, O(mlogn)


              range_begin = ::fsc::lower_bound<linear>(el_end, range_end, v, Base::less);

              el_end = ::fsc::upper_bound<linear>(range_begin, range_end, v, Base::less);

              // difference between the 2 iterators is the part that's equal.
              if (range_begin == range_end) return 0;

              // add the output entry, if found.
              output = ::std::copy(range_begin, el_end, output);

              return ::std::distance(range_begin, el_end);
          }
          template<bool linear, class DBIter, typename Query, class OutputIter, class Predicate = TruePredicate>
          size_t operator()(DBIter &range_begin, DBIter &el_end, DBIter const &range_end,
        		  Query const &v, OutputIter &output, Predicate const &pred) const {
              // TODO: LINEAR SEARCH, O(n).  vs BINARY SEARCH, O(mlogn)
        	  //OutputIter output_orig = output;

              range_begin = ::fsc::lower_bound<linear>(el_end, range_end, v, Base::less);

              el_end = ::fsc::upper_bound<linear>(range_begin, range_end, v, Base::less);

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

      using Base::redistribute;

      // default is for multimap
      virtual void redistribute() {

        BL_BENCH_INIT(rehash);


        bool balanced = this->is_balanced();
        bool gsorted = this->is_globally_sorted();

        if (balanced && gsorted) // already balanced and globally_sorted.
        {
          // so the splitters are correct as well.

          // ensure locally sorted.
          BL_BENCH_START(rehash);
          this->local_sort();
          BL_BENCH_END(rehash, "local_sort", this->c.size());

          // then return.
          BL_BENCH_REPORT_MPI_NAMED(rehash, "sorted_multimap:rehash", this->comm);
          return;
        }



        //printf("c size before: %lu\n", this->c.size());
        BL_BENCH_START(rehash);
        BL_BENCH_END(rehash, "begin", this->c.size());

        // stop if there are no data to rehash on any of the nodes.
        if (this->empty()) {
          this->key_to_rank.map.clear();
          BL_BENCH_REPORT_MPI_NAMED(rehash, "sorted_multimap:rehash", this->comm);

          this->sorted = true;
          this->set_balanced(true);
          this->set_globally_sorted(true);

          return;
        }


        // stop if there are no data to rehash on any of the nodes.
        if (this->empty()) {
          this->key_to_rank.map.clear();
          BL_BENCH_REPORT_MPI_NAMED(rehash, "sorted_multimap:rehash", this->comm);
          return;
        }

        if (this->comm.size() > 1) {
          // first balance

          // TODO: stable_block_decompose uses all2all internally.  is it better to move the deltas ourselves?
          if (!balanced) {
            BL_BENCH_START(rehash);
            ::mxx::stable_distribute(this->c, this->comm).swap(this->c);
            BL_BENCH_END(rehash, "block1", this->c.size());
          }

          // sort if needed
          if (!gsorted) {
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

          // note that key_to_rank.map needs to be unique.
          auto map_end = std::unique(this->key_to_rank.map.begin(), this->key_to_rank.map.end(),
                                     Base::Base::equal);
          this->key_to_rank.map.erase(map_end, this->key_to_rank.map.end());
          assert(this->key_to_rank.map.size() > 0);

          BL_BENCH_END(rehash, "splitter1", this->key_to_rank.map.size());


          //        mxx::datatype<::std::pair<Key, T> > dt;
          //        MPI_Datatype mpi_dt = dt.type();
          //        this->key_to_rank.map = ::mxx::sample_block_decomp(d.begin(), d.end(), Base::Base::less, this->comm.size() - 1, this->comm, mpi_dt);

          // redistribute using new pivots.  in trouble if we have a value that takes up a large part of the partition.
          BL_BENCH_START(rehash);
          ::std::vector<size_t> send_counts =
              ::fsc::sorted_bucketing(this->c, this->key_to_rank.map, Base::Base::less, this->comm.size());
          BL_BENCH_END(rehash, "bucket", this->c.size());

//          size_t tt = std::accumulate(send_counts.begin(), send_counts.end(), 0UL);
//          if (tt == 0) printf("rank %d has nothing.\n", this->comm.rank());
//          if (tt == send_counts[this->comm.rank()]) printf("rank %d nothing to send to others\n", this->comm.rank());

          BL_BENCH_COLLECTIVE_START(rehash, "a2a", this->comm);
          // TODO: readjust boundaries using all2all.  is it better to move the deltas ourselves?
          this->c = mxx::all2allv(this->c, send_counts, this->comm);
          BL_BENCH_END(rehash, "a2a", this->c.size());

        } else {

          BL_BENCH_START(rehash);
          this->local_sort();
          BL_BENCH_END(rehash, "local_sort", this->c.size());

        }

        // all redistributed.
        this->sorted = true;
        this->set_balanced(true);
        this->set_globally_sorted(true);

        // update the local unique_count.  only when redistributed.
        local_unique_count = 0;
        if (this->c.size() > 1) {
          // now count unique
          auto start = this->c.begin();
          while ((start = ::std::adjacent_find(start, this->c.end(), Base::Base::less)) != this->c.end()) {
            ++start;
            ++local_unique_count;
          }
          ++local_unique_count; // for the last one?
        } else {
          local_unique_count = this->c.size();
        }


        BL_BENCH_REPORT_MPI_NAMED(rehash, "sorted_multimap:rehash", this->comm);

      }

    public:


      sorted_multimap(const mxx::comm& _comm) : Base(_comm), local_unique_count(0) {}

      virtual ~sorted_multimap() {}


      //========= override find, so can pass in LocalFind object.
      template <class Predicate = TruePredicate>
      ::std::vector<::std::pair<Key, T> > find(::std::vector<Key>& keys, bool sorted_input = false,
    		  Predicate const& pred = Predicate()) const {
          return Base::find(find_element, keys, sorted_input, pred);
      }
      template <class Predicate = TruePredicate>
      ::std::vector<::std::pair<Key, T> > find_collective(::std::vector<Key>& keys, bool sorted_input = false,
    		  Predicate const& pred = Predicate()) const {
          return Base::find_a2a(find_element, keys, sorted_input, pred);
      }
      template <class Predicate = TruePredicate>
      ::std::vector<::std::pair<Key, T> > find(Predicate const& pred = Predicate()) const {
          return Base::find(find_element, pred);
      }

      // explicitly get the base class version of insert.
      using Base::insert;
      using Base::erase;
      using Base::count;

      /// update the multiplicity.  only multimap needs to do this.
      virtual float get_multiplicity() const {
          BL_BENCH_INIT(multiplicity);

        // one approach is to add up the number of repeats for the key of each entry, then divide by total count.
        //  sum(count per key) / c.size.
        // problem with this approach is that for unordered map, to get the count for a key is essentially O(count), so we get quadratic time.
        // The approach is VERY SLOW for large repeat count.  - (0.0078125 human: 52 sec, synth: FOREVER.)

        // a second approach is to count the number of unique key then divide the map size by that.
        //  c.size / #unique.  requires unique set
        // To find unique set, we take each bucket, copy to vector, sort it, and then count unique.
        // This is precise, and is faster than the approach above.  (0.0078125 human: 54 sec.  synth: 57sec.)
        // but the n log(n) sort still grows with the duplicate count
        // compute new multiplicity
        BL_BENCH_START(multiplicity);

        size_t n_unique = this->unique_size();
        float multiplicity = 1.0f;
        if (n_unique > 0) {
          // local unique
           multiplicity =
              static_cast<float>(this->size()) /
              static_cast<float>(n_unique);
        }
        BL_BENCH_END(multiplicity, "multiplicity", this->c.size());

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

        BL_BENCH_REPORT_MPI_NAMED(multiplicity, "sorted_multimap:multiplicity", this->comm);

        return multiplicity;
      }

      /**
       * @brief global unique size
       * @note: calculated differently than local_unique size -
       *        requires a redistribute first, which also allows caching of value if distribution is not changed.
       */
      virtual size_t unique_size() const {
        // side effect of redistribute is unique calc.
        this->redistribute();

        if (this->comm.size() == 1)
          return this->local_unique_count;
        else
          return ::mxx::allreduce(this->local_unique_count, this->comm);
      }

      /// get the size of unique keys in the current local container.
      virtual size_t local_unique_size() const {

        // first make sure that the global vector is sorted.
        this->local_sort();

        // now count unique
        size_t count = 0;
        if (this->c.size() > 1) {
          auto start = this->c.begin();
          while ((start = ::std::adjacent_find(start, this->c.end(), Base::less)) != this->c.end()) {
            ++start;
            ++count;
          }
          ++count; // for the last one?
        } else {
          count = this->c.size();
        }
        return count;
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


      virtual void local_reduction(::std::vector<::std::pair<Key, T> >& input, bool sorted_input = false) {
        if (input.size() == 0) return;

        ::fsc::sort(input, sorted_input, Base::Base::Base::less);

        auto b = input.begin();
        auto e = input.end();

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

      // explicitly get the base class version of insert.
      using Base::insert;
      using Base::erase;
      using Base::count;
      using Base::find;
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

      using key_type              = Key;
      using mapped_type           = T;
      using value_type            = ::std::pair<Key, T>;
      using iterator              = typename local_container_type::iterator;
      using const_iterator        = typename local_container_type::const_iterator;
      using size_type             = typename local_container_type::size_type;
      using difference_type       = typename local_container_type::difference_type;


      counting_sorted_map(const mxx::comm& _comm) : Base(_comm) {}

      virtual ~counting_sorted_map() {};

      // explicitly get the base class version of insert.
      using Base::insert;
      using Base::erase;
      using Base::count;
      using Base::find;

      /**
       * @brief insert new elements in the distributed sorted_multimap.  convert from Key to Key-count pair
       * @param first
       * @param last
       */
      template <class Predicate = TruePredicate>
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
        BL_BENCH_REPORT_MPI_NAMED(count_insert, "count_sorted_map:insert", this->comm);
        return count;
      }



  };



} /* namespace dsc */


#endif // BLISS_DISTRIBUTED_SORTED_MAP_HPP
