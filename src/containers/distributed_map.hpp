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
 * @file    distributed_map.hpp
 * @ingroup index
 * @author  Tony Pan <tpan7@gatech.edu>
 * @brief   Implements the distributed_multimap, distributed map, and distributed_reduction_map
 *          data structures.
 *
 *          implementation is hash-base (O(1) lookup). later will support sort-based (load balanced).
 *
 *          PERFORMANCE OF THIS CLASS IS VERY BAD BECAUSE OF THE TREE INSERT/SEARCH.  USE AT YOUR OWN PERIL.
 *
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

#ifndef BLISS_DISTRIBUTED_MAP_HPP
#define BLISS_DISTRIBUTED_MAP_HPP


#include <map>  // local storage hash table  // for multimap
#include <unordered_set>  // local storage hash table  // for multimap
#include <utility> 			  // for std::pair

//#include <sparsehash/dense_hash_map>  // not a multimap, where we need it most.
#include <functional> 		// for std::function and std::hash
#include <algorithm> 		// for sort, stable_sort, unique, is_sorted
#include <iterator>  // advance, distance

#include <cstdint>  // for uint8, etc.

#include <type_traits>

#include <mxx/collective.hpp>
#include <mxx/reduction.hpp>
#include "utils/benchmark_utils.hpp"  // for timing.
#include "utils/logging.h"
#include "containers/distributed_map_base.hpp"


namespace dsc  // distributed std container
{


  /**
   * @brief  distributed unordered map following std unordered map's interface.
   * @details   This class is modeled after the std::map.
   *         it has as much of the same methods of std::map as possible.  however, all methods consider the fact
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
   *  any operation that uses sort is not going to scale well.  this includes "sort_unique_key, sort_unique_tuple, local_reduction"...
   *  to reduce this, we can try by using a hash set instead.  http://www.vldb.org/pvldb/2/vldb09-257.pdf, http://www.vldb.org/pvldb/vol7/p85-balkesen.pdf
   *
   *  conditional version of insert/erase/find/count supports predicate that operate on INTERMEDIATE RESULTS.  input (key,key-value pair) can be pre-filtered.
   *    output (query result, e.g.) can be post filtered (and optionally reduce comm volume).
   *    intermediate results (such as counting in multimap only if T has certain value) can only be filtered at during local_operation.
   *
   *
   * key to proc assignment can be done as hash or splitters in sorted range.
   * tuples can be sotred in hash table or as sorted array.
   *   hash-hash combination works
   *  sort-sort combination works as well
   *  hash-sort combination can work.  advantage is in range query.
   *  sort-hash combination would be expensive for updating splitters
   *
   * This version is the hash-hash.
   *
   * @tparam Key
   * @tparam T
   * @tparam Container  default to map and unordered multimap, requiring 5 template params.
   * @tparam Comm   default to mpi_collective_communicator       communicator for global communication. may hash or sort.
   * @tparam KeyTransform   transform function for the key.  can supply identity.  requires a single template argument (Key).  useful for mapping kmolecule to kmer.
   * @tparam Hash   hash function for local and distribution.  requires a template arugment (Key), and a bool (prefix, chooses the MSBs of hash instead of LSBs)
   * @tparam Equal   default to ::std::equal_to<Key>   equal function for the local storage.
   * @tparam Alloc  default to ::std::allocator< ::std::pair<const Key, T> >    allocator for local storage.
   */
  template<typename Key, typename T,
    template <typename, typename, typename, typename> class Container,
    class Comm,
    template <typename> class KeyTransform,
    template <typename, bool> class Hash,
    class Less = ::std::less<Key>,
    class Alloc = ::std::allocator< ::std::pair<const Key, T> >
  > class distributed_map_base : public ::dsc::map_base<Key, T, Comm, KeyTransform, Less, ::std::equal_to<Key>, Alloc> {

    protected:
      using Base = ::dsc::map_base<Key, T, Comm, KeyTransform, Less, ::std::equal_to<Key>, Alloc>;


      struct KeyToRank {
          Hash<Key, true> proc_hash;
          const int p;

          // 2x comm size to allow more even distribution?
          KeyToRank(int comm_size) : proc_hash(ceilLog2(comm_size)), p(comm_size) {};

          inline int operator()(Key const & x) const {
            //            printf("KeyToRank operator. commsize %d  key.  hashed to %d, mapped to proc %d \n", p, proc_hash(Base::trans(x)), proc_hash(Base::trans(x)) % p);
            return proc_hash(Base::trans(x)) % p;
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
      struct QueryProcessor {  // assume unique, always.

          // assumes that container is sorted. and exact overlap region is provided.  do not filter output here since it's an output iterator.
          template <class DB, class QueryIter, class OutputIter, class Operator, class Predicate = TruePredicate>
          static size_t process(DB &db,
                                QueryIter query_begin, QueryIter query_end,
                                OutputIter &output, Operator & op,
                                bool sorted_query = false, Predicate const &pred = Predicate()) {

              if (query_begin == query_end) return 0;

              typename ::std::iterator_traits<QueryIter>::value_type v;
              size_t count = 0;  // before size.
              for (auto it = query_begin; it != query_end; ++it) {
                v = *it;
                if (!::std::is_same<Predicate, TruePredicate>::value)
                  count += op(db, v, output, pred);
                else
                  count += op(db, v, output);
              }
              return count;
          }

      };




    public:
      using local_container_type = Container<Key, T, typename Base::TransformedLess, Alloc>;

      // std::multimap public members.
      using key_type              = typename local_container_type::key_type;
      using mapped_type           = typename local_container_type::mapped_type;
      using value_type            = typename local_container_type::value_type;
      using iterator              = typename local_container_type::iterator;
      using const_iterator        = typename local_container_type::const_iterator;
      using size_type             = typename local_container_type::size_type;
      using difference_type       = typename local_container_type::difference_type;

    protected:
      local_container_type c;

      struct LocalCount {
          // unfiltered.
          template<class DB, typename Query, class OutputIter>
          size_t operator()(DB &db, Query const &v, OutputIter &output) const {
              *output = ::std::move(::std::make_pair(v, db.count(v)));
              ++output;
              return 1;
          }
          // filtered element-wise.
          template<class DB, typename Query, class OutputIter, class Predicate = TruePredicate>
          size_t operator()(DB &db, Query const &v, OutputIter &output,
                            Predicate const& pred) const {
              auto range = db.equal_range(v);

              // add the output entry.
              size_t count = 0;
              if (pred(range.first, range.second))  // operator to decide if range matches.
                count = ::std::count_if(range.first, range.second, pred);  // operator for each element in range.

              *output = ::std::move(::std::make_pair(v, count));
              ++output;
              return 1;
          }
          // no filter by range AND elemenet for now.
      } count_element;

      struct LocalErase {
          /// Return how much was KEPT.
          template<class DB, typename Query, class OutputIter>
          size_t operator()(DB &db, Query const &v, OutputIter &output) {
              size_t before = db.size();
              db.erase(v);
              return before - db.size();
          }
          /// Return how much was KEPT.
          template<class DB, typename Query, class OutputIter, class Predicate = TruePredicate>
          size_t operator()(DB &db, Query const &v, OutputIter &,
                            Predicate const & pred) {
              auto range = db.equal_range(v);

              // check range first.  then erase.  only removed iterators are invalidated.
              // order of remaining elements preserved.  true for map/multimap, unordered or not.
              size_t count = 0;
              auto tmp = range.first;
              if (pred(range.first, range.second)) { // operator to decide if range matches.
                for (auto it = range.first; it != range.second;) {
                  if (pred(*it)) {
                    // advance, then remove.
                    tmp = it;  ++it;
                    // remove.
                    db.erase(tmp);  // erase entry at 1 iterator position. (not by value).
                    ++count;
                  } else {
                    // keep.  so just advance
                    ++it;
                  }
                }
              }

              return count;
          }
          // no filter by range AND elemenet for now.
      } erase_element;

      /**
       * @brief insert new elements in the distributed multimap.
       * @param first
       * @param last
       */
      template <class InputIterator>
      size_t local_insert(InputIterator first, InputIterator last) {
          if (first == last) return 0;

          size_t before = c.size();
          c.insert(first, last);
          return c.size() - before;
      }

      /**
       * @brief insert new elements in the distributed multimap.  example use: stop inserting if more than x entries.
       * @param first
       * @param last
       */
      template <class InputIterator, class Predicate>
      size_t local_insert(InputIterator first, InputIterator last, Predicate const &pred) {
          if (first == last) return 0;

          size_t before = c.size();
          for (auto it = first; it != last; ++it) {
            if (pred(*it)) c.emplace(::std::move(*it));
          }
          return c.size() - before;
      }


      /// clears the map
      virtual void local_clear() noexcept {
          c.clear();
      }


      /**
       * @brief find elements with the specified keys in the distributed unordered_multimap.
       * @param first
       * @param last
       */
      template <class LocalFind, typename Predicate = TruePredicate>
      ::std::vector<::std::pair<Key, T> > find_a2a(LocalFind & find_element, ::std::vector<Key>& keys, bool sorted_input = false, Predicate const& pred = Predicate()) const {
        BL_BENCH_INIT(find);

        BL_BENCH_START(find);
        ::std::vector<::std::pair<Key, T> > results;
        ::fsc::back_emplace_iterator<::std::vector<::std::pair<Key, T> > > emplace_iter(results);
        // even if count is 0, still need to participate in mpi calls.  if (keys.size() == 0) return results;
        BL_BENCH_END(find, "begin", keys.size());

        if (this->comm.size() > 1) {
            BL_BENCH_COLLECTIVE_START(find, "dist_query", this->comm);
            // distribute (communication part)
            std::vector<size_t> recv_counts(::dsc::distribute_unique(keys, this->key_to_rank, sorted_input, this->comm,
                                                             typename Base::TransformedFarmHash(), typename Base::TransformedEqual()));

            BL_BENCH_END(find, "dist_query", keys.size());

          // local find. memory utilization a potential problem.
          // do for each src proc one at a time.

          BL_BENCH_START(find);
          results.reserve(keys.size() * this->key_multiplicity);                   // TODO:  should estimate coverage.
          //printf("reserving %lu\n", keys.size() * this->key_multiplicity);
          BL_BENCH_END(find, "reserve", results.capacity());

          BL_BENCH_START(find);
          std::vector<size_t> send_counts(this->comm.size(), 0);

          auto start = keys.begin();
          auto end = start;
          for (int i = 0; i < this->comm.size(); ++i) {
            ::std::advance(end, recv_counts[i]);

            // work on query from process i.
            send_counts[i] = QueryProcessor::process(c, start, end, emplace_iter, find_element, sorted_input, pred);
           // if (this->comm.rank() == 0) BL_DEBUGF("R %d added %d results for %d queries for process %d\n", this->comm.rank(), send_counts[i], recv_counts[i], i);

            start = end;
          }
          BL_BENCH_END(find, "local_find", results.size());

          BL_BENCH_COLLECTIVE_START(find, "a2a2", this->comm);
          // send back using the constructed recv count
          results = mxx::all2allv(results, send_counts, this->comm);
          BL_BENCH_END(find, "a2a2", results.size());

        } else {

          BL_BENCH_START(find);
          // keep unique keys
          ::fsc::unique(keys, sorted_input,
                            typename Base::TransformedFarmHash(), typename Base::TransformedEqual());
          BL_BENCH_END(find, "uniq1", keys.size());

          BL_BENCH_START(find);
          results.reserve(keys.size() * this->key_multiplicity);                   // TODO:  should estimate coverage.
          //printf("reserving %lu\n", keys.size() * this->key_multiplicity);
          BL_BENCH_END(find, "reserve", results.capacity() );

          BL_BENCH_START(find);
          QueryProcessor::process(c, keys.begin(), keys.end(), emplace_iter, find_element, sorted_input, pred);
          BL_BENCH_END(find, "local_find", results.size());
        }

        BL_BENCH_REPORT_MPI_NAMED(find, "base_map:find_a2a", this->comm);

        return results;

      }



      /**
       * @brief find elements with the specified keys in the distributed unordered_multimap.
       *
       * why this version that uses isend and irecv?  because all2all version requires all result data to be in memory.
       * this one can do it one source process at a time.
       *
       * @param first
       * @param last
       */
      template <class LocalFind, typename Predicate = TruePredicate>
      ::std::vector<::std::pair<Key, T> > find(LocalFind & find_element, ::std::vector<Key>& keys, bool sorted_input = false, Predicate const& pred = Predicate()) const {
        BL_BENCH_INIT(find);

        BL_BENCH_START(find);
        ::std::vector<::std::pair<Key, T> > results;
        ::fsc::back_emplace_iterator<::std::vector<::std::pair<Key, T> > > emplace_iter(results);
        // even if count is 0, still need to participate in mpi calls.  if (keys.size() == 0) return results;

        ::std::vector<::std::pair<Key, T> > local_results;
        ::fsc::back_emplace_iterator<::std::vector<::std::pair<Key, T> > > local_emplace_iter(local_results);
        BL_BENCH_END(find, "begin", keys.size());

        if (this->comm.size() > 1) {
            BL_BENCH_COLLECTIVE_START(find, "dist_query", this->comm);
            // distribute (communication part)
            std::vector<size_t> recv_counts(::dsc::distribute_unique(keys, this->key_to_rank, sorted_input, this->comm,
                                                                     typename Base::TransformedFarmHash(), typename Base::TransformedEqual()));
            BL_BENCH_END(find, "dist_query", keys.size());


          // local count to determine amount of memory to allocate at destination.
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
            QueryProcessor::process(c, start, end, count_emplace_iter, count_element, sorted_input, pred);
            send_counts[i] =
                ::std::accumulate(count_results.begin(), count_results.end(), static_cast<size_t>(0),
                                  [](size_t v, ::std::pair<Key, size_t> const & x) {
              return v + x.second;
            });
//            for (auto it = count_results.begin(), max = count_results.end(); it != max; ++it) {
//              send_counts[i] += it->second;
//            }
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
          //printf("reserving %lu\n", keys.size() * this->key_multiplicity);
          BL_BENCH_END(find, "reserve", resp_total);

          BL_BENCH_START(find);
          auto recv_displs = mxx::impl::get_displacements(recv_counts);  // compute response displacements.
          int recv_from, send_to;
          size_t found;
          total = 0;
          std::vector<MPI_Request> reqs(2 * this->comm.size());

          mxx::datatype dt = mxx::get_datatype<::std::pair<Key, T> >();
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
            found = QueryProcessor::process(c, start, end, local_emplace_iter, find_element, sorted_input, pred);
           // if (this->comm.rank() == 0) BL_DEBUGF("R %d added %d results for %d queries for process %d\n", this->comm.rank(), send_counts[i], recv_counts[i], i);
            total += found;
            //== now send the results immediately - minimizing data usage so we need to wait for both send and recv to complete right now.

            MPI_Isend(&(local_results[0]), found, dt.type(), send_to,
                      i, this->comm, &reqs[2 * i + 1]);

            // wait for both requests to complete.
            MPI_Waitall(2, &reqs[2 * i], MPI_STATUSES_IGNORE);


            // verify correct? done by comparing to previous code.


            //printf("Rank %d local find send to %d:  query %d result sent %d (%d).  recv from %d received %d\n", this->comm.rank(), send_to, recv_counts[send_to], found, send_counts[send_to], recv_from, resp_counts[recv_from]);
          }
          //printf("Rank %d total find %lu\n", this->comm.rank(), total);
          BL_BENCH_END(find, "find_send", results.size());

        } else {

          BL_BENCH_START(find);
          // keep unique keys
          ::fsc::unique(keys, sorted_input,
                            typename Base::TransformedFarmHash(), typename Base::TransformedEqual());
          BL_BENCH_END(find, "uniq1", keys.size());


          BL_BENCH_START(find);
          ::std::vector<::std::pair<Key, size_t> > count_results;
          count_results.reserve(keys.size());
          ::fsc::back_emplace_iterator<::std::vector<::std::pair<Key, size_t> > > count_emplace_iter(count_results);

          // count now.
          QueryProcessor::process(c, keys.begin(), keys.end(), count_emplace_iter, count_element, sorted_input, pred);
          size_t count = ::std::accumulate(count_results.begin(), count_results.end(), static_cast<size_t>(0),
                                            [](size_t v, ::std::pair<Key, size_t> const & x) {
                        return v + x.second;
                      });
//          for (auto it = count_results.begin(), max = count_results.end(); it != max; ++it) {
//            count += it->second;
//          }
          BL_BENCH_END(find, "local_count", count);

          BL_BENCH_START(find);
          results.reserve(count);                   // TODO:  should estimate coverage.
          //printf("reserving %lu\n", keys.size() * this->key_multiplicity);
          BL_BENCH_END(find, "reserve", count);

          BL_BENCH_START(find);
          QueryProcessor::process(c, keys.begin(), keys.end(), emplace_iter, find_element, sorted_input, pred);
          BL_BENCH_END(find, "local_find", results.size());
        }

        BL_BENCH_REPORT_MPI_NAMED(find, "base_map:find_isend", this->comm);

        return results;

      }


      template <class LocalFind, typename Predicate = TruePredicate>
      ::std::vector<::std::pair<Key, T> > find(LocalFind & find_element, Predicate const& pred = Predicate()) const {
        ::std::vector<::std::pair<Key, T> > results;
        ::fsc::back_emplace_iterator<::std::vector<::std::pair<Key, T> > > emplace_iter(results);

        auto keys = this->keys();
        results.reserve(keys.size() * this->key_multiplicity);                   // TODO:  should estimate coverage.

        QueryProcessor::process(c, keys.begin(), keys.end(), emplace_iter, find_element, true, pred);

        if (this->comm.size() > 1) this->comm.barrier();

        return results;
      }


      distributed_map_base(::mxx::comm const & _comm) : Base(_comm),
          key_to_rank(_comm.size()) {}


    public:

      virtual ~distributed_map_base() {};

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
      virtual size_t update_multiplicity() { return this->key_multiplicity; }


      /// convert the map to a vector
      virtual void to_vector(std::vector<std::pair<Key, T> > & result) const {
        result.clear();
        if (c.empty()) return;
        result.reserve(c.size());
        ::fsc::back_emplace_iterator<::std::vector<::std::pair<Key, T> > > emplace_iter(result);
        ::std::copy(c.begin(), c.end(), emplace_iter);
      }
      /// extract the unique keys of a map.
      virtual void keys(std::vector<Key> & result) const {
        result.clear();
        if (c.empty()) return;

        typename Base::template UniqueKeySetUtilityType<Key> temp(c.size());
        auto end = c.end();
        for (auto it = c.begin(); it != end; ++it) {
          temp.emplace((*it).first);
        }
        result.assign(temp.begin(), temp.end());
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
      void reserve( size_t n) {}


      /// rehash the local container.  n is the local container size.  this allows different processes to individually adjust its own size.
      void rehash( size_type n) {}


      /**
       * @brief count elements with the specified keys in the distributed unordered_multimap.
       * @param first
       * @param last
       */
      template <class Predicate = TruePredicate>
      ::std::vector<::std::pair<Key, size_type> > count(::std::vector<Key>& keys, bool sorted_input = false,
                                                        Predicate const& pred = Predicate() ) const {
        BL_BENCH_INIT(count);

        BL_BENCH_START(count);
        ::std::vector<::std::pair<Key, size_type> > results;
        ::fsc::back_emplace_iterator<::std::vector<::std::pair<Key, size_type> > > emplace_iter(results);
        // even if count is 0, still need to participate in mpi calls.  if (keys.size() == 0) return results;
        BL_BENCH_END(count, "begin", keys.size());


        if (this->comm.size() > 1) {

            BL_BENCH_START(count);
            std::vector<size_t> recv_counts(::dsc::distribute_unique(keys, this->key_to_rank, sorted_input, this->comm,
                                                                     typename Base::TransformedFarmHash(), typename Base::TransformedEqual()));
            BL_BENCH_END(count, "dist_query", keys.size());


          // local count. memory utilization a potential problem.
          // do for each src proc one at a time.
          BL_BENCH_START(count);
          results.reserve(keys.size());                   // TODO:  should estimate coverage.
          BL_BENCH_END(count, "reserve", results.capacity());

          BL_BENCH_START(count);
          auto start = keys.begin();
          auto end = start;
          for (int i = 0; i < this->comm.size(); ++i) {
            ::std::advance(end, recv_counts[i]);

            // within start-end, values are unique, so don't need to set unique to true.
            QueryProcessor::process(c, start, end, emplace_iter, count_element, sorted_input, pred);

            if (this->comm.rank() == 0)
            	BL_DEBUGF("R %d added %lu results for %lu queries for process %d\n", this->comm.rank(), send_counts[i], recv_counts[i], i);

            start = end;
          }
          BL_BENCH_END(count, "local_count", results.size());

          // send back using the constructed recv count
          BL_BENCH_COLLECTIVE_START(count, "a2a2", this->comm);
          results = mxx::all2allv(results, recv_counts, this->comm);
          BL_BENCH_END(count, "a2a2", results.size());
        } else {

          BL_BENCH_START(count);
          // keep unique keys
          ::fsc::unique(keys, sorted_input,
                            typename Base::TransformedFarmHash(), typename Base::TransformedEqual());
          BL_BENCH_END(count, "uniq1", keys.size());



          BL_BENCH_START(count);
          results.reserve(keys.size());                   // TODO:  should estimate coverage.
          BL_BENCH_END(count, "reserve", results.capacity());


          BL_BENCH_START(count);
          // within start-end, values are unique, so don't need to set unique to true.
          QueryProcessor::process(c, keys.begin(), keys.end(), emplace_iter, count_element, sorted_input, pred);
          BL_BENCH_END(count, "local_count", results.size());
        }

        BL_BENCH_REPORT_MPI_NAMED(count, "base_map:count", this->comm);

        return results;

      }



      template <typename Predicate = TruePredicate>
      ::std::vector<::std::pair<Key, size_type> > count(Predicate const & pred = Predicate()) const {
        ::std::vector<::std::pair<Key, size_type> > results;
        ::fsc::back_emplace_iterator<::std::vector<::std::pair<Key, size_t> > > emplace_iter(results);

        auto keys = this->keys();
        results.reserve(keys.size());

        QueryProcessor::process(c, keys.begin(), keys.end(), emplace_iter, count_element, true, pred);

        if (this->comm.size() > 1) this->comm.barrier();
        return results;
      }



      /**
       * @brief erase elements with the specified keys in the distributed unordered_multimap.
       * @param first
       * @param last
       */
      template <class Predicate = TruePredicate>
      size_t erase(::std::vector<Key>& keys, bool sorted_input = false, Predicate const& pred = Predicate() ) {
        // even if count is 0, still need to participate in mpi calls.  if (keys.size() == 0) return;
          size_t before = this->c.size();
          BL_BENCH_INIT(erase);

        if (this->comm.size() > 1) {

            BL_BENCH_START(erase);
            auto recv_counts(::dsc::distribute_unique(keys, this->key_to_rank, sorted_input, this->comm,
                                                      typename Base::TransformedFarmHash(), typename Base::TransformedEqual()));
            BLISS_UNUSED(recv_counts);
            BL_BENCH_END(erase, "dist_query", keys.size());

        }

        // this uses sort order to remove duplicates.  map handled sorted access better,
        // and log access is slow, so get unique first.
		BL_BENCH_START(erase);
		// remove duplicates
    ::fsc::unique(keys, sorted_input,
                      typename Base::TransformedFarmHash(), typename Base::TransformedEqual());
		BL_BENCH_END(erase, "unique", keys.size());

        BL_BENCH_START(erase);
        // then call local remove.
        auto dummy_iter = keys.end();  // process requires a reference.
        QueryProcessor::process(c, keys.begin(), keys.end(), dummy_iter, erase_element, sorted_input, pred);
        BL_BENCH_END(erase, "erase", keys.size());

        BL_BENCH_REPORT_MPI_NAMED(erase, "base_map:erase", this->comm);



        return before - this->c.size();
      }


      template <typename Predicate = TruePredicate>
      size_t erase(Predicate const & pred = Predicate()) {

    	  size_t count = 0;

          if (!::std::is_same<Predicate, TruePredicate>::value) {


        auto keys = this->keys();

        auto dummy_iter = keys.end();  // process requires a reference.

        // keys are ordered because this is based on a map.
        count = QueryProcessor::process(c, keys.begin(), keys.end(), dummy_iter, erase_element, true, pred);

          } else {  // identity matches all, so all are cleared.
        	  count = this->local_size();

        	  this->local_clear();
          }
        if (this->comm.size() > 1) this->comm.barrier();

        return count;
      }

      // update done via erase/insert.

  };


  /**
   * @brief  distributed unordered map following std unordered map's interface.
   * @details   This class is modeled after the std::map.
   *         it has as much of the same methods of std::map as possible.  however, all methods consider the fact
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
   * @tparam Alloc  default to ::std::allocator< ::std::pair<const Key, T> >    allocator for local storage.
   */
  template<typename Key, typename T,
  class Comm,
  template <typename> class KeyTransform,
  template <typename, bool> class Hash,
  class Less = ::std::less<Key>,
  class Alloc = ::std::allocator< ::std::pair<const Key, T> >
  >
  class map : public distributed_map_base<Key, T, ::std::map, Comm, KeyTransform, Hash, Less, Alloc> {
    protected:
      using Base = distributed_map_base<Key, T, ::std::map, Comm, KeyTransform, Hash, Less, Alloc>;


    public:
      using local_container_type = typename Base::local_container_type;

      // std::multimap public members.
      using key_type              = typename local_container_type::key_type;
      using mapped_type           = typename local_container_type::mapped_type;
      using value_type            = typename local_container_type::value_type;
      using iterator              = typename local_container_type::iterator;
      using const_iterator        = typename local_container_type::const_iterator;
      using size_type             = typename local_container_type::size_type;
      using difference_type       = typename local_container_type::difference_type;

    protected:

      struct LocalFind {
          // unfiltered.
          template<class DB, typename Query, class OutputIter>
          size_t operator()(DB &db, Query const &v, OutputIter &output) const {
              auto iter = db.find(v);

              if (iter != db.end()) {
                *output = *iter;
                ++output;
                return 1;
              }  // no insert if can't find it.
              return 0;
          }
          // filtered element-wise.
          template<class DB, typename Query, class OutputIter, class Predicate = TruePredicate>
          size_t operator()(DB &db, Query const &v, OutputIter &output,
                            Predicate const& pred) const {
              auto iter = db.find(v);

              // add the output entry.
              auto next = iter;  ++next;
              if (iter != db.end()) {
                if (pred(iter, next) && pred(*iter)) {
                  *output = *iter;
                  ++output;
                  return 1;
                }
              }
              return 0;
          }
          // no filter by range AND elemenet for now.
      } find_element;


      virtual void local_reduction(::std::vector<::std::pair<Key, T> > &input, bool & sorted_input) {
        ::fsc::sorted_unique(input, sorted_input, typename Base::TransformedLess(), typename Base::TransformedEqual());
      }


    public:

      map(::mxx::comm const & _comm) : Base(_comm) {}

      virtual ~map() {};

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


      /**
       * @brief insert new elements in the distributed unordered_multimap.
       * @param first
       * @param last
       */
      template <typename Predicate = TruePredicate>
      size_t insert(std::vector<::std::pair<Key, T> >& input, bool sorted_input = false, Predicate const & pred = Predicate()) {
        // even if count is 0, still need to participate in mpi calls.  if (input.size() == 0) return;
        BL_BENCH_INIT(insert);

        BL_BENCH_START(insert);
        BL_BENCH_END(insert, "start", input.size());


        // communication part
        if (this->comm.size() > 1) {

            BL_BENCH_START(insert);
            // get mapping to proc
            // keep unique only may not be needed - comm speed may be faster than we can compute unique.
            auto recv_counts(::dsc::distribute_unique(input, this->key_to_rank, sorted_input, this->comm,
                                                             typename Base::TransformedFarmHash(), typename Base::TransformedEqual()));

            BLISS_UNUSED(recv_counts);
            BL_BENCH_END(insert, "dist_data", input.size());
        }


        BL_BENCH_START(insert);
        // local compute part.  called by the communicator.
        size_t count = 0;
        if (!::std::is_same<Predicate, TruePredicate>::value)
          count = this->Base::local_insert(input.begin(), input.end(), pred);
        else
          count = this->Base::local_insert(input.begin(), input.end());
        BL_BENCH_END(insert, "insert", this->c.size());

        BL_BENCH_REPORT_MPI_NAMED(insert, "hashmap:insert", this->comm);

        return count;
      }


  };



  /**
   * @brief  distributed unordered multimap following std unordered multimap's interface.
   * @details   This class is modeled after the std::multimap.
   *         it does not have all the methods of std::multimap.  Whatever methods that are present considers the fact
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
   * @tparam Alloc  default to ::std::allocator< ::std::pair<const Key, T> >    allocator for local storage.
   */
  template<typename Key, typename T,
  class Comm,
  template <typename> class KeyTransform,
  template <typename, bool> class Hash,
  class Less = ::std::less<Key>,
  class Alloc = ::std::allocator< ::std::pair<const Key, T> >
  >
  class multimap : public distributed_map_base<Key, T, ::std::multimap, Comm, KeyTransform, Hash, Less, Alloc> {
    protected:
      using Base = distributed_map_base<Key, T, ::std::multimap, Comm, KeyTransform, Hash, Less, Alloc>;


    public:
      using local_container_type = typename Base::local_container_type;

      // std::multimap public members.
      using key_type              = typename local_container_type::key_type;
      using mapped_type           = typename local_container_type::mapped_type;
      using value_type            = typename local_container_type::value_type;
      using iterator              = typename local_container_type::iterator;
      using const_iterator        = typename local_container_type::const_iterator;
      using size_type             = typename local_container_type::size_type;
      using difference_type       = typename local_container_type::difference_type;

    protected:

      struct LocalFind {
          // unfiltered.
          template<class DB, typename Query, class OutputIter>
          size_t operator()(DB &db, Query const &v, OutputIter &output) const {
              auto range = db.equal_range(v);

              // range's iterators are not random access iterators, so insert calling distance uses ++, slowing down the process.
              // manually insert improves performance here.
              size_t count = 0;
              for (auto it2 = range.first; it2 != range.second; ++it2) {
                *output = *it2;
                ++output;
                ++count;
              }

              return count;
          }
          // filtered element-wise.
          template<class DB, typename Query, class OutputIter, class Predicate = TruePredicate>
          size_t operator()(DB &db, Query const &v, OutputIter &output,
                            Predicate const& pred) const {
              auto range = db.equal_range(v);

              // add the output entry.
              size_t count = 0;
              if (pred(range.first, range.second)) {
                for (auto it2 = range.first; it2 != range.second; ++it2) {
                  if (pred(*it2)) {
                    *output = *it2;
                    ++output;
                    ++count;
                  }
                }
              }
              return count;
          }
          // no filter by range AND elemenet for now.
      } find_element;



    public:


      multimap(::mxx::comm const & _comm) : Base(_comm) { }

      virtual ~multimap() {}

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
        size_t uniq_count = 0;

        //        ::std::vector< ::std::pair<Key, T> > temp;
        //        KeyTransform<Key> trans;
        //        for (int i = 0, max = this->c.bucket_count(); i < max; ++i) {
        //          if (this->c.bucket_size(i) == 0) continue;  // empty bucket. move on.
        //
        //          // copy and sort.
        //          temp.assign(this->c.begin(i), this->c.end(i));  // copy the bucket
        //          // sort the bucket
        //          ::std::sort(temp.begin(), temp.end(), [&] ( ::std::pair<Key, T> const & x,  ::std::pair<Key, T> const & y){
        //            return trans(x.first) < trans(y.first);
        //          });
        // //          auto end = ::std::unique(temp.begin(), temp.end(), this->key_equal_op);
        // //          uniq_count += ::std::distance(temp.begin(), end);
        //
        //          // count via linear scan..
        //          auto x = temp.begin();
        //          ++uniq_count;  // first entry.
        //          // compare pairwise.
        //          auto y = temp.begin();  ++y;
        //          while (y != temp.end()) {
        //            if (trans(x->first) != trans(y->first)) {
        //              ++uniq_count;
        //              x = y;
        //            }
        //            ++y;
        //          }
        //        }
        //        printf("%lu elements, %lu buckets, %lu unique\n", this->c.size(), this->c.bucket_count(), uniq_count);
        // alternative approach to get number of unique keys is to use an set.  this will take more memory but probably will be faster than sort for large buckets (high repeats).
        typename Base::Base::template UniqueKeySetUtilityType<Key> unique_set(this->c.size());
        auto max = this->c.end();
        for (auto it = this->c.begin(); it != max; ++it) {
          unique_set.emplace(it->first);
        }
        uniq_count = unique_set.size();
        if (uniq_count == 0) {
        	this->key_multiplicity = 1;
        } else {

        	this->key_multiplicity = (this->c.size() + uniq_count - 1) / uniq_count + 1;
        }
        //printf("%lu elements, %lu buckets, %lu unique, key multiplicity = %lu\n", this->c.size(), this->c.bucket_count(), uniq_count, this->key_multiplicity);


        //        // third approach is to assume each bucket contains only 1 kmer/kmolecule.
        //        // This is not generally true for all hash functions, so this is an over estimation of the repeat count.
        //        // we equate bucket size to the number of repeats for that key.
        //        // we can use mean, max, or mean+stdev.
        //        // max overestimates significantly with potentially value > 1000, so don't use max.  (0.0078125 human: 50 sec. synth  32 sec)
        //        // mean may be underestimating for well behaving hash function.   (0.0078125 human: 50 sec. synth  32 sec)
        //        // mean + 2 stdev gets 95% of all entries.  1 stdev covers 67% of all entries, which for high coverage genome is probably better.
        //        //    (1 stdev:  0.0078125 human: 49 sec. synth  32 sec;  2stdev: 0.0078125 human 49s synth: 33 sec)
        //        double nBuckets = 0.0;
        //        for (size_t i = 0, max = this->c.bucket_count(); i < max; ++i) {
        //          if (this->c.bucket_size(i) > 0) nBuckets += 1.0;
        //        }
        //        double mean = static_cast<double>(this->c.size()) / nBuckets;
        //        // do stdev = sqrt((1/nBuckets)  * sum((x - u)^2)).  value is more centered compared to summing the square of x.
        //        double stdev = 0.0;
        //        double entry = 0;
        //        for (size_t i = 0, max = this->c.bucket_count(); i < max; ++i) {
        //          if (this->c.bucket_size(i) == 0) continue;
        //          entry = static_cast<double>(this->c.bucket_size(i)) - mean;
        //          stdev += (entry * entry);
        //        }
        //        stdev = ::std::sqrt(stdev / nBuckets);
        //        this->key_multiplicity = ::std::ceil(mean + 1.0 * stdev);  // covers 95% of data.
        //        printf("%lu elements, %lu buckets, %f occupied, mean = %f, stdev = %f, key multiplicity = %lu\n", this->c.size(), this->c.bucket_count(), nBuckets, mean, stdev, this->key_multiplicity);

        // finally, hard coding.  (0.0078125 human:  50 sec.  synth:  32 s)
        // this->key_multiplicity = 50;

        return this->key_multiplicity;
      }


//
//
//      size_t count_unique(::std::vector<::std::pair<Key, T> > const & input) const {
//        // alternative approach to get number of unique keys is to use an unordered_set.  this will take more memory but probably will be faster than sort for large buckets (high repeats).
//        ::std::unordered_set<Key, typename Base::TransformedHash, typename Base::Base::TransformedEqual > unique_set(this->c.size());
//        for (auto it = input.begin(), max = input.end(); it != max; ++it) {
//          unique_set.insert(it->first);
//        }
//       // printf("r %d: %lu elements, %lu unique\n", this->comm_rank, input.size(), unique_set.size());
//        return unique_set.size();
//      }
//
//      template <typename _TargetP>
//      ::std::vector<::std::pair<Key, T> > bucketing(::std::vector<::std::pair<Key, T> > const & msgs, _TargetP target_p_fun, MPI_Comm comm) {
//
//        int p;
//        MPI_Comm_size(comm, &p);
//
//        // bucket input by their target processor
//        // TODO: in-place bucketing??
//        std::vector<int> send_counts(p, 0);
//        std::vector<int> pids(msgs.size());
//        for (int i = 0; i < msgs.size(); ++i)
//        {
//          pids[i] = target_p_fun(msgs[i]);
//          send_counts[pids[i]]++;
//        }
//
//        // get all2all params
//        std::vector<int> offset = mxx::get_displacements(send_counts);
//
//        // copy.  need to be able to track current position within each block.
//        ::std::vector<::std::pair<Key, T> > send_buffer;
//        if (msgs.size() > 0)
//          send_buffer.resize(msgs.size());
//        for (int i = 0; i < msgs.size(); ++i)
//        {
//          send_buffer[offset[pids[i]]++] = msgs[i];
//        }
//        return send_buffer;
//      }


      /**
       * @brief insert new elements in the distributed unordered_multimap.
       * @param first
       * @param last
       */
      template <typename Predicate = TruePredicate>
      size_t insert(std::vector<::std::pair<Key, T> >& input, bool sorted_input = false, Predicate const & pred = Predicate()) {
        // even if count is 0, still need to participate in mpi calls.  if (input.size() == 0) return;
        BL_BENCH_INIT(insert);


        //        printf("r %d key size %lu, val size %lu, pair size %lu, tuple size %lu\n", this->comm.rank(), sizeof(Key), sizeof(T), sizeof(::std::pair<Key, T>), sizeof(::std::tuple<Key, T>));
        //        count_unique(input);
        //        count_unique(bucketing(input, this->key_to_rank, this->comm));

        // communication part
        if (this->comm.size() > 1) {
            BL_BENCH_START(insert);
            // first remove duplicates.  sort, then get unique, finally remove the rest.  may not be needed
            auto recv_counts(::dsc::distribute(input, this->key_to_rank, sorted_input, this->comm));

            BLISS_UNUSED(recv_counts);
            BL_BENCH_END(insert, "dist_data", input.size());
        }

        //        count_unique(input);

        BL_BENCH_START(insert);
        // local compute part.  called by the communicator.
        size_t count = 0;
        if (!::std::is_same<Predicate, TruePredicate>::value)
          count = this->Base::local_insert(input.begin(), input.end(), pred);
        else
          count = this->Base::local_insert(input.begin(), input.end());
        BL_BENCH_END(insert, "insert", this->c.size());

        BL_BENCH_REPORT_MPI_NAMED(insert, "hash_multimap:insert", this->comm);
        return count;
      }

  };



  /**
   * @brief  distributed unordered reduction map following std unordered map's interface.  Insertion applies the binary reduction operator between the existing and inserted element (in that order).
   * @details   This class is modeled after the std::map, but allows a binary reduction operator to be used during insertion.
   *
   *         the reduction operator is not assumed to be associative.  The operator is called with parameters existing element, then new element to insert.
   *
   *         it has as much of the same methods of std::map as possible.  however, all methods consider the fact
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
   * @tparam Alloc  default to ::std::allocator< ::std::pair<const Key, T> >    allocator for local storage.
   */
  template<typename Key, typename T,
  class Comm,
  template <typename> class KeyTransform,
  template <typename, bool> class Hash,
  typename Reduc = ::std::plus<T>,
  class Less = ::std::less<Key>,
  class Alloc = ::std::allocator< ::std::pair<const Key, T> >
  >
  class reduction_map : public map<Key, T, Comm, KeyTransform, Hash, Less, Alloc> {
      static_assert(::std::is_arithmetic<T>::value, "mapped type has to be arithmetic");

    protected:
      using Base = map<Key, T, Comm, KeyTransform, Hash, Less, Alloc>;


    public:
      using local_container_type = typename Base::local_container_type;

      // std::multimap public members.
      using key_type              = typename local_container_type::key_type;
      using mapped_type           = typename local_container_type::mapped_type;
      using value_type            = typename local_container_type::value_type;
      using iterator              = typename local_container_type::iterator;
      using const_iterator        = typename local_container_type::const_iterator;
      using size_type             = typename local_container_type::size_type;
      using difference_type       = typename local_container_type::difference_type;

    protected:
      Reduc r;

      /**
       * @brief insert new elements in the distributed multimap.
       * @param first
       * @param last
       */
      template <class InputIterator>
      size_t local_insert(InputIterator first, InputIterator last) {
          size_t before = this->c.size();
          for (auto it = first; it != last; ++it) {
            this->c[it->first] = r(this->c[it->first], it->second);
          }
          return this->c.size() - before;
      }

      /**
       * @brief insert new elements in the distributed multimap.
       * @param first
       * @param last
       */
      template <class InputIterator, class Predicate>
      size_t local_insert(InputIterator first, InputIterator last, Predicate const & pred) {
          size_t before = this->c.size();

          for (auto it = first; it != last; ++it) {
            if (pred(*it)) this->c[it->first] = r(this->c[it->first], it->second);
          }
          return this->c.size() - before;

      }

      virtual void local_reduction(::std::vector<::std::pair<Key, T> >& input, bool& sorted_input) {

        if (input.size() == 0) return;

        // sort is slower.  use unordered map.
        BL_BENCH_INIT(reduce_tuple);

        BL_BENCH_START(reduce_tuple);
        local_container_type temp;
        BL_BENCH_END(reduce_tuple, "reserve", input.size());

        BL_BENCH_START(reduce_tuple);
        auto end = input.end();
        for (auto it = input.begin(); it != end; ++it) {
          if (temp.count(it->first) == 0) temp[it->first] = it->second;  // don't rely on initialization to set T to 0.
          else temp[it->first] = r(temp[it->first], it->second);
        }
        BL_BENCH_END(reduce_tuple, "reduce", temp.size());

        BL_BENCH_START(reduce_tuple);
        input.assign(temp.begin(), temp.end());
        BL_BENCH_END(reduce_tuple, "copy", input.size());

        BL_BENCH_REPORT_MPI_NAMED(reduce_tuple, "reduction_map:local_reduce", this->comm);
      }


    public:


      reduction_map(::mxx::comm const & _comm) : Base(_comm) {}

      virtual ~reduction_map() {};

      /**
       * @brief insert new elements in the distributed multimap.
       * @param first
       * @param last
       */
      template <typename Predicate = TruePredicate>
      size_t insert(std::vector<::std::pair<Key, T> >& input, bool sorted_input = false, Predicate const & pred = Predicate()) {
        // even if count is 0, still need to participate in mpi calls.  if (input.size() == 0) return;

        // even if count is 0, still need to participate in mpi calls.  if (input.size() == 0) return;
        BL_BENCH_INIT(insert);

        BL_BENCH_START(insert);
        BL_BENCH_END(insert, "begin", input.size());

        // communication part
        if (this->comm.size() > 1) {

          BL_BENCH_START(insert);
          // first remove duplicates.  sort, then get unique, finally remove the rest.  may not be needed
          this->local_reduction(input, sorted_input);
          BL_BENCH_END(insert, "local_reduce", input.size());

          BL_BENCH_START(insert);
          // first remove duplicates.  sort, then get unique, finally remove the rest.  may not be needed
          auto recv_counts(::dsc::distribute(input, this->key_to_rank, sorted_input, this->comm));
          BLISS_UNUSED(recv_counts);
          BL_BENCH_END(insert, "dist_data", input.size());
        }
        //
        //        // after communication, sort again to keep unique  - may not be needed
        //        local_reduction(input);

        // local compute part.  called by the communicator.
        BL_BENCH_START(insert);
        size_t count = 0;
        if (!::std::is_same<Predicate, TruePredicate>::value)
          count = this->local_insert(input.begin(), input.end(), pred);
        else
          count = this->local_insert(input.begin(), input.end());
        BL_BENCH_END(insert, "local_insert", this->local_size());

        BL_BENCH_REPORT_MPI_NAMED(insert, "reduction_map:insert", this->comm);

        return count;

      }



  };




  /**
   * @brief  distributed unordered counting map following std unordered map's interface.  Insertion applies the binary reduction operator between the existing and inserted element (in that order).
   * @details   This class is modeled after the std::map, but allows a binary reduction operator to be used during insertion.
   *
   *         the reduction operator is not assumed to be associative.  The operator is called with parameters existing element, then new element to insert.
   *
   *         it has as much of the same methods of std::map as possible.  however, all methods consider the fact
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
   * @tparam Alloc  default to ::std::allocator< ::std::pair<const Key, T> >    allocator for local storage.
   */
  template<typename Key, typename T,
  class Comm,
  template <typename> class KeyTransform,
  template <typename, bool> class Hash,
  class Less = ::std::less<Key>,
  class Alloc = ::std::allocator< ::std::pair<const Key, T> >
  >
  class counting_map : public reduction_map<Key, T, Comm, KeyTransform, Hash, ::std::plus<T>, Less,Alloc> {
      static_assert(::std::is_integral<T>::value, "count type has to be integral");

    protected:
      using Base = reduction_map<Key, T, Comm, KeyTransform, Hash, ::std::plus<T>, Less, Alloc>;

    public:
      using local_container_type = typename Base::local_container_type;

      // std::multimap public members.
      using key_type              = typename local_container_type::key_type;
      using mapped_type           = typename local_container_type::mapped_type;
      using value_type            = typename local_container_type::value_type;
      using iterator              = typename local_container_type::iterator;
      using const_iterator        = typename local_container_type::const_iterator;
      using size_type             = typename local_container_type::size_type;
      using difference_type       = typename local_container_type::difference_type;


      counting_map(::mxx::comm const &_comm) : Base(_comm) {}

      virtual ~counting_map() {};


      /**
       * @brief insert new elements in the distributed multimap.
       * @param first
       * @param last
       */
      template <typename Predicate = TruePredicate>
      size_t insert(std::vector< Key >& input, bool sorted_input = false, Predicate const &pred = Predicate()) {
        // even if count is 0, still need to participate in mpi calls.  if (input.size() == 0) return;
        BL_BENCH_INIT(count_insert);

        BL_BENCH_START(count_insert);
        ::std::vector<::std::pair<Key, T> > temp;
        temp.reserve(input.size());
        ::fsc::back_emplace_iterator<::std::vector<::std::pair<Key, T> > > emplace_iter(temp);
        ::std::transform(input.begin(), input.end(), emplace_iter, [](Key const & x) { return ::std::make_pair(x, T(1)); });
        BL_BENCH_END(count_insert, "convert", input.size());


        BL_BENCH_START(count_insert);
        // local compute part.  called by the communicator.
        size_t count = this->Base::insert(temp, sorted_input, pred);
        ::std::vector<::std::pair<Key, T> >().swap(temp);  // clear the temp.
        BL_BENCH_END(count_insert, "insert", this->c.size());


        // distribute
        BL_BENCH_REPORT_MPI_NAMED(count_insert, "count_map:insert_key", this->comm);

        return count;

      }



      /**
       * @brief insert new elements in the distributed multimap.
       * @param first
       * @param last
       */
      template <typename Predicate = TruePredicate>
      size_t insert(std::vector< std::pair< Key, T> >& input, bool sorted_input = false, Predicate const &pred = Predicate()) {
        // even if count is 0, still need to participate in mpi calls.  if (input.size() == 0) return;
        // local compute part.  called by the communicator.
        size_t count = this->Base::insert(input, sorted_input, pred);

        return count;

      }

  };

} /* namespace dsc */


#endif // BLISS_DISTRIBUTED_MAP_HPP
