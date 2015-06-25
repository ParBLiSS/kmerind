/**
 * @file    distributed_unordered_map.hpp
 * @ingroup index
 * @author  Tony Pan <tpan7@gatech.edu>
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   Implements the distributed_multimap, distributed map, and distributed_reduction_map
 *          data structures.
 *
 *          implementation is hash-base (O(1) lookup). later will support sort-based (load balanced).
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
 *
 * Copyright (c) 2014 Georgia Institute of Technology. All Rights Reserved.
 *
 * TODO add Licence
 */

#ifndef BLISS_DISTRIBUTED_UNORDERED_MAP_HPP
#define BLISS_DISTRIBUTED_UNORDERED_MAP_HPP


#include <unordered_map>  // local storage hash table  // for multimap
#include <unordered_set>  // local storage hash table  // for multimap
#include <utility> 			  // for std::pair

//#include <sparsehash/dense_hash_map>  // not a multimap, where we need it most.
#include <functional> 		// for std::function and std::hash
#include <algorithm> 		// for sort, stable_sort, unique, is_sorted
#include <iterator>  // advance, distance

#include <cstdint>  // for uint8, etc.

#include <type_traits>

#include "mxx/collective.hpp"
#include "mxx/reduction.hpp"
#include "io/mpi_utils.hpp"
#include "utils/timer.hpp"  // for timing.
#include "utils/logging.h"

#include "wip/distributed_map_base.hpp"
#include "wip/unordered_multimap.hpp"

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
   * @tparam Container  default to unordered_map and unordered multimap, requiring 5 template params.
   * @tparam Comm   default to mpi_collective_communicator       communicator for global communication. may hash or sort.
   * @tparam KeyTransform   transform function for the key.  can supply identity.  requires a single template argument (Key).  useful for mapping kmolecule to kmer.
   * @tparam Hash   hash function for local and distribution.  requires a template arugment (Key), and a bool (prefix, chooses the MSBs of hash instead of LSBs)
   * @tparam Equal   default to ::std::equal_to<Key>   equal function for the local storage.
   * @tparam Alloc  default to ::std::allocator< ::std::pair<const Key, T> >    allocator for local storage.
   */
  template<typename Key, typename T,
    template <typename, typename, typename, typename, typename> class Container,
    class Comm,
    template <typename> class KeyTransform,
    template <typename, bool> class Hash,
    class Equal = ::std::equal_to<Key>,
    class Alloc = ::std::allocator< ::std::pair<const Key, T> >
  > class unordered_map_base : public ::dsc::map_base<Key, T, Comm, KeyTransform, ::std::less<Key>, Equal, Alloc> {

    protected:
      using Base = ::dsc::map_base<Key, T, Comm, KeyTransform, ::std::less<Key>, Equal, Alloc>;

      struct TransformedHash {
          Hash<Key, false> h;

          inline uint64_t operator()(Key const& k) const {
            return h(Base::trans(k));
          }
          template<typename V>
          inline uint64_t operator()(::std::pair<Key, V> const& x) const {
            return this->operator()(x.first);
          }
          template<typename V>
          inline uint64_t operator()(::std::pair<const Key, V> const& x) const {
            return this->operator()(x.first);
          }
      } hash;

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
          template <class DB, class QueryIter, class OutputIter, class Operator, class Predicate = Identity>
          static size_t process(DB &db,
                                QueryIter query_begin, QueryIter query_end,
                                OutputIter &output, Operator const & op,
                                bool sorted_query = false, Predicate const &pred = Predicate()) {

              if (query_begin == query_end) return 0;

              typename ::std::iterator_traits<QueryIter>::value_type v;
              size_t count = 0;  // before size.
              for (auto it = query_begin; it != query_end; ++it) {
                v = *it;
                if (!::std::is_same<Predicate, Identity>::value)
                  count += op(db, v, output, pred);
                else
                  count += op(db, v, output);
              }
              return count;
          }

      };




    public:
      using local_container_type = Container<Key, T, TransformedHash, typename Base::TransformedEqual, Alloc>;

      // std::unordered_multimap public members.
      using key_type              = typename local_container_type::key_type;
      using mapped_type           = typename local_container_type::mapped_type;
      using value_type            = typename local_container_type::value_type;
      using hasher                = typename local_container_type::hasher;
      using key_equal             = typename local_container_type::key_equal;
      using allocator_type        = typename local_container_type::allocator_type;
      using iterator              = typename local_container_type::iterator;
      using const_iterator        = typename local_container_type::const_iterator;
      using size_type             = typename local_container_type::size_type;
      using difference_type       = typename local_container_type::difference_type;

    protected:
      local_container_type c;

      /// reserve space.  n is the local container size.  this allows different processes to individually adjust its own size.
      void local_reserve( size_t n) {
        c.reserve(n);
      }

      /// rehash the local container.  n is the local container size.  this allows different processes to individually adjust its own size.
      void local_rehash( size_type n) {
      }

      struct LocalCount {
          // unfiltered.
          template<class DB, typename Query, class OutputIter>
          size_t operator()(DB &db, Query const &v, OutputIter &output) const {
              *output = ::std::move(::std::make_pair(v, db.count(v)));
              ++output;
              return 1;
          }
          // filtered element-wise.
          template<class DB, typename Query, class OutputIter, class Predicate = Identity>
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
          template<class DB, typename Query, class OutputIter, class Predicate = Identity>
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
       * @brief insert new elements in the distributed unordered_multimap.
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
       * @brief insert new elements in the distributed unordered_multimap.  example use: stop inserting if more than x entries.
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


      /// clears the unordered_map
      virtual void local_clear() noexcept {
          c.clear();
      }


      ///  keep the unique keys in the input. primarily for reducing comm volume.   output is sorted.  equal operator forces comparison to Key
      template <typename V>
      void retain_unique(::std::vector< V >& input, bool sorted_input = false) const {
        if (input.size() == 0) return;
        if (sorted_input) {  // already sorted, then just get the unique stuff and remove rest.
          auto end = ::std::unique(input.begin(), input.end(), this->equal);
          input.erase(end, input.end());
        } else {  // not sorted, so use an unordered_set to keep the first occurence.

          // sorting is SLOW and not scalable.  use unordered set instead.

          ::std::unordered_set<V, TransformedHash, typename Base::TransformedEqual > temp(input.begin(), input.end(), input.size());
          input.assign(temp.begin(), temp.end());
        }
      }

      /**
       * @brief find elements with the specified keys in the distributed unordered_multimap.
       * @param first
       * @param last
       */
      template <class LocalFind, typename Predicate = Identity>
      ::std::vector<::std::pair<Key, T> > find(LocalFind const & find_element, ::std::vector<Key>& keys, bool sorted_input = false, Predicate const& pred = Predicate()) const {
        TIMER_INIT(find);

        TIMER_START(find);
        ::std::vector<::std::pair<Key, T> > results;
        back_emplace_iterator<::std::vector<::std::pair<Key, T> > > emplace_iter(results);
        // even if count is 0, still need to participate in mpi calls.  if (keys.size() == 0) return results;
        TIMER_END(find, "begin", keys.size());

        TIMER_START(find);
        // keep unique keys
        this->retain_unique(keys, sorted_input);
        TIMER_END(find, "uniq1", keys.size());

        if (this->comm_size > 1) {
          TIMER_START(find);
          // distribute (communication part)
          std::vector<int> recv_counts = mxx2::msgs_all2all(keys, this->key_to_rank, this->comm);
          TIMER_END(find, "a2a1", keys.size());

          // local find. memory utilization a potential problem.
          // do for each src proc one at a time.

          TIMER_START(find);
          std::vector<int> send_counts(this->comm_size, 0);
          results.reserve(keys.size() * this->key_multiplicity);                   // TODO:  should estimate coverage.
          //printf("reserving %lu\n", keys.size() * this->key_multiplicity);
          TIMER_END(find, "reserve", (keys.size() * this->key_multiplicity));

          TIMER_START(find);
          auto start = keys.begin();
          auto end = start;
          for (int i = 0; i < this->comm_size; ++i) {
            ::std::advance(end, recv_counts[i]);

            // work on query from process i.
            QueryProcessor::process(c, start, end, emplace_iter, find_element, true, pred);
           // if (this->comm_rank == 0) DEBUGF("R %d added %d results for %d queries for process %d\n", this->comm_rank, send_counts[i], recv_counts[i], i);

            start = end;
          }
          TIMER_END(find, "local_find", results.size());

          TIMER_START(find);
          // send back using the constructed recv count
          mxx::all2all(results, send_counts, this->comm);
          TIMER_END(find, "a2a2", results.size());

        } else {
          TIMER_START(find);
          results.reserve(keys.size() * this->key_multiplicity);                   // TODO:  should estimate coverage.
          //printf("reserving %lu\n", keys.size() * this->key_multiplicity);
          TIMER_END(find, "reserve", (keys.size() * this->key_multiplicity));

          TIMER_START(find);
          QueryProcessor::process(c, keys.begin(), keys.end(), emplace_iter, find_element, true, pred);
          TIMER_END(find, "local_find", results.size());
        }

        TIMER_REPORT_MPI(find, this->comm_rank, this->comm);

        return results;

      }

      template <class LocalFind, typename Predicate = Identity>
      ::std::vector<::std::pair<Key, T> > find(LocalFind const & find_element, Predicate const& pred = Predicate()) const {
        ::std::vector<::std::pair<Key, T> > results;
        back_emplace_iterator<::std::vector<::std::pair<Key, T> > > emplace_iter(results);

        auto keys = this->keys();
        results.reserve(keys.size() * this->key_multiplicity);                   // TODO:  should estimate coverage.

        QueryProcessor::process(c, keys.begin(), keys.end(), emplace_iter, find_element, true, pred);
        return results;
      }


      unordered_map_base(MPI_Comm _comm, int _comm_size) : Base(_comm, _comm_size),
          key_to_rank(_comm_size) {}


    public:

      virtual ~unordered_map_base() {};

      /// returns the local storage.  please use sparingly.
      local_container_type& get_local_container() { return c; }

      /// update the multiplicity.  only multimap needs to do this.
      virtual size_t update_multiplicity() { return this->key_multiplicity; }

      /// convert the map to a vector.
      virtual std::vector<std::pair<Key, T> > to_vector() const {
        std::vector<std::pair<Key, T> > result;
        this->to_vector(result);
        return result;
      }

      /// convert the map to a vector
      virtual void to_vector(std::vector<std::pair<Key, T> > & result) const {
        result.clear();
        if (c.empty()) return;
        result.reserve(c.size());
        ::dsc::back_emplace_iterator<::std::vector<::std::pair<Key, T> > > emplace_iter(result);
        ::std::copy(c.begin(), c.end(), emplace_iter);
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

        ::std::unordered_set<Key, TransformedHash, typename Base::TransformedEqual > temp;
        temp.reserve(c.size());
        for (auto it = c.begin(), end = c.end(); it != end; ++it) {
          temp.emplace(it->first);
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
      void reserve( size_t n) {
        // direct reserve + barrier
        this->local_reserve(n);
        if (this->comm_size > 1) MPI_Barrier(this->comm);
      }


      /// rehash the local container.  n is the local container size.  this allows different processes to individually adjust its own size.
      void rehash( size_type n) {
        // direct rehash + barrier
        local_rehash(n);
        if (this->comm_size > 1) MPI_Barrier(this->comm);
      }


      /**
       * @brief count elements with the specified keys in the distributed unordered_multimap.
       * @param first
       * @param last
       */
      template <class Predicate = Identity>
      ::std::vector<::std::pair<Key, size_type> > count(::std::vector<Key>& keys, bool sorted_input = false,
                                                        Predicate const& pred = Predicate() ) const {
        TIMER_INIT(count);

        TIMER_START(count);
        ::std::vector<::std::pair<Key, size_type> > results;
        back_emplace_iterator<::std::vector<::std::pair<Key, size_type> > > emplace_iter(results);
        // even if count is 0, still need to participate in mpi calls.  if (keys.size() == 0) return results;
        TIMER_END(count, "begin", keys.size());

        TIMER_START(count);
        // keep unique keys
        retain_unique(keys, sorted_input);
        TIMER_END(count, "uniq1", keys.size());


        if (this->comm_size > 1) {
          TIMER_START(count);
          // distribute (communication part)
          std::vector<int> recv_counts = mxx2::msgs_all2all(keys, this->key_to_rank, this->comm);
          TIMER_END(count, "a2a1", keys.size());

          // local count. memory utilization a potential problem.
          // do for each src proc one at a time.
          TIMER_START(count);
          results.reserve(keys.size());                   // TODO:  should estimate coverage.
          TIMER_END(count, "reserve", keys.size());

          TIMER_START(count);
          auto start = keys.begin();
          auto end = start;
          for (int i = 0; i < this->comm_size; ++i) {
            ::std::advance(end, recv_counts[i]);

            // within start-end, values are unique, so don't need to set unique to true.
            QueryProcessor::process(c, start, end, emplace_iter, count_element, true, pred);

            if (this->comm_rank == 0) DEBUGF("R %d added %d results for %d queries for process %d\n", this->comm_rank, send_counts[i], recv_counts[i], i);

            start = end;
          }
          TIMER_END(count, "local_count", results.size());

          // send back using the constructed recv count
          TIMER_START(count);
          mxx::all2all(results, recv_counts, this->comm);
          TIMER_END(count, "a2a2", results.size());
        } else {
          TIMER_START(count);
          results.reserve(keys.size());                   // TODO:  should estimate coverage.
          TIMER_END(count, "reserve", keys.size());


          TIMER_START(count);
          // within start-end, values are unique, so don't need to set unique to true.
          QueryProcessor::process(c, keys.begin(), keys.end(), emplace_iter, count_element, true, pred);
          TIMER_END(count, "local_count", results.size());
        }

        TIMER_REPORT_MPI(count, this->comm_rank, this->comm);

        return results;

      }


      template <typename Predicate = Identity>
      ::std::vector<::std::pair<Key, size_type> > count(Predicate const & pred = Predicate()) const {
        ::std::vector<::std::pair<Key, size_type> > results;
        back_emplace_iterator<::std::vector<::std::pair<Key, size_t> > > emplace_iter(results);

        auto keys = this->keys();
        results.reserve(keys.size());

        QueryProcessor::process(c, keys.begin(), keys.end(), emplace_iter, count_element, true, pred);

        if (this->comm_size > 1) MPI_Barrier(this->comm);
        return results;
      }



      /**
       * @brief erase elements with the specified keys in the distributed unordered_multimap.
       * @param first
       * @param last
       */
      template <class Predicate = Identity>
      size_t erase(::std::vector<Key>& keys, bool sorted_input = false, Predicate const& pred = Predicate() ) {
        // even if count is 0, still need to participate in mpi calls.  if (keys.size() == 0) return;
          size_t before = this->c.size();

          bool si = sorted_input;
        if (this->comm_size > 1) {

          // remove duplicates
          retain_unique(keys, si);

          // redistribute keys
          mxx2::msgs_all2all(keys, this->key_to_rank, this->comm);
          si = false;
        }

        // remove duplicates
        retain_unique(keys, si);

        // then call local remove.
        size_t count = QueryProcessor::process(c, keys.begin(), keys.end(), keys.end(), erase_element, true, pred);

        return before - this->c.size();
      }


      template <typename Predicate = Identity>
      size_t erase(Predicate const & pred = Predicate()) {
        auto keys = this->keys();

        size_t count = QueryProcessor::process(c, keys.begin(), keys.end(), keys.end(), erase_element, true, pred);

        if (this->comm_size > 1) MPI_Barrier(this->comm);

        return count;
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
   * @tparam Alloc  default to ::std::allocator< ::std::pair<const Key, T> >    allocator for local storage.
   */
  template<typename Key, typename T,
  class Comm,
  template <typename> class KeyTransform,
  template <typename, bool> class Hash,
  class Equal = ::std::equal_to<Key>,
  class Alloc = ::std::allocator< ::std::pair<const Key, T> >
  >
  class unordered_map : public unordered_map_base<Key, T, ::std::unordered_map, Comm, KeyTransform, Hash, Equal, Alloc> {
      using Base = unordered_map_base<Key, T, ::std::unordered_map, Comm, KeyTransform, Hash, Equal, Alloc>;


    public:
      using local_container_type = typename Base::local_container_type;

      // std::unordered_multimap public members.
      using key_type              = typename local_container_type::key_type;
      using mapped_type           = typename local_container_type::mapped_type;
      using value_type            = typename local_container_type::value_type;
      using hasher                = typename local_container_type::hasher;
      using key_equal             = typename local_container_type::key_equal;
      using allocator_type        = typename local_container_type::allocator_type;
      using reference             = typename local_container_type::reference;
      using const_reference       = typename local_container_type::const_reference;
      using pointer               = typename local_container_type::pointer;
      using const_pointer         = typename local_container_type::const_pointer;
      using iterator              = typename local_container_type::iterator;
      using const_iterator        = typename local_container_type::const_iterator;
      using size_type             = typename local_container_type::size_type;
      using difference_type       = typename local_container_type::difference_type;

    protected:

      // defined Communicator as a friend
      friend Comm;

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
          template<class DB, typename Query, class OutputIter, class Predicate = Identity>
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


      virtual void local_reduction(std::vector<::std::pair<Key, T> > &input, bool sorted_input = false) {
        this->Base::retain_unique(input, sorted_input);
      }


    public:

      unordered_map(MPI_Comm _comm, int _comm_size) : Base(_comm, _comm_size) {}

      virtual ~unordered_map() {};

      template <class Predicate = Identity>
      ::std::vector<::std::pair<Key, T> > find(::std::vector<Key>& keys, bool sorted_input = false,
          Predicate const& pred = Predicate()) const {
          return Base::find(find_element, keys, sorted_input, pred);
      }

      template <class Predicate = Identity>
      ::std::vector<::std::pair<Key, T> > find(Predicate const& pred = Predicate()) const {
          return Base::find(find_element, pred);
      }


      /**
       * @brief insert new elements in the distributed unordered_multimap.
       * @param first
       * @param last
       */
      template <typename Predicate = Identity>
      size_t insert(std::vector<::std::pair<Key, T> >& input, bool sorted_input = false, Predicate const & pred = Predicate()) {
        // even if count is 0, still need to participate in mpi calls.  if (input.size() == 0) return;
        TIMER_INIT(insert);

        TIMER_START(insert);
        TIMER_END(insert, "start", input.size());


        // communication part
        if (this->comm_size > 1) {
          TIMER_START(insert);
          // first remove duplicates.  sort, then get unique, finally remove the rest.  may not be needed
          local_reduction(input, sorted_input);
          TIMER_END(insert, "uniq1", input.size());

          TIMER_START(insert);
          mxx2::msgs_all2all(input, this->key_to_rank, this->comm);
          TIMER_END(insert, "a2a", input.size());
        }


        TIMER_START(insert);
        // local compute part.  called by the communicator.
        size_t count = 0;
        if (!::std::is_same<Predicate, Identity>::value)
          count = this->Base::local_insert(input.begin(), input.end(), pred);
        else
          count = this->Base::local_insert(input.begin(), input.end());
        TIMER_END(insert, "insert", this->c.size());

        TIMER_REPORT_MPI(insert, this->comm_rank, this->comm);

        return count;
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
   * @tparam Alloc  default to ::std::allocator< ::std::pair<const Key, T> >    allocator for local storage.
   */
  template<typename Key, typename T,
  class Comm,
  template <typename> class KeyTransform,
  template <typename, bool> class Hash,
  class Equal = ::std::equal_to<Key>,
  class Alloc = ::std::allocator< ::std::pair<const Key, T> >
  >
  class unordered_multimap : public unordered_map_base<Key, T, ::std::unordered_multimap, Comm, KeyTransform, Hash, Equal, Alloc> {
      using Base = unordered_map_base<Key, T, ::std::unordered_multimap, Comm, KeyTransform, Hash, Equal, Alloc>;


    public:
      using local_container_type = typename Base::local_container_type;

      // std::unordered_multimap public members.
      using key_type              = typename local_container_type::key_type;
      using mapped_type           = typename local_container_type::mapped_type;
      using value_type            = typename local_container_type::value_type;
      using hasher                = typename local_container_type::hasher;
      using key_equal             = typename local_container_type::key_equal;
      using allocator_type        = typename local_container_type::allocator_type;
      using reference             = typename local_container_type::reference;
      using const_reference       = typename local_container_type::const_reference;
      using pointer               = typename local_container_type::pointer;
      using const_pointer         = typename local_container_type::const_pointer;
      using iterator              = typename local_container_type::iterator;
      using const_iterator        = typename local_container_type::const_iterator;
      using size_type             = typename local_container_type::size_type;
      using difference_type       = typename local_container_type::difference_type;

    protected:

      // defined Communicator as a friend
      friend Comm;

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
          template<class DB, typename Query, class OutputIter, class Predicate = Identity>
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


      unordered_multimap(MPI_Comm _comm, int _comm_size) : Base(_comm, _comm_size) {
        this->key_multiplicity = 50;
      }

      virtual ~unordered_multimap() {}

      template <class Predicate = Identity>
      ::std::vector<::std::pair<Key, T> > find(::std::vector<Key>& keys, bool sorted_input = false,
          Predicate const& pred = Predicate()) const {
          return Base::find(find_element, keys, sorted_input, pred);
      }

      template <class Predicate = Identity>
      ::std::vector<::std::pair<Key, T> > find(Predicate const& pred = Predicate()) const {
          return Base::find(find_element, pred);
      }


      /// update the multiplicity.  only multimap needs to do this.
      virtual size_t update_multiplicity() const {
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
        // alternative approach to get number of unique keys is to use an unordered_set.  this will take more memory but probably will be faster than sort for large buckets (high repeats).
        ::std::unordered_set<Key, typename Base::TransformedHash, typename Base::Base::TransformedEqual > unique_set(this->c.size());
        for (auto it = this->c.begin(), max = this->c.end(); it != max; ++it) {
          unique_set.emplace(it->first);
        }
        uniq_count = unique_set.size();
        this->key_multiplicity = (this->c.size() + uniq_count - 1) / uniq_count + 1;
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
      template <typename Predicate = Identity>
      size_t insert(std::vector<::std::pair<Key, T> >& input, bool sorted_input = false, Predicate const & pred = Predicate()) {
        // even if count is 0, still need to participate in mpi calls.  if (input.size() == 0) return;
        TIMER_INIT(insert);

        TIMER_START(insert);

        //        printf("r %d key size %lu, val size %lu, pair size %lu, tuple size %lu\n", this->comm_rank, sizeof(Key), sizeof(T), sizeof(::std::pair<Key, T>), sizeof(::std::tuple<Key, T>));
        //        count_unique(input);
        //        count_unique(bucketing(input, this->key_to_rank, this->comm));

        // communication part
        if (this->comm_size > 1) {
          mxx2::msgs_all2all(input, this->key_to_rank, this->comm);
        }
        TIMER_END(insert, "a2a", input.size());

        //        count_unique(input);

        TIMER_START(insert);
        // local compute part.  called by the communicator.
        size_t count = 0;
        if (!::std::is_same<Predicate, Identity>::value)
          count = this->Base::local_insert(input.begin(), input.end(), pred);
        else
          count = this->Base::local_insert(input.begin(), input.end());
        TIMER_END(insert, "insert", this->c.size());

        TIMER_REPORT_MPI(insert, this->comm_rank, this->comm);
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
   * @tparam Alloc  default to ::std::allocator< ::std::pair<const Key, T> >    allocator for local storage.
   */
  template<typename Key, typename T,
  class Comm,
  template <typename> class KeyTransform,
  template <typename, bool> class Hash,
  typename Reduc = ::std::plus<T>,
  class Equal = ::std::equal_to<Key>,
  class Alloc = ::std::allocator< ::std::pair<const Key, T> >
  >
  class reduction_unordered_map : public unordered_map<Key, T, Comm, KeyTransform, Hash, Equal, Alloc> {
      using Base = unordered_map<Key, T, Comm, KeyTransform, Hash, Equal, Alloc>;

      static_assert(::std::is_arithmetic<T>::value, "mapped type has to be arithmetic");

    public:
      using local_container_type = typename Base::local_container_type;

      // std::unordered_multimap public members.
      using key_type              = typename local_container_type::key_type;
      using mapped_type           = typename local_container_type::mapped_type;
      using value_type            = typename local_container_type::value_type;
      using hasher                = typename local_container_type::hasher;
      using key_equal             = typename local_container_type::key_equal;
      using allocator_type        = typename local_container_type::allocator_type;
      using reference             = typename local_container_type::reference;
      using const_reference       = typename local_container_type::const_reference;
      using pointer               = typename local_container_type::pointer;
      using const_pointer         = typename local_container_type::const_pointer;
      using iterator              = typename local_container_type::iterator;
      using const_iterator        = typename local_container_type::const_iterator;
      using size_type             = typename local_container_type::size_type;
      using difference_type       = typename local_container_type::difference_type;

    protected:
      Reduc r;

      // defined Communicator as a friend
      friend Comm;

      /**
       * @brief insert new elements in the distributed unordered_multimap.
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
       * @brief insert new elements in the distributed unordered_multimap.
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

      virtual void local_reduction(::std::vector<::std::pair<Key, T> >& input) {

        if (input.size() == 0) return;

        // sort is slower.  use unordered map.
        TIMER_INIT(reduce_tuple);

        TIMER_START(reduce_tuple);
        local_container_type temp(input.size());
        TIMER_END(reduce_tuple, "reserve", input.size());

        TIMER_START(reduce_tuple);
        for (auto it = input.begin(), end = input.end(); it != end; ++it) {
          if (temp.count(it->first) == 0) temp[it->first] = it->second;  // don't rely on initialization to set T to 0.
          else temp[it->first] = r(temp[it->first], it->second);
        }
        TIMER_END(reduce_tuple, "reduce", temp.size());

        TIMER_START(reduce_tuple);
        input.assign(temp.begin(), temp.end());
        TIMER_END(reduce_tuple, "copy", input.size());

        TIMER_REPORT_MPI(reduce_tuple, this->comm_rank, this->comm);
      }


    public:


      reduction_unordered_map(MPI_Comm _comm, int _comm_size) : Base(_comm, _comm_size) {}

      virtual ~reduction_unordered_map() {};

      /**
       * @brief insert new elements in the distributed unordered_multimap.
       * @param first
       * @param last
       */
      template <typename Predicate = Identity>
      size_t insert(std::vector<::std::pair<Key, T> >& input, bool sorted_input = false, Predicate const & pred = Predicate()) {
        // even if count is 0, still need to participate in mpi calls.  if (input.size() == 0) return;


        // communication part
        if (this->comm_size > 1) {
          // first remove duplicates.  sort, then get unique, finally remove the rest.  may not be needed
          this->local_reduction(input);

          mxx2::msgs_all2all(input, this->key_to_rank, this->comm);
        }
        //
        //        // after communication, sort again to keep unique  - may not be needed
        //        local_reduction(input);

        // local compute part.  called by the communicator.
        size_t count = 0;
        if (!::std::is_same<Predicate, Identity>::value)
          count = this->local_insert(input.begin(), input.end(), pred);
        else
          count = this->local_insert(input.begin(), input.end());

        return count;
      }

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
   * @tparam Alloc  default to ::std::allocator< ::std::pair<const Key, T> >    allocator for local storage.
   */
  template<typename Key, typename T,
  class Comm,
  template <typename> class KeyTransform,
  template <typename, bool> class Hash,
  class Equal = ::std::equal_to<Key>,
  class Alloc = ::std::allocator< ::std::pair<const Key, T> >
  >
  class counting_unordered_map : public reduction_unordered_map<Key, T, Comm, KeyTransform, Hash, ::std::plus<T>, Equal,Alloc> {
      using Base = reduction_unordered_map<Key, T, Comm, KeyTransform, Hash, ::std::plus<T>, Equal, Alloc>;

      static_assert(::std::is_integral<T>::value, "count type has to be integral");

    public:
      using local_container_type = typename Base::local_container_type;

      // std::unordered_multimap public members.
      using key_type              = typename local_container_type::key_type;
      using mapped_type           = typename local_container_type::mapped_type;
      using value_type            = typename local_container_type::value_type;
      using hasher                = typename local_container_type::hasher;
      using key_equal             = typename local_container_type::key_equal;
      using allocator_type        = typename local_container_type::allocator_type;
      using reference             = typename local_container_type::reference;
      using const_reference       = typename local_container_type::const_reference;
      using pointer               = typename local_container_type::pointer;
      using const_pointer         = typename local_container_type::const_pointer;
      using iterator              = typename local_container_type::iterator;
      using const_iterator        = typename local_container_type::const_iterator;
      using size_type             = typename local_container_type::size_type;
      using difference_type       = typename local_container_type::difference_type;

    protected:

      // defined Communicator as a friend
      friend Comm;


    public:
      counting_unordered_map(MPI_Comm _comm, int _comm_size) : Base(_comm, _comm_size) {}

      virtual ~counting_unordered_map() {};


      /**
       * @brief insert new elements in the distributed unordered_multimap.
       * @param first
       * @param last
       */
      template <typename Predicate = Identity>
      size_t insert(std::vector< Key >& input, bool sorted_input = false, Predicate const &pred = Predicate()) {
        // even if count is 0, still need to participate in mpi calls.  if (input.size() == 0) return;
        TIMER_INIT(count_insert);

        TIMER_START(count_insert);
        ::std::vector<::std::pair<Key, T> > temp;
        temp.reserve(input.size());
        back_emplace_iterator<::std::vector<::std::pair<Key, T> > > emplace_iter(temp);
        ::std::transform(input.begin(), input.end(), emplace_iter, [](Key const & x) { return ::std::make_pair(x, T(1)); });
        TIMER_END(count_insert, "convert", input.size());


        TIMER_START(count_insert);
        // local compute part.  called by the communicator.
        size_t count = this->Base::insert(temp, sorted_input, pred);

        TIMER_END(count_insert, "insert", this->c.size());


        // distribute
        TIMER_REPORT_MPI(count_insert, this->comm_rank, this->comm);

        return count;

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
   * @tparam Alloc  default to ::std::allocator< ::std::pair<const Key, T> >    allocator for local storage.
   */
  template<typename Key, typename T,
  class Comm,
  template <typename> class KeyTransform,
  template <typename, bool> class Hash,
  class Equal = ::std::equal_to<Key>,
  class Alloc = ::std::allocator< ::std::pair<const Key, T> >
  >
  class unordered_multimap_vec : public unordered_map_base<Key, T, ::fsc::unordered_multimap, Comm, KeyTransform, Hash, Equal, Alloc> {
      using Base = unordered_map_base<Key, T, ::fsc::unordered_multimap, Comm, KeyTransform, Hash, Equal, Alloc>;


    public:
      using local_container_type = typename Base::local_container_type;

      // std::unordered_multimap public members.
      using key_type              = typename local_container_type::key_type;
      using mapped_type           = typename local_container_type::mapped_type;
      using value_type            = typename local_container_type::value_type;
      using hasher                = typename local_container_type::hasher;
      using key_equal             = typename local_container_type::key_equal;
      using allocator_type        = typename local_container_type::allocator_type;
      using reference             = typename local_container_type::reference;
      using const_reference       = typename local_container_type::const_reference;
      using pointer               = typename local_container_type::pointer;
      using const_pointer         = typename local_container_type::const_pointer;
      using iterator              = typename local_container_type::iterator;
      using const_iterator        = typename local_container_type::const_iterator;
      using size_type             = typename local_container_type::size_type;
      using difference_type       = typename local_container_type::difference_type;

    protected:

      // defined Communicator as a friend
      friend Comm;


      struct LocalFind {
          // unfiltered.
          template<class DB, typename Query, class OutputIter>
          size_t operator()(DB &db, Query const &v, OutputIter &output) const {
              auto range = db.equal_range(v);
              output = ::std::copy(range.first, range.second, output);  // tons faster to emplace - almost 3x faster
              return db.count(v);
          }
          // filtered element-wise.
          template<class DB, typename Query, class OutputIter, class Predicate = Identity>
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


      unordered_multimap_vec(MPI_Comm _comm, int _comm_size) : Base(_comm, _comm_size) {
        this->key_multiplicity = 50;
      }

      virtual ~unordered_multimap_vec() {}

      template <class Predicate = Identity>
      ::std::vector<::std::pair<Key, T> > find(::std::vector<Key>& keys, bool sorted_input = false,
          Predicate const& pred = Predicate()) const {
          return Base::find(find_element, keys, sorted_input, pred);
      }

      template <class Predicate = Identity>
      ::std::vector<::std::pair<Key, T> > find(Predicate const& pred = Predicate()) const {
          return Base::find(find_element, pred);
      }


      /// update the multiplicity.  only multimap needs to do this.
      virtual size_t update_multiplicity() const {
        // one approach is to add up the number of repeats for the key of each entry, then divide by total count.
        //  sum(count per key) / c.size.
        // problem with this approach is that for unordered map, to get the count for a key is essentially O(count), so we get quadratic time.
        // The approach is VERY SLOW for large repeat count.  - (0.0078125 human: 52 sec, synth: FOREVER.)

        // a second approach is to count the number of unique key then divide the map size by that.
        //  c.size / #unique.  requires unique set
        // To find unique set, we take each bucket, copy to vector, sort it, and then count unique.
        // This is precise, and is faster than the approach above.  (0.0078125 human: 54 sec.  synth: 57sec.)
        // but the n log(n) sort still grows with the duplicate count
        size_t uniq_count = this->c.unique_size();
        this->key_multiplicity = (this->c.size() + uniq_count - 1) / uniq_count + 1;
        //printf("%lu elements, %lu buckets, %lu unique, key multiplicity = %lu\n", this->c.size(), this->c.bucket_count(), uniq_count, this->key_multiplicity);

        return this->key_multiplicity;
      }

      /**
       * @brief insert new elements in the distributed unordered_multimap.
       * @param first
       * @param last
       */
      template <typename Predicate = Identity>
      size_t insert(std::vector<::std::pair<Key, T> >& input, bool sorted_input = false, Predicate const & pred = Predicate()) {
        // even if count is 0, still need to participate in mpi calls.  if (input.size() == 0) return;
        TIMER_INIT(insert);

        TIMER_START(insert);

        //        printf("r %d key size %lu, val size %lu, pair size %lu, tuple size %lu\n", this->comm_rank, sizeof(Key), sizeof(T), sizeof(::std::pair<Key, T>), sizeof(::std::tuple<Key, T>));
        //        count_unique(input);
        //        count_unique(bucketing(input, this->key_to_rank, this->comm));

        // communication part
        if (this->comm_size > 1) {
          mxx2::msgs_all2all(input, this->key_to_rank, this->comm);
        }
        TIMER_END(insert, "a2a", input.size());

        //        count_unique(input);

        TIMER_START(insert);
        // local compute part.  called by the communicator.
        size_t count = 0;
        if (!::std::is_same<Predicate, Identity>::value)
          count = this->Base::local_insert(input.begin(), input.end(), pred);
        else
          count = this->Base::local_insert(input.begin(), input.end());
        TIMER_END(insert, "insert", this->c.size());

        TIMER_REPORT_MPI(insert, this->comm_rank, this->comm);
        return count;

      }

  };


} /* namespace dsc */


#endif // BLISS_DISTRIBUTED_UNORDERED_MAP_HPP
