/**
 * @file    distributed_sorted_map.hpp
 * @ingroup index
 * @author  Tony Pan <tpan7@gatech.edu>
 * @author  Patrick Flick <patrick.flick@gmail.com>
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
 *
 * Copyright (c) 2014 Georgia Institute of Technology. All Rights Reserved.
 *
 * TODO add Licence
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

#include "mxx/collective.hpp"
#include "mxx/reduction.hpp"
#include "mxx/sort.hpp"
#include "io/mpi_utils.hpp"
#include "utils/timer.hpp"  // for timing.
#include "utils/logging.h"

namespace dsc  // distributed std container
{
  /// from http://stackoverflow.com/questions/18724999/why-no-emplacement-iterators-in-c11-or-c14
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

  struct Identity {
	  template <typename T>
      inline bool operator()(T const & x) const { return true; };
  };



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
   *       4.
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
  >  class sorted_map_base {

    protected:

      static KeyTransform<Key> trans;

      template <typename Comparator>
      struct TransformedComp {
          Comparator comp;
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


      struct Greater {
          Less lt;
          inline bool operator()(Key const &x, Key const &y) const {
            return lt(y, x);
          }
      };

      using TransformedEqual = TransformedComp<Equal>;
      using TransformedLess = TransformedComp<Less>;
      using TransformedGreater = TransformedComp<Greater>;
      static TransformedEqual equal;
      static TransformedLess less;
      static TransformedGreater greater;


      // splitters.  there should be p-1 entries.

      // uses sorted lookup table to map to target ranks.  assume we store in a vector a pair:  <first kmer on proc, proc id>.  then we can use lower_bound to search.
      // note that this makes in-between values go with the larger proc id.
      struct KeyToRank {
          ::std::vector<::std::pair<Key, unsigned int> > map;
          int p;
          KeyToRank(int _comm_size) : p(_comm_size) {};

          // return id of selected element
          inline unsigned int operator()(Key const & x) const {
            auto pos = ::std::lower_bound(map.begin(), map.end(), x, less);
            return (pos == map.end()) ? (p-1) : pos->second;
          }
          template<typename V>
          inline unsigned int operator()(::std::pair<Key, V> const & x) const {
            return this->operator()(x.first);
          }
          template<typename V>
          inline unsigned int operator()(::std::pair<const Key, V> const & x) const {
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
            if (!sorted_query) sort_ascending(query_begin, query_end);
//            if (!sorted_dest) sort_ascending(db_begin, db_end);  can't sort the container because its begin()/end() keeps returning const iterators if "this" is const.

            // find the bounds.
            auto range_begin = ::std::lower_bound(db_begin, db_end, *query_begin, less);
            auto query_last = query_begin; ::std::advance(query_last, dist_query - 1);  // last real entry.
            auto range_end = ::std::upper_bound(range_begin, db_end, *query_last, less);

            return ::std::make_pair(range_begin, range_end);
          }

          // assumes that container is sorted. and exact overlap region is provided.
        template <class DBIter, class QueryIter, class OutputIter, class Operator, typename Predicate = Identity>
        static size_t process(DBIter range_begin, DBIter range_end,
        		QueryIter query_begin, QueryIter query_end,
        		OutputIter output, Operator const & op,
        		bool sorted_query = false, Predicate const &pred = Predicate()) {

            // no matches in container.
            if (range_begin == range_end) return 0;
            if (query_begin == query_end) return 0;  // no input

            auto dist_range = ::std::distance(range_begin, range_end);
            auto dist_query = ::std::distance(query_begin, query_end);

            //if (!sorted_target) sort_ascending(range_begin, range_end);  range_begin and range_end often are const iterators.
            if (!sorted_query) sort_ascending(query_begin, query_end);

          auto el_end = range_begin;
          size_t count = 0;

          if ((static_cast<double>(dist_query) * ::std::log2(dist_range)) > dist_range) {  // based on number of input and search source, choose a method to search.

            // iterate through the input and search in output -
            for (auto it = query_begin; it != query_end;) {
              auto v = *it;

              if (::std::is_same<Predicate, Identity>::value || pred(v))
                count += op.template operator()<true>(range_begin, el_end, range_end, v, output);

              // compiler optimize out the conditional.
              if (unique) it = upper_bound<true>(it, query_end, v);
              else ++it;

            }
          } else {
            // use logarithmic search

            // iterate through the input and search in output -
            for (auto it = query_begin; it != query_end;) {
              auto v = *it;

              if (::std::is_same<Predicate, Identity>::value || pred(v))
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

      mutable size_t key_multiplicity;

      // communication stuff...
      MPI_Comm comm;
      int comm_size;
      int comm_rank;

      // defined Communicator as a friend
      friend Comm;

      /// reserve space.  n is the local container size.  this allows different processes to individually adjust its own size.
      void local_reserve( size_type n) {
        c.reserve(n);
      }

      template <bool linear, class DBIter, class Query = typename ::std::iterator_traits<DBIter>::value_type>
      static inline DBIter lower_bound(DBIter b, DBIter e, Query& v) {
          // compiler choose one.
          if (linear) while ((b != e) && less(*b, v)) ++b;
          else b = ::std::lower_bound(b, e, v, less);
          return b;
      }


      template <bool linear, class DBIter, class Query = typename ::std::iterator_traits<DBIter>::value_type>
      static inline DBIter upper_bound(DBIter b, DBIter e,  Query& v) {
          // compiler choose one.
          if (linear) while ((b != e) && !less(v, *b)) ++b;
          else b = ::std::upper_bound(b, e, v, less);
          return b;
      }


      template <class DBIter>
      static void sort_ascending(DBIter b, DBIter e) {
          if (!::std::is_sorted(b, e, less)) {
            ::std::sort(b, e, less);
          }
      }

      /// rehash the local container.  n is the local container size.  this allows different processes to individually adjust its own size.
      void local_sort() {
        if (sorted) return;
        sort_ascending(c.begin(), c.end());
        sorted = true;
      }

      void assert_sorted_locally() const { if (!sorted) throw ::std::logic_error("local_sort needed to be called to sort local vector"); }
      bool is_sorted_locally() const { return sorted; }

      struct LocalCount {
          template<bool linear, class DBIter, typename Query, class OutputIter>
          size_t operator()(DBIter &range_begin, DBIter &el_end, DBIter const &range_end, Query const &v, OutputIter &output) const {
            // TODO: LINEAR SEARCH, O(n).  vs BINARY SEARCH, O(mlogn)
            range_begin = lower_bound<linear>(el_end, range_end, v);  // range_begin at equal or greater than v.
            el_end = upper_bound<linear>(range_begin, range_end, v);  // el_end at greater than v.
            // difference between the 2 iterators is the part that's equal.

            // add the output entry.
            *output = ::std::move(::std::make_pair(v, ::std::distance(range_begin, el_end)));
            ++output;
            return 1;
          }
      } local_count;


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
      } local_erase;




      /// clears the sorted_map
      void local_clear() noexcept {
          c.clear();
      }

      ///  keep the unique keys in the input.   output is sorted.  equal operator forces comparison to Key
      template <typename V>
      void retain_unique(::std::vector< V >& input, bool sorted_input = false) const {
        if (input.size() == 0) return;
        if (!sorted_input) sort_ascending(input.begin(), input.end());
        auto end = ::std::unique(input.begin(), input.end(), this->equal);
        input.erase(end, input.end());
      }

      sorted_map_base(MPI_Comm _comm, int _comm_size) : key_to_rank(_comm_size), sorted(true), key_multiplicity(1), comm(_comm), comm_size(_comm_size) {
        MPI_Comm_rank(comm, &comm_rank);
      }

      /**
       * @brief find elements with the specified keys in the distributed sorted_multimap.
       * @param first
       * @param last
       */
      template <class LocalFind, class Predicate = Identity>
      ::std::vector<::std::pair<Key, T> > find(LocalFind const & local_find, ::std::vector<Key>& keys, bool sorted_input = false, Predicate const& pred = Predicate()) const {
        TIMER_INIT(find);

        TIMER_START(find);
        ::std::vector<::std::pair<Key, T> > results;
        // even if count is 0, still need to participate in mpi calls.  if (keys.size() == 0) return results;
        back_emplace_iterator<::std::vector<::std::pair<Key, T> > > emplace_iter(results);
        TIMER_END(find, "begin", keys.size());

        this->assert_sorted_locally();

        // keep unique keys
        TIMER_START(find);
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
          results.reserve(keys.size() * this->key_multiplicity);  // 1 result per key.


          TIMER_END(find, "reserve", (keys.size() * this->key_multiplicity));

          TIMER_START(find);
          auto start = keys.begin();
          auto end = start;
          for (int i = 0; i < this->comm_size; ++i) {
            ::std::advance(end, recv_counts[i]);

            // work on query from process i.
            auto overlap = Intersect<false>::intersect(this->c.begin(), this->c.end(), start, end, true);

            // within start-end, values are unique, so don't need to set unique to true.
            send_counts[i] = Intersect<false>::process(overlap.first, overlap.second, start, end, emplace_iter, local_find, true, pred);

            if (this->comm_rank == 0) DEBUGF("R %d added %d results for %d queries for process %d\n", this->comm_rank, send_counts[i], recv_counts[i], i);

            start = end;
          }
          TIMER_END(find, "local_find", results.size());

          // send back using the constructed recv count
          TIMER_START(find);
          mxx::all2all(results, send_counts, this->comm);
          TIMER_END(find, "a2a2", results.size());


        } else {

          TIMER_START(find);
          results.reserve(keys.size() * this->key_multiplicity);  // 1 result per key.
          TIMER_END(find, "reserve", (keys.size() * this->key_multiplicity));


          TIMER_START(find);
          auto overlap = Intersect<false>::intersect(this->c.begin(), this->c.end(), keys.begin(), keys.end(), true);

          // within start-end, values are unique, so don't need to set unique to true.
          Intersect<false>::process(overlap.first, overlap.second, keys.begin(), keys.end(), emplace_iter, local_find, true, pred);

          TIMER_END(find, "local_find", results.size());

        }

        TIMER_REPORT_MPI(find, this->comm_rank, this->comm);

        return results;
      }

      template <class LocalFind, class Predicate = Identity>
      ::std::vector<::std::pair<Key, T> > find(LocalFind const & local_find, Predicate const & pred = Predicate()) const {
          ::std::vector<::std::pair<Key, T> > results;
          back_emplace_iterator<::std::vector<::std::pair<Key, T> > > emplace_iter(results);

          auto keys = this->keys();
          results.reserve(keys.size() * this->key_multiplicity);  // 1 result per key.

          // within start-end, values are unique, so don't need to set unique to true.
          Intersect<false>::process(this->c.begin(), this->c.end(), keys.begin(), keys.end(), emplace_iter, local_find, true, pred);

          if (this->comm_size > 1) MPI_Barrier(this->comm);
          return results;
      }


    public:

      virtual ~sorted_map_base() {};

      /// returns the local storage.  please use sparingly.
      local_container_type& get_local_container() { return c; }

      /// update the multiplicity.  only multimap needs to do this.
      virtual size_t update_multiplicity() {
        this->rehash();
        return key_multiplicity;
      }

      /// convert the map to a vector.
      std::vector<std::pair<Key, T> > to_vector() const {
        std::vector<std::pair<Key, T> > result(c.begin(), c.end());
        return result;
      }

      /// convert the map to a vector
      void to_vector(std::vector<std::pair<Key, T> > & result) const {
        result.clear();
        if (c.empty()) return;

        result.assign(c.begin(), c.end());
      }

      /// extract the keys of a map.
      std::vector<Key> keys() const {
        std::vector<Key> result;
        this->keys(result);
        return result;
      }

      /// extract the unique keys of a map.
      void keys(std::vector<Key> & result) const {
        result.clear();
        if (c.empty()) return;

        // copy the keys
        for (auto it = c.begin(), end = c.end(); it != end; ++it) {
          result.emplace_back(it->first);
        }

        // and then find unique.
        retain_unique(result, this->sorted);
      }

      // note that for each method, there is a local version of the operartion.
      // this is for use by the asynchronous version of communicator as callback for any messages received.
      /// check if empty.
      bool local_empty() const noexcept {
        return c.empty();
      }

      /// get size of local container
      size_type local_size() const noexcept {
        return c.size();
      }

      /// check if empty.
      bool empty() const noexcept {
        if (comm_size == 1)
          return local_empty();
        else // all reduce
          return mxx::test_all(local_empty(), comm);
      }

      /// get size of distributed container
      size_type size() const noexcept {
        size_type s = local_size();
        if (comm_size == 1)
          return s;
        else
          return mxx::allreduce(s, comm);
      }


      /// reserve space.  n is the local container size.  this allows different processes to individually adjust its own size.
      void reserve( size_type n) {
        // direct reserve + barrier
        local_reserve(n);
        if (comm_size > 1) MPI_Barrier(comm);
      }

      // default is for multimap
      virtual void rehash() = 0;

//      /// rehash the local container.  n is the local container size.  this allows different processes to individually adjust its own size.
//      virtual void rehash( size_type n = 0) {
//        // direct rehash + barrier
//
//
//        this->key_to_rank.map.clear();
//        if (this->comm_size > 1) {
//
//        	// at this point, local data is sorted and unique locally for map and derived, and not sorted nor unique for multimap
//        	// content sizes may be uneven
//
//        	// first balance
//            printf("R %d container size b4 block decompose: %lu\n", this->comm_rank, c.size());
//            auto d = ::mxx::stable_block_decompose(c, this->comm);
//            printf("R %d container size after block decompose: %lu\n", this->comm_rank, d.size());
//
//        	// then sort
//            ::mxx::sort(d.begin(), d.end(), this->less, this->comm, false);
//
//            // get new splitter
//            // next compute the splitters.
//            if ((this->comm_rank < (this->comm_size - 1)) && (d.size() > 0)) {
//              // only send for the first p-1 proc, and only if they have a kmer to split with.
//              this->key_to_rank.map.emplace_back(d.back().first, this->comm_rank);
//            }
//            this->key_to_rank.map = ::mxx::allgatherv(this->key_to_rank.map, this->comm);
//
//            for (int i = 0; i < this->key_to_rank.map.size(); ++i) {
//          	  printf("R %d key to rank %s -> %d\n", this->comm_rank, this->key_to_rank.map[i].first.toAlphabetString().c_str(), this->key_to_rank.map[i].second);
//            }
//
//            // redistribute.
//            c.clear();
//            ::mxx2::msgs_all2all(d, this->key_to_rank, this->comm);
//
//            c.assign(d.begin(), d.end());
//            printf("R %d container size after redistribute: %lu\n", this->comm_rank, c.size());
//
//            this->sorted = true;
//
//        } else {
//			// locally sort
//			local_rehash(n);
//        }
//
//      }


      /**
       * @brief count elements with the specified keys in the distributed sorted_multimap.
       * @param first
       * @param last
       */
      template <typename Predicate = Identity>
      ::std::vector<::std::pair<Key, size_type> > count(::std::vector<Key>& keys, bool sorted_input = false, Predicate const & pred = Predicate()) const {

        ::std::vector<::std::pair<Key, size_type> > results;
        // even if count is 0, still need to participate in mpi calls.  if (keys.size() == 0) return results;
        back_emplace_iterator<::std::vector<::std::pair<Key, size_type> > > emplace_iter(results);

        this->assert_sorted_locally();

        if (comm_size > 1) {
          // keep unique keys
          retain_unique(keys, sorted_input);

          // distribute (communication part)
          std::vector<int> recv_counts = mxx2::msgs_all2all(keys, this->key_to_rank, comm);

          // local find. memory utilization a potential problem.
          // do for each src proc one at a time.

          results.reserve(keys.size());                   // TODO:  should estimate coverage.


          auto start = keys.begin();
          auto end = start;
          for (int i = 0; i < comm_size; ++i) {
            ::std::advance(end, recv_counts[i]);

            // work on query from process i.
            auto overlap = Intersect<false>::intersect(this->c.begin(), this->c.end(), start, end, true);

            // within start-end, values are unique, so don't need to set unique to true.
            Intersect<false>::process(overlap.first, overlap.second, start, end, emplace_iter, local_count, true, pred);

            if (comm_rank == 0) DEBUGF("R %d added %d results for %d queries for process %d\n", comm_rank, send_counts[i], recv_counts[i], i);

            start = end;
          }

          // send back using the constructed recv count
          mxx::all2all(results, recv_counts, comm);

        } else {
          results.reserve(keys.size());
          // work on query from process i.
          auto overlap = Intersect<true>::intersect(this->c.begin(), this->c.end(), keys.begin(), keys.end(), sorted_input);

          // within key, values may not be unique,
          Intersect<true>::process(overlap.first, overlap.second, keys.begin(), keys.end(), emplace_iter, local_count, true, pred);
        }

        return results;

      }


      template <typename Predicate = Identity>
      ::std::vector<::std::pair<Key, size_type> > count(Predicate const & pred = Predicate()) const {
        ::std::vector<::std::pair<Key, size_type> > results;
        back_emplace_iterator<::std::vector<::std::pair<Key, size_type> > > emplace_iter(results);

        auto keys = this->keys();
        results.reserve(keys.size());

        // keys already unique
        Intersect<false>::process(c.begin(), c.end(), keys.begin(), keys.end(), emplace_iter, local_count, true, pred);

        if (comm_size > 1) MPI_Barrier(comm);
        return results;
      }


      /**
       * @brief insert new elements in the distributed sorted_multimap.  example use: stop inserting if more than x entries.
       * @param src_begin
       * @param src_end
       */
      template <class InputIter, class Predicate = Identity>
      size_t insert(InputIter src_begin, InputIter src_end, bool sorted_input = false, Predicate const &pred = Predicate()) {
          if (src_begin == src_end) return 0;

          this->sorted = false; this->balanced = false; this->globally_sorted = false;
          back_emplace_iterator<local_container_type> emplace_iter(c);

          if (::std::is_same<Predicate, Identity>::value) ::std::copy(src_begin, src_end, emplace_iter);
          else ::std::copy_if(src_begin, src_end, emplace_iter, pred);

          return ::std::distance(src_begin, src_end);
      }

      /**
       * @brief insert new elements in the distributed sorted_multimap.  example use: stop inserting if more than x entries.
       * @param first
       * @param last
       */
      template <class Predicate = Identity>
      size_t insert(::std::vector<::std::pair<Key, T> > &input, bool sorted_input = false, Predicate const &pred = Predicate()) {
    	  return insert(input.begin(), input.end(), sorted_input, pred);
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
        if (comm_size > 1) {
          // remove duplicates
          retain_unique(keys, si);

          // redistribute keys
          mxx2::msgs_all2all(keys, this->key_to_rank, comm);

          si = false;
        }

        this->local_sort();

        //== now call local remove.
        // first get the range of intersection
        auto overlap = Intersect<true>::intersect(this->c.begin(), this->c.end(), keys.begin(), keys.end(), si);

        // do the work.  unique = true
        size_t kept = Intersect<true>::process(overlap.first, overlap.second, keys.begin(), keys.end(), overlap.first, local_erase, true, pred);

        // move the last part.
        auto end = overlap.first;
        ::std::advance(end, kept);
        end = ::std::move(overlap.second, c.end(), end);
        size_t before = c.size();
        c.erase(end, c.end());

        this->balanced = false;

        return before - c.size();
      }

      template <typename Predicate = Identity>
      size_t erase(Predicate const & pred = Predicate()) {
        size_t before = c.size();
        if (!::std::is_same<Predicate, Identity>::value) {

          local_sort();

          auto end = ::std::partition(c.begin(), c.end(), pred);
          c.erase(end, c.end());

          this->balanced = false;
        }

        if (comm_size > 1) MPI_Barrier(comm);

        return before - c.size();
      }

      // update done via erase/insert.

      /// clears the sorted_map
      void clear() noexcept {
        // clear + barrier.
        local_clear();
        if (comm_size > 1) MPI_Barrier(comm);
        this->sorted = true; this->balanced = true; this->globally_sorted = true;
      }


  };


  template<typename Key, typename T,
      template <typename, typename> class Container,
      class Comm,
      template <typename> class KeyTransform,
      class Less,
      class Equal,
      class Alloc>
  KeyTransform<Key> sorted_map_base<Key, T, Container, Comm, KeyTransform, Less, Equal, Alloc>::trans;

  template<typename Key, typename T,
      template <typename, typename> class Container,
      class Comm,
      template <typename> class KeyTransform,
      class Less,
      class Equal,
      class Alloc>
  typename sorted_map_base<Key, T, Container, Comm, KeyTransform, Less, Equal, Alloc>::TransformedLess
  	  sorted_map_base<Key, T, Container, Comm, KeyTransform, Less, Equal, Alloc>::less;

  template<typename Key, typename T,
      template <typename, typename> class Container,
      class Comm,
      template <typename> class KeyTransform,
      class Less,
      class Equal,
      class Alloc>
  typename sorted_map_base<Key, T, Container, Comm, KeyTransform, Less, Equal, Alloc>::TransformedEqual
  	  sorted_map_base<Key, T, Container, Comm, KeyTransform, Less, Equal, Alloc>::equal;

  template<typename Key, typename T,
      template <typename, typename> class Container,
      class Comm,
      template <typename> class KeyTransform,
      class Less,
      class Equal,
      class Alloc>
  typename sorted_map_base<Key, T, Container, Comm, KeyTransform, Less, Equal, Alloc>::TransformedGreater
  	  sorted_map_base<Key, T, Container, Comm, KeyTransform, Less, Equal, Alloc>::greater;


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
            }

            return 1;
          }
      } local_find;

//      /**
//       * @brief insert new elements in the distributed sorted_map.  result is not sorted.  input must not hvae duplicates;
//       * @param first
//       * @param last
//       */
//      template <class InputIterator>
//      void local_insert(InputIterator first, InputIterator last, bool sorted_input = false) {
//          if (first == last) return;
//
//          // sort both
//          this->local_rehash();
//
//          // also ensure that input is unique
//          if (!sorted_input) this->sort_ascending(first, last);
//
//          // walk through.  if exiting entry, remove from input.
//          // now walk through the input to count stuff.
//          auto t_start = ::std::lower_bound(this->c.begin(), this->c.end(), *first, this->less);
//          auto last2 = first; ::std::advance(last2, ::std::distance(first, last) - 1);
//          auto t_end = ::std::upper_bound(t_start, this->c.end(), *last2, this->less);
//
//          // go through the input and search the container.
//
//          // modify input so as to not invalidate the container's iterators, t_start, t_end;
//            // iterate through the input and search in output,  moving items up as we go.
//            auto it = first;
//            auto end = first;
//            for (it = first; it != last; ) {  // from large to small
//              auto v = *it;
//
//              // find entry in container.
//              t_start = this->template lower_bound<true>(t_start, t_end, v);
//
//              // if we found it, then t_start != t_end, and equal(v, *t_start)
//              // else
//              if ((t_start == t_end) || (!this->equal(v, *t_start))) {
//                // not matched.  copy the value, move pos up by 1.
//                if (end != it) *end = v;
//                ++end;
//              } // else matched.  so skip it.
//
//              it = this->upper_bound<true>(it, last, v);
//            }
//
//          // insert at the end.
//          if (first != end) {
//            this->c.insert(this->c.end(), first, end);
//            // after, not sorted
//            this->sorted = false;
//          }
//
//      }
//
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
//          if (!sorted_input) this->sort_ascending(first, last);
//
//          // walk through.  if exiting entry, remove from input.
//          // now walk through the input to count stuff.
//          auto t_start = ::std::lower_bound(this->c.begin(), this->c.end(), *first, this->less);
//          auto last2 = first; ::std::advance(last2, ::std::distance(first, last) - 1);
//          auto t_end = ::std::upper_bound(t_start, this->c.end(), *last2, this->less);
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

      virtual void local_reduction(std::vector<::std::pair<Key, T> > &input, bool sorted_input = false) {
        if (!sorted_input) this->sort_ascending(input.begin(), input.end());

        auto end = ::std::unique(input.begin(), input.end(), this->equal);
        input.erase(end, input.end());
      }

    public:

      sorted_map(MPI_Comm _comm, int _comm_size) : Base(_comm, _comm_size) {}

      virtual ~sorted_map() {};

      template <class Predicate = Identity>
      ::std::vector<::std::pair<Key, T> > find(::std::vector<Key>& keys, bool sorted_input = false, Predicate const& pred = Predicate()) const {
    	  return Base::find(local_find, keys, sorted_input, pred);
      }

      template <class Predicate = Identity>
      ::std::vector<::std::pair<Key, T> > find(Predicate const& pred = Predicate()) const {
    	  return Base::find(local_find, pred);
      }


      // default is for map
      virtual void rehash() {
printf("c size before: %lu\n", this->c.size());
        if (this->comm_size > 1) {
          // first balance

          if (!this->balanced) {
        	  this->c = ::std::move(::mxx::stable_block_decompose(this->c, this->comm));
          }

          // sort if needed
          if (!this->globally_sorted) {
            // kway merge / sort
            this->sort_ascending(this->c.begin(), this->c.end());

            // local unique
            this->local_reduction(this->c, true);

            // sample
            mxx::datatype<::std::pair<Key, T> > dt;
            MPI_Datatype mpi_dt = dt.type();
            this->key_to_rank.map = ::mxx::impl::sample_arbit_decomp(this->c.begin(), this->c.end(), this->less, this->comm_size -1, this->comm, mpi_dt);

            // distribute
            ::mxx2::msgs_all2all(this->c, this->key_to_rank, this->comm);

            // kway merge / sort
            this->sort_ascending(this->c.begin(), this->c.end());

            // local unique
            this->local_reduction(this->c, true);

            // rebalance
            this->c = ::std::move(::mxx::stable_block_decompose(this->c, this->comm));
          }

          // get new pivots
          // next compute the splitters.
          this->key_to_rank.map.clear();
          if ((this->comm_rank < (this->comm_size - 1)) && (this->c.size() > 0)) {
            // only send for the first p-1 proc, and only if they have a kmer to split with.
            this->key_to_rank.map.emplace_back(this->c.back().first, this->comm_rank);
          }
          this->key_to_rank.map = ::mxx::allgatherv(this->key_to_rank.map, this->comm);

          for (int i = 0; i < this->key_to_rank.map.size(); ++i) {
            printf("R %d key to rank %s -> %d\n", this->comm_rank, this->key_to_rank.map[i].first.toAlphabetString().c_str(), this->key_to_rank.map[i].second);
          }

          // no need to redistribute - each entry is unique so nothing is going to span processor boundaries.

        } else {
          // local unique
          this->retain_unique(this->c, this->sorted);
        }
        this->sorted = true; this->balanced = true; this->globally_sorted = true;
        printf("c size after: %lu\n", this->c.size());

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
      } local_find;


    public:


      sorted_multimap(MPI_Comm _comm, int _comm_size) : Base(_comm, _comm_size) {
        this->key_multiplicity = 50;
      }

      virtual ~sorted_multimap() {}


      // default is for multimap
      virtual void rehash() {
    	  if (this->comm_size > 1) {
			  // first balance

			  if (!this->balanced) {
				  this->c = std::move(::mxx::stable_block_decompose(this->c, this->comm));
			  }

			  // sort if needed
			  if (!this->globally_sorted) ::mxx::sort(this->c.begin(), this->c.end(), this->less, this->comm, false);

			  // get new pivots
			  // next compute the splitters.
			  if ((this->comm_rank < (this->comm_size - 1)) && (this->c.size() > 0)) {
				// only send for the first p-1 proc, and only if they have a kmer to split with.
				this->key_to_rank.map.emplace_back(this->c.back().first, this->comm_rank);
			  }
			  this->key_to_rank.map = ::mxx::allgatherv(this->key_to_rank.map, this->comm);

			  for (int i = 0; i < this->key_to_rank.map.size(); ++i) {
				  printf("R %d key to rank %s -> %d\n", this->comm_rank, this->key_to_rank.map[i].first.toAlphabetString().c_str(), this->key_to_rank.map[i].second);
			  }


//			  mxx::datatype<::std::pair<Key, T> > dt;
//			  MPI_Datatype mpi_dt = dt.type();
//			  this->key_to_rank.map = ::mxx::sample_block_decomp(d.begin(), d.end(), this->less, this->comm_size - 1, this->comm, mpi_dt);

			  // redistribute
			  ::mxx2::msgs_all2all(this->c, this->key_to_rank, this->comm);

    	  } else {
    		  this->local_sort();
    	  }
          this->sorted = true; this->balanced = true; this->globally_sorted = true;

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
        printf("%lu elements, %lu unique, key multiplicity = %lu\n", this->c.size(), uniq_count, this->key_multiplicity);


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
      ::std::vector<::std::pair<Key, T> > find(::std::vector<Key>& keys, bool sorted_input = false, Predicate const& pred = Predicate()) const {
    	  return Base::find(local_find, keys, sorted_input, pred);
      }

      template <class Predicate = Identity>
      ::std::vector<::std::pair<Key, T> > find(Predicate const& pred = Predicate()) const {
    	  return Base::find(local_find, pred);
      }

//      size_t count_unique(::std::vector<::std::pair<Key, T> > const & input) const {
//        // alternative approach to get number of unique keys is to use an unordered_set.  this will take more memory but probably will be faster than sort for large buckets (high repeats).
//        ::std::unordered_set<Key, typename Base::TransformedHash, typename Base::template TransformedComp<Equal> > unique_set(this->c.size());
//        for (auto it = input.begin(), max = input.end(); it != max; ++it) {
//          unique_set.insert(it->first);
//        }
//        printf("r %d: %lu elements, %lu unique\n", this->comm_rank, input.size(), unique_set.size());
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
      using Base = sorted_map<Key, T, Comm, KeyTransform, Less, Equal, Alloc>;

      static_assert(::std::is_arithmetic<T>::value, "mapped type has to be arithmetic");

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



//      /**
//       * @brief insert new elements in the distributed sorted_map.  result is not sorted.  input should not have duplicates
//       * @note  input CAN have duplicates
//       * @param first
//       * @param last
//       */
//      template <class InputIterator>
//      void local_insert(InputIterator first, InputIterator last, bool sorted_input = false) {
//          if (first == last) return;
//
//          this->local_rehash();
//
//          // has to be locally reduced first.
//          auto newlast = local_reduction(first, last, sorted_input);
//
//          // walk through.  if exiting entry, remove from input.
//          // now walk through the input to count stuff.
//          auto t_start = ::std::lower_bound(this->c.begin(), this->c.end(), *first, this->less);
//          auto t_end = ::std::upper_bound(t_start, this->c.end(), *(newlast - 1), this->less);
//
//          // go through the input and search the container.
//          // iterate through the input and search in output,  moving items up as we go.
//          auto end = first;
//          for (auto it = first; it != newlast; ++it) {  // walk though all input entries
//              auto v = *it;
//
//                // find entry in container.
//                t_start = this->lower_bound<true>(t_start, t_end, v);
//
//                if ((t_start == t_end) || (!this->equal(v, *t_start))) {
//                  // not matched.  copy the value, move pos up by 1.
//                  if (end != it) *end = v;
//                  ++end;
//                } else {
//                  // matched.  so need to reduce
//                  t_start->second = r(t_start->second, v.second);
//                }
//            }
//          // insert at the end.
//          if (first != end) {
//            this->c.insert(this->c.end(), first, end);
//            this->sorted = false;
//          }
//      }

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
//          auto t_start = ::std::lower_bound(this->c.begin(), this->c.end(), *first, this->less);
//          auto t_end = ::std::upper_bound(t_start, this->c.end(), *(newlast - 1), this->less);
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
          if (!sorted_input) this->sort_ascending(b, e);

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


      reduction_sorted_map(MPI_Comm _comm, int _comm_size) : Base(_comm, _comm_size) {}

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
      using Base = reduction_sorted_map<Key, T, Comm, KeyTransform, Less, ::std::plus<T>, Equal, Alloc>;

      static_assert(::std::is_integral<T>::value, "count type has to be integral");

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

      // convert key to a pair.
      ::std::vector<::std::pair<Key, T> > castToPair(::std::vector<Key>& input, bool sorted_input = false) {

        ::std::vector<::std::pair<Key, T> > output;
        output.reserve(input.size());

        for (int i = 0; i < input.size(); ++i) {
        	output.emplace_back(input[i], 1);
        }

//        if (input.size() == 0) return output;
//
//        output.reserve(input.size());
//
//
//        if (!sorted_input) this->sort_ascending(input.begin(), input.end());
//
//        // then do reduction
//        auto curr = input.begin();
//        auto v = *curr;
//        output.emplace_back(v, 1);
//        ++curr;
//        while (curr != input.end()) {
//          if (this->equal(v, *curr))  // if same, do reduction
//            ++(output.back().second);
//          else {  // else reset first.
//            v = *curr;
//            output.emplace_back(v, 1);
//          }
//          ++curr;  // increment second.
//        }
        return output;
      }


    public:
      counting_sorted_map(MPI_Comm _comm, int _comm_size) : Base(_comm, _comm_size) {}

      virtual ~counting_sorted_map() {};


      /**
       * @brief insert new elements in the distributed sorted_multimap.
       * @param first
       * @param last
       */
      void insert(std::vector< Key >& input, bool sorted_input = false) {
        // even if count is 0, still need to participate in mpi calls.  if (input.size() == 0) return;
        TIMER_INIT(count_insert);

        TIMER_START(count_insert);
        TIMER_END(count_insert, "start", input.size());

        // distribute
        TIMER_START(count_insert);

        // first remove duplicates.  sort, then get unique, finally remove the rest.  may not be needed
        auto temp = this->castToPair(input, sorted_input);
        printf("r %d count %lu unique %lu\n", this->comm_rank, input.size(), temp.size());
        TIMER_END(count_insert, "reduc1", temp.size());

        bool si = true;

//        if (this->comm_size > 1) {
//
//
//          // distribute
//          TIMER_START(count_insert);
//
//          // communication part
//          mxx2::msgs_all2all(temp, this->key_to_rank, this->comm);
//          TIMER_END(count_insert, "a2a", temp.size());
//
//          si = false;
//        }

          // distribute
          TIMER_START(count_insert);

          // local compute part.  called by the communicator.
          this->Base::insert(temp, si);
          TIMER_END(count_insert, "insert", this->c.size());

        // distribute
        TIMER_REPORT_MPI(count_insert, this->comm_rank, this->comm);

      }



  };



} /* namespace dsc */


#endif // BLISS_DISTRIBUTED_SORTED_MAP_HPP
