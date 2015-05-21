/**
 * @file    distributed_map.hpp
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

#ifndef BLISS_DISTRIBUTED_MAP_HPP
#define BLISS_DISTRIBUTED_MAP_HPP


#include <unordered_map>  // local storage hash table  // for multimap
#include <unordered_set>  // local storage hash table  // for multimap
#include <utility> 			  // for std::pair

//#include <sparsehash/dense_hash_map>  // local storage hash table, google dense hash map.
#include <functional> 		// for std::function and std::hash
#include <algorithm> 		// for sort, stable_sort, unique, is_sorted
#include <iterator>  // advance, distance

#include <cstdint>  // for uint8, etc.

#include "mxx/collective.hpp"
#include "mxx/reduction.hpp"
#include "io/mpi_utils.hpp"
#include "utils/timer.hpp"  // for timing.
#include "utils/logging.h"

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
  >
  class unordered_map_base {
    protected:
      template <typename K = Key>
      struct TransformedHash {
        KeyTransform<K> trans;
        Hash<K, false> hash;
        inline uint64_t operator()(K const& k) const {
          return hash(trans(k));
        }
      };

      template <typename K = Key>
      struct TransformedEqual {
          KeyTransform<K> trans;
          Equal equal;
          inline bool operator()(K const & x, K const & y) const {
            return equal(trans(x), trans(y));
          }
      };

      struct KeyToRank {
          Hash<Key, true> proc_hash;
          KeyTransform<Key> trans;

          int p;

          // 2x comm size to allow more even distribution?
          KeyToRank(int comm_size) : proc_hash(ceilLog2(2 * comm_size)), p(comm_size) {};

          inline int operator()(Key const & x) const {
//            printf("KeyToRank operator. commsize %d  key.  hashed to %d, mapped to proc %d \n", p, proc_hash(trans(x)), proc_hash(trans(x)) % p);
            return proc_hash(trans(x)) % p;
          }
          inline int operator()(::std::pair<Key, T> const & x) const {
//            printf("KeyToRank operator.  pair\n");
            return this->operator()(x.first);
          }
          inline int operator()(::std::pair<const Key, T> const & x) const {
            return this->operator()(x.first);
          }
      } key_to_rank;


    public:
      using local_container_type = Container<Key, T, TransformedHash<Key>, TransformedEqual<Key>, Alloc>;

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
      using local_iterator        = typename local_container_type::local_iterator;
      using const_local_iterator  = typename local_container_type::const_local_iterator;
      using size_type             = typename local_container_type::size_type;
      using difference_type       = typename local_container_type::difference_type;

    protected:
      local_container_type c;

      mutable size_t key_multiplicity;

      // communication stuff...
      MPI_Comm comm;
      int comm_size;
      int comm_rank;

      // defined Communicator as a friend
      friend Comm;

      KeyTransform<Key> trans;

      /// reserve space.  n is the local container size.  this allows different processes to individually adjust its own size.
      void local_reserve( size_type n) {
        c.reserve(n);
      }

      /// rehash the local container.  n is the local container size.  this allows different processes to individually adjust its own size.
      void local_rehash( size_type n) {
        c.rehash(n);
      }


      /**
       * @brief count elements with the specified keys in the distributed unordered_multimap.
       * @param first
       * @param last
       */
      template <class InputIterator>
      size_t local_count(InputIterator first, InputIterator last, ::std::vector<::std::pair<Key, size_type> > &output) const {
          if (first == last) return 0;

          size_t before = output.size();  // before size.
          for (auto it = first; it != last; ++it) {
             output.emplace_back(*it, c.count(*it));
          }
          return output.size() - before;
      }

      /**
       * @brief predicate version of local_count.  example use: count items with counts of 1 only.
       * @param first
       * @param last
       * @param output
       * @param pred        Predicate should account for key transform.
       * @return
       */
      template<class InputIterator, class Predicate>
      size_t local_count_if(InputIterator first, InputIterator last, ::std::vector<::std::pair<Key, size_type> > &output, Predicate const & pred) const {
          if (first == last) return 0;

          size_t before = output.size();  // before size.
          for (auto it = first; it != last; ++it) {
             if (pred(trans(*it))) output.emplace_back(*it, c.count(*it));
          }
          return output.size() - before;
      }


      /**
       * @brief insert new elements in the distributed unordered_multimap.
       * @param first
       * @param last
       */
      template <class InputIterator>
      void local_insert(InputIterator first, InputIterator last) {
          if (first == last) return;

          c.insert(first, last);
      }

      /**
       * @brief insert new elements in the distributed unordered_multimap.  example use: stop inserting if more than x entries.
       * @param first
       * @param last
       */
      template <class InputIterator, class Predicate>
      void local_insert_if(InputIterator first, InputIterator last, Predicate const &pred) {
          if (first == last) return;

          for (auto it = first; it != last; ++it) {
            if (pred(trans(*it))) c.insert(*it);
          }
      }


      /**
       * @brief erase elements with the specified keys in the distributed unordered_multimap.
       * @param first
       * @param last
       */
      template <class InputIterator>
      void local_erase(InputIterator first, InputIterator last) {
        if (first == last) return;

        for (auto it = first; it != last; ++it) {
          c.erase(*it);
        }      
      }

      /**
       * @brief erase elements with the specified keys in the distributed unordered_multimap.  example use:  remove all entries with low frequency
       * @param first
       * @param last
       */
      template <class InputIterator, class Predicate>
      void local_erase_if(InputIterator first, InputIterator last, Predicate const &pred) {
        if (first == last) return;

        for (auto it = first; it != last; ++it) {
          if (pred(trans(*it))) c.erase(*it);
        }
      }



      /// clears the unordered_map
      void local_clear() noexcept {
          c.clear();
      }

      void retain_unique_keys(::std::vector< Key >& input) const {
        if (input.size() == 0) return;

        // sorting is slow.
//        if (! ::std::is_sorted(input.begin(), input.end(), key_less_op))
//          ::std::stable_sort(input.begin(), input.end(), key_less_op);  // using stable sort to maintain same behavior as map - keeping first inserted.
//        auto end = ::std::unique(input.begin(), input.end(), key_equal_op);
//        input.erase(end, input.end());

        // this function primarily is used to reduce communication volume while keeping the behavior of unordered_map - i.e. first entry is kept.
        // simplest approach is to use the same type of data container to get the behavior,
        ::std::unordered_set<Key, TransformedHash<Key>, TransformedEqual<Key> > temp;
        temp.reserve(input.size());
        temp.insert(input.begin(), input.end());  // this process removes duplicates
        input.assign(temp.begin(), temp.end());   // copy back into original input.
      }


      unordered_map_base(MPI_Comm _comm, int _comm_size) : key_to_rank(_comm_size), key_multiplicity(1), comm(_comm), comm_size(_comm_size) {
        MPI_Comm_rank(comm, &comm_rank);
      }


    public:

      virtual ~unordered_map_base() {};

      /// returns the local storage.  please use sparingly.
      local_container_type& get_local_container() { return c; }

      /// update the multiplicity.  only multimap needs to do this.
      virtual size_t update_multiplicity() const { return key_multiplicity; }

      /// convert the map to a vector.
      std::vector<std::pair<Key, T> > to_vector() const {
        std::vector<std::pair<Key, T> > result;
        this->to_vector(result);
        return result;
      }

      /// convert the map to a vector
      void to_vector(std::vector<std::pair<Key, T> > & result) const {
        result.clear();
        if (c.empty()) return;
        result.reserve(c.size());
        result.insert(result.end(), c.begin(), c.end());
      }

      /// extract the keys of a map.
      std::vector<Key> keys() const {
        std::vector<Key> result;
        this->keys(result);
        return result;
      }

      /// extract the keys of a map.
      void keys(std::vector<Key> & result) const {
        result.clear();
        if (c.empty()) return;

        ::std::unordered_set<Key, TransformedHash<Key>, TransformedEqual<Key> > temp;
        temp.reserve(c.size());
        for (auto it = c.begin(), end = c.end(); it != end; ++it) {
          temp.insert(it->first);
        }
        result.assign(temp.begin(), temp.end());
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


      /// rehash the local container.  n is the local container size.  this allows different processes to individually adjust its own size.
      void rehash( size_type n) {
        // direct rehash + barrier
        local_rehash(n);
        if (comm_size > 1) MPI_Barrier(comm);
      }


      /**
       * @brief count elements with the specified keys in the distributed unordered_multimap.
       * @param first
       * @param last
       */
      ::std::vector<::std::pair<Key, size_type> > count(::std::vector<Key>& keys) const {
        TIMER_INIT(count);

        TIMER_START(count);
        ::std::vector<::std::pair<Key, size_type> > results;
        // even if count is 0, still need to participate in mpi calls.  if (keys.size() == 0) return results;
        TIMER_END(count, "begin", keys.size());

        TIMER_START(count);
        // keep unique keys
        retain_unique_keys(keys);
        TIMER_END(count, "uniq1", keys.size());


        if (comm_size > 1) {
          TIMER_START(count);
          // distribute (communication part)
          std::vector<int> recv_counts = mxx2::msgs_all2all(keys, this->key_to_rank, comm);
          TIMER_END(count, "a2a1", keys.size());

          // local count. memory utilization a potential problem.
          // do for each src proc one at a time.
          TIMER_START(count);
          results.reserve(keys.size());                   // TODO:  should estimate coverage.
          TIMER_END(count, "reserve", keys.size());

          TIMER_START(count);
          auto start = keys.begin();
          auto end = start;
          for (int i = 0; i < comm_size; ++i) {
            ::std::advance(end, recv_counts[i]);

            // work on query from process i.
            local_count( start, end, results );

            if (comm_rank == 0) DEBUGF("R %d added %d results for %d queries for process %d\n", comm_rank, send_counts[i], recv_counts[i], i);

            start = end;
          }
          TIMER_END(count, "local_count", results.size());

          // send back using the constructed recv count
          TIMER_START(count);
          mxx::all2all(results, recv_counts, comm);
          TIMER_END(count, "a2a2", results.size());
        } else {
          TIMER_START(count);
          local_count(keys.begin(), keys.end(), results);
          TIMER_END(count, "local_count", results.size());
        }

        TIMER_REPORT_MPI(count, this->comm_rank, this->comm);

        return results;

      }

      /**
       * @brief count elements with the specified keys in the distributed unordered_multimap.
       * @param first
       * @param last
       */
      template <typename Predicate>
      ::std::vector<::std::pair<Key, size_type> > count_if(::std::vector<Key>& keys, Predicate const & pred) const {
        ::std::vector<::std::pair<Key, size_type> > results;
        // even if count is 0, still need to participate in mpi calls.  if (keys.size() == 0) return results;

        // keep unique keys
        retain_unique_keys(keys);

        if (comm_size > 1) {
          // distribute (communication part)
          std::vector<int> recv_counts;

          recv_counts = mxx2::msgs_all2all(keys, this->key_to_rank, comm);

          // local find. memory utilization a potential problem.
          // do for each src proc one at a time.

          results.reserve(keys.size());                   // TODO:  should estimate coverage.
          auto start = keys.begin();
          auto end = start;
          for (int i = 0; i < comm_size; ++i) {
            ::std::advance(end, recv_counts[i]);

            // work on query from process i.
            local_count_if( start, end, results, pred);

            if (comm_rank == 0) DEBUGF("R %d added %d results for %d queries for process %d\n", comm_rank, send_counts[i], recv_counts[i], i);

            start = end;
          }

          // send back using the constructed recv count
          mxx::all2all(results, recv_counts, comm);

        } else {
          local_count_if(keys.begin(), keys.end(), results, pred);
        }

        return results;

      }


      template <typename Predicate>
      ::std::vector<::std::pair<Key, size_type> > count_if(Predicate const & pred) const {
        ::std::vector<::std::pair<Key, size_type> > results;
        auto keys = this->keys();

        size_t count = local_count_if(keys.begin(), keys.end(), results, pred);

        if (comm_size > 1) MPI_Barrier(comm);
        return results;
      }

      /**
       * @brief erase elements with the specified keys in the distributed unordered_multimap.
       * @param first
       * @param last
       */
      void erase(::std::vector<Key>& keys) {
        // even if count is 0, still need to participate in mpi calls.  if (keys.size() == 0) return;

          if (comm_size > 1) {

            // remove duplicates
            retain_unique_keys(keys);

            // redistribute keys
            mxx2::msgs_all2all(keys, this->key_to_rank, comm);
            
//
//            // remove duplicates again
//            retain_unique_keys(keys);
          }

          // then call local remove.
          local_erase(keys.begin(), keys.end());
      }

      /**
       * @brief erase elements with the specified keys in the distributed unordered_multimap.
       * @param first
       * @param last
       */
      template <typename Predicate>
      void erase_if(::std::vector<Key>& keys, Predicate const & pred) {
        // even if count is 0, still need to participate in mpi calls.  if (keys.size() == 0) return;


          if (comm_size > 1) {
            // remove duplicates
            retain_unique_keys(keys);

            // redistribute keys
            mxx2::msgs_all2all(keys, this->key_to_rank, comm);

//
//            // remove duplicates again
//            retain_unique_keys(keys);
          }

          // then call local remove.
          local_erase_if(keys.begin(), keys.end(), pred);
      }

      template <typename Predicate>
      void erase_if(Predicate const & pred) {
          auto keys = this->keys();

          size_t count = local_erase_if(keys.begin(), keys.end(), pred);

          if (comm_size > 1) MPI_Barrier(comm);
      }

      // update done via erase/insert.

      /// clears the unordered_map
      void clear() noexcept {
        // clear + barrier.
        local_clear();
        if (comm_size > 1) MPI_Barrier(comm);
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
      using local_iterator        = typename local_container_type::local_iterator;
      using const_local_iterator  = typename local_container_type::const_local_iterator;
      using size_type             = typename local_container_type::size_type;
      using difference_type       = typename local_container_type::difference_type;

    protected:

      // defined Communicator as a friend
      friend Comm;


      /**
       * @brief find elements with the specified keys in the distributed unordered_multimap.
       * @param first
       * @param last
       */
      template <class InputIterator>
      size_t local_find(InputIterator first, InputIterator last, ::std::vector<::std::pair<Key, T> > & output) const {
          if (first == last) return 0;

          size_t before = output.size();  // before size.
          for (auto it = first; it != last; ++it) {
             auto iter = this->c.find(*it);

            if (iter != this->c.end()) {
              output.push_back(*iter);
            }  // no insert if can't find it.
          }
          return output.size() - before;  // after size.
      }

      /**
       * @brief find elements with the specified keys in the distributed unordered_multimap.
       * @param first
       * @param last
       */
      template <class InputIterator, class Predicate>
      size_t local_find_if(InputIterator first, InputIterator last, ::std::vector<::std::pair<Key, T> > & output, Predicate const & pred) const {
          if (first == last) return 0;

          size_t before = output.size();  // before size.
          for (auto it = first; it != last; ++it) {
            if (!pred(this->trans(*it))) continue;

            auto iter = this->c.find(*it);

            if (iter != this->c.end()) {
              output.push_back(*iter);
            }  // no insert if can't find it.
          }
          return output.size() - before;  // after size.
      }


      void retain_unique_tuple(::std::vector<::std::pair<Key, T> >& input) const {
        if (input.size() == 0) return;

        // sorting can be slow.
//        if (! ::std::is_sorted(input.begin(), input.end(), this->key_less_op))
//          ::std::stable_sort(input.begin(), input.end(), this->key_less_op);
//        auto end = ::std::unique(input.begin(), input.end(), this->key_equal_op);
//        input.erase(end, input.end());

        // use the same data structure to do the uniqueness reduction
        local_container_type temp;
        temp.reserve(input.size());
        temp.insert(input.begin(), input.end());  // keep first entry only
        input.assign(temp.begin(), temp.end());   // copy back to input.
      }

    public:

      unordered_map(MPI_Comm _comm, int _comm_size) : Base(_comm, _comm_size) {}

      virtual ~unordered_map() {};

      /**
       * @brief find elements with the specified keys in the distributed unordered_multimap.
       * @param first
       * @param last
       */
      ::std::vector<::std::pair<Key, T> > find(::std::vector<Key>& keys) const {
        TIMER_INIT(find);

        TIMER_START(find);
        ::std::vector<::std::pair<Key, T> > results;
        // even if count is 0, still need to participate in mpi calls.  if (keys.size() == 0) return results;
        TIMER_END(find, "begin", keys.size());

        // keep unique keys
        TIMER_START(find);
        this->retain_unique_keys(keys);
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
            send_counts[i] = local_find( start, end, results);

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
          local_find(keys.begin(), keys.end(), results);
          TIMER_END(find, "local_find", results.size());

        }

        TIMER_REPORT_MPI(find, this->comm_rank, this->comm);

        return results;
      }

      /**
       * @brief find elements with the specified keys in the distributed unordered_multimap.
       * @param first
       * @param last
       */
      template <class Predicate>
      ::std::vector<::std::pair<Key, T> > find_if(::std::vector<Key>& keys, Predicate const & pred) const {
        ::std::vector<::std::pair<Key, T> > results;
        // even if count is 0, still need to participate in mpi calls.  if (keys.size() == 0) return results;

        // keep unique keys
        this->retain_unique_keys(keys);


        if (this->comm_size > 1) {
          // distribute (communication part)
          std::vector<int> recv_counts = mxx2::msgs_all2all(keys, this->key_to_rank, this->comm);

          // local find. memory utilization a potential problem.
          // do for each src proc one at a time.

          std::vector<int> send_counts(this->comm_size, 0);
          results.reserve(keys.size() * this->key_multiplicity);  // 1 result per key.
          auto start = keys.begin();
          auto end = start;
          for (int i = 0; i < this->comm_size; ++i) {
            ::std::advance(end, recv_counts[i]);

            // work on query from process i.
            send_counts[i] = local_find_if( start, end, results, pred);

            if (this->comm_rank == 0) DEBUGF("R %d added %d results for %d queries for process %d\n", this->comm_rank, send_counts[i], recv_counts[i], i);

            start = end;
          }

          // send back using the constructed recv count
          mxx::all2all(results, send_counts, this->comm);

        } else {
          local_find_if(keys.begin(), keys.end(), results, pred);
        }

        return results;
      }

      template <class Predicate>
      ::std::vector<::std::pair<Key, T> > find_if(Predicate const & pred) const {
          ::std::vector<::std::pair<Key, T> > results;

          auto keys = this->keys();
          local_find_if(keys.begin(), keys.end(), results, pred);

          return results;
      }


      /**
       * @brief insert new elements in the distributed unordered_multimap.
       * @param first
       * @param last
       */
      void insert(std::vector<::std::pair<Key, T> >& input) {
        // even if count is 0, still need to participate in mpi calls.  if (input.size() == 0) return;
        TIMER_INIT(insert);

        TIMER_START(insert);
        TIMER_END(insert, "start", input.size());


        TIMER_START(insert);
        // first remove duplicates.  sort, then get unique, finally remove the rest.  may not be needed
        retain_unique_tuple(input);
        TIMER_END(insert, "uniq1", input.size());

        TIMER_START(insert);
        // communication part
        if (this->comm_size > 1) {
          mxx2::msgs_all2all(input, this->key_to_rank, this->comm);
        }
        TIMER_END(insert, "a2a", input.size());
//
//        TIMER_START(insert);
//        // after communication, sort again to keep unique  - may not be needed
//        retain_unique_tuple(input);
//        TIMER_END(insert, "uniq2", input.size());

        TIMER_START(insert);
        // local compute part.  called by the communicator.
        this->Base::local_insert(input.begin(), input.end());
        TIMER_END(insert, "insert", this->csize());

        TIMER_REPORT_MPI(insert, this->comm_rank, this->comm);
      }


      /**
       * @brief insert new elements in the distributed unordered_multimap.
       * @param first
       * @param last
       */
      template <class Predicate>
      void insert_if(std::vector<::std::pair<Key, T> >& input, Predicate const & pred) {
          // even if count is 0, still need to participate in mpi calls.  if (input.size() == 0) return;

        // first remove duplicates.  sort, then get unique, finally remove the rest.  may not be needed
        retain_unique_tuple(input);

        // communication part
        if (this->comm_size > 1) {
          mxx2::msgs_all2all(input, this->key_to_rank, this->comm);
        }
//
//        // after communication, sort again to keep unique  - may not be needed
//        retain_unique_tuple(input);

        // local compute part.  called by the communicator.
        this->Base::local_insert_if(input.begin(), input.end(), pred);
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
      using local_iterator        = typename local_container_type::local_iterator;
      using const_local_iterator  = typename local_container_type::const_local_iterator;
      using size_type             = typename local_container_type::size_type;
      using difference_type       = typename local_container_type::difference_type;

    protected:

      // defined Communicator as a friend
      friend Comm;


      /**
       * @brief find elements with the specified keys in the distributed unordered_multimap.
       * @param first
       * @param last
       */
      template <class InputIterator>
      size_t local_find(InputIterator first, InputIterator last, ::std::vector<::std::pair<Key, T> > &output) const {
          if (first == last) return 0;

          size_t before = output.size();  // before size.
          for (auto it = first; it != last; ++it) {
             auto range = this->c.equal_range(*it);

             // range's iterators are not random access iterators, so insert needs to call distance repeatedly, slowing down the process.
             // manually insert improves performance here.
             for (auto it2 = range.first; it2 != range.second; ++it2) {
               output.push_back(*it2);
             }
//            if (range.first != range.second) {
//              output.insert(output.end(), range.first, range.second);
//            }  // no insert if can't find it.
          }
          return output.size() - before;  // after size.
      }

      /**
       * @brief find elements with the specified keys in the distributed unordered_multimap.
       * @param first
       * @param last
       */
      template <class InputIterator, class Predicate>
      size_t local_find_if(InputIterator first, InputIterator last, ::std::vector<::std::pair<Key, T> > &output, Predicate const& pred) const {
          if (first == last) return 0;

          size_t before = output.size();  // before size.
          for (auto it = first; it != last; ++it) {
            if (!pred(this->trans(*it))) continue;

             auto range = this->c.equal_range(*it);

             // range's iterators are not random access iterators, so insert needs to call distance repeatedly, slowing down the process.
             // manually insert improves performance here.
             for (auto it2 = range.first; it2 != range.second; ++it2) {
               output.push_back(*it2);
             }
//            if (range.first != range.second) {
//              output.insert(output.end(), range.first, range.second);
//            }  // no insert if can't find it.
          }
          return output.size() - before;  // after size.
      }

    public:


      unordered_multimap(MPI_Comm _comm, int _comm_size) : Base(_comm, _comm_size) {
        this->key_multiplicity = 50;
      }

      virtual ~unordered_multimap() {}

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
//        size_t uniq_count = 0;
//        ::std::vector< ::std::pair<Key, T> > temp;
//        for (int i = 0, max = this->c.bucket_count(); i < max; ++i) {
//          if (this->c.bucket_size(i) == 0) continue;  // empty bucket. move on.
//
//          // copy and sort.
//          temp.assign(this->c.begin(i), this->c.end(i));  // copy the bucket
//          // sort the bucket
//          ::std::sort(temp.begin(), temp.end(), this->key_less_op);
// //          auto end = ::std::unique(temp.begin(), temp.end(), this->key_equal_op);
// //          uniq_count += ::std::distance(temp.begin(), end);
//
//          // count via linear scan..
//          auto x = temp.begin();
//          ++uniq_count;  // first entry.
//          // compare pairwise.
//          auto y = temp.begin();  ++y;
//          while (y != temp.end()) {
//            if (! (this->key_equal_op(*x, *y))) {
//              ++uniq_count;
//              x = y;
//            }
//            ++y;
//          }
//        }
//        this->key_multiplicity = (this->c.size() + uniq_count - 1) / uniq_count;


        // third approach is to assume each bucket contains only 1 kmer/kmolecule.
        // This is not generally true for all hash functions, so this is an over estimation of the repeat count.
        // we equate bucket size to the number of repeats for that key.
        // we can use mean, max, or mean+stdev.
        // max overestimates significantly with potentially value > 1000, so don't use max.  (0.0078125 human: 50 sec. synth  32 sec)
        // mean may be underestimating for well behaving hash function.   (0.0078125 human: 50 sec. synth  32 sec)
        // mean + 2 stdev gets 95% of all entries and may be a reasonable compromise.
        //    (1 stdev:  0.0078125 human: 49 sec. synth  32 sec;  2stdev: 0.0078125 human 49s synth: 33 sec)
        size_t nbuckets = 0;
        size_t sum_sqr = 0;
        size_t entry = 0;
        for (int i = 0, max = this->c.bucket_count(); i < max; ++i) {
          entry = this->c.bucket_size(i);
          if (entry == 0) continue;  // empty bucket. move on.
          ++nbuckets;
          sum_sqr += entry * entry;
        }
        double mean = static_cast<double>(this->c.size()) / static_cast<double>(nbuckets);
        double meansqr = static_cast<double>(sum_sqr) / static_cast<double>(nbuckets);
        double stdev = ::std::sqrt(meansqr - mean * mean);
        this->key_multiplicity = ::std::ceil(mean + 2.0 * stdev);  // covers 95% of data.
        printf("stdev = %f\n", stdev);

        // finally, hard coding.  (0.0078125 human:  50 sec.  synth:  32 s)
        // this->key_multiplicity = 50;

        return this->key_multiplicity;
      }


      /**
       * @brief find elements with the specified keys in the distributed unordered_multimap.
       * @param first
       * @param last
       */
      ::std::vector<::std::pair<Key, T> > find(::std::vector<Key>& keys) const {
        TIMER_INIT(find);

        TIMER_START(find);
        ::std::vector<::std::pair<Key, T> > results;
        // even if count is 0, still need to participate in mpi calls.  if (keys.size() == 0) return results;
        TIMER_END(find, "begin", keys.size());

        TIMER_START(find);
        // keep unique keys
        this->retain_unique_keys(keys);
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
          TIMER_END(find, "reserve", (keys.size() * this->key_multiplicity));

          TIMER_START(find);
          auto start = keys.begin();
          auto end = start;
          for (int i = 0; i < this->comm_size; ++i) {
            ::std::advance(end, recv_counts[i]);

            // work on query from process i.
            send_counts[i] = local_find( start, end, results);

            if (this->comm_rank == 0) DEBUGF("R %d added %d results for %d queries for process %d\n", this->comm_rank, send_counts[i], recv_counts[i], i);

            start = end;
          }
          TIMER_END(find, "local_find", results.size());

          TIMER_START(find);
          // send back using the constructed recv count
          mxx::all2all(results, send_counts, this->comm);
          TIMER_END(find, "a2a2", results.size());

        } else {
          TIMER_START(find);
          local_find(keys.begin(), keys.end(), results);
          TIMER_END(find, "local_find", results.size());
        }

        TIMER_REPORT_MPI(find, this->comm_rank, this->comm);

        return results;

      }

      /**
       * @brief find elements with the specified keys in the distributed unordered_multimap.
       * @param first
       * @param last
       */
      template <typename Predicate>
      ::std::vector<::std::pair<Key, T> > find_if(::std::vector<Key>& keys, Predicate const& pred) const {
        ::std::vector<::std::pair<Key, T> > results;
        // even if count is 0, still need to participate in mpi calls.  if (keys.size() == 0) return results;

        // keep unique keys
        this->retain_unique_keys(keys);


        if (this->comm_size > 1) {
          // distribute (communication part)
          std::vector<int> recv_counts = mxx2::msgs_all2all(keys, this->key_to_rank, this->comm);

          // local find. memory utilization a potential problem.
          // do for each src proc one at a time.

          std::vector<int> send_counts(this->comm_size, 0);
          results.reserve(keys.size() * this->key_multiplicity);                   // TODO:  should estimate coverage.
          auto start = keys.begin();
          auto end = start;
          for (int i = 0; i < this->comm_size; ++i) {
            ::std::advance(end, recv_counts[i]);

            // work on query from process i.
            send_counts[i] = local_find_if( start, end, results, pred);

            if (this->comm_rank == 0) DEBUGF("R %d added %d results for %d queries for process %d\n", this->comm_rank, send_counts[i], recv_counts[i], i);

            start = end;
          }

          // send back using the constructed recv count
          mxx::all2all(results, send_counts, this->comm);

        } else {
          local_find_if(keys.begin(), keys.end(), results, pred);
        }

        return results;

      }

      template <typename Predicate>
      ::std::vector<::std::pair<Key, T> > find_if(Predicate const& pred) const {
        ::std::vector<::std::pair<Key, T> > results;

        auto keys = this->keys();
        local_find_if(keys.begin(), keys.end(), results, pred);

        return results;
      }

      /**
       * @brief insert new elements in the distributed unordered_multimap.
       * @param first
       * @param last
       */
      void insert(std::vector<::std::pair<Key, T> >& input) {
        // even if count is 0, still need to participate in mpi calls.  if (input.size() == 0) return;
        TIMER_INIT(insert);



        TIMER_START(insert);

        // communication part
        if (this->comm_size > 1) {
          mxx2::msgs_all2all(input, this->key_to_rank, this->comm);
        }
        TIMER_END(insert, "a2a", input.size());

        TIMER_START(insert);
        // local compute part.  called by the communicator.
        this->Base::local_insert(input.begin(), input.end());
        TIMER_END(insert, "insert", this->c.size());

        TIMER_REPORT_MPI(insert, this->comm_rank, this->comm);

      }

      /**
       * @brief insert new elements in the distributed unordered_multimap.
       * @param first
       * @param last
       */
      template <typename Predicate>
      void insert_if(std::vector<::std::pair<Key, T> >& input, Predicate const & pred) {
        // even if count is 0, still need to participate in mpi calls.  if (input.size() == 0) return;

        // communication part
        if (this->comm_size > 1) {
          mxx2::msgs_all2all(input, this->key_to_rank, this->comm);
        }

        // local compute part.  called by the communicator.
        this->Base::local_insert_if(input.begin(), input.end(), pred);
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
      using local_iterator        = typename local_container_type::local_iterator;
      using const_local_iterator  = typename local_container_type::const_local_iterator;
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
      void local_insert(InputIterator first, InputIterator last) {
          for (auto it = first; it != last; ++it) {
            this->c[it->first] = r(this->c[it->first], it->second);
          }
      }

      /**
       * @brief insert new elements in the distributed unordered_multimap.
       * @param first
       * @param last
       */
      template <class InputIterator, class Predicate>
      void local_insert_if(InputIterator first, InputIterator last, Predicate const & pred) {
          for (auto it = first; it != last; ++it) {
            if (pred(this->trans(*it))) this->c[it->first] = r(this->c[it->first], it->second);
          }
      }

      void local_reduction(::std::vector<::std::pair<Key, T> >& input) {
  
        if (input.size() == 0) return;

        // sort is slow
//        if (! ::std::is_sorted(input.begin(), input.end(), this->key_less_op))
//          ::std::sort(input.begin(), input.end(), this->key_less_op);
//
//        // then do reduction
//        auto first = input.begin();
//        auto curr = first;  ++curr;
//        while (curr != input.end()) {
//          if (this->key_equal_op(*first, *curr))  // if same, do reduction
//            first->second = r(first->second, curr->second);
//          else  // else reset first.
//            first = curr;
//          ++curr;  // increment second.
//        }
//
//        // keep the first element in each replicated range
//        auto end = ::std::unique(input.begin(), input.end(), this->key_equal_op);
//        input.erase(end, input.end());

        // use the target container type.

        TIMER_INIT(reduce_tuple);

        TIMER_START(reduce_tuple);
        local_container_type temp;
        temp.reserve(input.size());
        TIMER_END(reduce_tuple, "reserve", input.size());

        TIMER_START(reduce_tuple);
        for (auto it = input.begin(), end = input.end(); it != end; ++it) {
          temp[it->first] = r(temp[it->first], it->second);
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
      void insert(std::vector<::std::pair<Key, T> >& input) {
        // even if count is 0, still need to participate in mpi calls.  if (input.size() == 0) return;

        // first remove duplicates.  sort, then get unique, finally remove the rest.  may not be needed
        local_reduction(input);

        // communication part
        if (this->comm_size > 1) {
          mxx2::msgs_all2all(input, this->key_to_rank, this->comm);
        }
//
//        // after communication, sort again to keep unique  - may not be needed
//        local_reduction(input);

        // local compute part.  called by the communicator.
        local_insert(input.begin(), input.end());
      }

      /**
       * @brief insert new elements in the distributed unordered_multimap.
       * @param first
       * @param last
       */
      template <typename Predicate>
      void insert_if(std::vector<::std::pair<Key, T> >& input, Predicate const & pred) {
        // even if count is 0, still need to participate in mpi calls.  if (input.size() == 0) return;

        // first remove duplicates.  sort, then get unique, finally remove the rest.  may not be needed
        local_reduction(input);

        // communication part
        if (this->comm_size > 1) {
          mxx2::msgs_all2all(input, this->key_to_rank, this->comm);
        }
//
//        // after communication, sort again to keep unique  - may not be needed
//        local_reduction(input);

        // local compute part.  called by the communicator.
        local_insert_if(input.begin(), input.end(), pred);
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
      using local_iterator        = typename local_container_type::local_iterator;
      using const_local_iterator  = typename local_container_type::const_local_iterator;
      using size_type             = typename local_container_type::size_type;
      using difference_type       = typename local_container_type::difference_type;

    protected:

      // defined Communicator as a friend
      friend Comm;

      ::std::vector<::std::pair<Key, T> > local_reduction(::std::vector<Key>& input) {

        TIMER_INIT(local_reduc);

        TIMER_START(local_reduc);
        ::std::vector<::std::pair<Key, T> > output;
        if (input.size() == 0) return output;
        TIMER_END(local_reduc, "start", input.size());


        // sort is slow.
//        TIMER_START(local_reduc);
//        if (! ::std::is_sorted(input.begin(), input.end(), this->key_less_op))
//          ::std::sort(input.begin(), input.end(), this->key_less_op);
//
//        TIMER_END(local_reduc, "sort", input.size());
//
//
//        // reserve
//        TIMER_START(local_reduc);
//        // then do reduction
//        output.reserve(input.size());
//        TIMER_END(local_reduc, "reserve", input.size());
//
//
//        // reduce
//        TIMER_START(local_reduc);
//        auto first = input.begin();
//        output.emplace_back(*first, 1);
//
//        auto curr = first;  ++curr;
//        while (curr != input.end()) {
//          if (this->key_equal_op(*first, *curr))  // if same, do reduction
//            ++(output.back().second);
//          else { // else reset first.
//            first = curr;
//            output.emplace_back(*first, 1);
//          }
//          ++curr;  // increment second.
//        }
//        TIMER_END(local_reduc, "reduce", input.size());

        // use target container type instead
        TIMER_START(local_reduc);
        local_container_type temp;
        temp.reserve(input.size());
        for (auto it = input.begin(), end = input.end(); it != end; ++it) {
          temp[*it]++;
        }
        TIMER_END(local_reduc, "reduce", temp.size());

        // reserve
        TIMER_START(local_reduc);
        // then do reduction
        output.reserve(temp.size());
        TIMER_END(local_reduc, "reserve", temp.size());


        TIMER_START(local_reduc);
        output.assign(temp.begin(), temp.end());
        TIMER_END(local_reduc, "copy", output.size());

        TIMER_REPORT_MPI(local_reduc, this->comm_rank, this->comm);

        return output;
      }


    public:
      counting_unordered_map(MPI_Comm _comm, int _comm_size) : Base(_comm, _comm_size) {}

      virtual ~counting_unordered_map() {};


      /**
       * @brief insert new elements in the distributed unordered_multimap.
       * @param first
       * @param last
       */
      void insert(std::vector< Key >& input) {
        // even if count is 0, still need to participate in mpi calls.  if (input.size() == 0) return;
        TIMER_INIT(count_insert);

        TIMER_START(count_insert);
        TIMER_END(count_insert, "start", input.size());


        // distribute
        TIMER_START(count_insert);

        // first remove duplicates.  sort, then get unique, finally remove the rest.  may not be needed
        auto temp = local_reduction(input);
        TIMER_END(count_insert, "reduc1", temp.size());


        // distribute
        TIMER_START(count_insert);

        // communication part
        if (this->comm_size > 1) {
          mxx2::msgs_all2all(temp, this->key_to_rank, this->comm);
        }
        TIMER_END(count_insert, "a2a", temp.size());


//        // distribute
//        TIMER_START(count_insert);
//
//        // after communication, sort again to keep unique  - NOT be needed
//        this->Base::local_reduction(temp);
//        TIMER_END(count_insert, "reduc2", temp.size());


        // distribute
        TIMER_START(count_insert);

        // local compute part.  called by the communicator.
        this->Base::local_insert(temp.begin(), temp.end());
        TIMER_END(count_insert, "insert", this->c.size());


        // distribute
        TIMER_REPORT_MPI(count_insert, this->comm_rank, this->comm);

      }

      /**
       * @brief insert new elements in the distributed unordered_multimap.
       * @param first
       * @param last
       */
      template <typename Predicate>
      void insert_if(std::vector< Key >& input, Predicate &pred) {
        // even if count is 0, still need to participate in mpi calls.  if (input.size() == 0) return;

        // first remove duplicates.  sort, then get unique, finally remove the rest.  may not be needed
        auto temp = local_reduction(input);

        // communication part
        if (this->comm_size > 1) {
          mxx2::msgs_all2all(temp, this->key_to_rank, this->comm);
        }
//
//        // after communication, sort again to keep unique  - NOT be needed
//        this->Base::local_reduction(temp);

        // local compute part.  called by the communicator.
        this->Base::local_insert_if(temp.begin(), temp.end(), pred);
      }


  };



} /* namespace dsc */


#endif // BLISS_DISTRIBUTED_MAP_HPP
