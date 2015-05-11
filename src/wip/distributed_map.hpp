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
#include <utility> 			  // for std::pair

//#include <sparsehash/dense_hash_map>  // local storage hash table, google dense hash map.
//#include <vector>
#include <functional> 		// for std::function and std::hash
//#include <limits>
//#include <stdexcept>
#include <algorithm> 		// for sort, etc
//#include <atomic>
//#include <tuple>
//#include "iterators/concatenating_iterator.hpp"
#include <iterator>  // advance, distance

//// include MPI
//#include <mpi.h>

#include <cstdint>  // for uint8, etc.

#include "mxx/collective.hpp"
#include "mxx/reduction.hpp"
#include "io/mpi_utils.hpp"
#include "common/kmer_hash.hpp"

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
   * @tparam Key
   * @tparam T
   * @tparam Comm   default to mpi_collective_communicator       communicator for global communication. may hash or sort.
   * @tparam Hash   default to ::std::hash<Key>       hash function for the local storage.
   * @tparam Pred   default to ::std::equal_to<Key>   equal function for the local storage.
   * @tparam Alloc  default to ::std::allocator< ::std::pair<const Key, T> >    allocator for local storage.
   */
  template<typename Key, typename T,
      template <typename, typename, typename, typename, typename> class Container,
      class Comm,
      class Hash = ::std::hash<Key>,
      class Pred = ::std::equal_to<Key>,
      class Alloc = ::std::allocator< ::std::pair<const Key, T> >
  >
  class unordered_map_base {
    public:
      using local_container_type = Container<Key, T, Hash, Pred, Alloc>;

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
      bliss::hash::kmer::hash<Key, bliss::hash::kmer::detail::farm::hash, bliss::hash::kmer::LexicographicLessCombiner, true> hashf;
      MPI_Comm comm;
      int comm_size;
      int comm_rank;

      // defined Communicator as a friend
      friend Comm;



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
       * @param pred
       * @return
       */
      template<class InputIterator, class Predicate>
      size_t local_count_if(InputIterator first, InputIterator last, ::std::vector<::std::pair<Key, size_type> > &output, Predicate const & pred) const {
          if (first == last) return 0;

          size_t before = output.size();  // before size.
          for (auto it = first; it != last; ++it) {
             if (pred(*it)) output.emplace_back(*it, c.count(*it));
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
            if (pred(*it)) c.insert(*it);
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
          if (pred(*it)) c.erase(*it);
        }
      }

      /// clears the unordered_map
      void local_clear() noexcept {
          c.clear();
      }

      void retain_unique_keys(::std::vector< Key >& input) const {
        if (input.size() == 0) return;

        if (! ::std::is_sorted(input.begin(), input.end()))
          ::std::stable_sort(input.begin(), input.end());
        auto end = ::std::unique(input.begin(), input.end());
        input.erase(end, input.end());
      }



      unordered_map_base(MPI_Comm _comm, int _comm_size) : key_multiplicity(1), hashf(ceilLog2(comm_size)), comm(_comm), comm_size(_comm_size) {
        MPI_Comm_rank(comm, &comm_rank);
      }





    public:

      virtual ~unordered_map_base() {};

      /// returns the local storage.  please use sparingly.
      local_container_type& get_local_container() { return c; }

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
        ::std::vector<::std::pair<Key, size_type> > results;
        if (keys.size() == 0) return results;

        // keep unique keys
        retain_unique_keys(keys);


        if (comm_size > 1) {
          // distribute (communication part)
          std::vector<int> recv_counts;

          recv_counts = mxx2::msgs_all2all(keys, [&] ( Key const &x) {
             return (this->hashf(x) % comm_size);
           }, comm);

          // local find. memory utilization a potential problem.
          // do for each src proc one at a time.

          results.reserve(keys.size());                   // TODO:  should estimate coverage.
          auto start = keys.begin();
          auto end = start;
          for (int i = 0; i < comm_size; ++i) {
            ::std::advance(end, recv_counts[i]);

            // work on query from process i.
            local_count( start, end, results);

            if (comm_rank == 0) DEBUGF("R %d added %d results for %d queries for process %d\n", comm_rank, send_counts[i], recv_counts[i], i);

            start = end;
          }

          // send back using the constructed recv count
          mxx::all2all(results, recv_counts, comm);

        } else {
          local_count(keys.begin(), keys.end(), results);
        }

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
        if (keys.size() == 0) return results;

        // keep unique keys
        retain_unique_keys(keys);

        if (comm_size > 1) {
          // distribute (communication part)
          std::vector<int> recv_counts;

          recv_counts = mxx2::msgs_all2all(keys, [&] ( Key const &x) {
             return (this->hashf(x) % comm_size);
           }, comm);

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


      /**
       * @brief erase elements with the specified keys in the distributed unordered_multimap.
       * @param first
       * @param last
       */
      void erase(::std::vector<Key>& keys) {
          if (keys.size() == 0) return;

          // remove duplicates
          retain_unique_keys(keys);

          if (comm_size > 1) {
            // redistribute keys
            mxx2::msgs_all2all(keys, [&] ( Key const &x) {
               return (this->hashf(x) % comm_size);
             }, comm);
            
          
            // remove duplicates again
            retain_unique_keys(keys);
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
          if (keys.size() == 0) return;

          // remove duplicates
          retain_unique_keys(keys);

          if (comm_size > 1) {
            // redistribute keys
            mxx2::msgs_all2all(keys, [&] ( Key const &x) {
               return (this->hashf(x) % comm_size);
             }, comm);


            // remove duplicates again
            retain_unique_keys(keys);
          }

          // then call local remove.
          local_erase_if(keys.begin(), keys.end(), pred);
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
   * @tparam Hash   default to ::std::hash<Key>       hash function for the local storage.
   * @tparam Pred   default to ::std::equal_to<Key>   equal function for the local storage.
   * @tparam Alloc  default to ::std::allocator< ::std::pair<const Key, T> >    allocator for local storage.
   */
  template<typename Key, typename T,
      class Comm,
      class Hash = ::std::hash<Key>,
      class Pred = ::std::equal_to<Key>,
      class Alloc = ::std::allocator< ::std::pair<const Key, T> >
  >
  class unordered_map : public unordered_map_base<Key, T, ::std::unordered_map, Comm, Hash, Pred, Alloc> {
      using Base = unordered_map_base<Key, T, ::std::unordered_map, Comm, Hash, Pred, Alloc>;


    public:
      using local_container_type = ::std::unordered_map<Key, T, Hash, Pred, Alloc>;

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

      struct TupleLess {
        bool operator()(::std::pair<Key, T> const &x, ::std::pair<Key, T> const &y) {
          return x.first < y.first;
        }
      } tuple_less;

      struct TupleEqual {
        bool operator()(::std::pair<Key, T> const &x, ::std::pair<Key, T> const &y) {
          return x.first == y.first;
        }
      } tuple_equal;

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
              output.insert(output.end(), *iter);
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
            if (!pred(*it)) continue;

            auto iter = this->c.find(*it);

            if (iter != this->c.end()) {
              output.insert(output.end(), *iter);
            }  // no insert if can't find it.
          }
          return output.size() - before;  // after size.
      }

      void retain_unique(::std::vector<::std::pair<Key, T> >& input) const {
        if (input.size() == 0) return;

        if (! ::std::is_sorted(input.begin(), input.end(), tuple_less))
          ::std::stable_sort(input.begin(), input.end(), tuple_less);
        auto end = ::std::unique(input.begin(), input.end(), tuple_equal);
        input.erase(end, input.end());
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
        ::std::vector<::std::pair<Key, T> > results;
        if (keys.size() == 0) return results;

        // keep unique keys
        this->retain_unique_keys(keys);


        if (this->comm_size > 1) {
          // distribute (communication part)
          std::vector<int> recv_counts = mxx2::msgs_all2all(keys, [&] ( Key const &x) {
             return (this->hashf(x) % this->comm_size);
           }, this->comm);

          // local find. memory utilization a potential problem.
          // do for each src proc one at a time.

          std::vector<int> send_counts(this->comm_size, 0);
          results.reserve(keys.size() * this->key_multiplicity);  // 1 result per key.
          auto start = keys.begin();
          auto end = start;
          for (int i = 0; i < this->comm_size; ++i) {
            ::std::advance(end, recv_counts[i]);

            // work on query from process i.
            send_counts[i] = local_find( start, end, results);

            if (this->comm_rank == 0) DEBUGF("R %d added %d results for %d queries for process %d\n", this->comm_rank, send_counts[i], recv_counts[i], i);

            start = end;
          }

          // send back using the constructed recv count
          mxx::all2all(results, send_counts, this->comm);

        } else {
          local_find(keys.begin(), keys.end(), results);
        }

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
        if (keys.size() == 0) return results;

        // keep unique keys
        this->retain_unique_keys(keys);


        if (this->comm_size > 1) {
          // distribute (communication part)
          std::vector<int> recv_counts = mxx2::msgs_all2all(keys, [&] ( Key const &x) {
             return (this->hashf(x) % this->comm_size);
           }, this->comm);

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


      /**
       * @brief insert new elements in the distributed unordered_multimap.
       * @param first
       * @param last
       */
      void insert(std::vector<::std::pair<Key, T> >& input) {
        if (input.size() == 0) return;

        // first remove duplicates.  sort, then get unique, finally remove the rest.  may not be needed
        retain_unique(input);

        // communication part
        if (this->comm_size > 1) {
          mxx2::msgs_all2all(input, [&] ( ::std::pair<Key, T> const &x) {
             return (this->hashf(x.first) % this->comm_size);
           }, this->comm);
        }

        // after communication, sort again to keep unique  - may not be needed
        retain_unique(input);

        // local compute part.  called by the communicator.
        this->Base::local_insert(input.begin(), input.end());
      }


      /**
       * @brief insert new elements in the distributed unordered_multimap.
       * @param first
       * @param last
       */
      template <class Predicate>
      void insert_if(std::vector<::std::pair<Key, T> >& input, Predicate const & pred) {
        if (input.size() == 0) return;

        // first remove duplicates.  sort, then get unique, finally remove the rest.  may not be needed
        retain_unique(input);

        // communication part
        if (this->comm_size > 1) {
          mxx2::msgs_all2all(input, [&] ( ::std::pair<Key, T> const &x) {
             return (this->hashf(x.first) % this->comm_size);
           }, this->comm);
        }

        // after communication, sort again to keep unique  - may not be needed
        retain_unique(input);

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
   * @tparam Hash   default to ::std::hash<Key>       hash function for the local storage.
   * @tparam Pred   default to ::std::equal_to<Key>   equal function for the local storage.
   * @tparam Alloc  default to ::std::allocator< ::std::pair<const Key, T> >    allocator for local storage.
   */
  template<typename Key, typename T,
      class Comm,
      class Hash = ::std::hash<Key>,
      class Pred = ::std::equal_to<Key>,
      class Alloc = ::std::allocator< ::std::pair<const Key, T> >
  >
  class unordered_multimap : public unordered_map_base<Key, T, ::std::unordered_multimap, Comm, Hash, Pred, Alloc> {
      using Base = unordered_map_base<Key, T, ::std::unordered_multimap, Comm, Hash, Pred, Alloc>;

    public:
      using local_container_type = ::std::unordered_multimap<Key, T, Hash, Pred, Alloc>;

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

            if (range.first != range.second) {
              output.insert(output.end(), range.first, range.second);
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
      size_t local_find_if(InputIterator first, InputIterator last, ::std::vector<::std::pair<Key, T> > &output, Predicate const& pred) const {
          if (first == last) return 0;

          size_t before = output.size();  // before size.
          for (auto it = first; it != last; ++it) {
            if (!pred(*it)) continue;

             auto range = this->c.equal_range(*it);

            if (range.first != range.second) {
              output.insert(output.end(), range.first, range.second);
            }  // no insert if can't find it.
          }
          return output.size() - before;  // after size.
      }


    public:


      unordered_multimap(MPI_Comm _comm, int _comm_size) : Base(_comm, _comm_size) {
        this->key_multiplicity = 50;
      }

      virtual ~unordered_multimap() {}



      /**
       * @brief find elements with the specified keys in the distributed unordered_multimap.
       * @param first
       * @param last
       */
      ::std::vector<::std::pair<Key, T> > find(::std::vector<Key>& keys) const {
        ::std::vector<::std::pair<Key, T> > results;
        if (keys.size() == 0) return results;

        // keep unique keys
        this->retain_unique_keys(keys);


        if (this->comm_size > 1) {
          // distribute (communication part)
          std::vector<int> recv_counts = mxx2::msgs_all2all(keys, [&] ( Key const &x) {
             return (this->hashf(x) % this->comm_size);
           }, this->comm);

          // local find. memory utilization a potential problem.
          // do for each src proc one at a time.

          std::vector<int> send_counts(this->comm_size, 0);
          results.reserve(keys.size() * this->key_multiplicity);                   // TODO:  should estimate coverage.
          auto start = keys.begin();
          auto end = start;
          for (int i = 0; i < this->comm_size; ++i) {
            ::std::advance(end, recv_counts[i]);

            // work on query from process i.
            send_counts[i] = local_find( start, end, results);

            if (this->comm_rank == 0) DEBUGF("R %d added %d results for %d queries for process %d\n", this->comm_rank, send_counts[i], recv_counts[i], i);

            start = end;
          }

          // send back using the constructed recv count
          mxx::all2all(results, send_counts, this->comm);

        } else {
          local_find(keys.begin(), keys.end(), results);
        }

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
        if (keys.size() == 0) return results;

        // keep unique keys
        this->retain_unique_keys(keys);


        if (this->comm_size > 1) {
          // distribute (communication part)
          std::vector<int> recv_counts = mxx2::msgs_all2all(keys, [&] ( Key const &x) {
             return (this->hashf(x) % this->comm_size);
           }, this->comm);

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



      /**
       * @brief insert new elements in the distributed unordered_multimap.
       * @param first
       * @param last
       */
      void insert(std::vector<::std::pair<Key, T> >& input) {
        if (input.size() == 0) return;

        // communication part
        if (this->comm_size > 1) {
          mxx2::msgs_all2all(input, [&] ( ::std::pair<Key, T> const &x) {
             return (this->hashf(x.first) % this->comm_size);
           }, this->comm);
        }

        // local compute part.  called by the communicator.
        this->Base::local_insert(input.begin(), input.end());
      }

      /**
       * @brief insert new elements in the distributed unordered_multimap.
       * @param first
       * @param last
       */
      template <typename Predicate>
      void insert_if(std::vector<::std::pair<Key, T> >& input, Predicate const & pred) {
        if (input.size() == 0) return;

        // communication part
        if (this->comm_size > 1) {
          mxx2::msgs_all2all(input, [&] ( ::std::pair<Key, T> const &x) {
             return (this->hashf(x.first) % this->comm_size);
           }, this->comm);
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
   * @tparam Reduc  default to ::std::plus<T>         reduction operator for the content.
   * @tparam Hash   default to ::std::hash<Key>       hash function for the local storage.
   * @tparam Pred   default to ::std::equal_to<Key>   equal function for the local storage.
   * @tparam Alloc  default to ::std::allocator< ::std::pair<const Key, T> >    allocator for local storage.
   */
  template<typename Key, typename T,
      class Comm,
      typename Reduc = ::std::plus<T>,
      class Hash = ::std::hash<Key>,
      class Pred = ::std::equal_to<Key>,
      class Alloc = ::std::allocator< ::std::pair<const Key, T> >
  >
  class reduction_unordered_map :
    public unordered_map<Key, T, Comm, Hash, Pred, Alloc> {
      static_assert(::std::is_arithmetic<T>::value, "mapped type has to be arithmetic");

      using Base = unordered_map<Key, T, Comm, Hash, Pred, Alloc>;

    public:
      using local_container_type = ::std::unordered_map<Key, T, Hash, Pred, Alloc>;

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
            if (pred(*it)) this->c[it->first] = r(this->c[it->first], it->second);
          }
      }

      void local_reduction(::std::vector<::std::pair<Key, T> >& input) {
        if (input.size() == 0) return;

        // sort
        if (! ::std::is_sorted(input.begin(), input.end(), Base::tuple_less))
          ::std::sort(input.begin(), input.end(), Base::tuple_less);

        // then do reduction
        auto first = input.begin();
        auto curr = first;  ++curr;
        while (curr != input.end()) {
          if (Base::tuple_equal(*first, *curr))  // if same, do reduction
            first->second = r(first->second, curr->second);
          else  // else reset first.
            first = curr;
          ++curr;  // increment second.
        }

        // keep the first element in each replicated range
        auto end = ::std::unique(input.begin(), input.end(), Base::tuple_equal);
        input.erase(end, input.end());
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
        if (input.size() == 0) return;

        // first remove duplicates.  sort, then get unique, finally remove the rest.  may not be needed
        local_reduction(input);

        // communication part
        if (this->comm_size > 1) {
          mxx2::msgs_all2all(input, [&] ( ::std::pair<Key, T> const &x) {
             return (this->hashf(x.first) % this->comm_size);
           }, this->comm);
        }

        // after communication, sort again to keep unique  - may not be needed
        local_reduction(input);

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
        if (input.size() == 0) return;

        // first remove duplicates.  sort, then get unique, finally remove the rest.  may not be needed
        local_reduction(input);

        // communication part
        if (this->comm_size > 1) {
          mxx2::msgs_all2all(input, [&] ( ::std::pair<Key, T> const &x) {
             return (this->hashf(x.first) % this->comm_size);
           }, this->comm);
        }

        // after communication, sort again to keep unique  - may not be needed
        local_reduction(input);

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
   * @tparam Reduc  default to ::std::plus<T>         reduction operator for the content.
   * @tparam Hash   default to ::std::hash<Key>       hash function for the local storage.
   * @tparam Pred   default to ::std::equal_to<Key>   equal function for the local storage.
   * @tparam Alloc  default to ::std::allocator< ::std::pair<const Key, T> >    allocator for local storage.
   */
  template<typename Key, typename T,
      class Comm,
      class Hash = ::std::hash<Key>,
      class Pred = ::std::equal_to<Key>,
      class Alloc = ::std::allocator< ::std::pair<const Key, T> >
  >
  class counting_unordered_map :
      public reduction_unordered_map<Key, T, Comm, ::std::plus<T>, Hash, Pred, Alloc> {
      static_assert(::std::is_integral<T>::value, "count type has to be integral");

      using Base = reduction_unordered_map<Key, T, Comm, ::std::plus<T>, Hash, Pred, Alloc>;


    public:
      using local_container_type = ::std::unordered_map<Key, T, Hash, Pred, Alloc>;

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
        ::std::vector<::std::pair<Key, T> > output;
        if (input.size() == 0) return output;

        // sort
        if (! ::std::is_sorted(input.begin(), input.end()))
          ::std::sort(input.begin(), input.end());

        // then do reduction
        output.reserve(input.size());

        auto first = input.begin();
        output.emplace_back(*first, 1);

        auto curr = first;  ++curr;
        while (curr != input.end()) {
          if (*first == *curr)  // if same, do reduction
            ++(output.back().second);
          else { // else reset first.
            first = curr;
            output.emplace_back(*first, 1);
          }
          ++curr;  // increment second.
        }

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
        if (input.size() == 0) return;

        // first remove duplicates.  sort, then get unique, finally remove the rest.  may not be needed
        auto temp = local_reduction(input);

        // communication part
        if (this->comm_size > 1) {
          mxx2::msgs_all2all(temp, [&] ( ::std::pair<Key, T> const &x) {
             return (this->hashf(x.first) % this->comm_size);
           }, this->comm);
        }

        // after communication, sort again to keep unique  - may not be needed
        this->Base::local_reduction(temp);

        // local compute part.  called by the communicator.
        this->Base::local_insert(temp.begin(), temp.end());
      }

      /**
       * @brief insert new elements in the distributed unordered_multimap.
       * @param first
       * @param last
       */
      template <typename Predicate>
      void insert_if(std::vector< Key >& input, Predicate &pred) {
        if (input.size() == 0) return;

        // first remove duplicates.  sort, then get unique, finally remove the rest.  may not be needed
        auto temp = local_reduction(input);

        // communication part
        if (this->comm_size > 1) {
          mxx2::msgs_all2all(temp, [&] ( ::std::pair<Key, T> const &x) {
             return (this->hashf(x.first) % this->comm_size);
           }, this->comm);
        }

        // after communication, sort again to keep unique  - may not be needed
        this->Base::local_reduction(temp);

        // local compute part.  called by the communicator.
        this->Base::local_insert_if(temp.begin(), temp.end(), pred);
      }


  };



} /* namespace dsc */


#endif // BLISS_DISTRIBUTED_MAP_HPP
