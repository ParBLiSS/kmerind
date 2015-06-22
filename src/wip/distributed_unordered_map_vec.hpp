/**
 * @file    distributed_unordered_map_vec.hpp
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

#ifndef BLISS_DISTRIBUTED_UNORDERED_MAP_VEC_HPP
#define BLISS_DISTRIBUTED_UNORDERED_MAP_VEC_HPP


#include "wip/unordered_multimap.hpp"

#include "wip/distributed_unordered_map.hpp"

namespace dsc  // distributed std container
{



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


      /**
       * @brief find elements with the specified keys in the distributed unordered_multimap.
       * @param first
       * @param last
       */
      template <class InputIterator, class OutputIterator>
      size_t local_find(InputIterator first, InputIterator last, OutputIterator &output) const {
          if (first == last) return 0;

          auto count = 0;  // before size.
          for (auto it = first; it != last; ++it) {
            auto range = this->c.equal_range_subcontainer(*it);
            count += this->c.count(*it);
            // range's iterators are not random access iterators, so insert needs to call distance repeatedly, slowing down the process.
            // manually insert improves performance here.
//            for (auto it2 = range.first; it2 != range.second; ++it2) {
//              output.emplace_back(*it2);
//            }
            //output.insert(output.end(), range.first, range.second);
            output = ::std::copy(range.first, range.second, output);  // tons faster to emplace - almost 3x faster

            //            if (range.first != range.second) {
            //              output.insert(output.end(), range.first, range.second);
            //            }  // no insert if can't find it.
          }
          return count;  // after size.
      }

      /**
       * @brief find elements with the specified keys in the distributed unordered_multimap.
       * @param first
       * @param last
       */
      template <class InputIterator, class OutputIterator, class Predicate>
      size_t local_find_if(InputIterator first, InputIterator last, OutputIterator &output, Predicate const& pred) const {
          if (first == last) return 0;

          auto count = 0;  // before size.
          for (auto it = first; it != last; ++it) {

            auto range = this->c.equal_range_subcontainer(*it);
            // if (!pred(*it)) continue;  range predicate first.

            // range's iterators are not random access iterators, so insert needs to call distance repeatedly, slowing down the process.
            // manually insert improves performance here.
            for (auto it2 = range.first; it2 != range.second; ++it2) {
              if (pred(*it2)) {
                *output = *it2;
                ++count;
                ++output;
              }
            }
            //output.insert(output.end(), range.first, range.second);
            //::std::copy_if(range.first, range.second, output, pred);  // tons faster to emplace - almost 3x.
            //            if (range.first != range.second) {
            //              output.insert(output.end(), range.first, range.second);
            //            }  // no insert if can't find it.
          }
          return count;  // after size.
      }

    public:


      unordered_multimap_vec(MPI_Comm _comm, int _comm_size) : Base(_comm, _comm_size) {
        this->key_multiplicity = 50;
      }

      virtual ~unordered_multimap_vec() {}

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
       * @brief find elements with the specified keys in the distributed unordered_multimap.
       * @param first
       * @param last
       */
      ::std::vector<::std::pair<Key, T> > find(::std::vector<Key>& keys) const {
        TIMER_INIT(find);

        TIMER_START(find);
        ::std::vector<::std::pair<Key, T> > results;
        back_emplace_iterator<::std::vector<::std::pair<Key, T> > > emplace_iter(results);

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
          //printf("reserving %lu\n", keys.size() * this->key_multiplicity);
          TIMER_END(find, "reserve", (keys.size() * this->key_multiplicity));

          TIMER_START(find);
          auto start = keys.begin();
          auto end = start;
          for (int i = 0; i < this->comm_size; ++i) {
            ::std::advance(end, recv_counts[i]);

            // work on query from process i.
            send_counts[i] = local_find( start, end, emplace_iter);

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
          results.reserve(keys.size() * this->key_multiplicity);                   // TODO:  should estimate coverage.
          TIMER_END(find, "reserve", (keys.size() * this->key_multiplicity));


          TIMER_START(find);
          local_find(keys.begin(), keys.end(), emplace_iter);
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
        back_emplace_iterator<::std::vector<::std::pair<Key, T> > > emplace_iter(results);

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
            send_counts[i] = local_find_if( start, end, emplace_iter, pred);

            if (this->comm_rank == 0) DEBUGF("R %d added %d results for %d queries for process %d\n", this->comm_rank, send_counts[i], recv_counts[i], i);

            start = end;
          }

          // send back using the constructed recv count
          mxx::all2all(results, send_counts, this->comm);

        } else {
          results.reserve(keys.size() * this->key_multiplicity);                   // TODO:  should estimate coverage.

          local_find_if(keys.begin(), keys.end(), emplace_iter, pred);
        }

        return results;

      }

      template <typename Predicate>
      ::std::vector<::std::pair<Key, T> > find_if(Predicate const& pred) const {
        ::std::vector<::std::pair<Key, T> > results;
        back_emplace_iterator<::std::vector<::std::pair<Key, T> > > emplace_iter(results);

        auto keys = this->keys();
        results.reserve(keys.size() * this->key_multiplicity);                   // TODO:  should estimate coverage.

        local_find_if(keys.begin(), keys.end(), emplace_iter, pred);

        return results;
      }


      size_t count_unique(::std::vector<::std::pair<Key, T> > const & input) const {
        // alternative approach to get number of unique keys is to use an unordered_set.  this will take more memory but probably will be faster than sort for large buckets (high repeats).
        ::std::unordered_set<Key, typename Base::TransformedHash, typename Base::Base::TransformedEqual > unique_set(this->c.size());
        for (auto it = input.begin(), max = input.end(); it != max; ++it) {
          unique_set.insert(it->first);
        }
       // printf("r %d: %lu elements, %lu unique\n", this->comm_rank, input.size(), unique_set.size());
        return unique_set.size();
      }

      template <typename _TargetP>
      ::std::vector<::std::pair<Key, T> > bucketing(::std::vector<::std::pair<Key, T> > const & msgs, _TargetP target_p_fun, MPI_Comm comm) {

        int p;
        MPI_Comm_size(comm, &p);

        // bucket input by their target processor
        // TODO: in-place bucketing??
        std::vector<int> send_counts(p, 0);
        std::vector<int> pids(msgs.size());
        for (int i = 0; i < msgs.size(); ++i)
        {
          pids[i] = target_p_fun(msgs[i]);
          send_counts[pids[i]]++;
        }

        // get all2all params
        std::vector<int> offset = mxx::get_displacements(send_counts);

        // copy.  need to be able to track current position within each block.
        ::std::vector<::std::pair<Key, T> > send_buffer;
        if (msgs.size() > 0)
          send_buffer.resize(msgs.size());
        for (int i = 0; i < msgs.size(); ++i)
        {
          send_buffer[offset[pids[i]]++] = msgs[i];
        }
        return send_buffer;
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



} /* namespace dsc */


#endif // BLISS_DISTRIBUTED_UNORDERED_MAP_HPP
