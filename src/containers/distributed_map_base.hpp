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
 * @file    distributed_map_base.hpp
 * @ingroup
 * @author  tpan
 * @brief
 * @details
 *
 */
#ifndef SRC_CONTAINERS_DISTRIBUTED_MAP_BASE_HPP_
#define SRC_CONTAINERS_DISTRIBUTED_MAP_BASE_HPP_

#include <functional>
#include <algorithm>
#include <iterator>
#include <vector>
#include <unordered_set>
#include "containers/container_utils.hpp"
#include <mxx/collective.hpp>

#include "utils/benchmark_utils.hpp"

namespace dsc
{

  /**
   * @brief base class for a IteratorFilter.  uses CRTP
   * @details   goal of the class is to provide subclasses that adhere to the same interface without introducing dynamic polymorphism.
   *        4 types of subclasses are possible, the first is provided here.
   *          1. static so that every evaluation returns the same "true" or "false"
   *          2. depend on element only.  operator(iter, iter) always return 0, indicating dependence on element.
   *          3. depend on range only.  constructor evaluates, and the value is return for both operators.
   *          4. combination of both.  developer should implement the appropriate logic.
   */
  template<typename Iter, typename Derived>
  struct IteratorFilter {


      /// default constructor
      IteratorFilter() {};

      /// constructor with start and end iterator.  This is to allow for evaluation based on range alone, or range + element.
      IteratorFilter(Iter b, Iter e) {};

      /// return result for an element.  value may be based on the range.  subclass define eval function.
      bool operator()(typename ::std::iterator_traits<Iter>::value_type const & x) {
        return static_cast<Derived>(*this).eval(x);
      }

      /// return result for a range.  subclass define eval function
      int operator()() {
        return static_cast<Derived>(*this).eval();
      }
  };

  struct Identity {
      template <typename T>
      bool operator()(T const & x) const { return true; }
      template <typename Iter>
      bool operator()(Iter b, Iter e) const { return true; }
  };


  template<typename Key, typename T,
      class Comm,
      template <typename> class KeyTransform,
      class Less = ::std::less<Key>,
      class Equal = ::std::equal_to<Key>,
      class Alloc = ::std::allocator< ::std::pair<Key, T> >
  >  class map_base {

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


      struct TransformedFarmHash {
          ::bliss::kmer::hash::farm<Key, false> h;

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


//      struct Greater {
//          Less lt;
//          inline bool operator()(Key const &x, Key const &y) const {
//            return lt(y, x);
//          }
//      };

      using TransformedEqual = TransformedComp<Equal>;
      using TransformedLess = TransformedComp<Less>;
//      using TransformedGreater = TransformedComp<Greater>;
      static TransformedEqual equal;
      static TransformedLess less;
//      static TransformedGreater greater;

      mutable size_t key_multiplicity;

      // communication stuff...
      const mxx::comm& comm;
      int comm_size;

      // defined Communicator as a friend
      friend Comm;

      template <class Iter>
      static void sort_ascending(Iter b, Iter e) {
          if (!::std::is_sorted(b, e, less)) {
            ::std::sort(b, e, less);
          }
      }

      virtual void local_clear() noexcept = 0;
      //virtual void local_reserve(size_t n) = 0;
      virtual bool local_empty() const noexcept = 0;
      virtual size_t local_size() const noexcept = 0;


      map_base(const mxx::comm& _comm) :
          key_multiplicity(1), comm(_comm), comm_size(comm.size()) {
      }

      /**
       * @param[IN/OUT] keys  	keys to distribute. sortedness is NOT kept because of inplace bucketing.,
       * @param[IN/OUT] sorted_input  	indicates if input is sorted.  and whether each bucket is sorted.
       * @return received counts
       */
      template <typename V, typename ToRank>
      ::std::vector<size_t> distribute_unique(::std::vector<V>& vals, ToRank const & to_rank, bool & sorted_input) const {
    	  BL_BENCH_INIT(distribute);

          BL_BENCH_START(distribute);
          // distribute (communication part)
          std::vector<size_t> send_counts = mxx::bucketing_inplace(vals, to_rank, this->comm.size());
          sorted_input = false;
          BL_BENCH_END(distribute, "bucket", vals.size());

          BL_BENCH_START(distribute);
          // distribute (communication part)
          ::fsc::bucket_unique<::std::unordered_set<V,
                                                   TransformedFarmHash,
                                                   TransformedEqual >,
                                                   TransformedEqual >(vals, send_counts, sorted_input);
          BL_BENCH_END(distribute, "unique", vals.size());


          // distribute (communication part)
          BL_BENCH_COLLECTIVE_START(distribute, "a2a", this->comm);
          vals = mxx::all2allv(vals, send_counts, this->comm);
          BL_BENCH_END(distribute, "a2a", vals.size());

          BL_BENCH_START(distribute);
          std::vector<size_t> recv_counts= mxx::all2all(send_counts, this->comm);
          BL_BENCH_END(distribute, "a2a_counts", vals.size());

          BL_BENCH_REPORT_MPI_NAMED(distribute, "map_base:distribute", this->comm);


          return recv_counts;
      }

      /**
       * @param[IN/OUT] vals  	vals to distribute.  sortedness is NOT kept because of inplace bucketing.
       * @param[IN/OUT] sorted_input  	indicates if input is sorted.  and whether each bucket is sorted.
       *
       * @return received counts.
       */
      template <typename V, typename ToRank>
      ::std::vector<size_t> distribute(::std::vector<V>& vals, ToRank const & to_rank, bool & sorted_input) const {
    	  BL_BENCH_INIT(distribute);

          BL_BENCH_START(distribute);
          // distribute (communication part)
          std::vector<size_t> send_counts = mxx::bucketing_inplace(vals, to_rank, this->comm.size());
          sorted_input = false;
          BL_BENCH_END(distribute, "bucket", vals.size());

          // distribute (communication part)
          BL_BENCH_COLLECTIVE_START(distribute, "a2a", this->comm);
          vals = mxx::all2allv(vals, send_counts, this->comm);
          BL_BENCH_END(distribute, "a2a", vals.size());

          BL_BENCH_START(distribute);
          std::vector<size_t> recv_counts= mxx::all2all(send_counts, this->comm);
          BL_BENCH_END(distribute, "a2a_counts", vals.size());

          BL_BENCH_REPORT_MPI_NAMED(distribute, "map_base:distribute", this->comm);


          return recv_counts;
      }


      /**
       * @param[IN/OUT] keys  	keys to distribute. sortedness is KEPT.,
       * @param[IN/OUT] sorted_input  	indicates if input is sorted.  and whether each bucket is sorted.
       * @return received counts
       */
      template <typename V, typename ToRank>
      ::std::vector<size_t> stable_distribute_unique(::std::vector<V>& vals, ToRank const & to_rank, bool & sorted_input) const {
    	  BL_BENCH_INIT(distribute);

          BL_BENCH_START(distribute);
          // distribute (communication part)
          std::vector<size_t> send_counts = mxx::bucketing(vals, to_rank, this->comm.size());
          BL_BENCH_END(distribute, "bucket", vals.size());

          BL_BENCH_START(distribute);
          // distribute (communication part)
          ::fsc::bucket_unique<::std::unordered_set<V,
                                                   TransformedFarmHash,
                                                   TransformedEqual >,
                                                   TransformedEqual >(vals, send_counts, sorted_input);
          BL_BENCH_END(distribute, "unique", vals.size());


          // distribute (communication part)
          BL_BENCH_COLLECTIVE_START(distribute, "a2a", this->comm);
          vals = mxx::all2allv(vals, send_counts, this->comm);
          BL_BENCH_END(distribute, "a2a", vals.size());

          BL_BENCH_START(distribute);
          std::vector<size_t> recv_counts= mxx::all2all(send_counts, this->comm);
          BL_BENCH_END(distribute, "a2a_counts", vals.size());

          BL_BENCH_REPORT_MPI_NAMED(distribute, "map_base:distribute", this->comm);


          return recv_counts;
      }

      /**
       * @param[IN/OUT] vals  	vals to distribute.  sortedness is KEPT within each bucket
       * @param[IN/OUT] sorted_input  	indicates if input is sorted.  and whether each bucket is sorted.
       *
       * @return received counts.
       */
      template <typename V, typename ToRank>
      ::std::vector<size_t> stable_distribute(::std::vector<V>& vals, ToRank const & to_rank, bool & sorted_input) const {
    	  BL_BENCH_INIT(distribute);

          BL_BENCH_START(distribute);
          // distribute (communication part)
          std::vector<size_t> send_counts = mxx::bucketing(vals, to_rank, this->comm.size());
          BL_BENCH_END(distribute, "bucket", vals.size());

          // distribute (communication part)
          BL_BENCH_COLLECTIVE_START(distribute, "a2a", this->comm);
          vals = mxx::all2allv(vals, send_counts, this->comm);
          BL_BENCH_END(distribute, "a2a", vals.size());

          BL_BENCH_START(distribute);
          std::vector<size_t> recv_counts= mxx::all2all(send_counts, this->comm);
          BL_BENCH_END(distribute, "a2a_counts", vals.size());

          BL_BENCH_REPORT_MPI_NAMED(distribute, "map_base:distribute", this->comm);


          return recv_counts;
      }



      ///  keep the unique keys in the input. primarily for reducing comm volume.
      ///  sortedness is NOT changed.  equal operator forces comparison to Key
      template <typename V>
      void hash_unique(::std::vector< V >& input, bool & sorted_input) const {
        if (input.size() == 0) return;
        if (sorted_input) {  // already sorted, then just get the unique stuff and remove rest.
          auto end = ::std::unique(input.begin(), input.end(), this->equal);
          input.erase(end, input.end());
        } else {  // not sorted, so use a set to keep the first occurence.

          // sorting is SLOW and not scalable.  use unordered set instead.  memory use is higher.
          // unordered_set for large data is memory intensive.  depending on use, bucket per processor first.
          ::std::unordered_set<V, TransformedFarmHash, TransformedEqual> temp(input.begin(),
        		  input.end(), input.size());
          input.assign(temp.begin(), temp.end());
        }
      }

      ///  keep the unique keys in the input.   output is SORTED.  equal operator forces comparison to Key
      template <typename V>
      void sort_unique(::std::vector< V >& input, bool & sorted_input) const {
        if (input.size() == 0) return;
        if (!sorted_input) map_base::sort_ascending(input.begin(), input.end());
        auto end = ::std::unique(input.begin(), input.end(), this->equal);
        input.erase(end, input.end());

        sorted_input = true;
      }


    public:
      virtual ~map_base() {};

      virtual size_t update_multiplicity() = 0;
      virtual std::vector<std::pair<Key, T> > to_vector() const = 0;
      virtual void to_vector(std::vector<std::pair<Key, T> > & result) const  = 0;
      virtual std::vector<Key> keys() const = 0;
      virtual void keys(std::vector<Key> & result) const = 0;

      /// check if empty.
      virtual bool empty() const noexcept {
        if (comm_size == 1)
          return this->local_empty();
        else // all reduce
          return mxx::all_of(this->local_empty(), comm);
      }

      /// get size of distributed container
      virtual size_t size() const noexcept {
        size_t s = this->local_size();
        if (comm_size == 1)
          return s;
        else
          return mxx::allreduce(s, comm);
      }


      /// clears the unordered_map
      void clear() noexcept {
        // clear + barrier.
        this->local_clear();
        if (comm_size > 1)
          comm.barrier();
      }

  };



  template<typename Key, typename T,
      class Comm,
      template <typename> class KeyTransform,
      class Less,
      class Equal,
      class Alloc>
  KeyTransform<Key> map_base<Key, T, Comm, KeyTransform, Less, Equal, Alloc>::trans;

  template<typename Key, typename T,
      class Comm,
      template <typename> class KeyTransform,
      class Less,
      class Equal,
      class Alloc>
  typename map_base<Key, T, Comm, KeyTransform, Less, Equal, Alloc>::TransformedLess
      map_base<Key, T, Comm, KeyTransform, Less, Equal, Alloc>::less;

  template<typename Key, typename T,
      class Comm,
      template <typename> class KeyTransform,
      class Less,
      class Equal,
      class Alloc>
  typename map_base<Key, T, Comm, KeyTransform, Less, Equal, Alloc>::TransformedEqual
      map_base<Key, T, Comm, KeyTransform, Less, Equal, Alloc>::equal;

//  template<typename Key, typename T,
//      class Comm,
//      template <typename> class KeyTransform,
//      class Less,
//      class Equal,
//      class Alloc>
//  typename map_base<Key, T, Comm, KeyTransform, Less, Equal, Alloc>::TransformedGreater
//      map_base<Key, T, Comm, KeyTransform, Less, Equal, Alloc>::greater;

}



#endif /* SRC_CONTAINERS_DISTRIBUTED_MAP_BASE_HPP_ */
