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
//
//  /**
//   * @brief base class for a IteratorFilter.  uses CRTP  for find_if, erase_if, count_if
//   * @details   goal of the class is to provide subclasses that adhere to the same interface without introducing dynamic polymorphism.
//   *        4 types of subclasses are possible, the first is provided here.
//   *          1. static so that every evaluation returns the same "true" or "false"
//   *          2. depend on element only.  operator(iter, iter) always return 0, indicating dependence on element.
//   *          3. depend on range only.  constructor evaluates, and the value is return for both operators.
//   *          4. combination of both.  developer should implement the appropriate logic.
//   */
//  template<typename Iter, typename Derived>
//  struct IteratorFilter {
//
//      /// default constructor
//      IteratorFilter() {};
//
//      /// constructor with start and end iterator.  This is to allow for evaluation based on range alone, or range + element.
//      IteratorFilter(Iter b, Iter e) {};
//
//      /// return result for an element.  value may be based on the range.  subclass define eval function.
//      bool operator()(typename ::std::iterator_traits<Iter>::value_type const & x) {
//        return static_cast<Derived>(*this).eval(x);
//      }
//
//      /// return result for a range.  subclass define eval function
//      int operator()() {
//        return static_cast<Derived>(*this).eval();
//      }
//  };
//
  struct TruePredicate {
      template <typename T>
      bool operator()(T const & x) const { return true; }
      template <typename Iter>
      bool operator()(Iter b, Iter e) const { return true; }
  };

  /**
   * @brief parameter pack for distributed map
   * @tparam Key	Key to be transformed.
   * @tparam InputTrans		template template parameter for converting input key (to another key)
   * @tparam DistTrans		template template parameter for converting key (to another key) prior to computing distribution function, e.g. lex_less
   * @tparam DistFunc		template template parameter for computing mapping from key to rank. e.g. std::hash<>, or std::less
   * @tparam DistEqual
   * @tparam StoreTrans		template template parameter for converting key (to another key) prior to computing storage function, e.g. lex_less
   * @tparam StoreFunc		template template parameter for storing key. e.g. std::hash<> or std::less.
   * @tparam StoreEqual
   * @tparam DistTransFunc  fsc transformedHash or fsc transformedComparator
   * @tparam StoreTransFunc  fsc transformedHash or fsc transformedComparator
   */
  template <typename Key,
  	  	  	  template <typename> class InputTrans,
  	  	  	  template <typename> class DistTrans,
  	  	  	  template <typename> class DistFunc,
  	  	  	  template <typename> class DistEqual,
  	  	  	  template <typename> class StoreTrans,
  	  	  	  template <typename> class StoreFunc,
  	  	  	  template <typename> class StoreEqual,
  	  	  	  template <typename, template <typename> class, template <typename> class> class DistTransFunc,
  	  	  	  template <typename, template <typename> class, template <typename> class> class StoreTransFunc
  	  	  	  >
  struct DistributedMapParams {
	  using InputTransform = InputTrans<Key>;

	  template <typename K>
	  using DistFunction = DistFunc<K>;
	  template <typename K>
	  using DistTransform = DistTrans<K>;

	  using DistributionTransformedFunction = DistTransFunc<Key, DistFunc, DistTrans>;
	  using DistributionTransformedEqual = ::fsc::TransformedComparator<Key, DistEqual, DistTrans>;

	  template <typename K>
	  using StorageFunction = StoreFunc<K>;
	  template <typename K>
	  using StorageTransform = StoreTrans<K>;

	  using StorageTransformedFunction = StoreTransFunc<Key, StoreFunc, StoreTrans>;
	  using StorageTransformedEqual = ::fsc::TransformedComparator<Key, StoreEqual, StoreTrans>;
  };

  /**
   * KeyTransformParams should be an alias of a specialization of DistributedMapParams.  see subclass for example.
   */
  template<typename Key, typename T,
      template <typename> class MapParams,
      class Alloc = ::std::allocator< ::std::pair<Key, T> >
  >  class map_base {

    protected:

	  using InputTransform = typename MapParams<Key>::InputTransform;

	  using DistFunc = typename MapParams<Key>::template DistFunction<Key>;
	  using DistTrans = typename MapParams<Key>::template DistTransform<Key>;

	  using DistTransformedFunc  = typename MapParams<Key>::DistributionTransformedFunction;
	  using DistTransformedEqual = typename MapParams<Key>::DistributionTransformedEqual;

	  using StoreTransformedFunc = typename MapParams<Key>::StorageTransformedFunction;
	  using StoreTransformedEqual = typename MapParams<Key>::StorageTransformedEqual;

	  // primarily for use with distributed_map and sorted_map, where the TransformedFunction is
	  // a comparator, but we need a hash function.
	  template <typename K>
	  using StoreFarmHash = ::bliss::kmer::hash::farm<K, false>;
	  template <typename K>
	  using StoreTransform = typename MapParams<Key>::template StorageTransform<K>;
	  using StoreTransformedFarmHash = ::fsc::TransformedHash<Key, StoreFarmHash, StoreTransform>;

	  // primarily for use fpr sorting by hashmaps.
	  using StoreTransformedLess = ::fsc::TransformedComparator<Key, ::std::less, StoreTransform>;


      template <typename V>
      using UniqueKeySetUtilityType = ::std::unordered_set<V, StoreTransformedFarmHash, StoreTransformedEqual>;

      // communication stuff...
      const mxx::comm& comm;

      // ============= local modifiers.  not directly accessible publically.  meant to be called via collective calls.

      // abstract declarations - need to access the local containers, therefore override in subclases.
      virtual void local_clear() = 0;
      virtual void local_reserve(size_t n) = 0;

      map_base(const mxx::comm& _comm) : comm(_comm) {}

    public:
      virtual ~map_base() {};


      // ================ data access functions
      virtual void to_vector(std::vector<std::pair<Key, T> > & result) const  = 0;
      virtual void keys(std::vector<Key> & result) const = 0;

      /// convert the map to a vector.
      virtual std::vector<std::pair<Key, T> > to_vector() const {
        std::vector<std::pair<Key, T> > result;
        this->to_vector(result);
        return result;
      }

      /// extract the keys of a map.
      virtual std::vector<Key> keys() const {
        std::vector<Key> result;
        this->keys(result);
        return result;
      }

      // =========== local accessors.  abstract methods since they need access to local containers.
      virtual bool local_empty() const = 0;
      virtual size_t local_size() const = 0;
      virtual size_t local_unique_size() const = 0;

      // =========== collective accessors

      /// check if empty.
      bool empty() const {
        if (comm.size() == 1)
          return this->local_empty();
        else // all reduce
          return mxx::all_of(this->local_empty(), comm);
      }

      /// get TOTAL size of distributed container
      size_t size() const {
        size_t s = this->local_size();

        if (comm.rank() == 0) printf("rank %d map_base local_size in size() is %lu\n", comm.rank(), s);

        if (comm.size() == 1)
          return s;
        else
          return ::mxx::allreduce(s, comm);
      }

      virtual size_t unique_size() const {
          size_t s = this->local_unique_size();
          if (comm.size() == 1)
            return s;
          else
            return ::mxx::allreduce(s, comm);
      }

      /// access the current the multiplicity.  only multimap needs to override this.
      virtual float get_multiplicity() const {
        // multimaps would add a collective function to change the multiplicity
        if (this->comm.rank() == 0) printf("rank %d map_base get_multiplicity called\n", this->comm.rank());

        return 1.0f;
      }

      // ============= collective modifiers

      /// reserve space.  n is the local container size.  this allows different processes to individually adjust its own size.
      virtual void reserve( size_t n) {
        // direct reserve + barrier
        this->local_reserve(n);
        if (this->comm.size() > 1) comm.barrier();
      }

      /// clears the distributed container.
      virtual void clear() {
        // clear + barrier.
        this->local_clear();
        if (comm.size() > 1)
          comm.barrier();
      }

      template <typename V>
      void transform_input(std::vector<V> & input) const {
    	  std::for_each(input.begin(), input.end(), InputTransform());
      }

  };

}



#endif /* SRC_CONTAINERS_DISTRIBUTED_MAP_BASE_HPP_ */
