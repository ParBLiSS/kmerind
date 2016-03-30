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
      using TransformedComp = ::fsc::TransformedComparator<Key, Comparator, KeyTransform>;

      using TransformedFarmHash = ::fsc::TransformedHash<Key, ::bliss::kmer::hash::farm<Key, false>, KeyTransform>;
      static TransformedFarmHash local_hash;

      using TransformedEqual = TransformedComp<Equal>;
      static TransformedEqual equal;

      using TransformedLess = TransformedComp<Less>;
      static TransformedLess less;

      template <typename V>
      using UniqueKeySetUtilityType = ::std::unordered_set<V, TransformedFarmHash, TransformedEqual>;

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


  template<typename Key, typename T,
      class Comm,
      template <typename> class KeyTransform,
      class Less,
      class Equal,
      class Alloc>
  typename map_base<Key, T, Comm, KeyTransform, Less, Equal, Alloc>::TransformedFarmHash
      map_base<Key, T, Comm, KeyTransform, Less, Equal, Alloc>::local_hash;

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
