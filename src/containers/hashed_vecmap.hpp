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
 * @file    unordered_vecmap.hpp
 * @ingroup
 * @author  tpan
 * @brief
 * @details
 *
 */
#ifndef SRC_WIP_HASHED_VECMAP_HPP_
#define SRC_WIP_HASHED_VECMAP_HPP_

#include <unordered_map>
#include <vector>
#include <functional>  // hash, equal_to, etc
#include <tuple>   // pair
#include <scoped_allocator>
#include <algorithm>
#include <cmath>   // ceil

#include "utils/logging.h"

namespace fsc {  // fast standard container

  /**
   * @brief my version of hashed map.  because std::unordered_map does chaining for collision using a LINKED LIST.
   * @details std::unordered_map's use of linked list is fine for low collision.  for high collision, it's great for insertion but terrible for search or delete
   *          it is also not suitable for copying to a vector and sort/search.  in addition, the search (equal_range, count, and probably erase) are implemented
   *          using linear search.
   *
   *          This class attemps to address this shortcoming by maitaining the unordered map interface while replacing the LINKED_LIST with std::vector or std::multimap
   *
   *          Also note: google dense hash has strong requirements for the hash object to be move constructable/assignable (in our case, TransformedHash, and the underlying hash functions.
   *            because of the special constructors requirements, copy constructor and default constructors also need to be defined.
   *          Performance:  with high repeat input set (e.g. test.fastq, with 254K repeats per kmer, 40 unique kmers) std::unordered_map search time is approximately 1.5 sec for 1%, count is 0.56 sec.  build is 1.6 s.  pos+qual map)
   *          	hashmap:  build 1.6s, find 1.5s, count 0.56s
   *          	sorted vector:  build = 2 s, find 0.27s, count 0.01s
   *          with google dense_hash_map:  can't use it - google dense_hash_map requires uniqueness, and also requiers a NULL element.
   *
   *          tested using std::multimap and map.  VERY slow on insertion, especially for real data.  about 33% faster for find and count compared to unordered map for synthetic, high repeat datasets.  still about 4 times slower than sort based.
   *          for that dataset.  build = 3.6s, find 0.98s, count 0.73s
   *
   *          WE HAVE TO BUILD OUR OWN VERSION.
   *          	prefer using std::vector as internal - self growing,
   *
   *          Note that this class has an incomplete implementation of the multimap interface - only those interfaces that I needed are implemented here.
   *
   *          DESIGN:
   *          instead of requiring sorting, since unordered_map (without multimap) has reasonable performance, we can
   *          build a multimap using unordered map by setting the mapped_type to be ::std::vector<::std::pair<Key, T> >
   *
   *          instead of building this at the distributed map level, we build it at the local unordered map level,
   *          which is easier to debug and benchmark.
   *
   * 		  vector holds pair<Key, T> because then we can directly copy (fast)
   * 		  a version with pair<T> would be more space efficient but slower when large amount of data is present.
   *
   *          this class internally uses a hashmap of vectors.
   *
   *
   *      memory usage:
   *        supercontainer = map, stores ::std::pair<K, std::vector<> > in linked list for chaining.
   *          each "bucket" hash size_t + pointer to head of linked list. - 16 bytes, HU hash unique elements.
   *          each link list node has ::std::pair<K, std::vector> as payload, and at least a next ptr. - 24 bytes,  U unique elements
   *          each vector stores  std::pair<K, T>, so 16 or 24 bytes.  N elements
   *
   *          total: (16N or 24N) + 24U + 16 HU.  assume good hash function, HU = U.
   *
   *          bottomline - large amount of memory is needed.
   *
   *  // no compact vec - can't effectively sort by key when key is not there.
   */

  /**
   * uncompacted version of the vecmap
   *
   * internally has a single vector, sorted by (hash(key) % buckets) then by hash(key) then by key
   * insert into back of vector.   when insert, when map load_factor may trigger a rehash, we manually resort the vector and rebuild the map.
   *
   * note that as we are attempting to adhere to unordered_multimap template interface, we can't really add a comparacter parameter without causing
   * problems elsewhere, namely in distributed unordered hashvec map.  the transform that is inherently in the hash needs to be applied, else we might
   * get into a situation where std::less results in 2 separate keys while hash lumps the keys together.
   *
   * we also cannot sort by hash value, because of possibility of collision.
   */
  template <typename Key,
  typename T,
  typename Hash = ::std::hash<Key>,
  typename Comparator = ::std::less<Key>,
  typename Equal = ::std::equal_to<Key>,
  typename Allocator = ::std::allocator<::std::pair<Key, T> > >
  class hashed_vecmap {

    protected:
      struct Less {
        Comparator l;

        inline bool operator()(Key const &x, Key const &y ) {
          return l(x, y);
        }

        template <typename V>
        inline bool operator()(::std::pair<Key, V> const & x, ::std::pair<Key, V> const & y) {
          return l(x.first, y.first);
        }
        template <typename V>
        inline bool operator()(::std::pair<const Key, V> const & x, ::std::pair<const Key, V> const & y) {
          return l(x.first, y.first);
        }
      };


      // data container
      using subcontainer_type = ::std::vector<::std::pair<Key, T>, Allocator >;

      // index "pointers"
      using value_range_type = ::std::pair<typename subcontainer_type::iterator,
           typename subcontainer_type::iterator>;

      using superallocator_type = ::std::allocator<::std::pair<const Key, value_range_type > >;
      using supercontainer_type =
          ::std::unordered_map<Key, value_range_type,
                         Hash, Equal, superallocator_type >;


      // use of vector - increased cost during construction but query will be fast.
      // group all entries with the same key with vector - list chaing but already sorted.
      // group all entries with same hash with vector - list chaining but with randomly accessible container.


      subcontainer_type vec;
      supercontainer_type map;

      inline size_t dist(Key k) const {
        return std::distance(map.at(k).first, map.at(k).second);
      }
      inline size_t dist(typename supercontainer_type::iterator map_iter) const {
        return std::distance(map_iter->second.first, map_iter->second.second);
      }
      inline size_t dist(typename supercontainer_type::const_iterator map_iter) const {
        return std::distance(map_iter->second.first, map_iter->second.second);
      }
      inline size_t dist(value_range_type iters) const {
        return std::distance(iters.first, iters.second);
      }

      /// compact the vector after disjoint entries are deleted.
      void compact() {
        if (vec.size() == 0 || map.size() == 0) {
          map.clear();
          vec.clear();
          return;
        }

        // allocate a new vector
        subcontainer_type temp;
        temp.resize(vec.size());
        auto curr = temp.begin();
        auto prev = map.begin();
        // for each map entry, copy content over.  the range iterators in map are moved to temp
        for (auto it = map.begin(), max = map.end(); it != max; ++it) {
          // copy over the content
          it->second.second = std::copy(it->second.first, it->second.second, curr);

          // update the iterators
          it->second.first = curr;
          curr = it->second.second;

          prev = it;
        }
        // erase everything after.  curr will be pointing to temp.end() now.
        temp.erase(curr, temp.end());

        // swap
        vec.swap(temp);

        // now last iter's second (was point to the last valid element, which then became the end of temp, and now needs to be come vec's end.
        prev->second.second = vec.end();
        // (all other iterators are not invalided.  just one that does not point to a real element (i.e. end iterator)

      }

      /// compact the vector after disjoint entries are deleted.  ASSUMPTION: all entries with same key are contiguous in memory, and map points to a subrange of each of these ranges.
      void inplace_compact() {
        if (vec.size() == 0 || map.size() == 0) {
            map.clear();
            vec.clear();
            return;
          }

        Less less;

        auto compacted_it = vec.begin();
        auto key = compacted_it->first;
        auto map_it = map.find(key);


        size_t i = 0;

        // go through all elements in vec, copy to beginning,
        for (auto it = vec.begin(), max = vec.end(); it != max; ++i) {
          // get the current key
          key = it->first;
          // find the map entry
          map_it = map.find(key);


          // if in map, then compact
          if ((map_it != map.end()) && (dist(map_it) > 0)) {

              it = map_it->second.first;  // use original range to jump past stuff.
              std::advance(it, dist(map_it) - 1);
              it = std::adjacent_find(it, max, less);

              // in map.  we copy over the entries to keep.  advances the target it
              map_it->second.second = std::copy(map_it->second.first, map_it->second.second, compacted_it);
              // update the map
              map_it->second.first = compacted_it;
              compacted_it = map_it->second.second;

              //printf("iter %ld curr valid count = %ld, curr source count = %ld\n", i, std::distance(vec.begin(), compacted_it), std::distance(vec.begin(), it));

          } // else not in map.  don't update map, don't update current target it, but do update current it.

          // but advance it.  sorted, so != is same as <
          it = std::adjacent_find(it, max, less);
          // we are now at the last of the identical.  loop will increment to next.

          if (it != max) ++it;  // if nothing is found, max is returned, so need to check that.
        }
        // erase the extra
        vec.erase(compacted_it, vec.end());

        // the last map entry needs to have its second point to vec.end() now.  all other iterators are still valid.
        map_it->second.second = vec.end();

      }



      /// rehash to rebuild the hashmap index.
      void rebuild() {
        map.clear();

        if (vec.size() == 0) return;

        map.reserve(vec.size());
        Less less;

        auto first = vec.begin();
        auto key = first->first;
        for (auto it = vec.begin(), max = vec.end(); it != max;) {
          first = it;
          key = first->first;

          // find the last of the entries with same key
          it = std::adjacent_find(it, max, less);

          if (it != max) {
            // not last entry, so advance 1.
            ++it;
          }
          map.emplace(std::move(key), std::move(std::make_pair(first, it)));

        }

//        size_t bucket_max = 0;
//        for (size_t i = 0; i < map.bucket_count(); ++i) {
//          bucket_max = ::std::max(bucket_max, map.bucket_size(i));
//        }
//        printf("map size: %ld, map buckets %ld, max_bucket %ld, map loadfactor %f\n", map.size(), map.bucket_count(), bucket_max, map.load_factor());
      }




    public:
      using key_type              = Key;
      using mapped_type           = T;
      using value_type            = ::std::pair<const Key, T>;
      using hasher                = Hash;
      using key_equal             = Equal;
      using allocator_type        = Allocator;
      using reference             = value_type&;
      using const_reference       = const value_type&;
      using pointer               = typename std::allocator_traits<Allocator>::pointer;
      using const_pointer         = typename std::allocator_traits<Allocator>::const_pointer;
      using iterator              = typename subcontainer_type::iterator;
      using const_iterator        = typename subcontainer_type::const_iterator;
      using size_type             = typename subcontainer_type::size_type;
      using difference_type       = typename subcontainer_type::difference_type;


      //  if multiplicity of dataset is kind of known, can initialize data to that to avoid growing vector on average.
      hashed_vecmap(size_type load_factor = 1,
                   size_type bucket_count = 128,
                         const Hash& hash = Hash(),
                         const Less& less = Less(),
                         const Equal& equal = Equal(),
                         const Allocator& alloc = Allocator()) :
                           map(bucket_count, hash, equal, alloc) {};

      template<class InputIt>
      hashed_vecmap(InputIt first, InputIt last,
                         size_type load_factor = 1,
                         size_type bucket_count = 128,
                         const Hash& hash = Hash(),
                         const Less& less = Less(),
                         const Equal& equal = Equal(),
                         const Allocator& alloc = Allocator()) :
                         vec(first, last),
                         map(vec.size(), hash, equal, alloc) {

          std::sort(vec.begin(), vec.end(), less);

          this->rebuild();
      };

      virtual ~hashed_vecmap() {};


      iterator begin() {
        return vec.begin();
      }
      const_iterator begin() const {
        return vec.cbegin();
      }
      const_iterator cbegin() const {
        return vec.cbegin();
      }

      iterator end() {
        return vec.end();
      }
      const_iterator end() const {
        return vec.cend();
      }
      const_iterator cend() const {
        return vec.cend();
      }

      bool empty() const {
        return map.empty();
      }

      size_type size() const {
        return vec.size();
      }

      size_type unique_size() const {
        return map.size();
      }

      void clear() {
    	  vec.clear();
        map.clear();
      }


      /// rehash for new count number of BUCKETS.  iterators are invalidated.
      void rehash(size_type count) {
        // only rehash if new bucket count is greater than old bucket count
        if (count > map.bucket_count())
          map.rehash(count);
      }

      /// bucket count.  same as underlying buckets
      size_type bucket_count() { return map.bucket_count(); }

      /// max load factor.  this is the map's max load factor (vectors per bucket) x multiplicity = elements per bucket.  side effect is multiplicity is updated.
      float max_load_factor() {
        return (map.size() == 0) ? map.max_load_factor() : map.max_load_factor() * (static_cast<float>(vec.size()) / static_cast<float>(map.size()));
      }



      // choices:  sort first, then insert in ranges, or no sort, insert one by one.  second is O(n) but pays the random access and mem realloc cost
      template <class InputIt>
      void insert(InputIt first, InputIt last) {
          static_assert(::std::is_convertible<std::pair<Key, T>,
                        typename ::std::iterator_traits<InputIt>::value_type>::value,
                        "InputIt value type for insert cannot be converted to std::pair<Key, T> type");

          if (first == last) return;

          size_t prev_size = vec.size();

          // copy it in.
          vec.reserve(prev_size + std::distance(first, last));
          auto middle = vec.insert(vec.end(), first, last);

          Less less;
          // sort the new part
          std::sort(middle, vec.end(), less);

          // merge with previous
          if (prev_size > 0) ::std::inplace_merge(vec.begin(), middle, vec.end(), less);

          // rebuild index
          this->rebuild();
      }

      /// inserting sorted range
      void insert(::std::vector<::std::pair<Key, T> > & input) {

        if (input.size() == 0) return;

        if (vec.empty()) {
          vec.swap(input);

          Less less;
          // sort the new part
          std::sort(vec.begin(), vec.end(), less);

          this->rebuild();
        }
        else {
          this->insert(input.begin(), input.end());
        }
      }

      template <typename InputIt, typename Pred>
      size_t erase(InputIt first, InputIt last, Pred const & pred) {
        static_assert(::std::is_convertible<Key, typename ::std::iterator_traits<InputIt>::value_type>::value,
                      "InputIt value type for erase cannot be converted to key type");

        if (first == last) return 0;

        size_t before = this->size();

        bool erased = false;

        // mark for erasure
        for (; first != last; ++first) {
          auto iter = map.find(*(first));
          if (iter == map.end()) continue;

          iter->second.first = ::std::partition(iter->second.first, iter->second.second, pred);

          if (dist(iter) == 0) map.erase(iter);

          erased = true;
        }
        if (erased) this->inplace_compact();

    	  return before - this->size();
      }

      template <typename InputIt>
      size_t erase(InputIt first, InputIt last) {
        static_assert(::std::is_convertible<Key, typename ::std::iterator_traits<InputIt>::value_type>::value,
                      "InputIt value type for erase cannot be converted to key type");


        if (first == last) return 0;

        size_t before = this->size();

        bool erased = false;

        size_t j = 0;

        // mark for erasure
        for (; first != last; ++first, ++j) {
          auto iter = map.find(*first);
          if (iter == map.end()) continue;

          //iter->second.first = iter->second.second;
          map.erase(iter);

          erased = true;
        }

        if (erased) this->inplace_compact();

        return before - this->size();
      }

      template <typename Pred>
      size_t erase(Pred const & pred) {

        if (this->size() == 0) return 0;

        size_t before = this->size();

        bool erased = false;

        auto new_end = ::std::remove_if(vec.begin(), vec.end(), pred);
        vec.erase(new_end, vec.end());

        this->rebuild();

        return before - this->size();
      }

      size_type count(Key const & key) const {
        if (map.find(key) == map.end()) return 0;
        else return dist(key);
      }



      void report() {
          BL_INFOF("vecmap bucket count: %lu\n", map.bucket_count());
          BL_INFOF("vecmap load factor: %f\n", map.load_factor());
          BL_INFOF("vecmap unique entries: %lu\n", map.size());
          BL_INFOF("vecmap total size: %lu\n", s);
      }


      size_type get_max_multiplicity() const {
        size_type max_multiplicity = 0;
        auto max = map.cend();
        for (auto it = map.cbegin(); it != max; ++it) {
          max_multiplicity = ::std::max(max_multiplicity, dist(it));
        }
        return max_multiplicity;
      }

      size_type get_min_multiplicity() const {
        size_type min_multiplicity = ::std::numeric_limits<size_type>::max();
        auto max = map.cend();
        size_type ss = 0;
        for (auto it = map.cbegin(); it != max; ++it) {
          ss = dist(it);
          if (ss > 0)
            min_multiplicity = ::std::min(min_multiplicity, ss);
        }
        return min_multiplicity;
      }

      double get_mean_multiplicity() const {
        return static_cast<double>(vec.size()) / double(map.size());
      }
      double get_stdev_multiplicity() const {
        double stdev_multiplicity = 0;
        auto max = map.cend();
        double key_s;
        for (auto it = map.cbegin(); it != max; ++it) {
          key_s = dist(it);
          stdev_multiplicity += (key_s * key_s);
        }
        return stdev_multiplicity / double(map.size()) - get_mean_multiplicity();
      }


      ::std::pair<iterator, iterator> equal_range(Key const & key) {
        auto iter = map.find(key);

        if (iter == map.end()) return ::std::make_pair(vec.end(), vec.end());

        return iter->second;
      }
      ::std::pair<const_iterator, const_iterator> equal_range(Key const & key) const {
        auto iter = map.find(key);

        if (iter == map.cend()) return ::std::make_pair(vec.cend(), vec.cend());

        return iter->second;

      }
      // NO bucket interfaces

  };



} // end namespace fsc.





#endif /* SRC_WIP_HASHED_VECMAP_HPP_ */
