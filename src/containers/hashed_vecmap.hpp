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
 * @file    hashed_vecmap.hpp
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
   *          This class attempts to address this shortcoming by maintaining the unordered map interface while replacing the LINKED_LIST with std::vector or std::multimap
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
      using subiter_type = typename ::std::vector<::std::pair<Key, T>, Allocator >::iterator;
      using const_subiter_type = typename ::std::vector<::std::pair<Key, T>, Allocator >::const_iterator;

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

      inline size_t distance(value_range_type const & iters) const {
        return std::distance(iters.first, iters.second);
      }
      inline size_t distance(Key const & k) const {
        return distance(map.at(k));
      }
      inline size_t distance(typename supercontainer_type::iterator map_iter) const {
        return distance(map_iter->second);
      }
      inline size_t distance(typename supercontainer_type::const_iterator map_iter) const {
        return distance(map_iter->second);
      }

      /**
       * @class    bliss::iterator::ConcatenatingIterator
       * @brief    this class presents a single/sequential view of a series of underlying iterator ranges
       * @details  random access iterator is not supported. other iterator categories are okay.
       *
       */
      template<typename V>
      class concat_iter :
        public ::std::iterator<
          typename ::std::random_access_iterator_tag,
          V
        >
      {
        protected:
          using subiterator_type = typename ::std::conditional<::std::is_const<V>::value,
              typename subcontainer_type::const_iterator, typename subcontainer_type::iterator>::type;
    	  using superiterator_type = typename ::std::conditional<::std::is_const<V>::value,
              typename supercontainer_type::const_iterator, typename supercontainer_type::iterator>::type;
          using type = concat_iter<V>;

          using inner_value_type = typename ::std::iterator_traits<subiterator_type>::value_type;

        public:
          template <typename KK, typename TT, typename HH, typename EE, typename AA, typename OutputIterator>
          OutputIterator
          copy(typename ::fsc::hashed_vecmap<KK, TT, HH, EE, AA>::template concat_iter<V> first,
               typename ::fsc::hashed_vecmap<KK, TT, HH, EE, AA>::template concat_iter<V> last, OutputIterator result);


        protected:
          /// the current position in the ranges list
          superiterator_type curr_iter;

          /// the current iterator position in the range of interest.
          subiterator_type curr_pos;

          superiterator_type max_iter;

          bool at_max;

          /// enforce that iterator is at a dereferenceable position
          void ensure_dereferenceable() {
            if (at_max) return;

            if (curr_iter == max_iter) {
              at_max = true;
              return;
            }

            // check to see if we are at end of subcontainer.  if so, move to next dereferenceable position.
            // end of a subcontainer is treated as same position as beginning of next subcontainer.
            while (curr_pos == curr_iter->second.second) {
              ++curr_iter;
              if (curr_iter == max_iter) {
                at_max = true;
                break; // reached the very end, can't deref max_iter.
              }
              else curr_pos = curr_iter->second.first;
            }
          }


        public:
          using difference_type = typename ::std::iterator_traits<subiterator_type>::difference_type;


          /// constructor for end concat iterator.  _end refers to end of supercontainer.
          concat_iter(superiterator_type _end) : curr_iter(_end), max_iter(_end), at_max(true) {};


          /// constructor for start concatenating iterator.  general version
          concat_iter(superiterator_type _iter, superiterator_type _end, subiterator_type _pos) :
            curr_iter(_iter), curr_pos(_pos), max_iter(_end), at_max(_iter == _end) {
            ensure_dereferenceable();
          };

          /// constructor for start concatenating iterator.  general version with checking that _pos belongs to _iter subcontainer via distance check.
          concat_iter(superiterator_type _iter, superiterator_type _end, subiterator_type _pos, difference_type distance_check) :
            concat_iter<V>(_iter, _end, _pos) {
        	  subiterator_type temp = _iter->second.first;
            if (!at_max && (distance_check != ::std::distance(temp, _pos) ) )
              throw std::logic_error("unordered_compact_vecmap constructor failing distance check, suggesting that _pos is not from same subcontainer as what _iter points to");
            ensure_dereferenceable();
          };


          // note that explicit keyword cannot be on copy and move constructors else the constructors are not defined/found.

          // copy constructor, assignment operator, move constructor, assignment operator should
          // should be default since the member vars are simple.

          bool is_at_max() const {
            at_max = (curr_iter == max_iter);
            return at_max;
          }

          /**
           * @brief increment:  move to the next position in the concatenating iterator, which may cross range boundaries.
           * @note  side effect: set at_end variable.
           * @return
           */
          type& operator++() {
            // if at end, return
            if (!at_max) {
              // now increment.  since we are careful to leave iterator at a dereferenceable state, curr_pos is not at a subcontainer's end.
              // so just increment.
              ++curr_pos;

              // now make sure we don't end up at a subcontainer's end.
              ensure_dereferenceable();
            }
            return *this;
          }


          /**
           * post increment.  make a copy then increment that.
           */
          type operator++(int)
          {
            type output(*this);
            this->operator++();
            return output;
          }

          //=== input iterator specific

          /// comparison operator
          bool operator==(const type& rhs) const
            {
            if (max_iter != rhs.max_iter) throw std::logic_error("the iterators being compared do not have the same internal end iterators so they are not comparable.");

            if (at_max && rhs.at_max) return true;
            if (at_max || rhs.at_max) return false;

            return ((curr_iter == rhs.curr_iter) && (curr_pos == rhs.curr_pos));
            }

          /// comparison operator
          bool operator!=(const type& rhs) const
            {
            if (max_iter != rhs.max_iter) throw std::logic_error("the iterators being compared do not have the same internal end iterators so they are not comparable.");

            if (at_max && rhs.at_max) return false;
            if (at_max || rhs.at_max) return true;

            return ((curr_iter != rhs.curr_iter) || (curr_pos != rhs.curr_pos));
            }


          template <
              typename VV = V,
              typename IV = inner_value_type,
              typename = typename ::std::enable_if<::std::is_constructible<VV, IV>::value>::type>
          inline V operator*() const {
            return *curr_pos;
          }

          /*=== NOT output iterator.  this is a map, does not make sense to change the  */
          /* content via iterator.                                                      */


          //=== NOT full forward iterator - no default constructor

          //=== NOT bidirectional iterator - no decrement because map produces forward iterator only.

          //=== NOT full random access iterator - only have +, += but not -. -=.  no comparison operators. have offset dereference operator [].


          /**
           * @brief     Advances this iterator by `n` positions.
           *            used by std::advance when randomaccess iterator.
           * @param n   The number of positions to advance.
           * @return    A reference to this after advancing.
           */
          type& operator+=(difference_type n)
          {
            // ::std::advance will use this or the ++ operator.
            if (n < 0) throw ::std::logic_error("::fsc::hashed_vecmap::iterator does not support decrement.");
            if (n == 0) return *this;  // nothing to add.
            if (at_max) return *this;  // iterator at the end.

            auto orig_iter = curr_iter;

            // dereferenceable right now
            subiterator_type temp = curr_iter->second.second;
            auto curr_dist = ::std::distance(curr_pos, temp);
            while (n >= curr_dist) {
              // not at end, and n is larger than curr dist, so go to next subcontainer.
              n -= curr_dist;  // consume some entries
              ++curr_iter;     // go to next container.
              at_max = (curr_iter == max_iter);
              if (at_max) return *this;  // if we are at end right now, then we can just return.
              else curr_dist = ::std::distance(curr_iter->second.first, curr_iter->second.second);    // see how much the next container has.

            }  // when exiting here, we are at a subcontainer that has more entries than n.  n could be 0.

            // now reset the curr_pos if curr_iter has been moved.
            if (curr_iter != orig_iter) curr_pos = curr_iter->second.first;
            ::std::advance(curr_pos, n);

            return *this;
          }


          /**
           * @brief     Advances a copy of this iterator by `n` positions.
           *
           * @param n   The number of positions to advance.
           * @return    The advanced iterator.
           */
          type operator+(difference_type n)
          {
            // reduced to += operator
            type output(*this);
            output += n;
            return output;
          }

          /**
           * @brief     Advances a copy of the `right` iterator by `n` positions.
           *
           * @param n   The number of positions to advance.
           * @return    The advanced iterator.
           */
          friend type operator+(difference_type n, const type& right)
          {
            // reduced to + operator
            return right + n;
          }

          /**
           * @brief     Returns the n'th element as seen from the current iterator
           *            position.
           *
           * @param n   The offset.
           *
           * @return    The element at offset `n` from the current position.
           */
          V operator[](difference_type n)
          {
            // reduce to the following:
            return *(*this + n);
          }

          /// difference between 2 iterators.  used by std::distance.
          friend difference_type operator-(const type& last, const type& first) {
            if (last == first) return 0;  // if both are at end, then we say dist is 0.
            if (first.at_max) return ::std::numeric_limits<difference_type>::lowest();  // first is at end, can't get to last.


            // now try to increment first until we get to last, or fail at getting to last.
            difference_type n = 0;
            auto cit = first.curr_iter;

            // init distance.  only meaningful here when first and last are not on same subcontainer.
            subiterator_type temp = cit->second.second;
            auto dist = std::distance(first.curr_pos, temp);

            // walk until either we are in same subcontainer, or at end of first iterator.
            while (cit != last.curr_iter) {
              n += dist;  //
              ++cit;
              if (cit == first.max_iter) break;
              else dist = ::std::distance(cit->second.first, cit->second.second);
            }

            // at this point, we have cit == last.curr_iter, or cit == eit (cit == eit == last_curr_iter possible)
            if (cit == first.max_iter) {  // cit = eit
              // if at end of first, but not at curr of last, not reachable.
              // else if at end of first, and at curr of last (== end), then reached.  return n.
              return ((cit != last.curr_iter) ? ::std::numeric_limits<difference_type>::lowest() : n);
            }

            // else we have cit == last.curr_iter.  dist is either size of cit subcontainer, or from curr_pos to subcontainer's end, both not correct.
            // need to recalculate dist now.  first move cpos

            // recalc distance. if cit hasn't moved, use original.  else use beginning of current subcontainer.
            temp = cit->second.first;
            dist = std::distance(((cit == first.curr_iter) ? first.curr_pos : temp), last.curr_pos);
            // if dist is negative, then last.curr_pos is before cpos.  not reachable.
            // else add the distance to running total.
            return ((dist < 0) ? ::std::numeric_limits<difference_type>::lowest() : (n + dist));


            //            // iterate until both are in the same subcontainer, or first is at end.
            //            while (!it.at_end() && (it.curr_iter != last.curr_iter)) {
            //              n += ::std::distance(it.curr_pos, it.curr_iter->second.second);
            //              ++(it.curr_iter);
            //              if (it.curr_iter != it.max_iter) it.curr_pos = (it.curr_iter)->second.first;
            //            }
            //
            //            // both are at same place (including end)
            //            if (it == last) return n;
            //            // if first is at its end, and first != last, then last is not reachable.
            //            if (it.at_end) return ::std::numeric_limits<difference_type>::lowest();
            //
            //            // first is not at end, and first != last, then last is in same container as first.
            //            n += ::std::distance(it.curr_pos, (last.curr_iter == last.max_iter) ? last.curr_iter->second.first : last.curr_pos);
            //
            //            return n;
          }
      };



      subcontainer_type vec;
      supercontainer_type map;
      size_t s;

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
          if ((map_it != map.end()) && (distance(map_it) > 0)) {

              it = map_it->second.first;  // use original range to jump past stuff.
              std::advance(it, distance(map_it) - 1);
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
        s = 0UL;

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
          s += std::distance(first, it);

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
      using iterator              = concat_iter<value_type>;
      using const_iterator        = concat_iter<const value_type>;
      using size_type             = typename subcontainer_type::size_type;
      using difference_type       = typename subcontainer_type::difference_type;


      //  if multiplicity of dataset is kind of known, can initialize data to that to avoid growing vector on average.
      hashed_vecmap(size_type load_factor = 1,
                   size_type bucket_count = 128,
                         const Hash& hash = Hash(),
                         const Less& less = Less(),
                         const Equal& equal = Equal(),
                         const Allocator& alloc = Allocator()) :
                           map(bucket_count, hash, equal, alloc), s(0UL) {};

      template<class InputIt>
      hashed_vecmap(InputIt first, InputIt last,
                         size_type load_factor = 1,
                         size_type bucket_count = 128,
                         const Hash& hash = Hash(),
                         const Less& less = Less(),
                         const Equal& equal = Equal(),
                         const Allocator& alloc = Allocator()) :
                         vec(first, last),
                         map(vec.size(), hash, equal, alloc), s(0UL) {

          std::sort(vec.begin(), vec.end(), less);

          this->rebuild();
      };

      virtual ~hashed_vecmap() {};



      iterator begin() {
        return iterator(map.begin(), map.end(), map.begin()->second.first, 0);
      }
      const_iterator begin() const {
        return cbegin();
      }
      const_iterator cbegin() const {
        return const_iterator(map.cbegin(), map.cend(), map.cbegin()->second.first, 0);
      }



      iterator end() {
        return iterator(map.end());
      }
      const_iterator end() const {
        return cend();
      }
      const_iterator cend() const {
        return const_iterator(map.cend());
      }

      bool empty() const {
        return map.empty();
      }

      size_type size() const {
        return s;
      }

      size_type unique_size() const {
        return map.size();
      }

      void clear() {
    	  vec.clear();
        map.clear();
        s = 0UL;
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
          vec.insert(vec.end(), first, last);
          auto middle = vec.begin() + prev_size;   // get the middle this way because clang 3.5 and gcc 4.8's STL implementation (libc++?)
                                                   // on Travis is not c++11 compliant and does not return the insertion point iterator.

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

        size_t count = 0;
//        bool erased = false;

        // mark for erasure
        auto middle = map.begin()->second.first;
        for (; first != last; ++first) {
          auto iter = map.find(*(first));
          if (iter == map.end()) continue;

          middle = ::std::partition(iter->second.first, iter->second.second, pred);

          count += std::distance(iter->second.first, middle);
          iter->second.first = middle;

          if (distance(iter) == 0) map.erase(iter);
//          erased = true;
        }
//        if (erased) this->inplace_compact();
        s -= count;
        return count;
      }

      template <typename InputIt>
      size_t erase(InputIt first, InputIt last) {
        static_assert(::std::is_convertible<Key, typename ::std::iterator_traits<InputIt>::value_type>::value,
                      "InputIt value type for erase cannot be converted to key type");


        if (first == last) return 0;

        size_t count = 0;

//        bool erased = false;

        // mark for erasure
        for (; first != last; ++first) {
          auto iter = map.find(*first);
          if (iter == map.end()) continue;

          count += distance(iter);
          //iter->second.first = iter->second.second;
          map.erase(iter);

//          erased = true;
        }

//        if (erased) this->inplace_compact();
        s -= count;
        return count;
      }

      template <typename Pred>
      size_t erase(Pred const & pred) {

        if (this->size() == 0) return 0;

        size_t before = this->size();

        auto new_end = ::std::remove_if(vec.begin(), vec.end(), pred);
        vec.erase(new_end, vec.end());

        this->rebuild();

        return before - this->size();
      }

      size_type count(Key const & key) const {
        if (map.find(key) == map.end()) return 0;
        else return distance(key);
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
          max_multiplicity = ::std::max(max_multiplicity, distance(it));
        }
        return max_multiplicity;
      }

      size_type get_min_multiplicity() const {
        size_type min_multiplicity = ::std::numeric_limits<size_type>::max();
        auto max = map.cend();
        size_type ss = 0;
        for (auto it = map.cbegin(); it != max; ++it) {
          ss = distance(it);
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
          key_s = distance(it);
          stdev_multiplicity += (key_s * key_s);
        }
        return stdev_multiplicity / double(map.size()) - get_mean_multiplicity();
      }


      ::std::pair<subiter_type, subiter_type> equal_range_value_only(Key const & key) {
        auto iter = map.find(key);

        if (iter == map.end()) return ::std::make_pair(subiter_type(), subiter_type());

        return ::std::make_pair(iter->second.first, iter->second.second);
      }
      ::std::pair<const_subiter_type, const_subiter_type> equal_range_value_only(Key const & key) const {
        auto iter = map.find(key);

        if (iter == map.end()) return ::std::make_pair(const_subiter_type(), const_subiter_type());

        return ::std::make_pair(iter->second.first, iter->second.second);

      }


      ::std::pair<iterator, iterator> equal_range(Key const & key) {
        auto iter = map.find(key);

        if (iter == map.end()) return ::std::make_pair(iterator(map.end()), iterator(map.end()));

        return ::std::make_pair(iterator(iter, map.end(), iter->second.first, 0),
                                iterator(iter, map.end(), iter->second.second, distance(iter->second)));
      }
      ::std::pair<const_iterator, const_iterator> equal_range(Key const & key) const {
        auto iter = map.find(key);

        if (iter == map.cend()) return ::std::make_pair(const_iterator(map.cend()), const_iterator(map.cend()));

        return ::std::make_pair(const_iterator(iter, map.cend(), iter->second.first, 0),
                                const_iterator(iter, map.cend(), iter->second.second, distance(iter->second)));

      }
      // NO bucket interfaces

  };



} // end namespace fsc.

namespace std {

  template <typename Key, typename T, typename Hash, typename Equal, typename Allocator, typename OutputIterator>
  OutputIterator
  copy(typename ::fsc::hashed_vecmap<Key, T, Hash, Equal, Allocator>::iterator first,
       typename ::fsc::hashed_vecmap<Key, T, Hash, Equal, Allocator>::iterator last, OutputIterator result) {

    // can last be reach from first?
    if ((last - first) <= 0) return result;

    // reachable.  so now walk.  do not need to do as much checking.
    auto out_iter = result;

    // now try to increment first until we get to last, or fail at getting to last.
    auto cit = first.curr_iter;
    auto cpos = first.curr_pos;

    // walk until either we are in same subcontainer.
    // since last is reachable from first, we don't need to check first's end_iter.
    while (cit != last.curr_iter) {
      out_iter = ::std::copy(cpos, cit->second.second, out_iter);
      ++cit;
      cpos = cit->second.first;
    }

    // now we're in same container.  take care the rest.  cpos was updated if cit was moved.
    out_iter = ::std::copy(cpos, last.curr_pos, out_iter);

    return out_iter;
  }

  template <typename Key, typename T, typename Hash, typename Equal, typename Allocator, typename OutputIterator>
  OutputIterator
  copy(typename ::fsc::hashed_vecmap<Key, T, Hash, Equal, Allocator>::const_iterator first,
       typename ::fsc::hashed_vecmap<Key, T, Hash, Equal, Allocator>::const_iterator last, OutputIterator result) {

    // can last be reach from first?
    if ((last - first) <= 0) return result;

    // reachable.  so now walk.  do not need to do as much checking.
    auto out_iter = result;

    // now try to increment first until we get to last, or fail at getting to last.
    auto cit = first.curr_iter;
    auto cpos = first.curr_pos;

    // walk until either we are in same subcontainer.
    // since last is reachable from first, we don't need to check first's end_iter.
    while (cit != last.curr_iter) {
      out_iter = ::std::copy(cpos, cit->second.second, out_iter);
      ++cit;
      cpos = cit->second.first;
    }

    // now we're in same container.  take care the rest.  cpos was updated if cit was moved.
    out_iter = ::std::copy(cpos, last.curr_pos, out_iter);


    return out_iter;
  }

}




#endif /* SRC_WIP_HASHED_VECMAP_HPP_ */
