/**
 * @file    unordered_vecmap.hpp
 * @ingroup
 * @author  tpan
 * @brief
 * @details
 *
 * Copyright (c) 2015 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef SRC_WIP_UNORDERED_VECMAP_HPP_
#define SRC_WIP_UNORDERED_VECMAP_HPP_

#include <unordered_map>
#include <vector>
#include <functional>  // hash, equal_to, etc
#include <tuple>   // pair
#include <scoped_allocator>
#include <algorithm>
#include <cmath>   // ceil


namespace fsc {  // fast standard container

  /**
   * @brief my version of unordered map.  because std::unordered_map does chaining for collision using a LINKED LIST.
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
   */
  template <typename Key,
  typename T,
  typename Hash = ::std::hash<Key>,
  typename Equal = ::std::equal_to<Key>,
  typename Allocator =
      ::std::scoped_allocator_adaptor<::std::allocator<::std::pair<const Key, ::std::vector<T, ::std::allocator<T> > > >,
       ::std::allocator<T> >
  >
  class unordered_compact_vecmap {

    protected:
      using subcontainer_type = ::std::vector<T, ::std::allocator<T> >;
      using supercontainer_type =
          ::std::unordered_map<Key, subcontainer_type, Hash, Equal, Allocator >;

      using subiter_type = typename subcontainer_type::iterator;
      using const_subiter_type = typename subcontainer_type::const_iterator;

      // use of vector - increased cost during construction but query will be fast.
      // group all entries with the same key with vector - list chaing but already sorted.
      // group all entries with same hash with vector - list chaining but with randomly accessible container.

      //===== create special iterator for this class - dereferences iterator to internal_map, walks through
      // vector's iterator one by one via ++ and dereference operator.
      // acts like a concatenating iterator, but source list of component iterators is from the internal map's set of vectors.

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

        public:
          template <typename KK, typename TT, typename HH, typename EE, typename AA, typename OutputIterator>
          OutputIterator
          copy(typename ::fsc::unordered_compact_vecmap<KK, TT, HH, EE, AA>::template concat_iter<V> first,
               typename ::fsc::unordered_compact_vecmap<KK, TT, HH, EE, AA>::template concat_iter<V> last, OutputIterator result);


        protected:
          /// the current position in the range
          superiterator_type curr_iter;
          subiterator_type curr_pos;

          /// parent container to help with identifying end iterator.
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
            while (curr_pos == curr_iter->second.end()) {
              ++curr_iter;
              if (curr_iter == max_iter) {
                at_max = true;
                break; // reached the very end, can't deref max_iter.
              }
              else curr_pos = curr_iter->second.begin();
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
            if (!at_max && (distance_check != ::std::distance(_iter->second.begin(), _pos) ) )
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

          /// get pointer.  since we have to construct the pair during return, pointer does not make sense.
          //          inline V* operator->() const
          //          {
          //        	// TO FIX
          //            return &(*curr_pos);
          //          }

          /// dereference operator.  normal iterator: returns a rvalue
          inline V operator*() const
          {
            return V(curr_iter->first, *curr_pos);
          }

          //=== NOT output iterator.  this is a map, does not make sense to change the content via iterator.

          //          /// dereference operator.  return lvalue  for non-const item
          //          template<typename VV = V>
          //          inline typename std::enable_if<!std::is_const<VV>::value, VV& >::type operator*()
          //          {
          //            return *curr_pos;
          //          }

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
            if (n < 0) throw ::std::logic_error("::fsc::unordered_compact_vecmap::iterator does not support decrement.");
            if (n == 0) return *this;  // nothing to add.
            if (at_max) return *this;  // iterator at the end.

            auto orig_iter = curr_iter;

            // dereferenceable right now
            auto curr_dist = ::std::distance(curr_pos, curr_iter->second.end());
            while (n >= curr_dist) {
              // not at end, and n is larger than curr dist, so go to next subcontainer.
              n -= curr_dist;  // consume some entries
              ++curr_iter;  	 // go to next container.
              at_max = (curr_iter == max_iter);
              if (at_max) return *this;  // if we are at end right now, then we can just return.
              else curr_dist = curr_iter->second.size();	  // see how much the next container has.

            }  // when exiting here, we are at a subcontainer that has more entries than n.  n could be 0.

            // now reset the curr_pos if curr_iter has been moved.
            if (curr_iter != orig_iter) curr_pos = curr_iter->second.begin();
            ::std::advance(curr_pos, n);

            //
            //
            //            while ((!at_end) && (n > 0)) {
            //              curr_dist = ::std::distance(curr_pos, curr_iter->second.end());
            //
            //              if (curr_dist > n) {
            //                ::std::advance(curr_pos, n); // curr_pos would not be at subcontainer end, so derefernceable, and not at end.
            //                break;
            //              } else if (curr_dist == n) {
            //            	  // n goes to 0, curr_pos should go to next subcontainer's start if possible,
            //            	 ++curr_iter;
            //            	 at_end = (curr_iter == end_iter);
            //            	 if (!at_end) curr_pos = curr_iter->second.begin();
            //            	 break;
            //              } else { // if (curr_dist < n)
            //                n -= curr_dist;
            //                // go to next subcontainer's beginning.
            //				 ++curr_iter;
            //				 at_end = (curr_iter == end_iter);
            //				 if (!at_end) curr_pos = curr_iter->second.begin();
            //              }
            //            }
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
            auto dist = cit->second.end() - first.curr_pos;

            // walk until either we are in same subcontainer, or at end of first iterator.
            while (cit != last.curr_iter) {
              n += dist;  //
              ++cit;
              if (cit == first.max_iter) break;
              else dist = cit->second.size();
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
            dist = last.curr_pos - ((cit == first.curr_iter) ? first.curr_pos : cit->second.begin());
            // if dist is negative, then last.curr_pos is before cpos.  not reachable.
            // else add the distance to running total.
            return ((dist < 0) ? ::std::numeric_limits<difference_type>::lowest() : (n + dist));


            //            // iterate until both are in the same subcontainer, or first is at end.
            //            while (!it.at_end() && (it.curr_iter != last.curr_iter)) {
            //              n += ::std::distance(it.curr_pos, it.curr_iter->second.end());
            //              ++(it.curr_iter);
            //              if (it.curr_iter != it.end_iter) it.curr_pos = (it.curr_iter)->second.begin();
            //            }
            //
            //            // both are at same place (including end)
            //            if (it == last) return n;
            //            // if first is at its end, and first != last, then last is not reachable.
            //            if (it.at_end) return ::std::numeric_limits<difference_type>::lowest();
            //
            //            // first is not at end, and first != last, then last is in same container as first.
            //            n += ::std::distance(it.curr_pos, (last.curr_iter == last.end_iter) ? last.curr_iter->second.begin() : last.curr_pos);
            //
            //            return n;
          }
      };


      supercontainer_type map;
      size_t s;
      size_t multiplicity;


    public:
      using key_type              = Key;
      using mapped_type           = T;
      using value_type            = ::std::pair<const Key, T>;
      using hasher                = Hash;
      using key_equal             = Equal;
      using allocator_type        = Allocator;
      using reference 			  = value_type&;
      using const_reference	      = const value_type&;
      using pointer				  = typename std::allocator_traits<Allocator>::pointer;
      using const_pointer		  = typename std::allocator_traits<Allocator>::const_pointer;
      using iterator              = concat_iter<value_type>;
      using const_iterator        = concat_iter<const value_type>;
      using size_type             = typename subcontainer_type::size_type;
      using difference_type       = typename subcontainer_type::difference_type;


      //  if multiplicity of dataset is kind of known, can initialize data to that to avoid growing vector on average.
      unordered_compact_vecmap(size_type load_factor = 1,
                       size_type bucket_count = 128,
                       const Hash& hash = Hash(),
                       const Equal& equal = Equal(),
                       const Allocator& alloc = Allocator()) :
                         map(bucket_count, hash, equal, alloc),
                         s(0UL),
                         multiplicity(load_factor) {};

      template<class InputIt>
      unordered_compact_vecmap(InputIt first, InputIt last,
                       size_type load_factor = 1,
                       size_type bucket_count = 128,
                       const Hash& hash = Hash(),
                       const Equal& equal = Equal(),
                       const Allocator& alloc = Allocator()) :
                       ::fsc::unordered_compact_vecmap<Key, T, Hash, Equal, Allocator>(load_factor,
                                                                               ::std::max(bucket_count, static_cast<size_type>(::std::distance(first, last))),
                                                                                hash, equal, alloc) {
          this->insert(first, last);
      };

      virtual ~unordered_compact_vecmap() {};


      iterator begin() {
        return iterator(map.begin(), map.end(), map.begin()->second.begin(), 0);
      }
      const_iterator begin() const {
        return cbegin();
      }
      const_iterator cbegin() const {
        return const_iterator(map.cbegin(), map.cend(), map.cbegin()->second.cbegin(), 0);
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

      void clear() {
        map.clear();
      }

      /// rehash for new count number of BUCKETS.  iterators are invalidated.  side effect is multiplicity is updated.
      void rehash(size_type count) {
        // compute current average number of entries per vector.  this is side effect.
        multiplicity = (s + map.size() - 1) / map.size();

        // only rehash if new bucket count is greater than old bucket count
        if (count > map.bucket_count())
          map.rehash(count);
      }

      /// bucket count.  same as underlying buckets
      size_type bucket_count() { return map.bucket_count(); }

      /// max load factor.  this is the map's max load factor (vectors per bucket) x multiplicity = elements per bucket.  side effect is multiplicity is updated.
      float max_load_factor() {
        // compute current average number of entries per vector.  this is side effect.
        multiplicity = (s + map.size() - 1) / map.size();

        return map.max_load_factor() * (static_cast<float>(s) / static_cast<float>(map.size()));
      }


      /// reserve for new count of elements.  iterators may be invalidated.
      void reserve(size_type count) {
        // compute current average number of entries per vector.  this is side effect
        multiplicity = (s + map.size() - 1) / map.size();

        // compute number of buckets required.
        this->rehash(std::ceil(static_cast<float>(count) / this->max_load_factor()));
      }

      iterator insert(const value_type & value) {
        return emplace(::std::forward<value_type>(value));
      }
      iterator insert(value_type && value) {
        return emplace(::std::forward<value_type>(value));
      }
      iterator emplace(value_type && value) {
        auto key = value.first;
        auto iter = map.find(key);
        if (iter == map.end()) {
          iter = map.emplace(key, subcontainer_type()).first;
          //          map.at(key).reserve(multiplicity);
        }
        iter->second.emplace_back(::std::forward<T>(value.second));
        auto pos = iter->second.end();
        ++s;
        return iterator(iter, map.end(), --pos, iter->second.size() - 1);
      }
      iterator emplace(Key&& key, T&& value) {
        auto iter = map.find(key);
        if (iter == map.end()) {
          iter = map.emplace(key, subcontainer_type()).first;
          //          map[key].reserve(multiplicity);
        }
        iter->second.emplace_back(::std::forward<T>(value));
        auto pos = iter->second.end();
        ++s;
        return iterator(iter, map.end(), --pos, iter->second.size() - 1);
      }

      // choices:  sort first, then insert in ranges, or no sort, insert one by one.  second is O(n)
      template <class InputIt>
      void insert(InputIt first, InputIt last) {
          for (; first != last; ++first) {
            //            if (map.find(first->first) == map.end()) {
            //              // emplace, and take the result iterator, get the second component (subcontainer), reserve the size.
            //              map.emplace(first->first, subcontainer_type());
            // //              map[k].reserve(multiplicity);
            //            }
            //            map.at(first->first).emplace_back(::std::forward<value_type>(*first));
            map[first->first].emplace_back(::std::forward<T>(first->second));
            ++s;
          }
      }


      iterator erase(const key_type& key) {
        s -= count(key);
        auto iter = map.erase(key);
        return iterator(iter, map.end(), iter->second.begin(), 0);
      }

      size_type count(Key const & key) const {
        if (map.find(key) == map.end()) return 0;
        else return map.at(key).size();
      }

      void shrink_to_fit() {
        // update multiplicity.
        this->multiplicity = (s + map.size() - 1) / map.size();

        // for each vector, shrink it.
        for (auto it = map.begin(), max = map.end(); it != max; ++it) {
          it->second.shrink_to_fit();
        }
      }

      void report() {
        INFOF("compact vecmap bucket count: %lu", map.bucket_count());
        INFOF("compact vecmap load factor: %f", map.load_factor());
        INFOF("compact vecmap unique entries: %lu", map.size());
        INFOF("compact vecmap total entries: %lu", this->s);
      }


      size_type unique_size() const {
        return map.size();
      }
      size_type get_max_multiplicity() const {
        size_type max_multiplicity = 0;
        for (auto it = map.cbegin(), max = map.cend(); it != max; ++it) {
          max_multiplicity = ::std::max(max_multiplicity, it->second.size());
        }
        return max_multiplicity;
      }
      size_type get_min_multiplicity() const {
        size_type min_multiplicity = ::std::numeric_limits<size_type>::max();
        for (auto it = map.cbegin(), max = map.cend(); it != max; ++it) {
          min_multiplicity = ::std::min(min_multiplicity, it->second.size());
        }
        return min_multiplicity;
      }
      double get_mean_multiplicity() const {
        return static_cast<double>(s) / double(map.size());
      }
      double get_stdev_multiplicity() const {
        double stdev_multiplicity = 0;
        for (auto it = map.cbegin(), max = map.cend(); it != max; ++it) {
          stdev_multiplicity += (it->second.size() * it->second.size());
        }
        return stdev_multiplicity / double(map.size()) - get_mean_multiplicity();
      }



      ::std::pair<subiter_type, subiter_type> equal_range_value_only(Key const & key) {
        auto iter = map.find(key);

        if (iter == map.end()) return ::std::make_pair(subiter_type(), subiter_type());

        return ::std::make_pair(iter->second.begin(), iter->second.end());
      }
      ::std::pair<const_subiter_type, const_subiter_type> equal_range_value_only(Key const & key) const {
        auto iter = map.find(key);

        if (iter == map.end()) return ::std::make_pair(const_subiter_type(), const_subiter_type());

        return ::std::make_pair(iter->second.cbegin(), iter->second.cend());

      }


      ::std::pair<iterator, iterator> equal_range(Key const & key) {
        auto iter = map.find(key);

        if (iter == map.end()) return ::std::make_pair(iterator(map.end()), iterator(map.end()));

        return ::std::make_pair(iterator(iter, map.end(), iter->second.begin(), 0),
                                iterator(iter, map.end(), iter->second.end(), iter->second.size()));
      }
      ::std::pair<const_iterator, const_iterator> equal_range(Key const & key) const {
        auto iter = map.find(key);

        if (iter == map.cend()) return ::std::make_pair(const_iterator(map.cend()), const_iterator(map.cend()));

        return ::std::make_pair(const_iterator(iter, map.cend(), iter->second.cbegin(), 0),
                                const_iterator(iter, map.cend(), iter->second.cend(), iter->second.size()));

      }

  };


  /**
   * uncompacted version of the vecmap
   */
  template <typename Key,
  typename T,
  typename Hash = ::std::hash<Key>,
  typename Equal = ::std::equal_to<Key>,
  typename Allocator = ::std::allocator<::std::pair<const Key, T> > >
  class unordered_vecmap {

    protected:
      using subcontainer_type = ::std::vector<::std::pair<const Key, T>, Allocator >;
      using superallocator_type = ::std::allocator<::std::pair<const Key, subcontainer_type > >;
      using supercontainer_type =
          ::std::unordered_map<Key, subcontainer_type, Hash, Equal, superallocator_type >;

      using subiter_type = typename subcontainer_type::iterator;
      using const_subiter_type = typename subcontainer_type::const_iterator;

      // use of vector - increased cost during construction but query will be fast.
      // group all entries with the same key with vector - list chaing but already sorted.
      // group all entries with same hash with vector - list chaining but with randomly accessible container.

      //===== create special iterator for this class - dereferences iterator to internal_map, walks through
      // vector's iterator one by one via ++ and dereference operator.
      // acts like a concatenating iterator, but source list of component iterators is from the internal map's set of vectors.

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

        public:
          template <typename KK, typename TT, typename HH, typename EE, typename AA, typename OutputIterator>
          OutputIterator
          copy(typename ::fsc::unordered_vecmap<KK, TT, HH, EE, AA>::template concat_iter<V> first,
               typename ::fsc::unordered_vecmap<KK, TT, HH, EE, AA>::template concat_iter<V> last, OutputIterator result);


        protected:
          /// the current position in the ranges list
          superiterator_type curr_iter;
          superiterator_type end_iter;

          /// the current iterator position in the range of interest.
          subiterator_type curr_pos;


        public:
          using difference_type = typename ::std::iterator_traits<subiterator_type>::difference_type;


          /// constructor for start concatenating iterator using copy semantic
          concat_iter(superiterator_type _iter, superiterator_type _end) : curr_iter(_iter), end_iter(_end) {
            if (_iter != _end) curr_pos = _iter->second.begin();
          };
          concat_iter(superiterator_type _iter, superiterator_type _end, subiterator_type _pos) :
            curr_iter(_iter), end_iter(_end), curr_pos(_pos) {};

          // note that explicit keyword cannot be on copy and move constructors else the constructors are not defined/found.

          // copy constructor, assignment operator, move constructor, assignment operator should
          // should be default since the member vars are simple.

          bool at_end() const {
            if (curr_iter == end_iter) return true;

            auto next_iter = curr_iter; ++next_iter;
            if ((next_iter == end_iter) && (curr_pos == curr_iter->second.end())) {
              return true;
            }

            return false;
          }


          /**
           * @brief increment:  move to the next position in the concatenating iterator, which may cross range boundaries.
           * @return
           */
          type& operator++()
            {
            // if at end, return
            if (curr_iter == end_iter) return *this;

            // check to see if we are at end of subcontainer.  if so, move to next dereferenceable position.
            // end of a subcontainer is treated as same position as beginning of next subcontainer.
            while (curr_pos == curr_iter->second.end()) {
              ++curr_iter;
              // again check if we are at end.  else go to beginning of next
              if (curr_iter == end_iter) return *this;
              else curr_pos = curr_iter->second.begin();
            }

            // now increment.
            ++curr_pos;

            // and check again that we are at a dereferenceable position.
            // end of a subcontainer is treated as same position as beginning of next subcontainer.
            while (curr_pos == curr_iter->second.end()) {
              ++curr_iter;
              // again check if we are at end.  else go to beginning of next
              if (curr_iter == end_iter) return *this;
              else curr_pos = curr_iter->second.begin();
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
          inline bool operator==(const type& rhs) const
            {
            if (at_end() && rhs.at_end()) return true;
            if (at_end() || rhs.at_end()) return false;

            return ((curr_iter == rhs.curr_iter) && (curr_pos == rhs.curr_pos));
            }

          /// comparison operator
          inline bool operator!=(const type& rhs) const
            {
            return !(this->operator==(rhs));
            }

          /// get pointer
          inline V* operator->() const
          {
            return &(*curr_pos);
          }

          /// dereference operator.  normal iterator: returns a rvalue
          inline V operator*() const
          {
            return *curr_pos;
          }

          //== output iterator requirement.

          /// dereference operator.  return lvalue  for non-const item
          template<typename VV = V>
          inline typename std::enable_if<
          !std::is_const<VV>::value,
          VV&
          >::type operator*()
          {
            return *curr_pos;
          }

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
            if (n < 0) throw ::std::logic_error("::fsc::unordered_vecmap::iterator does not support decrement.");
            if (n == 0) return *this;  // nothing to add.

            difference_type curr_dist;
            while ((!at_end()) && (n > 0)) {
              curr_dist = ::std::distance(curr_pos, curr_iter->second.end());
              if (curr_dist > n) {
                ::std::advance(curr_pos, n);
                n = 0;
              } else { // if (curr_dist <= n)
                n -= curr_dist;
                ++(curr_iter);
                if (curr_iter != end_iter) curr_pos = (curr_iter)->second.begin();
              }
            }
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
            difference_type n = 0;
            type it = first;

            // iterate until both are in the same subcontainer, or first is at end.
            while (!it.at_end() && (it.curr_iter != last.curr_iter)) {
              n += ::std::distance(it.curr_pos, it.curr_iter->second.end());
              ++(it.curr_iter);
              if (it.curr_iter != it.end_iter) it.curr_pos = (it.curr_iter)->second.begin();
            }

            // both are at same place (including end)
            if (it == last) return n;

            // if first is at its end, and first != last, then last is not reachable.
            if (it.at_end()) return ::std::numeric_limits<difference_type>::lowest();

            // first is not at end, and first != last, then last is in same container as first.
            n += ::std::distance(it.curr_pos, (last.curr_iter == last.end_iter) ? last.curr_iter->second.begin() : last.curr_pos);

            return n;
          }
      };


      supercontainer_type map;
      size_t s;
      size_t multiplicity;


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
      unordered_vecmap(size_type load_factor = 1,
                   size_type bucket_count = 128,
                         const Hash& hash = Hash(),
                         const Equal& equal = Equal(),
                         const Allocator& alloc = Allocator()) :
                           map(bucket_count, hash, equal, alloc),
                           s(0UL),
                           multiplicity(load_factor) {};

      template<class InputIt>
      unordered_vecmap(InputIt first, InputIt last,
                         size_type load_factor = 1,
                         size_type bucket_count = 128,
                         const Hash& hash = Hash(),
                         const Equal& equal = Equal(),
                         const Allocator& alloc = Allocator()) :
                         ::fsc::unordered_vecmap<Key, T, Hash, Equal, Allocator>(load_factor, bucket_count, hash, equal, alloc) {
          this->insert(first, last);
      };

      virtual ~unordered_vecmap() {};


      iterator begin() {
        return iterator(map.begin(), map.end());
      }
      const_iterator begin() const {
        return cbegin();
      }
      const_iterator cbegin() const {
        return const_iterator(map.cbegin(), map.cend());
      }



      iterator end() {
        return iterator(map.end(), map.end());
      }
      const_iterator end() const {
        return cend();
      }
      const_iterator cend() const {
        return const_iterator(map.cend(), map.cend());
      }



      bool empty() const {
        return map.empty();
      }

      size_type size() const {
        return s;
      }

      void clear() {
        map.clear();
      }

      /// rehash for new count number of BUCKETS.  iterators are invalidated.  side effect is multiplicity is updated.
      void rehash(size_type count) {
        // compute current average number of entries per vector.  this is side effect.
        multiplicity = (s + map.size() - 1) / map.size();

        // only rehash if new bucket count is greater than old bucket count
        if (count > map.bucket_count())
          map.rehash(count);
      }

      /// bucket count.  same as underlying buckets
      size_type bucket_count() { return map.bucket_count(); }

      /// max load factor.  this is the map's max load factor (vectors per bucket) x multiplicity = elements per bucket.  side effect is multiplicity is updated.
      float max_load_factor() {
        // compute current average number of entries per vector.  this is side effect.
        multiplicity = (s + map.size() - 1) / map.size();

        return map.max_load_factor() * (static_cast<float>(s) / static_cast<float>(map.size()));
      }


      /// reserve for new count of elements.  iterators may be invalidated.
      void reserve(size_type count) {
        // compute current average number of entries per vector.  this is side effect
        multiplicity = (s + map.size() - 1) / map.size();

        // compute number of buckets required.
        this->rehash(std::ceil(static_cast<float>(count) / this->max_load_factor()));
      }

      iterator insert(const value_type & value) {
        return emplace(::std::forward<value_type>(value));
      }
      iterator insert(value_type && value) {
        return emplace(::std::forward<value_type>(value));
      }
      iterator emplace(value_type && value) {
        auto key = value.first;
        auto iter = map.find(key);
        if (iter == map.end()) {
          iter = map.emplace(key, subcontainer_type()).first;
//          map.at(key).reserve(multiplicity);
        }
        map.at(key).emplace_back(::std::forward<value_type>(value));
        auto pos = map.at(key).end();
        ++s;
        return iterator(iter, map.end(), --pos);
      }
      iterator emplace(Key&& key, T&& value) {
        auto iter = map.find(key);
        if (iter == map.end()) {
          iter = map.emplace(key, subcontainer_type()).first;
//          map[key].reserve(multiplicity);
        }
        map.at(key).emplace_back(::std::forward<Key>(key), ::std::forward<T>(value));
        auto pos = map.at(key).end();
        ++s;
        return iterator(iter, map.end(), --pos);
      }

      // choices:  sort first, then insert in ranges, or no sort, insert one by one.  second is O(n)
      template <class InputIt>
      void insert(InputIt first, InputIt last) {
          for (; first != last; ++first) {
            if (map.find(first->first) == map.end()) {
              // emplace, and take the result iterator, get the second component (subcontainer), reserve the size.
              map.emplace(first->first, subcontainer_type());
//              map[k].reserve(multiplicity);
            }
            map.at(first->first).emplace_back(::std::forward<value_type>(*first));

            ++s;
          }
      }


      iterator erase(const key_type& key) {
        s -= count(key);
        return iterator(map.erase(key), map.end());
      }

      size_type count(Key const & key) const {
        if (map.find(key) == map.end()) return 0;
        else return map.at(key).size();
      }

      void shrink_to_fit() {
        // update multiplicity.
        this->multiplicity = (s + map.size() - 1) / map.size();

        // for each vector, shrink it.
        for (auto it = map.begin(), max = map.end(); it != max; ++it) {
          it->second.shrink_to_fit();
        }
      }


      void report() {
        printf("vecmap bucket count: %lu\n", map.bucket_count());
        printf("vecmap load factor: %f\n", map.load_factor());
        printf("vecmap unique entries: %lu\n", map.size());
        printf("vecmap total size: %lu\n", s);

      }


      size_type unique_size() const {
        return map.size();
      }
      size_type get_max_multiplicity() const {
        size_type max_multiplicity = 0;
        for (auto it = map.cbegin(), max = map.cend(); it != max; ++it) {
          max_multiplicity = ::std::max(max_multiplicity, it->second.size());
        }
        return max_multiplicity;
      }
      size_type get_min_multiplicity() const {
        size_type min_multiplicity = ::std::numeric_limits<size_type>::max();
        for (auto it = map.cbegin(), max = map.cend(); it != max; ++it) {
          min_multiplicity = ::std::min(min_multiplicity, it->second.size());
        }
        return min_multiplicity;
      }
      double get_mean_multiplicity() const {
        return static_cast<double>(s) / double(map.size());
      }
      double get_stdev_multiplicity() const {
        double stdev_multiplicity = 0;
        for (auto it = map.cbegin(), max = map.cend(); it != max; ++it) {
          stdev_multiplicity += (it->second.size() * it->second.size());
        }
        return stdev_multiplicity / double(map.size()) - get_mean_multiplicity();
      }



      ::std::pair<subiter_type, subiter_type> equal_range(Key const & key) {
        auto iter = map.find(key);

        if (iter == map.end()) return ::std::make_pair(subiter_type(), subiter_type());

        return ::std::make_pair(iter->second.begin(), iter->second.end());
      }
      ::std::pair<const_subiter_type, const_subiter_type> equal_range(Key const & key) const {
        auto iter = map.find(key);

        if (iter == map.end()) return ::std::make_pair(const_subiter_type(), const_subiter_type());

        return ::std::make_pair(iter->second.cbegin(), iter->second.cend());

      }


      // NO bucket interfaces

  };



} // end namespace fsc.


namespace std {

  /// specialization for vecmap's iterator.  atomic.
  template <typename Key, typename T, typename Hash, typename Equal, typename Allocator, typename OutputIterator>
  OutputIterator
  copy(typename ::fsc::unordered_compact_vecmap<Key, T, Hash, Equal, Allocator>::iterator first,
       typename ::fsc::unordered_compact_vecmap<Key, T, Hash, Equal, Allocator>::iterator last, OutputIterator result) {

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
      out_iter = ::std::transform(cpos, cit->second.end(), out_iter, [&cit](T const & x) {return ::std::make_pair(cit->first, x); });
      ++cit;
      cpos = cit->second.begin();
    }

    // now we're in same container.  take care the rest.  cpos was updated if cit was moved.
    out_iter = ::std::transform(cpos, last.curr_pos, out_iter, [&cit](T const & x) {return ::std::make_pair(cit->first, x); });

    //    // iterate until both are in the same subcontainer, or first is at end.
    //    while (!first.at_end() && (first.curr_iter != last.curr_iter)) {
    //      out_iter = ::std::copy(first.curr_pos, (first.curr_iter)->second.end(), out_iter);
    //      ++(first.curr_iter);
    //      if (first.curr_iter != first.end_iter) first.curr_pos = (first.curr_iter)->second.begin();
    //    }
    //
    //    // both are at same place (including end)
    //    if (first == last) return out_iter;
    //
    //    // if first is at its end, and first != last, then last is not reachable.
    //    if (first.at_end()) return result;
    //
    //    // first is not at end, and first != last, then last is in same container as first.
    //    out_iter = ::std::copy(first.curr_pos, (last.curr_iter == last.end_iter) ? last.curr_iter->second.begin() : last.curr_pos, out_iter);

    return out_iter;
  }

  template <typename Key, typename T, typename Hash, typename Equal, typename Allocator, typename OutputIterator>
  OutputIterator
  copy(typename ::fsc::unordered_compact_vecmap<Key, T, Hash, Equal, Allocator>::const_iterator first,
       typename ::fsc::unordered_compact_vecmap<Key, T, Hash, Equal, Allocator>::const_iterator last, OutputIterator result) {
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
      out_iter = ::std::transform(cpos, cit->second.end(), out_iter, [&cit](T const & x) {return ::std::make_pair(cit->first, x); });
      ++cit;
      cpos = cit->second.begin();
    }

    // now we're in same container.  take care the rest.  cpos was updated if cit was moved.
    out_iter = ::std::transform(cpos, last.curr_pos, out_iter, [&cit](T const & x) {return ::std::make_pair(cit->first, x); });

    //    // iterate until both are in the same subcontainer, or first is at end.
    //    while (!first.at_end() && (first.curr_iter != last.curr_iter)) {
    //      out_iter = ::std::copy(first.curr_pos, (first.curr_iter)->second.end(), out_iter);
    //      ++(first.curr_iter);
    //      if (first.curr_iter != first.end_iter) first.curr_pos = (first.curr_iter)->second.begin();
    //    }
    //
    //    // both are at same place (including end)
    //    if (first == last) return out_iter;
    //
    //    // if first is at its end, and first != last, then last is not reachable.
    //    if (first.at_end()) return result;
    //
    //    // first is not at end, and first != last, then last is in same container as first.
    //    out_iter = ::std::copy(first.curr_pos, (last.curr_iter == last.end_iter) ? last.curr_iter->second.begin() : last.curr_pos, out_iter);

    return out_iter;
  }



  template <typename Key, typename T, typename Hash, typename Equal, typename Allocator, typename OutputIterator>
  OutputIterator
  copy(typename ::fsc::unordered_vecmap<Key, T, Hash, Equal, Allocator>::iterator first,
       typename ::fsc::unordered_vecmap<Key, T, Hash, Equal, Allocator>::iterator last, OutputIterator result) {
    auto out_iter = result;

    // iterate until both are in the same subcontainer, or first is at end.
    while (!first.at_end() && (first.curr_iter != last.curr_iter)) {
      out_iter = ::std::copy(first.curr_pos, (first.curr_iter)->second.end(), out_iter);
      ++(first.curr_iter);
      if (first.curr_iter != first.end_iter) first.curr_pos = (first.curr_iter)->second.begin();
    }

    // both are at same place (including end)
    if (first == last) return out_iter;

    // if first is at its end, and first != last, then last is not reachable.
    if (first.at_end()) return result;

    // first is not at end, and first != last, then last is in same container as first.
    out_iter = ::std::copy(first.curr_pos, (last.curr_iter == last.end_iter) ? last.curr_iter->second.begin() : last.curr_pos, out_iter);

    return out_iter;
  }

  template <typename Key, typename T, typename Hash, typename Equal, typename Allocator, typename OutputIterator>
  OutputIterator
  copy(typename ::fsc::unordered_vecmap<Key, T, Hash, Equal, Allocator>::const_iterator first,
       typename ::fsc::unordered_vecmap<Key, T, Hash, Equal, Allocator>::const_iterator last, OutputIterator result) {
    auto out_iter = result;

    // iterate until both are in the same subcontainer, or first is at end.
    while (!first.at_end() && (first.curr_iter != last.curr_iter)) {
      out_iter = ::std::copy(first.curr_pos, (first.curr_iter)->second.end(), out_iter);
      ++(first.curr_iter);
      if (first.curr_iter != first.end_iter) first.curr_pos = (first.curr_iter)->second.begin();
    }

    // both are at same place (including end)
    if (first == last) return out_iter;

    // if first is at its end, and first != last, then last is not reachable.
    if (first.at_end()) return result;

    // first is not at end, and first != last, then last is in same container as first.
    out_iter = ::std::copy(first.curr_pos, (last.curr_iter == last.end_iter) ? last.curr_iter->second.begin() : last.curr_pos, out_iter);

    return out_iter;
  }




}



#endif /* SRC_WIP_UNORDERED_VECMAP_HPP_ */
