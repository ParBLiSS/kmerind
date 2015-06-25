/**
 * @file    unordered_multimap.hpp
 * @ingroup
 * @author  tpan
 * @brief
 * @details
 *
 * Copyright (c) 2015 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef SRC_WIP_UNORDERED_MULTIMAP_HPP_
#define SRC_WIP_UNORDERED_MULTIMAP_HPP_

#include <unordered_map>
#include <vector>
#include <functional>  // hash, equal_to, etc
#include <tuple>   // pair


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
   */
  template <typename Key,
  typename T,
  typename Hash = ::std::hash<Key>,
  typename Equal = ::std::equal_to<Key>,
  typename Allocator = ::std::allocator<::std::pair<const Key, T> > >
  class unordered_multimap {

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
          copy(typename ::fsc::unordered_multimap<KK, TT, HH, EE, AA>::template concat_iter<V> first,
               typename ::fsc::unordered_multimap<KK, TT, HH, EE, AA>::template concat_iter<V> last, OutputIterator result);


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
            return ++output;
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
            if (n < 0) throw ::std::logic_error("::fsc::unordered_multimap::iterator does not support decrement.");
            if (n == 0) return *this;  // nothing to add.

            size_t curr_dist;
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
      using reference 			      = value_type&;
      using const_reference	      = const value_type&;
      using pointer				        = typename std::allocator_traits<Allocator>::pointer;
      using const_pointer		      = typename std::allocator_traits<Allocator>::const_pointer;
      using iterator              = concat_iter<value_type>;
      using const_iterator        = concat_iter<const value_type>;
      using size_type             = typename subcontainer_type::size_type;
      using difference_type       = typename subcontainer_type::difference_type;


      //  if multiplicity of dataset is kind of known, can initialize data to that to avoid growing vector on average.
      unordered_multimap(size_type load_factor = 1, size_type bucket_count = 128,
                         const Hash& hash = Hash(),
                         const Equal& equal = Equal(),
                         const Allocator& alloc = Allocator()) :
                           map(bucket_count, hash, equal, alloc),
                           s(0UL),
                           multiplicity(load_factor) {};

      template<class InputIt>
      unordered_multimap(InputIt first, InputIt last,
                         size_type load_factor = 1,
                         size_type bucket_count = 128,
                         const Hash& hash = Hash(),
                         const Equal& equal = Equal(),
                         const Allocator& alloc = Allocator()) :
                         ::fsc::unordered_multimap<Key, T, Hash, Equal, Allocator>(load_factor, bucket_count, hash, equal, alloc) {
          this->insert(first, last);
      };

      ~unordered_multimap() {};


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

      /// rehash for new count.  iterators are invalidated.
      void rehash(size_type count) {
        map.rehash((count + multiplicity - 1) / multiplicity);
      }

      /// reserve for new count.  iterators may be invalidated.
      void reserve(size_type count) {
        map.reserve((count + multiplicity - 1) / multiplicity);
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
          map[key].reserve(multiplicity);
        }
        auto pos = map[key].end();
        map[key].emplace_back(::std::forward<value_type>(value));
        ++s;
        return iterator(iter, map.end(), pos);
      }
      iterator emplace(Key&& key, T&& value) {
        auto iter = map.find(key);
        if (iter == map.end()) {
          iter = map.emplace(key, subcontainer_type()).first;
          map[key].reserve(multiplicity);
        }
        auto pos = map[key].end();
        map[key].emplace_back(::std::forward<Key>(key), ::std::forward<T>(value));
        ++s;
        return iterator(iter, map.end(), pos);
      }

      // choices:  sort first, then insert in ranges, or no sort, insert one by one.  second is O(n)
      template <class InputIt>
      void insert(InputIt first, InputIt last) {
          for (; first != last; ++first, ++s) {
            Key k = first->first;
            if (map.find(k) == map.end()) {
              map.emplace(k, subcontainer_type());
              map[k].reserve(multiplicity);
            }
            map[k].emplace_back(*first);
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

      size_type unique_size() const {
        return map.size();
      }


//      ::std::pair<iterator, iterator> equal_range(Key const & key) {
//        auto iter = map.find(key);
//
//        if (iter == map.end()) return ::std::make_pair(iterator(iter, iter),
//                                                       iterator(iter, iter));
//
//        auto iter2 = iter; ++iter2;
//        return ::std::make_pair(iterator(iter, iter2, iter->second.begin()),
//                                iterator(iter, iter2, iter->second.end()));
//      }
//      ::std::pair<const_iterator, const_iterator> equal_range(Key const & key) const {
//        auto iter = map.find(key);
//
//        if (iter == map.end()) return ::std::make_pair(const_iterator(iter, iter),
//                                                       const_iterator(iter, iter));
//
//        auto iter2 = iter; ++iter2;
//        return ::std::make_pair(const_iterator(iter, iter2, iter->second.cbegin()),
//                                const_iterator(iter, iter2, iter->second.cend()));
//
//      }

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



  template <typename Key, typename T, typename Hash, typename Equal, typename Allocator, typename OutputIterator>
  OutputIterator
  copy(typename ::fsc::unordered_multimap<Key, T, Hash, Equal, Allocator>::iterator first,
       typename ::fsc::unordered_multimap<Key, T, Hash, Equal, Allocator>::iterator last, OutputIterator result) {
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
  copy(typename ::fsc::unordered_multimap<Key, T, Hash, Equal, Allocator>::const_iterator first,
       typename ::fsc::unordered_multimap<Key, T, Hash, Equal, Allocator>::const_iterator last, OutputIterator result) {
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



#endif /* SRC_WIP_UNORDERED_MAP_HPP_ */
