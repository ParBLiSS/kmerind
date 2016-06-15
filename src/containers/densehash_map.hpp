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
 * @file    densehash_map.hpp
 * @ingroup
 * @author  tpan
 * @brief
 * @details
 *
 */
#ifndef SRC_WIP_DENSEHASH_VECMAP_HPP_
#define SRC_WIP_DENSEHASH_VECMAP_HPP_

#include <sparsehash/dense_hash_map>
#include <vector>
#include <functional>  // hash, equal_to, etc
#include <tuple>   // pair
#include <scoped_allocator>
#include <algorithm>
#include <cmath>   // ceil
#include <memory>  // allocator
#include <iostream>

#include "iterators/transform_iterator.hpp"
#include "iterators/filter_iterator.hpp"
#include "iterators/counting_iterator.hpp"
#include "iterators/zip_iterator.hpp"
#include "iterators/unzip_iterator.hpp"

#include "containers/fsc_container_utils.hpp"

#include "utils/logging.h"

namespace fsc {  // fast standard container


namespace sparsehash {

	// ==========
	// sparsehash specific functors, not specific to key type.

  template <typename Key, template <typename> class Comparator, template <typename> class Transform>
  struct threshold {

      Comparator<Key> comp;
      Transform<Key> trans;

      mutable Key upper_bound;

      bool operator()(Key const & x) const {
        return comp(trans(x), upper_bound); // lower is positive
      }
      template <typename V>
      bool operator()(::std::pair<Key, V> const & x) const {
        return this->operator()(x.first); // lower is positive.
      }
    };


  /// special comparison operator with knowledge of the deleted and empty keys.
  template <typename Key, template <typename> class Comparator, template <typename> class Transform>
  struct compare {
      Comparator<Key> comp;
      Transform<Key> trans;

      mutable Key empty;
      mutable Key deleted;

      //====  since dense hash table makes copies of equal operators left and right
      // we need these constructors and assignment operators.
      // also, note that there is no default constructor since the keys need to be set.

      compare() = delete;

      compare(Key const & em, Key const & del) : empty(em), deleted(del) {}

      compare(compare  const & other) : empty(other.empty), deleted(other.deleted) {}
      compare(compare && other) : empty(other.empty), deleted(other.deleted) {}

      compare & operator=(compare  const & other) {
    		  empty = other.empty;
    		  deleted = other.deleted;

    		  return *this;
      }
      compare & operator=(compare && other) {
    		  empty = other.empty;
    		  deleted = other.deleted;

    		  return *this;
      }


      inline bool is_key(Key const & x) const {
        return (x == empty) || (x == deleted);  // not using comp because it may not be equal comparator.
      }

      /// comparison.  if key, don't transform them
      inline bool operator()(Key const & x, Key const & y) const {
        return comp((is_key(x) ? x : trans(x)), (is_key(y) ? y : trans(y)));
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

  /// special comparison operator for sparsehash.  specialized for equal_to operator..
  template <typename Key, template <typename> class Transform>
  struct compare<Key, ::std::equal_to, Transform> {
      Transform<Key> trans;

      mutable Key empty;
      mutable Key deleted;

      //====  since dense hash table makes copies of equal operators left and right
      // we need these constructors and assignment operators.
      // also, note that there is no default constructor since the keys need to be set.

      compare() = delete;

      compare(Key const & em, Key const & del) : empty(em), deleted(del) {}

      compare(compare  const & other) : empty(other.empty), deleted(other.deleted) {}
      compare(compare && other) : empty(other.empty), deleted(other.deleted) {}

      compare & operator=(compare  const & other) {
    		  empty = other.empty;
    		  deleted = other.deleted;

    		  return *this;
      }
      compare & operator=(compare && other) {
    		  empty = other.empty;
    		  deleted = other.deleted;

    		  return *this;
      }


      inline bool is_key(Key const & x) const {
        return (x == empty) || (x == deleted);  // not using comp because it may not be equal comparator.
      }

      /// comparison operator.  if key, don't transform them.  also shortcuts some cases
      inline bool operator()(Key const & x, Key const & y) const {
        bool x_is_key = is_key(x);
        bool y_is_key = is_key(y);

//        std::cout << " x " << x << std::endl;
//        std::cout << " y " << y << std::endl;
//        std::cout << " empty " << empty << std::endl;
//        std::cout << " deleted " << deleted << std::endl;
//        std::cout << " x is empty ? " <<  ((x == empty) ? "y" : "n")  << " x is deleted ? " <<  (comp(x, deleted) ? "y" : "n") << std::endl;
//        std::cout << " y is empty ? " <<  (comp(y, empty) ? "y" : "n")  << " y is deleted ? " <<  (comp(y, deleted) ? "y" : "n") << std::endl;
//        std::cout << " x y same ? " <<  (comp(x, y) ? "y" : "n")  << " trans x y same ? " <<  (comp(trans(x), trans(y)) ? "y" : "n") << std::endl;

        return (x_is_key != y_is_key) ? false : (x_is_key) ? (x == y) : (trans(x) == trans(y));
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


  /// specialized sparsehash comparator for when there is no transform.
  template <typename Key, template <typename> class Comparator>
  struct compare<Key, Comparator, fsc::identity> {
      Comparator<Key> comp;

      // keys are not transformed.  so no need to treat them specially, and no need to store them.

      //====  since dense hash table makes copies of equal operators left and right
      // we need these constructors and assignment operators.
      // also, note that there is no default constructor since the keys need to be set.

      compare() = delete;

      compare(Key const & em, Key const & del) {}

      compare(compare  const & other)  {}
      compare(compare && other)  {}

      compare & operator=(compare  const & other) {
    		  return *this;
      }
      compare & operator=(compare && other) {
    		  return *this;
      }

      inline bool operator()(Key const & x, Key const & y) const {
        return comp(x, y);
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

  /// special comparison operator for sparsehash.  specialized for equal_to operator AND identity transform, here so that template is not ambiguous for previous 2 definitions.
  template <typename Key>
  struct compare<Key, ::std::equal_to, ::fsc::identity> {

      //====  since dense hash table makes copies of equal operators left and right
      // we need these constructors and assignment operators.
      // also, note that there is no default constructor since the keys need to be set.

      compare() = delete;

      compare(Key const & em, Key const & del) {}

      compare(compare  const & other) {}
      compare(compare && other)  {}

      compare & operator=(compare  const & other) {
    		  return *this;
      }
      compare & operator=(compare && other) {
    		  return *this;
      }


      /// comparison operator.  if key, don't transform them.  also shortcuts some cases
      inline bool operator()(Key const & x, Key const & y) const {
        return (x == y);
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



}


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
   *            hashmap:  build 1.6s, find 1.5s, count 0.56s
   *            sorted vector:  build = 2 s, find 0.27s, count 0.01s
   *          with google dense_hash_map:  can't use it - google dense_hash_map requires uniqueness, and also requiers a NULL element.
   *
   *          tested using std::multimap and map.  VERY slow on insertion, especially for real data.  about 33% faster for find and count compared to unordered map for synthetic, high repeat datasets.  still about 4 times slower than sort based.
   *          for that dataset.  build = 3.6s, find 0.98s, count 0.73s
   *
   *          WE HAVE TO BUILD OUR OWN VERSION.
   *            prefer using std::vector as internal - self growing,
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
   *      vector holds pair<Key, T> because then we can directly copy (fast)
   *      a version with pair<T> would be more space efficient but slower when large amount of data is present.
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
   *
   *
   *
   *  this file contains the GOOGLE DENSE HASH MAP version of map and multimap that are compatible with kmer indexing.
   */
// key values span entire key space.
template <typename Key,
typename T,
bool split = true,   // for full kmers, set this to true.
template<typename> class Transform = ::fsc::identity,
typename Hash =  ::fsc::TransformedHash<Key, ::std::hash, Transform>,
typename Equal = ::fsc::sparsehash::compare<Key, ::std::equal_to, Transform>,
typename Allocator = ::std::allocator<::std::pair<const Key, T> > >
class densehash_map {

protected:
    using supercontainer_type =
        ::google::dense_hash_map<Key, T,
                       Hash, Equal, Allocator >;

    using Splitter = ::fsc::sparsehash::threshold<Key, std::less, Transform>;

    Splitter splitter;

    supercontainer_type lower_map;
    supercontainer_type upper_map;


    template <typename InputIt>
    InputIt partition_input(InputIt first, InputIt last) {
    	return ::std::stable_partition(first, last, splitter);
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
//    using iterator              = concat_iter<::std::pair<const Key, T> >;
//    using const_iterator        = concat_iter<const ::std::pair<const Key, T> >;
    using iterator              = typename supercontainer_type::iterator;
    using const_iterator              = typename supercontainer_type::const_iterator;
    using size_type             = size_t;
    using difference_type       = ptrdiff_t;

//    densehash_map(size_type bucket_count = 128,
//                       const Hash& hash = Hash(),
//                       const Equal& equal = Equal(),
//                       const Allocator& alloc = Allocator()) :
//             lower_map(bucket_count / 2, hash, equal, alloc),
//						 upper_map(bucket_count / 2, hash, Equal(), alloc)
//    {
//    };

    //  if multiplicity of dataset is kind of known, can initialize data to that to avoid growing vector on average.
    densehash_map(Key const & empty_key, Key const & deleted_key,
    		 	 	   Key const & upper_empty_key, Key const & upper_deleted_key,
    		 	 	   Key const & split_key,
             size_type bucket_count = 128,
                       const Hash& hash = Hash(),
                       const Equal& equal = Equal(),
                       const Allocator& alloc = Allocator()) :
                    	   lower_map(bucket_count / 2, hash, Equal(empty_key, deleted_key), alloc),
						   upper_map(bucket_count / 2, hash, Equal(upper_empty_key, upper_deleted_key), alloc)
    {

        reserve_keys(empty_key, deleted_key);
    	reserve_upper_keys(upper_empty_key, upper_deleted_key);
    	reserve_split(split_key);
    };

    template<class InputIt>
    densehash_map(InputIt first, InputIt last,
    		Key const & empty_key, Key const & deleted_key,
    		    		 	 	   Key const & upper_empty_key, Key const & upper_deleted_key,
    		    		 	 	 Key const & split_key,
                       size_type bucket_count = 128,
                       const Hash& hash = Hash(),
                       const Equal& equal = Equal(),
                       const Allocator& alloc = Allocator()) :
					   densehash_map(empty_key, deleted_key, upper_empty_key, upper_deleted_key,
							   split_key, std::distance(first, last), hash, equal, alloc) {

    	this->insert(first, last);
    };

    virtual ~densehash_map() {};

protected:
    void reserve_keys(Key const & empty_key, Key const & deleted_key) {
      lower_map.set_empty_key(empty_key);
      lower_map.set_deleted_key(deleted_key);
    }

    void reserve_upper_keys(Key const & empty_key, Key const & deleted_key) {
      upper_map.set_empty_key(empty_key);
      upper_map.set_deleted_key(deleted_key);
    }

    void reserve_split(Key const & split_key) {
      splitter.upper_bound = split_key;
    }
public:

//
//    iterator begin() {
//    	return concat_iter<value_type>(lower_map, upper_map, lower_map.begin());
//    }
//    const_iterator begin() const {
//      return cbegin();
//    }
//    const_iterator cbegin() const {
//    	return concat_iter<const value_type>(lower_map, upper_map, lower_map.begin());
//    }
//
//
//
//    iterator end() {
//    	return concat_iter<value_type>(lower_map, upper_map);
//    }
//    const_iterator end() const {
//      return cend();
//    }
//    const_iterator cend() const {
//    	return concat_iter<const value_type>(lower_map, upper_map);
//    }
    std::vector<Key> keys() const  {
      std::vector<Key> ks;

      keys(ks);

      return ks;
    }
    void keys(std::vector<Key> & ks) const  {
      ks.clear();
      ks.reserve(size());

      for (auto it = lower_map.begin(); it != lower_map.end(); ++it) {
        ks.emplace_back(it->first);
      }
      for (auto it = upper_map.begin(); it != upper_map.end(); ++it) {
        ks.emplace_back(it->first);
      }
    }

    std::vector<std::pair<Key, T>> to_vector() const  {
      std::vector<std::pair<Key, T>> vs;

      to_vector(vs);

      return vs;
    }
    void to_vector(  std::vector<std::pair<Key, T>> & vs) const  {
      vs.clear();
      vs.reserve(size());

      for (auto it = lower_map.begin(); it != lower_map.end(); ++it) {
        vs.emplace_back(*it);
      }
      for (auto it = upper_map.begin(); it != upper_map.end(); ++it) {
        vs.emplace_back(*it);
      }

    }



    bool empty() const {
      return lower_map.empty() && upper_map.empty();
    }

    size_type size() const {
      return lower_map.size() + upper_map.size();
    }

    size_type unique_size() const {
      return lower_map.size() + upper_map.size();
    }

    void clear() {
      lower_map.clear_no_resize();
      upper_map.clear_no_resize();
    }

    void resize(size_t const n) {
      lower_map.resize(static_cast<float>(n) / 2.0 / lower_map.max_load_factor() );
      upper_map.resize(static_cast<float>(n) / 2.0 / upper_map.max_load_factor() );

    }

    /// rehash for new count number of BUCKETS.  iterators are invalidated.
    void rehash(size_type count) {
      this->resize(count);
    }

    /// bucket count.  same as underlying buckets
    size_type bucket_count() { return lower_map.bucket_count() + upper_map.bucket_count(); }

    /// max load factor.  this is the map's max load factor (vectors per bucket) x multiplicity = elements per bucket.  side effect is multiplicity is updated.
    float load_factor() {
      return  static_cast<float>(size()) / static_cast<float>(bucket_count());
    }



    // choices:  sort first, then insert in ranges, or no sort, insert one by one.  second is O(n) but pays the random access and mem realloc cost
    template <class InputIt>
    void insert(InputIt first, InputIt last) {

        InputIt middle = partition_input(first, last);

        lower_map.resize(static_cast<float>(lower_map.size() + ::std::distance(first, middle)) / lower_map.max_load_factor() ) ;
        lower_map.insert(first, middle);

        upper_map.resize(static_cast<float>(upper_map.size() + ::std::distance(middle, last)) / upper_map.max_load_factor()) ;
        upper_map.insert(middle, last);
    }

    /// inserting sorted range
    void insert(::std::vector<::std::pair<Key, T> > & input) {
    	insert(input.begin(), input.end());
    }

    std::pair<typename supercontainer_type::iterator, bool> insert(::std::pair<Key, T> const & x) {
      if (splitter(x.first)) {
        return lower_map.insert(x);
      }
      else {
        return upper_map.insert(x);
      }
    }

    template <typename InputIt, typename Pred>
    size_t erase(InputIt first, InputIt last, Pred const & pred) {
        static_assert(::std::is_convertible<Key, typename ::std::iterator_traits<InputIt>::value_type>::value,
                      "InputIt value type for erase cannot be converted to key type");


      if (first == last) return 0;

      size_t count = 0;

      InputIt middle = partition_input(first, last);

      // mark for erasure
      for (; first != middle; ++first) {
        auto iter = lower_map.find(*(first));
        if (iter == lower_map.end()) continue;

        if (pred(*iter)) {
        	lower_map.erase(iter);
        	++count;
        }
      }

      for (; first != last; ++first) {
        auto iter = upper_map.find(*(first));
        if (iter == upper_map.end()) continue;

        if (pred(*iter)) {
        	upper_map.erase(iter);
        	++count;
        }
      }

      return count;
    }

    template <typename InputIt>
    size_t erase(InputIt first, InputIt last) {

        static_assert(::std::is_convertible<Key, typename ::std::iterator_traits<InputIt>::value_type>::value,
                      "InputIt value type for erase cannot be converted to key type");

        if (first == last) return 0;

        size_t count = 0;


        InputIt middle = partition_input(first, last);


        // mark for erasure
        for (; first != middle; ++first) {
          auto iter = lower_map.find(*(first));
          if (iter == lower_map.end()) continue;

          lower_map.erase(iter);
          ++count;
        }

        for (; first != last; ++first) {
          auto iter = upper_map.find(*(first));
          if (iter == upper_map.end()) continue;

          	upper_map.erase(iter);
          	++count;
        }


        return count;
    }

    template <typename Pred>
    size_t erase(Pred const & pred) {
    	size_t before = size();

    	for (auto it = lower_map.begin(); it != lower_map.end(); ++it) {
    		if (pred(*it))
    				lower_map.erase(it);
    	}
    	for (auto it = upper_map.begin(); it != upper_map.end(); ++it) {
    		if (pred(*it))
    				upper_map.erase(it);
    	}

    	return before - size();
    }

    size_type count(Key const & key) const {
    	if (splitter(key))
    		return lower_map.count(key);
    	else
    		return upper_map.count(key);
    }


    ::std::pair<iterator, iterator> equal_range(Key const & key) {
    	if (splitter(key)) {
    		return lower_map.equal_range(key);
    	}
    	else {
    		return upper_map.equal_range(key);
    	}

    }
    ::std::pair<const_iterator, const_iterator> equal_range(Key const & key) const {
    	if (splitter(key)) {
    		return lower_map.equal_range(key);
    	}
    	else {
    		return upper_map.equal_range(key);
    	}
    }
    // NO bucket interfaces

};


// Key values does not span entire key space.
template <typename Key,
typename T,
template <typename > class Transform,
typename Hash,
typename Equal,
typename Allocator >
class densehash_map<Key, T, false, Transform, Hash, Equal, Allocator> {

  protected:

    using supercontainer_type =
        ::google::dense_hash_map<Key, T,
                       Hash, Equal, Allocator >;

    supercontainer_type map;



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
    using iterator              = typename supercontainer_type::iterator;
    using const_iterator        = typename supercontainer_type::const_iterator;
    using size_type             = size_t;
    using difference_type       = ptrdiff_t;


    //  if multiplicity of dataset is kind of known, can initialize data to that to avoid growing vector on average.
//    densehash_map(size_type bucket_count = 128,
//                       const Hash& hash = Hash(),
//                       const Equal& equal = Equal(),
//                       const Allocator& alloc = Allocator()) :
//                         map(bucket_count, hash, equal, alloc) {
//    };

    densehash_map(Key empty_key, Key deleted_key,
                   size_type bucket_count = 128,
                       const Hash& hash = Hash(),
                       const Equal& equal = Equal(),
                       const Allocator& alloc = Allocator()) :
                    	   map(bucket_count, hash, Equal(empty_key, deleted_key), alloc) {
    	reserve_keys(empty_key, deleted_key);
    };

    template<class InputIt>
    densehash_map(InputIt first, InputIt last,
                      Key empty_key, Key deleted_key,
                       size_type bucket_count = 128,
                       const Hash& hash = Hash(),
                       const Equal& equal = Equal(),
                       const Allocator& alloc = Allocator()) :
                       densehash_map(empty_key, deleted_key, std::distance(first, last), hash, equal, alloc) {

        this->insert(first, last);
    };

    virtual ~densehash_map() {};

protected:
    void reserve_keys(Key const & empty_key, Key const & deleted_key) {
    	map.set_empty_key(empty_key);
        map.set_deleted_key(deleted_key);
    }
    void reserve_upper_keys(Key const & empty_key, Key const & deleted_key) {
    	//printf("reserve_upper_keys no-op\n");
    }
    void reserve_split(Key const & split_key) {
    }

  public:

    iterator begin() {
      return map.begin();
    }
    const_iterator begin() const {
      return cbegin();
    }
    const_iterator cbegin() const {
      return map.begin();
    }



    iterator end() {
      return map.begin();
    }
    const_iterator end() const {
      return cend();
    }
    const_iterator cend() const {
      return map.begin();
    }


    std::vector<Key> keys() const {
      std::vector<Key> ks;

      keys(ks);

      return ks;
    }
    void keys(std::vector<Key> & ks) const {
      ks.clear();
      ks.reserve(size());

      for (auto it = map.begin(); it != map.end(); ++it) {
        ks.emplace_back(it->first);
      }
    }

    std::vector<std::pair<Key, T> > to_vector() const {
      std::vector<std::pair<Key, T>> vs;

      to_vector(vs);

      return vs;
    }
    void to_vector(  std::vector<std::pair<Key, T> > & vs) const {
      vs.clear();
      vs.reserve(size());

     for (auto it = map.begin(); it != map.end(); ++it) {
        vs.emplace_back(*it);
      }
    }


    bool empty() const {
      return map.empty();
    }

    size_type size() const {
      return map.size();
    }
    size_type unique_size() const {
      return map.size();
    }

    void clear() {
      map.clear_no_resize();
    }

    void resize(size_t const n) {
      map.resize(static_cast<float>(n) / map.max_load_factor());
    }

    /// rehash for new count number of BUCKETS.  iterators are invalidated.
    void rehash(size_type count) {
      this->resize(count);
    }

    /// bucket count.  same as underlying buckets
    size_type bucket_count() { return map.bucket_count(); }

    /// max load factor.  this is the map's max load factor (vectors per bucket) x multiplicity = elements per bucket.  side effect is multiplicity is updated.
    float load_factor() {
      return  static_cast<float>(map.size()) / static_cast<float>(map.bucket_count());
    }



    // choices:  sort first, then insert in ranges, or no sort, insert one by one.  second is O(n) but pays the random access and mem realloc cost
    template <class InputIt>
    void insert(InputIt first, InputIt last) {
      this->resize(map.size() + std::distance(first, last));

      map.insert(first, last);
    }

    /// inserting sorted range
    void insert(::std::vector<::std::pair<Key, T> > & input) {
      insert(input.begin(), input.end());
    }

    std::pair<typename supercontainer_type::iterator, bool> insert(::std::pair<Key, T> const & x) {
      return map.insert(x);
    }


    template <typename InputIt, typename Pred>
    size_t erase(InputIt first, InputIt last, Pred const & pred) {
      static_assert(::std::is_convertible<Key, typename ::std::iterator_traits<InputIt>::value_type>::value,
                    "InputIt value type for erase cannot be converted to key type");

      if (first == last) return 0;

      size_t count = 0;

      // mark for erasure
      for (; first != last; ++first) {
        auto iter = map.find(*(first));
        if (iter == map.end()) continue;

        if (pred(*iter)) {
          map.erase(iter);
          ++count;
        }
      }
      return count;
    }

    template <typename InputIt>
    size_t erase(InputIt first, InputIt last) {
        static_assert(::std::is_convertible<Key, typename ::std::iterator_traits<InputIt>::value_type>::value,
                      "InputIt value type for erase cannot be converted to key type");

        if (first == last) return 0;

        size_t count = 0;

        // mark for erasure
        for (; first != last; ++first) {
          auto iter = map.find(*(first));
          if (iter == map.end()) continue;

            map.erase(iter);
            ++count;
        }
        return count;
    }

    template <typename Pred>
    size_t erase(Pred const & pred) {
      size_t before = map.size();

      for (auto it = map.begin(); it != map.end(); ++it) {
        if (pred(*it))
            map.erase(it);
      }

      return before - map.size();
    }

    size_type count(Key const & key) const {
      return map.count(key);
    }


    ::std::pair<iterator, iterator> equal_range(Key const & key) {
      return map.equal_range(key);
    }
    ::std::pair<const_iterator, const_iterator> equal_range(Key const & key) const {
      return map.equal_range(key);
    }
    // NO bucket interfaces

};







/**
 * @brief multimap implemented using densehash.
 * @details:  earlier attempt used a sorted vector.  while still faster than unordered map, the sorting step does not scale well.
 *            a second attempt tried to use key+counter as key, hash only on key, compare on full key+counter, and query via bucket interface and compare only key portion.
 *              this did not work as google dense hash map bucket points to array element, and does not represent a collection of entries with same hash value like
 *              std::unordered_multimap.
 *
 *            current implementation uses 1 vector for singletons, and 1 vector of vectors for multiple entries.  memory utilization is probably suboptimal (up to 2x)
 *            map stores index in these vectors, with sign bit signifying that it is a multiple entry.
 *
 *            this trades off sorting with more vector allocations and copying, but appears to be faster at least for the case when there are large number of singleton entries.
 *
 */
template <typename Key,
typename T,
bool split = true,   // for full kmers, set this to true.
template<typename> class Transform = ::fsc::identity,
typename Hash = ::fsc::TransformedHash<Key, ::std::hash, Transform>,
typename Equal = ::fsc::sparsehash::compare<Key, std::equal_to, Transform>,
typename Allocator = ::std::allocator<::std::pair<const Key, T> > >
class densehash_multimap {

  protected:

    // data container
    using tuple_allocator_type = typename Allocator::template rebind< ::std::pair<Key, T> >::other;
    using subcontainer_type = ::std::vector<::std::pair<Key, T>, tuple_allocator_type >;
    using subiter_type = typename subcontainer_type::iterator;
    using const_subiter_type = typename subcontainer_type::const_iterator;

    // index in the vector - this is so we don't have to worry about pointer or iterators being invalidated when vector resizes
    // non-negative values indicate singletons.  negative values indicate multiple entries.
    // multiple entries is stored in a vector of vector.  the index is internal_val_type with sign bit removed.
    using internal_val_type = int64_t;

    using super_allocator_type = typename Allocator::template rebind<std::pair< const Key, internal_val_type> >::other;
    using supercontainer_type =
        ::google::dense_hash_map<Key, internal_val_type,
                       Hash, Equal, super_allocator_type >;

    using Splitter = ::fsc::sparsehash::threshold<Key, std::less, Transform>;

    Splitter splitter;

    supercontainer_type lower_map;
    supercontainer_type upper_map;

    subcontainer_type vec1;
    using vector_allocator_type = typename Allocator::template rebind< subcontainer_type >::other;
    std::vector<subcontainer_type, vector_allocator_type> vecX;
    size_t s;

    // TODO: provide iterator implementation for  begin/end.

    template <typename InputIt>
    InputIt partition_input(InputIt first, InputIt last) {
      return ::std::stable_partition(first, last, splitter);
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

  protected:
    // choices:  sort first, then insert in ranges, or no sort, insert one by one.  second is O(n) but pays the random access and mem realloc cost
    template <class InputIt>
    void insert_impl(InputIt first, InputIt last, supercontainer_type & map) {


        // get previous sizes so we know where to start from
        int64_t idx1 = vec1.size();
        int64_t idxX = vecX.size();

        // reserve but do not yet copy in.  we will have random access to check if already exist,
        // so don't copying in does not save a whole lot.
        vec1.reserve(vec1.size() + std::distance(first, last));

        // iterator over all and insert into map.
        std::pair<typename supercontainer_type::iterator, bool> insert_result;
        Key k;
        int64_t idx;
        for (InputIt it = first, max = last; it != max; ++it) {
          k = (*it).first;

          // try inserting
          insert_result = map.insert(std::make_pair(k, idx1));

          if (insert_result.second) {
            // successful map insertion, so now insert into vec1.
            vec1.emplace_back(*it);
            // new map entry already inserted. just need to increment.
            ++idx1;
          } else {
            // entry already there.  get the iterator and check the idx
            idx = insert_result.first->second;
            if (idx < 0) {
              // previously inserted multiple entries, so reuse the vector, and insert
              vecX[idx & ::std::numeric_limits<int64_t>::max()].emplace_back(*it);
              // no map update is needed.
            } else {
              // previously inserted 1 entry. so
              // update map with new vecX value
              insert_result.first->second = idxX | ~(::std::numeric_limits<int64_t>::max());
              // create a new entry in vecX
              vecX.emplace_back(subcontainer_type());
              // get previous value from vec1 and insert into vecX
              vecX[idxX].emplace_back(std::move(vec1[idx]));
              // insert new value into vecX
              vecX[idxX].emplace_back(*it);
              // update new vecX idx.
              ++idxX;
            }
          }
        }

        s += std::distance(first, last);

        // TODO: compact
    }



    template <typename InputIt, typename Pred>
    void erase_impl(InputIt first, InputIt last, Pred const & pred, supercontainer_type & map) {

      size_t dist = 0;
      int64_t idx;
      // mark for erasure
      for (; first != last; ++first) {
        auto iter = map.find(*first);
        if (iter == map.end()) continue;

        idx = iter->second;

        if (idx < 0) {

          subcontainer_type & vec = vecX[idx & ::std::numeric_limits<int64_t>::max()];

          // multi.
          auto new_end = ::std::remove_if(vec.begin(), vec.end(), pred);

          dist = std::distance(new_end, vec.end());

          if (dist == vec.size()) {
            // remove entry from map.
            map.erase(iter);

            // all deleted, so clean up.
            subcontainer_type().swap(vec);

          } else if (dist + 1 == vec.size()) {
            // update the map
            iter->second = vec1.size();

            // 1 left.  migrate back to vec1
            vec1.emplace_back(std::move(vec[0]));

            // clear the vecX entry
            subcontainer_type().swap(vec);
          } else {
            // clean up, no change to map
            vec.erase(new_end, vec.end());
          }

          s -= dist;

        } else {
          // else single.  if match, erase
          if (pred(vec1[idx])) {
            --s;
            // now clear the map entry.
            map.erase(iter);
          }
        }

      }

    }



    template <typename InputIt>
    void erase_impl(InputIt first, InputIt last, supercontainer_type & map) {

      int64_t idx;
      // mark for erasure
      for (; first != last; ++first) {
        auto iter = map.find(*first);
        if (iter == map.end()) continue;

        idx = iter->second;

        if (idx < 0) {
          // multi.
          subcontainer_type & vec = vecX[idx & ::std::numeric_limits<int64_t>::max()];

          s -= vec.size();
          subcontainer_type().swap(vec);
        } else {
          // else single.  do nothing.
          --s;
        }

        // now clear the map entry.
        map.erase(iter);
      }

    }


    template <typename Pred>
    void erase_impl(Pred const & pred, supercontainer_type & map) {

      size_t dist = 0;
      int64_t idx;


      // mark for erasure
      auto max = map.end();
      for (auto iter = map.begin(); iter != max; ++iter) {
        idx = iter->second;

        if (idx < 0) {

          subcontainer_type & vec = vecX[idx & ::std::numeric_limits<int64_t>::max()];

          // multi.
          auto new_end = ::std::remove_if(vec.begin(), vec.end(), pred);

          dist = std::distance(new_end, vec.end());

          if (dist == vec.size()) {
            // remove entry from map.
            map.erase(iter);

            // all deleted, so clean up.
            subcontainer_type().swap(vec);

          } else if (dist + 1 == vec.size()) {
            // update the map
            iter->second = vec1.size();

            // 1 left.  migrate back to vec1
            vec1.emplace_back(std::move(vec[0]));

            // clear the vecX entry
            subcontainer_type().swap(vec);
          } else {
            // clean up, no change to map
            vec.erase(new_end, vec.end());
          }

          s -= dist;

        } else {
          // else single.  if match, erase
          if (pred(vec1[idx])) {
            --s;
            // now clear the map entry.
            map.erase(iter);
          }
        }
      }
    }


    size_type count_impl(Key const & key, supercontainer_type const & map) const {
      auto iter = map.find(key);
      if (iter == map.end()) return 0;

      if (iter->second < 0) {
        // multiple entries
        return vecX[iter->second &  ::std::numeric_limits<int64_t>::max()].size();
      } else {
        return 1;
      }
    }


    ::std::pair<iterator, iterator> equal_range_impl(Key const & key, supercontainer_type & map) {


      auto iter = map.find(key);

      // not found
      if (iter == map.end()) return ::std::make_pair(vec1.end(), vec1.end());

      // found, has 1 value
      if (iter->second >= 0)  return ::std::make_pair(vec1.begin() + iter->second, vec1.begin() + iter->second + 1);

      // found, has multiple values
      subcontainer_type & vec = vecX[iter->second & ::std::numeric_limits<int64_t>::max()];

      return std::make_pair(vec.begin(), vec.end());

    }
    ::std::pair<const_iterator, const_iterator> equal_range_impl(Key const & key, supercontainer_type const & map) const {
      auto iter = map.find(key);

      // not found
      if (iter == map.end()) return ::std::make_pair(vec1.cend(), vec1.cend());

      // found, has 1 value
      if (iter->second >= 0)  return ::std::make_pair(vec1.cbegin() + iter->second, vec1.cbegin() + iter->second + 1);

      // found, has multiple values
      subcontainer_type & vec = vecX[iter->second & ::std::numeric_limits<int64_t>::max()];

      return std::make_pair(vec.cbegin(), vec.cend());

    }

  public:


//    //  if multiplicity of dataset is kind of known, can initialize data to that to avoid growing vector on average.
//    densehash_multimap(size_type bucket_count = 128,
//                       const Hash& hash = Hash(),
//                       const Equal& equal = Equal(),
//                       const Allocator& alloc = Allocator()) :
//                         lower_map(bucket_count / 2, hash, equal, alloc),
//                         upper_map(bucket_count / 2, hash, Equal(), alloc), s(0UL) {
//    };

    densehash_multimap(Key empty_key, Key deleted_key,
                             Key upper_empty_key, Key upper_deleted_key,
                             Key const & split_key,
                       size_type bucket_count = 128,
                       const Hash& hash = Hash(),
                       const Equal& equal = Equal(),
                       const Allocator& alloc = Allocator()) :
                    	   lower_map(bucket_count / 2, hash, Equal(empty_key, deleted_key), alloc),
						   upper_map(bucket_count / 2, hash, Equal(upper_empty_key, upper_deleted_key), alloc), s(0UL) {
    	reserve_keys(empty_key, deleted_key);
    	reserve_upper_keys(upper_empty_key, upper_deleted_key);
      reserve_split(split_key);

    };

    template<class InputIt>
    densehash_multimap(InputIt first, InputIt last,
                      Key empty_key, Key deleted_key,
                      Key upper_empty_key, Key upper_deleted_key,
                      Key const & split_key,
                       size_type bucket_count = 128,
                       const Hash& hash = Hash(),
                       const Equal& equal = Equal(),
                       const Allocator& alloc = Allocator()) :
                       densehash_multimap(empty_key, deleted_key, upper_empty_key, upper_deleted_key,
                                          split_key, std::distance(first, last), hash, equal, alloc) {
        this->insert(first, last);
    };

    virtual ~densehash_multimap() {
      BL_DEBUG(" split. singles: " << vec1.size() << " multiples: " << vecX.size());
    };

protected:
    void reserve_keys(Key const & empty_key, Key const & deleted_key) {
      lower_map.set_empty_key(empty_key);
      lower_map.set_deleted_key(deleted_key);
    }

    void reserve_upper_keys(Key const & empty_key, Key const & deleted_key) {
      upper_map.set_empty_key(empty_key);
      upper_map.set_deleted_key(deleted_key);
    }

    void reserve_split(Key const & split_key) {
      splitter.upper_bound = split_key;
    }

public:
    // TODO: begin and end iterator accessors.


    std::vector<Key> keys() const  {
      std::vector<Key> ks;

      keys(ks);

      return ks;
    }
    void keys(std::vector<Key> & ks) const  {
      ks.clear();
      ks.reserve(size());

      for (auto it = lower_map.begin(); it != lower_map.end(); ++it) {
        ks.emplace_back(it->first);
      }
      for (auto it = upper_map.begin(); it != upper_map.end(); ++it) {
        ks.emplace_back(it->first);
      }
    }

    std::vector<std::pair<Key, T>> to_vector() const  {
      std::vector<std::pair<Key, T>> vs;

      to_vector(vs);

      return vs;
    }
    void to_vector(  std::vector<std::pair<Key, T>> & vs) const  {
      vs.clear();
      vs.reserve(size());

      for (auto it = lower_map.begin(); it != lower_map.end(); ++it) {
        if (it->second < 0) {
          auto it2 = vecX[it->second & ::std::numeric_limits<int64_t>::max()].begin();
          auto max2 = vecX[it->second & ::std::numeric_limits<int64_t>::max()].end();
          for (; it2 != max2; ++it2) {
            vs.emplace_back(*it2);
          }
        } else {  // singleton
          vs.emplace_back(vec1[it->second]);
        }
      }
      for (auto it = upper_map.begin(); it != upper_map.end(); ++it) {
        if (it->second < 0) {
          auto it2 = vecX[it->second & ::std::numeric_limits<int64_t>::max()].begin();
          auto max2 = vecX[it->second & ::std::numeric_limits<int64_t>::max()].end();
          for (; it2 != max2; ++it2) {
            vs.emplace_back(*it2);
          }
        } else {  // singleton
          vs.emplace_back(vec1[it->second]);
        }
      }
    }


    bool empty() const {
      return lower_map.empty() && upper_map.empty();
    }

    size_type size() const {
      return s;
    }


    size_type unique_size() const {
      return lower_map.size() + upper_map.size();
    }

    void clear() {
      vec1.clear();
      vecX.clear();
      lower_map.clear_no_resize();
      upper_map.clear_no_resize();
      s = 0UL;
    }

    void resize(size_t const n) {
      lower_map.resize(static_cast<float>(n) / 2.0 / lower_map.max_load_factor());
      upper_map.resize(static_cast<float>(n) / 2.0 / lower_map.max_load_factor());
    }

    void rehash(size_type count) {
      this->resize(count);
    }

    /// bucket count.  same as underlying buckets
    size_type bucket_count() { return lower_map.bucket_count() + upper_map.bucket_count(); }

    /// max load factor.  this is the map's max load factor (vectors per bucket) x multiplicity = elements per bucket.  side effect is multiplicity is updated.
    float max_load_factor() {
      return static_cast<float>(size()) / static_cast<float>(bucket_count());
    }





    // choices:  sort first, then insert in ranges, or no sort, insert one by one.  second is O(n) but pays the random access and mem realloc cost
    template <class InputIt>
    void insert(InputIt first, InputIt last) {
        static_assert(::std::is_convertible<std::pair<Key, T>,
                      typename ::std::iterator_traits<InputIt>::value_type>::value,
                      "InputIt value type for insert cannot be converted to std::pair<Key, T> type");

        if (first == last) return;

        auto middle = partition_input(first, last);

        lower_map.resize(static_cast<float>(lower_map.size() + std::distance(first, middle)) / lower_map.max_load_factor());
        insert_impl(first, middle, lower_map);
        upper_map.resize(static_cast<float>(upper_map.size() + std::distance(first, middle)) / upper_map.max_load_factor());
        insert_impl(middle, last, upper_map);

        // TODO: compact
    }

    /// inserting sorted range
    void insert(::std::vector<::std::pair<Key, T> > & input) {
      // TODO: more memory efficient version of this.

      insert(input.begin(), input.end());
    }

    template <typename InputIt, typename Pred>
    size_t erase(InputIt first, InputIt last, Pred const & pred) {
      static_assert(::std::is_convertible<Key, typename ::std::iterator_traits<InputIt>::value_type>::value,
                    "InputIt value type for erase cannot be converted to key type");

      if (first == last) return 0;


      size_t before = s;

      auto middle = partition_input(first, last);
      erase_impl(first, middle, pred, lower_map);
      erase_impl(middle, last, pred, upper_map);

      return before - s;
    }

    template <typename InputIt>
    size_t erase(InputIt first, InputIt last) {
      static_assert(::std::is_convertible<Key, typename ::std::iterator_traits<InputIt>::value_type>::value,
                    "InputIt value type for erase cannot be converted to key type");

      if (first == last) return 0;


      size_t before = s;

      auto middle = partition_input(first, last);
      erase_impl(first, middle, lower_map);
      erase_impl(middle, last, upper_map);

      return before - s;
    }

    template <typename Pred>
    size_t erase(Pred const & pred) {

      if (this->size() == 0) return 0;

      size_t before = s;

      erase_impl(pred, lower_map);
      erase_impl(pred, upper_map);

      return before - s;
    }

    size_type count(Key const & key) const {

      if (splitter(key)) {
        return count_impl(key, lower_map);
      } else {
        return count_impl(key, upper_map);
      }

    }



    ::std::pair<iterator, iterator> equal_range(Key const & key) {
      if (splitter(key)) {
        return equal_range_impl(key, lower_map);
      } else {
        return equal_range_impl(key, upper_map);
      }
    }
    ::std::pair<const_iterator, const_iterator> equal_range(Key const & key) const {
      if (splitter(key)) {
        return equal_range_impl(key, lower_map);
      } else {
        return equal_range_impl(key, upper_map);
      }
    }
    // NO bucket interfaces

};




/**
 * specialization for densehash that does not require key space partitioning
 */
template <typename Key,
typename T,
template <typename > class Transform ,
typename Hash ,
typename Equal,
typename Allocator>
class densehash_multimap<Key, T, false, Transform, Hash, Equal, Allocator> {

  protected:

    // data container
    using tuple_allocator_type = typename Allocator::template rebind< ::std::pair<Key, T> >::other;
    using subcontainer_type = ::std::vector<::std::pair<Key, T>, tuple_allocator_type >;
    using subiter_type = typename subcontainer_type::iterator;
    using const_subiter_type = typename subcontainer_type::const_iterator;

    // index in the vector - this is so we don't have to worry about pointer or iterators being invalidated when vector resizes
    // non-negative values indicate singletons.  negative values indicate multiple entries.
    // multiple entries is stored in a vector of vector.  the index is internal_val_type with sign bit removed.
    using internal_val_type = int64_t;

    using super_allocator_type = typename Allocator::template rebind<std::pair< const Key, internal_val_type> >::other;
    using supercontainer_type =
        ::google::dense_hash_map<Key, internal_val_type,
                       Hash, Equal, super_allocator_type >;

    subcontainer_type vec1;
    using vector_allocator_type = typename Allocator::template rebind< subcontainer_type >::other;
    std::vector<subcontainer_type, vector_allocator_type> vecX;
    supercontainer_type map;
    size_t s;

    // TODO: provide iterator implementation for  begin/end.


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
//    densehash_multimap(size_type bucket_count = 128,
//                       const Hash& hash = Hash(),
//                       const Equal& equal = Equal(),
//                       const Allocator& alloc = Allocator()) :
//                         s(0UL) {
//    };

    densehash_multimap(Key empty_key, Key deleted_key,
                       size_type bucket_count = 128,
                       const Hash& hash = Hash(),
                       const Equal& equal = Equal(),
                       const Allocator& alloc = Allocator()) :
                       	map(bucket_count, hash, Equal(empty_key, deleted_key), alloc), s(0UL) {
    	reserve_keys(empty_key, deleted_key);
    };

    template<class InputIt>
    densehash_multimap(InputIt first, InputIt last,
                      Key empty_key, Key deleted_key,
                       size_type bucket_count = 128,
                       const Hash& hash = Hash(),
                       const Equal& equal = Equal(),
                       const Allocator& alloc = Allocator()) :
                       densehash_multimap(empty_key, deleted_key,
                                          std::distance(first, last), hash, equal, alloc) {
        this->insert(first, last);
    };

    virtual ~densehash_multimap() {
      BL_DEBUG(" unsplit. singles: " << vec1.size() << " multiples: " << vecX.size());
    };


protected:
    void reserve_keys(Key const & empty_key, Key const & deleted_key) {

//    	std::cout << "empty: " << empty_key << " deleted " << deleted_key << std::endl;
//    	std::cout << "empty: " << map.key_eq().empty << " deleted " << map.key_eq().deleted << std::endl;

        map.set_empty_key(empty_key);
        map.set_deleted_key(deleted_key);
    }

    void reserve_upper_keys(Key const & empty_key, Key const & deleted_key, const ::fsc::TruePredicate & _splitter = fsc::TruePredicate()) {
    	//printf("reserve_upper_keys no-op\n");
    }

    void reserve_split(Key const & split_key) {
    }

public:

    // TODO: begin and end iterator accessors.
    std::vector<Key> keys() const  {
      std::vector<Key> ks;

      keys(ks);

      return ks;
    }
    void keys(std::vector<Key> & ks) const  {
      ks.clear();
      ks.reserve(size());

      for (auto it = map.begin(); it != map.end(); ++it) {
        ks.emplace_back(it->first);
      }
    }

    std::vector<std::pair<Key, T>> to_vector() const  {
      std::vector<std::pair<Key, T>> vs;

      to_vector(vs);

      return vs;
    }
    void to_vector(  std::vector<std::pair<Key, T>> & vs) const  {
      vs.clear();
      vs.reserve(size());

      int64_t idx;
      for (auto it = map.begin(); it != map.end(); ++it) {
        if (it->second < 0) {
          idx = it->second & ::std::numeric_limits<int64_t>::max();
          auto max2 = vecX[idx].end();
          for (auto it2 = vecX[idx].begin(); it2 != max2; ++it2) {
            vs.emplace_back(*it2);
          }
        } else {  // singleton
          vs.emplace_back(vec1[it->second]);
        }
      }
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
      vec1.clear();
      vecX.clear();
      map.clear_no_resize();
      s = 0UL;
    }

    void resize(size_t const n) {
      map.resize(static_cast<float>(n) / map.max_load_factor());
    }

    void rehash(size_type count) {
      this->resize(count);
    }

    /// bucket count.  same as underlying buckets
    size_type bucket_count() { return map.bucket_count(); }

    /// max load factor.  this is the map's max load factor (vectors per bucket) x multiplicity = elements per bucket.  side effect is multiplicity is updated.
    float max_load_factor() {
      return static_cast<float>(map.size()) / static_cast<float>(map.bucket_count());
    }



    // choices:  sort first, then insert in ranges, or no sort, insert one by one.  second is O(n) but pays the random access and mem realloc cost
    template <class InputIt>
    void insert(InputIt first, InputIt last) {
        static_assert(::std::is_convertible<std::pair<Key, T>,
                      typename ::std::iterator_traits<InputIt>::value_type>::value,
                      "InputIt value type for insert cannot be converted to std::pair<Key, T> type");

        if (first == last) return;


        // get previous sizes so we know where to start from
        int64_t idx1 = vec1.size();
        int64_t idxX = vecX.size();

        // reserve but do not yet copy in.  we will have random access to check if already exist,
        // so don't copying in does not save a whole lot.
        vec1.reserve(vec1.size() + std::distance(first, last));
        this->resize(vec1.size() + std::distance(first, last));

        // iterator over all and insert into map.
        std::pair<typename supercontainer_type::iterator, bool> insert_result;
        Key k;
        int64_t idx;
        for (InputIt it = first, max = last; it != max; ++it) {
          k = (*it).first;

          // try inserting
          insert_result = map.insert(std::make_pair(k, idx1));

          if (insert_result.second) {
            // successful map insertion, so now insert into vec1.
            vec1.emplace_back(*it);
            // new map entry already inserted. just need to increment.
            ++idx1;
          } else {
            // entry already there.  get the iterator and check the idx
            idx = insert_result.first->second;
            if (idx < 0) {
              // previously inserted multiple entries, so reuse the vector, and insert
              vecX[idx & ::std::numeric_limits<int64_t>::max()].emplace_back(*it);
              // no map update is needed.
            } else {
              // previously inserted 1 entry. so
              // update map with new vecX value
              insert_result.first->second = idxX | ~(::std::numeric_limits<int64_t>::max());
              // create a new entry in vecX
              vecX.emplace_back(subcontainer_type());
              // get previous value from vec1 and insert into vecX
              vecX[idxX].emplace_back(std::move(vec1[idx]));
              // insert new value into vecX
              vecX[idxX].emplace_back(*it);
              // update new vecX idx.
              ++idxX;
            }
          }
        }

        s += std::distance(first, last);

        // TODO: compact
    }

    /// insert
    void insert(::std::vector<::std::pair<Key, T> > & input) {

      if (input.size() == 0) return;
        // there is data in there already.  so do normal insert
        this->insert(input.begin(), input.end());
    }

    template <typename InputIt, typename Pred>
    size_t erase(InputIt first, InputIt last, Pred const & pred) {
      static_assert(::std::is_convertible<Key, typename ::std::iterator_traits<InputIt>::value_type>::value,
                    "InputIt value type for erase cannot be converted to key type");

      if (first == last) return 0;

      size_t before = s;
      size_t dist = 0;
      int64_t idx;
      // mark for erasure
      for (; first != last; ++first) {
        auto iter = map.find(*first);
        if (iter == map.end()) continue;

        idx = iter->second;

        if (idx < 0) {

          subcontainer_type & vec = vecX[idx & ::std::numeric_limits<int64_t>::max()];

          // multi.
          auto new_end = ::std::remove_if(vec.begin(), vec.end(), pred);

          dist = std::distance(new_end, vec.end());

          if (dist == vec.size()) {
            // remove entry from map.
            map.erase(iter);

            // all deleted, so clean up.
            subcontainer_type().swap(vec);

          } else if (dist + 1 == vec.size()) {
            // update the map
            iter->second = vec1.size();

            // 1 left.  migrate back to vec1
            vec1.emplace_back(std::move(vec[0]));

            // clear the vecX entry
            subcontainer_type().swap(vec);
          } else {
            // clean up, no change to map
            vec.erase(new_end, vec.end());
          }

          s -= dist;

        } else {
          // else single.  if match, erase
          if (pred(vec1[idx])) {
            --s;
            // now clear the map entry.
            map.erase(iter);
          }
        }

      }

      return before - s;
    }

    template <typename InputIt>
    size_t erase(InputIt first, InputIt last) {
      static_assert(::std::is_convertible<Key, typename ::std::iterator_traits<InputIt>::value_type>::value,
                    "InputIt value type for erase cannot be converted to key type");

      if (first == last) return 0;

      size_t before = s;

      int64_t idx;
      // mark for erasure
      for (; first != last; ++first) {
        auto iter = map.find(*first);
        if (iter == map.end()) continue;

        idx = iter->second;

        if (idx < 0) {
          // multi.
          subcontainer_type & vec = vecX[idx & ::std::numeric_limits<int64_t>::max()];

          s -= vec.size();
          subcontainer_type().swap(vec);
        } else {
          // else single.  do nothing.
          --s;
        }

        // now clear the map entry.
        map.erase(iter);
      }

      return before - s;
    }

    template <typename Pred>
    size_t erase(Pred const & pred) {

      if (this->size() == 0) return 0;

      size_t before = s;
      size_t dist = 0;
      int64_t idx;


      // mark for erasure
      auto  max = map.end();
      for (auto iter = map.begin(); iter != max; ++iter) {
        idx = iter->second;

        if (idx < 0) {

          subcontainer_type & vec = vecX[idx & ::std::numeric_limits<int64_t>::max()];

          // multi.
          auto new_end = ::std::remove_if(vec.begin(), vec.end(), pred);

          dist = std::distance(new_end, vec.end());

          if (dist == vec.size()) {
            // remove entry from map.
            map.erase(iter);

            // all deleted, so clean up.
            subcontainer_type().swap(vec);

          } else if (dist + 1 == vec.size()) {
            // update the map
            iter->second = vec1.size();

            // 1 left.  migrate back to vec1
            vec1.emplace_back(std::move(vec[0]));

            // clear the vecX entry
            subcontainer_type().swap(vec);
          } else {
            // clean up, no change to map
            vec.erase(new_end, vec.end());
          }

          s -= dist;

        } else {
          // else single.  if match, erase
          if (pred(vec1[idx])) {
            --s;
            // now clear the map entry.
            map.erase(iter);
          }
        }
      }

      return before - s;
    }

    size_type count(Key const & key) const {
      auto iter = map.find(key);
      if (iter == map.end()) return 0;

      if (iter->second < 0) {
        // multiple entries
        return vecX[iter->second &  ::std::numeric_limits<int64_t>::max()].size();
      } else {
        return 1;
      }
    }



    ::std::pair<iterator, iterator> equal_range(Key const & key) {
      auto iter = map.find(key);

      // not found
      if (iter == map.end()) return ::std::make_pair(vec1.end(), vec1.end());

      // found, has 1 value
      if (iter->second >= 0)  return ::std::make_pair(vec1.begin() + iter->second, vec1.begin() + iter->second + 1);

      // found, has multiple values
      subcontainer_type & vec = vecX[iter->second & ::std::numeric_limits<int64_t>::max()];

      return std::make_pair(vec.begin(), vec.end());

    }
    ::std::pair<const_iterator, const_iterator> equal_range(Key const & key) const {
      auto iter = map.find(key);

      // not found
      if (iter == map.end()) return ::std::make_pair(vec1.cend(), vec1.cend());

      // found, has 1 value
      if (iter->second >= 0)  return ::std::make_pair(vec1.cbegin() + iter->second, vec1.cbegin() + iter->second + 1);

      // found, has multiple values
      int64_t idx = iter->second & ::std::numeric_limits<int64_t>::max();

      return std::make_pair(vecX[idx].cbegin(), vecX[idx].cend());


    }
    // NO bucket interfaces

};


} // end namespace fsc.




#endif /* SRC_WIP_DENSEHASH_VECMAP_HPP_ */
