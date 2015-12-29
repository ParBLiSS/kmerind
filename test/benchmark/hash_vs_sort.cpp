/**
 * @file    hashVsSort.cpp

 * @ingroup
 * @author  tpan
 * @brief   test the scalability of sort and search vs hashing.
 * @details  traditionally, hash is used for O(1) access.  sort is assume to take O(log(N)) time on sorted list.
 *          benefit of hashmap is that different hashes can be used to improve locality or load balance, drawback is load imbalance
 *          benefit of sort is that load is balanced for different processors. adjacent values are on same processor.
 *            draw back is log(N) access for each query
 *
 *          any transformation of values has additional cost.
 *
 *          when log(N) become greater than hash function time constant, hashmap wins.  otherwise (small N) sort wins.
 *          std::sort of uint64_t:  break even at just under 8M entries in vector.  this is a pretty large count.
 *
 *          what can be done:
 *          1. choose good hash function to improve load balance.   cost:  hash time per element during build and query
 *          2. 2 or more levels of indirection in sorting - break sorted array into buckets and search within buckets.
 *              does not change complexity.
 *          3. create buckets based on hash key, assign uneven number of buckets to processors to load balance.  use
 *              lookup table to map hash value to processor.  drawback:  could still have badly behaving buckets.
 *              can't do this progressively without some rebalancing.
 *
 *
 *
 *
 *
 * Copyright (c) 2015 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */


#include "farmhash/src/farmhash.cc"
#include <vector>
#include <unordered_map>
#include <algorithm>  // sort, equal_range,
#include <cstdlib> // atoi, srand, rand
#include <utility> // pair
#include <chrono>  // clock, time_point, duration.
#include <cstdio>  // printf.

#define TIMER_INIT(session)      std::chrono::steady_clock::time_point session##_t1, session##_t2; \
                                 std::chrono::duration<double> session##_time_span

#define TIMER_START(session)     session##_t1 = std::chrono::steady_clock::now()
#define TIMER_END(session)       session##_t2 = std::chrono::steady_clock::now(); \
                                 session##_time_span = std::chrono::duration_cast<std::chrono::duration<double>>(session##_t2 - session##_t1);
#define TIMER_REPORT(name, size, session) printf("%ld %s duration %f\n", size, name, session##_time_span.count())


template <typename K>
class FarmHash {
  public:
    /// operator to compute hash
    size_t operator()(const K & key) const {
      return ::util::Hash(reinterpret_cast<const char*>(&key), sizeof(K));
    }
};


template <typename K, typename V, typename H>
void reserve(std::unordered_multimap<K, V, H> &map, const size_t entries) {
  map.reserve(entries);
}

template <typename K, typename V, typename H>
void insert(std::unordered_multimap<K, V, H> &map, const size_t entries, unsigned int rand_seed = 1) {

  srand(rand_seed);
  int key = 0;
  for (size_t i = 0; i < entries; ++i) {
    key = rand();
    map.emplace(key, key);
  }
}

template <typename K, typename V, typename H>
void query(const std::unordered_multimap<K, V, H> &map, const size_t entries, unsigned int rand_seed = 1) {

  srand(rand_seed);
  volatile int val = 0;
  for (size_t i = 0; i < entries; ++i) {
    auto range = map.equal_range(rand());  // should be a hash to get bucket.  then what happens within bucket?
    for (auto it = range.first; it != range.second; ++it) {
      val = it->second;
    }
  }
  srand(val);
}



template <typename K, typename V>
struct key_less_than {
  bool operator()(const std::pair<K, V> & x, const std::pair<K, V> & y) {
    return x.first < y.first;
  }
  bool operator()(const std::pair<K, V> & x, const K & Ky) {
    return x.first < Ky;
  }
  bool operator()(const K & Kx, const std::pair<K, V> & y) {
    return Kx < y.first;
  }
};


template <typename K, typename V>
void reserve(std::vector<std::pair<K, V> > &vec, const size_t entries) {
  vec.reserve(entries);
}

template <typename K, typename V>
void insert(std::vector<std::pair<K, V> > &vec, const size_t entries, unsigned int rand_seed = 1) {

  srand(rand_seed);
  int key = 0;
  for (size_t i = 0; i < entries; ++i) {
    key = rand();
    vec.emplace_back(key, key);
  }
}

template <typename K, typename V>
void sort(std::vector<std::pair<K, V> > &vec) {
  std::sort(vec.begin(), vec.end(), key_less_than<K, V>());
}


template <typename K, typename V>
void query(const std::vector<std::pair<K, V> > &vec, const size_t entries, unsigned int rand_seed = 1) {

  key_less_than<K, V> comparator;

  srand(rand_seed);
  volatile int val = 0;
  for (size_t i = 0; i < entries; ++i) {
    auto range = std::equal_range(vec.begin(), vec.end(), rand(), comparator);  // a binary search each time.
    for (auto it = range.first; it != range.second; ++it) {
      val = it->second;
    }
  }
  srand(val);
}

/// batched query allows for query terms to be sorted, then searches can be sped up.
/// use middle element to split up things.
template <typename K, typename V>
void batched_query(const std::vector<std::pair<K, V> > &vec, const size_t entries, std::vector<K> &batch, const size_t batch_size, unsigned int rand_seed = 1) {

  key_less_than<K, V> comparator;

  srand(rand_seed);
  volatile int val = 0;
  auto iter = vec.begin();

  for (size_t i = 0; i < entries; i += batch_size) {

    for (size_t j = 0; j < batch_size; ++j)
    {
      batch[j] = rand();
    }
    std::sort(batch.begin(), batch.end());

    iter = vec.begin();
    for (size_t j = 0; j < batch_size; ++j)
    {
      auto range = std::equal_range(iter, vec.end(), batch[j], comparator);  // a binary search each time.
      for (auto it = range.first; it != range.second; ++it) {
        val = it->second;
      }
      iter = range.second;
    }
  }
  srand(val);
}



int main(int argc, char** argv) {
  using KeyType = uint64_t;
  using ValType = int64_t;
  using PairType = std::pair<KeyType, ValType>;  // typical tuple for kmers   16 bytes

  // param:  max exponent for 2^exp.  from 1 to max.  max default to space of 4GB, so 256M entries, so max = 28

  uint8_t max_exp = 27;
  if (argc > 1) {
    max_exp = atoi(argv[1]);
  }

  size_t max_size = 0x1 << max_exp;

  size_t size = 0x1;

  TIMER_INIT(hash);
  TIMER_INIT(sort);
  TIMER_INIT(hash_total);
  TIMER_INIT(sort_total);

  for (; size < max_size; size <<= 1)
  {
    {
      TIMER_START(hash_total);


      TIMER_START(hash);
      std::unordered_multimap<KeyType, ValType, FarmHash<KeyType> > container;
      //============ do the hash test
      // reserve size
      reserve(container, size);
      TIMER_END(hash);
      TIMER_REPORT("hash reserve", size, hash);

      // insert into hash
      TIMER_START(hash);
      insert(container, size);
      TIMER_END(hash);
      TIMER_REPORT("hash insert", size, hash);

      TIMER_END(hash_total);
      TIMER_REPORT("hash BUILD TOTAL", size, hash_total);

      TIMER_START(hash_total);


      // query hash, exactly the same items
      TIMER_START(hash);
      query(container, size);
      TIMER_END(hash);
      TIMER_REPORT("hash query all in map", size, hash);

      // query hash, different items
      TIMER_START(hash);
      query(container, size, 1234);
      TIMER_END(hash);
      TIMER_REPORT("hash query random", size, hash);


      TIMER_END(hash_total);
      TIMER_REPORT("hash QUERY TOTAL", size, hash_total);

    }  // clear the memory.

    //============ do the vector test
    {
      TIMER_START(sort_total);

      TIMER_START(sort);
      std::vector<PairType> container;
      //============ do the sort test
      // reserve size
      reserve(container, size);
      TIMER_END(sort);
      TIMER_REPORT("sort reserve", size, sort);

      // insert into hash
      TIMER_START(sort);
      insert(container, size);
      TIMER_END(sort);
      TIMER_REPORT("sort insert", size, sort);

      // sort the container
      TIMER_START(sort);
      sort(container);
      TIMER_END(sort);
      TIMER_REPORT("sort sorting", size, sort);

      TIMER_END(sort_total);
      TIMER_REPORT("sort BUILD TOTAL", size, sort_total);

      TIMER_START(sort_total);

      // query hash, exactly the same items
      TIMER_START(sort);
      query(container, size);
      TIMER_END(sort);
      TIMER_REPORT("sort query all in map", size, sort);

      // query hash, different items
      TIMER_START(sort);
      query(container, size, 1234);
      TIMER_END(sort);
      TIMER_REPORT("sort query random", size, sort);

      TIMER_END(sort_total);
      TIMER_REPORT("sort QUERY TOTAL", size, sort_total);

      TIMER_START(sort_total);

      size_t batch_size = 128*1024;
      std::vector<KeyType> keys(batch_size);

      // query hash, exactly the same items
      TIMER_START(sort);
      batched_query(container, size, keys, batch_size);
      TIMER_END(sort);
      TIMER_REPORT("batched sort query all in map", size, sort);

      // query hash, different items
      TIMER_START(sort);
      batched_query(container, size, keys, batch_size, 1234);
      TIMER_END(sort);
      TIMER_REPORT("batched sort query random", size, sort);

      TIMER_END(sort_total);
      TIMER_REPORT("batched sort QUERY TOTAL", size, sort_total);

    }  // clear the memory.


    //============ report timing


    //============ increment the size
  }


}
