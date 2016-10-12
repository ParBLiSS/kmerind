#include <unordered_map>
#include <vector>
#include <random>
#include <cstdint>

#include <tommyds/tommyalloc.h>
#include <tommyds/tommyalloc.c>
#include <tommyds/tommyhashdyn.h>
#include <tommyds/tommyhashdyn.c>
#include <tommyds/tommyhashlin.h>
#include <tommyds/tommyhashlin.c>
#include <tommyds/tommytrie.h>
#include <tommyds/tommytrie.c>


#include "containers/unordered_vecmap.hpp"
//#include "containers/hashed_vecmap.hpp"
#include "containers/densehash_map.hpp"

#include "common/kmer.hpp"
#include "common/kmer_transform.hpp"
#include "index/kmer_hash.hpp"


#include "mxx/env.hpp"
#include "mxx/comm.hpp"

#include "utils/benchmark_utils.hpp"
#include "utils/transform_utils.hpp"

// comparison of some hash tables.  note that this is not exhaustive and includes only the well tested ones and my own.  not so much
// the one-off ones people wrote.
// see http://preshing.com/20110603/hash-table-performance-tests/  - suggests google sparsehash dense, and Judy array
//      http://incise.org/hash-table-benchmarks.html  - suggests google dense hash map and glib ghashtable
//      https://attractivechaos.wordpress.com/2008/08/28/comparison-of-hash-table-libraries/  - suggests google sparsehash dense and khash (distant second)
//      http://preshing.com/20130107/this-hash-table-is-faster-than-a-judy-array/  - suggests judy array
//      http://www.tommyds.it/doc/index.html  - suggets Tommy_hashtable and google dense.  at the range we operate, Tommy and Google Densehash are competitive.
//      http://www.nothings.org/computer/judy/ - shows that judy performs poorly with random data insertion. sequential is better (so sorted kmers would be better)

// hashedvec is very slow because of repeated unordered map traversal.
// same with unordered vecmap
// unordered multimap is very slow relative to google dense

// google dense requires an empty key and a deleted key. it is definitely fast.
//      it does not support multimap...
// judy array is a map from integer to integer (treat as pointer?)  - essentially a base256 radix tree (hierarchical index).  64byte cacheline systems only
//      does NOT support multimap natively.
//      limit to 4B entries.
// Tommy provides mapping from integer to object pointers.  hash is computed by user (integer) as is the object pointer.
//          Tommy's limitation is there can be at most 2^32 -1 entries in the container - i.e. about 4Billion elements.  it probably is not an issue
//          Another of Tommy's issue is that it only index pointers so there is significant wrapper required, for example, hash value computation is done by user.
//          for smaller counts (< 1M entries) it's faster than google dense map - but experiments show that google dense hash is still faster.
//              can potentially avoid copying the input vector.
//      natively is a multimap.
//      tommytrie can potentially be used directly on k-mers - at least, on the Most Significant word.
//      tommyhashdyn or hashlin can be used for hashed version.
//      tommytrie segv deleting an already deleted entry, so need to be careful.

// results:  google dense hash is fastest.  tommy is pretty fast (about 2x faster), especially when using hashdyn.
// vecmap is not correct.  and tommytrie core dumps when n = 100M.  sparsehash is over 2x faster at insetion and count than tommy.
// question is how to make google dense hash support multimap style operations?  vector is expensive...



template <typename Kmer, typename Value>
void generate_input(std::vector<::std::pair<Kmer, Value> > & output, size_t const count, bool canonical = false) {
  output.resize(count);

  srand(23);
  for (size_t i = 0; i < count; ++i) {
    for (size_t j = 0; j < Kmer::nWords; ++j) {
      output[i].first.getDataRef()[j] = static_cast<typename Kmer::KmerWordType>(static_cast<long>(rand()) << 32) ^ static_cast<long>(rand());
    }
    output[i].first.sanitize();
    //output[i].second = static_cast<Value>(static_cast<long>(rand()) << 32) ^ static_cast<long>(rand());
    output[i].second = i;
    if (rand() < (1 << 24)) {
      // make a repeat
      ++i;
      if (i < count) {
        output[i] = output[i-1];
        //output[i].second = static_cast<Value>(static_cast<long>(rand()) << 32) ^ static_cast<long>(rand());
        output[i].second = i;
      }
    }
    if (rand() < (1 << 16)) {
      // make a repeat
      ++i;
      if (i < count) {

        output[i] = output[i-1];
        //output[i].second = static_cast<Value>(static_cast<long>(rand()) << 32) ^ static_cast<long>(rand());
        output[i].second = i;
      }
    }
    if (rand() < (1 << 8)) {
      // make a repeat
      ++i;
      if (i < count) {
        output[i] = output[i-1];
        //output[i].second = static_cast<Value>(static_cast<long>(rand()) << 32) ^ static_cast<long>(rand());
        output[i].second = i;
      }
    }
  }

  if (canonical) {
	  for (size_t i = 0; i < output.size(); ++i) {
		  Kmer revcomp = output[i].first.reverse_complement();
		  if (revcomp < output[i].first) output[i].first = revcomp;
	  }
  }
}

template <typename Kmer, typename Value>
void benchmark_unordered_map(size_t const count, size_t const query_frac, ::mxx::comm const & comm) {
  BL_BENCH_INIT(map);

  std::vector<Kmer> query;

  BL_BENCH_START(map);
  // no transform involved
  ::std::unordered_map<Kmer, Value, ::bliss::kmer::hash::farm<Kmer, false> > map(count);
  BL_BENCH_END(map, "reserve", count);

  {
//    BL_BENCH_START(map);
    std::vector<::std::pair<Kmer, Value> > input(count);
//    BL_BENCH_END(map, "reserve input", count);

//    BL_BENCH_START(map);
    generate_input(input, count);
    query.resize(count / query_frac);
    std::transform(input.begin(), input.begin() + input.size() / query_frac, query.begin(),
                   [](::std::pair<Kmer, Value> const & x){
      return x.first;
    });
//    BL_BENCH_END(map, "generate input", input.size());

    BL_BENCH_START(map);
    map.insert(input.begin(), input.end());
    BL_BENCH_END(map, "insert", map.size());
  }


  BL_BENCH_START(map);
  size_t result = 0;
  for (size_t i = 0, max = count / query_frac; i < max; ++i) {
    auto iters = map.equal_range(query[i]);
    for (auto it = iters.first; it != iters.second; ++it)
      result ^= it->second;
  }
  BL_BENCH_END(map, "find", result);

  BL_BENCH_START(map);
  result = 0;
  for (size_t i = 0, max = count / query_frac; i < max; ++i) {
    result += map.count(query[i]);
  }
  BL_BENCH_END(map, "count", result);


  BL_BENCH_START(map);
  result = 0;
  for (size_t i = 0, max = count / query_frac; i < max; ++i) {
    result += map.erase(query[i]);
  }
  BL_BENCH_END(map, "erase", result);

  BL_BENCH_REPORT_MPI_NAMED(map, "unordered_map", comm);
}



template <typename Kmer, typename Value>
void benchmark_unordered_multimap(size_t const count, size_t const query_frac, ::mxx::comm const & comm) {
  BL_BENCH_INIT(map);

  std::vector<Kmer> query;

  BL_BENCH_START(map);
  // no transform involved
  ::std::unordered_multimap<Kmer, Value, ::bliss::kmer::hash::farm<Kmer, false> > map(count);
  BL_BENCH_END(map, "reserve", count);


  {
//    BL_BENCH_START(map);
    std::vector<::std::pair<Kmer, Value> > input(count);
//    BL_BENCH_END(map, "reserve input", count);

//    BL_BENCH_START(map);
    generate_input(input, count);
    query.resize(count / query_frac);
    std::transform(input.begin(), input.begin() + input.size() / query_frac, query.begin(),
                   [](::std::pair<Kmer, Value> const & x){
      return x.first;
    });
//    BL_BENCH_END(map, "generate input", input.size());



    BL_BENCH_START(map);
    map.insert(input.begin(), input.end());
    BL_BENCH_END(map, "insert", map.size());
  }


  BL_BENCH_START(map);
  size_t result = 0;
  for (size_t i = 0, max = count / query_frac; i < max; ++i) {
    auto iters = map.equal_range(query[i]);
    for (auto it = iters.first; it != iters.second; ++it)
      result ^= it->second;
  }
  BL_BENCH_END(map, "find", result);

  BL_BENCH_START(map);
  result = 0;
  for (size_t i = 0, max = count / query_frac; i < max; ++i) {
    result += map.count(query[i]);
  }
  BL_BENCH_END(map, "count", result);


  BL_BENCH_START(map);
  result = 0;
  for (size_t i = 0, max = count / query_frac; i < max; ++i) {
    result += map.erase(query[i]);
  }
  BL_BENCH_END(map, "erase", result);

  BL_BENCH_REPORT_MPI_NAMED(map, "unordered_multimap", comm);
}



template <typename Kmer, typename Value>
void benchmark_unordered_vecmap(size_t const count, size_t const query_frac, ::mxx::comm const & comm) {
  BL_BENCH_INIT(map);

  std::vector<Kmer > query;
  BL_BENCH_START(map);
  // no transform involved.
  ::fsc::unordered_vecmap<Kmer, Value, ::bliss::kmer::hash::farm<Kmer, false> > map(1, count);
  BL_BENCH_END(map, "reserve", count);


  {

//    BL_BENCH_START(map);
    std::vector<::std::pair<Kmer, Value> > input(count);
//    BL_BENCH_END(map, "reserve input", count);

//    BL_BENCH_START(map);
    generate_input(input, count);
    query.resize(count / query_frac);
    std::transform(input.begin(), input.begin() + input.size() / query_frac, query.begin(),
                   [](::std::pair<Kmer, Value> const & x){
      return x.first;
    });
//    BL_BENCH_END(map, "generate input", input.size());

    BL_BENCH_START(map);
    map.insert_sorted(input.begin(), input.end());
    BL_BENCH_END(map, "insert", map.size());
  }

  BL_BENCH_START(map);
  size_t result = 0;
  for (size_t i = 0, max = count / query_frac; i < max; ++i) {
    auto iters = map.equal_range(query[i]);
    for (auto it = iters.first; it != iters.second; ++it)
      result ^= (*it).second;
  }
  BL_BENCH_END(map, "find", result);

  BL_BENCH_START(map);
  result = 0;
  for (size_t i = 0, max = count / query_frac; i < max; ++i) {
    result += map.count(query[i]);
  }
  BL_BENCH_END(map, "count", result);

  BL_BENCH_START(map);
  result = 0;
  for (size_t i = 0, max = count / query_frac; i < max; ++i) {
    result += map.erase(query[i]);
  }
  BL_BENCH_END(map, "erase", result);

  BL_BENCH_REPORT_MPI_NAMED(map, "unordered_vecmap", comm);
}

//template <typename Kmer, typename Value>
//void benchmark_hashed_vecmap(size_t const count, size_t const query_frac, ::mxx::comm const & comm) {
//  BL_BENCH_INIT(map);
//
//  std::vector<Kmer> query;
//
//  BL_BENCH_START(map);
//  // no transform involved.
//  ::fsc::hashed_vecmap<Kmer, Value, ::bliss::kmer::hash::farm<Kmer, false> > map(1, count);
//  BL_BENCH_END(map, "reserve", count);
//
//
//  {
////    BL_BENCH_START(map);
//    std::vector<::std::pair<Kmer, Value> > input(count);
////    BL_BENCH_END(map, "reserve input", count);
//
////    BL_BENCH_START(map);
//    generate_input(input, count);
//    query.resize(count / query_frac);
//    std::transform(input.begin(), input.begin() + input.size() / query_frac, query.begin(),
//                   [](::std::pair<Kmer, Value> const & x){
//      return x.first;
//    });
////    BL_BENCH_END(map, "generate input", input.size());
//
//
//    BL_BENCH_START(map);
//    map.insert(input);
//    BL_BENCH_END(map, "insert", map.size());
//  }
//
//  BL_BENCH_START(map);
//  size_t result = 0;
//  for (size_t i = 0, max = count / query_frac; i < max; ++i) {
//    auto iters = map.equal_range(query[i]);
//    for (auto it = iters.first; it != iters.second; ++it)
//      result ^= (*it).second;
//  }
//  BL_BENCH_END(map, "find", result);
//
//
//  BL_BENCH_START(map);
//  result = 0;
//  for (size_t i = 0, max = count / query_frac; i < max; ++i) {
//    result += map.count(query[i]);
//  }
//  BL_BENCH_END(map, "count", result);
//
//  BL_BENCH_START(map);
//  result = map.erase(query.begin(), query.end());
//  BL_BENCH_END(map, "erase", result);
//
//
//  BL_BENCH_REPORT_MPI_NAMED(map, "hashed_vecmap", comm);
//}

template <typename Kmer, typename Value>
void benchmark_densehash_map(size_t const count, size_t const query_frac, ::mxx::comm const & comm) {
  BL_BENCH_INIT(map);

  std::vector<Kmer> query;

  BL_BENCH_START(map);
  // no transform involved.
  ::fsc::densehash_map<Kmer, Value, 
	::bliss::kmer::hash::sparsehash::special_keys<Kmer, false>,
	::bliss::transform::identity,
	::bliss::kmer::hash::farm<Kmer, false> > map(count);
  BL_BENCH_END(map, "reserve", count);



  {
//    BL_BENCH_START(map);
    std::vector<::std::pair<Kmer, Value> > input(count);
//    BL_BENCH_END(map, "reserve input", count);

//    BL_BENCH_START(map);
    generate_input(input, count);
    query.resize(count / query_frac);
    std::transform(input.begin(), input.begin() + input.size() / query_frac, query.begin(),
                   [](::std::pair<Kmer, Value> const & x){
      return x.first;
    });
//    BL_BENCH_END(map, "generate input", input.size());

    BL_BENCH_START(map);
    map.insert(input.begin(), input.end());
    BL_BENCH_END(map, "insert", map.size());
  }

  BL_BENCH_START(map);
  size_t result = 0;
  for (size_t i = 0, max = count / query_frac; i < max; ++i) {
    auto iters = map.equal_range(query[i]);
    for (auto it = iters.first; it != iters.second; ++it)
      result ^= it->second;
  }
  BL_BENCH_END(map, "find", result);

  BL_BENCH_START(map);
  result = 0;
  for (size_t i = 0, max = count / query_frac; i < max; ++i) {
    result += map.count(query[i]);
  }
  BL_BENCH_END(map, "count", result);

  BL_BENCH_START(map);
  result = map.erase(query.begin(), query.end());

//  for (size_t i = 0, max = count / query_frac; i < max; ++i) {
//    result += map.erase(query[i]);
//  }
  map.resize(0);
  BL_BENCH_END(map, "erase", result);


  BL_BENCH_REPORT_MPI_NAMED(map, "densehash_map", comm);
}


template <typename Kmer, typename Value, bool canonical = false>
void benchmark_densehash_full_map(size_t const count, size_t const query_frac, ::mxx::comm const & comm) {
  BL_BENCH_INIT(map);

  std::vector<Kmer> query;

  BL_BENCH_START(map);
  // no transform involved.
  ::fsc::densehash_map<Kmer, Value, 
	::bliss::kmer::hash::sparsehash::special_keys<Kmer, canonical>,
	::bliss::transform::identity,
	::bliss::kmer::hash::farm<Kmer, false> > map(count);

  BL_BENCH_END(map, "reserve", count);



  {
//    BL_BENCH_START(map);
    std::vector<::std::pair<Kmer, Value> > input(count);
//    BL_BENCH_END(map, "reserve input", count);

//    BL_BENCH_START(map);
    generate_input(input, count, canonical);
    query.resize(count / query_frac);
    std::transform(input.begin(), input.begin() + input.size() / query_frac, query.begin(),
                   [](::std::pair<Kmer, Value> const & x){
      return x.first;
    });
//    BL_BENCH_END(map, "generate input", input.size());

    BL_BENCH_START(map);
    map.insert(input.begin(), input.end());
    BL_BENCH_END(map, "insert", map.size());
  }

  BL_BENCH_START(map);
  size_t result = 0;
  for (size_t i = 0, max = count / query_frac; i < max; ++i) {
    auto iters = map.equal_range(query[i]);
    for (auto it = iters.first; it != iters.second; ++it)
      result ^= (*it).second;
  }
  BL_BENCH_END(map, "find", result);

  BL_BENCH_START(map);
  result = 0;
  for (size_t i = 0, max = count / query_frac; i < max; ++i) {
    result += map.count(query[i]);
  }
  BL_BENCH_END(map, "count", result);

  BL_BENCH_START(map);
  result = map.erase(query.begin(), query.end());

//  for (size_t i = 0, max = count / query_frac; i < max; ++i) {
//    result += map.erase(query[i]);
//  }
  map.resize(0);
  BL_BENCH_END(map, "erase", result);


  BL_BENCH_REPORT_MPI_NAMED(map, "densehash_full_map", comm);
}

/*
template <typename Kmer, typename Value>
void benchmark_densehash_vecmap(size_t const count, size_t const query_frac, ::mxx::comm const & comm) {
  BL_BENCH_INIT(map);

  std::vector<Kmer> query;

  BL_BENCH_START(map);
  Kmer empty(true);
  empty.getDataRef()[Kmer::nWords - 1] = ~(~(static_cast<typename Kmer::KmerWordType>(0)) >> 1);  // for now, works for odd k values
  Kmer deleted(true);
  deleted.getDataRef()[Kmer::nWords - 1] = empty.getDataRef()[Kmer::nWords - 1] >> 1;  // for now, works for odd k values

  ::fsc::densehash_vecmap<Kmer, Value, ::bliss::kmer::hash::farm<Kmer, false> > map(empty, deleted, count);
  BL_BENCH_END(map, "reserve", count);


  {

//    BL_BENCH_START(map);
    std::vector<::std::pair<Kmer, Value> > input(count);
//    BL_BENCH_END(map, "reserve input", count);

//    BL_BENCH_START(map);
    generate_input(input, count);
    query.resize(count / query_frac);
    std::transform(input.begin(), input.begin() + input.size() / query_frac, query.begin(),
                   [](::std::pair<Kmer, Value> const & x){
      return x.first;
    });
//    BL_BENCH_END(map, "generate input", input.size());

    BL_BENCH_START(map);
    map.insert(input);
    BL_BENCH_END(map, "insert", map.size());
  }

  BL_BENCH_START(map);
  size_t result = 0;
  for (size_t i = 0, max = count / query_frac; i < max; ++i) {
    auto iters = map.equal_range(query[i]);
    for (auto it = iters.first; it != iters.second; ++it)
      result ^= (*it).second;
  }
  BL_BENCH_END(map, "find", result);


  BL_BENCH_START(map);
  result = 0;
  for (size_t i = 0, max = count / query_frac; i < max; ++i) {
    result += map.count(query[i]);
  }
  BL_BENCH_END(map, "count", result);



  BL_BENCH_START(map);
  result = map.erase(query.begin(), query.end());
  BL_BENCH_END(map, "erase", result);


  BL_BENCH_REPORT_MPI_NAMED(map, "densehash_vecmap", comm);
}
*/



template <typename Kmer, typename Value>
void benchmark_densehash_multimap(size_t const count, size_t const query_frac, ::mxx::comm const & comm) {
  BL_BENCH_INIT(map);

  std::vector<Kmer> query;

  BL_BENCH_START(map);
  // no transform involved.
  ::fsc::densehash_multimap<Kmer, Value, 
	::bliss::kmer::hash::sparsehash::special_keys<Kmer, false>,
	::bliss::transform::identity,
	::bliss::kmer::hash::farm<Kmer, false> > map(count);

  BL_BENCH_END(map, "reserve", count);


  {

//    BL_BENCH_START(map);
    std::vector<::std::pair<Kmer, Value> > input(count);
//    BL_BENCH_END(map, "reserve input", count);

//    BL_BENCH_START(map);
    generate_input(input, count);
    query.resize(count / query_frac);
    std::transform(input.begin(), input.begin() + input.size() / query_frac, query.begin(),
                   [](::std::pair<Kmer, Value> const & x){
      return x.first;
    });
//    BL_BENCH_END(map, "generate input", input.size());

    BL_BENCH_START(map);
    map.insert(input);
    BL_BENCH_END(map, "insert", map.size());
  }

  BL_BENCH_START(map);
  size_t result = 0;
  for (size_t i = 0, max = count / query_frac; i < max; ++i) {
    auto iters = map.equal_range(query[i]);
    for (auto it = iters.first; it != iters.second; ++it)
      result ^= (*it).second;
  }
  BL_BENCH_END(map, "find", result);


  BL_BENCH_START(map);
  result = 0;
  for (size_t i = 0, max = count / query_frac; i < max; ++i) {
    result += map.count(query[i]);
  }
  BL_BENCH_END(map, "count", result);



  BL_BENCH_START(map);
  result = map.erase(query.begin(), query.end());
  BL_BENCH_END(map, "erase", result);


  BL_BENCH_REPORT_MPI_NAMED(map, "densehash_multimap", comm);
}

template <typename Kmer, typename Value, bool canonical = false>
void benchmark_densehash_full_multimap(size_t const count, size_t const query_frac, ::mxx::comm const & comm) {
  BL_BENCH_INIT(map);

  std::vector<Kmer> query;

  BL_BENCH_START(map);
  // no transform involved.
  ::fsc::densehash_multimap<Kmer, Value, 
	::bliss::kmer::hash::sparsehash::special_keys<Kmer, canonical>,
	::bliss::transform::identity,
	::bliss::kmer::hash::farm<Kmer, false> > map(count);

  BL_BENCH_END(map, "reserve", count);


  {

//    BL_BENCH_START(map);
    std::vector<::std::pair<Kmer, Value> > input(count);
//    BL_BENCH_END(map, "reserve input", count);

//    BL_BENCH_START(map);
    generate_input(input, count, canonical);
    query.resize(count / query_frac);
    std::transform(input.begin(), input.begin() + input.size() / query_frac, query.begin(),
                   [](::std::pair<Kmer, Value> const & x){
      return x.first;
    });
//    BL_BENCH_END(map, "generate input", input.size());

    BL_BENCH_START(map);
    map.insert(input);
    BL_BENCH_END(map, "insert", map.size());
  }

  BL_BENCH_START(map);
  size_t result = 0;
  for (size_t i = 0, max = count / query_frac; i < max; ++i) {
    auto iters = map.equal_range(query[i]);
    for (auto it = iters.first; it != iters.second; ++it)
      result ^= (*it).second;
  }
  BL_BENCH_END(map, "find", result);


  BL_BENCH_START(map);
  result = 0;
  for (size_t i = 0, max = count / query_frac; i < max; ++i) {
    result += map.count(query[i]);
  }
  BL_BENCH_END(map, "count", result);



  BL_BENCH_START(map);
  result = map.erase(query.begin(), query.end());
  BL_BENCH_END(map, "erase", result);


  BL_BENCH_REPORT_MPI_NAMED(map, "densehash_full_multimap", comm);
}





template <typename DataType>
struct tommy_obj {
    tommy_node node;
    DataType * value;
};

template <typename Kmer, typename Value>
void benchmark_tommyhashdyn(size_t const count, size_t const query_frac, ::mxx::comm const & comm) {
  BL_BENCH_INIT(map);

  std::vector<Kmer> query;


  BL_BENCH_START(map);
  // init tommy hash
  tommy_hashdyn hashdyn;
  tommy_hashdyn_init(&hashdyn);
  // create array of tommy objects
  struct tommy_obj<std::pair<Kmer, Value> > *tommys = new tommy_obj<std::pair<Kmer, Value> >[count];
  // and a hasher
  ::bliss::kmer::hash::farm<Kmer, false> hasher;
  BL_BENCH_END(map, "init", count);


//  BL_BENCH_START(map);
  std::vector<::std::pair<Kmer, Value> > input(count);
//  BL_BENCH_END(map, "reserve input", count);

//  BL_BENCH_START(map);
  generate_input(input, count);
  query.resize(count / query_frac);
  std::transform(input.begin(), input.begin() + input.size() / query_frac, query.begin(),
                 [](::std::pair<Kmer, Value> const & x){
    return x.first;
  });
//  BL_BENCH_END(map, "generate input", count);


  BL_BENCH_START(map);
  // insert.
  uint key;
  for (size_t i = 0; i < count; ++i) {
    tommys[i].value = &(input[i]);
    key = hasher(input[i].first);
    tommy_hashdyn_insert(&hashdyn, &(tommys[i].node), &tommys[i], key);
  }
  BL_BENCH_END(map, "insert", tommy_hashdyn_count(&hashdyn));


  BL_BENCH_START(map);
  size_t result = 0;
  tommy_hashdyn_node * it;
  std::pair<Kmer, Value> * tmp;
  for (size_t i = 0, max = count / query_frac; i < max; ++i) {
    key = hasher(query[i]);
    it = tommy_hashdyn_bucket(&hashdyn, key);
    while (it) {
      tmp = static_cast<struct tommy_obj<std::pair<Kmer, Value> >*>(it->data)->value;  // data in a tommy_obj at the (tree) node.  pointer to a tuple of my data
      if (tmp->first == query[i]) result ^= tmp->second;  // bucket may have collided entries.
      it = it->next;
    }
  }
  BL_BENCH_END(map, "find", result);

  BL_BENCH_START(map);
  result = 0;
  for (size_t i = 0, max = count / query_frac; i < max; ++i) {
    key = hasher(query[i]);
    it = tommy_hashdyn_bucket(&hashdyn, key);
    while (it) {
      tmp = static_cast<struct tommy_obj<std::pair<Kmer, Value> >*>(it->data)->value;  // data in a tommy_obj at the (tree) node.  pointer to a tuple of my data
      if (tmp->first == query[i]) ++result;  // bucket may have collided entries.
      it = it->next;
    }
  }
  BL_BENCH_END(map, "count", result);

  BL_BENCH_START(map);
  result = 0;
//  std::vector<tommy_hashdyn_node *> todelete;
  tommy_hashdyn_node * next;
  for (size_t i = 0, max = count / query_frac; i < max; ++i) {
    key = hasher(query[i]);
    it = tommy_hashdyn_bucket(&hashdyn, key);
    while (it) {
      next = it->next;
      tmp = static_cast<struct tommy_obj<std::pair<Kmer, Value> >*>(it->data)->value;  // data in a tommy_obj at the (tree) node.  pointer to a tuple of my data
      if (tmp->first == query[i]) {
//        todelete.push_back(it);
        tommy_hashdyn_remove_existing(&hashdyn, it);
        ++result;
      }
      it = next;
    }
  }
//  for (size_t i = 0; i < todelete.size(); ++i) {
//    tommy_hashdyn_remove_existing(&hashdyn, todelete[i]);
//    ++result;  // bucket may have collided entries.
//  }
  BL_BENCH_END(map, "erase", result);

  BL_BENCH_START(map);
  tommy_hashdyn_done(&hashdyn);
  delete [] tommys;
  BL_BENCH_END(map, "cleanup", count);

  BL_BENCH_REPORT_MPI_NAMED(map, "tommyhashdyn", comm);
}


template <typename Kmer, typename Value>
void benchmark_tommyhashlin(size_t const count, size_t const query_frac, ::mxx::comm const & comm) {
  BL_BENCH_INIT(map);

  std::vector<Kmer> query;


  BL_BENCH_START(map);
  // init tommy hash
  tommy_hashlin hashlin;
  tommy_hashlin_init(&hashlin);
  // create array of tommy objects
  struct tommy_obj<std::pair<Kmer, Value> > *tommys = new tommy_obj<std::pair<Kmer, Value> >[count];
  // and a hasher
  ::bliss::kmer::hash::farm<Kmer, false> hasher;
  BL_BENCH_END(map, "init", count);


//  BL_BENCH_START(map);
  std::vector<::std::pair<Kmer, Value> > input(count);
//  BL_BENCH_END(map, "reserve input", count);

//  BL_BENCH_START(map);
  generate_input(input, count);
  query.resize(count / query_frac);
  std::transform(input.begin(), input.begin() + input.size() / query_frac, query.begin(),
                 [](::std::pair<Kmer, Value> const & x){
    return x.first;
  });
//  BL_BENCH_END(map, "generate input", count);


  BL_BENCH_START(map);
  // insert.
  uint key;
  for (size_t i = 0; i < count; ++i) {
    tommys[i].value = &(input[i]);
    key = hasher(input[i].first);
    tommy_hashlin_insert(&hashlin, &(tommys[i].node), &tommys[i], key);
  }
  BL_BENCH_END(map, "insert", tommy_hashlin_count(&hashlin));


  BL_BENCH_START(map);
  size_t result = 0;
  tommy_hashlin_node * it;
  std::pair<Kmer, Value> * tmp;
  for (size_t i = 0, max = count / query_frac; i < max; ++i) {
    key = hasher(query[i]);
    it = tommy_hashlin_bucket(&hashlin, key);
    while (it) {
      tmp = static_cast<struct tommy_obj<std::pair<Kmer, Value> >*>(it->data)->value;  // data in a tommy_obj at the (tree) node.  pointer to a tuple of my data
      if (tmp->first == query[i]) result ^= tmp->second;  // bucket may have collided entries.
      it = it->next;
    }
  }
  BL_BENCH_END(map, "find", result);

  BL_BENCH_START(map);
  result = 0;
  for (size_t i = 0, max = count / query_frac; i < max; ++i) {
    key = hasher(query[i]);
    it = tommy_hashlin_bucket(&hashlin, key);
    while (it) {
      tmp = static_cast<struct tommy_obj<std::pair<Kmer, Value> >*>(it->data)->value;  // data in a tommy_obj at the (tree) node.  pointer to a tuple of my data
      if (tmp->first == query[i]) ++result;  // bucket may have collided entries.
      it = it->next;
    }
  }
  BL_BENCH_END(map, "count", result);

  BL_BENCH_START(map);
  result = 0;
//  std::vector<tommy_hashlin_node *> todelete;
  tommy_hashlin_node * next;
  for (size_t i = 0, max = count / query_frac; i < max; ++i) {
    key = hasher(query[i]);
    it = tommy_hashlin_bucket(&hashlin, key);
    while (it) {
      next = it->next;
      tmp = static_cast<struct tommy_obj<std::pair<Kmer, Value> >*>(it->data)->value;  // data in a tommy_obj at the (tree) node.  pointer to a tuple of my data
      if (tmp->first == query[i]) {
        tommy_hashlin_remove_existing(&hashlin, it);
        ++result;
      }
      it = next;
    }
  }
//  for (size_t i = 0; i < todelete.size(); ++i) {
//    tommy_hashlin_remove_existing(&hashlin, todelete[i]);
//    ++result;  // bucket may have collided entries.
//  }
  BL_BENCH_END(map, "erase", result);

  BL_BENCH_START(map);
  tommy_hashlin_done(&hashlin);
  delete [] tommys;
  BL_BENCH_END(map, "cleanup", count);

  BL_BENCH_REPORT_MPI_NAMED(map, "tommyhashlin", comm);
}


template <typename Kmer, typename Value>
void benchmark_tommytrie(size_t const count, size_t const query_frac, ::mxx::comm const & comm) {
  BL_BENCH_INIT(map);

  std::vector<Kmer> query;


  BL_BENCH_START(map);
  // init tommy hash
  tommy_trie trie;
  tommy_allocator alloc;
  tommy_allocator_init(&alloc, TOMMY_TRIE_BLOCK_SIZE, TOMMY_TRIE_BLOCK_SIZE);

  tommy_trie_init(&trie, &alloc);
  // create array of tommy objects
  struct tommy_obj<std::pair<Kmer, Value> > *tommys = new tommy_obj<std::pair<Kmer, Value> >[count];
  // and a hasher
  BL_BENCH_END(map, "init", count);


//  BL_BENCH_START(map);
  std::vector<::std::pair<Kmer, Value> > input(count);
//  BL_BENCH_END(map, "reserve input", count);

//  BL_BENCH_START(map);
  generate_input(input, count);
  query.resize(count / query_frac);
  std::transform(input.begin(), input.begin() + input.size() / query_frac, query.begin(),
                 [](::std::pair<Kmer, Value> const & x){
    return x.first;
  });

//  BL_BENCH_END(map, "generate input", count);


  BL_BENCH_START(map);
  // insert.
  uint key;
  for (size_t i = 0; i < count; ++i) {
    tommys[i].value = &(input[i]);
    key = input[i].first.getDataRef()[Kmer::nWords - 1];
    tommy_trie_insert(&trie, &(tommys[i].node), &tommys[i], key);
  }
  BL_BENCH_END(map, "insert", tommy_trie_count(&trie));



  BL_BENCH_START(map);
  size_t result = 0;
  tommy_trie_node * it;
  std::pair<Kmer, Value> * tmp;
  for (size_t i = 0, max = count / query_frac; i < max; ++i) {
    key = query[i].getDataRef()[Kmer::nWords - 1];
    it = tommy_trie_bucket(&trie, key);
    while (it) {
      tmp = static_cast<struct tommy_obj<std::pair<Kmer, Value> >*>(it->data)->value;  // data in a tommy_obj at the (tree) node.  pointer to a tuple of my data
      if (tmp->first == query[i]) result ^= tmp->second;  // bucket may have collided entries.
      it = it->next;
    }
  }
  BL_BENCH_END(map, "find", result);

  BL_BENCH_START(map);
  result = 0;
  for (size_t i = 0, max = count / query_frac; i < max; ++i) {
    key = query[i].getDataRef()[Kmer::nWords - 1];
    it = tommy_trie_bucket(&trie, key);
    while (it) {
      tmp = static_cast<struct tommy_obj<std::pair<Kmer, Value> >*>(it->data)->value;  // data in a tommy_obj at the (tree) node.  pointer to a tuple of my data
      if (tmp->first == query[i]) ++result;  // bucket may have collided entries.
      it = it->next;
    }
  }
  BL_BENCH_END(map, "count", result);

  BL_BENCH_START(map);
  result = 0;
//  std::vector<tommy_trie_node *> todelete;
  tommy_trie_node * next;
  for (size_t i = 0, max = count / query_frac; i < max; ++i) {
    key = query[i].getDataRef()[Kmer::nWords - 1];
    it = tommy_trie_bucket(&trie, key);
    while (it) {
      next = it->next;
      tmp = static_cast<struct tommy_obj<std::pair<Kmer, Value> >*>(it->data)->value;  // data in a tommy_obj at the (tree) node.  pointer to a tuple of my data
      if (tmp->first == query[i]) {
        //todelete.push_back(it);
        tommy_trie_remove_existing(&trie, it);
        ++result;
      }
      it = next;
    }
  }
//  for (size_t i = 0; i < todelete.size(); ++i) {
//    tommy_trie_remove_existing(&trie, todelete[i]);
//    ++result;  // bucket may have collided entries.
//  }
  BL_BENCH_END(map, "erase", result);

  BL_BENCH_START(map);
  tommy_allocator_done(&alloc);
  delete [] tommys;
  BL_BENCH_END(map, "cleanup", count);

  BL_BENCH_REPORT_MPI_NAMED(map, "tommytrie", comm);
}


//========== has problems with casting to and from void*
// invalid conversion from ‘void*’ to ‘Word_t* {aka long unsigned int*}’ [-fpermissive]
//template <typename Kmer, typename V>
//void benchmark_judyl(size_t const count, size_t const query_frac, ::mxx::comm const & comm) {
//  BL_BENCH_INIT(map);
//
//
//  BL_BENCH_START(map);
//  std::vector<::std::pair<Kmer, V> > input(count);
//  BL_BENCH_END(map, "reserve input", count);
//
//  BL_BENCH_START(map);
//  generate_input(input, count);
//  BL_BENCH_END(map, "generate input", count);
//
//
//  BL_BENCH_START(map);
//
//  Word_t   Index;                     // array index
//  Word_t   Value;                     // array element value
//  Word_t  *PValue;                    // pointer to array element value
//  int      Rc_int;                    // return code
//  Word_t    Rc_word;
//  std::pair<Kmer, V> *tmp;
//  Word_t tmp_ptr;
//  Pvoid_t  PJLArray = (Pvoid_t) NULL; // initialize JudyL array
//
//  // and a hasher
//  ::bliss::kmer::hash::farm<Kmer, false> hasher;
//  BL_BENCH_END(map, "init", input.size());
//
//
//  BL_BENCH_START(map);
//  // insert.
//  Word_t key;
//  for (size_t i = 0; i < count; ++i) {
//    key = hasher(input[i].first);
//    JLI(PValue, PJLArray, key);
//    if (PValue == PJERR) continue;
//    tmp_ptr = reinterpret_cast<Word_t>(&input[i]);
//    *PValue = tmp_ptr;
//  }
//  Word_t count_total;
//  JLC(count_total, PJLArray, 0, -1);
//  BL_BENCH_END(map, "insert", count_total);
//
//
//  BL_BENCH_START(map);
//  size_t result = 0;
//  for (size_t i = 0, max = count / query_frac; i < max; ++i) {
//    key = hasher(input[i].first);
//    JLF(PValue, PJLArray, key);
//    tmp = reinterpret_cast<std::pair<Kmer, V> *>(*PValue);  // get pointer
//    if (tmp->first == input[i].first) result ^= tmp->second;  // bucket may have collided entries.
//  }
//  BL_BENCH_END(map, "find", result);
//
//  BL_BENCH_START(map);
//  result = 0;
//  for (size_t i = 0, max = count / query_frac; i < max; ++i) {
//    key = hasher(input[i].first);
//    JLF(PValue, PJLArray, key);
//    tmp = reinterpret_cast<std::pair<Kmer, V> *>(*PValue);  // get pointer
//    if (tmp->first == input[i].first) ++result;  // bucket may have collided entries.
//  }
//  BL_BENCH_END(map, "count", result);
//
//  BL_BENCH_START(map);
//  result = 0;
//  for (size_t i = 0, max = count / query_frac; i < max; ++i) {
//    key = hasher(input[i].first);
//    JLF(PValue, PJLArray, key);
//    tmp = reinterpret_cast<std::pair<Kmer, V> *>(*PValue);  // get pointer
//    if (tmp->first == input[i].first) {
//      JLD(Rc_int, PJLArray, key);
//      if (Rc_int == 1) ++result;  // bucket may have collided entries.
//    }
//  }
//  BL_BENCH_END(map, "erase", result);
//
//  BL_BENCH_START(map);
//  JLFA(Rc_word, PJLArray);
//  BL_BENCH_END(map, "cleanup", count);
//
//  BL_BENCH_REPORT_MPI_NAMED(map, "judyL", comm);
//}
//
//
//template <typename Kmer, typename Value>
//void benchmark_judyhs(size_t const count, size_t const query_frac, ::mxx::comm const & comm) {
//  BL_BENCH_INIT(map);
//
//
//  BL_BENCH_START(map);
//  std::vector<::std::pair<Kmer, size_t> > input(count);
//  BL_BENCH_END(map, "reserve input", count);
//
//  BL_BENCH_START(map);
//  generate_input(input, count);
//  BL_BENCH_END(map, "generate input", count);
//
//
//  BL_BENCH_START(map);
//
//  Word_t  * PValue;                           // JudyHS array element
//  int       Rc_int;                           // return flag
//  Word_t    Rc_word;                          // full word return value
//  Pvoid_t  PJLArray = (Pvoid_t) NULL; // initialize JudyL array
//  uint8_t * Index;                            // array-of-bytes pointer
//  Word_t    Length;                           // number of bytes in Index
//  std::pair<Kmer, Value> *tmp;
//  // and a hasher
//  ::bliss::kmer::hash::farm<Kmer, false> hasher;
//  BL_BENCH_END(map, "init", input.size());
//
//
//  BL_BENCH_START(map);
//  // insert.
//  Word_t key;
//  for (size_t i = 0; i < count; ++i) {
//    key = hasher(input[i].first);
//    JHSI(PValue, PJLArray, &key, sizeof(Word_t));
//    if (PValue == PJERR) continue;
//    if (*PValue == 0) continue; // duplicate
//    *PValue = reinterpret_cast<Word_t>(&input[i]);
//  }
//  BL_BENCH_END(map, "insert", count);
//
//
//  BL_BENCH_START(map);
//  size_t result = 0;
//  for (size_t i = 0, max = count / query_frac; i < max; ++i) {
//    key = hasher(input[i].first);
//    JHSG(PValue, PJLArray, &key, sizeof(Word_t));
//    tmp = *PValue;  // get pointer
//    if (tmp->first == input[i].first) result ^= tmp->second;  // bucket may have collided entries.
//  }
//  BL_BENCH_END(map, "find", result);
//
//  BL_BENCH_START(map);
//  result = 0;
//  for (size_t i = 0, max = count / query_frac; i < max; ++i) {
//    key = hasher(input[i].first);
//    JHSG(PValue, PJLArray, &key, sizeof(Word_t));
//    tmp = *PValue;  // get pointer
//    if (tmp->first == input[i].first) ++result;  // bucket may have collided entries.
//  }
//  BL_BENCH_END(map, "count", result);
//
//  BL_BENCH_START(map);
//  result = 0;
//  for (size_t i = 0, max = count / query_frac; i < max; ++i) {
//    key = hasher(input[i].first);
//    JHSD(Rc_int, PJLArray, &key, sizeof(Word_t));
//    if (Rc_int == 1) ++result;
//  }
//  BL_BENCH_END(map, "erase", result);
//
//  BL_BENCH_START(map);
//  JHSFA(Rc_word, PJLArray);
//  BL_BENCH_END(map, "cleanup", count);
//
//  BL_BENCH_REPORT_MPI_NAMED(map, "judyhs", comm);
//}

int main(int argc, char** argv) {

  size_t count = 10000000;
//  size_t count = 100;
  size_t query_frac = 10;


  mxx::env e(argc, argv);
  mxx::comm comm;

  if (comm.rank() == 0) printf("EXECUTING %s\n", argv[0]);

  comm.barrier();



  using Kmer = ::bliss::common::Kmer<31, ::bliss::common::DNA, uint64_t>;
  using DNA5Kmer = ::bliss::common::Kmer<21, ::bliss::common::DNA5, uint64_t>;
  using FullKmer = ::bliss::common::Kmer<32, ::bliss::common::DNA, uint64_t>;

  BL_BENCH_INIT(test);

  comm.barrier();

  BL_BENCH_START(test);
  benchmark_unordered_map<Kmer, size_t>(count, query_frac, comm);
  BL_BENCH_COLLECTIVE_END(test, "unordered_map", count, comm);

  BL_BENCH_START(test);
  benchmark_densehash_map<Kmer, size_t>(count, query_frac, comm);
  BL_BENCH_COLLECTIVE_END(test, "densehash_map_warmup", count, comm);

  BL_BENCH_START(test);
  benchmark_densehash_map<Kmer, size_t>(count, query_frac, comm);
  BL_BENCH_COLLECTIVE_END(test, "densehash_map", count, comm);

  BL_BENCH_START(test);
  benchmark_densehash_map<DNA5Kmer, size_t>(count, query_frac, comm);
  BL_BENCH_COLLECTIVE_END(test, "densehash_map_DNA5", count, comm);

  BL_BENCH_START(test);
  benchmark_densehash_full_map<FullKmer, size_t, true>(count, query_frac, comm);
  BL_BENCH_COLLECTIVE_END(test, "densehash_full_map_canonical", count, comm);

  BL_BENCH_START(test);
  benchmark_densehash_full_map<FullKmer, size_t, false>(count, query_frac, comm);
  BL_BENCH_COLLECTIVE_END(test, "densehash_full_map", count, comm);


  BL_BENCH_START(test);
  benchmark_unordered_multimap<Kmer, size_t>(count, query_frac, comm);
  BL_BENCH_COLLECTIVE_END(test, "unordered_multimap", count, comm);

//  BL_BENCH_START(test);
//  benchmark_densehash_vecmap<Kmer, size_t>(count, query_frac, comm);
//  BL_BENCH_COLLECTIVE_END(test, "densehash_vecmap", count, comm);

  BL_BENCH_START(test);
  benchmark_densehash_multimap<Kmer, size_t>(count, query_frac, comm);
  BL_BENCH_COLLECTIVE_END(test, "densehash_multimap_warmup", count, comm);

  BL_BENCH_START(test);
  benchmark_densehash_multimap<Kmer, size_t>(count, query_frac, comm);
  BL_BENCH_COLLECTIVE_END(test, "densehash_multimap", count, comm);

  BL_BENCH_START(test);
  benchmark_densehash_multimap<DNA5Kmer, size_t>(count, query_frac, comm);
  BL_BENCH_COLLECTIVE_END(test, "densehash_multimap_DNA5", count, comm);

  BL_BENCH_START(test);
  benchmark_densehash_full_multimap<FullKmer, size_t, false>(count, query_frac, comm);
  BL_BENCH_COLLECTIVE_END(test, "densehash_full_multimap", count, comm);

  BL_BENCH_START(test);
  benchmark_densehash_full_multimap<FullKmer, size_t, true>(count, query_frac, comm);
  BL_BENCH_COLLECTIVE_END(test, "densehash_full_multimap_canonical", count, comm);


  BL_BENCH_START(test);
  benchmark_tommyhashdyn<Kmer, size_t>(count, query_frac, comm);
  BL_BENCH_COLLECTIVE_END(test, "tommyhashdyn", count, comm);


  BL_BENCH_START(test);
  benchmark_tommyhashlin<Kmer, size_t>(count, query_frac, comm);
  BL_BENCH_COLLECTIVE_END(test, "tommyhashlin", count, comm);

  BL_BENCH_START(test);
  benchmark_tommytrie<Kmer, size_t>(count, query_frac, comm);
  BL_BENCH_COLLECTIVE_END(test, "tommytrie", count, comm);

  BL_BENCH_START(test);
  benchmark_unordered_vecmap<Kmer, size_t>(count, query_frac, comm);
  BL_BENCH_COLLECTIVE_END(test, "unordered_vecmap", count, comm);

//  BL_BENCH_START(test);
//  benchmark_hashed_vecmap<Kmer, size_t>(count, query_frac, comm);
//  BL_BENCH_COLLECTIVE_END(test, "hashed_vecmap", count, comm);


  BL_BENCH_REPORT_MPI_NAMED(test, "hashmaps", comm);

}

