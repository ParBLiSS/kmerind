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

#include <cstdio>
#include <stdint.h>

#include <tuple>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>


#include "common/kmer.hpp"
#include "index/kmer_hash.hpp"
#include "common/kmer_transform.hpp"

#include "boost/pool/pool_alloc.hpp"

#include "utils/benchmark_utils.hpp"


using Kmer = bliss::common::Kmer<21, bliss::common::DNA, uint64_t>;
using KmerPos = ::std::pair<Kmer, uint64_t>;
using KmerPosQual = ::std::pair<Kmer, ::std::pair<uint64_t, float> >;

typedef bliss::kmer::transform::lex_less<Kmer> KeyTransform;
typedef bliss::kmer::hash::farm<Kmer, false> KmerHash;

struct TransformedHash {
    KmerHash h;
    KeyTransform trans;

    inline uint64_t operator()(Kmer const& k) const {
      return h(trans(k));
    }
    template<typename V>
    inline uint64_t operator()(::std::pair<Kmer, V> const& x) const {
      return this->operator()(x.first);
    }
    template<typename V>
    inline uint64_t operator()(::std::pair<const Kmer, V> const& x) const {
      return this->operator()(x.first);
    }
};

template <typename Comparator>
struct TransformedComp {
    Comparator comp;
    KeyTransform trans;

    inline bool operator()(Kmer const & x, Kmer const & y) const {
      return comp(trans(x), trans(y));
    }
    template<typename V>
    inline bool operator()(::std::pair<Kmer, V> const & x, Kmer const & y) const {
      return this->operator()(x.first, y);
    }
    template<typename V>
    inline bool operator()(::std::pair<const Kmer, V> const & x, Kmer const & y) const {
      return this->operator()(x.first, y);
    }
    template<typename V>
    inline bool operator()(Kmer const & x, ::std::pair<Kmer, V> const & y) const {
      return this->operator()(x, y.first);
    }
    template<typename V>
    inline bool operator()(Kmer const & x, ::std::pair<const Kmer, V> const & y) const {
      return this->operator()(x, y.first);
    }
    template<typename V>
    inline bool operator()(::std::pair<Kmer, V> const & x, ::std::pair<Kmer, V> const & y) const {
      return this->operator()(x.first, y.first);
    }
    template<typename V>
    inline bool operator()(::std::pair<const Kmer, V> const & x, ::std::pair<const Kmer, V> const & y) const {
      return this->operator()(x.first, y.first);
    }
};

typedef TransformedComp<::std::less<Kmer> > TransformedLess;
typedef TransformedComp<::std::equal_to<Kmer> > TransformedEqual;

template <typename T>
using vector_stl_alloc = ::std::vector<T, ::std::allocator<T > >;

template <typename T>
using vector_boost_pool = ::std::vector<T, ::boost::pool_allocator<T > >;

template <typename K, typename T>
using umap_stl_alloc = ::std::unordered_map<K, T, TransformedHash, TransformedEqual, ::std::allocator<::std::pair<const K, T> > >;

template <typename K, typename T>
using umap_boost_pool = ::std::unordered_map<K, T, TransformedHash, TransformedEqual, ::boost::fast_pool_allocator<::std::pair<const K, T> > >;




class KmerHelper
{
public:

    std::default_random_engine generator;
    std::uniform_int_distribution<uint64_t> distribution;
    uint64_t v[1];

    Kmer random_kmer() {
        *v = distribution(generator);
        return Kmer(v);
    }

};




int main(int argc, char** argv) {

	int iterations = 1000000;
	if (argc > 1) {
		iterations = atoi(argv[1]);
	}


	printf("SIZES:\n");
	printf("kmer size: %lu\n", sizeof(Kmer));
	printf("Kmer Pos size: %lu\n", sizeof(KmerPos));
	printf("Kmer PosQual size: %lu\n", sizeof(KmerPosQual));

	printf("Vector of Kmer, base size: %lu\n", sizeof(::std::vector<Kmer>));
	printf("Vector of KmerPos, base size: %lu\n", sizeof(::std::vector<KmerPos >));
	printf("Vector of KmerPosQual, base size: %lu\n", sizeof(::std::vector<KmerPosQual >));

	printf("Kmer Vector1 pair, base size: %lu\n", sizeof(::std::pair<Kmer, ::std::vector<Kmer> >));
	printf("Kmer Vector2 pair, base size: %lu\n", sizeof(::std::pair<Kmer, ::std::vector<KmerPos> >));
	printf("Kmer Vector3 pair, base size: %lu\n", sizeof(::std::pair<Kmer, ::std::vector<KmerPosQual > >));



	KmerHelper helper;

	BL_BENCH_INIT(vector);

	Kmer result;

	BL_BENCH_START(vector);
	{
		vector_stl_alloc<KmerPos> stlvec;
		for (int i = 0; i < iterations; ++i) {
			stlvec.emplace_back(helper.random_kmer(), 1UL);
		}
		result = helper.random_kmer();
		for (int i = 0; i < iterations; ++i) {
			result ^= stlvec[i].first.reverse_complement();
		}
	}
	BL_BENCH_END(vector, "stl", iterations);
	printf("result : %s\n", result.toAlphabetString().c_str());


	BL_BENCH_START(vector);
	{
		vector_boost_pool<KmerPos> boostvec;
		for (int i = 0; i < iterations; ++i) {
			boostvec.emplace_back(helper.random_kmer(), 1UL);
		}
		result = helper.random_kmer();
		for (int i = 0; i < iterations; ++i) {
			result ^= boostvec[i].first.reverse_complement();
		}
		boost::singleton_pool<boost::pool_allocator_tag, sizeof(KmerPos)>::release_memory();

		// could be 4x slower than stl allocator.
	}
	BL_BENCH_END(vector, "boost", iterations);
	printf("result : %s\n", result.toAlphabetString().c_str());


	BL_BENCH_REPORT(vector);



	BL_BENCH_INIT(map);

	BL_BENCH_START(map);
	{
		umap_stl_alloc<Kmer, uint64_t> stlumap;
		for (int i = 0; i < iterations; ++i) {
			stlumap.emplace(helper.random_kmer(), 1UL);
		}
		result = helper.random_kmer();
		for (auto it = stlumap.begin(), max = stlumap.end(); it != max; ++it) {
			result ^= it->first.reverse_complement();
		}
	}
	BL_BENCH_END(map, "stl", iterations);
	printf("result : %s\n", result.toAlphabetString().c_str());


	// cannot get this to work...
//	BL_BENCH_START(map);
//	{
//		umap_boost_pool<Kmer, uint64_t> boostumap;
//		for (int i = 0; i < iterations; ++i) {
//			boostumap.emplace(helper.random_kmer(), 1UL);
//		}
//		result = helper.random_kmer();
//		for (auto it = boostumap.begin(), max = boostumap.end(); it != max; ++it) {
//			result ^= it->first.reverse_complement();
//		}
//
//		boost::singleton_pool<boost::pool_allocator_tag, sizeof(KmerPos)>::release_memory();
//
//		// could be 4x slower than stl allocator.
//	}
//	BL_BENCH_END(map, "boost", iterations);
//	printf("result : %s\n", result.toAlphabetString().c_str());


	BL_BENCH_REPORT(map);


}
