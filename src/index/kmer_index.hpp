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
 * @file    kmer_index.hpp
 * @ingroup index
 * @author  Tony Pan <tpan7@gatech.edu>
 * @brief	High level kmer indexing API
 * @details	3 primary Kmer Index classes are currently provided:
 * 			Kmer Count Index,
 * 			Kmer Position Index, and
 * 			Kmer Position + Quality score index.
 *
 *
 *		Currently, KmerIndex classes only supports FASTQ files, but the intent is to
 *		  allow replaceable module for file reading.
 *
 *		Data distribution is via replaceable Hash functions.  currently,
 *		  PrefixHasher is used to distribute to MPI processes based on the first log_2(P) bits of the kmer
 *		  InfixHasher is used to distribute to threads within a process, based on the next log_2(p) bits of the kmer
 *		  SuffixHasher is used for threadlocal unordered_map hashing
 *
 *		All are deterministic to allow simple lookup processes.
 *
 */
#ifndef KMER_INDEX_HPP_
#define KMER_INDEX_HPP_

#include "bliss-config.hpp"

#if defined(USE_MPI)
#include "mpi.h"
#endif

//#if defined(USE_OPENMP)
//#include "omp.h"
//#endif

#include <tuple>        // tuple and utility functions
#include <utility>      // pair and utility functions.
#include <type_traits>
#include <cctype>       // tolower.

#include "io/file.hpp"
#include "io/fastq_loader.hpp"
#include "io/fasta_loader.hpp"

#include "utils/logging.h"
#include "common/kmer.hpp"
#include "common/base_types.hpp"
#include "common/sequence.hpp"
#include "utils/kmer_utils.hpp"
#include "index/kmer_hash.hpp"
#include "common/kmer_transform.hpp"

#include "io/kmer_file_helper.hpp"
#include "io/kmer_parser.hpp"
#include "io/mxx_support.hpp"
#include "containers/distributed_unordered_map.hpp"
#include "containers/distributed_sorted_map.hpp"
//#include "containers/distributed_hashed_vec.hpp"
//#include "containers/distributed_map.hpp"
#include "containers/distributed_densehash_map.hpp"

#include "utils/benchmark_utils.hpp"
#include "utils/file_utils.hpp"
#include "utils/transform_utils.hpp"


#include <fstream> // debug only
#include <iostream>  // debug only
#include <sstream>  // debug only

namespace bliss
{
namespace index
{
namespace kmer
{

/**
 * @tparam MapType  	container type
 * @tparam KmerParser		functor to generate kmer (tuple) from input.  specified here so we specialize for different index.  note KmerParser needs to be supplied with a data type.
 */
template <typename MapType, typename KmerParser>
class Index {
protected:

	MapType map;

	const mxx::comm& comm;

public:
	using KmerType = typename MapType::key_type;
	// TODO: make this consistent with map data type conventions?
	using ValueType   = typename MapType::mapped_type;
	using TupleType =  std::pair<KmerType, ValueType>;
	using Alphabet = typename KmerType::KmerAlphabet;

	using KmerParserType = KmerParser;

	Index(const mxx::comm& _comm) : map(_comm), comm(_comm) {
	}

	virtual ~Index() {};

	MapType & get_map() {
		return map;
	}
	MapType const & get_map() const {
		return map;
	}



//	std::vector<TupleType> find_overlap(std::vector<KmerType> &query) const {
//		return map.find_overlap(query);
//	}
	auto find(std::vector<KmerType> &query) const
		-> decltype(::std::declval<MapType>().find(::std::declval<std::vector<KmerType> &>())) {
		return map.find(query);
	}
//	std::vector<TupleType> find_collective(std::vector<KmerType> &query) const {
//		return map.find_collective(query);
//	}
//  std::vector<TupleType> find_sendrecv(std::vector<KmerType> &query) const {
//    return map.find_sendrecv(query);
//  }
	auto count(std::vector<KmerType> &query) const
	-> decltype(::std::declval<MapType>().count(::std::declval<std::vector<KmerType> &>())){
		return map.count(query);
	}

	void erase(std::vector<KmerType> &query) {
		map.erase(query);
	}


//	template <typename Predicate>
//	std::vector<TupleType> find_if_overlap(std::vector<KmerType> &query, Predicate const &pred) const {
//		return map.find_overlap(query, false, pred);
//	}
	template <typename Predicate>
	auto find_if(std::vector<KmerType> &query, Predicate const &pred) const
	-> decltype(::std::declval<MapType>().find(::std::declval<std::vector<KmerType> &>())) {
		return map.find(query, false, pred);
	}
//	template <typename Predicate>
//	std::vector<TupleType> find_if_collective(std::vector<KmerType> &query, Predicate const &pred) const {
//		return map.find_collective(query, false, pred);
//	}
//  template <typename Predicate>
//  std::vector<TupleType> find_if_sendrecv(std::vector<KmerType> &query, Predicate const &pred) const {
//    return map.find_sendrecv(query, false, pred);
//  }
	template <typename Predicate>
	std::vector<TupleType> find_if(Predicate const &pred) const {
		return map.find(pred);
	}

	template <typename Predicate>
	auto count_if(std::vector<KmerType> &query, Predicate const &pred) const
	-> decltype(::std::declval<MapType>().find(::std::declval<std::vector<KmerType> &>())) {
		return map.count(query, false, pred);
	}

	template <typename Predicate>
	std::vector<std::pair<KmerType, size_t> > count_if(Predicate const &pred) const {
		return map.count(pred);
	}


	template <typename Predicate>
	void erase_if(std::vector<KmerType> &query, Predicate const &pred) {
		map.erase(query, false, pred);
	}

	template <typename Predicate>
	void erase_if(Predicate const &pred) {
		map.erase(pred);
	}


	/**
	 * @tparam T 	input type may not be same as map's value types, so map need to provide overloads (and potentially with transform operators)
	 */
	 template <typename T>
	void insert(std::vector<T> &temp) {
		BL_BENCH_INIT(insert);

		// do not reserve until insertion - less transient memory used.
//		BL_BENCH_START(build);
//		this->map.reserve(this->map.size() + temp.size());
//		BL_BENCH_END(build, "reserve", temp.size());

		// distribute
		BL_BENCH_START(insert);
		this->map.insert(temp);  // COLLECTIVE CALL...
		BL_BENCH_END(insert, "map_insert", this->map.local_size());

#if (BL_BENCHMARK == 1)
		BL_BENCH_START(insert);
		size_t m = 0;  // here because sortmap needs it.
		m = this->map.get_multiplicity();
		BL_BENCH_END(insert, "multiplicity", m);
#else
		auto result = this->map.get_multiplicity();
		BLISS_UNUSED(result);
#endif
		BL_BENCH_REPORT_MPI_NAMED(insert, "index:insert", this->comm);

	 }

	 // Note that KmerParserType may depend on knowing the Sequence Parser Type (e.g. provide quality score iterators)
	 //	Output type of KmerParserType may not match Map value type, in which case the map needs to do its own transform.
	 //     since Kmer template parameter is not explicitly known, we can't hard code the return types of KmerParserType.

	 //============= THESE ARE TO BE DEPRECATED


	 // Note that KmerParserType may depend on knowing the Sequence Parser Type (e.g. provide quality score iterators)
	 //	Output type of KmerParserType may not match Map value type, in which case the map needs to do its own transform.
	 //     since Kmer template parameter is not explicitly known, we can't hard code the return types of KmerParserType.

	 /// convenience function for building index.
	 template <template <typename> class SeqParser, template <typename, template <typename> class> class SeqIterType>
	 void build_mpiio(const std::string & filename, MPI_Comm comm) {

		 // file extension determines SeqParserType
		 std::string extension = ::bliss::utils::file::get_file_extension(filename);
		 std::transform(extension.begin(), extension.end(), extension.begin(), ::tolower);
		 if ((extension.compare("fastq") != 0) && (extension.compare("fasta") != 0) && (extension.compare("fa") != 0)) {
			 throw std::invalid_argument("input filename extension is not supported.");
		 }

		 // check to make sure that the file parser will work
		 if ((extension.compare("fastq") == 0) && (!std::is_same<SeqParser<char*>, ::bliss::io::FASTQParser<char*> >::value)) {
			 throw std::invalid_argument("Specified File Parser template parameter does not support files with fastq extension.");
		 } else if (((extension.compare("fasta") == 0) || (extension.compare("fa") == 0)) && (!std::is_same<SeqParser<char*>, ::bliss::io::FASTAParser<char*> >::value)) {
			 throw std::invalid_argument("Specified File Parser template parameter does not support files with fasta extension.");
		 }
     BL_BENCH_INIT(build);

		 // proceed
     BL_BENCH_START(build);
		 ::std::vector<typename KmerParser::value_type> temp;
		 bliss::io::KmerFileHelper::template read_file_mpiio<KmerParser, SeqParser, SeqIterType>(filename, temp, comm);
     BL_BENCH_END(build, "read", temp.size());


		 //        // dump the generated kmers to see if they look okay.
		 //         std::stringstream ss;
		 //         ss << "test." << commRank << ".log";
		 //         std::ofstream ofs(ss.str());
		 //         for (int i = 0; i < temp.size(); ++i) {
		 //          ofs << "item: " << i << " value: " << temp[i] << std::endl;
		 //         }
		 //         ofs.close();

     BL_BENCH_START(build);
		 this->insert(temp);
     BL_BENCH_END(build, "insert", temp.size());


     BL_BENCH_REPORT_MPI_NAMED(build, "index:build_mpiio", this->comm);

	 }


	  /// convenience function for building index.
	   template <template <typename> class SeqParser, template <typename,  template <typename> class> class SeqIterType>
	   void build_mmap(const std::string & filename, MPI_Comm comm) {

	     // file extension determines SeqParserType
	     std::string extension = ::bliss::utils::file::get_file_extension(filename);
	     std::transform(extension.begin(), extension.end(), extension.begin(), ::tolower);
	     if ((extension.compare("fastq") != 0) && (extension.compare("fasta") != 0) && (extension.compare("fa") != 0)) {
	       throw std::invalid_argument("input filename extension is not supported.");
	     }

	     // check to make sure that the file parser will work
	     if ((extension.compare("fastq") == 0) && (!std::is_same<SeqParser<char*>, ::bliss::io::FASTQParser<char*> >::value)) {
	       throw std::invalid_argument("Specified File Parser template parameter does not support files with fastq extension.");
	     } else if (((extension.compare("fasta") == 0) || (extension.compare("fa") == 0)) && (!std::is_same<SeqParser<char*>, ::bliss::io::FASTAParser<char*> >::value)) {
	       throw std::invalid_argument("Specified File Parser template parameter does not support files with fasta extension.");
	     }
	     BL_BENCH_INIT(build);

	     // proceed
	     BL_BENCH_START(build);
	     ::std::vector<typename KmerParser::value_type> temp;
	     bliss::io::KmerFileHelper::template read_file_mmap<KmerParser, SeqParser, SeqIterType>(filename, temp, comm);
	      BL_BENCH_END(build, "read", temp.size());


	     //        // dump the generated kmers to see if they look okay.
	     //         std::stringstream ss;
	     //         ss << "test." << commRank << ".log";
	     //         std::ofstream ofs(ss.str());
	     //         for (int i = 0; i < temp.size(); ++i) {
	     //          ofs << "item: " << i << " value: " << temp[i] << std::endl;
	     //         }
	     //         ofs.close();

	     BL_BENCH_START(build);
	     this->insert(temp);
	      BL_BENCH_END(build, "insert", temp.size());


	     BL_BENCH_REPORT_MPI_NAMED(build, "index:build_mmap", this->comm);

	   }




		 /// convenience function for building index.
		 template <template <typename> class SeqParser, template <typename,  template <typename> class> class SeqIterType>
		 void build_posix(const std::string & filename, MPI_Comm comm) {

			 // file extension determines SeqParserType
			 std::string extension = ::bliss::utils::file::get_file_extension(filename);
			 std::transform(extension.begin(), extension.end(), extension.begin(), ::tolower);
			 if ((extension.compare("fastq") != 0) && (extension.compare("fasta") != 0) && (extension.compare("fa") != 0)) {
				 throw std::invalid_argument("input filename extension is not supported.");
			 }

			 // check to make sure that the file parser will work
			 if ((extension.compare("fastq") == 0) && (!std::is_same<SeqParser<char*>, ::bliss::io::FASTQParser<char*> >::value)) {
				 throw std::invalid_argument("Specified File Parser template parameter does not support files with fastq extension.");
			 } else if (((extension.compare("fasta") == 0) || (extension.compare("fa") == 0)) && (!std::is_same<SeqParser<char*>, ::bliss::io::FASTAParser<char*> >::value)) {
				 throw std::invalid_argument("Specified File Parser template parameter does not support files with fasta extension.");
			 }
	     BL_BENCH_INIT(build);

			 // proceed
	     BL_BENCH_START(build);
			 ::std::vector<typename KmerParser::value_type> temp;
			 bliss::io::KmerFileHelper::template read_file_posix<KmerParser, SeqParser, SeqIterType>(filename, temp, comm);
	     BL_BENCH_END(build, "read", temp.size());


			 //        // dump the generated kmers to see if they look okay.
			 //         std::stringstream ss;
			 //         ss << "test." << commRank << ".log";
			 //         std::ofstream ofs(ss.str());
			 //         for (int i = 0; i < temp.size(); ++i) {
			 //          ofs << "item: " << i << " value: " << temp[i] << std::endl;
			 //         }
			 //         ofs.close();

	     BL_BENCH_START(build);
			 this->insert(temp);
	     BL_BENCH_END(build, "insert", temp.size());


	     BL_BENCH_REPORT_MPI_NAMED(build, "index:build_posix", this->comm);

		 }




   typename MapType::const_iterator cbegin() const
   {
     return map.cbegin();
   }

   typename MapType::const_iterator cend() const {
     return map.cend();
   }

   size_t size() const {
     return map.size();
   }

   size_t local_size() const {
     return map.local_size();
   }

};



// TODO: the types of Map that is used should be restricted.  (perhaps via map traits)
template <typename MapType>
using KmerIndex = Index<MapType, KmerParser<typename MapType::key_type> >;

template <typename MapType>
using PositionIndex = Index<MapType, KmerPositionTupleParser<std::pair<typename MapType::key_type, typename MapType::mapped_type> > >;

template <typename MapType>
using PositionQualityIndex = Index<MapType, KmerPositionQualityTupleParser<std::pair<typename MapType::key_type, typename MapType::mapped_type> > >;

template <typename MapType>
using CountIndex = Index<MapType, KmerCountTupleParser<std::pair<typename MapType::key_type, typename MapType::mapped_type> > >;
template <typename MapType>
using CountIndex2 = Index<MapType, KmerParser<typename MapType::key_type> >;

// template aliases for hash to be used as distribution hash
template <typename Key>
using DistHashFarm = ::bliss::kmer::hash::farm<Key, true>;
template <typename Key>
using DistHashMurmur = ::bliss::kmer::hash::murmur<Key, true>;
template <typename Key>
using DistHashStd = ::bliss::kmer::hash::cpp_std<Key, true>;
template <typename Key>
using DistHashIdentity = ::bliss::kmer::hash::identity<Key, true>;


template <typename Key>
using StoreHashFarm = ::bliss::kmer::hash::farm<Key, false>;
template <typename Key>
using StoreHashMurmur = ::bliss::kmer::hash::murmur<Key, false>;
template <typename Key>
using StoreHashStd = ::bliss::kmer::hash::cpp_std<Key, false>;
template <typename Key>
using StoreHashIdentity = ::bliss::kmer::hash::identity<Key, false>;

// =================  Partially defined aliases for MapParams, for distributed_xxx_maps.
// NOTE: when using this, need to further alias so that only Key param remains.
// =================
template <typename Key,
			template <typename> class DistHash  = DistHashMurmur,
			template <typename> class StoreHash = StoreHashMurmur,
			template <typename> class DistTrans = ::bliss::transform::identity
			>
using SingleStrandHashMapParams = ::dsc::HashMapParams<
		Key,
		::bliss::transform::identity,  // precanonalizer
		 DistTrans,  				// could be iden, xor, lex_less
		  DistHash,
		  ::std::equal_to,
		   ::bliss::transform::identity,  // only one that makes sense given InputTransform
		    StoreHash,
		    ::std::equal_to
		  >;


template <typename Key,
	template <typename> class DistHash  = DistHashMurmur,
	template <typename> class StoreHash = StoreHashMurmur
>
using CanonicalHashMapParams = ::dsc::HashMapParams<
		Key,
		::bliss::kmer::transform::lex_less,  // precanonalizer
		 ::bliss::transform::identity,  // only one that makes sense given InputTransform
		  DistHash,
		  ::std::equal_to,
		   ::bliss::transform::identity,
		    StoreHash,
		    ::std::equal_to
		  >;

template <typename Key,
	template <typename> class DistHash  = DistHashMurmur,
	template <typename> class StoreHash = StoreHashMurmur
	>
using BimoleculeHashMapParams = ::dsc::HashMapParams<
		Key,
		::bliss::transform::identity,  // precanonalizer - only one that makes sense for bimole
		 ::bliss::kmer::transform::lex_less,  // only one that makes sense for bimole
		  DistHash,
		  ::std::equal_to,
		   ::bliss::kmer::transform::lex_less,
		    StoreHash,
		    ::std::equal_to
		  >;

//template <typename Key,
//			template <typename> class DistHash = DistHashMurmur,
//			template <typename> class StoreLess = ::std::less,
//			template <typename> class DistTrans = ::bliss::transform::identity >
//using SingleStrandOrderedMapParams = ::dsc::OrderedMapParams<
//		Key,
//		::bliss::transform::identity,  // precanonalizer
//		 DistTrans,  				// could be iden, xor, lex_less
//		  DistHash,
//		  ::std::equal_to,
//		   ::bliss::transform::identity,  // only one that makes sense given InputTransform
//		    StoreLess,
//		    ::std::equal_to
//		  >;
//
//
//template <typename Key,
//			template <typename> class DistHash = DistHashMurmur,
//			template <typename> class StoreLess = ::std::less
//			>
//using CanonicalOrderedMapParams = ::dsc::OrderedMapParams<
//		Key,
//		::bliss::kmer::transform::lex_less,  // precanonalizer
//		 ::bliss::transform::identity,  // only one that makes sense given InputTransform
//		  DistHash,
//		  ::std::equal_to,
//		   ::bliss::transform::identity,
//		    StoreLess,
//		    ::std::equal_to
//		  >;
//
//template <typename Key,
//			template <typename> class DistHash = DistHashMurmur,
//			template <typename> class StoreLess = ::std::less
//			>
//using BimoleculeOrderedMapParams = ::dsc::OrderedMapParams<
//		Key,
//		::bliss::transform::identity,  // precanonalizer - only one that makes sense for bimole
//		 ::bliss::kmer::transform::lex_less,  // only one that makes sense for bimole
//		  DistHash,
//		  ::std::equal_to,
//		   ::bliss::kmer::transform::lex_less,
//		    StoreLess,
//		    ::std::equal_to
//		  >;


template <typename Key,
  template <typename> class Less = ::std::less
>
using SingleStrandSortedMapParams = ::dsc::SortedMapParams<
		Key,
		::bliss::transform::identity,  // precanonalizer
		   ::bliss::transform::identity,  // only one that makes sense given InputTransform
		    Less,
		    ::std::equal_to
		  >;


template <typename Key,
			template <typename> class Less = ::std::less
>
using CanonicalSortedMapParams = ::dsc::SortedMapParams<
		Key,
		::bliss::kmer::transform::lex_less,  // precanonalizer
		 ::bliss::transform::identity,  // only one that makes sense given InputTransform
		    Less,
		    ::std::equal_to
		  >;

template <typename Key,
			template <typename> class Less = ::std::less
>
using BimoleculeSortedMapParams = ::dsc::SortedMapParams<
		Key,
		::bliss::transform::identity,  // precanonalizer - only one that makes sense for bimole
		 ::bliss::kmer::transform::lex_less,  // only one that makes sense for bimole
		    Less,
		    ::std::equal_to
		  >;
} /* namespace kmer */

} /* namespace index */
} /* namespace bliss */

#endif /* KMER_INDEX_HPP_ */
