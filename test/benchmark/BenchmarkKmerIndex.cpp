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
 * @file    test_threads.cpp
 * @ingroup
 * @author  Tony Pan <tpan7@gatech.edu>
 * @brief
 * @details
 *

 */

#include "bliss-config.hpp"



#include <functional>
#include <random>
#include <algorithm>
#include <string>
#include <sstream>
#include <chrono>
#include <iostream>  // for system("pause");

#include "utils/logging.h"

#include "common/alphabets.hpp"
#include "common/kmer.hpp"
#include "common/base_types.hpp"
#include "utils/kmer_utils.hpp"

#include "io/mxx_support.hpp"

#include "io/sequence_iterator.hpp"
#include "io/sequence_id_iterator.hpp"

#include "iterators/transform_iterator.hpp"

#include "common/kmer_iterators.hpp"

#include "iterators/zip_iterator.hpp"

#include "index/quality_score_iterator.hpp"

#include "index/kmer_index.hpp"

#include "utils/benchmark_utils.hpp"
#include "utils/exception_handling.hpp"

#include "tclap/CmdLine.h"

#include "mxx/env.hpp"
#include "mxx/comm.hpp"

// ================ define preproc macro constants
// needed as #if can only calculate constant int expressions
#define FASTA 1
#define FASTQ 0

#define IDEN 10
#define LEX 11
#define XOR 12

#define STD 21
#define MURMUR 22
#define FARM 23

#define POS 31
#define POSQUAL 32
#define COUNT 33

#define SORTED 41
#define ORDERED 42
#define VEC 43
#define COMPACTVEC 44
#define HASHEDVEC 45
#define UNORDERED 46
#define DENSEHASH 47

#define SINGLE 51
#define CANONICAL 52
#define BIMOLECULE 53



//================= define types - changeable here...

//========   Kmer parameters
#if (pDNA == 16)
using Alphabet = bliss::common::DNA16;
#elif (pDNA == 5)
using Alphabet = bliss::common::DNA5;
#elif (pDNA == 4)
using Alphabet = bliss::common::DNA;
#endif

#if defined(pK)
using KmerType = bliss::common::Kmer<pK, Alphabet, WordType>;
#else
using KmerType = bliss::common::Kmer<21, Alphabet, WordType>;
#endif

//============== index input file format
#if (pPARSER == FASTA)
	using IdType = bliss::common::LongSequenceKmerId;
#define PARSER_TYPE ::bliss::io::FASTAParser
#elif (pPARSER == FASTQ)
	using IdType = bliss::common::ShortSequenceKmerId;
#define PARSER_TYPE ::bliss::io::FASTQParser
#endif



// ============  index value type
using QualType = float;
using KmerInfoType = std::pair<IdType, QualType>;
using CountType = uint32_t;

#if (pINDEX == POS)
	using ValType = IdType;
#elif (pINDEX == POSQUAL)
	using ValType = KmerInfoType;
#elif (pINDEX == COUNT)
	using ValType = CountType;
#endif




//============== MAP properties


//----- get them all. may not use subsequently.

// distribution transforms
#if (pDistTrans == LEX)
	template <typename KM>
	using DistTrans = bliss::kmer::transform::lex_less<KM>;
#elif (pDistTrans == XOR)
	template <typename KM>
	using DistTrans = bliss::kmer::transform::xor_rev_comp<KM>;
#else //if (pDistTrans == IDEN)
	template <typename KM>
	using DistTrans = bliss::kmer::transform::identity<KM>;
#endif

// distribution hash
#if (pDistHash == STD)
	template <typename KM>
	using DistHash = bliss::kmer::hash::cpp_std<KM, true>;
#elif (pDistHash == IDEN)
	template <typename KM>
	using DistHash = bliss::kmer::hash::identity<KM, true>;
#elif (pDistHash == MURMUR)
	template <typename KM>
	using DistHash = bliss::kmer::hash::murmur<KM, true>;
#else // if (pDistHash == FARM)
	template <typename KM>
	using DistHash = bliss::kmer::hash::farm<KM, true>;
#endif


// storage hash type
#if (pStoreHash == STD)
	template <typename KM>
	using StoreHash = bliss::kmer::hash::cpp_std<KM, false>;
#elif (pStoreHash == IDEN)
	template <typename KM>
	using StoreHash = bliss::kmer::hash::identity<KM, false>;
#elif (pStoreHash == MURMUR)
	template <typename KM>
	using StoreHash = bliss::kmer::hash::murmur<KM, false>;
#else //if (pStoreHash == FARM)
	template <typename KM>
	using StoreHash = bliss::kmer::hash::farm<KM, false>;
#endif



// ==== define Map parameter
#if (pMAP == SORTED)
	// choose a MapParam based on type of map and kmer model (canonical, original, bimolecule)
	#if (pKmerStore == SINGLE)  // single stranded
		template <typename Key>
		using MapParams = ::bliss::index::kmer::SingleStrandSortedMapParams<Key>;
	#elif (pKmerStore == CANONICAL)
		template <typename Key>
		using MapParams = ::bliss::index::kmer::CanonicalSortedMapParams<Key>;
	#elif (pKmerStore == BIMOLECULE)  // bimolecule
		template <typename Key>
		using MapParams = ::bliss::index::kmer::BimoleculeSortedMapParams<Key>;
	#endif

	// DEFINE THE MAP TYPE base on the type of data to be stored.
	#if (pINDEX == POS) || (pINDEX == POSQUAL)  // multimap
		using MapType = ::dsc::sorted_multimap<
				KmerType, ValType, MapParams>;
	#elif (pINDEX == COUNT)  // map
		using MapType = ::dsc::counting_sorted_map<
				KmerType, ValType, MapParams>;
	#endif


#elif (pMAP == ORDERED)
	// choose a MapParam based on type of map and kmer model (canonical, original, bimolecule)
	#if (pKmerStore == SINGLE)  // single stranded
		template <typename Key>
		using MapParams = ::bliss::index::kmer::SingleStrandOrderedMapParams<Key, DistHash, ::std::less, DistTrans>;
	#elif (pKmerStore == CANONICAL)
		template <typename Key>
		using MapParams = ::bliss::index::kmer::CanonicalOrderedMapParams<Key, DistHash>;
	#elif (pKmerStore == BIMOLECULE)  // bimolecule
		template <typename Key>
		using MapParams = ::bliss::index::kmer::BimoleculeOrderedMapParams<Key, DistHash>;
	#endif

	// DEFINE THE MAP TYPE base on the type of data to be stored.
	#if (pINDEX == POS) || (pINDEX == POSQUAL)  // multimap
		using MapType = ::dsc::multimap<
				KmerType, ValType, MapParams>;
	#elif (pINDEX == COUNT)  // map
		using MapType = ::dsc::counting_map<
				KmerType, ValType, MapParams>;
	#endif


#else  // hashmap
		struct MSBSplitter {
		    template <typename Kmer>
		    bool operator()(Kmer const & kmer) const {
		      return (kmer.getData()[Kmer::nWords - 1] & ~(~(static_cast<typename Kmer::KmerWordType>(0)) >> 2)) > 0;
		    }

		    template <typename Kmer, typename V>
		    bool operator()(std::pair<Kmer, V> const & x) const {
		      return (x.first.getData()[Kmer::nWords - 1] & ~(~(static_cast<typename Kmer::KmerWordType>(0)) >> 2)) > 0;
		    }
		};


  // choose a MapParam based on type of map and kmer model (canonical, original, bimolecule)
  #if (pKmerStore == SINGLE)  // single stranded
    template <typename Key>
    using MapParams = ::bliss::index::kmer::SingleStrandHashMapParams<Key, DistHash, StoreHash, DistTrans>;
  #elif (pKmerStore == CANONICAL)
    template <typename Key>
    using MapParams = ::bliss::index::kmer::CanonicalHashMapParams<Key, DistHash, StoreHash>;
  #elif (pKmerStore == BIMOLECULE)  // bimolecule
    template <typename Key>
    using MapParams = ::bliss::index::kmer::BimoleculeHashMapParams<Key, DistHash, StoreHash>;
  #endif

	#if (pDNA == 5) || (pKmerStore == CANONICAL)
	using Splitter = ::fsc::TruePredicate;
	#else
	using Splitter = typename ::std::conditional<(KmerType::nBits == (KmerType::nWords * sizeof(typename KmerType::KmerWordType) * 8)),
			MSBSplitter, ::fsc::TruePredicate>::type;
	#endif


  // DEFINE THE MAP TYPE base on the type of data to be stored.
  #if (pINDEX == POS) || (pINDEX == POSQUAL)  // multimap
    #if (pMAP == VEC)
      using MapType = ::dsc::unordered_multimap_vec<
          KmerType, ValType, MapParams>;

    #elif (pMAP == UNORDERED)
      using MapType = ::dsc::unordered_multimap<
          KmerType, ValType, MapParams>;
    #elif (pMAP == COMPACTVEC)
      using MapType = ::dsc::unordered_multimap_compact_vec<
          KmerType, ValType, MapParams>;
    #elif (pMAP == HASHEDVEC)
      using MapType = ::dsc::unordered_multimap_hashvec<
          KmerType, ValType, MapParams>;
    #elif (pMAP == DENSEHASH)
      using MapType = ::dsc::densehash_multimap<
          KmerType, ValType, MapParams, Splitter>;
    #endif
  #elif (pINDEX == COUNT)  // map
    #if (pMAP == DENSEHASH)
      using MapType = ::dsc::counting_densehash_map<
        KmerType, ValType, MapParams, Splitter>;
    #else
      using MapType = ::dsc::counting_unordered_map<
        KmerType, ValType, MapParams>;
    #endif
  #endif


#endif




//================ FINALLY, the actual index type.

#if (pINDEX == POS)
	using IndexType = bliss::index::kmer::PositionIndex<MapType>;

#elif (pINDEX == POSQUAL)
  using IndexType = bliss::index::kmer::PositionQualityIndex<MapType>;

#elif (pINDEX == COUNT)  // map
	using IndexType = bliss::index::kmer::CountIndex<MapType>;
#endif



/*
 * BENCHMARK Kmer index building.
 *
 * variables:   hash function   (std, murmur, farm)
 *              store canonical or not (canonicalize on the fly, canonicalize on build/query, no op on build, query doubled)
 *              k (15, 21, 31, 63)
 *              backing container type (hashmap, unordered vecmap, sorted array)
 *
 *              file reader type - mpi-io, mmap, fileloader without prefetch.  mpi-io performs really well when file's cached.
 *
 */

template <typename IndexType, typename KmerType = typename IndexType::KmerType>
std::vector<KmerType> readForQuery(const std::string & filename, MPI_Comm comm) {

  ::std::vector<KmerType> query;

  IndexType::template read_file<PARSER_TYPE, ::bliss::index::kmer::KmerParser<KmerType> >(filename, query, comm);

  return query;
}


template <typename IndexType, typename KmerType = typename IndexType::KmerType>
std::vector<KmerType> readForQuery_mpiio(const std::string & filename, MPI_Comm comm) {

  ::std::vector<KmerType> query;

  IndexType::template read_file_mpiio<PARSER_TYPE, ::bliss::index::kmer::KmerParser<KmerType> >(filename, query, comm);

  return query;
}

template <typename IndexType, typename KmerType = typename IndexType::KmerType>
std::vector<KmerType> readForQuery_mmap(const std::string & filename, MPI_Comm comm) {

  ::std::vector<KmerType> query;

  // default to including quality score iterators.
  IndexType::template read_file_mmap<PARSER_TYPE, ::bliss::index::kmer::KmerParser<KmerType> >(filename, query, comm);

  return query;
}


template <typename IndexType, typename KmerType = typename IndexType::KmerType>
std::vector<KmerType> readForQuery_posix(const std::string & filename, MPI_Comm comm) {

  ::std::vector<KmerType> query;

  // default to including quality score iterators.
  IndexType::template read_file_posix<PARSER_TYPE, ::bliss::index::kmer::KmerParser<KmerType> >(filename, query, comm);

  return query;
}


template<typename KmerType>
void sample(std::vector<KmerType> &query, size_t n, unsigned int seed, mxx::comm const & comm) {
  std::shuffle(query.begin(), query.end(), std::default_random_engine(seed));

  size_t n_p = (n / comm.size());
  std::vector<size_t> send_counts(comm.size(), n_p);

  if (n < static_cast<size_t>(comm.size())) {
    n_p = 1;

    for (size_t i = 0; i < n; ++i) {
      send_counts[(i + comm.rank()) % comm.size()] = 1;
    }
    for (int i = n; i < comm.size(); ++i) {
      send_counts[(i + comm.rank()) % comm.size()] = 0;
    }
  }

  std::vector<KmerType> out = ::mxx::all2allv(query, send_counts, comm);
  query.swap(out);
}


/**
 *
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char** argv) {

  //////////////// init logging
  LOG_INIT();


  //////////////// initialize MPI and openMP

  mxx::env e(argc, argv);
  mxx::comm comm;

  if (comm.rank() == 0) printf("EXECUTING %s\n", argv[0]);

  comm.barrier();


  //////////////// parse parameters

  std::string filename;
  filename.assign(PROJ_SRC_DIR);
#if (pPARSER == FASTA)
      filename.append("/test/data/test.fasta");
#elif (pPARSER == FASTQ)
      filename.append("/test/data/test.fastq");
#endif
  std::string queryname(filename);

  int sample_ratio = 100;

  int reader_algo = -1;
  // Wrap everything in a try block.  Do this every time,
  // because exceptions will be thrown for problems.
  try {

    // Define the command line object, and insert a message
    // that describes the program. The "Command description message"
    // is printed last in the help text. The second argument is the
    // delimiter (usually space) and the last one is the version number.
    // The CmdLine object parses the argv array based on the Arg objects
    // that it contains.
    TCLAP::CmdLine cmd("Benchmark parallel kmer index building", ' ', "0.1");

    // Define a value argument and add it to the command line.
    // A value arg defines a flag and a type of value that it expects,
    // such as "-n Bishop".
    TCLAP::ValueArg<std::string> fileArg("F", "file", "FASTQ file path", false, filename, "string", cmd);
    TCLAP::ValueArg<std::string> queryArg("Q", "query", "FASTQ file path for query. default to same file as index file", false, "", "string", cmd);

    TCLAP::ValueArg<int> algoArg("A",
                                 "algo", "Reader Algorithm id. Fileloader w/o preload = 2, mmap = 5, posix=7, piio = 10. default is 7.",
                                 false, 7, "int", cmd);

    TCLAP::ValueArg<int> sampleArg("S",
                                 "query-sample", "sampling ratio for the query kmers. default=100",
                                 false, sample_ratio, "int", cmd);


    // Parse the argv array.
    cmd.parse( argc, argv );

    // Get the value parsed by each arg.
    queryname = queryArg.getValue();   // get this first
    if (queryname.empty()) // at default  set to same as input.
    	queryname = fileArg.getValue();

    filename = fileArg.getValue();
    reader_algo = algoArg.getValue();
    sample_ratio = sampleArg.getValue();

    // set the default for query to filename, and reparse



    // Do what you intend.

  } catch (TCLAP::ArgException &e)  // catch any exceptions
  {
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    exit(-1);
  }





  // ================  read and get file

  IndexType idx(comm);

#if (pMAP == DENSEHASH)
  KmerType empty_key = ::bliss::kmer::hash::sparsehash::empty_key<KmerType>::generate();
  KmerType deleted_key = ::bliss::kmer::hash::sparsehash::deleted_key<KmerType>::generate();

  	idx.get_map().reserve_keys(empty_key, deleted_key);

  	// upper key is negation of lower keys
  	KmerType upper_empty_key = empty_key;
  	KmerType upper_deleted_key = deleted_key;
  	for (size_t i = 0; i < KmerType::nWords; ++i) {
  		upper_empty_key.getDataRef()[i] = ~(upper_empty_key.getDataRef()[i]);
  		upper_deleted_key.getDataRef()[i] = ~(upper_deleted_key.getDataRef()[i]);
  	}

  	idx.get_map().reserve_upper_keys(upper_empty_key, upper_deleted_key, Splitter());


#endif

  BL_BENCH_INIT(test);

  if (comm.rank() == 0) printf("reading query %s via posix\n", queryname.c_str());
  BL_BENCH_START(test);
  auto query = readForQuery_posix<IndexType>(queryname, comm);
  BL_BENCH_COLLECTIVE_END(test, "read_query", query.size(), comm);

  BL_BENCH_START(test);
  sample(query, query.size() / sample_ratio, comm.rank(), comm);
  BL_BENCH_COLLECTIVE_END(test, "sample", query.size(), comm);


  {
	  ::std::vector<typename IndexType::TupleType> temp;

	  BL_BENCH_START(test);
	  if (reader_algo == 2)
	  {
		if (comm.rank() == 0) printf("reading %s via fileloader\n", filename.c_str());

		idx.read_file<PARSER_TYPE, typename IndexType::KmerParserType>(filename, temp, comm);

	  } else if (reader_algo == 5) {
		if (comm.rank() == 0) printf("reading %s via mmap\n", filename.c_str());
		idx.read_file_mmap<PARSER_TYPE, typename IndexType::KmerParserType>(filename, temp, comm);

	  } else if (reader_algo == 7) {
		if (comm.rank() == 0) printf("reading %s via posix\n", filename.c_str());
		idx.read_file_posix<PARSER_TYPE, typename IndexType::KmerParserType>(filename, temp, comm);

	  } else if (reader_algo == 10){
		if (comm.rank() == 0) printf("reading %s via mpiio\n", filename.c_str());
		idx.read_file_mpiio<PARSER_TYPE, typename IndexType::KmerParserType>(filename, temp, comm);
	  } else {
		throw std::invalid_argument("missing file reader type");
	  }
	  BL_BENCH_COLLECTIVE_END(test, "read", temp.size(), comm);

	  size_t total = mxx::allreduce(temp.size(), comm);
	  if (comm.rank() == 0) printf("total size is %lu\n", total);

	  BL_BENCH_START(test);
	  idx.insert(temp);
	  BL_BENCH_COLLECTIVE_END(test, "insert", idx.local_size(), comm);

    total = idx.size();
    if (comm.rank() == 0) printf("total size after insert/rehash is %lu\n", total);
  }

  {

	  {
		  auto lquery = query;
		  BL_BENCH_START(test);
		  auto counts = idx.count(lquery);
		  BL_BENCH_COLLECTIVE_END(test, "count", counts.size(), comm);
	  }
#ifndef pCompare
	  {
		  auto lquery = query;
		  BL_BENCH_START(test);
		  auto found = idx.find(lquery);
		  BL_BENCH_COLLECTIVE_END(test, "find", found.size(), comm);
	  }
#endif
#if 0
	  // separate test because of it being potentially very slow depending on imbalance.
	  {
		  auto lquery = query;

	  BL_BENCH_START(test);
	  auto found = idx.find_collective(lquery);
	  BL_BENCH_COLLECTIVE_END(test, "find_collective", found.size(), comm);
	  }
#endif
	    {
	      auto lquery = query;

	    BL_BENCH_START(test);
	    auto found = idx.find_overlap(lquery);
	    BL_BENCH_COLLECTIVE_END(test, "find_overlap", found.size(), comm);
	    }
    // separate test because of it being potentially very slow depending on imbalance.
#if 0
    {
      auto lquery = query;

    BL_BENCH_START(test);
    auto found = idx.find_sendrecv(lquery);
    BL_BENCH_COLLECTIVE_END(test, "find_sendrecv", found.size(), comm);
    }
#endif

	  BL_BENCH_START(test);
	  idx.erase(query);
	  BL_BENCH_COLLECTIVE_END(test, "erase", idx.local_size(), comm);

  }

  
  BL_BENCH_REPORT_MPI_NAMED(test, "app", comm);


  // mpi cleanup is automatic
  comm.barrier();

  return 0;
}
