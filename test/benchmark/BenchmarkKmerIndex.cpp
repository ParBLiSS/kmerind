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

#define LEX 1
#define XOR 2
#define PRELEX 3

#define STD 1
#define MURMUR 2
#define FARM 3

#define POS 1
#define POSQUAL 2
#define COUNT 3

#define SORTED 1
#define ORDERED 2
#define VEC 3
#define COMPACTVEC 4
#define UNORDERED 5

//================= define types - changeable here...


#if (pPARSER == FASTA)
#define PARSER_TYPE ::bliss::io::FASTAParser
#else
#define PARSER_TYPE ::bliss::io::FASTQParser
#endif

#if (pDNA == 16)
using Alphabet = bliss::common::DNA16;
#elif (pDNA == 5)
using Alphabet = bliss::common::DNA5;
#else
using Alphabet = bliss::common::DNA;
#endif

#if defined(pK)
using KmerType = bliss::common::Kmer<pK, Alphabet, WordType>;
#else
using KmerType = bliss::common::Kmer<21, Alphabet, WordType>;
#endif



#if (pTRANS == LEX)
	template <typename KM>
	using KmerTrans = bliss::kmer::transform::lex_less<KM>;

	template <typename KM>
	using PreCanonicalizer = bliss::kmer::transform::identity<KM>;
#elif (pTRANS == XOR)
	template <typename KM>
	using KmerTrans = bliss::kmer::transform::xor_rev_comp<KM>;

	template <typename KM>
	using PreCanonicalizer = bliss::kmer::transform::identity<KM>;
#elif (pTRANS == PRELEX)
	template <typename KM>
	using KmerTrans = bliss::kmer::transform::identity<KM>;

	template <typename KM>
	using PreCanonicalizer = bliss::kmer::transform::lex_less<KM>;
#else
	template <typename KM>
	using KmerTrans = bliss::kmer::transform::identity<KM>;

	template <typename KM>
	using PreCanonicalizer = bliss::kmer::transform::identity<KM>;
#endif


#if (pHASH == STD)
	template <typename KM, bool pre=false>
	using KmerHash = bliss::kmer::hash::cpp_std<KM, pre>;
#elif (pHASH == IDEN)
	template <typename KM, bool pre=false>
	using KmerHash = bliss::kmer::hash::identity<KM, pre>;
#elif (pHASH == MURMUR)
	template <typename KM, bool pre=false>
	using KmerHash = bliss::kmer::hash::murmur<KM, pre>;
#else
	template <typename KM, bool pre=false>
	using KmerHash = bliss::kmer::hash::farm<KM, pre>;
#endif


#if (pPARSER == FASTA)
	using IdType = bliss::common::LongSequenceKmerId;
#else
	using IdType = bliss::common::ShortSequenceKmerId;
#endif

using QualType = float;
using KmerInfoType = std::pair<IdType, QualType>;
using CountType = uint32_t;

#if (pINDEX == POS)
	using ValType = IdType;
#elif (pINDEX == POSQUAL)
	using ValType = KmerInfoType;
#else
	using ValType = CountType;
#endif

#if (pINDEX == POS) || (pINDEX == POSQUAL)
	#if (pMAP == SORTED)
		using MapType = ::dsc::sorted_multimap<
			  KmerType, ValType, int,
			  KmerTrans>;

	#elif (pMAP == ORDERED)
		using MapType = ::dsc::multimap<
			  KmerType, ValType, int,
			  KmerTrans,
			  KmerHash>;

	#elif (pMAP == VEC)
		using MapType = ::dsc::unordered_multimap_vec<
			  KmerType, ValType, int,
			  KmerTrans,
			  KmerHash>;

	#elif (pMAP == COMPACTVEC)
		using MapType = ::dsc::unordered_multimap_compact_vec<
			  KmerType, ValType, int,
			  KmerTrans,
			  KmerHash>;
	#else
		using MapType = ::dsc::unordered_multimap<
			  KmerType, ValType, int,
			  KmerTrans,
			  KmerHash>;
	#endif

#else
	#if (pMAP == SORTED)
	  using MapType = ::dsc::counting_sorted_map<
		  KmerType, ValType, int,
		  KmerTrans>;
	#elif (pMAP == ORDERED)
	  using MapType = ::dsc::counting_map<
		  KmerType, ValType, int,
		  KmerTrans,
		  KmerHash>;
	#else
	  using MapType = ::dsc::counting_unordered_map<
		  KmerType, ValType, int,
		  KmerTrans,
		  KmerHash>;
	#endif

#endif


#if (pINDEX == POS)
	using IndexType = bliss::index::kmer::PositionIndex<MapType>;

#elif (pINDEX == POSQUAL)
  using IndexType = bliss::index::kmer::PositionQualityIndex<MapType>;

#else
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

  IndexType::template read_file<PARSER_TYPE, PreCanonicalizer, ::bliss::index::kmer::KmerParser<KmerType> >(filename, query, comm);

  return query;
}


template <typename IndexType, typename KmerType = typename IndexType::KmerType>
std::vector<KmerType> readForQuery_mpiio(const std::string & filename, MPI_Comm comm) {

  ::std::vector<KmerType> query;

  IndexType::template read_file_mpiio<PARSER_TYPE, PreCanonicalizer, ::bliss::index::kmer::KmerParser<KmerType> >(filename, query, comm);

  return query;
}

template <typename IndexType, typename KmerType = typename IndexType::KmerType>
std::vector<KmerType> readForQuery_mmap(const std::string & filename, MPI_Comm comm) {

  ::std::vector<KmerType> query;

  // default to including quality score iterators.
  IndexType::template read_file_mmap<PARSER_TYPE, PreCanonicalizer, ::bliss::index::kmer::KmerParser<KmerType> >(filename, query, comm);

  return query;
}


template <typename IndexType, typename KmerType = typename IndexType::KmerType>
std::vector<KmerType> readForQuery_posix(const std::string & filename, MPI_Comm comm) {

  ::std::vector<KmerType> query;

  // default to including quality score iterators.
  IndexType::template read_file_posix<PARSER_TYPE, PreCanonicalizer, ::bliss::index::kmer::KmerParser<KmerType> >(filename, query, comm);

  return query;
}


template<typename KmerType>
void sample(std::vector<KmerType> &query, size_t n, unsigned int seed) {
  std::shuffle(query.begin(), query.end(), std::default_random_engine(seed));
  query.erase(query.begin() + n, query.end());
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
#else
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
                                 false, 100, "int", cmd);


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

  BL_BENCH_INIT(test);

  {
	  ::std::vector<typename IndexType::TupleType> temp;

	  BL_BENCH_START(test);
	  if (reader_algo == 2)
	  {
		if (comm.rank() == 0) printf("reading %s via fileloader\n", filename.c_str());

		idx.read_file<PARSER_TYPE, PreCanonicalizer, typename IndexType::KmerParserType>(filename, temp, comm);

	  } else if (reader_algo == 5) {
		if (comm.rank() == 0) printf("reading %s via mmap\n", filename.c_str());
		idx.read_file_mmap<PARSER_TYPE, PreCanonicalizer, typename IndexType::KmerParserType>(filename, temp, comm);

	  } else if (reader_algo == 7) {
		if (comm.rank() == 0) printf("reading %s via posix\n", filename.c_str());
		idx.read_file_posix<PARSER_TYPE, PreCanonicalizer, typename IndexType::KmerParserType>(filename, temp, comm);

	  } else if (reader_algo == 10){
		if (comm.rank() == 0) printf("reading %s via mpiio\n", filename.c_str());
		idx.read_file_mpiio<PARSER_TYPE, PreCanonicalizer, typename IndexType::KmerParserType>(filename, temp, comm);
	  } else {
		throw std::invalid_argument("missing file reader type");
	  }
	  BL_BENCH_COLLECTIVE_END(test, "read", temp.size(), comm);

	  BL_BENCH_START(test);
	  if (temp.size() > 0)
		idx.insert(temp);
	  BL_BENCH_COLLECTIVE_END(test, "insert", temp.size(), comm);
  }

  {
	  if (comm.rank() == 0) printf("reading query %s via posix\n", queryname.c_str());
	  BL_BENCH_START(test);
	  auto query = readForQuery_posix<IndexType>(queryname, comm);
	  BL_BENCH_COLLECTIVE_END(test, "read_query", query.size(), comm);

	  BL_BENCH_START(test);
	  sample(query, query.size() / sample_ratio, comm.rank());
	  BL_BENCH_COLLECTIVE_END(test, "sample", query.size(), comm);

	  {
	  BL_BENCH_START(test);
	  auto counts = idx.count(query);
	  BL_BENCH_COLLECTIVE_END(test, "count", counts.size(), comm);
	  }

	  {
	  BL_BENCH_START(test);
	  auto found = idx.find(query);
	  BL_BENCH_COLLECTIVE_END(test, "find", found.size(), comm);
	  }

	  {
	  BL_BENCH_START(test);
	  auto found = idx.find_collective(query);
	  BL_BENCH_COLLECTIVE_END(test, "find_collective", found.size(), comm);
	  }

	  BL_BENCH_START(test);
	  idx.erase(query);
	  printf("new local size = %lu\n", idx.local_size());
	  BL_BENCH_COLLECTIVE_END(test, "erase", idx.local_size(), comm);
  }

  
  BL_BENCH_REPORT_MPI_NAMED(test, "app", comm);


  // mpi cleanup is automatic
  comm.barrier();

  return 0;
}
