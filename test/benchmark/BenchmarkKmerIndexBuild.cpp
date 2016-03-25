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

#if defined(USE_FASTQ_PARSER)
#define PARSER_TYPE ::bliss::io::FASTQParser
#elif defined(USE_FASTA_PARSER)
#define PARSER_TYPE ::bliss::io::FASTAParser
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

  IndexType::template read_file<PARSER_TYPE, bliss::kmer::transform::identity, ::bliss::index::kmer::KmerParser<KmerType> >(filename, query, comm);

  return query;
}


template <typename IndexType, typename KmerType = typename IndexType::KmerType>
std::vector<KmerType> readForQuery_mpiio(const std::string & filename, MPI_Comm comm) {

  ::std::vector<KmerType> query;

  IndexType::template read_file_mpiio<PARSER_TYPE, bliss::kmer::transform::identity, ::bliss::index::kmer::KmerParser<KmerType> >(filename, query, comm);

  return query;
}

template <typename IndexType, typename KmerType = typename IndexType::KmerType>
std::vector<KmerType> readForQuery_mmap(const std::string & filename, MPI_Comm comm) {

  ::std::vector<KmerType> query;

  // default to including quality score iterators.
  IndexType::template read_file_mmap<PARSER_TYPE, bliss::kmer::transform::identity, ::bliss::index::kmer::KmerParser<KmerType> >(filename, query, comm);

  return query;
}



template<typename KmerType>
void sample(std::vector<KmerType> &query, size_t n, unsigned int seed) {
  std::shuffle(query.begin(), query.begin() + ::std::min(4 * n, query.size()), std::default_random_engine(seed));
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
  comm.barrier();

  //////////////// parse parameters

  std::string filename;
  filename.assign(PROJ_SRC_DIR);
#if defined(USE_FASTQ_PARSER)
      filename.append("/test/data/test.fastq");
#elif defined(USE_FASTA_PARSER)
      filename.append("/test/data/test.fasta");
#endif

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

    TCLAP::ValueArg<int> algoArg("A",
                                 "algo", "Reader Algorithm id. Fileloader w/o preload = 2, mmap = 5, mpiio = 10. If absent, fileloader is used (2).",
                                 false, 2, "int", cmd);

    // Parse the argv array.
    cmd.parse( argc, argv );

    // Get the value parsed by each arg.
    filename = fileArg.getValue();
    reader_algo = algoArg.getValue();

    // Do what you intend.

  } catch (TCLAP::ArgException &e)  // catch any exceptions
  {
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    exit(-1);
  }

  //================= define types

  using Alphabet = bliss::common::DNA;
  using KmerType = bliss::common::Kmer<21, Alphabet, WordType>;



  using MapType = ::dsc::counting_unordered_map<
      KmerType, uint64_t, int,
      bliss::kmer::transform::identity,
      bliss::kmer::hash::farm >;

  using IndexType = bliss::index::kmer::CountIndex<MapType>;

  // ================  read and get file

  ::std::vector<typename IndexType::TupleType> temp;
  IndexType idx(comm);

  BL_BENCH_INIT(test);

  BL_BENCH_START(test);
  if (reader_algo == 2)
  {
    if (comm.rank() == 0) printf("reading via fileloader\n");

    idx.read_file<PARSER_TYPE, bliss::kmer::transform::identity, typename IndexType::KmerParserType>(filename, temp, comm);

  } else if (reader_algo == 5) {
    if (comm.rank() == 0) printf("reading via mmap\n");
    idx.read_file_mmap<PARSER_TYPE, bliss::kmer::transform::identity, typename IndexType::KmerParserType>(filename, temp, comm);

  } else if (reader_algo == 10){
    if (comm.rank() == 0) printf("reading via mpiio\n");
    idx.read_file_mpiio<PARSER_TYPE, bliss::kmer::transform::identity, typename IndexType::KmerParserType>(filename, temp, comm);
  } else {
    throw std::invalid_argument("missing file reader type");
  }
  BL_BENCH_END(test, "read", temp.size());

  BL_BENCH_START(test);
  if (temp.size() > 0)
    idx.insert(temp);
  BL_BENCH_END(test, "insert", temp.size());


  BL_BENCH_REPORT_MPI_NAMED(test, "app:build", comm);


  // mpi cleanup is automatic
  comm.barrier();

  return 0;
}
