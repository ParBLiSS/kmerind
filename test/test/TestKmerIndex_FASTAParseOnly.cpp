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

#include <unistd.h>  // get hostname


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


template <typename IndexType, typename KmerType = typename IndexType::KmerType>
std::vector<KmerType> readForQuery(const std::string & filename, MPI_Comm comm) {

  ::std::vector<KmerType> query;

  std::string extension = ::bliss::utils::file::get_file_extension(filename);
  std::transform(extension.begin(), extension.end(), extension.begin(), ::tolower);
  if (extension.compare("fasta") == 0) {
    // default to including quality score iterators.
    IndexType::template read_file<::bliss::io::FASTAParser, bliss::kmer::transform::identity, ::bliss::index::kmer::KmerParser<KmerType> >(filename, query, comm);
  } else {
    throw std::invalid_argument("input filename extension is not supported.");
  }

  return query;
}


template <typename IndexType, typename KmerType = typename IndexType::KmerType>
std::vector<KmerType> readForQuery_direct(const std::string & filename, MPI_Comm comm) {

  ::std::vector<KmerType> query;

  std::string extension = ::bliss::utils::file::get_file_extension(filename);
  std::transform(extension.begin(), extension.end(), extension.begin(), ::tolower);
  if (extension.compare("fasta") == 0) {
    // default to including quality score iterators.
    IndexType::template read_file_direct<::bliss::io::FASTAParser, bliss::kmer::transform::identity, ::bliss::index::kmer::KmerParser<KmerType> >(filename, query, comm);
  } else {
    throw std::invalid_argument("input filename extension is not supported.");
  }

  return query;
}


//template <typename KmerType>
//std::vector<KmerType> readForQuery(const std::string & filename, MPI_Comm comm) {
//  using FileLoaderType = bliss::io::FASTQLoader<CharType, true, false, false>; // raw data type :  use CharType
//
//  //====  now process the file, one L1 block (block partition by MPI Rank) at a time
//  // from FileLoader type, get the block iter type and range type
//  using FileBlockIterType = typename FileLoaderType::L1BlockType::iterator;
//
//  using ParserType = bliss::io::FASTQParser<FileBlockIterType, void>;
//  using SeqType = typename ParserType::SequenceType;
//  using SeqIterType = bliss::io::SequencesIterator<ParserType>;
//
//  using Alphabet = typename KmerType::KmerAlphabet;
//
//  /// converter from ascii to alphabet values
//  using BaseCharIterator = bliss::iterator::transform_iterator<typename SeqType::IteratorType, bliss::common::ASCII2<Alphabet> >;
//
//  /// kmer generation iterator
//  using KmerIterType = bliss::common::KmerGenerationIterator<BaseCharIterator, KmerType>;
//  std::vector< KmerType > query;
//
//  int commSize, rank;
//  MPI_Comm_size(comm, &commSize);
//  MPI_Comm_rank(comm, &rank);
//
//
//  {
//    //==== create file Loader
//    FileLoaderType loader(comm, filename, 1, sysconf(_SC_PAGE_SIZE));  // this handle is alive through the entire building process.
//    typename FileLoaderType::L1BlockType partition = loader.getNextL1Block();
//    size_t est_size = (loader.getKmerCountEstimate(KmerType::size) + commSize - 1) / commSize;
//
//    // == create kmer iterator
//    //            kmer_iter start(data, range);
//    //            kmer_iter end(range.second,range.second);
//
//    query.reserve(est_size);
//
//    ParserType parser;
//    //=== copy into array
//    while (partition.getRange().size() > 0) {
//      //== process the chunk of data
//      SeqType read;
//
//      //==  and wrap the chunk inside an iterator that emits Reads.
//      SeqIterType seqs_start(parser, partition.begin(), partition.end(), partition.getRange().start);
//      SeqIterType seqs_end(partition.end());
//
//
//      //== loop over the reads
//      for (; seqs_start != seqs_end; ++seqs_start)
//      {
//        // first get read
//        read = *seqs_start;
//
//        // then compute and store into index (this will generate kmers and insert into index)
//        if (read.seq_begin == read.seq_end) continue;
//
//        //== set up the kmer generating iterators.
//        KmerIterType start(BaseCharIterator(read.seq_begin, bliss::common::ASCII2<Alphabet>()), true);
//        KmerIterType end(BaseCharIterator(read.seq_end, bliss::common::ASCII2<Alphabet>()), false);
//
//
//        query.insert(query.end(), start, end);
//        //        for (auto it = index_start; it != index_end; ++it) {
//        //          temp.push_back(*it);
//        //        }
//        //std::copy(index_start, index_end, temp.end());
//        //        BL_INFOF("R %d inserted.  new temp size = %lu", rank, temp.size());
//      }
//
//      partition = loader.getNextL1Block();
//    }
//  }
//  return query;
//}



template<typename KmerType>
void sample(std::vector<KmerType> &query, size_t n, unsigned int seed) {
  std::shuffle(query.begin(), query.begin() + ::std::min(4 * n, query.size()), std::default_random_engine(seed));
  query.erase(query.begin() + n, query.end());
}




template <typename IndexType, template <typename> class SeqParser>
void testIndex(const mxx::comm& comm, const std::string & filename, std::string test ) {

  IndexType idx(comm);

  BL_BENCH_INIT(test);

  if (comm.rank() == 0) BL_INFOF("RANK %d / %d: Testing %s", comm.rank(), comm.size(), test.c_str());

  BL_BENCH_START(test);
  idx.template build<SeqParser>(filename, comm);
  BL_BENCH_END(test, "build", idx.local_size());



  BL_BENCH_START(test);
  auto query = readForQuery<IndexType>(filename, comm);
  BL_BENCH_END(test, "read query", query.size());


  // for testing, query 1% (else could run out of memory.  if a kmer exists r times, then we may need r^2/p total storage.
  BL_BENCH_START(test);
  unsigned seed = comm.rank() * 23;
  sample(query, query.size() / 100, seed);
  BL_BENCH_END(test, "select 1%", query.size());

  BL_BENCH_REPORT_MPI(test, comm.rank(), comm);

}




template <typename IndexType, template <typename> class SeqParser>
void testIndexDirect(const mxx::comm& comm, const std::string & filename, std::string test ) {

  IndexType idx(comm);

  BL_BENCH_INIT(test);

  if (comm.rank() == 0) BL_INFOF("RANK %d / %d: Testing %s", comm.rank(), comm.size(), test.c_str());

  BL_BENCH_START(test);
  idx.template build_direct<SeqParser>(filename, comm);
  BL_BENCH_END(test, "build", idx.local_size());



  BL_BENCH_START(test);
  auto query = readForQuery_direct<IndexType>(filename, comm);
  BL_BENCH_END(test, "read query", query.size());


  // for testing, query 1% (else could run out of memory.  if a kmer exists r times, then we may need r^2/p total storage.
  BL_BENCH_START(test);
  unsigned seed = comm.rank() * 23;
  sample(query, query.size() / 100, seed);
  BL_BENCH_END(test, "select 1%", query.size());

  BL_BENCH_REPORT_MPI(test, comm.rank(), comm);

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

  //////////////// parse parameters

  std::string filename;
      filename.assign(PROJ_SRC_DIR);
      filename.append("/test/data/test.fasta");

  if (argc > 1)
  {
    filename.assign(argv[1]);
  }

  int which = -1;
  if (argc > 2)
	which = atoi(argv[2]);


  int rank = 0;
	int nthreads = 1;
  //////////////// initialize MPI and openMP
#ifdef USE_MPI

  if (nthreads > 1) {

    int provided;

    // one thread will be making all MPI calls.
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);

    if (provided < MPI_THREAD_FUNNELED) {
      BL_ERRORF("The MPI Library Does not have thread support.");
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  } else {
    MPI_Init(&argc, &argv);
  }

  MPI_Comm comm = MPI_COMM_WORLD;

  MPI_Comm_rank(comm, &rank);

  MPI_Barrier(comm);

  if (rank == 0) BL_INFOF("USE_MPI is set");
#else
  static_assert(false, "MPI used although compilation is not set to use MPI");
#endif
  
  if (which != -1) std::cin.get();


  using Alphabet = bliss::common::DNA;
  using KmerType = bliss::common::Kmer<21, Alphabet, WordType>;



  if (which == -1 || which == 1)
  {
  using MapType = ::dsc::counting_unordered_map<
      KmerType, uint32_t, int,
      bliss::kmer::transform::lex_less,
      bliss::kmer::hash::farm >;
  testIndex<bliss::index::kmer::CountIndex<MapType>, bliss::io::FASTAParser > (comm, filename, "ST, hash, all read, count index.");
    MPI_Barrier(comm);
}

  if (which == -1 || which == 3)
  {
  using MapType = ::dsc::counting_unordered_map<
      KmerType, uint32_t, int,
      bliss::kmer::transform::lex_less,
      bliss::kmer::hash::farm >;
  testIndexDirect<bliss::index::kmer::CountIndex<MapType>, bliss::io::FASTAParser > (comm, filename, "ST, hash, read_direct, count index.");
    MPI_Barrier(comm);
}

  //////////////  clean up MPI.
  MPI_Finalize();

  //BL_INFOF("M Rank %d called MPI_Finalize", rank);



  return 0;
}
