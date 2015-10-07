/**
 * @file    test_threads.cpp
 * @ingroup
 * @author  Tony Pan <tpan7@gatech.edu>
 * @brief
 * @details
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
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

#include "wip/kmer_index.hpp"

#include "utils/timer.hpp"


template <typename IndexType, typename KmerType = typename IndexType::KmerType, template<typename> class PreCanonicalizer=bliss::kmer::transform::identity>
std::vector<KmerType> readForQuery(const std::string & filename, MPI_Comm comm) {

  ::std::vector<KmerType> query;

  std::string extension = ::bliss::utils::file::get_file_extension(filename);
  std::transform(extension.begin(), extension.end(), extension.begin(), ::tolower);
  if (extension.compare("fastq") == 0) {
    // default to including quality score iterators.
    IndexType::template read_file<::bliss::io::FASTQParser, ::bliss::index::kmer::KmerParser<KmerType>, PreCanonicalizer >(filename, query, comm);
  } else {
    throw std::invalid_argument("input filename extension is not supported.");
  }

  return query;
}


template <typename IndexType, typename KmerType = typename IndexType::KmerType, template<typename> class PreCanonicalizer=bliss::kmer::transform::identity>
std::vector<KmerType> readForQuery_subcomm(const std::string & filename, MPI_Comm comm) {

  ::std::vector<KmerType> query;

  std::string extension = ::bliss::utils::file::get_file_extension(filename);
  std::transform(extension.begin(), extension.end(), extension.begin(), ::tolower);
  if (extension.compare("fastq") == 0) {
    // default to including quality score iterators.
    IndexType::template read_file_mpi_subcomm<::bliss::io::FASTQParser, ::bliss::index::kmer::KmerParser<KmerType>, PreCanonicalizer >(filename, query, comm);
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
//        if (read.seqBegin == read.seqEnd) continue;
//
//        //== set up the kmer generating iterators.
//        KmerIterType start(BaseCharIterator(read.seqBegin, bliss::common::ASCII2<Alphabet>()), true);
//        KmerIterType end(BaseCharIterator(read.seqEnd, bliss::common::ASCII2<Alphabet>()), false);
//
//
//        query.insert(query.end(), start, end);
//        //        for (auto it = index_start; it != index_end; ++it) {
//        //          temp.push_back(*it);
//        //        }
//        //std::copy(index_start, index_end, temp.end());
//        //        INFOF("R %d inserted.  new temp size = %lu", rank, temp.size());
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
void testIndex(MPI_Comm comm, const std::string & filename, std::string test ) {

  int nprocs = 1;
  int rank = 0;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &rank);

  IndexType idx(comm, nprocs);

  TIMER_INIT(test);

  if (rank == 0) printf("RANK %d / %d: Testing %s", rank, nprocs, test.c_str());

  TIMER_START(test);
  idx.template build<SeqParser>(filename, comm);
  TIMER_END(test, "build", idx.local_size());



  TIMER_START(test);
  auto query = readForQuery<IndexType>(filename, comm);
  TIMER_END(test, "read query", query.size());


  // for testing, query 1% (else could run out of memory.  if a kmer exists r times, then we may need r^2/p total storage.
  TIMER_START(test);
  unsigned seed = rank * 23;
  sample(query, query.size() / 100, seed);
  TIMER_END(test, "select 1%", query.size());

  auto query_orig = query;

  auto query1 = query_orig;
  query1.resize(1);

  // query 1
  TIMER_START(test);
  auto results3 = idx.find(query1);
  TIMER_END(test, "query 1", results3.size());

  query1 = query_orig;
  query1.resize(1);

  // query 1
  TIMER_START(test);
  auto results4 = idx.count(query1);
  TIMER_END(test, "count 1", results4.size());



  query = query_orig;

  // process query
  // query
  TIMER_START(test);
  auto results = idx.find(query);
  TIMER_END(test, "query 1%", results.size());

  query = query_orig;

  // count
  TIMER_START(test);
  auto results2 = idx.count(query);
  TIMER_END(test, "count 1%", results2.size());




  TIMER_REPORT_MPI(test, rank, comm);

}




template <typename IndexType, template <typename> class SeqParser, template<typename> class PreCanonicalizer=bliss::kmer::transform::lex_less>
void testIndexPrecomputeCanonical(MPI_Comm comm, const std::string & filename, std::string test ) {

  int nprocs = 1;
  int rank = 0;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &rank);

  IndexType idx(comm, nprocs);

  TIMER_INIT(test);

  if (rank == 0) printf("RANK %d / %d: Testing %s", rank, nprocs, test.c_str());

  TIMER_START(test);
  idx.template build<SeqParser, PreCanonicalizer>(filename, comm);
  TIMER_END(test, "build", idx.local_size());



  TIMER_START(test);
  auto query = readForQuery<IndexType, typename IndexType::KmerType, PreCanonicalizer>(filename, comm);
  TIMER_END(test, "read query", query.size());


  // for testing, query 1% (else could run out of memory.  if a kmer exists r times, then we may need r^2/p total storage.
  TIMER_START(test);
  unsigned seed = rank * 23;
  sample(query, query.size() / 100, seed);
  TIMER_END(test, "select 1%", query.size());

  auto query_orig = query;

  auto query1 = query_orig;
  query1.resize(1);

  // query 1
  TIMER_START(test);
  auto results3 = idx.find(query1);
  TIMER_END(test, "query 1", results3.size());

  query1 = query_orig;
  query1.resize(1);

  // query 1
  TIMER_START(test);
  auto results4 = idx.count(query1);
  TIMER_END(test, "count 1", results4.size());



  query = query_orig;

  // process query
  // query
  TIMER_START(test);
  auto results = idx.find(query);
  TIMER_END(test, "query 1%", results.size());

  query = query_orig;

  // count
  TIMER_START(test);
  auto results2 = idx.count(query);
  TIMER_END(test, "count 1%", results2.size());




  TIMER_REPORT_MPI(test, rank, comm);

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

  //std::string filename("/home/tpan/src/bliss/test/data/test.medium.fastq");
  std::string filename("/home/tpan/src/bliss/test/data/test.fastq");

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
      ERRORF("The MPI Library Does not have thread support.");
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  } else {
    MPI_Init(&argc, &argv);
  }

  MPI_Comm comm = MPI_COMM_WORLD;

  MPI_Comm_rank(comm, &rank);

  {
    char hostname[256];
    memset(hostname, 0, 256);
    gethostname(hostname, 256);
    //INFOF("Rank %d hostname [%s]", rank, hostname);
  }
  MPI_Barrier(comm);

  if (rank == 0) INFOF("USE_MPI is set");
#else
  static_assert(false, "MPI used although compilation is not set to use MPI");
#endif
  
//  if (which != -1) std::cin.get();


  using Alphabet = bliss::common::DNA;
  using KmerType = bliss::common::Kmer<21, Alphabet, WordType>;

  using IdType = bliss::common::ShortSequenceKmerId;
  using QualType = float;
  using KmerInfoType = std::pair<IdType, QualType>;


  if (which == -1 || which == 2)
  {
  using MapType = ::dsc::counting_unordered_map<
      KmerType, uint32_t, int,
      bliss::kmer::transform::lex_less,
      bliss::kmer::hash::farm >;
  testIndex<bliss::index::kmer::CountIndex<MapType>, bliss::io::FASTQParser > (comm, filename, "ST, hash, count index.");
    MPI_Barrier(comm);
}

  if (which == -1 || which == 4)
  {
  using MapType = ::dsc::unordered_multimap<
      KmerType, IdType, int,
      bliss::kmer::transform::lex_less,
      bliss::kmer::hash::farm >;
  testIndex<bliss::index::kmer::PositionIndex<MapType>, bliss::io::FASTQParser >(comm, filename, "ST, hash, position index.");
    MPI_Barrier(comm);
}
  if (which == -1 || which == 7)
  {
  using MapType = ::dsc::unordered_multimap<
      KmerType, KmerInfoType, int,
      bliss::kmer::transform::lex_less,
      bliss::kmer::hash::farm >;
  testIndex<bliss::index::kmer::PositionQualityIndex<MapType>, bliss::io::FASTQParser >(comm, filename , "ST, hash, pos+qual index");
    MPI_Barrier(comm);
}

  if (which == -1 || which == 1)
  {
  using MapType = ::dsc::unordered_multimap_vec<
      KmerType, IdType, int,
      bliss::kmer::transform::lex_less,
      bliss::kmer::hash::farm >;
  testIndex<bliss::index::kmer::PositionIndex<MapType>, bliss::io::FASTQParser >(comm, filename, "ST, hashvec, position index.");

  MPI_Barrier(comm);
  }

  if (which == -1 || which == 6)
  {
  using MapType = ::dsc::unordered_multimap_vec<
      KmerType, KmerInfoType, int,
      bliss::kmer::transform::lex_less,
      bliss::kmer::hash::farm >;
  testIndex<bliss::index::kmer::PositionQualityIndex<MapType>, bliss::io::FASTQParser >(comm, filename , "ST, hashvec, pos+qual index");
    MPI_Barrier(comm);
}


  if (which == -1 || which == 3)
  {
  using MapType = ::dsc::counting_sorted_map<
      KmerType, uint32_t, int,
      bliss::kmer::transform::lex_less>;
  testIndex<bliss::index::kmer::CountIndex<MapType>, bliss::io::FASTQParser > (comm, filename, "ST, sort, count index.");
    MPI_Barrier(comm);
}

    if (which == -1 || which == 5)
{
  using MapType = ::dsc::sorted_multimap<
      KmerType, IdType, int,
      bliss::kmer::transform::lex_less>;
  testIndex<bliss::index::kmer::PositionIndex<MapType>, bliss::io::FASTQParser >(comm, filename, "ST, sort, position index.");
    MPI_Barrier(comm);
}

  if (which == -1 || which == 8)
  {
  using MapType = ::dsc::sorted_multimap<
      KmerType, KmerInfoType, int,
      bliss::kmer::transform::lex_less>;
  testIndex<bliss::index::kmer::PositionQualityIndex<MapType>, bliss::io::FASTQParser >(comm, filename , "ST, sort, pos+qual index");
    MPI_Barrier(comm);
}

  if (which == -1 || which == 9)
  {
  using MapType = ::dsc::unordered_multimap_compact_vec<
      KmerType, IdType, int,
      bliss::kmer::transform::lex_less,
      bliss::kmer::hash::farm >;
  testIndex<bliss::index::kmer::PositionIndex<MapType>, bliss::io::FASTQParser >(comm, filename, "ST, hashvec2, position index.");

  MPI_Barrier(comm);
  }

  if (which == -1 || which == 10)
  {
  using MapType = ::dsc::unordered_multimap_compact_vec<
      KmerType, KmerInfoType, int,
      bliss::kmer::transform::lex_less,
      bliss::kmer::hash::farm >;
  testIndex<bliss::index::kmer::PositionQualityIndex<MapType>, bliss::io::FASTQParser >(comm, filename , "ST, hashvec2, pos+qual index");
    MPI_Barrier(comm);
}


  if (which == -1 || which == 11)
  {
  using MapType = ::dsc::counting_unordered_map<
      KmerType, uint32_t, int,
      bliss::kmer::transform::identity,
      bliss::kmer::hash::cpp_std >;
      testIndexPrecomputeCanonical<bliss::index::kmer::CountIndex<MapType>, bliss::io::FASTQParser, bliss::kmer::transform::lex_less > (comm, filename, "ST, hash, count index precanonical, std hash.");
      MPI_Barrier(comm);
}


  if (which == -1 || which == 12)
  {
  using MapType = ::dsc::counting_unordered_map<
      KmerType, uint32_t, int,
      bliss::kmer::transform::identity,
      bliss::kmer::hash::identity >;
      testIndexPrecomputeCanonical<bliss::index::kmer::CountIndex<MapType>, bliss::io::FASTQParser, bliss::kmer::transform::lex_less > (comm, filename, "ST, hash, count index, precononical, iden hash.");
      MPI_Barrier(comm);
}

  if (which == -1 || which == 13)
  {
  using MapType = ::dsc::counting_unordered_map<
      KmerType, uint32_t, int,
      bliss::kmer::transform::identity,
      bliss::kmer::hash::murmur >;
      testIndexPrecomputeCanonical<bliss::index::kmer::CountIndex<MapType>, bliss::io::FASTQParser, bliss::kmer::transform::lex_less > (comm, filename, "ST, hash, count index, precononical, murmur hash.");
      MPI_Barrier(comm);
}

  if (which == -1 || which == 14)
  {
  using MapType = ::dsc::counting_unordered_map<
      KmerType, uint32_t, int,
      bliss::kmer::transform::xor_rev_comp,
      bliss::kmer::hash::farm >;
      testIndexPrecomputeCanonical<bliss::index::kmer::CountIndex<MapType>, bliss::io::FASTQParser, bliss::kmer::transform::identity > (comm, filename, "ST, hash, count index, xor in flight, farm hash.");
      MPI_Barrier(comm);
}




//
//
//  if (which == 9)
//  {
//  using MapType = ::dsc::counting_map<
//      KmerType, uint32_t, int,
//      bliss::kmer::transform::lex_less,
//      bliss::kmer::hash::farm >;
//  testIndex<bliss::index::kmer::CountIndex<MapType>, bliss::io::FASTQParser > (comm, filename, "ST, map, count index.");
//    MPI_Barrier(comm);
//}
//
//  if (which == 10)
//  {
//  using MapType = ::dsc::multimap<
//      KmerType, IdType, int,
//      bliss::kmer::transform::lex_less,
//      bliss::kmer::hash::farm >;
//  testIndex<bliss::index::kmer::PositionIndex<MapType>, bliss::io::FASTQParser >(comm, filename, "ST, map, position index.");
//    MPI_Barrier(comm);
//}
//
//  if (which == 11)
//  {
//  using MapType = ::dsc::multimap<
//      KmerType, KmerInfoType, int,
//      bliss::kmer::transform::lex_less,
//      bliss::kmer::hash::farm >;
//  testIndex<bliss::index::kmer::PositionQualityIndex<MapType>, bliss::io::FASTQParser >(comm, filename , "ST, map, pos+qual index");
//    MPI_Barrier(comm);
//}

  //////////////  clean up MPI.
  MPI_Finalize();

  //INFOF("M Rank %d called MPI_Finalize", rank);



  return 0;
}
