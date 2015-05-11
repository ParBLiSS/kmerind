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
#include "config.hpp"

#include <unistd.h>  // get hostname


#include <functional>
#include <random>
#include <algorithm>
#include "utils/logging.h"

#include "common/alphabets.hpp"


#include "common/kmer.hpp"
#include "common/base_types.hpp"
#include <string>
#include <sstream>
#include "utils/kmer_utils.hpp"
#include <chrono>

#include "io/mxx_support.hpp"
#include "wip/distributed_map.hpp"
#include "io/sequence_iterator.hpp"
#include "io/sequence_id_iterator.hpp"
#include "iterators/transform_iterator.hpp"
#include "common/kmer_iterators.hpp"
#include "iterators/zip_iterator.hpp"
#include "index/quality_score_iterator.hpp"

#include "wip/kmer_index.hpp"


template <typename KmerType>
std::vector<KmerType> readForQuery(const std::string & filename, MPI_Comm comm) {
  using FileLoaderType = bliss::io::FASTQLoader<CharType, true, false>; // raw data type :  use CharType

  //====  now process the file, one L1 block (block partition by MPI Rank) at a time
  // from FileLoader type, get the block iter type and range type
  using FileBlockIterType = typename FileLoaderType::L1BlockType::iterator;

  using ParserType = bliss::io::FASTQParser<FileBlockIterType, void>;
  using SeqType = typename ParserType::SequenceType;
  using SeqIterType = bliss::io::SequencesIterator<ParserType>;

  using Alphabet = typename KmerType::KmerAlphabet;

  /// converter from ascii to alphabet values
  using BaseCharIterator = bliss::iterator::transform_iterator<typename SeqType::IteratorType, bliss::common::ASCII2<Alphabet> >;

  /// kmer generation iterator
  using KmerIterType = bliss::common::KmerGenerationIterator<BaseCharIterator, KmerType>;
  std::vector< KmerType > query;

  int commSize, rank;
  MPI_Comm_size(comm, &commSize);
  MPI_Comm_rank(comm, &rank);


  {
    std::chrono::high_resolution_clock::time_point t1, t2;
    std::chrono::duration<double> time_span;

    t1 = std::chrono::high_resolution_clock::now();
    //==== create file Loader
    FileLoaderType loader(comm, filename, 1, sysconf(_SC_PAGE_SIZE));  // this handle is alive through the entire building process.
    typename FileLoaderType::L1BlockType partition = loader.getNextL1Block();
    size_t est_size = (loader.getKmerCountEstimate(KmerType::size) + commSize - 1) / commSize;

    t2 = std::chrono::high_resolution_clock::now();
    time_span =
        std::chrono::duration_cast<std::chrono::duration<double>>(
            t2 - t1);
    INFOF("R %d file open time: %f", rank, time_span.count());


    // == create kmer iterator
    //            kmer_iter start(data, range);
    //            kmer_iter end(range.second,range.second);

    t1 = std::chrono::high_resolution_clock::now();
    query.reserve(est_size);

    ParserType parser;
    //=== copy into array
    while (partition.getRange().size() > 0) {
      //== process the chunk of data
      SeqType read;

      //==  and wrap the chunk inside an iterator that emits Reads.
      SeqIterType seqs_start(parser, partition.begin(), partition.end(), partition.getRange().start);
      SeqIterType seqs_end(partition.end());


      //== loop over the reads
      for (; seqs_start != seqs_end; ++seqs_start)
      {
        // first get read
        read = *seqs_start;

        // then compute and store into index (this will generate kmers and insert into index)
        if (read.seqBegin == read.seqEnd) continue;

        //== set up the kmer generating iterators.
        KmerIterType start(BaseCharIterator(read.seqBegin, bliss::common::ASCII2<Alphabet>()), true);
        KmerIterType end(BaseCharIterator(read.seqEnd, bliss::common::ASCII2<Alphabet>()), false);


        query.insert(query.end(), start, end);
        //        for (auto it = index_start; it != index_end; ++it) {
        //          temp.push_back(*it);
        //        }
        //std::copy(index_start, index_end, temp.end());
        //        INFOF("R %d inserted.  new temp size = %lu", rank, temp.size());
      }

      partition = loader.getNextL1Block();
    }
    t2 = std::chrono::high_resolution_clock::now();
    time_span =
        std::chrono::duration_cast<std::chrono::duration<double>>(
            t2 - t1);
    INFOF("R %d local kmer array time: %f. query size = %lu", rank, time_span.count(), query.size());



    //            //distributed_map m(element_count);
    //            m.reserve(element_count, communicator);
    //            m.insert(start, end, communicator);
    //            //m.local_rehash();

    //df.close();


  }
  return query;
}



template<typename KmerType>
void sample(std::vector<KmerType> &query, size_t n, unsigned int seed) {
  std::shuffle(query.begin(), query.end(), std::default_random_engine(seed));
  query.erase(query.begin() + n, query.end());
}




template <typename IndexType>
void testIndex(MPI_Comm comm, const std::string & filename, std::string testname ) {

  int nprocs = 1;
  int rank = 0;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &rank);

  DEBUGF("test nthreads is %d", nthreads);


  IndexType idx(comm, nprocs);

  using KmerType = typename IndexType::map_type::key_type;

  std::chrono::high_resolution_clock::time_point t1, t2;
  std::chrono::duration<double> time_span;
  std::vector<std::string> timespan_names;
  std::vector<double> timespans;
  std::vector<size_t> elements;

  size_t entries = 0;

  INFOF("RANK %d: Testing %s", rank, testname.c_str());

  t1 = std::chrono::high_resolution_clock::now();
  // initialize index
  // start processing.  enclosing with braces to make sure loader is destroyed before MPI finalize.
  DEBUGF("RANK %d: ***** building index first pass.  %d threads, callback_time %f ", rank, nthreads, callback_time);

  idx.build(filename, comm);
  INFO("RANK " << rank << " Index Building 1 for " << filename );

  t2 = std::chrono::high_resolution_clock::now();
  time_span =
      std::chrono::duration_cast<std::chrono::duration<double>>(
          t2 - t1);
  timespan_names.push_back("build");
  timespans.push_back(time_span.count());
  elements.push_back(idx.local_size());


  t1 = std::chrono::high_resolution_clock::now();
  auto query = readForQuery<KmerType>(filename, comm);
  t2 = std::chrono::high_resolution_clock::now();
  time_span =
      std::chrono::duration_cast<std::chrono::duration<double>>(
          t2 - t1);
  timespan_names.push_back("read query");
  timespans.push_back(time_span.count());
  elements.push_back(query.size());


  // for testing, query 1% (else could run out of memory.  if a kmer exists r times, then we may need r^2/p total storage.
  t1 = std::chrono::high_resolution_clock::now();
  unsigned seed = rank * 23;
  sample(query, query.size() / 100, seed);
  t2 = std::chrono::high_resolution_clock::now();
  time_span =
      std::chrono::duration_cast<std::chrono::duration<double>>(
          t2 - t1);
  timespan_names.push_back("select 1%%");
  timespans.push_back(time_span.count());
  elements.push_back(query.size());


  // process query

  // query
  t1 = std::chrono::high_resolution_clock::now();

  auto results = idx.find(query);

 t2 = std::chrono::high_resolution_clock::now();
  time_span =
      std::chrono::duration_cast<std::chrono::duration<double> >(
          t2 - t1);
  timespan_names.push_back("query 1%%");
  timespans.push_back(time_span.count());
  elements.push_back(results.size());

  // count
  t1 = std::chrono::high_resolution_clock::now();

  auto results2 = idx.count(query);

 t2 = std::chrono::high_resolution_clock::now();
  time_span =
      std::chrono::duration_cast<std::chrono::duration<double> >(
          t2 - t1);
  timespan_names.push_back("count 1%%");
  timespans.push_back(time_span.count());
  elements.push_back(results2.size());


  // select 1
  // for testing, query 1% (else could run out of memory.  if a kmer exists r times, then we may need r^2/p total storage.
  t1 = std::chrono::high_resolution_clock::now();
  seed = rank * 57;
  sample(query, 1, seed);
  t2 = std::chrono::high_resolution_clock::now();
  time_span =
      std::chrono::duration_cast<std::chrono::duration<double>>(
          t2 - t1);
  timespan_names.push_back("select 1");
  timespans.push_back(time_span.count());
  elements.push_back(query.size());


  // query 1
  t1 = std::chrono::high_resolution_clock::now();

  auto results3 = idx.find(query);

 t2 = std::chrono::high_resolution_clock::now();
  time_span =
      std::chrono::duration_cast<std::chrono::duration<double> >(
          t2 - t1);
  timespan_names.push_back("query 1");
  timespans.push_back(time_span.count());
  elements.push_back(results3.size());


  // query 1
  t1 = std::chrono::high_resolution_clock::now();

  auto results4 = idx.count(query);

 t2 = std::chrono::high_resolution_clock::now();
  time_span =
      std::chrono::duration_cast<std::chrono::duration<double> >(
          t2 - t1);
  timespan_names.push_back("count 1");
  timespans.push_back(time_span.count());
  elements.push_back(results4.size());

  std::stringstream ss;
  std::copy(timespan_names.begin(), timespan_names.end(), std::ostream_iterator<std::string>(ss, ","));
  std::stringstream ss2;
  std::copy(timespans.begin(), timespans.end(), std::ostream_iterator<double>(ss2, ","));
  std::stringstream ss3;
  std::copy(elements.begin(), elements.end(), std::ostream_iterator<size_t>(ss3, ","));


  INFOF("Rank %d Test %s\n\tphases [%s]\n\ttimes [%s]\n\tcount [%s]", rank, testname.c_str(), ss.str().c_str(), ss2.str().c_str(), ss3.str().c_str());


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
  //std::string filename("/home/tpan/src/bliss/test/data/test.fastq");
  std::string filename("/mnt/data/1000genome/HG00096/sequence_read/SRR077487_1.filt.0_0078235.fastq");
  if (argc > 1)
  {
    filename.assign(argv[1]);
  }


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
    INFOF("Rank %d hostname [%s]\n", rank, hostname);
  }
  MPI_Barrier(comm);

  if (rank == 0)
    INFOF("USE_MPI is set");
#else
  static_assert(false, "MPI used although compilation is not set to use MPI");
#endif


  using Alphabet = bliss::common::DNA;
  using KmerType = bliss::common::Kmer<21, Alphabet, WordType>;

  using IdType = bliss::io::FASTQ::SequenceId;
  using MapType = ::dsc::unordered_multimap<KmerType, IdType, int, 
	bliss::hash::kmer::hash<KmerType, bliss::hash::kmer::detail::farm::hash, 
				bliss::hash::kmer::LexicographicLessCombiner> >;
  testIndex<bliss::index::kmer::PositionIndex<MapType> >(comm, filename, "single thread, position index.");

  MPI_Barrier(comm);


  using MapType2 = ::dsc::counting_unordered_map<KmerType, uint32_t, int, 
	bliss::hash::kmer::hash<KmerType, bliss::hash::kmer::detail::farm::hash,
				bliss::hash::kmer::LexicographicLessCombiner> >;
  testIndex<bliss::index::kmer::CountIndex<MapType2> > (comm, filename, "single thread, count index.");

  MPI_Barrier(comm);

  using QualType = float;
  using KmerInfoType = std::pair<IdType, QualType>;
  using MapType3 = ::dsc::unordered_multimap<KmerType, KmerInfoType, int,
	bliss::hash::kmer::hash<KmerType, bliss::hash::kmer::detail::farm::hash,
				bliss::hash::kmer::LexicographicLessCombiner> >;
  testIndex<bliss::index::kmer::PositionQualityIndex<MapType3> >(comm, filename , "single thread, pos+qual index");

  MPI_Barrier(comm);

  //////////////  clean up MPI.
  MPI_Finalize();

  INFOF("M Rank %d called MPI_Finalize", rank);



  return 0;
}
