/**
 * @file    test_threads.cpp
 * @ingroup
 * @author  tpan
 * @brief
 * @details
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#include "config.hpp"

#include <functional>
#include "utils/logging.h"

#include "common/alphabets.hpp"


#include "index/kmer_index.hpp"
#include "common/kmer.hpp"
#include "common/base_types.hpp"
#include <string>
#include <sstream>
#include "utils/kmer_utils.hpp"
/*
 * TYPE DEFINITIONS
 */

using namespace std::placeholders;   // for _1, _2, _3

template<typename KmerIndexType>
void testQueryOld(MPI_Comm comm, const std::string & filename, const int nthreads, const int chunkSize, KmerIndexType& kmer_index) {
  {
//            t1 = std::chrono::high_resolution_clock::now();
    // create FASTQ Loader

    typename KmerIndexType::FileLoaderType loader(comm, filename, nthreads, chunkSize);  // this handle is alive through the entire execution.
    typename KmerIndexType::FileLoaderType::L1BlockType partition;

    partition = loader.getNextL1Block();

    //===  repeatedly load the next L1 Block.
    if (partition.getRange().size() > 0) {


      // get L2Block from L1Block, spread work to all threads.
      //== instantiate a local parser in each thread
      typename KmerIndexType::ParserType parser;

      int tid = 0;
      //== local variables for loop
      typename KmerIndexType::FileLoaderType::L2BlockType chunk;


      //== process L2Blocks in a loop.  loader.getNextL2Block is atomic.  the work is shared by all threads.
      chunk = loader.getNextL2Block(tid);
      if (chunk.getRange().size() > 0) {


        //== process the chunk of data
        typename KmerIndexType::SeqType read;

        //==  and wrap the chunk inside an iterator that emits Reads.
        typename KmerIndexType::SeqIterType fastq_start(parser, chunk.begin(), chunk.end(), chunk.getRange().start);

          // first get read
          read = *fastq_start;

          if (read.seqBegin != read.seqEnd) {

            typename KmerIndexType::KmoleculeOpType kmer_op;
            bliss::utils::KmoleculeToCanonicalKmerFunctor<typename KmerIndexType::KmerType> transform;
            typename KmerIndexType::KmerIterType start(typename KmerIndexType::KmoleculeIterType(read.seqBegin, kmer_op), transform);
            typename KmerIndexType::KmerIterType end(typename KmerIndexType::KmoleculeIterType(read.seqEnd, kmer_op), transform);

            unsigned int k = KmerIndexType::KmerType::size;
//            int i = 0;
//            for (; i < k; ++i, ++start) {};
//            DEBUGF(" skip over first %d", i);
//            for (; start != end; ++start) {};
//            {
//              DEBUGF(" kmer: %s", bliss::utils::KmerUtils::toASCIIString(*start).c_str());
//              //kmer = *start;
//              kmer_index.sendQuery(*start);  // right side either RVO assignment (if iterator/functor created the object) or copy (if returning a cached object)
//                                    // then copy assign
//            }
            unsigned int i = 0;
  //          KmerType idx;
            for (; i < (k - 1) && start != end; ++start, ++i);  // compute but discard the first K - 1.

            for (; start != end; ++start)
            {
//              DEBUGF(" kmer: %s", bliss::utils::KmerUtils::toASCIIString(*start).c_str());
              //kmer = *start;
              kmer_index.sendQuery(*start);  // right side either RVO assignment (if iterator/functor created the object) or copy (if returning a cached object)
                                    // then copy assign
            }
          }
      }
    }
    kmer_index.flushQuery();
  }  // scope to ensure file loader is destroyed.
}

template<typename KmerIndexType>
void testQuery(MPI_Comm comm, const std::string & filename, const int nthreads, const int chunkSize, KmerIndexType& kmer_index) {
  {
//            t1 = std::chrono::high_resolution_clock::now();
    // create FASTQ Loader

    typename KmerIndexType::FileLoaderType loader(comm, filename, nthreads, chunkSize);  // this handle is alive through the entire execution.
    typename KmerIndexType::FileLoaderType::L1BlockType partition;

    partition = loader.getNextL1Block();

    //===  repeatedly load the next L1 Block.
    if (partition.getRange().size() > 0) {


      // get L2Block from L1Block, spread work to all threads.
      //== instantiate a local parser in each thread
      typename KmerIndexType::ParserType parser;

      int tid = 0;
      //== local variables for loop
      typename KmerIndexType::FileLoaderType::L2BlockType chunk;


      //== process L2Blocks in a loop.  loader.getNextL2Block is atomic.  the work is shared by all threads.
      chunk = loader.getNextL2Block(tid);
      if (chunk.getRange().size() > 0) {


        //== process the chunk of data
        typename KmerIndexType::SeqType read;

        //==  and wrap the chunk inside an iterator that emits Reads.
        typename KmerIndexType::SeqIterType fastq_start(parser, chunk.begin(), chunk.end(), chunk.getRange().start);

          // first get read
          read = *fastq_start;

          if (read.seqBegin != read.seqEnd) {

          //== set up the kmer generating iterators.
          typename KmerIndexType::KmerIterType start(typename KmerIndexType::BaseCharIterator(read.seqBegin, bliss::ASCII2<DNA>()), true);
          typename KmerIndexType::KmerIterType end(typename KmerIndexType::BaseCharIterator(read.seqEnd, bliss::ASCII2<DNA>()), false);

  //          int kmerCount = 0;

  //          int tid = omp_get_thread_num();

            // no need to discard first K-1

            for (; start != end; ++start)
            {
              //kmer = *start;
              kmer_index.sendQuery(*start);  // right side either RVO assignment (if iterator/functor created the object) or copy (if returning a cached object)
                                    // then copy assign
            }

          }
      }


    }
    kmer_index.flushQuery();
  }  // scope to ensure file loader is destroyed.

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

  int nthreads = 1;
#ifdef USE_OPENMP
  if (argc > 1)
  {
    nthreads = atoi(argv[1]);
    if (nthreads == -1)
      nthreads = omp_get_max_threads();
  }
  omp_set_nested(1);
  omp_set_dynamic(0);
#else
  FATALF("NOT compiled with OPENMP")
#endif

  int chunkSize = sysconf(_SC_PAGE_SIZE);
  if (argc > 2)
  {
    chunkSize = atoi(argv[2]);
  }

  //std::string filename("/home/tpan/src/bliss/test/data/test.medium.fastq");
  std::string filename("/home/tpan/src/bliss/test/data/test.fastq");
  //std::string filename("/mnt/data/1000genome/HG00096/sequence_read/SRR077487_1.filt.fastq");
  if (argc > 3)
  {
    filename.assign(argv[3]);
  }

  int nprocs = 1;
  int rank = 0;
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

  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &rank);


  if (rank == 0)
    INFOF("USE_MPI is set");
#else
  static_assert(false, "MPI used although compilation is not set to use MPI");
#endif


  std::chrono::high_resolution_clock::time_point t1, t2;
  std::chrono::duration<double> time_span;

  size_t result = 0;
  size_t entries = 0;
  double score = 0;

#if defined(KMOLECULEINDEX)
  {
    t1 = std::chrono::high_resolution_clock::now();
    // initialize index
    INFOF("***** initializing index.");

    using IndexType = bliss::index::retired::KmerPositionIndexOld<21, DNA, bliss::io::FASTQ, true>;
    IndexType kmer_index(comm, nprocs,
                                       std::bind(IndexType::defaultReceivePositionAnswer,
                                                 _1, _2, nthreads, std::ref(result), std::ref(entries)),
                                       nthreads);

    // start processing.  enclosing with braces to make sure loader is destroyed before MPI finalize.
    INFOF("***** building index first pass.");

    kmer_index.build(filename, nthreads, chunkSize);
    //kmer_index.flush();

  //  MPI_Barrier(comm);
    //INFO("COUNT " << rank << " Index Building 1 for " << filename << " using " << nthreads << " threads, index size " << kmer_index.local_size());
    INFO("RANK " << rank << " Index Building 1 for " << filename << " using " << nthreads << " threads, index size " << kmer_index.local_size());

    t2 = std::chrono::high_resolution_clock::now();
    time_span =
        std::chrono::duration_cast<std::chrono::duration<double>>(
            t2 - t1);
    INFO("new position index: time " << time_span.count());

    // === test query.  make each node read 1.
    // get the file ready for read
    INFOF("******* query index");
    t1 = std::chrono::high_resolution_clock::now();

    testQueryOld(comm, filename, nthreads, chunkSize, kmer_index);

    t2 = std::chrono::high_resolution_clock::now();
    time_span =
        std::chrono::duration_cast<std::chrono::duration<double>>(
            t2 - t1);
    INFO("new position index query: time " << time_span.count());

    kmer_index.finalize();

  }
#elif defined(KMERINDEX)
  {
    t1 = std::chrono::high_resolution_clock::now();


  // initialize index
  INFOF("***** initializing index.");
  using IndexType = bliss::index::KmerPositionIndex<21, DNA, bliss::io::FASTQ, true>;
  IndexType kmer_index(comm, nprocs,
                                    std::bind(IndexType::defaultReceivePositionAnswer, // <bliss::Kmer<21, DNA>, bliss::io::FASTQSequenceId >
                                              _1, _2, nthreads, std::ref(result), std::ref(entries)),
                                              nthreads);

  // start processing.  enclosing with braces to make sure loader is destroyed before MPI finalize.
  INFOF("***** building index first pass.");

  kmer_index.build(filename, nthreads, chunkSize);
  //kmer_index.flush();

//  MPI_Barrier(comm);
  //INFO("COUNT " << rank << " Index Building 1 for " << filename << " using " << nthreads << " threads, index size " << kmer_index.local_size());
  INFO("RANK " << rank << " Index Building 1 for " << filename << " using " << nthreads << " threads, index size " << kmer_index.local_size());

  t2 = std::chrono::high_resolution_clock::now();
  time_span =
      std::chrono::duration_cast<std::chrono::duration<double>>(
          t2 - t1);
  INFO("new position index: time " << time_span.count());

  // === test query.  make each node read 1.
  // get the file ready for read
  INFOF("******* query index");
  t1 = std::chrono::high_resolution_clock::now();

  testQuery(comm, filename, nthreads, chunkSize, kmer_index);

  t2 = std::chrono::high_resolution_clock::now();
  time_span =
      std::chrono::duration_cast<std::chrono::duration<double>>(
          t2 - t1);
  INFO("new position index query: time " << time_span.count());
  kmer_index.finalize();

  }


#endif
  MPI_Barrier(comm);


  result = 0;
  entries = 0;
  score = 0;

#if defined(KMOLECULEINDEX)
  {
    t1 = std::chrono::high_resolution_clock::now();
  // initialize index
  INFOF("***** initializing index.");
  using IndexType = bliss::index::retired::KmerCountIndexOld<21, DNA, bliss::io::FASTQ, true>;
  IndexType kmer_index(comm, nprocs,
                           std::bind(&IndexType::defaultReceiveCountAnswer,
                                     _1, _2, nthreads, std::ref(result), std::ref(entries)),
                                     nthreads);

  // start processing.  enclosing with braces to make sure loader is destroyed before MPI finalize.
  INFOF("***** building index first pass.");

  kmer_index.build(filename, nthreads, chunkSize);
  //kmer_index.flush();

//  MPI_Barrier(comm);
  //INFO("COUNT " << rank << " Index Building 1 for " << filename << " using " << nthreads << " threads, index size " << kmer_index.local_size());
  INFO("RANK " << rank << " Index Building 1 for " << filename << " using " << nthreads << " threads, index size " << kmer_index.local_size());


  for (auto it = kmer_index.getLocalIndex().cbegin(); it != kmer_index.getLocalIndex().cend(); ++it) {
    INFOF("Entry: Key=%s (%s), val = %d", it->first.toString().c_str(), bliss::utils::KmerUtils::toASCIIString(it->first).c_str(), it->second);
  }
  t2 = std::chrono::high_resolution_clock::now();
  time_span =
      std::chrono::duration_cast<std::chrono::duration<double>>(
          t2 - t1);
  INFO("old count index: time " << time_span.count());

  //=== query
  // get the file ready for read

  INFOF("***** query index.");
  t2 = std::chrono::high_resolution_clock::now();

  testQueryOld(comm, filename, nthreads, chunkSize, kmer_index);

  t2 = std::chrono::high_resolution_clock::now();
  time_span =
      std::chrono::duration_cast<std::chrono::duration<double>>(
          t2 - t1);
  INFO("old count index query: time " << time_span.count());
  kmer_index.finalize();

  }

#elif defined(KMERINDEX)

  {
    t1 = std::chrono::high_resolution_clock::now();
  // initialize index
  INFOF("***** initializing index.");
  using IndexType = bliss::index::KmerCountIndex<21, DNA, bliss::io::FASTQ, true>;
  IndexType  kmer_index(comm, nprocs,
                              std::bind(&IndexType::defaultReceiveCountAnswer,
                                        _1, _2, nthreads, std::ref(result), std::ref(entries)),
                                        nthreads);

  // start processing.  enclosing with braces to make sure loader is destroyed before MPI finalize.
  INFOF("***** building index first pass.");

  kmer_index.build(filename, nthreads, chunkSize);
  //kmer_index.flush();

//  MPI_Barrier(comm);
  //INFO("COUNT " << rank << " Index Building 1 for " << filename << " using " << nthreads << " threads, index size " << kmer_index.local_size());
  INFO("RANK " << rank << " Index Building 1 for " << filename << " using " << nthreads << " threads, index size " << kmer_index.local_size());

  for (auto it = kmer_index.getLocalIndex().cbegin(); it != kmer_index.getLocalIndex().cend(); ++it) {
    INFOF("Entry: Key=%s (%s), val = %d", it->first.toString().c_str(), bliss::utils::KmerUtils::toASCIIString(it->first).c_str(), it->second);
  }

  t2 = std::chrono::high_resolution_clock::now();
  time_span =
      std::chrono::duration_cast<std::chrono::duration<double>>(
          t2 - t1);
  INFO("new count index: time " << time_span.count());

  // === test query.  make each node read 1.
  // get the file ready for read
  INFOF("******* query index");
  t1 = std::chrono::high_resolution_clock::now();

  testQuery(comm, filename, nthreads, chunkSize, kmer_index);

  t2 = std::chrono::high_resolution_clock::now();
  time_span =
      std::chrono::duration_cast<std::chrono::duration<double>>(
          t2 - t1);
  INFO("new count index query: time " << time_span.count());
  kmer_index.finalize();


  }

#endif
  MPI_Barrier(comm);


entries = 0;
result = 0;
score = 0;

#if defined(KMOLECULEINDEX)

  // with quality score....
  {
    t1 = std::chrono::high_resolution_clock::now();
  // initialize index
  INFOF("***** initializing index.");
  using IndexType = bliss::index::retired::KmerPositionAndQualityIndexOld<21, DNA, bliss::io::FASTQ, true>;
  IndexType kmer_index(comm, nprocs,
                                     std::bind(&IndexType::defaultReceivePositionAndQualityAnswer,
                                               _1, _2, nthreads, std::ref(result), std::ref(score), std::ref(entries)),
                                              nthreads);

  // start processing.  enclosing with braces to make sure loader is destroyed before MPI finalize.
  INFOF("***** building index first pass.");

  kmer_index.build(filename, nthreads, chunkSize);
  //kmer_index.flush();

//  MPI_Barrier(comm);
  //INFO("COUNT " << rank << " Index Building 1 for " << filename << " using " << nthreads << " threads, index size " << kmer_index.local_size());
  INFO("RANK " << rank << " Index Building 1 for " << filename << " using " << nthreads << " threads, index size " << kmer_index.local_size());

  t2 = std::chrono::high_resolution_clock::now();
  time_span =
      std::chrono::duration_cast<std::chrono::duration<double>>(
          t2 - t1);
  INFO("old pos + quality index: time " << time_span.count());

  // === test query.  make each node read 1.
  // get the file ready for read
  INFOF("******* query index");
  t1 = std::chrono::high_resolution_clock::now();

  testQueryOld(comm, filename, nthreads, chunkSize, kmer_index);

  t2 = std::chrono::high_resolution_clock::now();
  time_span =
      std::chrono::duration_cast<std::chrono::duration<double>>(
          t2 - t1);
  INFO("old pos + quality index query: time " << time_span.count());

  kmer_index.finalize();


  }
#elif defined(KMERINDEX)

  {
    t1 = std::chrono::high_resolution_clock::now();
  // initialize index
  INFOF("***** initializing index.");
  using IndexType = bliss::index::KmerPositionAndQualityIndex<21, DNA, bliss::io::FASTQ, true>;
  IndexType kmer_index(comm, nprocs,
                                    std::bind(&IndexType::defaultReceivePositionAndQualityAnswer,
                                               _1, _2, nthreads, std::ref(result), std::ref(score), std::ref(entries)),
                                               nthreads);

  // start processing.  enclosing with braces to make sure loader is destroyed before MPI finalize.
  INFOF("***** building index first pass.");

  kmer_index.build(filename, nthreads, chunkSize);
  //kmer_index.flush();

//  MPI_Barrier(comm);
  //INFO("COUNT " << rank << " Index Building 1 for " << filename << " using " << nthreads << " threads, index size " << kmer_index.local_size());
  INFO("RANK " << rank << " Index Building 1 for " << filename << " using " << nthreads << " threads, index size " << kmer_index.local_size());

  t2 = std::chrono::high_resolution_clock::now();
  time_span =
      std::chrono::duration_cast<std::chrono::duration<double>>(
          t2 - t1);
  INFO("new pos + quality index: time " << time_span.count());

  // === test query.  make each node read 1.
  // get the file ready for read
  INFOF("******* query index");
  t1 = std::chrono::high_resolution_clock::now();

  testQuery(comm, filename, nthreads, chunkSize, kmer_index);

  t2 = std::chrono::high_resolution_clock::now();
  time_span =
      std::chrono::duration_cast<std::chrono::duration<double>>(
          t2 - t1);
  INFO("new pos + quality index query: time " << time_span.count());

  kmer_index.finalize();


  }

#endif
  MPI_Barrier(comm);

  //////////////  clean up MPI.
  MPI_Finalize();

  DEBUGF("M Rank %d called MPI_Finalize", rank);



  return 0;
}
