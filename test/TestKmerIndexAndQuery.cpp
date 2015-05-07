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

    DEBUGF("testQueryOld nthreads is %d", nthreads);
    int rank;
        MPI_Comm_rank(comm, &rank);
    std::chrono::high_resolution_clock::time_point t1, t2;
    std::chrono::duration<double> time_span;

    t1 = std::chrono::high_resolution_clock::now();

    typename KmerIndexType::FileLoaderType loader(comm, filename, nthreads, chunkSize);  // this handle is alive through the entire execution.

    t2 = std::chrono::high_resolution_clock::now();
    time_span =
         std::chrono::duration_cast<std::chrono::duration<double> >(
             t2 - t1);
     INFOF("R %d query file open time: %f", rank, time_span.count());

     t1 = std::chrono::high_resolution_clock::now();


    typename KmerIndexType::FileLoaderType::L1BlockType partition = loader.getNextL1Block();

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
    t2 = std::chrono::high_resolution_clock::now();
     time_span =
         std::chrono::duration_cast<std::chrono::duration<double>>(
             t2 - t1);
     INFOF("R %d query file load/kmer gen time: %f", rank, time_span.count());

     t1 = std::chrono::high_resolution_clock::now();

    kmer_index.flushQuery();

    t2 = std::chrono::high_resolution_clock::now();
     time_span =
         std::chrono::duration_cast<std::chrono::duration<double>>(
             t2 - t1);
     INFOF("R %d query flush time: %f", rank, time_span.count());

  }  // scope to ensure file loader is destroyed.
}

// TODO: currently first read only.
template<typename KmerIndexType>
void testQuery(MPI_Comm comm, const std::string & filename, const int nthreads, const int chunkSize, KmerIndexType& kmer_index)
{
  {
//            t1 = std::chrono::high_resolution_clock::now();
    // create FASTQ Loader
    DEBUGF("testQuery nthreads is %d", nthreads);
    int rank;
    MPI_Comm_rank(comm, &rank);


    std::chrono::high_resolution_clock::time_point t1, t2;
    std::chrono::duration<double> time_span;

    t1 = std::chrono::high_resolution_clock::now();
    typename KmerIndexType::FileLoaderType loader(comm, filename, nthreads, chunkSize);  // this handle is alive through the entire execution.
    t2 = std::chrono::high_resolution_clock::now();
     time_span =
         std::chrono::duration_cast<std::chrono::duration<double> >(
             t2 - t1);
     INFOF("R %d query file open time: %f", rank, time_span.count());

     t1 = std::chrono::high_resolution_clock::now();
    typename KmerIndexType::FileLoaderType::L1BlockType partition = loader.getNextL1Block();

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
          typename KmerIndexType::KmerIterType start(typename KmerIndexType::BaseCharIterator(read.seqBegin, bliss::common::ASCII2<bliss::common::DNA>()), true);
          typename KmerIndexType::KmerIterType end(typename KmerIndexType::BaseCharIterator(read.seqEnd, bliss::common::ASCII2<bliss::common::DNA>()), false);

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

    t2 = std::chrono::high_resolution_clock::now();
     time_span =
         std::chrono::duration_cast<std::chrono::duration<double>>(
             t2 - t1);
     INFOF("R %d query file load/kmer gen time: %f", rank, time_span.count());

     t1 = std::chrono::high_resolution_clock::now();
    kmer_index.flushQuery();
    t2 = std::chrono::high_resolution_clock::now();
     time_span =
         std::chrono::duration_cast<std::chrono::duration<double>>(
             t2 - t1);
     INFOF("R %d query flush time: %f", rank, time_span.count());


  }  // scope to ensure file loader is destroyed.

}



template<typename IndexType>
void test(MPI_Comm comm, const std::string & filename, const int nthreads, const int chunkSize,
             const std::function<void(MPI_Comm comm, const std::string & , const int , const int, IndexType&)> &queryFunction,
             std::string testname) {

  int nprocs = 1;
  int rank = 0;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &rank);

  DEBUGF("test nthreads is %d", nthreads);



  std::chrono::high_resolution_clock::time_point t1, t2;
  std::chrono::duration<double> time_span;
  std::vector<std::string> timespan_names;
  std::vector<double> timespans;

  size_t entries = 0;

  INFOF("RANK %d: Testing %s", rank, testname.c_str());

  t1 = std::chrono::high_resolution_clock::now();
  // initialize index
  DEBUGF("RANK %d: ***** initializing %s.", rank, testname.c_str());

  double callback_time = 0;

  IndexType kmer_index(comm, nprocs,
                                     std::bind(IndexType::defaultReceiveAnswerCallback,
                                               _1, _2, nthreads, std::ref(entries), std::ref(callback_time)),
                                     nthreads);

  // start processing.  enclosing with braces to make sure loader is destroyed before MPI finalize.
  DEBUGF("RANK %d: ***** building index first pass.  %d threads, callback_time %f ", rank, nthreads, callback_time);

  kmer_index.build(filename, nthreads, chunkSize);
  //kmer_index.flush();

//  MPI_Barrier(comm);
  INFO("RANK " << rank << " Index Building 1 for " << filename << " using " << nthreads << " threads, index size " << kmer_index.local_size());

  t2 = std::chrono::high_resolution_clock::now();
  time_span =
      std::chrono::duration_cast<std::chrono::duration<double>>(
          t2 - t1);
  timespan_names.push_back("build");
  timespans.push_back(time_span.count());


  // === test query.  make each node read 1.
  // get the file ready for read
  DEBUGF("RANK %d: ******* query index", rank);
  callback_time = 0;
  t1 = std::chrono::high_resolution_clock::now();

  queryFunction(comm, filename, nthreads, chunkSize, kmer_index);

  t2 = std::chrono::high_resolution_clock::now();
  time_span =
      std::chrono::duration_cast<std::chrono::duration<double>>(
          t2 - t1);
  timespan_names.push_back("query");
  timespans.push_back(time_span.count());

  timespan_names.push_back("response portion");
  timespans.push_back(callback_time);


  // get the file ready for read
  DEBUGF("RANK %d: ******* query index 2", rank);
  callback_time = 0;
  t1 = std::chrono::high_resolution_clock::now();

  queryFunction(comm, filename, nthreads, chunkSize, kmer_index);

  t2 = std::chrono::high_resolution_clock::now();
  time_span =
      std::chrono::duration_cast<std::chrono::duration<double>>(
          t2 - t1);
  timespan_names.push_back("query");
  timespans.push_back(time_span.count());

  timespan_names.push_back("response portion");
  timespans.push_back(callback_time);


  std::stringstream ss;
  std::copy(timespan_names.begin(), timespan_names.end(), std::ostream_iterator<std::string>(ss, ","));
  std::stringstream ss2;
  std::copy(timespans.begin(), timespans.end(), std::ostream_iterator<double>(ss2, ","));

  std::stringstream ss3;
  std::vector<size_t> sizes = kmer_index.local_sizes();
  std::copy(sizes.begin(), sizes.end(), std::ostream_iterator<size_t>(ss3, ","));


  INFOF("Rank %d Test %s phases [%s] times [%s] size %lu (by thread [%s])", rank, testname.c_str(), ss.str().c_str(), ss2.str().c_str(), kmer_index.local_size(), ss3.str().c_str() );

  kmer_index.finalize();


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



#if defined(KMOLECULEINDEX)
  {
    using IndexType = bliss::index::retired::KmerPositionIndexOld<21, bliss::common::DNA, bliss::io::FASTQ>;
    test<IndexType>(comm, filename, nthreads, chunkSize, testQueryOld<IndexType>, "Kmolecule pos index");

  }
#elif defined(KMERINDEX)
  {
    using IndexType = bliss::index::KmerPositionIndex<21, bliss::common::DNA, bliss::io::FASTQ>;
    test<IndexType>(comm, filename, nthreads, chunkSize, testQuery<IndexType>, "Kmer pos index");

  }


#endif
  MPI_Barrier(comm);



#if defined(KMOLECULEINDEX)
  {
  using IndexType = bliss::index::retired::KmerCountIndexOld<21, bliss::common::DNA, bliss::io::FASTQ>;
  test<IndexType>(comm, filename, nthreads, chunkSize, testQueryOld<IndexType>, "Kmolecule count index");

  }

#elif defined(KMERINDEX)

  {
  using IndexType = bliss::index::KmerCountIndex<bliss::common::Kmer<21, bliss::common::DNA>, bliss::io::FASTQ>;
  test<IndexType>(comm, filename, nthreads, chunkSize, testQuery<IndexType>, "Kmer count index");

  }

#endif
  MPI_Barrier(comm);


#if defined(KMOLECULEINDEX)

  // with quality score....
  {
  using IndexType = bliss::index::retired::KmerPositionAndQualityIndexOld<21, bliss::common::DNA, bliss::io::FASTQ>;

  test<IndexType>(comm, filename, nthreads, chunkSize, testQueryOld<IndexType>, "Kmolecule pos+qual index");

  }
#elif defined(KMERINDEX)

  {
  using IndexType = bliss::index::KmerPositionAndQualityIndex<21, bliss::common::DNA, bliss::io::FASTQ>;

  test<IndexType>(comm, filename, nthreads, chunkSize, testQuery<IndexType>, "Kmer pos+qual index");

  }

#endif
  MPI_Barrier(comm);

  //////////////  clean up MPI.
  MPI_Finalize();

  INFOF("M Rank %d called MPI_Finalize", rank);



  return 0;
}
