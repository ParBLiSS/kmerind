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

#if defined(USE_MPI)
#include "mpi.h"
#endif

#if defined(USE_OPENMP)
#include "omp.h"
#endif

#include "utils/logging.h"


//#include <unistd.h>
//#include <string.h>
//#include <cstdio>
//#include <cmath>
//
//#include <iostream>
////#include <thread>
//#include <vector>
//#include <queue>
//#include <unordered_map>
//#include <chrono>

#include <common/base_types.hpp>
#include <common/alphabets.hpp>
#include <common/AlphabetTraits.hpp>
#include <partition/range.hpp>
#include <io/fastq_loader.hpp>
#include <index/kmer_index_element.hpp>
#include <index/kmer_index_functors.hpp>
#include <index/kmer_index_generator.hpp>
#include <io/CommunicationLayer.hpp>
#include <index/distributed_map.hpp>


/*
 * TYPE DEFINITIONS
 */

////////// COMMON TYPES

/// define kmer index type
typedef bliss::index::KmerSize<21>                                      KmerSize;
typedef uint64_t                                                        KmerType;
typedef float                                                           QualityType;
typedef DNA                                                             Alphabet;
typedef bliss::io::FASTQSequenceId                                    IdType;

/// define Range type
typedef bliss::partition::range<size_t>                                 RangeType;

/// define the input file type
// raw data type :  use CharType
typedef bliss::io::FASTQLoader<CharType, false, true>                   FileLoaderType;
typedef typename FileLoaderType::L2BlockType::iterator                      FileBlockIterType;

/// define read type
typedef bliss::io::SequenceWithQuality<FileBlockIterType, Alphabet, QualityType>  SeqType;

/// define the transform iterator type
typedef bliss::io::FASTQParser<FileBlockIterType, Alphabet, QualityType>    ParserType;
typedef bliss::io::SequencesIterator<ParserType, FileBlockIterType>             SeqIterType;

/// define kmer quality GENERATOR types
typedef bliss::index::SangerToLogProbCorrect<QualityType>               QualityEncoderType;
typedef bliss::index::generate_qual<SeqType, KmerSize, QualityType, QualityEncoderType > QualOpType;


///////////// INDEX TYPES
/// kmer index element type and corresponding generator type
typedef bliss::index::KmerIndexElementWithIdAndQuality<KmerSize, KmerType, IdType, QualityType>   KmerIndexValueType;
typedef bliss::index::generate_kmer<SeqType, KmerIndexValueType>                                  KmerOpType;

/// define the index storage type
typedef bliss::index::distributed_multimap<KmerType, KmerIndexValueType, bliss::io::CommunicationLayer> IndexType;
/// define the KmerIndexElement Generator type.
typedef bliss::index::KmerIndexGeneratorWithQuality<KmerOpType, IndexType, QualOpType>            KmerIndexComputeType;



////////////// COUNT TYPES
/// kmer index element type and corresponding generator type
typedef bliss::index::KmerIndexElement<KmerSize, KmerType>                                        KmerCountType;
typedef bliss::index::generate_kmer<SeqType, KmerCountType>                                       KmerCountOpType;

/// count indices
typedef bliss::index::distributed_counting_map<KmerType, bliss::io::CommunicationLayer>           CountIndexType;
/// define the KmerIndexElement Generator for counting maptype.
typedef bliss::index::KmerIndexGenerator<KmerCountOpType, CountIndexType>                         KmerCountComputeType;




// TODO: make these function with/without MPI, with/without OpenMP.
// TODO: 1. make a communicator that encapsulates the nprocs, rank, nthread, tid, etc - this refactors KmerIndexGenerator so there is no dependency on nprocs, rank, and tid.
// TODO: 2. refactor these "compute" methods so that the loop iterations are handled by a ThreadModel.
// TODO: 3. make the communicator short circuit when copying to local, and lock when multithreading.
// TODO: 4. communicator should have send buffer and receive buffer (actual data structure)

/*
 * first pattern for threading - has a master thread  - does not work well.  need each thread to get its own chunk range based on tid, so that means
 * each task needs to call the getNextChunkRange method.  but then the master thread cannot check the length of the chunk received.
 * the alternative is for master thread to use a counter. but if I have a counter, then might as well use parfor.
 *
 * function removed.  for this pattern, see omp_patterns.hpp.
 */

/*
 * peer to peer
 */
template <typename Compute, typename Index>
void buildIndex(FileLoaderType &loader, Index &index, const int& rank,
                const int& nthreads, const int& chunkSize)
{
  std::chrono::high_resolution_clock::time_point t1, t2;
  std::chrono::duration<double> time_span;
  t1 = std::chrono::high_resolution_clock::now();

  INFO("buildIndex rank.tid = " << rank << "." << omp_get_thread_num());

  // local variables for information only.
  int nReads = 0;    // read count
  int nChunks = 0;    // chunk count

  // uses the fastq iterator as the queue itself, instead of master/slave.
  //   at this point, no strong difference.
#pragma omp parallel num_threads(nthreads) shared(loader, nthreads, index, rank) reduction(+:nReads, nChunks) OMP_SHARE_DEFAULT
  {
    /// instantiate a local parser
    ParserType parser;

    int tid = 0;
#ifdef USE_OPENMP
    tid = omp_get_thread_num();
#endif

    /// instantiate a local kmer generator
    Compute op;

    /// local variables for loop
    typename FileLoaderType::L2BlockType chunk;
    SeqType read;

    /// initialize the loop by getting the first chunk
    chunk = loader.getNextL2Block(tid);
    while (chunk.getRange().size() > 0) {
      /// get the chunk of data

      ///  and wrap the chunk inside an iterator that emits Reads.
      SeqIterType fastq_start(parser, chunk.begin(), chunk.end(), chunk.getRange().start);
      SeqIterType fastq_end(chunk.end());

      /// loop over the reads
      for (; fastq_start != fastq_end; ++fastq_start)
      {
        // first get read
        read = *fastq_start;

        // then compute and store into index (this will generate kmers and insert into index)
        op(read, index);
        ++nReads;

        // do a little status print.
        if (nReads % 20000 == 0)
          INFO("buildIndex rank.tid=" << rank << "." <<  tid << " nReads=" << nReads);
      }

      /// get read for next loop iteration.
      chunk = loader.getNextL2Block(tid);
      ++nChunks;
    }

    INFO("buildIndex rank.tid=" << rank << "." << tid << " nChunks=" << nChunks);

  }  // compute threads parallel

  t2 = std::chrono::high_resolution_clock::now();
  time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
  INFO("buildIndex rank=" << rank << " nReads=" << nReads << " nChunks=" << nChunks << " elapsed time: " << time_span.count() << "s.");

  /// this MPI process is done.  now flush the index to all other nodes.
  t1 = std::chrono::high_resolution_clock::now();
  index.flush();
  t2 = std::chrono::high_resolution_clock::now();
  time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
  INFO("flushIndex rank=" << rank << " elapsed time: " << time_span.count() << "s.");

}



template<typename ComputeType, typename Index>
struct RunTask {
    void operator()(const std::string &filename, Index &index, MPI_Comm &comm, const int &nthreads, const int &chunkSize) {
      /////////////// initialize output variables

      int rank;
      int nprocs;

      MPI_Comm_size(comm, &nprocs);
      MPI_Comm_rank(comm, &rank);


      std::chrono::high_resolution_clock::time_point t1, t2;
      std::chrono::duration<double> time_span;


      //////////////// now partition and open the file.
      // create FASTQ Loader
      // call get NextPartition Range to block partition,
      // call "load" with the partition range.


      // get the file ready for read
      {
        t1 = std::chrono::high_resolution_clock::now();

        FileLoaderType loader(filename, comm, nthreads, chunkSize);  // this handle is alive through the entire execution.

        loader.getNextL1Block();


      t2 = std::chrono::high_resolution_clock::now();
      time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
      std::cout << "MMap rank " << rank << " elapsed time: " << time_span.count() << "s." << std::endl;

      std::cout << rank << " file partition: " << loader.getCurrentL1Block().getRange() << std::endl;



      /////////////// now process the file .
      // do some work using openmp  // version without master
      t1 = std::chrono::high_resolution_clock::now();
      buildIndex<ComputeType, Index>(loader, index, rank, nthreads, chunkSize);

      t2 = std::chrono::high_resolution_clock::now();
      time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
      std::cout << "Compute rank " << rank << " elapsed time: " << time_span.count() << "s." << std::endl;


      }  // scope to ensure file loader is destroyed.
    }
};





/**
 *
 * @param argc
 * @param argv
 * @return
 *
 *
 *
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
  printf("NOT compiled with OPENMP")
#endif

  int chunkSize = 4096;
  if (argc > 2)
  {
    chunkSize = atoi(argv[2]);
  }

  std::string filename("/home/tpan/src/bliss/test/data/test.fastq");
  //  std::string filename("/mnt/data/1000genome/HG00096/sequence_read/SRR077487_1.filt.fastq");
  if (argc > 3)
  {
    filename.assign(argv[3]);
  }

  int groupSize = 1;
  int id = 0;
  //////////////// initialize MPI and openMP
#ifdef USE_MPI


  if (nthreads > 1) {

    int provided;

    // one thread will be making all MPI calls.
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);

    if (provided < MPI_THREAD_FUNNELED) {
      printf("ERROR: The MPI Library Does not have thread support.\n");
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  } else {
    MPI_Init(&argc, &argv);
  }

  MPI_Comm comm = MPI_COMM_WORLD;

  MPI_Comm_size(comm, &groupSize);
  MPI_Comm_rank(comm, &id);


  if (id == 0)
    std::cout << "USE_MPI is set" << std::endl;
#else
  static_assert(false, "MPI used although compilation is not set to use MPI");
#endif
  // replace with MPIRunner


  {  // scoped to ensure index is deleted before MPI_Finalize()
    IndexType index(comm, groupSize);
    //index.setLooupAnswerCallback(std::function<void(std::pair<KmerType, std::vector<KmerIndexValueType> >&)>(&callback));


    /////////////// start processing.  enclosing with braces to make sure loader is destroyed before MPI finalize.
    RunTask<KmerIndexComputeType, IndexType> t;
    t(filename, index, comm, nthreads, chunkSize);

    printf("MPI number of entries in index for rank %d is %lu\n", id, index.local_size());



    //// query:  use the same file as input.  walk through and generate kmers as before.  send query

  }
  MPI_Barrier(comm);

////  {  // scoped to ensure index is deleted before MPI_Finalize()
////    CountIndexType index(comm, groupSize);
////    //index.setLooupAnswerCallback(std::function<void(std::pair<KmerType, std::vector<sKmerIndexValueType> >&)>(&callback));
////
////
////    /////////////// start processing.  enclosing with braces to make sure loader is destroyed before MPI finalize.
////    RunTask<KmerCountComputeType, CountIndexType> t;
////    t(filename, index, comm, nthreads, chunkSize);
////
////    printf("MPI number of entries in index for rank %d is %lu\n", id, index.local_size());
////
////    // TODO:  need to change how kmer indices returned from iterator.  do NOT compute the key as revcomp ^ kmer - let map do it.
////    // TODO:  need consistent interface for returning from iterator when there are additional information besides kmer.
////
////    //// query:  use the same file as input.  walk through and generate kmers as before.  send query
////
////  }
//
//  MPI_Barrier(comm);

  //////////////  clean up MPI.
  MPI_Finalize();

  return 0;
}
