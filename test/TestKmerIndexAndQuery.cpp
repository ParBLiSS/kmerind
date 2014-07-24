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
/// define kmer index type
typedef bliss::index::KmerSize<21>                                KmerSize;
typedef uint64_t                                                  KmerType;
typedef float                                                     QualityType;
typedef DNA                                                       Alphabet;
typedef bliss::io::fastq_sequence_id                              IdType;
typedef bliss::index::KmerIndexElementWithIdAndQuality<KmerSize, KmerType, IdType, QualityType>
                                                                  KmerIndexType;

/// define Range type
typedef bliss::partition::range<size_t>                           RangeType;

/// define the input file type
// raw data type :  use CharType
typedef bliss::io::FASTQLoader<CharType, true, false>             FileLoaderType;
typedef typename FileLoaderType::BlockIteratorType                BaseIterType;
typedef typename std::iterator_traits<BaseIterType>::value_type   BaseValueType;
/// define the transform iterator type
typedef bliss::io::fastq_parser<BaseIterType, Alphabet, QualityType>
                                                                  ParserType;
typedef bliss::io::fastq_iterator<ParserType, BaseIterType>       ReadIteratorType;

/// define read type
typedef bliss::io::fastq_sequence_quality<BaseIterType, Alphabet, QualityType>
                                                                  SequenceType;


/// define kmer GENERATOR types
typedef bliss::index::generate_kmer<SequenceType, KmerIndexType>  KmerOpType;

/// define kmer quality GENERATOR types
typedef bliss::index::SangerToLogProbCorrect<QualityType>         QualityEncodeType;
typedef bliss::index::generate_qual<SequenceType, KmerSize, QualityType, QualityEncodeType >
                                                                  QualOpType;

/// define the index storage type
typedef bliss::index::distributed_multimap<KmerType, KmerIndexType, bliss::io::CommunicationLayer>
                                                                  IndexType;
/// define the KmerIndexElement Generator type.
typedef bliss::index::KmerIndexGeneratorWithQuality<KmerOpType, IndexType, QualOpType>
                                                                  KmerIndexComputeType;






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
template <typename Compute>
void buildIndex(FileLoaderType &loader, IndexType &index, const int &rank,
                const int &nthreads, const int &chunkSize)
{

  ///  get the start and end iterator position.

  // compute threads
      std::chrono::high_resolution_clock::time_point t1, t2;
      std::chrono::duration<double> time_span;

      INFO("buildIndex rank.tid = " << rank << "." << omp_get_thread_num());

      t1 = std::chrono::high_resolution_clock::now();

      // uses the fastq iterator as the queue itself, instead of master/slave.
      //   at this point, no strong difference.
#pragma omp parallel num_threads(nthreads) shared(loader, nthreads, index) OMP_SHARE_DEFAULT
      {
        ParserType parser;

        int tid = 0;
#ifdef USE_OPENMP
        tid = omp_get_thread_num();
#endif
        Compute op;

        int i = 0;
        int j = 0;
        RangeType r = loader.getNextChunkRange(tid);

        while (r.size() > 0) {

          auto chunk = loader.getChunk(tid, r);

          ReadIteratorType fastq_start(parser, chunk.begin(), chunk.end(), r);
          ReadIteratorType fastq_end(parser, chunk.end(), r);

          SequenceType read;

          for (; fastq_start != fastq_end; ++fastq_start)
          {
            // get data, and move to next for the other threads

            // first get read
            read = *fastq_start;

            // then compute
            op(read, index);
            ++i;

            if (i % 20000 == 0)
              INFO("buildIndex rank.tid = " << rank << "." <<  tid << " processed seq " << i);
          }

          r = loader.getNextChunkRange(tid);
          ++j;

        }

        INFO("buildIndex rank.tid = " << rank << "." << tid << " processed total of " << j << " chunks");

      }  // compute threads parallel


      t2 = std::chrono::high_resolution_clock::now();
      time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
      INFO("buildIndex rank.tid = " << rank << "." << omp_get_thread_num() << " elapsed time: " << time_span.count() << "s.");



      // send the last part out.
      t1 = std::chrono::high_resolution_clock::now();
      index.flush();
      t2 = std::chrono::high_resolution_clock::now();
      time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
      INFO("buildIndex rank.tid = " << rank << "." << omp_get_thread_num() << " elapsed time: " << time_span.count() << "s.");
      //printf("flush rank %d thread %d elapsed time %f\n", rank, omp_get_thread_num(), time_span.count());


      //printf("rank %d done compute\n", rank);
}



template<typename ComputeType>
struct RunTask {
    void operator()(const std::string &filename, IndexType &index, MPI_Comm &comm, const int &nthreads, const int &chunkSize) {
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

        t1 = std::chrono::high_resolution_clock::now();

        // get the file ready for read
        FileLoaderType loader(filename, comm, nthreads, chunkSize);  // this handle is alive through the entire execution.
        RangeType pr = loader.getNextPartitionRange(rank);
        loader.load(pr);


        t2 = std::chrono::high_resolution_clock::now();
        time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
        std::cout << "MMap rank " << rank << " elapsed time: " << time_span.count() << "s." << std::endl;

        std::cout << rank << " file partition: " << loader.getMMapRange() << std::endl;



        /////////////// now process the file using version with master.
        // do some work using openmp  // version without master
        t1 = std::chrono::high_resolution_clock::now();
        buildIndex<ComputeType>(loader, index, rank, nthreads, chunkSize);

        t2 = std::chrono::high_resolution_clock::now();
        time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
        std::cout << "Compute rank " << rank << " elapsed time: " << time_span.count() << "s." << std::endl;



        /////////////// clean up
        t1 = std::chrono::high_resolution_clock::now();
        loader.unload();
        t2 = std::chrono::high_resolution_clock::now();
        time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
        std::cout << "Unload rank " << rank << " elapsed time: " << time_span.count() << "s." << std::endl;


      }  // scope to ensure file loader is destroyed.

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


          {
        IndexType index(comm, groupSize);
          //index.setLooupAnswerCallback(std::function<void(std::pair<KmerType, std::vector<KmerIndexType> >&)>(&callback));


        /////////////// start processing.  enclosing with braces to make sure loader is destroyed before MPI finalize.
        RunTask<KmerIndexComputeType> t;
        t(filename, index, comm, nthreads, chunkSize);

        printf("MPI number of entries in index for rank %d is %lu\n", id, index.local_size());

        MPI_Barrier(comm);


        //// query:  use the same file as input.  walk through and generate kmers as before.  send query

          }



  //////////////  clean up MPI.
  MPI_Finalize();

  return 0;
}
