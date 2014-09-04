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


#include "utils/logging.h"

#include "common/alphabets.hpp"


#include "index/KmerIndex.hpp"
/*
 * TYPE DEFINITIONS
 */


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
  printf("NOT compiled with OPENMP")
#endif

  int chunkSize = sysconf(_SC_PAGE_SIZE);
  if (argc > 2)
  {
    chunkSize = atoi(argv[2]);
  }

  std::string filename("/home/tpan/src/bliss/test/data/test.fastq");
  //std::string filename("/mnt/data/1000genome/HG00096/sequence_read/SRR077487_1.filt.fastq");
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

  // initialize index
  bliss::index::KmerPositionIndex<21, DNA> kmer_index(comm, groupSize);

  // start processing.  enclosing with braces to make sure loader is destroyed before MPI finalize.
  kmer_index.build(filename, nthreads);

      //// query:  use the same file as input.  walk through and generate kmers as before.  send query


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
