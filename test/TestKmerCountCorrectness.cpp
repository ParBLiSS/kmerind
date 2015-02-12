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
#include "index/distributed_map.hpp"
#include "io/CommunicationLayer.hpp"
/*
 * TYPE DEFINITIONS
 */

//Code snippet to execute bash command and return the result
std::string exec(char* cmd) {
  FILE* pipe = popen(cmd, "r");
  if (!pipe) return "ERROR";
  char buffer[8];
  std::string result = "";
  while(!feof(pipe)) {
    if(fgets(buffer, 8, pipe) != NULL)
      result += buffer;
  }
  pclose(pipe);
  return result;
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
    nthreads = omp_get_max_threads();
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

  std::string filename("/home/chirag/Documents/GRA/BigData/Code/bliss/test/data/natural.fastq");
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
      printf("ERROR: The MPI Library Does not have thread support.\n");
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  } else {
    MPI_Init(&argc, &argv);
  }

  MPI_Comm comm = MPI_COMM_WORLD;

  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &rank);


  if (rank == 0)
    std::cout << "USE_MPI is set" << std::endl;
#else
  static_assert(false, "MPI used although compilation is not set to use MPI");
#endif

  {
  // initialize index
  printf("***** initializing index.\n");
  typedef typename bliss::index::KmerCountIndex<21, DNA, bliss::io::FASTQ, true> KmerIndexType;

  KmerIndexType kmer_index(comm, nprocs, nthreads);

  // start processing.  enclosing with braces to make sure loader is destroyed before MPI finalize.
  printf("***** building index first pass.\n");

  kmer_index.build(filename, nthreads, chunkSize);

  std::cout << "COUNT " << rank << " Index Building 1 for " << filename << " using " << nthreads << " threads, index size " << kmer_index.local_size() << std::endl;
  std::cout << std::flush;

  fprintf(stderr, "COUNT %d index built pass 1 with index size: %ld\n", rank, kmer_index.local_size());

  auto& localIndex = kmer_index.getLocalIndex();

  //Run Jellyfish to build its index
  std::string jellyFishExec = "/home/chirag/Documents/GRA/BigData/Code/JelleyFish/jellyfish-2.2.0/bin/jellyfish";

  std::string jellyFishRunCommand = jellyFishExec + " count -m 21 -s 100000 -o jellyFishoutput -C /home/chirag/Documents/GRA/BigData/Code/bliss/test/data/natural.fastq"; 
  exec(&jellyFishRunCommand[0]); 

  //Iterate over all the kmers
  //Get and compare the frequency with Jellyfish's output
  for (auto iter=localIndex.cbegin(); iter!=localIndex.cend(); iter++)
  {
    //std::cout << (iter->first).toAlphabets() << " ";
    std::string jellyFishQueryCommand = jellyFishExec + " query jellyFishoutput " + (iter->first).toAlphabets() + " |  sed 's/[^0-9]*//g'";
    std::string jellyFishResult = exec(&jellyFishQueryCommand[0]);
    //std::cout << iter->second << ", " << std::stoi(jellyFishResult) << "\n";
    assert(iter->second == (unsigned)std::stoi(jellyFishResult));

  }



  //TODO: Uncomment this part after merging with Tony's code
  //kmer_index.finalize();
  MPI_Barrier(comm);

  }


  //////////////  clean up MPI.
  MPI_Finalize();

  DEBUGF("M Rank %d called MPI_Finalize", rank);

  std::cout << std::flush;

  return 0;
}
