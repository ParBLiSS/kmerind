/**
 * @file    TestKmerCountCorrectness.cpp
 * @ingroup
 * @author  cjain7
 * @brief   Tests the kmer frequency count output of the library against the jellyFish's histogram
 * @details
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#include "config.hpp"

#include <functional>
#include <fstream>

#include "utils/logging.h"

#include "common/alphabets.hpp"


#include "index/kmer_index.hpp"
#include "index/distributed_map.hpp"
#include "io/communication_layer.hpp"
/*
 * TYPE DEFINITIONS
 */

using namespace std::placeholders;   // for _1, _2, _3

//Code snippet to execute bash command and return the result
/*std::string exec(char* cmd) {*/
  //FILE* pipe = popen(cmd, "r");
  //if (!pipe) return "ERROR";
  //char buffer[8];
  //std::string result = "";
  //while(!feof(pipe)) {
    //if(fgets(buffer, 8, pipe) != NULL)
      //result += buffer;
  //}
  //pclose(pipe);
  //return result;
/*}*/

/**
 * @param filename  File with space seperated frequency and counts
 */
std::map<int, int> buildHistogramfromFile(std::string filename)
{
  std::map<int, int> histoMap;

  std::ifstream infile(filename);
  int freq, count;

  while (infile >> freq >> count)
  {
    histoMap.insert(std::make_pair(freq, count));
  }

  return histoMap;
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

  std::string filename;
  filename.assign(PROJ_SRC_DIR);
  filename.append("/test/data/natural.fastq");

  std::string solutionFileName;
  solutionFileName.assign(PROJ_SRC_DIR);
  solutionFileName.append("/test/data/natural.fastq.jellyFish.histo");

  std::cout << "DEBUGGING : " << filename << "\n";

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
  typedef bliss::index::KmerCountIndex<21, bliss::common::DNA, bliss::io::FASTQ, true> KmerIndexType;

  size_t result = 0, entries = 0;
  KmerIndexType kmer_index(comm, nprocs,
                                     std::bind(KmerIndexType::defaultReceiveAnswerCallback,
                                               _1, _2, nthreads, std::ref(result), std::ref(entries)), 
                                     nthreads);

  // start processing.  enclosing with braces to make sure loader is destroyed before MPI finalize.
  printf("***** building index first pass.\n");

  kmer_index.build(filename, nthreads, chunkSize);

  std::cout << "COUNT " << rank << " Index Building 1 for " << filename << " using " << nthreads << " threads, index size " << kmer_index.local_size() << std::endl;
  std::cout << std::flush;

  fprintf(stderr, "COUNT %d index built pass 1 with index size: %ld\n", rank, kmer_index.local_size());

  //Make sure every proc/ thread is done 
  kmer_index.finalize();

  auto& localIndex = kmer_index.getLocalIndex();

  //If n is the highest frequency, vector should contain n entries
  //ith element denotes the count of kmers with (i+1) frequency
  auto overallHistogram = localIndex.countHistogram();

  //Get a map with frequency as key
  auto solutionHistogram = buildHistogramfromFile(solutionFileName);

  //Start checking
  for(auto& item: solutionHistogram)
  {
    assert(overallHistogram[item.first -1]  == item.second);
    cout << item.first << ":" << item.second << " matches \n";
  }
  
  MPI_Barrier(comm);

  }


  //////////////  clean up MPI.
  MPI_Finalize();

  DEBUGF("M Rank %d called MPI_Finalize", rank);

  std::cout << std::flush;

  return 0;
}
