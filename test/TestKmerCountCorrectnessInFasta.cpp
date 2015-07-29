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
#include "bliss-config.hpp"

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

/**
 * @brief returns the histogram of kmer frequencies without using Bliss
 * @details If v is the vector returned, v[i] denotes the count of unique
 *          kmers which have frequence i. v[0] = 0
 */
std::vector<int> sequentialbuildHistogram(std::string filename, int kmerLength)
{
  //Assuming the given FASTA file is well-structured & (all reads are just a single line)
  std::ifstream infile(filename);
  std::string ignore, read;
  std::vector<std::string> kmerVectorWithDuplicates;

  //Process the file
  while (infile >> ignore >> read)
  {
    //Process this read
    for(unsigned int i= 0; i < read.size() - (unsigned)kmerLength + 1; i++)
      kmerVectorWithDuplicates.push_back(read.substr(i, kmerLength));
  }

  std::sort(begin(kmerVectorWithDuplicates), end(kmerVectorWithDuplicates));
  std::vector<int> allFrequenciesUniqueKmers;

  for(auto i= begin(kmerVectorWithDuplicates) ; i != end(kmerVectorWithDuplicates);)
  {
    auto r = std::equal_range(i, end(kmerVectorWithDuplicates), *i);
    allFrequenciesUniqueKmers.emplace_back( std::distance(r.first, r.second) );
    i = r.second;
  }

  //Build histogram
  std::sort(begin(allFrequenciesUniqueKmers), end(allFrequenciesUniqueKmers));
  std::vector<int> histogram;

  histogram.emplace_back(0);
  for(auto i= begin(allFrequenciesUniqueKmers) ; i != end(allFrequenciesUniqueKmers);)
  {
    auto r = std::equal_range(i, end(allFrequenciesUniqueKmers), *i);
    histogram.emplace_back( std::distance(r.first, r.second) );
    i = r.second;
  }

  return histogram;
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
  printf("NOT compiled with OPENMP");
#endif

  int chunkSize = sysconf(_SC_PAGE_SIZE);
  if (argc > 2)
  {
    chunkSize = atoi(argv[2]);
  }

  std::string filename;
  filename.assign(PROJ_SRC_DIR);
  filename.append("/test/data/natural.fasta");

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

  const int kmerLength = 21;
  typedef bliss::index::KmerCountIndex<bliss::common::Kmer<kmerLength, bliss::common::DNA>, bliss::io::FASTA> KmerIndexType;

  size_t entries = 0;
  double time = 0;
  KmerIndexType kmer_index(comm, nprocs,
                                     std::bind(KmerIndexType::defaultReceiveAnswerCallback,
                                               _1, _2, nthreads, std::ref(entries), std::ref(time)),
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
  //auto solutionHistogram = buildHistogramfromFile(solutionFileName);
  auto ourOwnHistogram = sequentialbuildHistogram(filename, kmerLength);

  for(auto const& e : ourOwnHistogram) std::cout << e << "; "; std::cout << "\n";
  for(auto const& e : overallHistogram) std::cout << e << "; "; std::cout << "\n";

  //assert(ourOwnHistogram.size() == overallHistogram.size());
  //bool is_equal = std::equal(overallHistogram.begin(), overallHistogram.end(), ourOwnHistogram.begin());
  
  //assert(is_equal == true);

  MPI_Barrier(comm);

  }


  //////////////  clean up MPI.
  MPI_Finalize();

  DEBUGF("M Rank %d called MPI_Finalize", rank);

  std::cout << std::flush;

  return 0;
}
