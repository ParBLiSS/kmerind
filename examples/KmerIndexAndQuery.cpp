/**
 *
 * Copyright (c) 2014 Georgia Institute of Technology
 *
 * TODO add Licence
 */

//! [example]


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

/**
 *
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char** argv) {

  int chunkSize = sysconf(_SC_PAGE_SIZE);

  std::string filename("/home/tpan/src/bliss/test/data/test.fastq");

  } else {
    MPI_Init(&argc, &argv);
  }

  MPI_Comm comm = MPI_COMM_WORLD;

  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &rank);

  // initialize index
  using IndexType =  bliss::index::KmerCountIndex<21, bliss::common::DNA, bliss::io::FASTQ, true>;
  IndexType kmer_index(comm, nprocs, nthreads);

  //Build the index
  kmer_index.build(filename, nthreads, chunkSize);

  //Finalize the process
  kmer_index.finalize();
  MPI_Barrier(comm);

  //////////////  clean up MPI.
  MPI_Finalize();
  return 0;
}
