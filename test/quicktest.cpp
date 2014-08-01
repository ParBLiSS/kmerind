/**
 * @file    quicktest.cpp
 * @ingroup
 * @author  tpan
 * @brief
 * @details
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */


#include "mpi.h"
#include <sys/stat.h>   // block size.
#include <iostream>

size_t getFileSize(const std::string& filename) throw (bliss::io::IOException) {
  struct stat64 filestat;

  // get the file state
  int ret = stat64(filename.c_str(), &filestat);

  // handle any error
  if (ret < 0) {
    std::cerr << "ERROR in file size detection: ["  << filename << "] error " << std::endl;
  }

  // return file size.
  return static_cast<size_t>(filestat.st_size);
}


int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);

  MPI_Comm comm = MPI_COMM_WORLD;

  int groupSize;
  int id;

  MPI_Comm_size(comm, &groupSize);
  MPI_Comm_rank(comm, &id);

  std::cout << id <<  " file size: " << getFileSize("/mnt/data/1000genome/HG00096/sequence_read/SRR077487.filt.fastq.gz") << std::endl;

  MPI_Barrier(comm);
  MPI_Finalize();
}
