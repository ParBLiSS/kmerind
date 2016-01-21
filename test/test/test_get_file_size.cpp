/*
 * Copyright 2015 Georgia Institute of Technology
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/**
 * @file    quicktest.cpp
 * @ingroup
 * @author  Tony Pan <tpan7@gatech.edu>
 * @brief
 * @details

 */


#include "mpi.h"

#include <sys/stat.h>   // block size.
#include <iostream>
#include "io/io_exception.hpp"
#include "utils/logging.h"


size_t getFileSize(const std::string& filename) throw (bliss::io::IOException) {
  struct stat64 filestat;

  // get the file state
  int ret = stat64(filename.c_str(), &filestat);

  // handle any error
  if (ret < 0) {
    BL_ERROR( "ERROR in file size detection: ["  << filename << "] error " );
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

  BL_INFO( id <<  " file size: " << getFileSize("/mnt/data/1000genome/HG00096/sequence_read/SRR077487.filt.fastq.gz") );

  MPI_Barrier(comm);
  MPI_Finalize();
}
