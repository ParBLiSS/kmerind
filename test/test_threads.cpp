/**
 * @file		test_threads.cpp
 * @ingroup
 * @author	tpan
 * @brief
 * @details
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */

#include <iostream>
#include <thread>

#include "mpi.h"
#include "omp.h"

#include "utils/logging.h"
#include "config.hpp"

#include "iterators/range.hpp"
#include "io/fastq_loader.hpp"
#include "io/fastq_iterator.hpp"


void fileio() {
  // open the file and create a fastq iterator
    // each dereference has pointers.  choices:
      // a. directly use.  mmap may be jumping around
      // b. preload in fastqloader.  jumping around in memory
      // c. get the pointers and copy the data locally.

  // for each compute thread, accessing next element on file iterator returns a preloaded
  // data object.  with locking to update pointer.

  // return data object pointing to the local copy.
  printf("file io\n");

}

void compute(int i, bliss::iterator::fastq_sequence<char*> &read) {
   // do something.

  printf("compute %d, %ld\n", i, read.id);

}

void networkwrite() {
  // instantiate bins (m per proc)

  // atomic add to bin. update bin count

  // if bin count is at limit, MPI send length, then send data
  //
  // final write
    // flush: sends length
      // if length > 0, send data
      // else length = 0, don't send anything.  return to waiting.

  // send termination signal (send empty message).
  printf("network write\n");

}

void networkread() {
  // initiate an array of senders, mark all as 1.

  // iprobe for a message.  if empty message, then don't need to listen on that node anymore

  // if has length,
    // if length > 0, then receive data and do something
  // else return to probe.

  printf("network read\n");
}


int main(int argc, char** argv) {

  LOG_INIT();

  std::string filename("/home/tpan/src/bliss/test/data/test.fastq");
//  std::string filename("/mnt/data/1000genome/HG00096/sequence_read/SRR077487_1.filt.fastq");

  if (argc > 1)
  {
    filename.assign(argv[1]);
  }

  int rank = 0, nprocs = 1;
#ifdef USE_MPI
  // initialize MPI
  MPI_Init(&argc, &argv);

  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == 0)
    std::cout << "USE_MPI is set" << std::endl;
#endif

  // first thread gets the file size.
  uint64_t file_size = 0;
  if (rank == 0)
  {
    struct stat filestat;
    stat(filename.c_str(), &filestat);
    file_size = static_cast<uint64_t>(filestat.st_size);
    fprintf(stderr, "block size is %ld\n", filestat.st_blksize);
    fprintf(stderr, "sysconf block size is %ld\n", sysconf(_SC_PAGE_SIZE));
  }

#ifdef USE_MPI
  if (nprocs > 1) {
    // broadcast to all
    MPI_Bcast(&file_size, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
  }
#endif

  if (rank == nprocs - 1)
  {
    fprintf(stderr, "file size is %ld\n", file_size);
  }
  /////////////// now try to open the file

  // real data:  mmap is better for large files and limited memory.
  //             preloading is better for smaller files and/or large amount of memory.
  // stream processing means data does not need to be buffered in memory - more efficient.

  // file access:  better to work with a few pages at a time, or to work with large block?

  // first generate an approximate partition.
  bliss::io::file_loader::range_type r =
      bliss::io::file_loader::range_type::block_partition(nprocs,
                                                          rank, 0, file_size);
  std::cout << rank << " equipart: " << r << std::endl;
  std::cout << rank << " test block aligned: "
            << r.align_to_page(sysconf(_SC_PAGE_SIZE)) << std::endl;

  // now open the file and begin reading
  std::chrono::high_resolution_clock::time_point t1, t2;
  std::chrono::duration<double> time_span;

  {
    t1 = std::chrono::high_resolution_clock::now();
    // now open the file
    bliss::io::fastq_loader loader(filename, r, file_size, true);
//    bliss::io::file_loader loader(filename, file_size);
    t2 = std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(
        t2 - t1);
    INFO("MMap rank " << rank << " elapsed time: " << time_span.count() << "s.");

    std::cout << rank << " record adjusted " << loader.getRange() << std::endl;

    t1 = std::chrono::high_resolution_clock::now();
    bliss::io::fastq_loader::iterator fastq_start = loader.begin();
    bliss::io::fastq_loader::iterator fastq_end = loader.end();





    std::cout << "starting threads" << std::endl;

  //  std::thread fileio_t(fileio);
    std::thread networkwrite_t(networkwrite);
    std::thread networkread_t(networkread);

    std::cout << "threads started" << std::endl;
    // do some work using openmp

    int threadcount = 3;

    bliss::iterator::fastq_sequence<char*> read;
  #pragma omp parallel shared(fastq_start, fastq_end) private(read) num_threads(threadcount)
    {
      while (fastq_start != fastq_end) {
        // get data, and move to next for the other threads

  #pragma omp task
        {

          // first get read (assume copy already made)
  #pragma omp critical
          {
            read = *fastq_start;

            ++fastq_start;
          }

          // then compute
          compute(omp_get_thread_num(), read);
        }

        // next iteration will check to see if the iterator is at end,
        // and if not, get and compute.
      }
    }
    std::cout << "computation done" << std::endl;


    //fileio_t.join();
    networkwrite_t.join();
    networkread_t.join();
  }
  std::cout << "threads done" << std::endl;

  return 0;
}
