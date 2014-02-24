/**
 *
 * test generating the indices
 *
 * 1,2,3,4 MPI procs, each reading a different part.
 *
 * step 1.  parse fastq into reads
 * setp 2.  for each read, iterate and generate kmers (and reverse complement)
 * step 3.  for each read, concurrently generate quality score
 * step 4.  bin the kmers for redistribution
 * step 5.  MPI send kmers.  (threading?)
 *
 */

#include <string>     // for std::string
#include <sys/stat.h>  // for stat
#include <cassert>
#include <iostream>
#include <chrono>

#include "mpi.h"

#include "utils/logging.h"
#include "config.hpp"

#include "iterators/range.hpp"
#include "io/file_loader.hpp"
#include "io/fastq_loader.hpp"





template<typename T>
struct generate kmers




int main(int argc, char* argv[]) {
  LOG_INIT();


  std::string filename("/home/tpan/src/bliss/test/test.fastq");

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
  }

#ifdef USE_MPI
  // broadcast to all
  MPI_Bcast(&file_size, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
#endif

  if (rank == nprocs - 1)
    fprintf(stderr, "file size is %ld\n", file_size);

  /////////////// now try to open the file


  // first generate an approximate partition.
  bliss::io::file_loader::range_type r = bliss::io::file_loader::range_type::block_partition(file_size, nprocs, rank);
  std::cout << rank << " equipart: " << r << std::endl;
  std::cout << rank << " test block aligned: " << r.align_to_page(sysconf(_SC_PAGE_SIZE)) << std::endl;

  // now open the file and begin reading
  std::chrono::high_resolution_clock::time_point t1, t2;
  std::chrono::duration<double> time_span1;
  std::chrono::duration<double> time_span2;
  std::chrono::duration<double> time_span3;
  std::chrono::duration<double> time_span;


  size_t len;

  {
    t1 = std::chrono::high_resolution_clock::now();
    // now open the file
    bliss::io::fastq_loader loader(filename, r, file_size);
    t2 = std::chrono::high_resolution_clock::now();
    time_span3 = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    INFO("MMap " << "elapsed time: " << time_span3.count() << "s.");

    std::cout << rank << " record adjusted " << loader.getRange() << std::endl;

    len = loader.getRange().end - loader.getRange().start;

//    // test parsing the sequences.
    t1 = std::chrono::high_resolution_clock::now();
    loader.compute_sequence_positions();
    t2 = std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1) + time_span3;
    INFO("get reads " << "elapsed time: " << time_span.count() << "s.");


    // test assign seq ids
    t1 = std::chrono::high_resolution_clock::now();
    loader.assign_sequence_ids();
    t2 = std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1) + time_span3;
    INFO("assign Ids " << "elapsed time: " << time_span.count() << "s.");


    // transform.

  }



#ifdef USE_MPI
  MPI_Finalize();

#endif


  return 0;


}
