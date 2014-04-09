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

#include "mpi.h"

#include <iostream>
//#include <thread>
#include <vector>
#include <unordered_map>
#include <chrono>

#include <unistd.h>
#include <string.h>
#include <cstdio>
#include <cmath>

#include "omp.h"

#include "utils/logging.h"
#include "utils/constexpr_array.hpp"
#include "config.hpp"

#include "iterators/range.hpp"
#include "io/fastq_loader.hpp"

#include "common/alphabets.hpp"
#include "common/AlphabetTraits.hpp"
#include "iterators/buffered_transform_iterator.hpp"

#include "io/MPISendBuffer.hpp"
#include "index/kmer_index_element.hpp"
#include "index/kmer_functors.hpp"


// index type
typedef bliss::index::KmerIndexElement<uint64_t, float> kmer_struct_type;

// MPI buffer type
typedef bliss::io::MPISendBuffer<kmer_struct_type, true> buffer_type;

#define K 21

// kmer generator type
typedef bliss::index::generate_kmer<DNA, char*, kmer_struct_type, K> kmer_op_type;

// can't use auto keyword.  declare and initialize in class declaration
// then "define" but not initialize outside class declaration, again.
template<typename ENCODING, typename Iterator, typename TO, int K>
constexpr std::array<typename ENCODING::value_type, ENCODING::size> bliss::index::generate_qual<
    ENCODING, Iterator, TO, K>::lut;


// quality score generator type
typedef bliss::index::generate_qual<bliss::index::SangerToLogProbCorrect<float>, char*, float, K> qual_op_type;


//void fileio() {
//  // open the file and create a fastq iterator
//    // each dereference has pointers.  choices:
//      // a. directly use.  mmap may be jumping around
//      // b. preload in fastqloader.  jumping around in memory
//      // c. get the pointers and copy the data locally.
//
//  // for each compute thread, accessing next element on file iterator returns a preloaded
//  // data object.  with locking to update pointer.
//
//  // return data object pointing to the local copy.
//  printf("file io\n");
//
//}
//
//

// use a vector of MPIBuffers to manage process'
//


typedef std::unordered_multimap<kmer_struct_type::kmer_type, kmer_struct_type> kmer_map_type;
void networkread(MPI_Comm comm, const int nprocs, const int rank,
                 const size_t buf_size, const int senders, kmer_map_type &kmers)
{

  // track how many active senders remain.
  int n_senders = senders;  // hack.  each proc has multiple threads.
  MPI_Status status;

  size_t capacity = buf_size / sizeof(kmer_struct_type);
  kmer_struct_type *array = new kmer_struct_type[capacity];
  //printf("created temp storage for read, capacity = %ld\n", capacity); fflush(stdout);
  memset(array, 0, capacity * sizeof(kmer_struct_type));
  int received = 0;
  int count = 0;
  int src;
  int tag;

  int hasMessage = 0;
  int usleep_duration = (1000 + nprocs - 1) / nprocs; // more procs, more frequent receive.

  while (n_senders > 0)
  {
    // probe for a message.  if empty message, then don't need to listen on that node anymore
    //printf("probing...\n");
    // NOTE: MPI_Probe does busy-wait (polling) - high CPU utilization.  reduce with --mca mpi_yield_when_idle 1
    // NOTE: that was still kind of slow.
    // NOTE: using Iprobe with usleep is even less CPU utilization.
    //        length between probing needs to be tuned.
    while (!hasMessage)
    {
      MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, comm, &hasMessage, &status);
      usleep(usleep_duration);
    }
    hasMessage = 0;

    src = status.MPI_SOURCE;
    tag = status.MPI_TAG;
    MPI_Get_count(&status, MPI_UNSIGNED_CHAR, &received); // length is always > 0?

    if (tag == buffer_type::END_TAG)
    {
      // end of messaging.
      //printf("RECV %d receiving END signal %d from %d\n", rank, received, src);
      fflush(stdout);
      MPI_Recv(reinterpret_cast<unsigned char*>(array), received,
               MPI_UNSIGNED_CHAR, src, tag, comm, MPI_STATUS_IGNORE);
      --n_senders;
      continue;
    }

    //printf("RECV %d receiving %d bytes from %d.%d\n", rank, received, src, tag - 1); fflush(stdout);

    // TODO: FIX segfault here at iter 2.
    MPI_Recv(reinterpret_cast<unsigned char*>(array), received,
             MPI_UNSIGNED_CHAR, src, tag, comm, MPI_STATUS_IGNORE);

    // if message size is 0, then no work to do.
    if (received == 0)
      continue;

    assert(received % sizeof(kmer_struct_type) == 0);  // no partial messages.
    count = received / sizeof(kmer_struct_type);
    assert(count > 0);
    //printf("RECV+ %d receiving %d bytes or %d records from %d.%d\n", rank, received, count, src, tag - 1); fflush(stdout);

    // now process the array
    //DEBUG: temp comment out.
    //    for (int i = 0; i < count; ++i) {
    //      kmers.insert(kmer_map_type::value_type(array[i].kmer, array[i]));
    //    }

    //    memset(array, 0, capacity * sizeof(kmer_struct_type));
  }

  delete[] array;

  //printf("network read done\n");
  fflush(stdout);
}

int main(int argc, char** argv)
{

  LOG_INIT();

  // TODO: need to fail if there is no OpenMP support or thread count is less than 2.

  std::string filename("/home/tpan/src/bliss/test/data/test.fastq");
  //  std::string filename("/mnt/data/1000genome/HG00096/sequence_read/SRR077487_1.filt.fastq");

  if (argc > 1)
  {
    filename.assign(argv[1]);
  }

  int rank = 0, nprocs = 1;
  int tid = 0, nthreads = 1;
#ifdef USE_MPI
  // initialize MPI
  int provided;

#ifdef USE_OPENMP
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

  nthreads = omp_get_max_threads();
  if (provided < MPI_THREAD_MULTIPLE)
  {
    printf("ERROR: The MPI Library Does not have full thread support.\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
#else
  MPI_Init(&argc, &argv);
#endif

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
  if (nprocs > 1)
  {
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
      bliss::io::file_loader::range_type::block_partition(nprocs, rank, 0,
                                                          file_size);
  std::cout << rank << " equipart: " << r << std::endl;
  std::cout << rank << " test block aligned: "
            << r.align_to_page(sysconf(_SC_PAGE_SIZE)) << std::endl;

  // now open the file and begin reading
  std::chrono::high_resolution_clock::time_point t1, t2;
  std::chrono::duration<double> time_span;

  {
    t1 = std::chrono::high_resolution_clock::now();
    // now open the file
    bliss::io::fastq_loader loader(filename, r, file_size);
    //    bliss::io::file_loader loader(filename, file_size);
    t2 = std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(
        t2 - t1);
    INFO("MMap rank " << rank << " elapsed time: " << time_span.count() << "s.");

    std::cout << rank << " record adjusted " << loader.getRange() << std::endl;

    std::unordered_multimap<kmer_struct_type::kmer_type, kmer_struct_type> index;

    //  std::thread fileio_t(fileio);
    //    std::thread networkwrite_t(networkwrite);
    //std::thread networkread_t(networkread);

    // do some work using openmp
    t1 = std::chrono::high_resolution_clock::now();
    bliss::io::fastq_loader::iterator fastq_start = loader.begin();
    bliss::io::fastq_loader::iterator fastq_end = loader.end();

#ifdef USE_OPENMP
    fprintf(stderr, "rank %d MAX THRADS = %d\n", rank, nthreads);
    omp_set_nested(1);

#endif

    int senders = 0;
    MPI_Allreduce(&nthreads, &senders, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    int buf_size = nthreads * nprocs;

    bliss::iterator::fastq_sequence<char*> read;
    int i = 0;

#pragma omp parallel sections num_threads(2) private(t2, time_span)
    {

#pragma omp section
      {
        networkread(MPI_COMM_WORLD, nprocs, rank, 8192 * 1024, senders, index);
        t2 = std::chrono::high_resolution_clock::now();
        time_span =
            std::chrono::duration_cast<std::chrono::duration<double>>(
                t2 - t1);
        INFO(
            "read receive rank " << rank << " elapsed time: " << time_span.count() << "s, " << index.size() << "records.");

      } // network read thread

#pragma omp section
      {  // compute threads

        std::vector<buffer_type> buffers;
        for (int j = 0; j < buf_size; ++j)
        {

          //printf("insert buffer %d for thread %d to rank %d\n", j, j / nprocs, j % nprocs);
          buffers.push_back(
              std::move(buffer_type(MPI_COMM_WORLD, j % nprocs, 8192 * 1024)));
        }


        std::vector<size_t> counts;
        for (int j = 0; j < nthreads; ++j) {
          counts.push_back(0);
        }

#pragma omp parallel num_threads(nthreads) private(tid)
        {
          tid = omp_get_thread_num();

          int li = 0;

#pragma omp single nowait
          {
            for (; fastq_start != fastq_end; ++fastq_start, ++i)
            {
              // get data, and move to next for the other threads

              // first get read
              read = *fastq_start;
              li = i;
              //
              //            ++fastq_start;
              //            ++i;
#pragma omp task firstprivate(read, li, tid)
              {
                // copy read.  not doing that right now.
                int tid2 = omp_get_thread_num();
                // then compute
                compute(read, nprocs, rank, tid2, li, buffers, counts);

                if (li % 100000 == 0)
                  printf("rank %d thread %d: %d\n", rank, tid2, li);
                // release resource
              }

              // next iteration will check to see if the iterator is at end,
              // and if not, get and compute.
            }
          } // omp single
#pragma omp taskwait

        } // compute threads parallel

        t2 = std::chrono::high_resolution_clock::now();
        time_span =
            std::chrono::duration_cast<std::chrono::duration<double>>(
                t2 - t1);
        INFO(
            "computation rank " << rank << " elapsed time: " << time_span.count() << "s.");


        auto t3 = std::chrono::high_resolution_clock::now();
        // send the last part out.
        for (int j = 0; j < buf_size; ++j)
        {
          buffers[(j + rank) % buf_size].flush();
        }
        auto t4  = std::chrono::high_resolution_clock::now();
        time_span =
            std::chrono::duration_cast<std::chrono::duration<double>>(
                t4 - t3);
        INFO(
            "flush rank " << rank << " thread " << tid << " elapsed time: " << time_span.count() << "s.");


        for (size_t j = 0; j < counts.size(); ++j) {
          printf("rank %d COUNTS by thread %lud: %lud\n", rank, j, counts[j]);
        }

      }  // compute threads section

    } // outer parallel sections

    //fileio_t.join();
    //    networkwrite_t.join();
    //networkread_t.join();
  }

#ifdef USE_MPI
  MPI_Finalize();
#endif

  return 0;
}
