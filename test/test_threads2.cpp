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


typedef bliss::index::kmer_index_element<uint64_t, float> kmer_struct_type;
typedef  bliss::io::MPISendBuffer<kmer_struct_type, false> buffer_type;



// can't use auto keyword.  declare and initialize in class declaration
// then "define" but not initialize outside class declaration, again.
template<typename ENCODING, typename Iterator, typename TO, int K>
constexpr std::array<typename ENCODING::value_type, ENCODING::size> bliss::index::generate_qual<ENCODING, Iterator, TO, K>::lut;


#define K 21


typedef bliss::index::generate_kmer<DNA, char*, kmer_struct_type, K> kmer_op_type;
typedef bliss::iterator::buffered_transform_iterator<kmer_op_type, char*> read_iter_type;


typedef bliss::index::generate_qual<bliss::index::SangerToLogProbCorrect<double>, char*, double, K> qual_op_type;
qual_op_type qual_op;
typedef bliss::iterator::buffered_transform_iterator<qual_op_type, char*> qual_iter_type;



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


// this can become a 1 to n transformer???
void compute(bliss::iterator::fastq_sequence<char*> &read, int nprocs, int rank, int nthreads, int pid, int j, std::vector< buffer_type > &buffers ) {

  kmer_op_type kmer_op(read.id);
  read_iter_type start(read.seq, kmer_op);
  read_iter_type end(read.seq_end, kmer_op);

  qual_iter_type qstart(read.qual, qual_op);
  qual_iter_type qend(read.qual_end, qual_op);


  std::pair<kmer_struct_type::kmer_type, kmer_struct_type> index_kmer;

  uint64_t kmerCount = 0;

  int i = -1;
  // NOTE: need to call *start to actually evaluate.  question is whether ++ should be doing computation.
  for (; (start != end) && (qstart != qend); ++start, ++qstart)
  {
    ++i;

    if (i < (K - 1))
      continue;
//        kmers.push_back(*start);
//    kmer ^= *start;
//    qual = *qstart - qual;
    index_kmer = *start;
    index_kmer.second.qual = *qstart;


    // some debugging output
   // printf("kmer send to %lx, key %lx, pos %d, qual %f\n", index_kmer.first, index_kmer.second.kmer, index_kmer.second.id.components.pos, index_kmer.second.qual);

    if (fabs(index_kmer.second.qual) > std::numeric_limits<typename kmer_struct_type::qual_type>::epsilon() ) {
      // sending the kmer.
//      buffers[index_kmer.first % (nprocs * nthreads)].buffer(index_kmer.second);
      buffers[index_kmer.first % nprocs].buffer(index_kmer.second);
//      printf("kmer send to %lx, key %lx, pos %d, qual %f\n", index_kmer.first, index_kmer.second.kmer, index_kmer.second.id.components.pos, index_kmer.second.qual);
      ++kmerCount;
    } else {
//      printf("BAD kmer quality.  key %lx, pos %d, qual %f\n", index_kmer.second.kmer, index_kmer.second.id.components.pos, index_kmer.second.qual);
    }
  }

}


typedef std::unordered_multimap<kmer_struct_type::kmer_type, kmer_struct_type> kmer_map_type;
void networkread(MPI_Comm comm, const int nprocs, const int rank, const size_t buf_size, const int senders, kmer_map_type &kmers) {

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
  int usleep_duration = (1000 + nprocs - 1) / nprocs;

  while (n_senders > 0) {
    // TODO:  have this thread handle all network IO.


    // probe for a message.  if empty message, then don't need to listen on that node anymore
    //printf("probing...\n");
    // NOTE: MPI_Probe does busy-wait (polling) - high CPU utilization.  reduce with --mca mpi_yield_when_idle 1
    // NOTE: that was still kind of slow.
    // NOTE: using Iprobe with usleep is even less CPU utilization.
    //        length between probing needs to be tuned.
    while (!hasMessage) {
      MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, comm, &hasMessage, &status);
      usleep(usleep_duration);
    }
    hasMessage = 0;

    src = status.MPI_SOURCE;
    tag = status.MPI_TAG;
    MPI_Get_count(&status, MPI_UNSIGNED_CHAR, &received);   // length is always > 0?

    if (tag ==  buffer_type::END_TAG) {
      // end of messaging.
      //printf("RECV %d receiving END signal %d from %d\n", rank, received, src); fflush(stdout);
      MPI_Recv(reinterpret_cast<unsigned char*>(array), received, MPI_UNSIGNED_CHAR, src, tag, comm, MPI_STATUS_IGNORE);
      --n_senders;
      continue;
    }

    //printf("RECV %d receiving %d bytes from %d.%d\n", rank, received, src, tag - 1); fflush(stdout);

    MPI_Recv(reinterpret_cast<unsigned char*>(array), received, MPI_UNSIGNED_CHAR, src, tag, comm, MPI_STATUS_IGNORE);

    // if message size is 0, then no work to do.
    if (received == 0)
      continue;

    assert(received % sizeof(kmer_struct_type) == 0);  // no partial messages.
    count = received / sizeof(kmer_struct_type);
    assert(count > 0);
    //printf("RECV+ %d receiving %d bytes or %d records from %d.%d\n", rank, received, count, src, tag - 1); fflush(stdout);

    // TODO:  if we change to this one thread handles all MPI comm,
    //  then need to have another thread handle the insert.

    // now process the array
    //TODO:  DEBUG: temp comment out.
    for (int i = 0; i < count; ++i) {
      kmers.insert(kmer_map_type::value_type(array[i].kmer, array[i]));
    }



//    memset(array, 0, capacity * sizeof(kmer_struct_type));
  }

  delete [] array;

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
  int provided;

#ifdef USE_OPENMP
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

  if (provided < MPI_THREAD_MULTIPLE) {
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

  int nthreads = 1;
#if defined(USE_OPENMP)
  nthreads = omp_get_max_threads();
#endif

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
    assert(nthreads >= 3);
    fprintf(stderr, "rank %d MAX THRADS = %d\n", rank, nthreads);
    omp_set_nested(1);
    omp_set_dynamic(0);
#endif

    int senders = 0;
    MPI_Allreduce(&nthreads, &senders, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

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
    } // omp network read section

#pragma omp section
    {  // compute threads

//      std::vector<  buffer_type > buffers;
//      for (int i = 0; i < nprocs; ++i) {
//        // over provision by the number of threads as well.  (this is okay for smaller number of procs)
//        for (int j = 0; j < nthreads; ++j) {
//          buffers.push_back(std::move( buffer_type(MPI_COMM_WORLD, i, 8192*1024)));
//        }
//      }
//

      // if atEnd, done.
      bool atEnd = false;
      int i = 0;

      // VERSION 2.  uses the fastq iterator as the queue itself, instead of master/slave.
      //   at this point, no strong difference.
#pragma omp parallel num_threads(nthreads)
      {
              std::vector<  buffer_type > buffers;
              for (int j = 0; j < nprocs; ++j) {
                // over provision by the number of threads as well.  (this is okay for smaller number of procs)
                  buffers.push_back(std::move( buffer_type(MPI_COMM_WORLD, j, 8192*1024)));
              }




        bliss::iterator::fastq_sequence<char*> read;
        bool hasData = false;
        int li = 0;
        int tid2 = omp_get_thread_num();

        do {

          // single thread at a time getting data.
#pragma omp critical
          {
            // get data, and move to next for the other threads
            if (!(atEnd = (fastq_start == fastq_end)))
            {
              read = *fastq_start;
              hasData = true;
              li = i;
              ++fastq_start;
              ++i;
            }
          }

          // now do computation.
          if (hasData) {
            compute(read, nprocs, rank, nthreads, tid2, li, buffers);

            if (li % 100000 == 0)
              printf("rank %d thread %d: %d\n", rank, tid2, li);
          }
        } while (!atEnd);

#pragma omp barrier

        // send the last part out.
        for (int j = 0; j < nprocs; ++j) {
          buffers[(j+ rank) % (nprocs)].flush();
        }
      }


      t2 = std::chrono::high_resolution_clock::now();
      time_span = std::chrono::duration_cast<std::chrono::duration<double>>(
          t2 - t1);
      INFO("computation rank " << rank << " elapsed time: " << time_span.count() << "s.");


//      // send the last part out.
//      for (int i = 0; i < nprocs * nthreads; ++i) {
//        buffers[(i+ rank) % (nprocs * nthreads)].flush();
//      }
//      printf("compute completed\n");

    } // omp compute section


    }  // outer parallel sections




    //fileio_t.join();
//    networkwrite_t.join();
    //networkread_t.join();

  }  // scope to ensure file loader is destroyed.


#ifdef USE_MPI
    MPI_Finalize();
#endif

  return 0;
}
