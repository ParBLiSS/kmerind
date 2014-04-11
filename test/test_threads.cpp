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

#include <unistd.h>
#include <string.h>
#include <cstdio>
#include <cmath>

#include <iostream>
//#include <thread>
#include <vector>
#include <unordered_map>
#include <chrono>

#include "omp.h"

#include "config.hpp"

#include "utils/logging.h"
#include "utils/constexpr_array.hpp"
#include "common/base_types.hpp"
#include "common/alphabets.hpp"
#include "common/AlphabetTraits.hpp"
#include "iterators/range.hpp"
#include "iterators/buffered_transform_iterator.hpp"
#include "io/fastq_loader.hpp"
#include "io/MPISendBuffer.hpp"
#include "index/kmer_index_element.hpp"
#include "index/kmer_index_functors.hpp"
#include "index/kmer_index_generation.hpp"


/*
 * TYPE DEFINITIONS
 */
// define kmer index type
typedef bliss::index::KmerSize<21> KmerSize;
typedef uint64_t KmerType;
typedef float QualityType;
typedef bliss::index::KmerIndexElementWithIdAndQuality<KmerSize, KmerType, bliss::io::fastq_sequence_id, QualityType> KmerIndexType;

// define buffer where to put the kmers
constexpr bool thread_safe = false;
typedef bliss::io::MPISendBuffer<KmerIndexType, thread_safe> BufferType;

// define raw data type :  use CharType
typedef CharType* BaseIterType;

// define read type
typedef DNA Alphabet;
typedef bliss::io::fastq_sequence_quality<BaseIterType, Alphabet, float>  SequenceType;


typedef bliss::index::generate_kmer<SequenceType, KmerIndexType> kmer_op_type;

typedef bliss::index::generate_qual<SequenceType, KmerSize, bliss::index::SangerToLogProbCorrect<float>, QualityType> qual_op_type;

typedef std::unordered_multimap<KmerType, KmerIndexType> IndexType;

typedef bliss::io::fastq_loader<Alphabet, QualityType> FileLoaderType;

typedef bliss::index::KmerIndexGeneratorWithQuality<std::vector<BufferType>, kmer_op_type, qual_op_type, bliss::index::XorModulus<KmerType> > ComputeType;


/**
 * receiver for MPI communications.
 */
void networkread(MPI_Comm comm, const int nprocs, const int rank, const size_t buf_size, const int senders, IndexType &kmers) {

  // track how many active senders remain.
  int n_senders = senders;  // hack.  each proc has multiple threads.
  MPI_Status status;

  size_t capacity = buf_size / sizeof(KmerIndexType);
  KmerIndexType *array = new KmerIndexType[capacity];
  //printf("created temp storage for read, capacity = %ld\n", capacity); fflush(stdout);
  memset(array, 0, capacity * sizeof(KmerIndexType));
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

    if (tag ==  BufferType::END_TAG) {
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

    assert(received % sizeof(KmerIndexType) == 0);  // no partial messages.
    count = received / sizeof(KmerIndexType);
    assert(count > 0);
    //printf("RECV+ %d receiving %d bytes or %d records from %d.%d\n", rank, received, count, src, tag - 1); fflush(stdout);

    // TODO:  if we change to this one thread handles all MPI comm,
    //  then need to have another thread handle the insert.

    // now process the array
    //TODO:  DEBUG: temp comment out.
    for (int i = 0; i < count; ++i) {
      kmers.insert(IndexType::value_type(array[i].kmer, array[i]));
    }



//    memset(array, 0, capacity * sizeof(KmerIndexType));
  }

  delete [] array;

}

/**
 * initializes MPI and OpenMP.
 *
 * @param argc
 * @param argv
 * @param comm
 * @param nproc
 * @param rank
 * @param nthreads
 * @param tid
 */
void init(int argc, char** argv, MPI_Comm &comm, int &nprocs, int &rank, int &nthreads) {

#if defined(USE_OPENMP)
  nthreads = omp_get_max_threads();
#endif

#ifdef USE_MPI
  // initialize MPI

#ifdef USE_OPENMP
  int provided;
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
#else
  //TODO:  need to support no MPI.
  fprintf(stderr, "ERROR: need MPI support\n");
  exit(1);
#endif

}

/**
 * finalizes MPI (and OpenMP)
 */
void finalize(MPI_Comm &comm) {
  #ifdef USE_MPI
      MPI_Finalize();
  #endif
}

FileLoaderType openFile(const std::string &filename, const int &nprocs, const int &rank) {
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

  // real data:  mmap is better for large files and limited memory.
  //             preloading is better for smaller files and/or large amount of memory.
  // stream processing means data does not need to be buffered in memory - more efficient.

  // file access:  better to work with a few pages at a time, or to work with large block?

  // first generate an approximate partition.
  FileLoaderType::RangeType r =
      FileLoaderType::RangeType::block_partition(nprocs, rank, 0, file_size);

  FileLoaderType loader = FileLoaderType(filename, r, file_size);

  return loader;
}

// TODO: make these function with/without MPI, with/without OpenMP.
// TODO: 1. make a communicator that encapsulates the nprocs, rank, nthread, tid, etc - this refactors KmerIndexGenerator so there is no dependency on nprocs, rank, and tid.
// TODO: 2. refactor these "compute" methods so that the loop iterations are handled by a ThreadModel.
// TODO: 3. make the communicator short circuit when copying to local, and lock when multithreading.
// TODO: 4. communicator should have send buffer and receive buffer (actual data structure)

/*
 * first pattern for threading - has a master thread
 */
template <typename Compute>
void compute_MPI_OMP_WithMaster(FileLoaderType &loader,
                                int nthreads, IndexType &index, int nprocs, int rank) {
  std::chrono::high_resolution_clock::time_point t1, t2;
  std::chrono::duration<double> time_span;

  // do some work using openmp
  t1 = std::chrono::high_resolution_clock::now();

  ///  get the start and end iterator position.
  FileLoaderType::IteratorType fastq_start = loader.begin();
  FileLoaderType::IteratorType fastq_end = loader.end();

#ifdef USE_OPENMP
  assert(nthreads >= 2);
  fprintf(stderr, "rank %d MAX THRADS = %d\n", rank, nthreads);
  omp_set_nested(1);
  omp_set_dynamic(0);
#endif

  int senders = nthreads;
#ifdef USE_MPI
  senders = 0;
  MPI_Allreduce(&nthreads, &senders, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif



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


      int buf_size = nthreads * nprocs;

      std::vector<BufferType> buffers;
      for (int j = 0; j < buf_size; ++j)
      {

        //printf("insert buffer %d for thread %d to rank %d\n", j, j / nprocs, j % nprocs);
        buffers.push_back(
            std::move(BufferType(MPI_COMM_WORLD, j % nprocs, 8192 * 1024)));
      }


      std::vector<size_t> counts;
      for (int j = 0; j < nthreads; ++j) {
        counts.push_back(0);
      }

#pragma omp parallel num_threads(nthreads)
      {
        Compute op(nprocs, rank, omp_get_thread_num());

        SequenceType read;
        int tid = omp_get_thread_num();

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
              op(read, li, buffers, counts);

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
      INFO("flush rank " << rank << " elapsed time: " << time_span.count() << "s.");


      for (size_t j = 0; j < counts.size(); ++j) {
        printf("rank %d COUNTS by thread %lud: %lud\n", rank, j, counts[j]);
      }

    }  // compute threads section

  } // outer parallel sections

}


/*
 * second pattern for threading
 */
template <typename Compute>
void compute_MPI_OMP_NoMaster(FileLoaderType &loader,
                              int nthreads, IndexType &index, int nprocs, int rank)
{
  std::chrono::high_resolution_clock::time_point t1, t2;
  std::chrono::duration<double> time_span;


  ///  get the start and end iterator position.
  FileLoaderType::IteratorType fastq_start = loader.begin();
  FileLoaderType::IteratorType fastq_end = loader.end();



#ifdef USE_OPENMP
    assert(nthreads >= 2);
    fprintf(stderr, "rank %d MAX THRADS = %d\n", rank, nthreads);
    omp_set_nested(1);
    omp_set_dynamic(0);
#endif

    int senders = nthreads;
#ifdef USE_MPI
    senders = 0;
    MPI_Allreduce(&nthreads, &senders, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif

#pragma omp parallel sections num_threads(2) private(t1, t2, time_span)
    {


#pragma omp section
    {
      t1 = std::chrono::high_resolution_clock::now();
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

//      std::vector<  BufferType > buffers;
//      for (int i = 0; i < nprocs; ++i) {
//        // over provision by the number of threads as well.  (this is okay for smaller number of procs)
//        for (int j = 0; j < nthreads; ++j) {
//          buffers.push_back(std::move( BufferType(MPI_COMM_WORLD, i, 8192*1024)));
//        }
//      }
//
      t1 = std::chrono::high_resolution_clock::now();

      // if atEnd, done.
      bool atEnd = false;
      int i = 0;

      std::vector<size_t> counts;
      for (int j = 0; j < nthreads; ++j) {
        counts.push_back(0);
      }


      // VERSION 2.  uses the fastq iterator as the queue itself, instead of master/slave.
      //   at this point, no strong difference.
#pragma omp parallel num_threads(nthreads)
      {

        Compute op(nprocs, rank, omp_get_thread_num());

        std::vector<  BufferType > buffers;
              for (int j = 0; j < nprocs; ++j) {
                // over provision by the number of threads as well.  (this is okay for smaller number of procs)
                  buffers.push_back(std::move( BufferType(MPI_COMM_WORLD, j, 8192*1024)));
              }




        SequenceType read;
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
            op(read, li, buffers, counts);

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

}

/**
 * source:              source data,  e.g. fastq_loader
 * Compute:             transformation to index
 * ThreadModel:         threading for compute
 * Communicator:        move output of compute to destination.  e.g. MPI send/recv.
 * Destination:         target container, e.g. index hash table
 *
 * compute is single threaded.
 * threadModel determines whether we are using master/slave or demand driven bag of tasks
 * communicator contains buffer, and is responsible for moving data around.
 * Destination is the target data structure.
 *
 * flow is source -> compute under thread model -> buffer via communicator --> send --> recv -> Destination.
 */
template <typename Source, typename Compute, typename ThreadModel, typename Communicator, typename Destination>
void compute(Source src, Destination dest) {



}


/**
 *
 * @param argc
 * @param argv
 * @return
 *
 *
 *
 */
int main(int argc, char** argv) {

  //////////////// init logging
  LOG_INIT();

  //////////////// parse parameters
  std::string filename("/home/tpan/src/bliss/test/data/test.fastq");

  bool WithMaster = true;
  if (argc > 1)
  {
    if (strcmp(argv[1], "0") == 0) {
      WithMaster = false;
    }
  }

//  std::string filename("/mnt/data/1000genome/HG00096/sequence_read/SRR077487_1.filt.fastq");
  if (argc > 2)
  {
    filename.assign(argv[2]);
  }



  //////////////// initialize MPI and openMP
  int rank = 0, nprocs = 1;
  int tid = 0, nthreads = 1;
  MPI_Comm comm;
  init(argc, argv, comm, nprocs, rank, nthreads);


  /////////////// initialize output variables
  IndexType index;


  /////////////// start processing.  enclosing with braces to make sure loader is destroyed before MPI finalize.
  {
    std::chrono::high_resolution_clock::time_point t1, t2;
    std::chrono::duration<double> time_span;

    //////////////// now partition and open the file.
    t1 = std::chrono::high_resolution_clock::now();

    FileLoaderType loader = openFile(filename, nprocs, rank);  // this handle is alive through the entire execution.

    t2 = std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    INFO("MMap rank " << rank << " elapsed time: " << time_span.count() << "s.");

    std::cout << rank << " file partition: " << loader.getRange() << std::endl;



    /////////////// now process the file using version with master.
    if (WithMaster) {
      // do some work using openmp  // version without master
      compute_MPI_OMP_WithMaster<ComputeType>(loader, nthreads, index, nprocs, rank);
    } else {
      /////////////// now process the file using version with master.
      // do some work using openmp  // version without master
      compute_MPI_OMP_NoMaster<ComputeType>(loader, nthreads, index, nprocs, rank);
    }

  }  // scope to ensure file loader is destroyed.


  //////////////  clean up MPI.
  finalize(comm);

  return 0;
}
