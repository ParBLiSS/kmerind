/**
 * @file		test_threads.cpp
 * @ingroup
 * @author	Tony Pan <tpan7@gatech.edu>
 * @brief
 * @details
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */

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
#define OLD_USE_MPI (USE_MPI)
#undef USE_MPI

#include "utils/logging.h"
#include "utils/constexpr_array.hpp"
#include "common/base_types.hpp"
#include "common/alphabets.hpp"
#include "common/alphabet_traits.hpp"
#include "iterators/range.hpp"
#include "iterators/buffered_transform_iterator.hpp"
#include "io/fastq_partition_helper.hpp"
#include "io/MPISendBuffer.hpp"
#include "index/kmer_index_element.hpp"
#include "index/kmer_index_functors.hpp"
#include "index/kmer_index_generation_no_mpi.hpp"


/*
 * TYPE DEFINITIONS
 */
// define kmer index type
typedef bliss::index::KmerSize<21> KmerSize;
typedef uint64_t KmerType;
typedef DNA Alphabet;
typedef bliss::index::KmerIndexElementWithId<KmerSize, KmerType, bliss::io::FASTQSequenceId> KmerIndexType;

// define buffer where to put the kmers
constexpr bool thread_safe = false;
typedef bliss::io::MPISendBuffer<KmerIndexType, thread_safe> BufferType;

// define raw data type :  use CharType
typedef bliss::io::FileLoader<CharType>                       FileLoaderType;
typedef bliss::io::FASTQPartitionHelper<CharType, size_t> PartitionHelperType;
typedef typename FileLoaderType::IteratorType                  BaseIterType;
typedef typename std::iterator_traits<BaseIterType>::value_type BaseValueType;


// define read type
typedef bliss::io::Sequence<BaseIterType>  SequenceType;


typedef bliss::index::generate_kmer<SequenceType, KmerIndexType> kmer_op_type;

typedef std::unordered_multimap<KmerType, KmerIndexType> IndexType;

typedef bliss::io::FASTQParser<BaseIterType>               ParserType;
typedef bliss::io::SequencesIterator<ParserType, BaseIterType>           IteratorType;
typedef bliss::index::KmerIndexGenerator<kmer_op_type, BufferType, bliss::index::XorModulus<KmerType>> ComputeType;









// TODO: make these function with/without MPI, with/without OpenMP.
// TODO: 1. make a communicator that encapsulates the nprocs, rank, nthread, tid, etc - this refactors KmerIndexGenerator so there is no dependency on nprocs, rank, and tid.
// TODO: 2. refactor these "compute" methods so that the loop iterations are handled by a ThreadModel.
// TODO: 3. make the communicator short circuit when copying to local, and lock when multithreading.
// TODO: 4. communicator should have send buffer and receive buffer (actual data structure)

/*
 * first pattern for threading - has a master thread
 */
template <typename Compute>
void compute_OMP_WithMaster(FileLoaderType &loader, PartitionHelperType &ph,
                                const int &nthreads, IndexType &index,
                                const int &nprocs, const int &rank, const int chunkSize) {

  INFO("HAS MASTER");

  // do some work using openmp

  ///  get the start and end iterator position.
  ParserType parser;


#ifdef USE_OPENMP
  INFO("rank " << rank << " MAX THRADS = " << nthreads);
  omp_set_nested(1);
  omp_set_dynamic(0);
#endif

  int senders = nthreads;
  int nt = nthreads;


int buffer_size = 8192*1024;


      INFO("Level 0: compute tid = " << omp_get_thread_num());
      std::chrono::high_resolution_clock::time_point t1, t2;
      std::chrono::duration<double> time_span;

      t1 = std::chrono::high_resolution_clock::now();

      int buf_size = nthreads * nprocs;
      int i = 0;


      std::vector<bliss::index::countType> counts;
      for (int j = 0; j < nthreads; ++j) {
        bliss::index::countType c;
        c.c = 0;
        counts.push_back(c);
      }


#pragma omp parallel num_threads(nthreads) firstprivate(t1, t2, time_span) shared(parser, loader, ph, i, counts, nprocs, rank, nthreads) default(none)
      {
        int tid = omp_get_thread_num();
        INFO("Level 1: compute: thread id = " << tid);

        Compute op(nprocs, rank, nthreads);
        bool copying = false;

#pragma omp single
        {
          INFO("Level 1: MASTER thread id = " << tid);

          // private vars, singly updated.
          SequenceType read;
          int li = 0;

#pragma omp task untied
          {
            BaseIterType begin = nullptr, end = nullptr;
            size_t chunkRead;

            while ((chunkRead = length(loader.getNextChunkAtomic(ph, begin, end, chunkSize, copying))) > 0) {
              li = i;
              ++i;

  #pragma omp task firstprivate(read, li, begin, end, chunkRead, copying) shared(op, rank, parser, counts) default(none)
              {

                IteratorType fastq_start(parser, begin, end, loader.getRange());
                IteratorType fastq_end(parser, end, loader.getRange());

                for (; fastq_start != fastq_end; ++fastq_start)
                {
                  // get data, and move to next for the other threads

                  // first get read
                  read = *fastq_start;
      //            SequenceType::allocCopy(*fastq_start, read);

                //
                //            ++fastq_start;
                //            ++i;
                  // copy read.  not doing that right now.
                  // then compute
                  op(read, li, counts);
    //            SequenceType::deleteCopy(read);


                  // release resource
                }

                if (li % 100000 == 0)
                  INFO("Level 1: rank " << rank << " thread " << omp_get_thread_num() << " processed chunk " << li);

                if (copying) delete [] begin;
              }
              // next iteration will check to see if the iterator is at end,
              // and if not, get and compute.

            }
          }  // untied task
        } // omp single
#pragma omp taskwait

      } // compute threads parallel
      INFO("Level 1: rank " << rank << " thread " << omp_get_thread_num() << " processed " << i);



      t2 = std::chrono::high_resolution_clock::now();
      time_span =
          std::chrono::duration_cast<std::chrono::duration<double>>(
              t2 - t1);
      INFO(
          "Level 0: Computation: rank " << rank << " thread " << omp_get_thread_num() << " elapsed time: " << time_span.count() << "s.");


      for (size_t j = 0; j < counts.size(); ++j) {
        INFO("rank " << rank << " COUNTS by thread "<< j << ": "<< counts[j].c);
      }


}


/*
 * second pattern for threading
 */
template <typename Compute>
void compute_OMP_NoMaster(FileLoaderType &loader, PartitionHelperType &ph,
                              int &nthreads, IndexType &index,
                              const int &nprocs, const int &rank, const int chunkSize)
{


  INFO("NO MASTER");

  ///  get the start and end iterator position.
  ParserType parser;



#ifdef USE_OPENMP
    INFO("rank " << rank << " MAX THRADS = " << nthreads);
    omp_set_nested(1);
    omp_set_dynamic(0);
#endif

    int senders = nthreads;
//    INFOF("senders = %d\n", senders);

      std::chrono::high_resolution_clock::time_point t1, t2;
      std::chrono::duration<double> time_span;

      INFO("Level 0: compute tid = " << omp_get_thread_num());

      t1 = std::chrono::high_resolution_clock::now();

      // if atEnd, done.
//      bool atEnd = false;
      int i = 0;



      std::vector<bliss::index::countType> counts;
      for (int j = 0; j < nthreads; ++j) {
        bliss::index::countType c;
        c.c = 0;
        counts.push_back(c);
      }


      // VERSION 2.  uses the fastq iterator as the queue itself, instead of master/slave.
      //   at this point, no strong difference.
#pragma omp parallel num_threads(nthreads) shared(loader, parser, ph, i, counts, nprocs, nthreads, rank) default(none)
      {
        Compute op(nprocs, rank, nthreads);


//        std::vector<  BufferType > buffers;
//        for (int j = 0; j < nprocs; ++j) {
//          // over provision by the number of threads as well.  (this is okay for smaller number of procs)
//            buffers.push_back(std::move( BufferType(comm, j, 8192*1024)));
//        }

        BaseIterType begin, end;
        size_t chunkRead = 0;
        int j = 0;
        bool copying = false;

        while ((chunkRead = loader.getNextChunkAtomic(ph, begin, end, chunkSize, copying)) > 0) {
          ++j;

#pragma omp atomic
          ++i;

          IteratorType fastq_start(parser, begin, end);
          IteratorType fastq_end(parser, end);

          SequenceType read;

          for (; fastq_start != fastq_end; ++fastq_start)
          {
            // get data, and move to next for the other threads

            // first get read
            read = *fastq_start;
//            SequenceType::allocCopy(*fastq_start, read);


            // then compute
            op(read, i, counts);
//            SequenceType::deleteCopy(read);


            // release resource
          }
          if (i % 100000 == 0)
            INFO("Level 1: rank " << rank << " thread " << omp_get_thread_num() << " processed chunk " << i);


          if (copying) delete [] begin;
        }
//        //private variables
//        SequenceType read;
//        bool hasData = false;
//        int li = 0;
//        int tid2 = omp_get_thread_num();
//        INFO("Level 1: compute: thread id = " << tid2);
//        int j = 0;
//        do {
//
//          // single thread at a time getting data.
//#pragma omp critical
//          {
//            // get data, and move to next for the other threads
//            atEnd = (fastq_start == fastq_end);
//#pragma omp flush(atEnd)
//            hasData = !atEnd;
//
//            if (hasData)
//            {
//              SequenceType::allocCopy(*fastq_start, read);
//              li = i;
//              ++fastq_start;
//              ++i;
//#pragma omp flush(i, fastq_start)
//            }
//
//          }
//
//          // now do computation.
//          if (hasData) {
//            op(read, li, buffers, counts);
//            ++j;
//
//            SequenceType::deleteCopy(read);
//
//            if (li % 1000000 == 0)
//              INFO("Level 1: rank " << rank << " thread " << tid2 << " processed " << li << " reads");
//          }
//
//
//        } while (!atEnd);


        INFO("Level 1: rank " << rank << " thread " << omp_get_thread_num() << " processed total of " << j << " chunks");

//        // send the last part out.
//        for (int j = 0; j < nprocs; ++j) {
//          buffers[(j+ rank) % (nprocs)].flush();
//        }
      }  // compute threads parallel


      t2 = std::chrono::high_resolution_clock::now();
      time_span = std::chrono::duration_cast<std::chrono::duration<double>>(
          t2 - t1);
      INFO(
          "Level 0: Computation: rank " << rank << " thread " << omp_get_thread_num() << " elapsed time: " << time_span.count() << "s.");


      for (size_t j = 0; j < counts.size(); ++j) {
        INFO("rank " << rank << " COUNTS by thread "<< j << ": "<< counts[j].c);
      }


}


/*
 * second pattern for threading
 */
template <typename Compute>
void compute_OMP_ParFor(FileLoaderType &loader, PartitionHelperType &ph,
                              int &nthreads, IndexType &index,
                              const int &nprocs, const int &rank, const int chunkSize)
{


  INFO("NO MASTER");

  ///  get the start and end iterator position.
  ParserType parser;



#ifdef USE_OPENMP
    INFO("rank " << rank << " MAX THRADS = " << nthreads);
    omp_set_nested(1);
    omp_set_dynamic(0);
#endif

    int senders = nthreads;

//    INFOF("senders = %d\n", senders);

      std::chrono::high_resolution_clock::time_point t1, t2;
      std::chrono::duration<double> time_span;

      INFO("Level 0: compute tid = " << omp_get_thread_num());

      t1 = std::chrono::high_resolution_clock::now();

      // if atEnd, done.
//      bool atEnd = false;
      int i = 0;

      std::vector<bliss::index::countType> counts;
      for (int j = 0; j < nthreads; ++j) {
        bliss::index::countType c;
        c.c = 0;
        counts.push_back(c);
      }


      // VERSION 2.  uses the fastq iterator as the queue itself, instead of master/slave.
      //   at this point, no strong difference.
#pragma omp parallel num_threads(nthreads) shared(loader, parser, ph, i, counts, nprocs, nthreads, rank) default(none)
      {
        Compute op(nprocs, rank, nthreads);

//        std::vector<  BufferType > buffers;
//        for (int j = 0; j < nprocs; ++j) {
//          // over provision by the number of threads as well.  (this is okay for smaller number of procs)
//            buffers.push_back(std::move( BufferType(comm, j, 8192*1024)));
//        }

        BaseIterType begin, end;
        size_t chunkRead = 0;
        bool copying = false;
        auto r = loader.getRange();
        SequenceType read;

#pragma omp for schedule(dynamic, 10) private(chunkRead, begin, end, read)
        for (int c = r.start; c < r.end; c += chunkSize) {
          chunkRead = loader.getNextChunkAtomic(ph, begin, end, chunkSize, copying);

          if (chunkRead > 0) {

            IteratorType fastq_start(parser, begin, end);
            IteratorType fastq_end(parser, end);
            ++i;

            for (; fastq_start != fastq_end; ++fastq_start)
            {
              // get data, and move to next for the other threads

              // first get read
              read = *fastq_start;
  //            SequenceType::allocCopy(*fastq_start, read);


              // then compute
              op(read, i, counts);
  //            SequenceType::deleteCopy(read);


              // release resource
            }

            if (copying) delete [] begin;
          }
        }

        if (i % 100000 == 0)
          INFO("Level 1: rank " << rank << " thread " << omp_get_thread_num() << " processed chunk " << i);

//        //private variables
//        SequenceType read;
//        bool hasData = false;
//        int li = 0;
//        int tid2 = omp_get_thread_num();
//        INFO("Level 1: compute: thread id = " << tid2);
//        int j = 0;
//        do {
//
//          // single thread at a time getting data.
//#pragma omp critical
//          {
//            // get data, and move to next for the other threads
//            atEnd = (fastq_start == fastq_end);
//#pragma omp flush(atEnd)
//            hasData = !atEnd;
//
//            if (hasData)
//            {
//              SequenceType::allocCopy(*fastq_start, read);
//              li = i;
//              ++fastq_start;
//              ++i;
//#pragma omp flush(i, fastq_start)
//            }
//
//          }
//
//          // now do computation.
//          if (hasData) {
//            op(read, li, buffers, counts);
//            ++j;
//
//            SequenceType::deleteCopy(read);
//
//            if (li % 1000000 == 0)
//              INFO("Level 1: rank " << rank << " thread " << tid2 << " processed " << li << " reads");
//          }
//
//
//        } while (!atEnd);


        INFO("Level 1: rank " << rank << " thread " << omp_get_thread_num() << " processed total of " << i << " chunks");

//        // send the last part out.
//        for (int j = 0; j < nprocs; ++j) {
//          buffers[(j+ rank) % (nprocs)].flush();
//        }
      }  // compute threads parallel


      t2 = std::chrono::high_resolution_clock::now();
      time_span = std::chrono::duration_cast<std::chrono::duration<double>>(
          t2 - t1);
      INFO(
          "Level 0: Computation: rank " << rank << " thread " << omp_get_thread_num() << " elapsed time: " << time_span.count() << "s.");

      for (size_t j = 0; j < counts.size(); ++j) {
        INFO("rank " << rank << " COUNTS by thread : "<< counts[j].c);
      }

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
  std::string filename("/home/Tony Pan <tpan7@gatech.edu>/src/bliss/test/data/test.fastq");

  int WithMaster = 1;
  if (argc > 1)
  {
    WithMaster = atoi(argv[1]);
  }

  int nthreads = omp_get_max_threads();
  if (argc > 2)
  {
    nthreads = atoi(argv[2]);
    if (nthreads == -1)
      nthreads = omp_get_max_threads();

  }

  int chunkSize = 4000;
  if (argc > 3)
  {
    chunkSize = atoi(argv[3]);
  }


//  std::string filename("/mnt/data/1000genome/HG00096/sequence_read/SRR077487_1.filt.fastq");
  if (argc > 4)
  {
    filename.assign(argv[4]);
  }


  //////////////// initialize MPI and openMP
  int rank = 0, nprocs = 1;


  /////////////// initialize output variables
  IndexType index;


  /////////////// start processing.  enclosing with braces to make sure loader is destroyed before MPI finalize.
  {
    std::chrono::high_resolution_clock::time_point t1, t2;
    std::chrono::duration<double> time_span;

    //////////////// now partition and open the file.
    t1 = std::chrono::high_resolution_clock::now();

    FileLoaderType loader(filename);  // this handle is alive through the entire execution.
    // adjust the partition of the file
    PartitionHelperType ph;
    loader.adjustRange(ph);
    loader.load();


    t2 = std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    INFO("MMap rank " << rank << " elapsed time: " << time_span.count() << "s.");

    INFO( rank << " file partition: " << loader.getRange() );



    /////////////// now process the file using version with master.
    if (WithMaster == 1) {
      // do some work using openmp  // version without master
      compute_OMP_WithMaster<ComputeType>(loader, ph, nthreads, index, nprocs, rank, chunkSize);
    } else if (WithMaster == 2) {
      /////////////// now process the file using version with master.
      // do some work using openmp  // version without master
      compute_OMP_ParFor<ComputeType>(loader, ph, nthreads, index, nprocs, rank, chunkSize);

    } else {
      /////////////// now process the file using version with master.
      // do some work using openmp  // version without master
      compute_OMP_NoMaster<ComputeType>(loader, ph, nthreads, index, nprocs, rank, chunkSize);
    }

    loader.unload();

  }  // scope to ensure file loader is destroyed.

  return 0;
}

#define USE_MPI (OLD_USE_MPI)
#undef OLD_USE_MPI
