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
#include "config.hpp"

#if defined(USE_MPI)
#include "mpi.h"
#endif

#include <unistd.h>
#include <string.h>
#include <cstdio>
#include <cmath>

#include <iostream>
//#include <thread>
#include <vector>
#include <queue>
#include <unordered_map>
#include <chrono>

#include "omp.h"

#include <xmmintrin.h>

#include "utils/logging.h"
#include "utils/constexpr_array.hpp"
#include "common/base_types.hpp"
#include "common/alphabets.hpp"
#include "common/AlphabetTraits.hpp"
#include "partition/range.hpp"
#include "iterators/buffered_transform_iterator.hpp"
#include "io/fastq_loader.hpp"
#include "io/MPISendBuffer.hpp"
#include "io/SendBuffer.hpp"
#include "index/kmer_index_element.hpp"
#include "index/kmer_index_functors.hpp"
#include "index/kmer_index_generation.hpp"
#include "concurrent/lockfree_queue.hpp"

/*
 * TYPE DEFINITIONS
 */
// define kmer index type
typedef bliss::index::KmerSize<21>                                KmerSize;
typedef uint64_t                                                  KmerType;
typedef float                                                     QualityType;
typedef DNA                                                       Alphabet;
typedef bliss::index::KmerIndexElementWithIdAndQuality<KmerSize, KmerType, bliss::io::FASTQSequenceId, QualityType>
                                                                  KmerIndexType;

// Range type
typedef bliss::partition::range<size_t>                           RangeType;

// define buffers and queue element types.
constexpr bool thread_safe = false;
typedef bliss::io::MPISendBuffer<KmerIndexType, thread_safe>      MPIBufferType;
typedef bliss::io::SendBuffer<KmerIndexType, false>               BufferType;
typedef typename BufferType::ExportType                           SendQueueElementType;
typedef std::pair<int, KmerIndexType*>                            RecvQueueElementType;

// define the input file type
// raw data type :  use CharType
typedef bliss::io::FASTQLoader<CharType, true, false>             FileLoaderType;
typedef typename FileLoaderType::BlockIteratorType                BaseIterType;
typedef typename std::iterator_traits<BaseIterType>::value_type   BaseValueType;


// define read type
typedef bliss::io::SequenceWithQuality<BaseIterType, QualityType>
                                                                  SequenceType;


// define kmer and quality operator types
typedef bliss::index::generate_kmer<SequenceType, KmerIndexType>  KmerOpType;

typedef bliss::index::SangerToLogProbCorrect<QualityType>         EncodeType;

typedef bliss::index::generate_qual<SequenceType, KmerSize, QualityType, EncodeType >
                                                                  QualOpType;

// define the index storage type
typedef std::unordered_multimap<KmerType, KmerIndexType>          IndexType;

// define the transform iterator type
typedef bliss::io::FASTQParser<BaseIterType, QualityType>
                                                                  ParserType;
typedef bliss::io::SequencesIterator<ParserType, BaseIterType>       IteratorType;


// define the computation type.
typedef bliss::index::KmerIndexGeneratorWithQuality<KmerOpType, MPIBufferType, bliss::index::XorModulus<KmerType>, QualOpType>
                                                                  MPIComputeType;
typedef bliss::index::KmerIndexGeneratorWithQuality<KmerOpType, BufferType, bliss::index::XorModulus<KmerType>, QualOpType>
                                                                  ComputeType;




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
    while (hasMessage == 0) {
      MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, comm, &hasMessage, &status);
      //usleep(usleep_duration);
      _mm_pause();
      //printf("%dw",rank); fflush(stdout);
    }
    //printf("\n");
    hasMessage = 0;

    src = status.MPI_SOURCE;
    tag = status.MPI_TAG;
    MPI_Get_count(&status, MPI_UNSIGNED_CHAR, &received);   // length is always > 0?

    if (tag == MPIBufferType::END_TAG) {
      // end of messaging.
      MPI_Recv(reinterpret_cast<unsigned char*>(array), received, MPI_UNSIGNED_CHAR, src, tag, comm, MPI_STATUS_IGNORE);
      --n_senders;
      //printf("RECV rank %d receiving END signal %d from %d, num sanders remaining is %d\n", rank, received, src, n_senders); fflush(stdout);

      continue;
    }

    //printf("RECV %d receiving %d bytes from %d.%d\n", rank, received, src, tag - 1); fflush(stdout);
    MPI_Recv(reinterpret_cast<unsigned char*>(array), received, MPI_UNSIGNED_CHAR, src, tag, comm, MPI_STATUS_IGNORE);

    // if message size is 0, then no work to do.
    if (received == 0)
      continue;

    //////// logic for storing index.

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







// TODO: make these function with/without MPI, with/without OpenMP.
// TODO: 1. make a communicator that encapsulates the nprocs, rank, nthread, tid, etc - this refactors KmerIndexGenerator so there is no dependency on nprocs, rank, and tid.
// TODO: 2. refactor these "compute" methods so that the loop iterations are handled by a ThreadModel.
// TODO: 3. make the communicator short circuit when copying to local, and lock when multithreading.
// TODO: 4. communicator should have send buffer and receive buffer (actual data structure)

/*
 * first pattern for threading - has a master thread  - does not work well.  need each thread to get its own chunk range based on tid, so that means
 * each task needs to call the getNextChunkRange method.  but then the master thread cannot check the length of the chunk received.
 * the alternative is for master thread to use a counter. but if I have a counter, then might as well use parfor.
 *
 * function removed.  for this pattern, see omp_patterns.hpp.
 */

/*
 * peer to peer
 */
template <typename Compute>
void compute_MPI_OMP_P2P(FileLoaderType &loader,
                              const int &nthreads, IndexType &index,
                              MPI_Comm &comm, const int &nprocs, const int &rank, const int &chunkSize)
{

  ///  get the start and end iterator position.

#ifdef USE_OPENMP
    INFO("rank " << rank << " MAX THRADS = " << nthreads);
    omp_set_nested(1);
    omp_set_dynamic(0);
#endif


    int senders = nthreads;
#ifdef USE_MPI
    senders = 0;
    int nt = nthreads;
    MPI_Allreduce(&nt, &senders, 1, MPI_INT, MPI_SUM, comm);
#endif

    //printf("senders = %d\n", senders);

#pragma omp parallel sections num_threads(2) shared(comm, nprocs, rank, index, nthreads, loader, senders) OMP_SHARE_DEFAULT
    {


#pragma omp section
    {
      std::chrono::high_resolution_clock::time_point t1, t2;
      std::chrono::duration<double> time_span;




      INFO("Level 0: index storage tid = " << omp_get_thread_num());
      t1 = std::chrono::high_resolution_clock::now();
      networkread(comm, nprocs, rank, 8192 * 1024, senders, index);
      t2 = std::chrono::high_resolution_clock::now();
      time_span =
          std::chrono::duration_cast<std::chrono::duration<double>>(
              t2 - t1);
      INFO(
          "Level 0: Index storage: rank " << rank << " elapsed time: " << time_span.count() << "s, " << index.size() << " records.");

      printf("rank %d done network receive\n", rank);
    } // omp network read section

#pragma omp section
    {  // compute threads
      std::chrono::high_resolution_clock::time_point t1, t2;
      std::chrono::duration<double> time_span;

      INFO("Level 0: compute tid = " << omp_get_thread_num());

      t1 = std::chrono::high_resolution_clock::now();

      // if atEnd, done.
//      bool atEnd = false;

      // create local buffers
      std::vector<  std::vector<MPIBufferType > > buffers;
      // buffers.reserve(nprocs); don't reserve.  no need.
      for (int j = 0; j < nthreads; ++j) {

        std::vector< MPIBufferType > temp;
        for (int i = 0; i < nprocs; ++i) {
          temp.push_back(std::move(MPIBufferType(comm, i, 8192*1024)));
        }

        buffers.push_back(std::move(temp));
      }
      INFO("Level 1: compute num buffers = " << buffers.size());



      // VERSION 2.  uses the fastq iterator as the queue itself, instead of master/slave.
      //   at this point, no strong difference.
#pragma omp parallel num_threads(nthreads) shared(loader, buffers, nprocs, nthreads, rank) OMP_SHARE_DEFAULT
      {
        ParserType parser;

        int tid = 0;
#ifdef USE_OPENMP
        tid = omp_get_thread_num();
#endif
        Compute op(nprocs, rank, nthreads);


//        std::vector<  MPIBufferType > buffers;
//        for (int j = 0; j < nprocs; ++j) {
//          // over provision by the number of threads as well.  (this is okay for smaller number of procs)
//            buffers.push_back(std::move( MPIBufferType(comm, j, 8192*1024)));
//        }

        int i = 0;
        int j = 0;
        RangeType r = loader.getNextChunkRange(tid);

        while (r.size() > 0) {

          auto chunk = loader.getChunk(tid, r);

          IteratorType fastq_start(parser, chunk.begin(), chunk.end(), r);
          IteratorType fastq_end(parser, chunk.end(), r);

          SequenceType read;

          for (; fastq_start != fastq_end; ++fastq_start)
          {
            // get data, and move to next for the other threads

            // first get read
            read = *fastq_start;

            // then compute
            op(read, buffers[tid]);
            ++i;

            if (i % 20000 == 0)
              INFO("Level 1: rank " << rank << " thread " << tid << " processed seq " << i);
          }

          r = loader.getNextChunkRange(tid);
          ++j;

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

        INFO("Level 1: rank " << rank << " thread " << tid << " processed total of " << j << " chunks");

//        // send the last part out.
//        for (int j = 0; j < nprocs; ++j) {
//          buffers[(j+ rank) % (nprocs)].flush();
//        }
      }  // compute threads parallel


      t2 = std::chrono::high_resolution_clock::now();
      time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
      INFO("Level 0: Computation: rank " << rank << " thread " << omp_get_thread_num() << " elapsed time: " << time_span.count() << "s.");



      // send the last part out.
//#pragma omp parallel for OMP_SHARE_DEFAULT  shared(nprocs, nthreads, rank, buffers)
      for (int j = 0; j < nthreads; ++j) {
        for (int i = 0; i < nprocs; ++i) {
        //printf("rank %d thread %d flushing buffer for buffer %d\n", rank, omp_get_thread_num(), (i+rank) % (nprocs * nthreads));
          buffers[j][(i + rank) % nprocs].flush();   /// + rank to avoid all threads sending to the same node.
        }
      }

      INFO("flush rank " << rank << " thread " << omp_get_thread_num() << " elapsed time: " << time_span.count() << "s.");
      //printf("flush rank %d thread %d elapsed time %f\n", rank, omp_get_thread_num(), time_span.count());


      //printf("rank %d done compute\n", rank);
    } // omp compute section


    }  // outer parallel sections
    //printf("rank %d done\n", rank);
}


/*
 * parfor is awkward as well.  the check for chunk done is empty range.  using the parallel for creates another layer of atomic updates.
 */


/**
 * receiver for MPI communications.
 */
struct RecvProcessor {

  void operator()(RecvQueueElementType &buf, IndexType &kmers) {
    if (buf.first == 0) return;

    for (int i = 0; i < buf.first; ++i) {
      kmers.insert(std::move(IndexType::value_type(buf.second[i].kmer, buf.second[i])));
    }
  }
};


/*
 * single thread for mpi messages.
 */
template <typename Compute>
void computeP2P(FileLoaderType &loader,
                              const int &nthreads, IndexType &index,
                              MPI_Comm &comm, const int &nprocs, const int &rank, const int &chunkSize)
{

  ///  get the start and end iterator position.

#ifdef USE_OPENMP
    INFO("rank " << rank << " MAX THRADS = " << nthreads);
    omp_set_nested(1);
    omp_set_dynamic(0);
#endif


    // shared flag to indicate receiving is done.
    int mpi_senders = nprocs;
    // shared flag to indicate compute is done.
    int computes = nthreads;  // for now, oversubscribe by 2.

    // outbound message structure
    // create shared queue for MPI outbound messages. produced by compute threads via MPISendBuffer. consumed by MPI comm thread.
    bliss::concurrent::ThreadSafeQueue<SendQueueElementType> sendQueue(nthreads * 2);

    // define inbound message structure
    // create shared queue for MPI inbound messages. produced by MPI comm thread.  consumed by index hashing thread.
    bliss::concurrent::ThreadSafeQueue<RecvQueueElementType> recvQueue(nprocs * 2);



#pragma omp parallel sections num_threads(3) shared(comm, nprocs, rank, index, nthreads, loader, mpi_senders, computes, sendQueue, recvQueue, ompi_mpi_unsigned_char, ompi_mpi_int, ompi_request_null) OMP_SHARE_DEFAULT
    {

#pragma omp section
      {  // section for handling the received MPI messages

        RecvQueueElementType re;
        RecvProcessor rp;

        /// process the inbound messages until MPI comm thread indicates no more messages.
        while (mpi_senders > 0) {
          /// check for messages
          auto result = std::move(recvQueue.tryPop());
          if (result.first) {
            re = result.second;
            //printf("%d received %d. rec queue %lu\n", rank, re.first, recvQueue.size());

            RecvQueueElementType lre(std::move(result.second));

            /// process messages;
            rp(lre, index);

            // clean up
            delete [] lre.second;

          } else {
            //usleep(50);
            _mm_pause();
          }



#pragma omp flush(mpi_senders)
        }

        /// received done message.  now flush the remainder
        auto result = std::move(recvQueue.tryPop());
        while (result.first) {

//          printf("flushing recv queue. %lu\n", recvQueue.size());

          RecvQueueElementType lre(std::move(result.second));

          /// process messages;
          rp(lre, index);

          // clean up.
          delete [] lre.second;

          /// check for next.
          result = std::move(recvQueue.tryPop());
        }

        printf("finally: %lu -> %lu\n", sendQueue.getSize(), recvQueue.getSize());
      }  // section for handling the received MPI messages



#pragma omp section
    {  // section for handling the MPI communications
      std::chrono::high_resolution_clock::time_point t1, t2;
      std::chrono::duration<double> time_span;

      /// loop until : 1. no more srcs.  2. no more compute threads.
      bool computeDone = false;  // 2 is computing, 1 is flushing, 0 is done
      bool recvDone = false;
      SendQueueElementType se;


      MPI_Status status;

      // for receiving.
      int hasMessage = 0;
      int src;
      int tag;
      int received = 0;

      // for sending
      std::queue<std::pair<MPI_Request, RecvQueueElementType> > recvInProgress;
      std::queue<std::pair<MPI_Request, SendQueueElementType> > sendInProgress;

      while (true) {

        /// clear out any completed requests.
        // MPI_Test - not using immediate mode Send/Recv.
        size_t waiting = recvInProgress.size();  // getting the size ahead of time.
        int fin = 0;
        for (int i = 0; i < waiting; ++i) {
          std::pair<MPI_Request, RecvQueueElementType> el = std::move(recvInProgress.front());
          recvInProgress.pop();

          MPI_Test(&(el.first), &fin, MPI_STATUS_IGNORE);

          if (fin == 1) {
            // finished receiving.  insert into recvQueue
            while (!recvQueue.tryPush(std::move(el.second))) {
              //usleep(50);
              _mm_pause();

            }
          } else {
            recvInProgress.push(std::move(el));
          }
        }
        // now check send
        waiting = sendInProgress.size();
        for (int i = 0; i < waiting; ++i) {
          std::pair<MPI_Request, SendQueueElementType> el = std::move(sendInProgress.front());
          sendInProgress.pop();

          MPI_Test(&(el.first), &fin, MPI_STATUS_IGNORE);

          if (fin == 0) {
            // finished sending.  insert into recvQueue
            sendInProgress.push(std::move(el));
          }
        }


        // try to receive 1 message.
        if (mpi_senders > 0) {
        /// probe for messages
          hasMessage = 0;
          MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, comm, &hasMessage, &status);

          if (hasMessage > 0) {
            //// if message, check if end message.  if so, decrement mpi_sender.
            //// if data message, put into queue

            src = status.MPI_SOURCE;
            tag = status.MPI_TAG;
            MPI_Get_count(&status, MPI_UNSIGNED_CHAR, &received);   // length is always > 0?

            if (tag == MPIBufferType::END_TAG) {
              // end of messaging.
              MPI_Recv(nullptr, received, MPI_UNSIGNED_CHAR, src, tag, comm, MPI_STATUS_IGNORE);
#pragma omp atomic
              --mpi_senders;
              printf("RECV rank %d receiving END signal %d from %d, num sanders remaining is %d\n", rank, received, src, mpi_senders);
#pragma omp flush(mpi_senders)
            } else {
              KmerIndexType* array = new KmerIndexType[received / sizeof(KmerIndexType)];
              //printf("created temp storage for read, capacity = %ld\n", capacity); fflush(stdout);
              MPI_Request req;

              printf("RECV %d->%d receiving %d bytes with tag %d\n", src, rank, received, tag);
              MPI_Irecv(reinterpret_cast<unsigned char*>(array), received, MPI_UNSIGNED_CHAR, src, tag, comm, &req);

              // put into queue. - std::move version.  this part has to complete.
              recvInProgress.push(std::move(std::pair<MPI_Request, RecvQueueElementType>(
                  req,
                  std::move(RecvQueueElementType(received/ sizeof(KmerIndexType), array))
              )));
              //printf("%d received %d. in progress: %lu\n", rank, received, recvInProgress.size());

            }
          }
        } else {
          if ((recvDone == false) && recvInProgress.empty())
            recvDone = true;
        }


#pragma omp flush(computes)
       if ((computes == 0) && sendQueue.isEmpty()) {

         if (sendInProgress.empty() && (computeDone == false)) {
           /// no more generators.  and queue is empty.
           computeDone = true;
           printf("%d compute done\n", rank);

           // now send end message to everyone except self.


           // next send to all others.
           printf("%d terminating\n", rank);
           MPI_Request endRequests[nprocs];
           for (int i = 0; i < nprocs; ++i) {
             if (i == rank) {
               // first deal with sending to self.
      #pragma omp atomic
               // local (same rank), just update.
               --mpi_senders;
      #pragma omp flush(mpi_senders)
               endRequests[i] = MPI_REQUEST_NULL;

             } else {
               // send end message.
               MPI_Isend(nullptr, 0, MPI_INT, i, MPIBufferType::END_TAG, comm, &(endRequests[i]));
             }
           }

           // wait for them all to finish
           printf("%d wait for everyone to finish\n", rank);
           MPI_Waitall(nprocs, endRequests, MPI_STATUSES_IGNORE);
         }
       } else {
         /// try dequeue and iSend.
         auto result = std::move(sendQueue.tryPop());
         if (result.first) {
           se = result.second;

           /// iSend.  put request into a vector..

           if (se.first == rank) {
             // local, directly handle by creating an output object and directly insert into the recv queue
             KmerIndexType* array = new KmerIndexType[se.second.size()];
             memcpy(array, se.second.data(), se.second.size() * sizeof(KmerIndexType));
             RecvQueueElementType temp(se.second.size(), array);
             while(!recvQueue.tryPush(std::move(temp))) {
               //usleep(50);
               _mm_pause();

             }

             //printf("%d->%d Sent %lu locally. in queue: sendQueue size %lu\n", rank, rank, se.second.size(), sendQueue.size());

           } else {
             MPI_Request req;

             SendQueueElementType lse(std::move(se));
             MPI_Isend(lse.second.data(), lse.second.size() * sizeof(KmerIndexType), MPI_UNSIGNED_CHAR, lse.first, rank+1, comm, &req);
             // lse is automatically cleaned up.

             //printf("%d->%d Sent %lu remotely. in queue: sendQueue size %lu\n", rank, lse.first, lse.second.size(), sendQueue.size());

             sendInProgress.push(std::move(std::pair<MPI_Request, SendQueueElementType>(req, std::move(lse))));
           }

         }

       }

        if (computeDone && recvDone) {
          //printf("inprogress:  send queue %lu -> %lu -> %lu -> %lu recv queue\n", sendQueue.size(), sendInProgress.size(), recvInProgress.size(), recvQueue.size());
          break;
        }
//        usleep(50);
        _mm_pause();

      }
    }  // section for handling the MPI communications

#pragma omp section
    {    // section for processing the data and generating the kmers.
      std::chrono::high_resolution_clock::time_point t1, t2;
      std::chrono::duration<double> time_span;

      INFO("Level 0: compute tid = " << omp_get_thread_num());

      t1 = std::chrono::high_resolution_clock::now();



      // create local buffers
      std::vector<  std::vector<BufferType > > buffers;
      // buffers.reserve(nprocs); don't reserve.  no need.
      for (int j = 0; j < nthreads; ++j) {

        std::vector< BufferType > temp;
        for (int i = 0; i < nprocs; ++i) {
          temp.push_back(std::move(BufferType(i, 8192*1024)));
        }

        buffers.push_back(std::move(temp));
      }
      INFO("Level 1: compute num buffers = " << buffers.size());

      int tcount = 0;


      // VERSION 2.  uses the fastq iterator as the queue itself, instead of master/slave.
      //   at this point, no strong difference.
#pragma omp parallel num_threads(nthreads) shared(loader, nprocs, nthreads, rank, sendQueue, computes, buffers) OMP_SHARE_DEFAULT reduction(+:tcount)
      {
        /// initialize variables.
        ParserType parser;

        int tid = 0;
#ifdef USE_OPENMP
        tid = omp_get_thread_num();
#endif
        Compute op(nprocs, rank, nthreads);


        printf("thread %d in %d.  requested %d\n", tid, omp_get_num_threads(), nthreads);


        tcount = 0;
        int j = 0;
        SequenceType read;
        typename FileLoaderType::DataBlockType chunk;
        /// repeatedly get from source, and process.

        RangeType r = loader.getNextChunkRange(tid);
        int targetRank;

        while (r.size() > 0) {
          /// get the actual chunk
           chunk = loader.getChunk(tid, r);

          /// processing
          IteratorType fastq_start(parser, chunk.begin(), chunk.end(), r);
          IteratorType fastq_end(parser, chunk.end(), r);

          for (; fastq_start != fastq_end; ++fastq_start)
          {
            // get data, and move to next for the other threads

            // first get read
            read = *fastq_start;

            // compute and buffer results
            op(read, buffers[tid]);
            ++tcount;

            if (tcount % 20000 == 0)
              INFO("Level 1: rank " << rank << " thread " << tid << " processed seq " << tcount);
          }

          /// now check the buffers and send if needed.
          for (int i = 0; i < buffers[tid].size(); ++i) {
            targetRank =(i + rank) % nprocs;
            if (buffers[tid][targetRank].isFull()) {
              //printf("%d pushing %d->%d into queue (size %lu). entries %lu\n", tid, rank, targetRank, buffers[tid].size(), buffers[tid][targetRank].size());
              sendQueue.waitAndPush(std::move(buffers[tid][targetRank].exportData()));
            }
          }

          r = loader.getNextChunkRange(tid);
          ++j;

        }

        INFO("Level 1: rank " << rank << " thread " << tid << " processed total of " << j << " chunks");

      }  // compute threads parallel


      /// send the last part out.
      int targetRank;
      for (int j = 0; j < buffers.size(); ++j) {
        for (int i = 0; i < buffers[j].size(); ++i) {
          //printf("rank %d thread %d flushing buffer for buffer %d\n", rank, omp_get_thread_num(), (i+rank) % (nprocs * nthreads));
          targetRank =(i + rank) % nprocs;
          //printf("%d flushing %d, %lu\n", omp_get_thread_num(), targetRank, buffers[j][targetRank].size());
          if (buffers[j][targetRank].size() > 0)
            sendQueue.waitAndPush(std::move(buffers[j][targetRank].exportData()));   /// + rank to avoid all threads sending to the same node.
        }
      }
      INFO("flush rank " << rank << " thread " << omp_get_thread_num() << "s.");
      //printf("flush rank %d thread %d elapsed time %f\n", rank, omp_get_thread_num(), time_span.count());


      /// signaling done for this thread.
      computes = 0;
#pragma omp flush(computes)



      t2 = std::chrono::high_resolution_clock::now();
      time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
      INFO("Level 0: Computation: rank " << rank << " thread " << omp_get_thread_num() << " elapsed time: " << time_span.count() << "s, " << tcount << "reads");

      INFO("rank " << rank << "COUNTS total: " << tcount);

      //printf("rank %d done compute\n", rank);
    }     // section for processing the data and generating the kmers.


  }  // outer parallel sections
    //printf("rank %d done\n", rank);
}






///**
// * source:              source data,  e.g. fastq_loader
// * Compute:             transformation to index
// * ThreadModel:         threading for compute
// * Communicator:        move output of compute to destination.  e.g. MPI send/recv.
// * Destination:         target container, e.g. index hash table
// *
// * compute is single threaded.
// * threadModel determines whether we are using master/slave or demand driven bag of tasks
// * communicator contains buffer, and is responsible for moving data around.
// * Destination is the target data structure.
// *
// * flow is source -> compute under thread model -> buffer via communicator --> send --> recv -> Destination.
// */
//template <typename Source, typename Compute, typename ThreadModel, typename Communicator, typename Destination>
//void compute(Source src, Destination dest) {
//
//
//
//}

template<typename ComputeType>
struct RunTask {
    void operator()(const std::string &filename, IndexType &index, MPI_Comm &comm, const int &nprocs, const int &rank, const int &nthreads, const int &chunkSize) {
      /////////////// initialize output variables


        std::chrono::high_resolution_clock::time_point t1, t2;
        std::chrono::duration<double> time_span;


        /// open file
        // create FASTQ Loader
        // call get NextPartition Range to block partition,
        // call "load" with the partition range.

        /// processing
        //  get Next Chunk Range
        // then call getChunk with the range.
        // call compute with the DataBlock.

        // set up OMP call:  input file_loader, compute op, and output buffer
        // set up MPI receiver. - capture output
        // (set up MPI sender)

        // query.







        //////////////// now partition and open the file.
        t1 = std::chrono::high_resolution_clock::now();

        // get the file ready for read
        FileLoaderType loader(filename, comm, nthreads, chunkSize);  // this handle is alive through the entire execution.
        RangeType pr = loader.getNextPartitionRange(rank);
        loader.load(pr);


        t2 = std::chrono::high_resolution_clock::now();
        time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
        std::cout << "MMap rank " << rank << " elapsed time: " << time_span.count() << "s." << std::endl;

        std::cout << rank << " file partition: " << loader.getMMapRange() << std::endl;

        t1 = std::chrono::high_resolution_clock::now();
        /////////////// now process the file using version with master.
        // do some work using openmp  // version without master
        compute_MPI_OMP_P2P<ComputeType>(loader, nthreads, index, comm, nprocs, rank, chunkSize);

        t2 = std::chrono::high_resolution_clock::now();
        time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
        std::cout << "Compute rank " << rank << " elapsed time: " << time_span.count() << "s." << std::endl;

        t1 = std::chrono::high_resolution_clock::now();
        loader.unload();
        t2 = std::chrono::high_resolution_clock::now();
        time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
        std::cout << "Unload rank " << rank << " elapsed time: " << time_span.count() << "s." << std::endl;

      }  // scope to ensure file loader is destroyed.

};


template<typename ComputeType>
struct RunTask2 {
    void operator()(const std::string &filename, IndexType &index, MPI_Comm &comm, const int &nprocs, const int &rank, const int &nthreads, const int &chunkSize) {
      /////////////// initialize output variables


        std::chrono::high_resolution_clock::time_point t1, t2;
        std::chrono::duration<double> time_span;


        /// open file
        // create FASTQ Loader
        // call get NextPartition Range to block partition,
        // call "load" with the partition range.

        /// processing
        //  get Next Chunk Range
        // then call getChunk with the range.
        // call compute with the DataBlock.

        // set up OMP call:  input file_loader, compute op, and output buffer
        // set up MPI receiver. - capture output
        // (set up MPI sender)

        // query.







        //////////////// now partition and open the file.
        t1 = std::chrono::high_resolution_clock::now();

        // get the file ready for read
        FileLoaderType loader(filename, comm, nthreads, chunkSize);  // this handle is alive through the entire execution.
        RangeType pr = loader.getNextPartitionRange(rank);
        loader.load(pr);


        t2 = std::chrono::high_resolution_clock::now();
        time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
        std::cout << "MMap rank " << rank << " elapsed time: " << time_span.count() << "s." << std::endl;

        std::cout << rank << " file partition: " << loader.getMMapRange() << std::endl;

        t1 = std::chrono::high_resolution_clock::now();
        /////////////// now process the file using version with master.
        // do some work using openmp  // version without master
        computeP2P<ComputeType>(loader, nthreads, index, comm, nprocs, rank, chunkSize);

        t2 = std::chrono::high_resolution_clock::now();
        time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
        std::cout << "Compute rank " << rank << " elapsed time: " << time_span.count() << "s." << std::endl;

        t1 = std::chrono::high_resolution_clock::now();
        loader.unload();
        t2 = std::chrono::high_resolution_clock::now();
        time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
        std::cout << "Unload rank " << rank << " elapsed time: " << time_span.count() << "s." << std::endl;

      }  // scope to ensure file loader is destroyed.

};



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

  int nthreads = 1;
#ifdef USE_OPENMP
  if (argc > 1)
  {
    nthreads = atoi(argv[1]);
    if (nthreads == -1)
      nthreads = omp_get_max_threads();
  }
#else
  printf("NOT compiled with OPENMP")
#endif

  int chunkSize = 4096;
  if (argc > 2)
  {
    chunkSize = atoi(argv[2]);
  }

  std::string filename("/home/tpan/src/bliss/test/data/test.fastq");
//  std::string filename("/mnt/data/1000genome/HG00096/sequence_read/SRR077487_1.filt.fastq");
  if (argc > 3)
  {
    filename.assign(argv[3]);
  }

  int groupSize = 1;
  int id = 0;
  //////////////// initialize MPI and openMP
#ifdef USE_MPI


          if (nthreads > 1) {

            int provided;
            MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

            if (provided < MPI_THREAD_MULTIPLE) {
              printf("ERROR: The MPI Library Does not have full thread support.\n");
              MPI_Abort(MPI_COMM_WORLD, 1);
            }
          } else {
            MPI_Init(&argc, &argv);
          }

          MPI_Comm comm = MPI_COMM_WORLD;

          MPI_Comm_size(comm, &groupSize);
          MPI_Comm_rank(comm, &id);


          if (id == 0)
            std::cout << "USE_MPI is set" << std::endl;
#else
          static_assert(false, "MPIRunner Used when compilation is not set to use MPI");
#endif
  // replace with MPIRunner



MPI_Barrier(comm);
          IndexType index;





          /////////////// start processing.  enclosing with braces to make sure loader is destroyed before MPI finalize.          index.clear();

          RunTask2<ComputeType> t2;
          t2(filename, index, comm, groupSize, id, nthreads, chunkSize);

          printf("number of entries in index for rank %d is %lu\n", id, index.size());


          MPI_Barrier(comm);



        /////////////// start processing.  enclosing with braces to make sure loader is destroyed before MPI finalize.
        RunTask<MPIComputeType> t;
        t(filename, index, comm, groupSize, id, nthreads, chunkSize);

        printf("MPI number of entries in index for rank %d is %lu\n", id, index.size());

        MPI_Barrier(comm);


        //// query:  use the same file as input.  walk through and generate kmers as before.  send query





  //////////////  clean up MPI.
  MPI_Finalize();

  return 0;
}
