/**
 * @file		MPISendBuffer.hpp
 * @ingroup
 * @author	tpan
 * @brief
 * @details
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef MPISENDBUFFER_HPP_
#define MPISENDBUFFER_HPP_

#include "mpi.h"

#include <cassert>

#include <vector>

#include "config.hpp"

namespace bliss
{
  namespace io
  {


    /**
     * @class			bliss::io::MPISendBuffer
     * @brief     a data structure to buffer and send MPI messages
     * @details   THREAD_SAFE enables omp_lock for OMP thread safety.  this likely introduces a good amount of contention.
     *
     *
     * TODO: 1. abstract the MPI requirement.  in fact, abstract the send/flush calls - make it a generic (thread-safe) buffer.
     */
    template<typename T, bool THREAD_SAFE=true>
    class MPISendBuffer {
        // need some default MPI buffer size, then pack in sizeof(T) blocks as many times as possible.

      public:
        static const int END_TAG = 0;


        MPISendBuffer(MPI_Comm _comm, const int _targetRank, const size_t nbytes)
           : accepting(true), targetRank(_targetRank), comm(_comm), send_request(MPI_REQUEST_NULL), total(0)  {
          // initialize internal buffers (double buffering)
          capacity = nbytes / sizeof(T);
          //printf("created buffer for send, capacity = %ld\n", capacity);

          active_buf.reserve(capacity);
          inactive_buf.reserve(capacity);

         // printf("targetrank is %d \n", targetRank);

          MPI_Comm_rank(comm, &rank);

    #ifdef USE_OPENMP
          if (THREAD_SAFE)
            omp_init_lock(&lock);
          tid = omp_get_thread_num();
    #endif


        }

        virtual ~MPISendBuffer() {

          // buffer should be empty at this point.
          assert(active_buf.size() == 0);
          assert(inactive_buf.size() == 0);

          // all previous messages should be done.

    #ifdef USE_OPENMP
          if (THREAD_SAFE)
            omp_destroy_lock(&lock);
    #endif

        }

        void buffer(const T & val) {
          assert(accepting);   // if this assert fails, the calling thread's logic is faulty.

          // locking.
    #ifdef USE_OPENMP
          if (THREAD_SAFE)
            omp_set_lock(&lock);
    #endif
          // store value
          active_buf.push_back(val);
          // if full, call send (block if other buffer is sending)
          if (active_buf.size() == capacity) {
            send();
          }

    #ifdef USE_OPENMP
          if (THREAD_SAFE)
            omp_unset_lock(&lock);
    #endif
        }

        void flush() {

          accepting = false;

          // no more entries to buffer. send all remaining.
          send();

          // one more time to wait for prev send to finish.
          send();

          // send termination signal (send empty message).
          //printf("%d done sending to %d\n", rank, targetRank); fflush(stdout);
          MPI_Send(nullptr, 0, MPI_INT, targetRank, END_TAG, comm);
          //printf("%d.%d  %ld flushed to %d.\n", rank, tid, total, targetRank); fflush(stdout);
        }

      protected:


        // 2 buffers per target
        std::vector<T> active_buf;
        std::vector<T> inactive_buf;

        // if buffer is still accepting stuff
        bool accepting;

        // rank of target entry.
        const int targetRank;
        size_t capacity;

        MPI_Comm comm;
        // inactive buffer send status
        MPI_Request send_request;
        int nthreads;
        int tid;
        int nprocs;
        int rank;

        size_t total;

    #ifdef USE_OPENMP
          omp_lock_t lock;
    #endif


        void send() {

          // only called from within locked code block
          std::chrono::high_resolution_clock::time_point time1, time2;
          std::chrono::duration<double> time_span;

          // check inactive buffer is sending?
          if (send_request != MPI_REQUEST_NULL) {
            // if being sent, wait for that to complete
            //printf("Waiting for prev send %d -> %d to finish\n", rank, targetRank); fflush(stdout);
            //time1 = std::chrono::high_resolution_clock::now();Traffic

            MPI_Wait(&send_request, MPI_STATUS_IGNORE);
            //time2 = std::chrono::high_resolution_clock::now();
            //time_span = std::chrono::duration_cast<std::chrono::duration<double>>(
            //    time2 - time1);

            //printf("Waiting for prev send %d -> %d finished in %f\n", rank, targetRank, time_span.count()); fflush(stdout);
          }
          // clear the inactive buf
          inactive_buf.clear();
          assert(inactive_buf.size() == 0);
          int s = active_buf.size();
          // swap active and inactive buffer  (this swaps the contents)
          inactive_buf.swap(active_buf);
          assert(active_buf.size() == 0);
          assert(inactive_buf.size() == s);

          // async send inactive buffer if it's not empty.  Requires that remote side uses probe to get the size of data first.

          // TODO: send the data as bytes.  This assumes the receive end knows how to parse the elements
          if (inactive_buf.size() > 0) {
            total += inactive_buf.size();
            //printf("%d sending to %d\n", rank, targetRank);
            MPI_Isend(inactive_buf.data(), inactive_buf.size() * sizeof(T), MPI_UNSIGNED_CHAR, targetRank, tid + 1, comm, &send_request);
            //printf("SEND %d.%d send to %d, number of entries %ld, total %ld\n", rank, pid, targetRank, inactive_buf.size(), total); fflush(stdout);

          }

        }
    };

  } /* namespace io */
} /* namespace bliss */

#endif /* MPISENDBUFFER_HPP_ */
