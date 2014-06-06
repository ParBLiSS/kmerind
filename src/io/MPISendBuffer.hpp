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
#include <unistd.h>  // for usleep

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
        static constexpr int END_TAG = 0;
        typedef T ValueType;

        MPISendBuffer(MPI_Comm _comm, const int _targetRank, const size_t nbytes)
           : accepting(true), targetRank(_targetRank), comm(_comm), send_request(MPI_REQUEST_NULL), tid(0), rank(0), total(0)  {
          // initialize internal buffers (double buffering)
          capacity = nbytes / sizeof(T);

          active_buf.reserve(capacity);
          inactive_buf.reserve(capacity);

         // printf("targetrank is %d \n", targetRank);

          MPI_Comm_rank(comm, &rank);
          int nprocs;
          MPI_Comm_size(comm, &nprocs);
          usleep_duration = (1000 + nprocs - 1) / nprocs;

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

        void buffer(T val) {
          assert(accepting);   // if this assert fails, the calling thread's logic is faulty.



          // locking.
    #ifdef USE_OPENMP
          if (THREAD_SAFE)
            omp_set_lock(&lock);
    #endif
          // store value
          //printf("val: %lu %lu %f\n", val.id.composite, val.kmer, val.qual);

          // careful.  when growing, the content is automatically copied to new and destroyed in old.
          //   the destruction could cause double free error or segv.

          active_buf.push_back(std::forward(val));
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
          //printf("flushing %d target %d.\n", rank, targetRank);

          accepting = false;

          // no more entries to buffer. send all remaining.
          send();

          // one more time to clear double buffer and wait for previous send to finish.
          send();

          // send termination signal (send empty message).
          //printf("%d done sending to %d\n", rank, targetRank); fflush(stdout);
          MPI_Send(nullptr, 0, MPI_INT, targetRank, END_TAG, comm);
          printf("%d.%d  %ld flushed to %d.\n", rank, tid, total, targetRank); fflush(stdout);
        }

        size_t size() {
          return active_buf.size();
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
        int tid;
        int rank;

        size_t total;
        int usleep_duration;

    #ifdef USE_OPENMP
          omp_lock_t lock;
    #endif


        void send() {

          // only called from within locked code block

          //printf("SEND %d -> %d,  inactive size: %lu, active size: %lu\n", rank, targetRank, inactive_buf.size(), active_buf.size());

          // check inactive buffer is sending?
          int finished = 0;
          do {
            // repeatedly check to see if request was completed.
            // okay to test against a NULL request.

            // if being sent, wait for that to complete
            MPI_Test(&send_request, &finished, MPI_STATUS_IGNORE );
            usleep(usleep_duration);

            //printf("SEND %d -> %d Waiting for prev send \n", rank, targetRank); fflush(stdout);

          } while(finished == 0);

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
            MPI_Isend(inactive_buf.data(), inactive_buf.size() * sizeof(T), MPI_UNSIGNED_CHAR, targetRank, tid + 1, comm, &send_request);
            //printf("SEND %d -> %d, number of entries %ld, total %ld\n", rank, targetRank, inactive_buf.size(), total); fflush(stdout);

          }

        }
    };

  } /* namespace io */
} /* namespace bliss */

#endif /* MPISENDBUFFER_HPP_ */
