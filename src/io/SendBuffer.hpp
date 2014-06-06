/**
 * @file		SendBuffer.hpp
 * @ingroup
 * @author	tpan
 * @brief
 * @details MPISendBuffer uses double buffering.  we can't use that here as the goal is to have a compute thread populate the buffer but have an MPI thread
 *            consume this and other buffers.  reusing would require the MPI thread to track where the buffer came from, to update when the MPI send request
 *            is complete, and for compute thread to wait for all previously queued buffers.
 *          so instead, create a new buffer each time.
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef SENDBUFFER_HPP_
#define SENDBUFFER_HPP_

#include <cassert>
#include <unistd.h>  // for usleep

#include <thread>
#include <mutex>
#include <condition_variable>
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
     *            approach:
     *                have reference to downstream queue.
     *                when full, create a new vector, swap, enqueue std::pair with target rank and stdvector.
     *
     */
    template<typename T, bool THREAD_SAFE=true>
    class SendBuffer {
        // need some default MPI buffer size, then pack in sizeof(T) blocks as many times as possible.

      protected:
        // 2 buffers per target
        std::vector<T> buf;

        // if buffer is still accepting stuff
        bool accepting;

        // rank of target entry.
        const int targetRank;
        size_t capacity;

        mutable std::mutex mutex;

      public:
        typedef T ValueType;

       SendBuffer(const int _targetRank, const size_t nbytes)
           : accepting(true), targetRank(_targetRank) {
          // initialize internal buffers (double buffering)
          capacity = nbytes / sizeof(T);

          buf.reserve(capacity);
        }

        virtual ~SendBuffer() {}

        bool buffer(T val) {
          if (THREAD_SAFE)
            std::unique_lock<std::mutex> lock(mutex);

          assert(accepting);   // if this assert fails, the calling thread's logic is faulty.
          if (isFull()) return false;

          // careful.  when growing, the content is automatically copied to new and destroyed in old.
          //   the destruction could cause double free error or segv.

          buf.push_back(std::move(val));
          return true;
        }

        size_t size() {
          if (THREAD_SAFE)
            std::unique_lock<std::mutex> lock(mutex);
          return buf.size();
        }

        bool isFull() {
          if (THREAD_SAFE)
            std::unique_lock<std::mutex> lock(mutex);
          return buf.size() == capacity;
        }

      protected:

        void send() {
          if (buf.size() == 0) return;

          // TODO: send the data as bytes.  This assumes the receive end knows how to parse the elements

          // clear the inactive buf
          inactive_buf.clear();
          assert(inactive_buf.size() == 0);
          int s = active_buf.size();
          // swap active and inactive buffer  (this swaps the contents)
          inactive_buf.swap(active_buf);
          assert(active_buf.size() == 0);
          assert(inactive_buf.size() == s);

          buf.data(),

            MPI_Isend(inactive_buf.data(), inactive_buf.size() * sizeof(T), MPI_UNSIGNED_CHAR, targetRank, tid + 1, comm, &send_request);
            //printf("SEND %d -> %d, number of entries %ld, total %ld\n", rank, targetRank, inactive_buf.size(), total); fflush(stdout);


        }
    };

  } /* namespace io */
} /* namespace bliss */

#endif /* SENDBUFFER_HPP_ */
