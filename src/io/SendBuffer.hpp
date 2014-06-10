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
#include <utility>

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
        // rank of target entry.
        int targetRank;

        // 2 buffers per target
        std::vector<T> buf;
        size_t capacity;

        mutable std::mutex mutex;

      public:
        typedef T ValueType;
        typedef std::pair<int, std::vector<T> > ExportType;

       SendBuffer(const int _targetRank, const size_t nbytes)
           : targetRank(_targetRank) {
          // initialize internal buffers (double buffering)
          capacity = 3 * nbytes / (2 * sizeof(T));

          buf.reserve(capacity);

          printf("construct size, capacity: %lu, %lu. datastructure size %lu\n", buf.size(), buf.capacity(), sizeof(*this));
          if (capacity == 0) printf("capacity is 0!!!\n");
        }

       // copy constructor.
       SendBuffer(const SendBuffer<T, THREAD_SAFE>& other) :
         targetRank(other.targetRank), buf(other.buf), capacity(other.capacity) {
         //printf("copy size: %lu\n", buf.size());
         if (capacity == 0) printf("capacity is 0!!!\n");

       }

       // copy assign operator
       SendBuffer<T, THREAD_SAFE>& operator=(const SendBuffer<T, THREAD_SAFE>& other) {
         targetRank = other.targetRank;
         buf = other.buf;  // copy
         capacity = other.capacity;
         //printf("copy assign size: %lu\n", buf.size());
         if (capacity == 0) printf("capacity is 0!!!\n");
         return *this;
       }

       // move constructor.
       SendBuffer(SendBuffer<T, THREAD_SAFE>&& other) :
         targetRank(other.targetRank), buf(std::move(other.buf)), capacity(other.capacity) {
         //printf("move size: %lu\n", buf.size());
         if (capacity == 0) printf("capacity is 0!!!\n");

       }

       // copy assign operator
       SendBuffer<T, THREAD_SAFE>& operator=(SendBuffer<T, THREAD_SAFE>&& other) {
         targetRank = other.targetRank;
         buf = std::move(other.buf);  // copy
         capacity = other.capacity;
         //printf("move assign size: %lu\n", buf.size());
         if (capacity == 0) printf("capacity is 0!!!\n");
         return *this;
       }


       virtual ~SendBuffer() {}

        void buffer(const T& val) {
          if (THREAD_SAFE)
            std::unique_lock<std::mutex> lock(mutex);

          // careful.  when growing, the content is automatically copied to new and destroyed in old.
          //   the destruction could cause double free error or segv if T does not have the move constructor/assignemnt operator.
//          if (capacity < 1) printf("capacity 0!!!\n");
//          if ((buf.size() >= capacity) && ((buf.size() - capacity) % 100000 == 0)) printf("size:  %lu\n", buf.size());
          buf.push_back(val);
        }

        void buffer(T&& val) {
          if (THREAD_SAFE)
            std::unique_lock<std::mutex> lock(mutex);

          // careful.  when growing, the content is automatically copied to new and destroyed in old.
          //   the destruction could cause double free error or segv if T does not have the move constructor/assignemnt operator.
//          if (capacity < 1) printf("capacity 0!!!\n");
//          if ((buf.size() >= capacity) && ((buf.size() - capacity) % 100000 == 0)) printf("size:  %lu\n", buf.size());
          buf.push_back(val);
        }

        size_t size() {
          if (THREAD_SAFE)
            std::unique_lock<std::mutex> lock(mutex);
          return buf.size();
        }

        bool isFull() {
          if (THREAD_SAFE)
            std::unique_lock<std::mutex> lock(mutex);
          return buf.size() >= capacity;
        }


        /**
         * moves the internal data to a new vector, and return that vector in the ExportType instance.
         * this readies the buffer for new data.
         * @return
         */
        ExportType exportData() {
          if (THREAD_SAFE)
            std::unique_lock<std::mutex> lock(mutex);

          //printf("exporting data. %lu\n", buf.size());
          std::vector<T> temp;
          temp.reserve(buf.size());
          std::swap(buf, temp);
          assert(buf.size() == 0);
          int trank = targetRank;
          return ExportType(std::move(trank), std::move(temp));
        }
    };

  } /* namespace io */
} /* namespace bliss */

#endif /* SENDBUFFER_HPP_ */
