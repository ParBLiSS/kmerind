/**
 * @file		MessageBuffer.hpp
 * @ingroup
 * @author	tpan
 * @brief   buffering data for MPI send/receive
 * @details MessageBuffer provides buffering for MPI Send/Receive messages.  provides functions for storing single and collection of data into the buffer.
 *          assumption is also that MPI receive will populate a whole buffer, so retrieve is for the entire buffer.
 *
 *          each element is destined for a particular mpi target rank, along with a tag that denotes its type.  The add method uses a byte array along with size.
 *
 *          buffer has a limited size.
 *
 *          buffer's lifetime is from the first call to the buffer function to the completion of the MPI send, or
 *              from the start of the MPI recv to the consumtion of the buffer.
 *
 *          QUESTIONS: reuse buffer?  std::move or copy content?
 *
 *          send/recv have different semantics?
 *
 *          event driven with observer pattern on full/receive with CommLayer as handler? - avoids looping to check if full, before enqueue.
 *            comm layer can do this work since it has to have a send method - can check for full and send.
 *
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef MESSAGEBUFFER_HPP_
#define MESSAGEBUFFER_HPP_

#include <cassert>
#include <unistd.h>  // for usleep

#include <thread>
#include <mutex>
#include <vector>
#include <utility>

#include "config.hpp"

namespace bliss
{
  namespace io
  {

    /**
     * @class			bliss::io::MessageBuffer
     * @brief     a data structure to buffer and send MPI messages
     * @details   THREAD_SAFE enables std::mutex for OMP thread safety.  this likely introduces a good amount of contention.
     *            approach:
     *                have reference to downstream queue.
     *                when full, create a new vector, swap, enqueue std::pair with target rank and stdvector.
     *
     */
    template<typename T, bool THREAD_SAFE=true>
    class MessageBuffer {
        // need some default MPI buffer size, then pack in sizeof(T) blocks as many times as possible.

      protected:
        // rank of target entry.
        int targetId;
        int typeId;

        // 2 buffers per target
        std::vector<T> buf;
        size_t capacity;

        mutable std::mutex mutex;

      public:
        typedef T ValueType;
        typedef std::pair<int, std::vector<T> > ExportType;

       MessageBuffer(const int _targetId, const int _typeid, const size_t nbytes)
           : targetId(_targetId), typeId(_typeid) {
          // initialize internal buffers (double buffering)
          capacity = 3 * nbytes / (2 * sizeof(T));

          buf.reserve(capacity);

          printf("construct size, capacity: %lu, %lu. datastructure size %lu\n", buf.size(), buf.capacity(), sizeof(*this));
          if (capacity == 0) printf("capacity is 0!!!\n");
        }

       // copy constructor.
       MessageBuffer(const MessageBuffer<T, THREAD_SAFE>& other) :
         targetId(other.targetId), typeId(other.typeId), buf(other.buf), capacity(other.capacity) {
         //printf("copy size: %lu\n", buf.size());
         if (capacity == 0) printf("capacity is 0!!!\n");

       }

       // copy assign operator
       MessageBuffer<T, THREAD_SAFE>& operator=(const MessageBuffer<T, THREAD_SAFE>& other) {
         targetId = other.targetId;
         typeId = other.typeId;
         buf = other.buf;  // copy
         capacity = other.capacity;
         //printf("copy assign size: %lu\n", buf.size());
         if (capacity == 0) printf("capacity is 0!!!\n");
         return *this;
       }

       // move constructor.
       MessageBuffer(MessageBuffer<T, THREAD_SAFE>&& other) :
         targetId(other.targetId), typeId(other.typeId),  buf(std::move(other.buf)), capacity(other.capacity) {
         //printf("move size: %lu\n", buf.size());
         if (capacity == 0) printf("capacity is 0!!!\n");

       }

       // copy assign operator
       MessageBuffer<T, THREAD_SAFE>& operator=(MessageBuffer<T, THREAD_SAFE>&& other) {
         targetId = other.targetId;
         typeId = other.typeId;
         buf = std::move(other.buf);  // copy
         capacity = other.capacity;
         //printf("move assign size: %lu\n", buf.size());
         if (capacity == 0) printf("capacity is 0!!!\n");
         return *this;
       }


       virtual ~MessageBuffer() {}

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
          int trank = targetId;
          return ExportType(std::move(trank), std::move(temp));
        }
    };

  } /* namespace io */
} /* namespace bliss */

#endif /* MESSAGEBUFFER_HPP_ */
