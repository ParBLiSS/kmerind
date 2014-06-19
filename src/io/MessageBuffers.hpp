/**
 * @file		MessageBuffers.hpp
 * @ingroup
 * @author	tpan
 * @brief   buffering data for MPI send/receive
 * @details MessageBuffers provides buffering for MPI Send/Receive messages.
 *          Terminology:  Message is an MPI compatible message, serialized from buffer, including the appropriate metadata.
 *                        Element is a single piece of data from the application code.
 *                        Buffer is a collection of elements.
 *
 *
 *          provides functions for storing single and collection of elements into the buffer, and function to serialize
 *          the buffer into a message.
 *
 *          2 modes:  send:  application code populates the buffer either 1 element at a time or multiple elements at a time
 *                            then application code will request the message
 *                    recv:  application code gets a buffer region and receive from comm channel
 *                            application code will request the entire buffer.
 *
 *          each element is destined for a particular mpi target rank, along with a tag that denotes its type.  The add method uses a byte array along with size.
 *          buffer has a limited size.
 *
 *          This class contains a vector of buffers, each with a different id.  The tag is set as template parameter, since we need to know type corresponding to the
 *          tag at compile time anyways.
 *
 *          The class allows sharing between thread but it can certainly be used within a single thread.
 *
 *          send/recv have different semantics?  single class would be nice - easier to support different patterns
 *            buffer -> send = send
 *            recv -> buffer = recv
 *            recv -> send   = forward
 *            buffer -> buffer = local, same rank send.
 *
 *          event driven with observer pattern on full/receive with CommLayer as handler? - avoids looping to check if full, before enqueue.
 *            comm layer can do this work without being an observer since it has to have a send method - can check for full and send.
 *
 *
 *
 *
 *          buffer's lifetime is from the first call to the buffer function to the completion of the MPI send, or
 *              from the start of the MPI recv to the consumption of the buffer.
 *
 *          QUESTIONS: reuse buffer?  std::move or copy content?
 *          MAY NEED A BUFFER OF THESE.  Do that later.
 *
 *
 *
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef MESSAGEBUFFERS_HPP_
#define MESSAGEBUFFERS_HPP_

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
     * @class			bliss::io::MessageBuffers
     * @brief     a data structure to buffer and send MPI messages
     * @details   THREAD_SAFE enables std::mutex for OMP thread safety.  this likely introduces a good amount of contention.
     *            approach:
     *                have reference to downstream queue.
     *                when full, create a new vector, swap, enqueue std::pair with target rank and stdvector.
     *
     */
    template<int TAG, bool THREAD_SAFE=true>
    class MessageBuffers {
        // need some default MPI buffer size, then pack in sizeof(T) blocks as many times as possible.

      protected:

        // shared, available buffers
        std::queue<  > BufferPool;


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

        static MessageBuffers<TAG, THREAD_SAFE>& getBuffer() {

        }


       MessageBuffers(const int _targetId, const int _typeid, const size_t nbytes)
           : targetId(_targetId), typeId(_typeid) {
          // initialize internal buffers (double buffering)
          capacity = 3 * nbytes / (2 * sizeof(T));

          buf.reserve(capacity);

          printf("construct size, capacity: %lu, %lu. datastructure size %lu\n", buf.size(), buf.capacity(), sizeof(*this));
          if (capacity == 0) printf("capacity is 0!!!\n");
        }

       // copy constructor.
       MessageBuffers(const MessageBuffers<TAG, THREAD_SAFE>& other) :
         targetId(other.targetId), typeId(other.typeId), buf(other.buf), capacity(other.capacity) {
         //printf("copy size: %lu\n", buf.size());
         if (capacity == 0) printf("capacity is 0!!!\n");

       }

       // copy assign operator
       MessageBuffers<TAG, THREAD_SAFE>& operator=(const MessageBuffers<TAG, THREAD_SAFE>& other) {
         targetId = other.targetId;
         typeId = other.typeId;
         buf = other.buf;  // copy
         capacity = other.capacity;
         //printf("copy assign size: %lu\n", buf.size());
         if (capacity == 0) printf("capacity is 0!!!\n");
         return *this;
       }

       // move constructor.
       MessageBuffers(MessageBuffers<TAG, THREAD_SAFE>&& other) :
         targetId(other.targetId), typeId(other.typeId),  buf(std::move(other.buf)), capacity(other.capacity) {
         //printf("move size: %lu\n", buf.size());
         if (capacity == 0) printf("capacity is 0!!!\n");

       }

       // copy assign operator
       MessageBuffers<TAG, THREAD_SAFE>& operator=(MessageBuffers<TAG, THREAD_SAFE>&& other) {
         targetId = other.targetId;
         typeId = other.typeId;
         buf = std::move(other.buf);  // copy
         capacity = other.capacity;
         //printf("move assign size: %lu\n", buf.size());
         if (capacity == 0) printf("capacity is 0!!!\n");
         return *this;
       }


       virtual ~MessageBuffers() {}


        void buffer(const T& val) {
          if (THREAD_SAFE) {
            std::unique_lock<std::mutex> lock(mutex);

          // careful.  when growing, the content is automatically copied to new and destroyed in old.
          //   the destruction could cause double free error or segv if T does not have the move constructor/assignemnt operator.
//          if (capacity < 1) printf("capacity 0!!!\n");
//          if ((buf.size() >= capacity) && ((buf.size() - capacity) % 100000 == 0)) printf("size:  %lu\n", buf.size());
            buf.push_back(val);
          } else
            buf.push_back(val);
        }

        void buffer(T&& val) {
          if (THREAD_SAFE) {
            std::unique_lock<std::mutex> lock(mutex);

          // careful.  when growing, the content is automatically copied to new and destroyed in old.
          //   the destruction could cause double free error or segv if T does not have the move constructor/assignemnt operator.
//          if (capacity < 1) printf("capacity 0!!!\n");
//          if ((buf.size() >= capacity) && ((buf.size() - capacity) % 100000 == 0)) printf("size:  %lu\n", buf.size());
            buf.push_back(val);
          } else
            buf.push_back(val);
        }

        size_t size() {
          if (THREAD_SAFE) {
            std::unique_lock<std::mutex> lock(mutex);
            return buf.size();
          } else
            return buf.size();

        }

        bool isFull() {
          if (THREAD_SAFE) {
            std::unique_lock<std::mutex> lock(mutex);
            return buf.size() >= capacity;
          } else
            return buf.size() >= capacity;
        }


        /**
         * moves the internal data to a new vector, and return that vector in the ExportType instance.
         * this readies the buffer for new data.
         * @return
         */
        ExportType exportData() {
          std::vector<T> temp;
          int trank = targetId;
          if (THREAD_SAFE) {
            std::unique_lock<std::mutex> lock(mutex);

            //printf("exporting data. %lu\n", buf.size());
            temp.reserve(buf.size());
            std::swap(buf, temp);
            assert(buf.size() == 0);
          } else {
            //printf("exporting data. %lu\n", buf.size());
            temp.reserve(buf.size());
            std::swap(buf, temp);
            assert(buf.size() == 0);
          }
          return ExportType(std::move(trank), std::move(temp));
        }
    };


  } /* namespace io */
} /* namespace bliss */

#endif /* MESSAGEBUFFERS_HPP_ */
