/**
 * @file		buffer_pool.hpp
 * @ingroup
 * @author	tpan
 * @brief   a pool of memory buffers.
 * @details provides a reusable pool of buffers.  Each buffer is a block of preallocated memory that can be copied into.
 *          this dynamically grows if capacity is not specified.
 *
 *          the pools RETAINS OWNERSHIP of the buffers.
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef BUFFERPOOL_HPP_
#define BUFFERPOOL_HPP_


#include <mutex>
#include <deque>
#include <vector>

#include "io/buffer.hpp"
#include "concurrent/concurrent.hpp"
#include "concurrent/threadsafe_queue.hpp"

namespace bliss
{
  namespace io
  {

    template<bliss::concurrent::ThreadSafety PoolTS, bliss::concurrent::ThreadSafety BufferThreadSafety>
    class BufferPool;


    /**
     * @class			bliss::io::BufferPool
     * @brief     thread-safe version of buffer pool
     * @details   this can contain thread-safe buffers or non-thread safe buffers.
     *            user of the pool will have to access the pool via index.
     */
    template<bliss::concurrent::ThreadSafety BufferThreadSafety>
    class BufferPool<bliss::concurrent::THREAD_SAFE, BufferThreadSafety>
    {
      public:
        typedef bliss::io::Buffer<BufferThreadSafety>     BufferType;

      protected:

        size_t capacity;
        size_t buffer_capacity;

        bliss::concurrent::ThreadSafeQueue<int> available;
        std::vector<BufferType> buffers;

        mutable std::mutex mutex;

        bool fixedSize;

        int createNewBuffer() {
          int bufferId;
          std::unique_lock<std::mutex> lock(mutex);
          bufferId = buffers.size();
          buffers.push_back(BufferType(buffer_capacity));
          lock.unlock();

          available.tryPush(bufferId);    // available is not fixedSized, so will just succeed.
          return bufferId;
        }


      public:
        /**
         * if capacity is set to size_t max value then we populate the array with 100.  else the _capacity number of entries are created.
         * @param _capacity
         */
        BufferPool(const size_t _buffer_capacity, const size_t pool_capacity = std::numeric_limits<size_t>::max()) :
          capacity(pool_capacity), buffer_capacity(_buffer_capacity), available(),   // ThreadSafeQueue is not size bound.
          buffers(), fixedSize(pool_capacity == std::numeric_limits<size_t>::max()) {

          // get an estimated size first, so don't have to keep growing the vector
          size_t size_hint = fixedSize ? 128 : capacity;

          // reserve the buffer, and configure.
          buffers.reserve(size_hint);
          for (int i = 0; i < size_hint; ++i) {
            buffers.push_back(BufferType(buffer_capacity));
            available.waitAndPush(i);
          }
        };

        /**
         * all auto deleted.
         */
        virtual ~BufferPool() {};


        /**
         *
         * @param bufferId
         * @return
         */
        bool tryAcquireBuffer(int & bufferId) {

          if (fixedSize) {
            // if there is something available in "available", get it.  else get false.
            return available.tryPop(bufferId);
          } else {
            // can grow, but first check for available.
            if (!available.tryPop(bufferId)) {
              // non available.  so allocate a new one and return.
              bufferId = createNewBuffer();
              return true;
            } else {
              //got one
              return true;
            }
          }
        }

        /**
         *
         * @param bufferId
         */
        void waitAndAcquireBuffer(int & bufferId) {
          if (fixedSize) {
            // have to wait for available
            available.waitAndPop(bufferId);

          } else {  // can grow, but first see if we can reuse one.
            // check if there are some available.
            if (!available.tryPop(bufferId)) {
              bufferId = createNewBuffer();
            }
          }
        }

        void releaseBuffer(const int& bufferId) {
          // insert the buffer into "available"
          available.tryPush(bufferId);
        }

        const BufferType& getBuffer(const int& bufferId) {
          return buffers[bufferId];
        }

    };



    /**
     * @class     bliss::io::BufferPool
     * @brief     non-thread-safe version of buffer pool
     * @details   this can only contain non-thread-safe version of buffers. (non-thread-safe buffer pool with thread-safe buffers does not make sense.
     *
     */
    template<>
    class BufferPool<bliss::concurrent::THREAD_UNSAFE, bliss::concurrent::THREAD_UNSAFE>
    {
      public:
        typedef bliss::io::Buffer<bliss::concurrent::THREAD_UNSAFE>     BufferType;

      protected:

        size_t capacity;
        size_t buffer_capacity;

        std::deque<int> available;
        std::vector<BufferType> buffers;

        bool fixedSize;

        int createNewBuffer() {
          int bufferId = buffers.size();
          buffers.push_back(BufferType(buffer_capacity));

          available.push_back(bufferId);    // available is not fixedSized, so will just succeed.
          return bufferId;
        }


      public:
        /**
         * if capacity is set to size_t max value then we populate the array with 100.  else the _capacity number of entries are created.
         * @param _capacity
         */
        BufferPool(const size_t _buffer_capacity, const size_t pool_capacity = std::numeric_limits<size_t>::max()) :
          capacity(pool_capacity), buffer_capacity(_buffer_capacity), available(),   // ThreadSafeQueue is not size bound.
          buffers(), fixedSize(pool_capacity == std::numeric_limits<size_t>::max()) {

          // get an estimated size first, so don't have to keep growing the vector
          size_t size_hint = fixedSize ? 128 : capacity;

          // reserve the buffer, and configure.
          buffers.reserve(size_hint);
          for (int i = 0; i < size_hint; ++i) {
            buffers.push_back(BufferType(buffer_capacity));
            available.push_back(i);
          }
        };

        /**
         * all auto deleted.
         */
        virtual ~BufferPool() {};


        /**
         *
         * @param bufferId
         * @return
         */
        bool tryAcquireBuffer(int & bufferId) {

          if (available.empty()) {
            // if empty,...
            if (fixedSize) {
              // none available and fixed size.  return false.
              return false;
            } else {
              // if there is no capacity limit, create a new buffer, insert into vector, and return it.
              bufferId = createNewBuffer();
              return true;
            }
          } else {
            // if there is something available in "available", get it
            bufferId = available.front();
            available.pop_front();
            return true;
          }
        }

        /**
         *  waitAndAcquireBuffer does not really make sense - if there is nothing available in a single thread, waiting
         *  for the same thread to update does not make sense.
         * @param bufferId
         */
//        void waitAndAcquireBuffer(int & bufferId) {
//          // if empty, and can grow, then grow.
//          if (available.empty()) {
//            if (!fixedSize) {
//              bufferId = createNewBuffer();
//            }
//          }
//
//          // now wait until available is n
//          if (fixedSize) {
//            // have to wait for available
//            available.waitAndPop(bufferId);
//
//          } else {
//            // check if there are some available.
//            if (!available.tryPop(bufferId)) {
//
//            }
//          }
//        }

        void releaseBuffer(const int& bufferId) {
          // insert the buffer into "available"
          available.push_back(bufferId);
        }

        const BufferType& getBuffer(const int& bufferId) {
          return buffers[bufferId];
        }
    };



  } /* namespace io */
} /* namespace bliss */

#endif /* BUFFERPOOL_HPP_ */
