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
#include <sstream>

#include "io/io_exception.hpp"
#include "io/buffer.hpp"
#include "concurrent/concurrent.hpp"
#include "concurrent/threadsafe_queue.hpp"

namespace bliss
{
  namespace io
  {

    template<bliss::concurrent::ThreadSafety PoolTS, bliss::concurrent::ThreadSafety BufferTS = PoolTS>
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
        typedef size_t                                    IdType;

      protected:

        size_t capacity;
        size_t buffer_capacity;

        bliss::concurrent::ThreadSafeQueue<size_t> available;

        std::vector<BufferType> buffers;

        mutable std::mutex mutex;

        bool fixedSize;

        /**
         * should be called from within a mutex locked block
         * @return
         */
        size_t createNewBuffer() {
          size_t bufferId = buffers.size();
          buffers.push_back(BufferType(buffer_capacity));

          return bufferId;
        }


      public:
        /**
         *  construct a Buffer Pool with buffers of capacity _buffer_capacity.  the number of buffers in the pool is set to _pool_capacity.
         *
         * @param _pool_capacity       number of buffers in the pool
         * @param _buffer_capacity    size of the individual buffers
         */
        BufferPool(const size_t _pool_capacity, const size_t _buffer_capacity) :
          capacity(_pool_capacity), buffer_capacity(_buffer_capacity), available(),   // ThreadSafeQueue is not size bound.
          buffers(), fixedSize(_pool_capacity != std::numeric_limits<size_t>::max()) {

          // get an estimated size first, so don't have to keep growing the vector
          size_t size_hint = fixedSize ? capacity : 128;

          // reserve the buffer, and configure.
          buffers.reserve(size_hint);
          for (size_t i = 0; i < size_hint; ++i) {
            buffers.push_back(BufferType(buffer_capacity));
            available.waitAndPush(i);
          }
        };

        /**
         * construct a Buffer pool with buffers of capacity _buffer_capacity.  the number of buffers in the pool is unbounded.
         *
         * @param _buffer_capacity    size of the individual buffers
         */
        explicit BufferPool(const size_t _buffer_capacity) :
            BufferPool<bliss::concurrent::THREAD_SAFE, BufferThreadSafety>(std::numeric_limits<size_t>::max(), _buffer_capacity) {};


        /**
         * all auto deleted.
         */
        virtual ~BufferPool() {};

        size_t getSize() {
          return buffers.size();
        }
        size_t getCapacity() {
          return (fixedSize ? capacity : getSize());
        }
        bool isFixedSize() {
          return fixedSize;
        }

        void reset() {
          std::lock_guard<std::mutex> lock(mutex);
          available.clear();
          for (size_t i = 0; i < buffers.size(); ++i) {
            available.tryPush(i);
            buffers[i].clear();
          }
        }

        /**
         *
         * @param bufferId
         * @return
         */
        bool tryAcquireBuffer(size_t & bufferId) {

          if (fixedSize) {
            // if there is something available in "available", get it.  else get false.
            return available.tryPop(bufferId);
          } else {
            // can grow, but first check for available.
            if (!available.tryPop(bufferId)) {
              // non available.  so allocate a new one and return.
              std::lock_guard<std::mutex> lock(mutex);
              bufferId = createNewBuffer();
            } // else already got one.
            return true;
          }
        }

        /**
         *
         * @param bufferId
         */
        void waitAndAcquireBuffer(size_t & bufferId) {
          if (fixedSize) {
            // have to wait for available
            available.waitAndPop(bufferId);

          } else {  // can grow, but first see if we can reuse one.
            // check if there are some available.
            if (!available.tryPop(bufferId)) {
              std::lock_guard<std::mutex> lock(mutex);
              bufferId = createNewBuffer();
            }
          }
        }

        /**
         * release a buffer back to pool, by id.  if id is incorrect, throw exception.
         * @param bufferId
         */
        void releaseBuffer(const size_t& bufferId) throw (bliss::io::IOException) {
          std::unique_lock<std::mutex> lock(mutex);
          if (bufferId >= buffers.size()) {
            std::stringstream ss;
            ss << "ERROR: BufferPool releasing buffer with id " << bufferId << ", buffer does not exist.";
            throw IOException(ss.str());
          }
          // insert the buffer into "available"
          lock.unlock();

          available.tryPush(bufferId);
        }

        /**
         * get the buffer by id.  if id is incorrect, throw exception.
         * @param bufferId
         * @return
         */
        const BufferType& operator[](const size_t& bufferId) throw (bliss::io::IOException) {
          if (bufferId >= buffers.size()) {
            if (fixedSize) {
              std::stringstream ss;
              ss << "ERROR: BufferPool get buffer with id " << bufferId << ", buffer does not exist.";
              throw IOException(ss.str());
            } else {
              std::unique_lock<std::mutex> lock(mutex);
              for (size_t i = buffers.size(); i <= bufferId; ++i) {
                createNewBuffer();
              }
            }
          }
          // insert the buffer into "available"

          return buffers[bufferId];
        }

        const BufferType& at(const size_t& bufferId) throw (bliss::io::IOException) {
          return this->operator[](bufferId);
        }
    };



    /**
     * @class     bliss::io::BufferPool
     * @brief     non-thread-safe version of buffer pool
     * @details   this can only contain non-thread-safe version of buffers. (non-thread-safe buffer pool with thread-safe buffers does not make sense.
     *
     */
    template<>
    class BufferPool<bliss::concurrent::THREAD_UNSAFE>
    {
      public:
        typedef bliss::io::Buffer<bliss::concurrent::THREAD_UNSAFE>     BufferType;
        typedef size_t                                                  IdType;
      protected:

        size_t capacity;
        size_t buffer_capacity;

        std::deque<size_t> available;
        std::vector<BufferType> buffers;

        bool fixedSize;

        size_t createNewBuffer() {
          size_t bufferId = buffers.size();
          buffers.push_back(BufferType(buffer_capacity));

          return bufferId;
        }


      public:
        /**
         *  construct a Buffer Pool with buffers of capacity _buffer_capacity.  the number of buffers in the pool is set to _pool_capacity.
         *
         * @param _pool_capacity       number of buffers in the pool
         * @param _buffer_capacity    size of the individual buffers
         */
        BufferPool(const size_t pool_capacity, const size_t _buffer_capacity) :
          capacity(pool_capacity), buffer_capacity(_buffer_capacity), available(),   // ThreadSafeQueue is not size bound.
          buffers(), fixedSize(pool_capacity != std::numeric_limits<size_t>::max()) {

          // get an estimated size first, so don't have to keep growing the vector
          size_t size_hint = fixedSize ? capacity : 128;

          // reserve the buffer, and configure.
          buffers.reserve(size_hint);
          for (size_t i = 0; i < size_hint; ++i) {
            buffers.push_back(BufferType(buffer_capacity));
            available.push_back(i);
          }
        };

        /**
         * construct a Buffer pool with buffers of capacity _buffer_capacity.  the number of buffers in the pool is unbounded.
         *
         * @param _buffer_capacity    size of the individual buffers
         */
        explicit BufferPool(const size_t _buffer_capacity) :
            BufferPool<bliss::concurrent::THREAD_UNSAFE>(std::numeric_limits<size_t>::max(), _buffer_capacity) {};

        /**
         * all auto deleted.
         */
        virtual ~BufferPool() {};

        size_t getSize() {
          return buffers.size();
        }
        size_t getCapacity() {
          return (fixedSize ? capacity : getSize());
        }
        bool isFixedSize() {
          return fixedSize;
        }


        void reset() {
          available.clear();
          for (size_t i = 0; i < buffers.size(); ++i) {
            available.push_back(i);
            buffers[i].clear();
          }
        }

        /**
         *
         * @param bufferId
         * @return
         */
        bool tryAcquireBuffer(size_t & bufferId) {

          if (available.empty()) {
            // if empty,...
            if (fixedSize) {
              // none available and fixed size.  return false.
              return false;
            } else {
              // if there is no capacity limit, create a new buffer, insert into vector, and return it.
              bufferId = createNewBuffer();
            }
          } else {

            // if there is something available in "available", get it
            bufferId = available.front();
            available.pop_front();
          }
          return true;

        }

        /**
         *  waitAndAcquireBuffer does not really make sense - if there is nothing available in a single thread, waiting
         *  for the same thread to update does not make sense.
         * @param bufferId
         */
//        void waitAndAcquireBuffer(size_t & bufferId) {
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

        void releaseBuffer(const size_t& bufferId) throw (bliss::io::IOException) {
          if (bufferId >= buffers.size()) {
            std::stringstream ss;
            ss << "ERROR: BufferPool releasing buffer with id " << bufferId << ", buffer does not exist.";
            throw IOException(ss.str());
          }

          // insert the buffer into "available"
          available.push_back(bufferId);
        }

        const BufferType& operator[](const size_t& bufferId) throw (bliss::io::IOException) {
          if (bufferId >= buffers.size()) {
            if (fixedSize) {
              std::stringstream ss;
              ss << "ERROR: BufferPool get buffer with id " << bufferId << ", buffer does not exist.";
              throw IOException(ss.str());
            } else {
              for (size_t i = buffers.size(); i <= bufferId; ++i) {
                createNewBuffer();
              }
            }
          }

          return buffers[bufferId];
        }

        const BufferType& at(const size_t& bufferId) throw (bliss::io::IOException) {
          return this->operator[](bufferId);
        }
    };



  } /* namespace io */
} /* namespace bliss */

#endif /* BUFFERPOOL_HPP_ */