/**
 * @file		buffer_pool.hpp
 * @ingroup bliss::io
 * @author	tpan
 * @brief   this file defines in-memory Buffer Pool.
 * @details
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef BUFFERPOOL_HPP_
#define BUFFERPOOL_HPP_

#include <cassert>

#include <mutex>
#include <deque>
#include <stdexcept>

#include "utils/logging.h"
#include "io/buffer.hpp"
#include "concurrent/concurrent.hpp"



namespace bliss
{
  namespace io
  {
    // TODO: move constructor and assignment operators between BufferPools of different thread safeties.
    // TODO: use threadsafe_queue where appropriate.

    /**
     * @class     BufferPool
     * @brief     BufferPool manages a set of reusable, in-memory buffers.  In particular, it can limit the amount of memory usage.
     * @details   The class is templated to provide thread-safe or unsafe BufferPools, managing thread-safe or unsafe Buffers
     *
     *            Each buffer is a block of preallocated memory that can be appended into.
     *
     *            The caller can acquire a new buffer from the pool and release an old buffer back into the pool, by Id.
     *            Released buffers are marked as empty and available.
     *
     *            When the pool is exhausted, Acquire function returns false.
     *
     *      life cycle:  pool acquire():  buffer unblock().  -> allow rw.
     *              application:  buffer block() -> allow ro
     *              pool release():  buffer clear().
     *
     *
     * @note      Each buffer should be acquired by a single thread  and released by a single thread.   (NOT REASONABLE)
     *            Note that Buffer does not actualy clear the memory, just resets the size, which is also the pointer for insertion.
     *            While the release method handles duplicate releases, there may be race conditions when multithreading, including
     *              likely loss of data (thread 1 appends data and thread 2 releases  it.
     *
     *            the pools passes ownership of the buffers to callers.  This is so to reduce race conditions when multiple
     *            threads are using and acquiring buffers.
     *
     *
     *
     *
     *  @tparam PoolTS    The thread safety property for the pool
     *  @tparam BufferTS  The thread safety property of each Buffer.
     */
    template<bliss::concurrent::ThreadSafety PoolTS, bliss::concurrent::ThreadSafety BufferTS = PoolTS>
    class BufferPool
    {
      public:
        /**
         * @brief     Type of Buffer in the BufferPool (thread safe or not)
         */
        using BufferType = bliss::io::Buffer<BufferTS>;
        /**
         * @brief     Index Type used to reference the Buffers in the BufferPool.
         */
        using BufferPtrType = std::unique_ptr<BufferType>;

      protected:
        /**
         * @brief     capacity of the BufferPool (in number of Buffers)
         */
        mutable int64_t capacity;

        /**
         * @brief     capacity of the individual Buffers (in bytes)
         */
        mutable size_t buffer_capacity;

        /**
         * @brief     current pool size
         */
        mutable typename std::conditional<PoolTS, std::atomic<int64_t>, int64_t>::type numBuffersAvailable;

        /**
         * @brief     Internal Set of available Buffers for immediate use.
         * @details   Using set instead of ThreadSafeQueue to ensure uniqueness of Buffer Ids.
         */
        std::deque<BufferPtrType>                          available;

        /**
         * @brief     mutex to control access.
         */
        mutable std::mutex mutex;


        /**
         * @brief     Private move constructor with mutex lock.
         * @param other   Source BufferPool object to move
         * @param l       Mutex lock on the Source BufferPool object
         */
        BufferPool(BufferPool<PoolTS, BufferTS>&& other, const std::lock_guard<std::mutex>&) :
          capacity(other.capacity),
          buffer_capacity(other.buffer_capacity),
          numBuffersAvailable(other.numBuffersAvailable),
          available(std::move(other.available)) {

          other.capacity = 0;
          other.numBuffersAvailable = 0;
          other.available.clear();
        };


      public:
        /**
         * @brief     construct a Buffer Pool with buffers of capacity _buffer_capacity.  the number of buffers in the pool is set to _pool_capacity.
         *
         * @param _pool_capacity       number of buffers in the pool
         * @param _buffer_capacity     size of the individual buffers
         */
        BufferPool(const int64_t _pool_capacity, const size_t _buffer_capacity) :
          capacity(_pool_capacity),
          buffer_capacity(_buffer_capacity),
          numBuffersAvailable(_pool_capacity), available()
          {};

        /**
         * @brief     default constructor is deleted.
         */
        BufferPool() = delete;

        /**
         * @brief     construct a Buffer pool with buffers of capacity _buffer_capacity.  the number of buffers in the pool is unbounded.
         *
         * @param _buffer_capacity    size of the individual buffers
         */
        explicit BufferPool(const size_t _buffer_capacity) :
            BufferPool<PoolTS, BufferTS>(std::numeric_limits<int64_t>::max(), _buffer_capacity) {};

        /**
         * @brief     default copy constructor is deleted.
         * @param other   source BufferPool object to copy from.
         */
        explicit BufferPool(const BufferPool<PoolTS, BufferTS>& other) = delete;

        /**
         * @brief      Move constructor.  Delegates to the Move constructor that locks access on the source.
         * @param other    source BufferPool object to move.
         */
        explicit BufferPool(BufferPool<PoolTS, BufferTS>&& other) :
          BufferPool<PoolTS, BufferTS>(std::move(other), std::lock_guard<std::mutex>(other.mutex)) {};

        /**
         * @brief     default copy assignment operator is deleted.
         * @param other   source BufferPool object to copy from.
         * @return        self, with member variables copied from other.
         */
        BufferPool<PoolTS, BufferTS>& operator=(const BufferPool<PoolTS, BufferTS>& other) = delete;

        /**
         * @brief     move assignment operator.
         * @param other   source BufferPool object to move from.
         * @return        self, with member variables moved from other.
         */
        BufferPool<PoolTS, BufferTS>& operator=(BufferPool<PoolTS, BufferTS>&& other) {
          std::unique_lock<std::mutex> mylock(mutex, std::defer_lock),
                                        otherlock(other.mutex, std::defer_lock);
          std::lock(mylock, otherlock);

          capacity = other.capacity; other.capacity = 0;
          buffer_capacity = other.buffer_capacity; other.buffer_capacity = 0;
          numBuffersAvailable = int64_t(other.numBuffersAvailable); other.numBuffersAvailable = 0;
          available = std::move(other.available);
          return *this;
        };


        /**
         * @brief     default destructor
         */
        virtual ~BufferPool() {};

        /**
         * @brief     Current size of the BufferPool
         * @return  size, type IdType (aka int).
         */
        const int64_t getAvailableCount() const  {
          return int64_t(numBuffersAvailable);
        }

        /**
         * @brief     Current capacity of the BufferPool
         * @return    capacity, type IdTyype (aka int)
         */
        const int64_t getCapacity() const {
          return capacity;
        }

        inline const bool isUnlimited() const {
          return capacity == std::numeric_limits<int64_t>::max();
        }

        /**
         * @brief     Each buffer's maximum capacity.
         * @return  each buffer's maximum capacity.
         */
        const size_t getBufferCapacity() const {
          return buffer_capacity;
        }


        /**
         * @brief     Resets the entire BufferPool: all Buffers in pool are  marked as released and available.
         *
         * @note      This is not entirely thread safe.  the available set is cleared by a single thread, but other
         *            threads may be acquiring buffers while they are being released.
         *
         *            It is envisioned that this function should be called from a single thread.
         *            however, care must be taken to ensure that all other threads have completed their work, else data loss is likely.
         */
        void reset() {
          numBuffersAvailable = capacity;  // can only produce up to capacity, ever.  some will be stored back.
        }


        /**
         * @brief     Get the next available Buffer by id.  if none available, return false in the first argument of the pair.
         * @return    std::pair with bool and Id,  first indicate whether acquire was successful, second indicate the BufferId if successful.
         */
        BufferPtrType tryAcquireBuffer() {
          BufferPtrType ptr;  // default is a null ptr.
          if (!isUnlimited()) {  // not fix size

            // if multiple threads, the predecrement is atomic, so each thread has a valid value to compare to empty
            // the total number of threads to continue will be same as total number of values that pass the empty test.
            if (--numBuffersAvailable < 0) {  // predecrement is a.fetch_sub(1)-1.  compare result to 0.
              // not unlimited, and already >= capacity.  bufferPtr is now lost, but that's okay.
              // undo and return.
              numBuffersAvailable++;    // post increment is a.fetch_add(1)
              return ptr;
            }  // else not empty and already decremented it.
          } // else don't need to touch numBuffersAvailable.

          // now get or create
          std::lock_guard<std::mutex> lock(mutex);
          if (available.empty()) {
            // none available for reuse
            // but has room to allocate, so do it.
            ptr = std::move(BufferPtrType(new BufferType(buffer_capacity)));
          } else {
            // has available for reuse.
            ptr = std::move(available.front());
            available.pop_front();
          }
          ptr->clear();
          ptr->unblock();

          return std::move(ptr);
        }

        /**
         * @brief     Release a buffer back to pool, by id.  if id is incorrect, throw std exception.  clears buffer before release back into available.
         * @note      THREAD SAFE
         * @param bufferId    The id of the Buffer to be released.
         */
        template<bliss::concurrent::ThreadSafety TS = PoolTS>
        typename std::enable_if<TS == bliss::concurrent::THREAD_SAFE, bool>::type releaseBuffer(BufferPtrType&& ptr) {
          if (!isUnlimited()) {  // check against size_t max - thread safe, and ensures there is no overflow.

            // if multiple threads, the preincrement is atomic, so each thread has a valid value to compare to capacity.
            // the total number of threads to continue will be  same as total number of values that pass the capacity test.
            if (++numBuffersAvailable > capacity) {  // preincrement is a.fetch_add(1)+1.  compare result to capacity.
              // not unlimited, and already >= capacity.  bufferPtr is now lost, but that's okay.
              // undo and return.
              numBuffersAvailable--;    // post decrement is a.fetch_sub(1)
              return false;
            }  // else not at capacity and already incremented it.
          } // else don't need to touch numBuffersAvailable.

          assert(ptr->isBlocked());
          ptr->clear();

          // now store the buffer.  make sure push_back is done one thread at a time.
          std::lock_guard<std::mutex> lock(mutex);
          available.push_back(std::move(ptr));

          return true;

        }


        /**
         * @brief     Release a buffer back to pool, by id.  if id is incorrect, throw std exception.  clears buffer before release back into available.
         * @note      THREAD UNSAFE
         * @param bufferId    The id of the Buffer to be released.
         */
        template<bliss::concurrent::ThreadSafety TS = PoolTS>
        typename std::enable_if<TS == bliss::concurrent::THREAD_UNSAFE, bool>::type releaseBuffer(BufferPtrType&& ptr) {
          if (!isUnlimited()) {
            if (numBuffersAvailable >= capacity)
              // not unlimited, and already at capacity.  bufferPtr is now lost, but that's okay.
              return false;
            else {
              // not at capacity, so increment.
              ++numBuffersAvailable;
            }
          } // else don't need to touch numBuffersAvailable.

          assert(ptr->isBlocked());
          ptr->clear();

          available.push_back(std::move(ptr));

          return true;

          //DEBUGF("buffer id %d is available for reuse.", id);
        }

    };


  } /* namespace io */
} /* namespace bliss */

#endif /* BUFFERPOOL_HPP_ */
