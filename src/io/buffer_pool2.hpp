/**
 * @file		buffer_pool2.hpp
 * @ingroup bliss::io
 * @author	tpan
 * @brief   this file defines in-memory Buffer Pool.
 * @details
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef BUFFERPOOL2_HPP_
#define BUFFERPOOL2_HPP_

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
    // TODO: move constructor and assignment operators between BufferPool2s of different thread safeties.

    /**
     * @class     BufferPool2
     * @brief     BufferPool2 manages a set of reusable, in-memory buffers.  In particular, it can limit the amount of memory usage.
     * @details   The class is templated to provide thread-safe or unsafe BufferPool2s, managing thread-safe or unsafe Buffers
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
    class BufferPool2
    {
      public:
        /**
         * @brief     Type of Buffer in the BufferPool2 (thread safe or not)
         */
        using BufferType = bliss::io::Buffer<BufferTS>;
        /**
         * @brief     Index Type used to reference the Buffers in the BufferPool2.
         */
        using BufferPtrType = std::unique_ptr<BufferType>;

      protected:
        /**
         * @brief     capacity of the BufferPool2 (in number of Buffers)
         */
        mutable size_t capacity;

        /**
         * @brief     capacity of the individual Buffers (in bytes)
         */
        const size_t buffer_capacity;

        /**
         * @brief     current pool size
         */
        mutable typename std::conditional<PoolTS, std::atomic<size_t>, size_t>::type numBuffersAvailable;

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
         * @param other   Source BufferPool2 object to move
         * @param l       Mutex lock on the Source BufferPool2 object
         */
        BufferPool2(BufferPool2<PoolTS, BufferTS>&& other, const std::lock_guard<std::mutex>& l) :
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
        BufferPool2(const size_t _pool_capacity, const size_t _buffer_capacity) :
          capacity(_pool_capacity),
          buffer_capacity(_buffer_capacity),
          numBuffersAvailable(_pool_capacity), available()
          {};

        /**
         * @brief     default constructor is deleted.
         */
        BufferPool2() = delete;

        /**
         * @brief     construct a Buffer pool with buffers of capacity _buffer_capacity.  the number of buffers in the pool is unbounded.
         *
         * @param _buffer_capacity    size of the individual buffers
         */
        explicit BufferPool2(const size_t _buffer_capacity) :
            BufferPool2<PoolTS, BufferTS>(std::numeric_limits<size_t>::max(), _buffer_capacity) {};

        /**
         * @brief     default copy constructor is deleted.
         * @param other   source BufferPool2 object to copy from.
         */
        explicit BufferPool2(const BufferPool2<PoolTS, BufferTS>& other) = delete;

        /**
         * @brief      Move constructor.  Delegates to the Move constructor that locks access on the source.
         * @param other    source BufferPool2 object to move.
         */
        explicit BufferPool2(BufferPool2<PoolTS, BufferTS>&& other) :
          BufferPool2<PoolTS, BufferTS>(std::move(other), std::lock_guard<std::mutex>(other.mutex)) {};

        /**
         * @brief     default copy assignment operator is deleted.
         * @param other   source BufferPool2 object to copy from.
         * @return        self, with member variables copied from other.
         */
        BufferPool2<PoolTS, BufferTS>& operator=(const BufferPool2<PoolTS, BufferTS>& other) = delete;

        /**
         * @brief     move assignment operator.
         * @param other   source BufferPool2 object to move from.
         * @return        self, with member variables moved from other.
         */
        BufferPool2<PoolTS, BufferTS>& operator=(BufferPool2<PoolTS, BufferTS>&& other) {
          std::unique_lock<std::mutex> mylock(mutex, std::defer_lock),
                                        otherlock(other.mutex, std::defer_lock);
          std::lock(mylock, otherlock);

          capacity = other.capacity; other.capacity = 0;
          buffer_capacity = other.buffer_capacity; other.buffer_capacity = 0;
          numBuffersAvailable = other.numBuffersAvailable; other.numBuffersAvailable = 0;
          available = std::move(other.available);
          return *this;
        };


        /**
         * @brief     default destructor
         */
        virtual ~BufferPool2() {};

        /**
         * @brief     Current size of the BufferPool2
         * @return  size, type IdType (aka int).
         */
        const size_t getAvailableCount() const  {
        	if (isFixedSize())
        		return size_t(numBuffersAvailable);
        	else
        		return capacity;
        }

        /**
         * @brief     Current capacity of the BufferPool2
         * @return    capacity, type IdTyype (aka int)
         */
        const size_t getCapacity() const {
          return capacity;
        }

        /**
         * @brief     Each buffer's maximum capacity.
         * @return  each buffer's maximum capacity.
         */
        const size_t getBufferCapacity() const {
          return buffer_capacity;
        }

        /**
         * @brief     whether the BufferPool2 is growable
         * @return    bool indicating if BufferPool2 is growable.
         */
        inline const bool isFixedSize() const {
          return capacity < std::numeric_limits<size_t>::max();
        }

        /**
         * @brief     Resets the entire BufferPool2: all Buffers in pool are  marked as released and available.
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
          if (this->getAvailableCount() > 0) {
              if (this->isFixedSize()) numBuffersAvailable--;

        	  std::lock_guard lock(mutex);
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
          } // else, got all buffers allowed.
          return std::move(ptr);
        }

        /**
         * @brief     Release a buffer back to pool, by id.  if id is incorrect, throw std exception.  clears buffer before release back into available.
         * @note      THREAD SAFE
         * @param bufferId    The id of the Buffer to be released.
         */
        template<bliss::concurrent::ThreadSafety TS = PoolTS>
        typename std::enable_if<TS == bliss::concurrent::THREAD_SAFE, bool>::type releaseBuffer(BufferPtrType&& ptr) {
          std::lock_guard<std::mutex> lock(mutex);
          if (numBuffersAvailable < capacity) {
			  assert(ptr->isBlocked());
			  ptr->clear();
			  available.push_back(std::move(ptr));
			  numBuffersAvailable++;  // post increment on atomic is 1 instruction less, but has an extract alloc for non-atomic.
			  return true;
          } else return false;

        }


        /**
         * @brief     Release a buffer back to pool, by id.  if id is incorrect, throw std exception.  clears buffer before release back into available.
         * @note      THREAD UNSAFE
         * @param bufferId    The id of the Buffer to be released.
         */
        template<bliss::concurrent::ThreadSafety TS = PoolTS>
        typename std::enable_if<TS == bliss::concurrent::THREAD_UNSAFE, bool>::type releaseBuffer(BufferPtrType&& ptr) {
        	if (numBuffersAvailable < capacity) {
			  assert(ptr->isBlocked());
			  ptr->clear();
			  available.push_back(std::move(ptr));
			  ++numBuffersAvailable;
			  return true;
        	} else {
        		return false;
        	}

          //DEBUGF("buffer id %d is available for reuse.", id);
        }

    };


  } /* namespace io */
} /* namespace bliss */

#endif /* BUFFERPOOL2_HPP_ */
