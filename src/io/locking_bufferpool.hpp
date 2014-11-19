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
#include <unordered_set>
#include <stdexcept>

#include "utils/logging.h"
#include "io/locking_buffer.hpp"
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
     *  @tparam PoolLT    The thread safety property for the pool
     *  @tparam BufferLT  The thread safety property of each Buffer.
     */
    template<bliss::concurrent::LockType PoolLT, bliss::concurrent::LockType BufferLT = PoolLT>
    class BufferPool
    {
      public:
        /**
         * @brief     Type of Buffer in the BufferPool (thread safe or not)
         */
        using BufferType = bliss::io::Buffer<BufferLT>;
        /**
         * @brief     Index Type used to reference the Buffers in the BufferPool.
         */
        using BufferPtrType = std::shared_ptr<BufferType>;   // shared pointer allows atomic operations as well as check for expired pointers.

        static const bliss::concurrent::LockType poolLT = PoolLT;
        static const bliss::concurrent::LockType bufferLT = BufferLT;


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
        mutable typename std::conditional<PoolLT == bliss::concurrent::LockType::LOCKFREE, std::atomic<int64_t>, int64_t>::type numBuffersAvailable;

        /**
         * @brief     Internal Set of available Buffers for immediate use.
         * @details   Using set instead of ThreadSafeQueue to ensure uniqueness of Buffer Ids.
         */
        std::deque<BufferPtrType>                          available;
        std::unordered_set<BufferPtrType>                  in_use;

        /**
         * @brief     mutex to control access.
         */
        mutable std::mutex mutex;
        mutable std::atomic_flag spinlock = ATOMIC_FLAG_INIT;


        /**
         * @brief     Private move constructor with mutex lock.
         * @param other   Source BufferPool object to move
         * @param l       Mutex lock on the Source BufferPool object
         */
        BufferPool(BufferPool<PoolLT, BufferLT>&& other, const std::lock_guard<std::mutex>&) :
          capacity(other.capacity),
          buffer_capacity(other.buffer_capacity),
          numBuffersAvailable((int64_t)(other.numBuffersAvailable)),
          available(std::move(other.available)),
          in_use(std::move(other.in_use())){

          other.capacity = 0;
          other.numBuffersAvailable = 0;
          other.available.clear();
          other.in_use.clear();
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
          numBuffersAvailable(_pool_capacity),
          available(), in_use()
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
            BufferPool<PoolLT, BufferLT>(std::numeric_limits<int64_t>::max(), _buffer_capacity) {};

        /**
         * @brief     default copy constructor is deleted.
         * @param other   source BufferPool object to copy from.
         */
        explicit BufferPool(const BufferPool<PoolLT, BufferLT>& other) = delete;

        /**
         * @brief      Move constructor.  Delegates to the Move constructor that locks access on the source.
         * @param other    source BufferPool object to move.
         */
        explicit BufferPool(BufferPool<PoolLT, BufferLT>&& other) :
          BufferPool<PoolLT, BufferLT>(std::move(other), std::lock_guard<std::mutex>(other.mutex)) {};

        /**
         * @brief     default copy assignment operator is deleted.
         * @param other   source BufferPool object to copy from.
         * @return        self, with member variables copied from other.
         */
        BufferPool<PoolLT, BufferLT>& operator=(const BufferPool<PoolLT, BufferLT>& other) = delete;

        /**
         * @brief     move assignment operator.
         * @param other   source BufferPool object to move from.
         * @return        self, with member variables moved from other.
         */
        BufferPool<PoolLT, BufferLT>& operator=(BufferPool<PoolLT, BufferLT>&& other) {
          std::unique_lock<std::mutex> mylock(mutex, std::defer_lock),
                                        otherlock(other.mutex, std::defer_lock);
          std::lock(mylock, otherlock);

          capacity = other.capacity; other.capacity = 0;
          buffer_capacity = other.buffer_capacity; other.buffer_capacity = 0;
          numBuffersAvailable = (int64_t)(other.numBuffersAvailable); other.numBuffersAvailable = 0;
          available.clear();  available.swap(other.available);
          in_use.clear();  in_use.swap(other.in_use);

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
          return (int64_t)(numBuffersAvailable);   // load is atomic
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
          numBuffersAvailable = capacity;  // store is atomic .
        }


        /**
         * @brief     Get the next available Buffer by id.  if none available, return false in the first argument of the pair.
         * @return    std::pair with bool and Id,  first indicate whether acquire was successful, second indicate the BufferId if successful.
         */
        template<bliss::concurrent::LockType LT = PoolLT>
                typename std::enable_if<LT == bliss::concurrent::LockType::MUTEX, BufferPtrType>::type tryAcquireBuffer() {
          BufferPtrType sptr;  // default is a null ptr.

          std::unique_lock<std::mutex> lock(mutex);

          if (!isUnlimited()) {
            if (numBuffersAvailable <= 0) {  // not fix size, and no more available
              lock.unlock();
              return sptr;
            } else {
              numBuffersAvailable--;
            }
          }

          // now get or create
          if (available.empty()) {
            // none available for reuse
            // but has room to allocate, so do it.
            sptr.reset(new BufferType(buffer_capacity));
          } else {
            // has available for reuse.
            sptr = available.front();
            available.pop_front();
          }

          // get the buffer ready for use.
          sptr->clear();
          sptr->unblock();
          in_use.insert(sptr);  // store the shared pointer.

          lock.unlock();

          return std::move(sptr);
        }

        /**
         * @brief     Get the next available Buffer by id.  if none available, return false in the first argument of the pair.
         * @return    std::pair with bool and Id,  first indicate whether acquire was successful, second indicate the BufferId if successful.
         */
        template<bliss::concurrent::LockType LT = PoolLT>
                typename std::enable_if<LT == bliss::concurrent::LockType::SPINLOCK, BufferPtrType>::type tryAcquireBuffer() {
          BufferPtrType sptr;  // default is a null ptr.

          while (spinlock.test_and_set());

          if (!isUnlimited()) {
            if (numBuffersAvailable <= 0) {  // not fix size, and no more available
              spinlock.clear();
              return sptr;
            } else {
              numBuffersAvailable--;
            }
          }

          // now get or create
          if (available.empty()) {
            // none available for reuse
            // but has room to allocate, so do it.
            sptr.reset(new BufferType(buffer_capacity));
          } else {
            // has available for reuse.
            sptr = available.front();
            available.pop_front();
          }


          // get the buffer ready for use.
          sptr->clear();
          sptr->unblock();
          in_use.insert(sptr);  // store the shared pointer.
          spinlock.clear();

          return std::move(sptr);
        }
//
//        template<bliss::concurrent::LockType LT = PoolLT>
//                typename std::enable_if<LT == bliss::concurrent::LockType::LOCKFREE, BufferPtrRefType>::type tryAcquireBuffer() {
//          static_assert(false, "BufferPool currently does not support LockFree operations due to internal use of deque");
//        }
//


        template<bliss::concurrent::LockType LT = PoolLT>
               typename std::enable_if<LT == bliss::concurrent::LockType::NONE, BufferPtrType>::type tryAcquireBuffer() {
          BufferPtrType sptr;  // default is a null ptr.

          if (!isUnlimited()) {
            if (numBuffersAvailable <= 0) {  // not fix size, and no more available
              return sptr;
            } else {
              numBuffersAvailable--;
            }
          }

          // now get or create
          if (available.empty()) {
            // none available for reuse
            // but has room to allocate, so do it.
            sptr.reset(new BufferType(buffer_capacity));
          } else {
            // has available for reuse.
            sptr = available.front();
            available.pop_front();
          }

          // get the buffer ready for use.
          sptr->clear();
          sptr->unblock();

          in_use.insert(sptr);  // store the shared pointer.

          return std::move(sptr);
        }

        /**
         * @brief     Release a buffer back to pool, by id.  if id is incorrect, throw std exception.  clears buffer before release back into available.
         * @note      THREAD SAFE
         * @param ptr weak_ptr to buffer.
         */
        template<bliss::concurrent::LockType LT = PoolLT>
        typename std::enable_if<LT == bliss::concurrent::LockType::MUTEX, bool>::type releaseBuffer(BufferPtrType&& ptr) {

          if (! ptr) return false;

          std::unique_lock<std::mutex> lock(mutex);
          if (!isUnlimited()) {  // check against size_t max - thread safe, and ensures there is no overflow.

            // if multiple threads, the preincrement is atomic, so each thread has a valid value to compare to capacity.
            // the total number of threads to continue will be  same as total number of values that pass the capacity test.
            if (numBuffersAvailable >= capacity) {  // preincrement is a.fetch_add(1)+1.  compare result to capacity.
              // not unlimited, and already >= capacity.  bufferPtr is now lost, but that's okay.
              // undo and return.
              lock.unlock();
              return false;
            } else {
              // else not at capacity and already incremented it.
              numBuffersAvailable++;
            }
          } // else don't need to touch numBuffersAvailable.

          // now store the buffer.  make sure push_back is done one thread at a time.
          in_use.erase(ptr);
          available.push_back(std::move(ptr));
          lock.unlock();

          return true;

        }


        /**
         * @brief     Release a buffer back to pool, by id.  if id is incorrect, throw std exception.  clears buffer before release back into available.
         * @note      THREAD SAFE
         * @param bufferId    The id of the Buffer to be released.
         */
        template<bliss::concurrent::LockType LT = PoolLT>
        typename std::enable_if<LT == bliss::concurrent::LockType::SPINLOCK, bool>::type releaseBuffer(BufferPtrType&& ptr) {
          if (! ptr) return false;

          while (spinlock.test_and_set());
          if (!isUnlimited()) {  // check against size_t max - thread safe, and ensures there is no overflow.

            // if multiple threads, the preincrement is atomic, so each thread has a valid value to compare to capacity.
            // the total number of threads to continue will be  same as total number of values that pass the capacity test.
            if (numBuffersAvailable >= capacity) {  // preincrement is a.fetch_add(1)+1.  compare result to capacity.
              // not unlimited, and already >= capacity.  bufferPtr is now lost, but that's okay.
              // undo and return.
              spinlock.clear();
              return false;
            } else {
              // else not at capacity and already incremented it.
              numBuffersAvailable++;
            }
          } // else don't need to touch numBuffersAvailable.

          // now store the buffer.  make sure push_back is done one thread at a time.
          in_use.erase(ptr);
          available.push_back(std::move(ptr));
          spinlock.clear();

          return true;

        }
//
//        /**
//         * @brief     Release a buffer back to pool, by id.  if id is incorrect, throw std exception.  clears buffer before release back into available.
//         * @note      THREAD SAFE
//         * @param bufferId    The id of the Buffer to be released.
//         */
//        template<bliss::concurrent::LockType LT = PoolLT>
//        typename std::enable_if<LT == bliss::concurrent::LockType::LOCKFREE, bool>::type releaseBuffer(BufferPtrRefType&& ptr) {
//          static_assert(false, "BufferPool currently does not support LockFree operations due to internal use of deque");
//        }


        /**
         * @brief     Release a buffer back to pool, by id.  if id is incorrect, throw std exception.  clears buffer before release back into available.
         * @note      THREAD UNSAFE
         * @param bufferId    The id of the Buffer to be released.
         */
        template<bliss::concurrent::LockType LT = PoolLT>
        typename std::enable_if<LT == bliss::concurrent::LockType::NONE, bool>::type releaseBuffer(BufferPtrType&& ptr) {
           if (! ptr) return false;

          if (!isUnlimited()) {  // check against size_t max - thread safe, and ensures there is no overflow.

            // if multiple threads, the preincrement is atomic, so each thread has a valid value to compare to capacity.
            // the total number of threads to continue will be  same as total number of values that pass the capacity test.
            if (numBuffersAvailable >= capacity) {  // preincrement is a.fetch_add(1)+1.  compare result to capacity.
              // not unlimited, and already >= capacity.  bufferPtr is now lost, but that's okay.
              // undo and return.
              return false;
            } else {
              // else not at capacity and already incremented it.
              numBuffersAvailable++;
            }
          } // else don't need to touch numBuffersAvailable.

          // now store the buffer.  make sure push_back is done one thread at a time.
          in_use.erase(ptr);
          available.push_back(std::move(ptr));

          return true;
        }

    };

    template<bliss::concurrent::LockType PoolLT, bliss::concurrent::LockType BufferLT>
    const bliss::concurrent::LockType BufferPool<PoolLT, BufferLT>::bufferLT;

    template<bliss::concurrent::LockType PoolLT, bliss::concurrent::LockType BufferLT>
    const bliss::concurrent::LockType BufferPool<PoolLT, BufferLT>::poolLT;


  } /* namespace io */
} /* namespace bliss */

#endif /* BUFFERPOOL_HPP_ */
