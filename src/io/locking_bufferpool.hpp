/**
 * @file		object_pool.hpp
 * @ingroup bliss::io
 * @author	tpan
 * @brief   this file defines in-memory Object Pool.
 * @details
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef OBJECTPOOL_HPP_
#define OBJECTPOOL_HPP_

#include <cassert>

#include <mutex>
#include <deque>
#include <unordered_set>
#include <stdexcept>

#include "utils/logging.h"
#include "io/locking_buffer.hpp"
#include "concurrent/concurrent.hpp"

// TODO: change to an object pool, not just a buffer pool.

namespace bliss
{
  namespace io
  {
    // TODO: move constructor and assignment operators between ObjectPools of different thread safeties.
    // TODO: use threadsafe_queue where appropriate.

    /**
     * @class     ObjectPool
     * @brief     ObjectPool manages a set of reusable, in-memory buffers.  In particular, it can limit the amount of memory usage.
     * @details   The class is templated to provide thread-safe or unsafe ObjectPools, managing thread-safe or unsafe Objects
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
     *            Note that Object does not actualy clear the memory, just resets the size, which is also the pointer for insertion.
     *            While the release method handles duplicate releases, there may be race conditions when multithreading, including
     *              likely loss of data (thread 1 appends data and thread 2 releases  it.
     *
     *            the pools passes ownership of the buffers to callers.  This is so to reduce race conditions when multiple
     *            threads are using and acquiring buffers.
     *
     *
     *
     *
     *  @tparam LockType    The thread safety property for the pool
     *  @tparam ObjectLT  The thread safety property of each Object.
     */
    template<typename ObjectType, typename ObjectAllocatorType, bliss::concurrent::LockType LockType>
    class ObjectPool
    {
      public:
        /**
         * @brief     Index Type used to reference the Objects in the ObjectPool.
         */
        using ObjectPtrType = ObjectType*;   // shared pointer allows atomic operations as well as check for expired pointers.

        static const bliss::concurrent::LockType poolLT = LockType;


      protected:
        /**
         * @brief     capacity of the ObjectPool (in number of Objects)
         */
        mutable int64_t capacity;

        /**
         * @brief     current pool size
         */
        mutable typename std::conditional<LockType == bliss::concurrent::LockType::LOCKFREE, std::atomic<int64_t>, int64_t>::type numAvailable;

        /**
         * @brief     Internal Set of available Objects for immediate use.
         * @details   Using set instead of ThreadSafeQueue to ensure uniqueness of Object Ids.
         */
        std::deque<ObjectPtrType>                          available;
        std::unordered_set<ObjectPtrType>                  in_use;

        ObjectAllocatorType alloc;

        /**
         * @brief     mutex to control access.
         */
        mutable std::mutex mutex;
        mutable std::atomic_flag spinlock = ATOMIC_FLAG_INIT;


        /**
         * @brief     Private move constructor with mutex lock.
         * @param other   Source ObjectPool object to move
         * @param l       Mutex lock on the Source ObjectPool object
         */
        ObjectPool(ObjectPool<ObjectType, ObjectAllocatorType, LockType>&& other, const std::lock_guard<std::mutex>&) :
          capacity(other.capacity),
          numAvailable((int64_t)(other.numAvailable)),
          alloc(std::move(other.alloc))

        {
          other.capacity = 0;
          other.numAvailable = 0;

          ObjectPtrType ptr;
          while (!available.empty()) {
            ptr = available.front();
            delete ptr;
            available.pop_front();
          }
          available.clear();  available.swap(other.available);

          for (auto ptr : in_use) {
            delete ptr;
          }
          in_use.clear();  in_use.swap(other.in_use);
        };


      public:
        /**
         * @brief     construct a Object Pool with buffers of capacity _buffer_capacity.  the number of buffers in the pool is set to _pool_capacity.
         *
         * @param _pool_capacity       number of buffers in the pool
         * @param _buffer_capacity     size of the individual buffers
         */
        ObjectPool(const int64_t _pool_capacity, ObjectAllocatorType & _alloc) :
          capacity(_pool_capacity),
          numAvailable(_pool_capacity),
          available(), in_use(), alloc(_alloc)
          {};

        /**
         * @brief     default constructor is deleted.
         */
        ObjectPool() = delete;

        /**
         * @brief     construct a Object pool with buffers of capacity _buffer_capacity.  the number of buffers in the pool is unbounded.
         *
         * @param _buffer_capacity    size of the individual buffers
         */
        explicit ObjectPool(const size_t _buffer_capacity) :
            ObjectPool<ObjectType, ObjectAllocatorType, LockType>(std::numeric_limits<int64_t>::max(), _buffer_capacity) {};

        /**
         * @brief     default copy constructor is deleted.
         * @param other   source ObjectPool object to copy from.
         */
        explicit ObjectPool(const ObjectPool<ObjectType, ObjectAllocatorType, LockType>& other) = delete;

        /**
         * @brief      Move constructor.  Delegates to the Move constructor that locks access on the source.
         * @param other    source ObjectPool object to move.
         */
        explicit ObjectPool(ObjectPool<ObjectType, ObjectAllocatorType, LockType>&& other) :
          ObjectPool<ObjectType, ObjectAllocatorType, LockType>(std::move(other), std::lock_guard<std::mutex>(other.mutex)) {};

        /**
         * @brief     default copy assignment operator is deleted.
         * @param other   source ObjectPool object to copy from.
         * @return        self, with member variables copied from other.
         */
        ObjectPool<ObjectType, ObjectAllocatorType, LockType>& operator=(const ObjectPool<ObjectType, ObjectAllocatorType, LockType>& other) = delete;

        /**
         * @brief     move assignment operator.
         * @param other   source ObjectPool object to move from.
         * @return        self, with member variables moved from other.
         */
        ObjectPool<ObjectType, ObjectAllocatorType, LockType>& operator=(ObjectPool<ObjectType, ObjectAllocatorType, LockType>&& other) {
          std::unique_lock<std::mutex> mylock(mutex, std::defer_lock),
                                        otherlock(other.mutex, std::defer_lock);
          std::lock(mylock, otherlock);

          capacity = other.capacity; other.capacity = 0;
          numAvailable = (int64_t)(other.numAvailable); other.numAvailable = 0;

          ObjectPtrType ptr;
          while (!available.empty()) {
            ptr = available.front();
            delete ptr;
            available.pop_front();
          }
          available.clear();  available.swap(other.available);

          for (auto ptr : in_use) {
            delete ptr;
          }
          in_use.clear();  in_use.swap(other.in_use);
          alloc = std::move(other.alloc);


          return *this;
        };


        /**
         * @brief     default destructor
         */
        virtual ~ObjectPool() {
          // delete all the objects
          std::lock_guard<std::mutex> lock(mutex);
          numAvailable = 0;
          capacity = 0;
          ObjectPtrType ptr;
          while (!available.empty()) {
            ptr = available.front();
            delete ptr;
            available.pop_front();
          }
          available.clear();

          for (auto ptr : in_use) {
            delete ptr;
          }
          in_use.clear();
        };

        /**
         * @brief     Current size of the ObjectPool
         * @return  size, type IdType (aka int).
         */
        const int64_t getAvailableCount() const  {
          return (int64_t)(numAvailable);   // load is atomic
        }

        /**
         * @brief     Current capacity of the ObjectPool
         * @return    capacity, type IdTyype (aka int)
         */
        const int64_t getCapacity() const {
          return capacity;
        }

        inline const bool isUnlimited() const {
          return capacity == std::numeric_limits<int64_t>::max();
        }


        /**
         * @brief     Resets the entire ObjectPool: all Objects in pool are  marked as released and available.
         *
         * @note      This is not entirely thread safe.  the available set is cleared by a single thread, but other
         *            threads may be acquiring objects while they are being released.
         *
         *            It is envisioned that this function should be called from a single thread.
         *            however, care must be taken to ensure that all other threads have completed their work, else data loss is likely.
         */
        void reset() {
          numAvailable = capacity;  // store is atomic .
        }


        /**
         * @brief     Get the next available Object by id.  if none available, return false in the first argument of the pair.
         * @return    std::pair with bool and Id,  first indicate whether acquire was successful, second indicate the ObjectId if successful.
         */
        template<bliss::concurrent::LockType LT = LockType>
                typename std::enable_if<LT == bliss::concurrent::LockType::MUTEX, ObjectPtrType>::type tryAcquireObject() {
          ObjectPtrType sptr = nullptr;  // default is a null ptr.

          std::unique_lock<std::mutex> lock(mutex);

          if (!isUnlimited()) {
            if (numAvailable <= 0) {  // not fix size, and no more available
              lock.unlock();
              return sptr;
            } else {
              numAvailable--;
            }
          }

          // now get or create
          if (available.empty()) {
            // none available for reuse
            // but has room to allocate, so do it.
            sptr = alloc();
          } else {
            // has available for reuse.
            sptr = available.front();
            available.pop_front();
          }

          // get the object ready for use.
          in_use.insert(sptr);  // store the shared pointer.

          lock.unlock();

          return std::move(sptr);
        }

        /**
         * @brief     Get the next available Object by id.  if none available, return false in the first argument of the pair.
         * @return    std::pair with bool and Id,  first indicate whether acquire was successful, second indicate the ObjectId if successful.
         */
        template<bliss::concurrent::LockType LT = LockType>
                typename std::enable_if<LT == bliss::concurrent::LockType::SPINLOCK, ObjectPtrType>::type tryAcquireObject() {
          ObjectPtrType sptr = nullptr;  // default is a null ptr.

          while (spinlock.test_and_set());

          if (!isUnlimited()) {
            if (numAvailable <= 0) {  // not fix size, and no more available
              spinlock.clear();
              return sptr;
            } else {
              numAvailable--;
            }
          }

          // now get or create
          if (available.empty()) {
            // none available for reuse
            // but has room to allocate, so do it.
            sptr = alloc();
          } else {
            // has available for reuse.
            sptr = available.front();
            available.pop_front();
          }


          // get the object ready for use.
          in_use.insert(sptr);  // store the shared pointer.
          spinlock.clear();

          return std::move(sptr);
        }
//
//        template<bliss::concurrent::LockType LT = LockType>
//                typename std::enable_if<LT == bliss::concurrent::LockType::LOCKFREE, ObjectPtrRefType>::type tryAcquireObject() {
//          static_assert(false, "ObjectPool currently does not support LockFree operations due to internal use of deque");
//        }
//


        template<bliss::concurrent::LockType LT = LockType>
               typename std::enable_if<LT == bliss::concurrent::LockType::NONE, ObjectPtrType>::type tryAcquireObject() {
          ObjectPtrType sptr = nullptr;  // default is a null ptr.

          if (!isUnlimited()) {
            if (numAvailable <= 0) {  // not fix size, and no more available
              return sptr;
            } else {
              numAvailable--;
            }
          }

          // now get or create
          if (available.empty()) {
            // none available for reuse
            // but has room to allocate, so do it.
            sptr = alloc();
          } else {
            // has available for reuse.
            sptr = available.front();
            available.pop_front();
          }

          // get the object ready for use.

          in_use.insert(sptr);  // store the shared pointer.

          return std::move(sptr);
        }

        /**
         * @brief     Release a buffer back to pool, by id.  if id is incorrect, throw std exception.  clears object before release back into available.
         * @note      THREAD SAFE
         * @param ptr weak_ptr to object.
         */
        template<bliss::concurrent::LockType LT = LockType>
        typename std::enable_if<LT == bliss::concurrent::LockType::MUTEX, bool>::type releaseObject(ObjectPtrType ptr) {

          if (! ptr) return false;

          std::unique_lock<std::mutex> lock(mutex);
          if (!isUnlimited()) {  // check against size_t max - thread safe, and ensures there is no overflow.

            // if multiple threads, the preincrement is atomic, so each thread has a valid value to compare to capacity.
            // the total number of threads to continue will be  same as total number of values that pass the capacity test.
            if (numAvailable >= capacity) {  // preincrement is a.fetch_add(1)+1.  compare result to capacity.
              // not unlimited, and already >= capacity.  objectPtr is now lost, but that's okay.
              // undo and return.
              lock.unlock();
              return false;
            } else {
              // else not at capacity and already incremented it.
              numAvailable++;
            }
          } // else don't need to touch numAvailable.

          // now store the object.  make sure push_back is done one thread at a time.
          in_use.erase(ptr);
          available.push_back(ptr);
          lock.unlock();

          return true;

        }


        /**
         * @brief     Release a object back to pool, by id.  if id is incorrect, throw std exception.  clears object before release back into available.
         * @note      THREAD SAFE
         * @param objectId    The id of the Object to be released.
         */
        template<bliss::concurrent::LockType LT = LockType>
        typename std::enable_if<LT == bliss::concurrent::LockType::SPINLOCK, bool>::type releaseObject(ObjectPtrType ptr) {
          if (! ptr) return false;

          while (spinlock.test_and_set());
          if (!isUnlimited()) {  // check against size_t max - thread safe, and ensures there is no overflow.

            // if multiple threads, the preincrement is atomic, so each thread has a valid value to compare to capacity.
            // the total number of threads to continue will be  same as total number of values that pass the capacity test.
            if (numAvailable >= capacity) {  // preincrement is a.fetch_add(1)+1.  compare result to capacity.
              // not unlimited, and already >= capacity.  objectPtr is now lost, but that's okay.
              // undo and return.
              spinlock.clear();
              return false;
            } else {
              // else not at capacity and already incremented it.
              numAvailable++;
            }
          } // else don't need to touch numAvailable.

          // now store the object.  make sure push_back is done one thread at a time.
          in_use.erase(ptr);
          available.push_back(ptr);
          spinlock.clear();

          return true;

        }
//
//        /**
//         * @brief     Release a object back to pool, by id.  if id is incorrect, throw std exception.  clears object before release back into available.
//         * @note      THREAD SAFE
//         * @param objectId    The id of the Object to be released.
//         */
//        template<bliss::concurrent::LockType LT = LockType>
//        typename std::enable_if<LT == bliss::concurrent::LockType::LOCKFREE, bool>::type releaseObject(ObjectPtrRefType&& ptr) {
//          static_assert(false, "ObjectPool currently does not support LockFree operations due to internal use of deque");
//        }


        /**
         * @brief     Release a object back to pool, by id.  if id is incorrect, throw std exception.  clears object before release back into available.
         * @note      THREAD UNSAFE
         * @param objectId    The id of the Object to be released.
         */
        template<bliss::concurrent::LockType LT = LockType>
        typename std::enable_if<LT == bliss::concurrent::LockType::NONE, bool>::type releaseObject(ObjectPtrType ptr) {
           if (! ptr) return false;

          if (!isUnlimited()) {  // check against size_t max - thread safe, and ensures there is no overflow.

            // if multiple threads, the preincrement is atomic, so each thread has a valid value to compare to capacity.
            // the total number of threads to continue will be  same as total number of values that pass the capacity test.
            if (numAvailable >= capacity) {  // preincrement is a.fetch_add(1)+1.  compare result to capacity.
              // not unlimited, and already >= capacity.  objectPtr is now lost, but that's okay.
              // undo and return.
              return false;
            } else {
              // else not at capacity and already incremented it.
              numAvailable++;
            }
          } // else don't need to touch numAvailable.

          // now store the object.  make sure push_back is done one thread at a time.
          in_use.erase(ptr);
          available.push_back(ptr);

          return true;
        }

    };

    template<typename ObjectType, typename ObjectAllocatorType, bliss::concurrent::LockType LockType>
    const bliss::concurrent::LockType ObjectPool<ObjectType, ObjectAllocatorType, LockType>::poolLT;


  } /* namespace io */
} /* namespace bliss */

#endif /* OBJECTPOOL_HPP_ */
