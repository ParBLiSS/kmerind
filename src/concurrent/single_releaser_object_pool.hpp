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

#include <atomic>
#include <mutex>
#include <deque>
#include <unordered_set>
#include <stdexcept>

#include "utils/logging.h"
#include "concurrent/concurrent.hpp"
#include "concurrent/lockfree_queue.hpp"

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
     *      life cycle:  pool acquire() by a single thread.
     *              application:  uses object, potentially by multiple threads
     *              pool release():  returns object to pool.  this class supports releasing an object by multiple threads
     *                only first one succeeds.
     *
     *
     *           note that the use of "in_use" set serves 2 purposes:  ensure there is no memory leak, and prevent multiple releases for the same object.
     *
     * @note      Each object should be acquired by a single thread.
     *            object may be used by multiple threads
     *            the object, when using this class, MAY be released by multiple threads.
     *              this is NOT the case for CommunicationLayer and SendMessageBuffer classes
     *                these 2 releases an object via a single thread
     *                CommLayer will also guarantee that all in-use buffers are released back to buffer.
     *
     *            It is important to address race conditions when multithreading, including
     *              likely loss of data (thread 1 appends data while thread 2 releases it;
     *              thread 2 and 3 release the same object, thread 1 successfully acquires before thread 3 releases object).
     *
     *  @note:  for now, rely on default constructor of T.  in the future, custom allocator may be better
     *
     *  @tparam LockType    The thread safety property for the pool
     *  @tparam T           The object instance type.  object pool stores pointers to T instances.
     */
    template<bliss::concurrent::LockType LockType, class T >
    class ObjectPool
    {
      public:
        /**
         * @brief     Index Type used to reference the Objects in the ObjectPool.
         */
        using ObjectType = T;
        using ObjectPtrType = T*;   // shared pointer allows atomic operations as well as check for expired pointers
                                    // however, GCC does not support it.

        static const bliss::concurrent::LockType poolLT = LockType;


      protected:
        /**
         * @brief     capacity of the ObjectPool (in number of Objects)
         */
        mutable int64_t capacity;

        /**
         * @brief     Internal Set of available Objects for immediate use.
         * @details   Using set instead of ThreadSafeQueue to ensure uniqueness of Object Ids.
         */
        bliss::concurrent::ThreadSafeQueue<ObjectPtrType>  available;

        std::atomic<int64_t>  size_in_use;

        /**
         * @brief     mutex to control access.
         */
        mutable std::mutex mutex;

        /**
         * NOT thread safe, so need to be wrapped in synchronized calls.
         */
        void clear_storage() {
          ObjectPtrType ptr;
          auto entry = available.tryPop();
          while (entry.first) {
             if (entry.second) delete entry.second;
             entry = available.tryPop();
          }
        }



        /**
         * @brief     Private move constructor with mutex lock.
         * @param other   Source ObjectPool object to move
         * @param l       Mutex lock on the Source ObjectPool object
         */
        ObjectPool(ObjectPool<LockType, T>&& other, const std::lock_guard<std::mutex>&) :
          capacity(other.capacity)

        {
          other.capacity = 0;
          size_in_use.store(other.size_in_use.exchange(0, std::memory_order_relaxed), std::memory_order_relaxed);

          available.clear();
          auto entry = other.available.tryPop();
          while (entry.first) {
             available.tryPush(entry.second);
          }
        };


        /**
         * @brief     default copy constructor is deleted.
         * @param other   source ObjectPool object to copy from.
         */
        explicit ObjectPool(const ObjectPool<LockType, T>& other) = delete;

        /**
         * @brief     default copy assignment operator is deleted.
         * @param other   source ObjectPool object to copy from.
         * @return        self, with member variables copied from other.
         */
        ObjectPool<LockType, T>& operator=(const ObjectPool<LockType, T>& other) = delete;


      public:
        /**
         * @brief     construct a Object Pool with buffers of capacity _buffer_capacity.  the number of buffers in the pool is set to _pool_capacity.
         *
         * @param _pool_capacity       number of buffers in the pool
         * @param _buffer_capacity     size of the individual buffers
         */
        explicit ObjectPool(const int64_t _pool_capacity = std::numeric_limits<int64_t>::max()) :
          capacity(_pool_capacity),
          available(_pool_capacity), size_in_use(0)
          {};


        /**
         * @brief      Move constructor.  Delegates to the Move constructor that locks access on the source.
         * @param other    source ObjectPool object to move.
         */
        explicit ObjectPool(ObjectPool<LockType, T>&& other) :
          ObjectPool<LockType, T>(std::move(other), std::lock_guard<std::mutex>(other.mutex)) {};

        /**
         * @brief     move assignment operator.
         * @param other   source ObjectPool object to move from.
         * @return        self, with member variables moved from other.
         */
        ObjectPool<LockType, T>& operator=(ObjectPool<LockType, T>&& other) {
          std::unique_lock<std::mutex> mylock(mutex, std::defer_lock),
                                        otherlock(other.mutex, std::defer_lock);
          std::lock(mylock, otherlock);

          capacity = other.capacity; other.capacity = 0;
          size_in_use.store(other.size_in_use.exchange(0, std::memory_order_relaxed), std::memory_order_relaxed);

          clear_storage();

          auto entry = other.available.tryPop();
          while (entry.first) {
             available.tryPush(entry.second);
          }

          return *this;
        };
//        explicit ObjectPool(ObjectPool<LockType, T>&& other) = delete;
//        ObjectPool<LockType, T>& operator=(ObjectPool<LockType, T>&& other) = delete;



        /**
         * @brief     default destructor
         */
        virtual ~ObjectPool() {
          // delete all the objects
          capacity = 0;
          size_in_use.store(0, std::memory_order_relaxed);

          clear_storage();
        };

        /**
         * @brief     Current size of the ObjectPool.  For debugging only.
         * @return  size, type IdType (aka int).
         */
        int64_t getAvailableCount() const  {

          return (isUnlimited() ? capacity : (capacity - size_in_use.load(std::memory_order_relaxed)));
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
          size_in_use.store(0, std::memory_order_release);
        }


        /**
         * @brief     Get the next available Object by id.  if none available, return false in the first argument of the pair.
         * @return    std::pair with bool and Id,  first indicate whether acquire was successful, second indicate the ObjectId if successful.
         */
        ObjectPtrType tryAcquireObject() {


          ObjectPtrType sptr = nullptr;  // default is a null ptr.

          int64_t prev_size = size_in_use.fetch_add(1, std::memory_order_relaxed);
          if (prev_size >= capacity) {
            size_in_use.fetch_sub(1, std::memory_order_relaxed);
            // leave ptr as nullptr.
          } else {

            // now get or create
            if (available.isEmpty()) {

                // none available for reuse
                // but has room to allocate, so do it.
                sptr = new T();

            } else {
              // has available for reuse.
              sptr = available.tryPop().second;
            }
          }

          return sptr;

        }




        /**
         * @brief     Release a object back to pool, by id.  if id is incorrect, throw std exception.  clears object before release back into available.
         * @note      THREAD SAFE
         * @param objectId    The id of the Object to be released.
         */
        bool releaseObject(ObjectPtrType ptr) {
          if (!ptr) return false;

          // now make object available.  make sure push_back is done one thread at a time.
          int64_t prev = size_in_use.fetch_sub(1, std::memory_order_relaxed);
          if (prev <= 0) {
            //delete ptr;  // this is a ptr that is beyond capacity of the internal queue, delete it and don't decrement.
            size_in_use.fetch_add(1, std::memory_order_release);
            return false;
          } else {
            return available.tryPush(ptr);    // push it onto the queue.
          }

        }
    };



    template<class T >
    class ObjectPool<bliss::concurrent::LockType::NONE, T>
    {
      public:
        /**
         * @brief     Index Type used to reference the Objects in the ObjectPool.
         */
        using ObjectType = T;
        using ObjectPtrType = T*;   // shared pointer allows atomic operations as well as check for expired pointers
                                    // however, GCC does not support it.

        static const bliss::concurrent::LockType poolLT = bliss::concurrent::LockType::NONE;


      protected:


        /**
         * @brief     capacity of the ObjectPool (in number of Objects)
         */
        mutable int64_t capacity;

        /**
         * @brief     Internal Set of available Objects for immediate use.
         * @details   Using set instead of ThreadSafeQueue to ensure uniqueness of Object Ids.
         */
        std::deque<ObjectPtrType>                          available;

        int64_t size_in_use;

        mutable std::mutex mutex;

        /**
         * NOT thread safe, so need to be wrapped in synchronized calls.
         */
        void clear_storage() {
          ObjectPtrType ptr;
          while (!available.empty()) {
           ptr = available.front();
           available.pop_front();
           delete ptr;
          }
          available.clear();

        }



        /**
         * @brief     Private move constructor with mutex lock.
         * @param other   Source ObjectPool object to move
         * @param l       Mutex lock on the Source ObjectPool object
         */
        ObjectPool(ObjectPool<poolLT, T>&& other, const std::lock_guard<std::mutex>&) :
          capacity(other.capacity), size_in_use(other.size_in_use)

        {
          other.capacity = 0;
          other.size_in_use = 0;

          available.swap(other.available);
        };


        /**
         * @brief     default copy constructor is deleted.
         * @param other   source ObjectPool object to copy from.
         */
        explicit ObjectPool(const ObjectPool<poolLT, T>& other) = delete;

        /**
         * @brief     default copy assignment operator is deleted.
         * @param other   source ObjectPool object to copy from.
         * @return        self, with member variables copied from other.
         */
        ObjectPool<poolLT, T>& operator=(const ObjectPool<poolLT, T>& other) = delete;


      public:
        /**
         * @brief     construct a Object Pool with buffers of capacity _buffer_capacity.  the number of buffers in the pool is set to _pool_capacity.
         *
         * @param _pool_capacity       number of buffers in the pool
         * @param _buffer_capacity     size of the individual buffers
         */
        explicit ObjectPool(const int64_t _pool_capacity = std::numeric_limits<int64_t>::max()) :
          capacity(_pool_capacity),
          available(), size_in_use(0)
          {};


        /**
         * @brief      Move constructor.  Delegates to the Move constructor that locks access on the source.
         * @param other    source ObjectPool object to move.
         */
        explicit ObjectPool(ObjectPool<poolLT, T>&& other) :
          ObjectPool<poolLT, T>(std::move(other), std::lock_guard<std::mutex>(other.mutex)) {};

        /**
         * @brief     move assignment operator.
         * @param other   source ObjectPool object to move from.
         * @return        self, with member variables moved from other.
         */
        ObjectPool<poolLT, T>& operator=(ObjectPool<poolLT, T>&& other) {
          std::unique_lock<std::mutex> mylock(mutex, std::defer_lock),
                                        otherlock(other.mutex, std::defer_lock);
          std::lock(mylock, otherlock);

          capacity = other.capacity; other.capacity = 0;
          size_in_use = other.size_in_use; other.size_in_use = 0;
          clear_storage();

          available.swap(other.available);

          return *this;
        };
//        explicit ObjectPool(ObjectPool<poolLT, T>&& other) = delete;
//        ObjectPool<poolLT, T>& operator=(ObjectPool<poolLT, T>&& other) = delete;

        /**
         * @brief     default destructor
         */
        virtual ~ObjectPool() {
          // delete all the objects
          capacity = 0;
          size_in_use = 0;

          clear_storage();
        };


        /**
         * @brief     Current size of the ObjectPool
         * @return  size, type IdType (aka int).
         */
        int64_t getAvailableCount() const  {
          return (isUnlimited() ? capacity : (capacity - size_in_use));
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

          size_in_use = 0;


        }


        /**
         * @brief     Get the next available Object by id.  if none available, return false in the first argument of the pair.
         * @return    std::pair with bool and Id,  first indicate whether acquire was successful, second indicate the ObjectId if successful.
         */
        ObjectPtrType tryAcquireObject() {


          ObjectPtrType sptr = nullptr;  // default is a null ptr.

          if (size_in_use >= capacity) {
            // leave ptr as nullptr.
          } else {
            ++size_in_use;
            // now get or create
            if (available.empty()) {

                // none available for reuse
                // but has room to allocate, so do it.
                sptr = new T();

            } else {
              // has available for reuse.
              sptr = available.front();
              available.pop_front();
            }
          }

          return sptr;

        }




        /**
         * @brief     Release a object back to pool, by id.  if id is incorrect, throw std exception.  clears object before release back into available.
         * @note      THREAD SAFE
         * @param objectId    The id of the Object to be released.
         */
        bool releaseObject(ObjectPtrType ptr) {
          if (!ptr) return false;

          // now make object available.  make sure push_back is done one thread at a time.
          if (size_in_use <= 0) {
            // this is a ptr that is beyond capacity of the internal queue, don't decrement.
            return false;
          } else {
            --size_in_use;
            available.push_back(ptr);    // push it onto the queue.
            return true;
          }

        }

    };

    template<bliss::concurrent::LockType LockType, class T>
    const bliss::concurrent::LockType ObjectPool<LockType, T>::poolLT;


  } /* namespace io */
} /* namespace bliss */

#endif /* OBJECTPOOL_HPP_ */
