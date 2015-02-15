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

#include <omp.h>

// TODO: convert to an allocator.

namespace bliss
{
  namespace io
  {

    /**
     * @class     ObjectPool
     * @brief     ObjectPool manages a set of reusable, in-memory buffers.  In particular, it can limit the amount of memory usage.
     * @details   The class is templated to provide thread-safe or unsafe ObjectPools, managing thread-safe or unsafe Objects
     *
     *            When the pool is exhausted, Acquire function returns false.
     *
     *      supports mutex or spinlock for thread safety or no safety.
     *
     *
     *      object life cycle:
     *              a thread calls pool's acquire() to get pointer to an allocated object
     *              application:  uses object, potentially by multiple threads
     *              when done, 1 thread releases the object to the pool.  repeated release of the same object pointer does nothing.
     *
     *
     *      internally, there is a queue of available objects and a set of in use objects.
     *           the use of "in_use" set serves 2 purposes:  ensure there is no memory leak, and prevent multiple releases for the same object.
     *
     * @note      Each object should be acquired by a single thread.
     *            object may be used by multiple threads
     *            multiple threads may acquire concurrently.
     *            the object MAY be released by multiple threads.
     *
     *            It is important to address race conditions when multithreading, including
     *              likely loss of data (thread 1 appends data while thread 2 releases it;
     *              thread 2 and 3 release the same object, thread 1 successfully acquires before thread 3 releases object).
     *
     *            MessageBuffers class internally ensures that this does not happen.
     *
     *
     *  @tparam LockType    The thread safety property for the pool
     *  @tparam T           The object instance type.  object pool stores pointers to T instances.
     */
    template<bliss::concurrent::LockType LockType, class T >
    class ObjectPool
    {
      public:
        /// object type
        using ObjectType = T;

        /// object pointer type  (prefer shared pointer as c++11 allows atomic operations, but gcc does not support it.
        using ObjectPtrType = T*;

        /// type of lock for the pool. (mutex, spinlock,  or none).
        static const bliss::concurrent::LockType poolLT = LockType;

      protected:
        /**
         * @brief     capacity of the ObjectPool (in number of Objects)
         */
        mutable int64_t capacity;

        /**
         * @brief     Internal queue of available Objects for immediate use.
         * @details   ThreadSafeQueue to ensure thread safety.
         */
        bliss::concurrent::ThreadSafeQueue<ObjectPtrType>  available;

        /**
         * @brief     Internal Set of objects in use..
         * @details   locked for threadsafe access.
         */
        std::unordered_set<ObjectPtrType>                  in_use;

        /// num objects in active use.
        std::atomic<int64_t>  size_in_use;

        /// mutex for thread safety
        mutable std::mutex mutex;

        /// spinlock for threadsafety.
        mutable std::atomic_flag spinlock = ATOMIC_FLAG_INIT;


        /// clear all allocated objects.
        void clear_storage() {
          ObjectPtrType ptr;
          auto entry = available.tryPop();
          while (entry.first) {
             delete entry.second;
             entry = available.tryPop();
          }

          std::unique_lock<std::mutex> lock(mutex);
          for (auto ptr : in_use) {
            delete ptr;
          }
          in_use.clear();
          lock.unlock();

          size_in_use.store(0, std::memory_order_release);
        }

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
         * @brief     construct a Object Pool with buffers of capacity _buffer_capacity.
         *   the number of buffers in the pool is set to _pool_capacity.
         *
         * @param _pool_capacity       number of buffers in the pool
         * @param _buffer_capacity     size of the individual buffers
         */
        explicit ObjectPool(const int64_t _pool_capacity = std::numeric_limits<int64_t>::max()) :
          capacity(_pool_capacity),
          available(_pool_capacity), in_use(), size_in_use(0)
          {};


        /// move constructor is deleted.
        explicit ObjectPool(ObjectPool<LockType, T>&& other) = delete;
        /// move assignment operator is deleted.
        ObjectPool<LockType, T>& operator=(ObjectPool<LockType, T>&& other) = delete;


        /**
         * @brief     default destructor
         */
        virtual ~ObjectPool() {
          // delete all the objects
          capacity = 0;

          clear_storage();
        };

        /**
         * @brief     Current size of the ObjectPool.  For debugging only.
         * @return  size,
         */
        int64_t getAvailableCount() const  {
          return (isUnlimited() ? capacity : (capacity - size_in_use.load(std::memory_order_relaxed)));
        }


        /**
         * @brief     Current capacity of the ObjectPool
         * @return    capacity
         */
        const int64_t getCapacity() const {
          return capacity;
        }

        /// check if the pool has unlimited capacity
        inline const bool isUnlimited() const {
          return capacity == std::numeric_limits<int64_t>::max();
        }


        /**
         * @brief     Resets the entire ObjectPool: all Objects in pool are  marked as released and available.
         *
         * @note      This may not be entirely thread safe.  the available set is cleared by a single thread, but other
         *            threads may be acquiring objects while they are being released.
         *
         *            It is envisioned that this function should be called from a single thread.
         *            however, care must be taken to ensure that all other threads have completed their work, else data loss is likely.
         */
        void reset() {
          std::lock_guard<std::mutex> lock(mutex);

          // move all from in_use to available
          for (auto iter = in_use.begin(); iter != in_use.end(); ++iter) {
            available.tryPush(*iter);
          }

          // TODO: clear the objects somehow?
          in_use.clear();
          size_in_use.store(0, std::memory_order_release);
        }


        /**
         * @brief     Get the next available Object.  Mutex Locked.
         * @details   if none available, and size is below capacity or pool is unlimited,
         *            a new object is created.
         * @return    pointer to acquired object, nullptr if failed.
         */
        template<bliss::concurrent::LockType LT = LockType>
                typename std::enable_if<LT == bliss::concurrent::LockType::MUTEX, ObjectPtrType>::type
        tryAcquireObject() {


          ObjectPtrType sptr = nullptr;  // default is a null ptr.

          // okay to add before compare, since object pool is long lasting.
          int64_t prev_size = size_in_use.fetch_add(1, std::memory_order_relaxed);
          if (prev_size >= capacity) {
            size_in_use.fetch_sub(1, std::memory_order_relaxed);
            // leave ptr as nullptr.
            if (this->isUnlimited()) {
              ERRORF("ERROR: pool is full but should be unlimited. prev size %lu.", prev_size);
            }
          } else {

            // now get or create
            if (available.isEmpty()) {

                // none available for reuse
                // but has room to allocate, so do it.
                sptr = new T();

            } else {
              // has available for reuse.
              sptr = available.tryPop().second;
              // available may return null because it's concurrent queue.
              // cannot waitAndPop, because available could be empty at this point.
            }

            if (sptr) {

              std::unique_lock<std::mutex> lock(mutex);
	      // save in in_use set.
              in_use.insert(sptr);  // store the shared pointer.

              lock.unlock();

            }
          }

          return sptr;

        }

        /**
         * @brief     Get the next available Object.  spinlocked.
         * @details   if none available, and size is below capacity or pool is unlimited,
         *            a new object is created.
         * @return    pointer to acquired object, nullptr if failed.
         */
        template<bliss::concurrent::LockType LT = LockType>
                typename std::enable_if<LT == bliss::concurrent::LockType::SPINLOCK, ObjectPtrType>::type 
	tryAcquireObject() {

          ObjectPtrType sptr = nullptr;  // default is a null ptr.

          // okay to add before compare, since object pool is long lasting.
          int64_t prev_size = size_in_use.fetch_add(1, std::memory_order_relaxed);
          if (prev_size >= capacity) {

            size_in_use.fetch_sub(1, std::memory_order_relaxed);

            if (this->isUnlimited()) {
              ERRORF("ERROR: pool is full but should be unlimited. prev size %lu.", prev_size);
            }

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

              // available may return null because it's concurrent queue.
              // cannot wait for available to return non-null, because available could be empty at this point.
            }

            if (sptr) {
              while (spinlock.test_and_set());
		
	      // save in in-use set.
              in_use.insert(sptr);  // store the shared pointer.

              spinlock.clear();
            }
          }

          return sptr;
        }


        /**
         * @brief     Release a buffer back to pool.  Threadsafe via mutex lock
	 * @details   clears object before release back into available.
         * @param ptr weak_ptr to object.
         */
        template<bliss::concurrent::LockType LT = LockType>
        typename std::enable_if<LT == bliss::concurrent::LockType::MUTEX, bool>::type
	releaseObject(ObjectPtrType ptr) {

          if (!ptr)
          {
            WARNINGF("WARNING: pool releasing a nullptr.");
            return true;
          }

          std::unique_lock<std::mutex> lock(mutex);
          int count = in_use.erase(ptr);
          lock.unlock();

          bool res = false;      

          if (count > 0) { // only put back in available queue if it was in use.
            // now make object available.  make sure push_back is done one thread at a time.
            res = available.tryPush(ptr);
            size_in_use.fetch_sub(count, std::memory_order_release);
          } // else return false.

          return res;
        }


        /**
         * @brief     Release a buffer back to pool.  Threadsafe via spinlock
	 * @details   clears object before release back into available.
         * @param ptr weak_ptr to object.
         */
        template<bliss::concurrent::LockType LT = LockType>
        typename std::enable_if<LT == bliss::concurrent::LockType::SPINLOCK, bool>::type 
	releaseObject(ObjectPtrType ptr) {

          if (!ptr)
          {
            WARNINGF("WARNING: pool releasing a nullptr.");
            return true;
          }

          while (spinlock.test_and_set());
          int count = in_use.erase(ptr);
          spinlock.clear();

          bool res = false;           

          if (count > 0) { // only put object back in available queue if it was in use.
            // now make object available.  make sure push_back is done one thread at a time.
            res = available.tryPush(ptr);
            size_in_use.fetch_sub(count, std::memory_order_release);

          } // else return false.

          return res;

        }
    };

    /**
     * 
     */
    template<class T >
    class ObjectPool<bliss::concurrent::LockType::THREADLOCAL, T>
    {
      public:
        /**
         * @brief     Index Type used to reference the Objects in the ObjectPool.
         */
        using ObjectType = T;
        using ObjectPtrType = T*;   // shared pointer allows atomic operations as well as check for expired pointers
                                    // however, GCC does not support it.

        static const bliss::concurrent::LockType poolLT = bliss::concurrent::LockType::THREADLOCAL;


      protected:


        /**
         * @brief     capacity of the ObjectPool (in number of Objects)
         */
        mutable int64_t capacity;

        /**
         * @brief     Internal Set of available Objects for immediate use.
         * @details   Using set instead of ThreadSafeQueue to ensure uniqueness of Object Ids.
         */
        std::vector< std::deque<ObjectPtrType> >                        available;
        std::vector< std::unordered_set<ObjectPtrType> >                 in_use;

        /**
         * NOT thread safe, so need to be wrapped in synchronized calls.
         */
        void clear_storage() {
          ObjectPtrType ptr;
          for (int64_t t = available.size()-1 ; t >= 0; --t) {
            while (!available[t].empty()) {
             ptr = available[t].front();
             available[t].pop_front();
             delete ptr;
            }
            available[t].clear();

            for (auto ptr : in_use[t]) {
              delete ptr;
            }
            in_use[t].clear();
          }
        }



//        /**
//         * @brief     Private move constructor with mutex lock.
//         * @param other   Source ObjectPool object to move
//         * @param l       Mutex lock on the Source ObjectPool object
//         */
//        ObjectPool(ObjectPool<lt, T>&& other, const std::lock_guard<std::mutex>&) :
//          capacity(other.capacity)
//
//        {
//          other.capacity = 0;
//
//          available.swap(other.available);
//          in_use.swap(other.in_use);
//        };


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
        explicit ObjectPool(int _num_threads, const int64_t _pool_capacity = std::numeric_limits<int64_t>::max()) :
          capacity(_num_threads, _pool_capacity), available(_num_threads, std::deque<ObjectPtrType>()),
          in_use(_num_threads, std::unordered_set<ObjectPtrType>())  {
        };


//        /**
//         * @brief      Move constructor.  Delegates to the Move constructor that locks access on the source.
//         * @param other    source ObjectPool object to move.
//         */
//        explicit ObjectPool(ObjectPool<LockType, T>&& other) :
//          ObjectPool<LockType, T>(std::move(other), std::lock_guard<std::mutex>(other.mutex)) {};
//
//        /**
//         * @brief     move assignment operator.
//         * @param other   source ObjectPool object to move from.
//         * @return        self, with member variables moved from other.
//         */
//        ObjectPool<LockType, T>& operator=(ObjectPool<LockType, T>&& other) {
//          std::unique_lock<std::mutex> mylock(mutex, std::defer_lock),
//                                        otherlock(other.mutex, std::defer_lock);
//          std::lock(mylock, otherlock);
//
//          capacity = other.capacity; other.capacity = 0;
//
//          clear_storage();
//
//          available.swap(other.available);
//          in_use.swap(other.in_use);
//
//          return *this;
//        };
        explicit ObjectPool(ObjectPool<poolLT, T>&& other) = delete;
        ObjectPool<poolLT, T>& operator=(ObjectPool<poolLT, T>&& other) = delete;

        /**
         * @brief     default destructor
         */
        virtual ~ObjectPool() {
          // delete all the objects
          clear_storage();
          capacity = 0;
        };


        /**
         * @brief     Current size of the ObjectPool
         * @return  size, type IdType (aka int).
         */
        int64_t getAvailableCount(int thread_id = -1) const  {
          int tid = thread_id < 0 ? omp_get_thread_num() : thread_id;
          return (isUnlimited() ? capacity : (capacity - in_use[tid].size()));
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
        void resetForThread(int thread_id = -1) {
          int tid = thread_id < 0 ? omp_get_thread_num() : thread_id;
          // move all from in_use to available
          available[tid].insert(available[tid].end(), in_use[tid].begin(), in_use[tid].end());
          // TODO: clear the objects somehow?
          in_use[tid].clear();
        }
        void reset() {
          for (int tid = available.size() - 1; tid >= 0; --tid) {
            // move all from in_use to available
            available[tid].insert(available[tid].end(), in_use[tid].begin(), in_use[tid].end());
            // TODO: clear the objects somehow?
            in_use[tid].clear();
          }
        }


        ObjectPtrType tryAcquireObject(int thread_id = -1) {
          int tid = thread_id < 0 ? omp_get_thread_num() : thread_id;

          ObjectPtrType sptr = nullptr;  // default is a null ptr.

          size_t size = in_use[tid].size();

            // now get or create
            if (available[tid].empty()) {
              if (size < capacity) {

                // none available for reuse
                // but has room to allocate, so do it.
                sptr = new T();
              }  else {
                // else already nullptr.
                if (this->isUnlimited()) {
                  ERRORF("ERROR: pool is full but should be unlimited. prev size %lu.", size);
//                } else {
//                  WARNINGF("WARNING: nonConcurrent pool is full. prev size %lu.", size);
                }
              }
            } else {
              // has available for reuse.
              sptr = available[tid].front();
              available[tid].pop_front();


              // available may return null because it's concurrent queue.
              // cannot wait for available to return non-null, because available could be empty at this point.
            }

            // get the object ready for use.
            if (sptr) in_use[tid].insert(sptr);  // store the shared pointer.
          return sptr;
        }

        /**
         * @brief     Release a object back to pool, by id.  if id is incorrect, throw std exception.  clears object before release back into available.
         * @note      THREAD UNSAFE
         * @param objectId    The id of the Object to be released.
         */
        bool releaseObject(ObjectPtrType ptr, int thread_id = -1) {

          if (!ptr)
          {
            WARNINGF("WARNING: pool releasing a nullptr.");
            //throw std::logic_error("ERROR pool releasing a nullptr.");
            return true;
          }

          bool res = false;            // nullptr would not be in in_use.
          int tid = thread_id < 0 ? omp_get_thread_num() : thread_id;

          if (in_use[tid].erase(ptr) > 0) { // object is in in_use.  if not in use, release could double count the entry.
            // now make object available.  make sure push_back is done one thread at a time.
            available[tid].push_back(ptr);
            res = true;
          } // else return false.

          return res;
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
        std::unordered_set<ObjectPtrType>                  in_use;

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

          for (auto ptr : in_use) {
            delete ptr;
          }
          in_use.clear();
        }



//        /**
//         * @brief     Private move constructor with mutex lock.
//         * @param other   Source ObjectPool object to move
//         * @param l       Mutex lock on the Source ObjectPool object
//         */
//        ObjectPool(ObjectPool<lt, T>&& other, const std::lock_guard<std::mutex>&) :
//          capacity(other.capacity)
//
//        {
//          other.capacity = 0;
//
//          available.swap(other.available);
//          in_use.swap(other.in_use);
//        };


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
          available(), in_use()
          {};


//        /**
//         * @brief      Move constructor.  Delegates to the Move constructor that locks access on the source.
//         * @param other    source ObjectPool object to move.
//         */
//        explicit ObjectPool(ObjectPool<LockType, T>&& other) :
//          ObjectPool<LockType, T>(std::move(other), std::lock_guard<std::mutex>(other.mutex)) {};
//
//        /**
//         * @brief     move assignment operator.
//         * @param other   source ObjectPool object to move from.
//         * @return        self, with member variables moved from other.
//         */
//        ObjectPool<LockType, T>& operator=(ObjectPool<LockType, T>&& other) {
//          std::unique_lock<std::mutex> mylock(mutex, std::defer_lock),
//                                        otherlock(other.mutex, std::defer_lock);
//          std::lock(mylock, otherlock);
//
//          capacity = other.capacity; other.capacity = 0;
//
//          clear_storage();
//
//          available.swap(other.available);
//          in_use.swap(other.in_use);
//
//          return *this;
//        };
        explicit ObjectPool(ObjectPool<poolLT, T>&& other) = delete;
        ObjectPool<poolLT, T>& operator=(ObjectPool<poolLT, T>&& other) = delete;

        /**
         * @brief     default destructor
         */
        virtual ~ObjectPool() {
          // delete all the objects
          capacity = 0;

          clear_storage();
        };


        /**
         * @brief     Current size of the ObjectPool
         * @return  size, type IdType (aka int).
         */
        int64_t getAvailableCount() const  {
          return (isUnlimited() ? capacity : (capacity - in_use.size()));
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

          // move all from in_use to available
          available.insert(available.end(), in_use.begin(), in_use.end());
          // TODO: clear the objects somehow?
          in_use.clear();

        }


        ObjectPtrType tryAcquireObject() {


          ObjectPtrType sptr = nullptr;  // default is a null ptr.

          size_t size = in_use.size();

            // now get or create
            if (available.empty()) {
              if (size < capacity) {

                // none available for reuse
                // but has room to allocate, so do it.
                sptr = new T();
              }  else {
                // else already nullptr.
                if (this->isUnlimited()) {
                  ERRORF("ERROR: pool is full but should be unlimited. prev size %lu.", size);
//                } else {
//                  WARNINGF("WARNING: nonConcurrent pool is full. prev size %lu.", size);
                }
              }
            } else {
              // has available for reuse.
              sptr = available.front();
              available.pop_front();


              // available may return null because it's concurrent queue.
              // cannot wait for available to return non-null, because available could be empty at this point.
            }

            // get the object ready for use.
            if (sptr) in_use.insert(sptr);  // store the shared pointer.
          return sptr;
        }

        /**
         * @brief     Release a object back to pool, by id.  if id is incorrect, throw std exception.  clears object before release back into available.
         * @note      THREAD UNSAFE
         * @param objectId    The id of the Object to be released.
         */
        bool releaseObject(ObjectPtrType ptr) {

          if (!ptr)
          {
            WARNINGF("WARNING: pool releasing a nullptr.");
            //throw std::logic_error("ERROR pool releasing a nullptr.");
            return true;
          }

          bool res = false;            // nullptr would not be in in_use.

          if (in_use.erase(ptr) > 0) { // object is in in_use.  if not in use, release could double count the entry.
            // now make object available.  make sure push_back is done one thread at a time.
            available.push_back(ptr);
            res = true;
          } // else return false.

          return res;
        }
    };

    template<bliss::concurrent::LockType LockType, class T>
    const bliss::concurrent::LockType ObjectPool<LockType, T>::poolLT;


  } /* namespace io */
} /* namespace bliss */

#endif /* OBJECTPOOL_HPP_ */
