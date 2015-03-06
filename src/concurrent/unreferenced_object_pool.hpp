/**
 * @fiule		unreferenced_object_pool.hpp
 * @ingroup concurrent
 * @author	Tony Pan <tpan7@gatech.edu>
 * @brief   this file defines in-memory Object Pool.
 * @details
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef UNREFERENCED_OBJECTPOOL_HPP_
#define UNREFERENCED_OBJECTPOOL_HPP_

#include <cassert>

#include <atomic>
#include <mutex>
#include <deque>
#include <unordered_set>
#include <stdexcept>

#include "utils/logging.h"
#include "concurrent/concurrent.hpp"
#include "concurrent/lockfree_queue.hpp"
#include "concurrent/object_pool.hpp"

#include <omp.h>

// TODO: convert to an allocator.

namespace bliss
{
  namespace concurrent
  {


    /**
     * @class     ObjectPool
     * @brief     ObjectPool manages a set of reusable, in-memory buffers.  In particular, it can limit the amount of memory usage.
     * @details   The class is templated to provide thread-safe or unsafe ObjectPools, managing thread-safe or unsafe Objects
     *
     *            When the pool is exhausted, Acquire function returns false.
     *
     *      supports mutex for thread safety.
     *
     *     this class assumes that while there may be multiple threads acquiring objects,
     *     but there is a SINGLE thread to release the object back.
     *        with this assumption, the pool can be lockfree, since the "in-use" set is unnecessary,
     *        and the available queue is already lockfree.
     *
     *        The releasing thread essentially acts as an "in-use" set.
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
    template<class T >
    class ObjectPool<T, bliss::concurrent::LockType::MUTEX, false> :
    public ObjectPoolBase<ObjectPool<T, bliss::concurrent::LockType::MUTEX, false> >
    {

      protected:

        using PoolType = ObjectPool<T, bliss::concurrent::LockType::MUTEX, false>;
        using BaseType = ObjectPoolBase<PoolType>;

        friend BaseType;

      public:
        /// object type
        using ObjectType = typename BaseType::ObjectType;

        /// object pointer type  (prefer shared pointer as c++11 allows atomic operations, however gcc does not support it.)
        using ObjectPtrType = typename BaseType::ObjectPtrType;


      protected:

        /// mutex for thread safety
        mutable std::mutex mutex;

        /**
         * @brief     default copy constructor is deleted.
         * @param other   source ObjectPool object to copy from.
         */
        explicit ObjectPool(const PoolType& other) = delete;

        /**
         * @brief     default copy assignment operator is deleted.
         * @param other   source ObjectPool object to copy from.
         * @return        self, with member variables copied from other.
         */
        PoolType& operator=(const PoolType& other) = delete;

        /// move constructor is deleted.
        explicit ObjectPool(PoolType&& other) = delete;
        /// move assignment operator is deleted.
        PoolType& operator=(PoolType&& other) = delete;


        /// clear all allocated objects.
        void clearStorageImpl() {
          std::lock_guard<std::mutex> lock(mutex);

          this->size_in_use = 0;

          // clear the available queue.
          while (!this->available.empty()) {
            delete this->available.front();
            this->available.pop_front();
          }
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
        void resetImpl() {
          std::lock_guard<std::mutex> lock(mutex);

          this->size_in_use = 0;
        }


        /**
         * @brief     Get the next available Object.  Mutex or Spin Locked.
         * @details   if none available, and size is below capacity or pool is unlimited,
         *            a new object is created.
         * @return    pointer to acquired object, nullptr if failed.
         */
        ObjectPtrType tryAcquireObjectImpl() {

          ObjectPtrType sptr = nullptr;  // default is a null ptr.

          std::unique_lock<std::mutex> lock(mutex);

          // first check if we are exceeding capacity.
          if (this->size_in_use >= this->capacity) {
            if (this->isUnlimited()) {
              WARNINGF("ERROR: pool is full but should be unlimited. size %lu.", this->size_in_use);
            }  // else limited size, so return nullptr.
          } else {  // has room.  get one.


            // get an object
            if (this->available.empty()) {  // empty. allocate one.
              sptr = new ObjectType();
            } else {

              sptr = this->available.front();
              this->available.pop_front();
            }

            assert(sptr != nullptr);

            ++this->size_in_use;

          }
          return sptr;
        }


        /**
         * @brief     Release a buffer back to pool.  Threadsafe via mutex lock
         * @details   clears object before release back into available.
         * @param ptr weak_ptr to object.
         */
        bool releaseObjectImpl(ObjectPtrType ptr) {

          if (!ptr) {
            WARNINGF("WARNING: pool releasing a nullptr.");
            return true;
          }

          std::lock_guard<std::mutex> lock(mutex);
          if (this->size_in_use == 0) {
            // trying to release when we are maxed out.
            return false;
          }

          // only put back in available queue if it was in use.
          // now make object available.  make sure push_back is done one thread at a time.
          this->available.emplace_back(ptr);
          this->size_in_use -= 1;

          return true;
        }


      public:
        /**
         * @brief     construct a Object Pool with buffers of capacity _buffer_capacity.
         *   the number of buffers in the pool is set to _pool_capacity.
         *
         * @param _pool_capacity       number of buffers in the pool
         */
        explicit ObjectPool(const int64_t _pool_capacity = std::numeric_limits<int64_t>::max()) :
        BaseType(_pool_capacity) {};


        /**
         * @brief     default destructor
         */
        virtual ~ObjectPool() {};


    };


    /**
     * @class     ObjectPool
     * @brief     ObjectPool manages a set of reusable, in-memory buffers.  In particular, it can limit the amount of memory usage.
     * @details   The class is templated to provide thread-safe or unsafe ObjectPools, managing thread-safe or unsafe Objects
     *
     *            When the pool is exhausted, Acquire function returns false.
     *
     *      supports spinlock for thread safety.
     *
     *     this class assumes that while there may be multiple threads acquiring objects,
     *     but there is a SINGLE thread to release the object back.
     *        with this assumption, the pool can be lockfree, since the "in-use" set is unnecessary,
     *        and the available queue is already lockfree.
     *
     *        The releasing thread essentially acts as an "in-use" set.
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
    template<class T >
    class ObjectPool<T, bliss::concurrent::LockType::SPINLOCK, false> :
    public ObjectPoolBase<ObjectPool<T, bliss::concurrent::LockType::SPINLOCK, false> >
    {

      protected:

        using PoolType = ObjectPool<T, bliss::concurrent::LockType::SPINLOCK, false>;
        using BaseType = ObjectPoolBase<PoolType>;

        friend BaseType;

      public:
        /// object type
        using ObjectType = typename BaseType::ObjectType;

        /// object pointer type  (prefer shared pointer as c++11 allows atomic operations, however gcc does not support it.)
        using ObjectPtrType = typename BaseType::ObjectPtrType;


      protected:

        /// mutex for thread safety
        mutable std::atomic_flag spinlock = ATOMIC_FLAG_INIT;

        /**
         * @brief     default copy constructor is deleted.
         * @param other   source ObjectPool object to copy from.
         */
        explicit ObjectPool(const PoolType& other) = delete;

        /**
         * @brief     default copy assignment operator is deleted.
         * @param other   source ObjectPool object to copy from.
         * @return        self, with member variables copied from other.
         */
        PoolType& operator=(const PoolType& other) = delete;

        /// move constructor is deleted.
        explicit ObjectPool(PoolType&& other) = delete;
        /// move assignment operator is deleted.
        PoolType& operator=(PoolType&& other) = delete;


        /// clear all allocated objects.
        void clearStorageImpl() {
          while (spinlock.test_and_set(std::memory_order_acq_rel));

          this->size_in_use = 0;

          // clear the available queue.
          while (!this->available.empty()) {
            delete this->available.front();
            this->available.pop_front();
          }

          spinlock.clear(std::memory_order_release);
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
        void resetImpl() {
          while (spinlock.test_and_set(std::memory_order_acq_rel));

          this->size_in_use = 0;

          spinlock.clear(std::memory_order_release);
        }


        /**
         * @brief     Get the next available Object.  Mutex or Spin Locked.
         * @details   if none available, and size is below capacity or pool is unlimited,
         *            a new object is created.
         * @return    pointer to acquired object, nullptr if failed.
         */
        ObjectPtrType tryAcquireObjectImpl() {

          ObjectPtrType sptr = nullptr;  // default is a null ptr.

          while (spinlock.test_and_set(std::memory_order_acq_rel));

          // first check if we are exceeding capacity.
          if (this->size_in_use >= this->capacity) {
            if (this->isUnlimited()) {
              WARNINGF("ERROR: pool is full but should be unlimited. size %lu.", this->size_in_use);
            }  // else limited size, so return nullptr.
          } else {  // has room.  get one.


            // get an object
            if (this->available.empty()) {  // empty. allocate one.
              sptr = new ObjectType();
            } else {

              sptr = this->available.front();
              this->available.pop_front();
            }

            assert(sptr != nullptr);

            ++this->size_in_use;

          }

          spinlock.clear(std::memory_order_release);
          return sptr;
        }


        /**
         * @brief     Release a buffer back to pool.  Threadsafe via mutex lock
         * @details   clears object before release back into available.
         * @param ptr weak_ptr to object.
         */
        bool releaseObjectImpl(ObjectPtrType ptr) {

          if (!ptr) {
            WARNINGF("WARNING: pool releasing a nullptr.");
            return true;
          }

          while (spinlock.test_and_set(std::memory_order_acq_rel));
          if (this->size_in_use == 0) {
            spinlock.clear(std::memory_order_release);

            // trying to release when we are maxed out.
            return false;
          }


          // only put back in available queue if it was in use.
          // now make object available.  make sure push_back is done one thread at a time.
          this->available.emplace_back(ptr);
          this->size_in_use -= 1;

          spinlock.clear(std::memory_order_release);
          return true;
        }


      public:
        /**
         * @brief     construct a Object Pool with buffers of capacity _buffer_capacity.
         *   the number of buffers in the pool is set to _pool_capacity.
         *
         * @param _pool_capacity       number of buffers in the pool
         */
        explicit ObjectPool(const int64_t _pool_capacity = std::numeric_limits<int64_t>::max()) :
        BaseType(_pool_capacity) {};


        /**
         * @brief     default destructor
         */
        virtual ~ObjectPool() {};


    };


    /**
     * @class     ObjectPool
     * @brief     ObjectPool manages a set of reusable, in-memory buffers.  In particular, it can limit the amount of memory usage.
     * @details   The class is templated to provide thread-safe ObjectPools, managing thread-safe or unsafe Objects
     *
     *            When the pool is exhausted, Acquire function returns false.
     *
     *      This template specialization is meant to support the use case where each thread has its own
     *      	in-use set, separate from other threads, yet the threads still share
     *      	the same available queue, which allows multiple threads to reuse the same pool of objects.
     *
     *      Note that this means we cannot acquire in one thread and release in another - this would cause a memory leak.
     *      This is generally lockfree, especially given the use of LockFree ThreadSafeQueue
     *
     *     this class assumes that while there may be multiple threads acquiring objects,
     *     but there is a SINGLE thread to release the object back.
     *        with this assumption, the pool can be lockfree, since the "in-use" set is unnecessary,
     *        and the available queue is already lockfree.
     *
     *        The releasing thread essentially acts as an "in-use" set.
     *
     *      object life cycle:
     *              a thread calls pool's acquire() to get pointer to an allocated object
     *              application:  uses object, potentially by multiple threads
     *              when done, 1 thread releases the object to the pool.  repeated release of the same object pointer does nothing.
     *
     *
     *      internally, there is a queue of available objects and a vector of sets of in use objects.
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
     *  @tparam LockType    The thread safety property for the pool.  This class defaults to TRHEAD Local storage
     *  @tparam T           The object instance type.  object pool stores pointers to T instances.
     */
    template<class T >
    class ObjectPool<T, bliss::concurrent::LockType::LOCKFREE, false> :
    public ObjectPoolBase<ObjectPool<T, bliss::concurrent::LockType::LOCKFREE, false> >
    {

      protected:

        using PoolType = ObjectPool<T, bliss::concurrent::LockType::LOCKFREE, false>;
        using BaseType = ObjectPoolBase<PoolType>;

        friend BaseType;

      public:
        /// object type
        using ObjectType = typename BaseType::ObjectType;

        /// object pointer type  (prefer shared pointer as c++11 allows atomic operations, however gcc does not support it.)
        using ObjectPtrType = typename BaseType::ObjectPtrType;


      protected:

        /**
         * @brief     default copy constructor is deleted.
         * @param other   source ObjectPool object to copy from.
         */
        explicit ObjectPool(const PoolType& other) = delete;

        /**
         * @brief     default copy assignment operator is deleted.
         * @param other   source ObjectPool object to copy from.
         * @return        self, with member variables copied from other.
         */
        PoolType& operator=(const PoolType& other) = delete;

        /// move constructor is deleted.
        explicit ObjectPool(PoolType&& other) = delete;
        /// move assignment operator is deleted.
        PoolType& operator=(PoolType&& other) = delete;


        /**
         * @brief clear all allocated objects.  Single Threaded.
         *
         */
        void clearStorageImpl() {
          // access is not locked on insertion side, so no point to lock here.

          std::atomic_thread_fence(std::memory_order_acq_rel);

          this->size_in_use.store(0, std::memory_order_relaxed);

          auto entry = this->available.tryPop();
          while (entry.first) {
            if (entry.second) delete entry.second;
            else WARNINGF("object pool contains nullptr in available queue!");
            entry = this->available.tryPop();
          }

          std::atomic_thread_fence(std::memory_order_release);
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
        void resetImpl() {
          // access is not locked on insertion side, so no point to lock here.
          this->size_in_use.store(0, std::memory_order_release);
        }

        /**
         * @brief     Get the next available Object.  Thread Local
         * @details   if none available, and size is below capacity or pool is unlimited,
         *            a new object is created.
         * @param thread_id   id of thread acquiring.  -1 means the current omp thread.
         * @return    pointer to acquired object, nullptr if failed.
         */
        ObjectPtrType tryAcquireObjectImpl() {

          ObjectPtrType sptr = nullptr;  // default is a null ptr.

          // okay to add before compare, since object pool is long lasting.
          int64_t prev_size = this->size_in_use.fetch_add(1, std::memory_order_acq_rel);  //reserve
          if (prev_size >= this->capacity) {
            this->size_in_use.fetch_sub(1, std::memory_order_release);
            // leave ptr as nullptr.
            if (this->isUnlimited()) {
              ERRORF("ERROR: pool is full but should be unlimited. prev size %lu.", prev_size);
            }
          } else {
            // try pop one.
            auto reuse = this->available.tryPop();
            std::atomic_thread_fence(std::memory_order_release);

            // if successful pop
            if (reuse.first) {
              // use it
              sptr = reuse.second;
            } else {
              // available is likely empty. allocate a new one.
              sptr = new ObjectType();
            }

          }

          return sptr;
        }


        /**
         * @brief     Release a buffer back to pool.  Thread local.
         * @details   clears object before release back into available.
         *        the ptr is checked against the thread's in-use set only.
         *        if the wrong thread releases the ptr, it will fail with false.
         *
         *        a buffer may be returned by a different thread.
         * @param ptr  ptr to object.
         * @param thread_id  id of thread to release object into.  -1 implies current omp thread.
         * @return    true if sucessful release. note if wrong thread releases then false is returned.
         */
        bool releaseObjectImpl(ObjectPtrType ptr) {

          if (!ptr)
          {
            WARNINGF("WARNING: pool releasing a nullptr.");
            return true;
          }

          bool res = false;            // nullptr would not be in in_use.

          auto v = this->size_in_use.fetch_sub(1, std::memory_order_acq_rel);
          if (v <= 0) {
            this->size_in_use.fetch_add(1, std::memory_order_release);
            return false;
          }

          // only put back in available queue if it was in use.
          // now make object available.  make sure push_back is done one thread at a time.
          res = this->available.tryPush(ptr);
          if (!res) {
            std::atomic_thread_fence(std::memory_order_acquire);
            delete ptr;
          }
          std::atomic_thread_fence(std::memory_order_release);

          return res;
        }


      public:
        /**
         * @brief     construct a Object Pool with buffers of capacity _buffer_capacity.
         *   the number of buffers in the pool is set to _pool_capacity, not _pool_capacity * _num_threads.
         *
         * @param _num_threasd		   the number of threads using this pool.
         * @param _pool_capacity       number of buffers in the pool
         * @param _buffer_capacity     size of the individual buffers
         */
        explicit ObjectPool(const int64_t _pool_capacity = std::numeric_limits<int64_t>::max()) :
        BaseType(_pool_capacity) {}

        /**
         * @brief     default destructor
         */
        virtual ~ObjectPool() {};

    };



    /**
     * @class     ObjectPool
     * @brief     ObjectPool manages a set of reusable, in-memory buffers.  In particular, it can limit the amount of memory usage.
     * @details   The class is templated to provide thread-safe ObjectPools, managing thread-safe or unsafe Objects
     *
     *            When the pool is exhausted, Acquire function returns false.
     *
     *      This template specialization is meant to support the use case where each thread has its own
     *      	pool that is not shared with any other thread.
     *
     *     this class assumes that while there may be multiple threads acquiring objects,
     *     but there is a SINGLE thread to release the object back.
     *        with this assumption, the pool can be lockfree, since the "in-use" set is unnecessary,
     *        and the available queue is already lockfree.
     *
     *        The releasing thread essentially acts as an "in-use" set.
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
    template<class T >
    class ObjectPool<T, bliss::concurrent::LockType::NONE, false> :
    public ObjectPoolBase<ObjectPool<T, bliss::concurrent::LockType::NONE, false> >
    {

      protected:

        using PoolType = ObjectPool<T, bliss::concurrent::LockType::NONE, false>;
        using BaseType = ObjectPoolBase<PoolType>;

        friend BaseType;

      public:
        /// object type
        using ObjectType = typename BaseType::ObjectType;

        /// object pointer type  (prefer shared pointer as c++11 allows atomic operations, however gcc does not support it.)
        using ObjectPtrType = typename BaseType::ObjectPtrType;


      protected:

        /**
         * @brief     default copy constructor is deleted.
         * @param other   source ObjectPool object to copy from.
         */
        explicit ObjectPool(const PoolType& other) = delete;

        /**
         * @brief     default copy assignment operator is deleted.
         * @param other   source ObjectPool object to copy from.
         * @return        self, with member variables copied from other.
         */
        PoolType& operator=(const PoolType& other) = delete;

        /// move constructor is deleted.
        explicit ObjectPool(PoolType&& other) = delete;
        /// move assignment operator is deleted.
        PoolType& operator=(PoolType&& other) = delete;


        /// clear all allocated objects.
        void clearStorageImpl() {
          this->size_in_use = 0;

          // clear the available queue.
          while (!this->available.empty()) {
            delete this->available.front();
            this->available.pop_front();
          }
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
        void resetImpl() {
          this->size_in_use = 0;

        }


        /**
         * @brief     Get the next available Object.  Mutex or Spin Locked.
         * @details   if none available, and size is below capacity or pool is unlimited,
         *            a new object is created.
         * @return    pointer to acquired object, nullptr if failed.
         */
        ObjectPtrType tryAcquireObjectImpl() {

          ObjectPtrType sptr = nullptr;  // default is a null ptr.

          // first check if we are exceeding capacity.
          if (this->size_in_use >= this->capacity) {
            if (this->isUnlimited()) {
              WARNINGF("ERROR: pool is full but should be unlimited. size %lu.", this->size_in_use);
            }  // else limited size, so return nullptr.
          } else {  // has room.  get one.


            // get an object
            if (this->available.empty()) {  // empty. allocate one.
              sptr = new ObjectType();
            } else {

              sptr = this->available.front();
              this->available.pop_front();
            }

            assert(sptr != nullptr);

            ++this->size_in_use;

          }

          return sptr;
        }


        /**
         * @brief     Release a buffer back to pool.  Threadsafe via mutex lock
         * @details   clears object before release back into available.
         * @param ptr weak_ptr to object.
         */
        bool releaseObjectImpl(ObjectPtrType ptr) {

          if (!ptr) {
            WARNINGF("WARNING: pool releasing a nullptr.");
            return true;
          }

          if (this->size_in_use == 0)
            return false;

          // only put back in available queue if it was in use.
          // now make object available.  make sure push_back is done one thread at a time.
          this->available.emplace_back(ptr);
          this->size_in_use -= 1;

          return true;
        }


      public:
        /**
         * @brief     construct a Object Pool with buffers of capacity _buffer_capacity.
         *   the number of buffers in the pool is set to _pool_capacity.
         *
         * @param _pool_capacity       number of buffers in the pool
         */
        explicit ObjectPool(const int64_t _pool_capacity = std::numeric_limits<int64_t>::max()) :
        BaseType(_pool_capacity) {};


        /**
         * @brief     default destructor
         */
        virtual ~ObjectPool() {};


    };


  } /* namespace concurrent */
} /* namespace bliss */

#endif /* UNREFERENCED_OBJECTPOOL_HPP_ */
