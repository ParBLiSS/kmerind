/**
 * @file		single_releaser_object_pool.hpp
 * @ingroup concurrent
 * @author	Tony Pan <tpan7@gatech.edu>
 * @brief   this file defines an in-memory Object Pool for the scenario where there is only 1 thread performing the obj release.
 * @details this is essentially a higher level allocator for large objects.
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

// TODO: convert this to a custom, threadsafe allocator?

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
     * @note      Each object should be acquired by a single thread.
     *            object may be used by multiple threads
     *            multiple threads can concurrently acquire
     *            the object MAY be released a single thread only.
     *
     *            It is important to address race conditions when multithreading, including
     *              likely loss of data (thread 1 appends data while thread 2 releases it;
     *              thread 2 and 3 release the same object, thread 1 successfully acquires before thread 3 releases object).
     *
     *            MessageBuffers class internally ensures that this does not happen.
     *
     *
     *  @tparam LockType    The thread safety property for the pool. Default to mutex.
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

        /// type of lock for the pool. (mutex or spinlock,or lockfree, or none).
        static const bliss::concurrent::LockType poolLT = LockType;

      protected:
        /**
         * @brief     capacity of the ObjectPool (in number of Objects)
         */
        mutable int64_t capacity;

        /**
         * @brief     Internal queue of available Objects for immediate use.
         * @details   ThreadSafeQueue to ensure thread safety.
         *            uniqueness is responsibility of releasing thread.
         */
        bliss::concurrent::ThreadSafeQueue<ObjectPtrType>  available;

        /// num objects in active use.
        std::atomic<int64_t>  size_in_use;

        /// mutex for thread safety
        mutable std::mutex mutex;

        /// clear all allocated objects.  NOT Thread Safe
        void clear_storage() {
          ObjectPtrType ptr;
          auto entry = available.tryPop();
          while (entry.first) {
             if (entry.second) {
               delete entry.second;
             }
             entry = available.tryPop();
          }

          // note that in-use information is held by releaser thread.
          size_in_use.store(0, std::memory_order_relaxed);

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
         *  the number of buffers in the pool is set to _pool_capacity.
         *
         * @param _pool_capacity       number of buffers in the pool
         * @param _buffer_capacity     size of the individual buffers
         */
        explicit ObjectPool(const int64_t _pool_capacity = std::numeric_limits<int64_t>::max()) :
          capacity(_pool_capacity),
          available(_pool_capacity), size_in_use(0)
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
         * @return    capacity,
         */
        const int64_t getCapacity() const {
          return capacity;
        }

        /// check if the pool has unlimited capacity
        inline const bool isUnlimited() const {
          return capacity == std::numeric_limits<int64_t>::max();
        }


        /**
         * @brief     Resets the entire ObjectPool: all Objects in pool are marked as released and available.
         *
         * @note      this method just sets in-use to 0.  the releasing thread may have objects to release.
         */
        void reset() {
          size_in_use.store(0, std::memory_order_release);
        }


        /**
         * @brief     Get the next available Object.  thread safe
         * @details   if none available, and size is below capacity or pool is unlimited,
         *            a new object is created.
         * @return    pointer to acquired object, nullptr if failed.
         */
        ObjectPtrType tryAcquireObject() {


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


            // try pop one.
            auto reuse = available.tryPop();

            // if successful pop
            if (reuse.first) {
              // use it
              sptr = reuse.second;
            } else {
              // available is likely empty. allocate a new one.
              sptr = new T();
            }

          }

          return sptr;

        }




    	/**
    	 * @brief     Release a buffer back to pool.  Single Thread Only..
    	 * @details   clears object before release back into available.
    	 * 				the ptr is NOT checked against any in-use set only.
    	 * 				if the ptr was not created by this pool, it may encounter
    	 * 				double free issues.
    	 * @param ptr  ptr to object.
    	 * @return		true if sucessful release. note if wrong thread releases then false is returned.
    	 */
        bool releaseObject(ObjectPtrType ptr) {

    		if (!ptr)
    		{
    			WARNINGF("WARNING: pool releasing a nullptr.");
    			return true;
    		}
          // now make object available.  make sure push_back is done one thread at a time.
          int64_t prev = size_in_use.fetch_sub(1, std::memory_order_relaxed);
          if (prev <= 0) {
        	  // failed release, but do not delete.
            size_in_use.fetch_add(1, std::memory_order_release);
            return false;
          } else {
            return available.tryPush(ptr);    // push it onto the queue.
          }

        }
    };


    /**
     * @class     ObjectPool
     * @brief     ObjectPool manages a set of reusable, in-memory buffers.  In particular, it can limit the amount of memory usage.
     * @details   The class is templated to provide thread-safe or unsafe ObjectPools, managing thread-safe or unsafe Objects
     *
     *            When the pool is exhausted, Acquire function returns false.
     *
     *     this class assumes single acquiring thread and single releasing thread.
     *        It is NOT thread safe.
     *
     *        The releasing thread essentially acts as an "in-use" set.
     *
     *      object life cycle:
     *              a thread calls pool's acquire() to get pointer to an allocated object
     *              application:  uses object, potentially by multiple threads
     *              when done, 1 thread releases the object to the pool.  repeated release of the same object pointer does nothing.
     *
     *
     * @note      Each object should be acquired by a single thread.
     *            object may be used by multiple threads
     *            multiple threads can concurrently acquire
     *            the object MAY be released a single thread only.
     *
     *  @tparam LockType    The thread safety property for the pool, default to NONE
     *  @tparam T           The object instance type.  object pool stores pointers to T instances.
     */
    template<class T >
    class ObjectPool<bliss::concurrent::LockType::NONE, T>
    {
      public:
        /// object type
        using ObjectType = T;
        /// object pointer type  (prefer shared pointer as c++11 allows atomic operations, but gcc does not support it.
        using ObjectPtrType = T*;

        /// type of lock for the pool. set to None
        static const bliss::concurrent::LockType poolLT = bliss::concurrent::LockType::NONE;


      protected:
        /**
         * @brief     capacity of the ObjectPool (in number of Objects)
         */
        mutable int64_t capacity;

        /**
         * @brief     Internal Set of available Objects for immediate use.
         * @details   Using std deque
         * 				uniqueness is responsibility of releasing thread.
         */
        std::deque<ObjectPtrType>                          available;

        /// num objects in active use.
        int64_t size_in_use;

        /// mutex for thread safety
        mutable std::mutex mutex;

        /// clear all allocated objects.  NOT Thread Safe
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
         * @brief     construct a Object Pool with buffers of capacity _buffer_capacity.
         *   the number of buffers in the pool is set to _pool_capacity.
         *
         * @param _pool_capacity       number of buffers in the pool
         * @param _buffer_capacity     size of the individual buffers
         */
        explicit ObjectPool(const int64_t _pool_capacity = std::numeric_limits<int64_t>::max()) :
          capacity(_pool_capacity),
          available(), size_in_use(0)
          {};


        /// move constructor is deleted.
        explicit ObjectPool(ObjectPool<poolLT, T>&& other) = delete;
        /// move assignment operator is deleted.
        ObjectPool<poolLT, T>& operator=(ObjectPool<poolLT, T>&& other) = delete;

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
         * @return  size
         */
        int64_t getAvailableCount() const  {
          return (isUnlimited() ? capacity : (capacity - size_in_use));
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
         * @brief     Resets the entire ObjectPool: all Objects in pool are marked as released and available.
         *
         * @note      this method just sets in-use to 0.  the releasing thread may have objects to release.
         */
        void reset() {
          size_in_use = 0;
        }


        /**
         * @brief     Get the next available Object.  single thread
         * @details   if none available, and size is below capacity or pool is unlimited,
         *            a new object is created.
         * @return    pointer to acquired object, nullptr if failed.
         */
        ObjectPtrType tryAcquireObject() {

          ObjectPtrType sptr = nullptr;  // default is a null ptr.

          if (size_in_use >= capacity) {
            // leave ptr as nullptr.
  			if (this->isUnlimited()) {
  				ERRORF("ERROR: pool is full but should be unlimited. prev size %lu.", size_in_use);
  			}
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
    	 * @brief     Release a buffer back to pool.  Single Thread Only..
    	 * @details   clears object before release back into available.
    	 * 				the ptr is NOT checked against any in-use set only.
    	 * 				if the ptr was not created by this pool, it may encounter
    	 * 				double free issues.
    	 * @param ptr  ptr to object.
    	 * @return		true if sucessful release. note if wrong thread releases then false is returned.
    	 */
        bool releaseObject(ObjectPtrType ptr) {
    		if (!ptr)
    		{
    			WARNINGF("WARNING: pool releasing a nullptr.");
    			return true;
    		}

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

    /// static variable for Pool LockType.
    template<bliss::concurrent::LockType LockType, class T>
    const bliss::concurrent::LockType ObjectPool<LockType, T>::poolLT;


  } /* namespace concurrent */
} /* namespace bliss */

#endif /* OBJECTPOOL_HPP_ */
