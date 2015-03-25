/**
 * @file		object_pool.hpp
 * @ingroup concurrent
 * @author	Tony Pan <tpan7@gatech.edu>
 * @brief   this file defines in-memory Object Pool.
 * @details
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef OBJECTPOOL_HPP_
#define OBJECTPOOL_HPP_

#include <deque>
#include <unordered_set>
#include <stdexcept>

#include "utils/logging.h"
#include "concurrent/concurrent.hpp"

#include "config/relacy_config.hpp"

#include "concurrent/lockfree_queue.hpp"


// TODO: convert to an allocator.

namespace bliss
{
  namespace concurrent
  {
    /// generic template class for ObjectPoolBase class.
    template<class Derived>
    class ObjectPoolBase;

    /**
     * @brief ObjectPool class.  Allows reuse of previously allocated objects.  Also tracks objects that are currently in use.
     * @detail   supports various Locking strategies.  Acquire and release locking strategies do not need to be the same under some conditions.
     *          Supported LockTypes are MUTEX, SPINLOCK, THREADLOCAL, and NONE.
     *
     *          For MUTEX, SPINLOCK, and NONE, the Acquire and Release LockTypes must be the same.
     *
     *          For THREADLOCAL, Release LockType can be either THREADLOCAL or NONE.   In the first case, each thread tracks in-use objects in its own set.
     *            In the second case, the releasing thread is assume to be tracking the in-use objects so that an object is released only once.
     *
     *          There is currently no support for LOCKFREE ObjectPool because the requirement to track objects that are in use translates to needing a
     *            lockfree set, which we do not have at the moment.  one possibility is to punt the requirement to the calling threads.
     *
     *
     *          These requirements should be enforced by the template specializations.
     *
     */
    template<class T, bliss::concurrent::LockType LockType, bool Referencing = true>
    class ObjectPool;

    /// generic template class for ObjectPoolBase class, specialized for Derived type == ObjectPool and used for extracting ObjectPool's template parameters.
    template<class T, bliss::concurrent::LockType LockType, bool Referencing>
    class ObjectPoolBase< ObjectPool<T, LockType, Referencing> >
    {

      protected:

        using Derived = ObjectPool<T, LockType, Referencing>;

        /// object type
        using ObjectType = T;

        /// object pointer type  (prefer shared pointer as c++11 allows atomic operations, however gcc does not support it.)
        using ObjectPtrType = T*;


        //================== member variables
        /**
         * @brief     capacity of the ObjectPool (in number of Objects)
         */
        mutable int64_t                                                   capacity;

        /**
         * @brief     Internal queue of available Objects for immediate use.
         * @details   ThreadSafeQueue to ensure thread safety.
         */
        typename std::conditional<LockType == bliss::concurrent::LockType::LOCKFREE ||
            LockType == bliss::concurrent::LockType::THREADLOCAL,
            bliss::concurrent::ThreadSafeQueue<ObjectPtrType, bliss::concurrent::LockType::LOCKFREE>,
            std::deque<ObjectPtrType> >::type                             available;

        /**
         * @brief     Internal Set of objects in use..
         * @details   locked for threadsafe access.
         */
        typename std::conditional<LockType == bliss::concurrent::LockType::THREADLOCAL,
            std::vector<std::unordered_set<ObjectPtrType> >,
            std::unordered_set<ObjectPtrType> >::type                      in_use;

        /// number of objects currently in use.
        typename std::conditional<LockType == bliss::concurrent::LockType::LOCKFREE ||
            LockType == bliss::concurrent::LockType::THREADLOCAL,
            std::atomic<int64_t>,
            VAR_T(int64_t)>::type                                                size_in_use;


        //================ member methods

        /**
         * @brief     construct a Object Pool with buffers of capacity _buffer_capacity.
         *   the number of buffers in the pool is set to _pool_capacity.
         *
         * @param _pool_capacity       number of buffers in the pool
         * @param _buffer_capacity     size of the individual buffers
         */
        explicit ObjectPoolBase(const int64_t _pool_capacity = std::numeric_limits<int64_t>::max()) :
          capacity(_pool_capacity), available(), in_use(), size_in_use(0)
        {};

        /// destructor
        virtual ~ObjectPoolBase() {
          // delete all the objects
          this->clear_storage();
        }

        /// deleted default constructor, move and copy constructors and assignment operators.
        explicit DELETED_FUNC_DECL(ObjectPoolBase(ObjectPoolBase const &));
        explicit DELETED_FUNC_DECL(ObjectPoolBase(ObjectPoolBase &&));
        ObjectPoolBase& DELETED_FUNC_DECL(operator=(ObjectPoolBase const &));
        ObjectPoolBase& DELETED_FUNC_DECL(operator=(ObjectPoolBase &&));


        //-----------  these are support functions.  they should not lock.

        /// clear all allocated objects.
        void clear_storage() {
          static_cast<Derived*>(this)->clearStorageImpl();
        }

        template<bliss::concurrent::LockType LT = LockType>
        inline typename std::enable_if<LT != bliss::concurrent::LockType::LOCKFREE &&
        LT != bliss::concurrent::LockType::THREADLOCAL, int64_t>::type getSizeInUse() const {
          return VAR(size_in_use);
        }

        template<bliss::concurrent::LockType LT = LockType>
        inline typename std::enable_if<LT == bliss::concurrent::LockType::LOCKFREE ||
        LT == bliss::concurrent::LockType::THREADLOCAL, int64_t>::type getSizeInUse() const {
          return size_in_use.load(std::memory_order_relaxed);
        }

      public:

        /// type of lock for the pool.
        static const bliss::concurrent::LockType poolLT = LockType;

        /// check if the pool has unlimited capacity
        inline bool isUnlimited() const {
          return capacity == std::numeric_limits<int64_t>::max();
        }


        /**
         * @brief     Current size of the ObjectPool.  For debugging only.
         * @return  size,
         */
        inline int64_t getAvailableCount() const  {
          return (isUnlimited() ? capacity : (capacity - this->getSizeInUse<LockType>()));
        }


        /**
         * @brief     Current capacity of the ObjectPool
         * @return    capacity
         */
        inline const int64_t getCapacity() const {
          return capacity;
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
          static_cast<Derived*>(this)->resetImpl();
        }

        /**
         * @brief     Get the next available Object.  Mutex or Spin Locked.
         * @details   if none available, and size is below capacity or pool is unlimited,
         *            a new object is created.
         * @return    pointer to acquired object, nullptr if failed.
         */
        ObjectPtrType tryAcquireObject() {
          return static_cast<Derived*>(this)->tryAcquireObjectImpl();
        }

        /**
         * @brief     Release a buffer back to pool.  Threadsafe via mutex lock
         * @details   clears object before release back into available.
         * @param ptr weak_ptr to object.
         */
        bool releaseObject(ObjectPtrType ptr) {
          return static_cast<Derived*>(this)->releaseObjectImpl(ptr);
        }

    };


    // static const pool LockType variable.
    template<class T, bliss::concurrent::LockType LockType, bool Referencing>
    const bliss::concurrent::LockType ObjectPoolBase<ObjectPool<T, LockType, Referencing> >::poolLT;


  } /* namespace concurrent */
} /* namespace bliss */

#endif /* OBJECTPOOL_HPP_ */
