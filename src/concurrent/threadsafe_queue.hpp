/**
 * @file		threadsafe_queue.hpp
 * @ingroup concurrent
 * @author	Tony Pan <tpan7@gatech.edu>
 * @brief   Thread Safe Queue Implementation
 * @details this header file contains the templated implementation of a thread-safe queue.  this class is used for MPI buffer management.
 *
 *    This class uses moodycamel's concurrent queue.
 *      can't use boost's lockfree queue:  boost expect T to have copy constructor, trivial assignment operator, and trivial destructor.
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef THREADSAFE_QUEUE_HPP_
#define THREADSAFE_QUEUE_HPP_

#include "config/relacy_config.hpp"

#include <limits>
#include <stdexcept>
#include <xmmintrin.h>   // for _mm_paush

#include "concurrent/concurrent.hpp"

namespace bliss
{
  namespace concurrent
  {
    /// default template type for ThreadSafeQueue, user facing.
    template <typename T, bliss::concurrent::LockType LT>
    class ThreadSafeQueue;


    /**
     * @class     bliss::concurrent::ThreadSafeQueue
     * @brief     a multi-producer, multi-consumer thread safe queue with an optional capacity
     * @details   This class is a wrapper for the moodycamel::concurrentqueue::ConcurrentQueue class.
     *            https://github.com/cameron314/concurrentqueue
     *
     *          This class inherits the following properties:
     *            growable
     *            supports move semantics
     *            c++11 compatible
     *            threadsafe
     *            multiple producer, multiple consumer
     *
     *          The wrapper provides additional capabilities including:
     *            capacity limit
     *            disabling and enabling push operation (to "pause" the queue)
     *
     *          DUE TO LOCKS, this is NOT fast.  but may be fast enough for MPI buffer management.
     */

    namespace detail
    {
      /// default base type for ThreadSafeQueue Base
      template<typename Derived>
      class ThreadSafeQueueBase;

    template <typename T, bliss::concurrent::LockType LT>
    class ThreadSafeQueueBase< ThreadSafeQueue<T, LT> >
    {
    	//static_assert(false, "ConcurrentQueue is currently having data race problem resulting in missing entry sometimes, manifesting in deadlock when control message is missing.");
      public:

        /// maximum possible size of a thread safe queue.  initialized to maximum size_t value.
        static constexpr int64_t MAX_SIZE = std::numeric_limits<int64_t>::max();
        /// flag indicating the queue is disabled.
        static constexpr int64_t DISABLED = std::numeric_limits<int64_t>::lowest();


      protected:

        using Derived = ThreadSafeQueue<T, LT>;
        using Base = ThreadSafeQueueBase< Derived >;


        /// capacity of the queue.  if set to std::numeric_limits<int64_t>::max() indicates unlimited size queue
        mutable int64_t capacity;

        /// size encodes 2 things:  sign bit encodes whether a calling thread can push into this queue.  use when suspending or terminating a queue.  rest is size of current queue.
        typename std::conditional<LT == bliss::concurrent::LockType::NONE,
            VAR_T(int64_t),
            std::atomic<int64_t> >::type size;


        /**
         * normal constructor allowing the caller to specify an optional capacity parameter.
         * @param _capacity   The maximum capacity for the thread safe queue.
         */
        explicit ThreadSafeQueueBase(const size_t _capacity = static_cast<size_t>(MAX_SIZE)) :
              capacity(static_cast<int64_t>(_capacity))
        {
          assert(_capacity <= static_cast<size_t>(MAX_SIZE));
          if (capacity == 0)
            throw std::invalid_argument("ThreadSafeQueueBase constructor parameter capacity is given as 0");

          if (LT == bliss::concurrent::LockType::NONE)
            VAR(size) = 0;
          else
            size.store(0, std::memory_order_relaxed);
        };


        /**
         * move constructor, DISABLED
         * @param other   the source ThreadSafeQueueBase from which to copy.
         */
        DELETED_FUNC_DECL(explicit ThreadSafeQueueBase(Base&& other));

        /**
         * move assignment operator.  disabled.
         * @param other   the source ThreadSafeQueueBase from which to copy.
         * @return
         */
        DELETED_FUNC_DECL(Base& operator=(Base&& other));

        /**
         * copy constructor, DISABLED
         * @param other   the source ThreadSafeQueueBase from which to copy.
         */
        DELETED_FUNC_DECL(explicit ThreadSafeQueueBase(const Base& other));

        /**
         * copy assignment operator.  disabled.
         * @param other   the source ThreadSafeQueueBase from which to copy.
         * @return
         */
        DELETED_FUNC_DECL(Base& operator=(const Base& other));

        /// internal, get value in size variable.  not thread safe
        template<bliss::concurrent::LockType L = LT>
        inline const typename std::enable_if<L != bliss::concurrent::LockType::NONE, int64_t>::type
        getRawSize() const {
          return size.load(std::memory_order_acquire);
        }

        /// internal, get value in size variable.  not thread safe.
        template<bliss::concurrent::LockType L = LT>
        inline const typename std::enable_if<L == bliss::concurrent::LockType::NONE, int64_t>::type
        getRawSize() const {
          std::atomic_thread_fence(std::memory_order_acquire);
          return VAR(size);
        }
        
      public:

        /// default destructor
        virtual ~ThreadSafeQueueBase() {};

        /**
         * get the capacity of the thread safe queue
         * @return    capacity of the queue
         */
        inline const size_t getCapacity() const {
          return static_cast<size_t>(capacity);
        }

        /// check if queue is fixed size or growable.
        inline bool isFixedSize() const {
        	return capacity < MAX_SIZE;
        }

        /**
         * check if the thread safe queue is full.
         * @return    boolean - whether the queue is full.
         */
        inline bool isFull() const {
          return (capacity < MAX_SIZE) && (getSize() >= getCapacity());
        }

        /**
         * check if the thread safe queue is empty
         * @return    boolean - whether the queue is empty.
         */
        inline bool isEmpty() const
        {
          return getSize() == 0UL;
        }

        /**
         * get the current size of the queue
         * @return    the current size of the queue
         */
        inline const size_t getSize() const
        {
          return static_cast<size_t>(getRawSize<LT>() & MAX_SIZE);   // size is atomic, so don't need strong memory ordering itself.
        }

        /**
         * clears the queue of all contents (discarding).
         */
        void clear() {
          static_cast<Derived*>(this)->clearImpl();
        }


        /**
         * set the queue to accept new elements
         */
        inline void enablePush() {
          static_cast<Derived*>(this)->enablePushImpl();
        }

        /**
         * set the queue to disallow insertion of new elements.
         */
        inline void disablePush() {
          static_cast<Derived*>(this)->disablePushImpl();
        }

        /**
         * check if the thread safe queue can accept new elements. (full or not)
         * @return    boolean - queue insertion allowed or not.
         */
        inline bool canPush() {
          // pushing thread does not immediately care what other threads did before this
          // size is >= 0 (not disabled), and less than capacity.
          // note that if we reinterpret_cast size to size_t, then disabled (< 0) will have MSB set to 1, so > max(int64_t).
          return getRawSize<LT>() >= 0L;   // int highest bit set means negative, and means cannot push
        }

        /**
         * @brief    check if the thread safe queue can produce an element now or in near future.
         * @details  ThreadSafeQueueBase can pop only if it has elements in the base queue. or
         *            if additional items can be pushed in.
         * @return    boolean - queue pop is allowed or not.
         */
        inline bool canPop() {
          // canPop == first bit is 0, OR has some elements (not 0 for remaining bits).  so basically, not 1000000000...
          // popping thread should have visibility of all changes to the queue
          return getRawSize<LT>() != DISABLED;
        }

      protected:
        /**
         * check if the thread safe queue can accept new elements. (full or not)
         * @return    boolean - queue insertion allowed or not.
         */
        inline bool canPushAndHasRoom() {
          // if we reinterpret this number as a uint64_t, then we only need to check less than capacity, since now highest bits are all way higher.
          int64_t v = getRawSize<LT>();
          return reinterpret_cast<uint64_t&>(v) < getCapacity();
        }

      public:

        /**
         * @brief     Non-blocking - pushes an element by constant reference (copy).
         * Returns true only if push was successful.
         * If queue is full, or if queue is not accepting new element inserts, return false.
         *    Does not modify the element if cannot push onto queue.
         *
         * @param data    data element to be pushed onto the thread safe queue
         * @return        whether push was successful.
         */
        bool tryPush (T const& data) {
          return static_cast<Derived*>(this)->tryPushImpl(data);
        }

        /**
         * @brief Non-blocking - pushes an element by constant reference (move).
         * Returns true only if push was successful.
         * If queue is full, or if queue is not accepting new element inserts, return false.
         *    Does not modify the element if cannot push onto queue.
         *
         * @details  if move fails, concurrentqueue does NOT touch data (also uses std::forward, not move),
         *          so can return the results.
         * @note     requires move constuctor CLEAR OUT OLD, since that is semantically correct.  else destuction
         *            of old object could invalidate the moved to object.
         *
         * @param data    data element to be pushed onto the thread safe queue
         * @return        whether push was successful.
         */
        std::pair<bool, T> tryPush (T && data) {
          return static_cast<Derived*>(this)->tryPushImpl(std::forward<T>(data));
        }

        /**
         * @brief Semi-blocking - pushes an element by constant reference (copy).
         * Returns true only if push was successful.
         * If queue is full, wait until space becomes available.
         * If queue is not accepting new element inserts, return false immediately.
         *    Does not modify the element if cannot push onto queue.
         *
         * this function signature is useful where the application would like to have guaranteed insertion without busy-waiting,
         * yet needs a way to lift the block when the queue has been marked by another thread as terminating (not accepting new elements)
         *
         * @param data    data element to be pushed onto the thread safe queue
         * @return        whether push was successful.
         */
        bool waitAndPush (T const& data) {
          return static_cast<Derived*>(this)->waitAndPushImpl(data);
        }

        /**
         * @brief Semi-blocking - pushes an element by constant reference (move).
         * Returns true only if push was successful.
         * If queue is full, wait until space becomes available.
         * If queue is not accepting new element inserts, return false.
         *    Does not modify the element if cannot push onto queue.
         *
         * this function signature is useful where the application would like to have guaranteed insertion without busy-waiting,
         * yet needs a way to lift the block when the queue has been marked by another thread as terminating (not accepting new elements)
         *
         * @param data    data element to be pushed onto the thread safe queue
         * @return        whether push was successful.
         */
        std::pair<bool, T> waitAndPush (T && data) {
          return static_cast<Derived*>(this)->waitAndPushImpl(std::forward<T>(data));
        }




        /**
         * @brief	 Non-blocking - remove the first element in the queue and return it to the calling thread.
         * returns a std::pair with first element indicating the success/failure of the Pop operation,
         * and second is the popped queue element, if pop were successful.
         *
         * Function fails when the queue is empty, returning false.
         * Function returns success and retrieved element, regardless of whether the thread safe queue is accepting new inserts.
         *
         * @details uses move assignment operator internally.  requires Move assignment operator to clear the source, else destructor could invalidate current one.
         *
         * @return    std::pair with boolean (successful pop?) and an element from the queue (if successful)
         */
        std::pair<bool, T> tryPop() {
          return static_cast<Derived*>(this)->tryPopImpl();
        }

        /**
         * @brief Semi-blocking - remove the first element in the queue and return it to the calling thread.
         * Returns a std::pair with first element indicating the success/failure of the Pop operation,
         *  and second is the popped queue element, if pop were successful.
         *
         * Function will wait for some element to be available for pop from the queue.
         * Function will return when it has retrieved some data from queue, or when if it is notified that the thread safe queue has terminated.
         * Returns false if queue is terminated (no new elements) and flushed.
         *
         * This function signature is useful where the application would like to have guaranteed data retrieval from the queue without busy-waiting,
         * yet needs a way to lift the block when the queue has been marked by another thread as terminated and flushed
         *   (no more new elements to be inserted and queue is empty.)
         *
         * @return    std::pair with boolean (successful pop?) and an element from the queue (if successful)
         */
        std::pair<bool, T> waitAndPop() {
          return static_cast<Derived*>(this)->waitAndPopImpl();
        }

    };

    /**
     * static templated MAX_SIZE definition.
     */
    template<typename T, bliss::concurrent::LockType LT> constexpr int64_t ThreadSafeQueueBase<ThreadSafeQueue<T, LT> >::MAX_SIZE;

    }

  } /* namespace concurrent */
} /* namespace bliss */

#endif /* THREADSAFE_QUEUE_HPP_ */
