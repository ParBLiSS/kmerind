/**
 * @file		lockfree_queue.hpp
 * @ingroup concurrent
 * @author	Tony Pan <tpan7@gatech.edu>
 * @brief   Thread Safe Queue Implementation
 * @details this header file contains the templated implementation of a thread-safe queue.  this class is used for MPI buffer management.
 *
 *    This class uses moodycamel's concurrent queue.
 *      can't use boost's lockfree queue:  boost expect T to have copy constuctor, tribial assignment operator, and trivial destuctor.
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef LOCKFREE_QUEUE_HPP_
#define LOCKFREE_QUEUE_HPP_

#include "concurrentqueue/concurrentqueue.h"

#include "concurrent/threadsafe_queue.hpp"

namespace bliss
{
  namespace concurrent
  {
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
     *            lockfree
     *            multiple producer, multiple consumer
     *
     *          The wrapper provides additional capabilities including:
     *            capacity limit
     *            disabling and enabling push operation (to "pause" the queue)
     *
     *          DUE TO LOCKS, this is NOT fast.  but may be fast enough for MPI buffer management.
     */
    template <typename T>
    class ThreadSafeQueue<T, bliss::concurrent::LockType::LOCKFREE> :
    public detail::ThreadSafeQueueBase<ThreadSafeQueue<T, bliss::concurrent::LockType::LOCKFREE> >
    {
    	//static_assert(false, "ConcurrentQueue is currently having data race problem resulting in missing entry sometimes, manifesting in deadlock when control message is missing.");
    
      protected:

        using Derived = ThreadSafeQueue<T, bliss::concurrent::LockType::LOCKFREE>;
        using Base = detail::ThreadSafeQueueBase<Derived>;

        /// underlying lockfree queue
        moodycamel::ConcurrentQueue<T> q;

        /**
         * copy constructor, DISABLED
         * @param other   the source ThreadSafeQueue from which to copy.
         */
        DELETED_FUNC_DECL(explicit ThreadSafeQueue(const Derived& other));

        /**
         * copy assignment operator.  disabled.
         * @param other   the source ThreadSafeQueue from which to copy.
         * @return
         */
        DELETED_FUNC_DECL(Derived& operator=(const Derived& other));

      public:

        /**
         * normal constructor allowing the caller to specify an optional capacity parameter.
         * @param _capacity   The maximum capacity for the thread safe queue.
         */
        explicit ThreadSafeQueue(const size_t _capacity = static_cast<size_t>(Base::MAX_SIZE)) :
          Base(_capacity), q( _capacity == static_cast<size_t>(Base::MAX_SIZE) ? 128 : _capacity)
        {};

        /**
         * move constructor.  mutex locks the src ThreadSafeQueue first before delegating to the private constructor.
         * @param other   the source ThreadSafeQueue from which to move.
         */
        explicit ThreadSafeQueue(Derived&& other) {
          //std::atomic_thread_fence(std::memory_order_acquire);
          this->size.exchange(other.size.exchange(Base::DISABLED, std::memory_order_acq_rel), std::memory_order_relaxed);
          this->capacity = other.capacity;  other.capacity = 0;
          q = std::move(other.q);
          std::atomic_thread_fence(std::memory_order_release);
        }

        /**
         * move assignment operator.  locks both the src and destination ThreadSafeQueue before performing the move.
         * @param other   the source ThreadSafeQueue from which to move.
         * @return
         */
        Derived& operator=(Derived&& other) {

          //std::atomic_thread_fence(std::memory_order_acquire);
          this->size.exchange(other.size.exchange(Base::DISABLED, std::memory_order_acq_rel), std::memory_order_relaxed);
          this->capacity = other.capacity; other.capacity = 0;
          q = std::move(other.q);
          std::atomic_thread_fence(std::memory_order_release);

          return *this;
        }

        /**
         * clears the queue of all contents (discarding).
         */
        void clearImpl()
        {
          // first clear the size, so threads' don't try to dequeue
         this->size.fetch_and(Base::DISABLED, std::memory_order_acq_rel);  // keep the push bit, and set size to 0

          // dequeue
          T val;
          while (q.try_dequeue(val)) ;
          std::atomic_thread_fence(std::memory_order_release);
        }

        /**
         * set the queue to accept new elements
         */
        inline void enablePushImpl() {
          // before enablePush, queue is probably empty or no one is writing to it.  so prior side effects don't need to be visible right away.
          // size itself is atomic
          this->size.fetch_and(Base::MAX_SIZE, std::memory_order_relaxed);  // clear the push bit, and leave size as is
        }

        /**
         * set the queue to disallow insertion of new elements.
         */
        inline void disablePushImpl() {
          // before disable push, should make sure that all writes are visible to all other threads, so release here.
          this->size.fetch_or(Base::DISABLED, std::memory_order_relaxed);   // set the push bit, and leave size as is.
        }


        /**
         * @brief     Non-blocking - pushes an element by constant reference (copy).
         * Returns true only if push was successful.
         * If queue is full, or if queue is not accepting new element inserts, return false.
         *    Does not modify the element if cannot push onto queue.
         *
         * @param data    data element to be pushed onto the thread safe queue
         * @return        whether push was successful.
         */
        bool tryPushImpl (T const& data) {

//          // code below produces race condition causes extra items to be pushed.
//          int64_t v = size.load(std::memory_order_relaxed);
//          if ( reinterpret_cast<uint64_t&>(v) < static_cast<uint64_t>(this->capacity) ) {
//            if (q.enqueue(data)) {
//              size.fetch_add(1, std::memory_order_relaxed);
//              return true;
//            }
//          }

          bool res = false;
          int64_t v = this->size.fetch_add(1, std::memory_order_relaxed);  // reserve
          if (reinterpret_cast<uint64_t&>(v) < static_cast<uint64_t>(this->capacity)) {

            // between 0 and capacity - 1.  enqueue.
            std::atomic_thread_fence(std::memory_order_acq_rel);
            res = q.enqueue(data);
            std::atomic_thread_fence(std::memory_order_release);

            if (res) {
              return true;
            } // else  enqueue failed.
          } // else disabled or full
			    this->size.fetch_sub(1, std::memory_order_relaxed);

          // at this point, if success, data should be empty.  else data should be untouched.
          return false;

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
        std::pair<bool, T> tryPushImpl (T && data) {
//          // code below race condition causes extra items to be pushed.
//          int64_t v = this->size.load(std::memory_order_relaxed);
//          if ( reinterpret_cast<uint64_t&>(v) < static_cast<uint64_t>(this->capacity) ) {
//            if (q.enqueue(std::forward<T>(data))) {   // depend on q not messing up data if enqueue fails
//              this->size.fetch_add(1, std::memory_order_relaxed);
//              return std::move(std::make_pair(true, std::forward<T>(data)));  // data already moved, so okay.
//            }
//          }
          std::pair<bool, T> res(false, T());

          int64_t v = this->size.fetch_add(1, std::memory_order_relaxed);
          if (reinterpret_cast<uint64_t&>(v) < static_cast<uint64_t>(this->capacity)) {

            // else between 0 and capacity - 1.  enqueue.
            std::atomic_thread_fence(std::memory_order_acq_rel);
            res.first = q.enqueue(std::forward<T>(data));
            std::atomic_thread_fence(std::memory_order_release);

            if (res.first) {
              return std::move(res);
            } // else  enqueue failed. decrement and try again.
              
          } // else disabled or full
					this->size.fetch_sub(1, std::memory_order_relaxed);
					
          // at this point, if success, data should be empty.  else data should be untouched.
          res.second = std::move(data);
          return std::move(res);
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
        bool waitAndPushImpl (T const& data) {

          int64_t v = 0L;
          bool res = false;

          do {  // loop forever, unless queue is disabled, or insert succeeded.

            v = this->size.load(std::memory_order_relaxed);

            if (v < 0L) { // disabled so return false
              res = false;
              break;
            } else if (v < this->capacity) { // else between 0 and capacity - 1.  enqueue.
              res = this->size.compare_exchange_strong(v, v+1, std::memory_order_relaxed);
              if (res) {
                std::atomic_thread_fence(std::memory_order_acq_rel);
                res = q.enqueue(data);
                std::atomic_thread_fence(std::memory_order_release);

                if (res) {
                  break;
                } else {// else  enqueue failed. decrement and try again.
                  this->size.fetch_sub(1, std::memory_order_relaxed);
                }
              } // else reservation failed.
            }  // else over capacity

            // pause if over capacity, or failed to reserve, or failed to insert.
            _mm_pause();
          } while (true);

          return res;  // disabled so return false

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
        std::pair<bool, T> waitAndPushImpl (T && data) {
          std::pair<bool, T> res(false, T());

          int64_t v = 0L;

          do {  // loop forever, unless queue is disabled, or insert succeeded.

            v = this->size.load(std::memory_order_relaxed);

            if (v < 0L) { // disabled so return false
              res.first = false;
              res.second = std::move(data);
              break;
            } else if (v < this->capacity) { // else between 0 and capacity - 1.  enqueue.
              res.first = this->size.compare_exchange_strong(v, v+1, std::memory_order_relaxed);
              if (res.first) {
                std::atomic_thread_fence(std::memory_order_acq_rel);
                res.first = q.enqueue(std::forward<T>(data));
                std::atomic_thread_fence(std::memory_order_release);

                if (res.first) {
                  break;
                } else {// else  enqueue failed. decrement and try again.
                  this->size.fetch_sub(1, std::memory_order_relaxed);
                }
              } // else reservation failed.
            }  // else over capacity

            // pause if over capacity, or failed to reserve, or failed to insert.
            _mm_pause();
          } while (true);

          return std::move(res);  // disabled so return false

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
        std::pair<bool, T> tryPopImpl() {
          // pop dequeues first, then decrement count
          std::pair<bool, T> res(false, T());

          int64_t v = this->size.load(std::memory_order_relaxed);

					if ((v & Base::MAX_SIZE) == 0L) return res;  // nothing to pop

					// use compare_exchange_strong to ensure we are not going into negative territory.
					res.first = this->size.compare_exchange_strong(v, v-1, std::memory_order_relaxed);

					if (res.first) {
            std::atomic_thread_fence(std::memory_order_acq_rel);
            res.first = q.try_dequeue(res.second);
            std::atomic_thread_fence(std::memory_order_release);

            if (!res.first) {
              this->size.fetch_add(1, std::memory_order_relaxed);
            }

					}
          return res;
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
        std::pair<bool, T> waitAndPopImpl() {
          std::pair<bool, T> res(false, T());

          int64_t v = 0L;

          do {
            v = this->size.load(std::memory_order_relaxed);

            if (v == Base::DISABLED) {  // disabled and empty. return.
              res.first = false;
              break;
            } else if (v != 0L) {

              res.first = this->size.compare_exchange_strong(v, v-1, std::memory_order_relaxed);
              if (res.first)
              {
                //  else has some entries
                std::atomic_thread_fence(std::memory_order_acq_rel);
                res.first = q.try_dequeue(res.second);
                std::atomic_thread_fence(std::memory_order_release);

                if (res.first) {
                  break;
                } else {
                  // failed dequeue. need to increment and try again
                  v = this->size.fetch_add(1, std::memory_order_relaxed);
                }
              } // failed reserve.  try again
            } // else empty, so wait

            _mm_pause();
          } while (true);

          return res;

        }


    };


    template<typename T>
    using LockfreeQueue = bliss::concurrent::ThreadSafeQueue<T, bliss::concurrent::LockType::LOCKFREE>;



  } /* namespace concurrent */
} /* namespace bliss */

#endif /* LOCKFREE_QUEUE_HPP_ */
