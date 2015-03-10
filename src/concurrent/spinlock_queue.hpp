/**
 * @file		spinlock_queue.hpp
 * @ingroup concurrent
 * @author	Tony Pan <tpan7@gatech.edu>
 * @brief   Thread Safe Queue Implementation
 * @details this header file contains the templated implementation of a thread safe queue.  this class is used for MPI buffer management.
 *
 *	this class is adapted from  http://www.justsoftwaresolutions.co.uk/threading/implementing-a-thread-safe-queue-using-condition-variables.html
 *            with SIGNIFICANT modifications
 *
 * Copyright (c) 2014 Georgia InOAstitute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef SPINLOCK_QUEUE_HPP_
#define SPINLOCK_QUEUE_HPP_

#include <deque>

#include "concurrent/threadsafe_queue.hpp"
#include "utils/logging.h"

namespace bliss
{
  namespace concurrent
  {
  /**
   * @class     bliss::concurrent::ThreadSafeQueue
   * @brief     a multi-producer, multi-consumer thread safe queue with an optional capacity
   * @details   this class is adapted from  http://www.justsoftwaresolutions.co.uk/threading/implementing-a-thread-safe-queue-using-condition-variables.html
   *            with SIGNIFICANT modifications
   *            1. use c++11 std::thread and atomic constructs
   *            2. incorporating move semantics.
   *            3. support a capacity limit.
   *            4. supports multiple producer, multiple consumer.
   *            5. allows blocking of the enqueue function, (useful when draining queue or when finishing the use of the queue).
   *
   *          note that this is NOT truly concurrent.  the class serializes parallel access.  move semantic minimizes the copy operations needed.
   *
   *          DUE TO SPINLOCKS, this is NOT fast.  but may be fast enough for MPI buffer management.
   *			this is faster than Mutex based version.
   */
    template <typename T>
    class ThreadSafeQueue<T, bliss::concurrent::LockType::SPINLOCK> :
      public detail::ThreadSafeQueueBase<ThreadSafeQueue<T, bliss::concurrent::LockType::SPINLOCK> >

    {
      protected:
        using Derived = ThreadSafeQueue<T, bliss::concurrent::LockType::SPINLOCK>;
        using Base = detail::ThreadSafeQueueBase< Derived >;

        /// spinlock for locking during enqueue, dequeue.
        mutable INIT_ATOMIC_FLAG(spinlock);

        /// internal queue, protected by spinlock.
        std::deque<T> q;

        /**
         * copy constructor, DISABLED
         * @param other   the source ThreadSafeQueue from which to copy.
         */
        DELETED_FUNC_DECL(explicit ThreadSafeQueue(const Derived& other));

        /**
         * copy assignment operator
         * @param other   the source ThreadSafeQueue from which to copy.
         * @return
         */
        DELETED_FUNC_DECL(Derived& operator=(const Derived& other));

      public:

        /**
         * normal constructor allowing the caller to specify an optional capacity parameter.
         * @param _capacity   The maximum capacity for the thread safe queue.
         */
        explicit ThreadSafeQueue(const size_t &_capacity = static_cast<size_t>(Base::MAX_SIZE)) :
          Base(_capacity), q() {};

        /**
         * move constructor.  mutex locks the src ThreadSafeQueue first before delegating to the private constructor.
         * @param other   the source ThreadSafeQueue from which to move.
         */
        explicit ThreadSafeQueue(Derived&& other) {
        	while (other.spinlock.test_and_set(std::memory_order_acq_rel));

    			std::swap(q, other.q);
		    	this->capacity = other.capacity; other.capacity = 0;
			    // ordering may seem strange, but see http://en.cppreference.com/w/cpp/atomic/memory_order re. operator ordering within same thread.
			    VAR(this->size) = VAR(other.size); VAR(other.size) = Base::DISABLED;
    			other.spinlock.clear(std::memory_order_release);
        };

        /**
         * move assignment operator.  locks both the src and destination ThreadSafeQueue before performing the move.
         * @param other   the source ThreadSafeQueue from which to move.
         * @return
         */
        Derived& operator=(Derived&& other) {

        	while (spinlock.test_and_set(std::memory_order_acq_rel));
        	while (other.spinlock.test_and_set(std::memory_order_acq_rel));

        	std::swap(q, other.q);
        	this->capacity = other.capacity; other.capacity = 0;
        	VAR(this->size) = VAR(other.size); VAR(other.size) = Base::DISABLED;

        	other.spinlock.clear(std::memory_order_release);
        	spinlock.clear(std::memory_order_release);

           return *this;
        }

        /**
         * clears the queue of all contents (discarding).
         */
        void clearImpl()
        {
          while(spinlock.test_and_set(std::memory_order_acq_rel));
          q.clear();
          // clear size
          VAR(this->size) &= Base::DISABLED;  // keep the push bit, and set size to 0
          spinlock.clear(std::memory_order_release);
        }


        /**
         * set the queue to accept new elements
         */
        inline void enablePushImpl() {
          while (spinlock.test_and_set(std::memory_order_acquire));
          VAR(this->size) &= Base::MAX_SIZE;  // clear the push bit, and leave size as is
          spinlock.clear(std::memory_order_release);
        }

        /**
         * set the queue to disallow insertion of new elements.
         */
        inline void disablePushImpl() {
          while (spinlock.test_and_set(std::memory_order_acquire));
          VAR(this->size) |= Base::DISABLED;  // clear the push bit, and leave size as is
          spinlock.clear(std::memory_order_release);
        }


        /**
         * @brief Non-blocking - pushes an element by constant reference (copy).
         * Returns true only if push was successful.
         * If queue is full, or if queue is not accepting new element inserts, return false.
         *    Does not modify the element if cannot push onto queue.
         *
         * @param data    data element to be pushed onto the thread safe queue
         * @return        whether push was successful.
         */
        bool tryPushImpl (T const& data) {
          bool res = false;
          while(spinlock.test_and_set(std::memory_order_acq_rel));    // critical block instructions confined here via acquire/release

          if (reinterpret_cast<uint64_t&>(VAR(this->size)) < static_cast<uint64_t>(this->capacity)) {
            // following 2 instructions cannot be reordered to above this line,
            // and the if check cannot be reordered to below this line
           // std::atomic_thread_fence(std::memory_order_acq_rel);

      			q.push_back(data);  // insert using predefined copy version of dequeue's push function
            res = true;
            ++VAR(this->size);
          } // else at capacity.
          spinlock.clear(std::memory_order_release);
          return res;
        }

        /**
         * @brief 	Non-blocking - pushes an element by constant reference (move).
         * Returns true only if push was successful.
         * If queue is full, or if queue is not accepting new element inserts, return false.
         *    Does not modify the element if cannot push onto queue.
         *
         * @param data    data element to be pushed onto the thread safe queue
         * @return        whether push was successful.
         */
        std::pair<bool, T> tryPushImpl (T && data) {
            std::pair<bool, T> res(false, T());

            while(spinlock.test_and_set(std::memory_order_acq_rel));   // see const ref version for memory order info.

            if (reinterpret_cast<uint64_t&>(VAR(this->size)) < static_cast<uint64_t>(this->capacity)) {
              //std::atomic_thread_fence(std::memory_order_acq_rel);

          	  q.emplace_back(std::forward<T>(data));  // insert using predefined copy version of dequeue's push function
              res.first = true;
              ++VAR(this->size);
            } else { // else at capacity.
              res.second = std::move(data);
            }

            spinlock.clear(std::memory_order_release);
            return std::move(res);
        }

        /**
         * Semi-blocking - pushes an element by constant reference (copy).
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
        bool waitAndPushImpl (T const& data) {

          int64_t v = 0L;
          bool res = false;

          do {  // loop forever, unless queue is disabled, or insert succeeded.
        	  while(spinlock.test_and_set(std::memory_order_acq_rel));

            v = VAR(this->size);

            if (v < 0L) { // disabled so return false
             	spinlock.clear(std::memory_order_release);
              break;
            } else if (v >= this->capacity) {   // over capacity.  wait a little
             	spinlock.clear(std::memory_order_release);
            	_mm_pause();
              continue;
            } else {
              // else between 0 and capacity - 1.  enqueue.
              q.push_back(data);  // successfully enqueued.
              res = true;
              ++VAR(this->size);

              spinlock.clear(std::memory_order_release);
              break;
            }
          } while (!res);

          return res;
        }

        /**
         * Semi-blocking - pushes an element by constant reference (move).
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
          	while(spinlock.test_and_set(std::memory_order_acq_rel));

            v = VAR(this->size);

            if (v < 0L) { // disabled so return false
              res.second = std::move(data);
            	spinlock.clear(std::memory_order_release);
              break;
            } else if (v >= this->capacity) {   // over capacity.  wait.
             	spinlock.clear(std::memory_order_release);
            	_mm_pause();
              continue;
            } else {

              // else between 0 and capacity - 1.  enqueue.
              q.emplace_back(std::forward<T>(data));  // successfully enqueued.
              res.first = true;
              ++VAR(this->size);

              spinlock.clear(std::memory_order_release);
              break;
            }
          } while (!res.first);

          return std::move(res);  // disabled so return false
        }


        /**
         * @brief	Non-blocking - remove the first element in the queue and return it to the calling thread.
         * returns a std::pair with first element indicating the success/failure of the Pop operation,
         * and second is the popped queue element, if pop were successful.
         *
         * Function fails when the queue is empty, returning false.
         * Function returns success and retrieved element, regardless of whether the thread safe queue is accepting new inserts.
         *
         * @return    std::pair with boolean (successful pop?) and an element from the queue (if successful)
         */
        std::pair<bool, T> tryPopImpl() {
          std::pair<bool, T> res(false, T());

          if (this->isEmpty()) return res;

          while (spinlock.test_and_set(std::memory_order_acq_rel));   // establish critical section

          if (!this->isEmpty()) {

        	  res.second = std::move(q.front());  // convert to movable reference and move-assign.
        	  q.pop_front();
            res.first = true;
            --VAR(this->size);
          }

          spinlock.clear(std::memory_order_release);

          return res;

        }

        /**
         * @brief	Semi-blocking - remove the first element in the queue and return it to the calling thread.
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
          bool breaking = false;

          do {

            while (spinlock.test_and_set(std::memory_order_acq_rel));

            v = VAR(this->size);

//            if (v < 0) DEBUGF("-");
            if (v == Base::DISABLED) {  // disabled and empty. return.
//              DEBUGF("-");
              breaking = true;
            } else if (v == 0L) {  // empty but not disabled. so wait
//              DEBUGF("%ld", v);
            } else {

//            DEBUGF("popping...");

              //  else has some entries
              res.second = std::move(q.front());  // convert to movable reference and move-assign.
              q.pop_front();
              
              res.first = true;
              --VAR(this->size);
              breaking = true;
            }
            spinlock.clear(std::memory_order_release);

            if (breaking) break;
          
            _mm_pause();
          } while (!res.first);

          return res;
        }



    };

    template<typename T>
    using SpinLockQueue = ThreadSafeQueue<T, bliss::concurrent::LockType::SPINLOCK>;



  } /* namespace concurrent */
} /* namespace bliss */

#endif /* SPINLOCK_QUEUE_HPP_ */
