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
        static_assert(false, "TSAN reports race conditions in access to deque (front() and push_back())");

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
        explicit ThreadSafeQueue(const size_t _capacity = Base::MAX_SIZE) :
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
		    	this->size.exchange(other.size.exchange(Base::DISABLED, std::memory_order_relaxed), std::memory_order_relaxed);
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
        	this->size.exchange(other.size.exchange(Base::DISABLED, std::memory_order_relaxed), std::memory_order_relaxed);

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
          this->size.fetch_and(Base::DISABLED, std::memory_order_relaxed);  // keep the push bit, and set size to 0
          q.clear();
          spinlock.clear(std::memory_order_release);
        }


        /**
         * set the queue to accept new elements
         */
        inline void enablePushImpl() {
          this->size.fetch_and(Base::MAX_SIZE, std::memory_order_relaxed);  // clear the push bit, and leave size as is        }
        }

        /**
         * set the queue to disallow insertion of new elements.
         */
        inline void disablePushImpl() {
          // before disable push, should make sure that all writes are visible to all other threads, so release here.
          this->size.fetch_or(Base::DISABLED, std::memory_order_relaxed);   // set the push bit, and leave size as is.
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
          size_t v = this->size.fetch_add(1, std::memory_order_relaxed);



          if (v < this->capacity) {
            // following 2 instructions cannot be reordered to above this line,
            // and the if check cannot be reordered to below this line
           // std::atomic_thread_fence(std::memory_order_acq_rel);

      			q.push_back(data);  // insert using predefined copy version of dequeue's push function
            res = true;
          } else {
            // else at capacity or disabled.
            this->size.fetch_sub(1, std::memory_order_relaxed); // failed enqueue, decrement size.
          }
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
            size_t v = this->size.fetch_add(1, std::memory_order_relaxed);
            if (v < this->capacity) {

          	  q.emplace_back(std::forward<T>(data));  // insert using predefined copy version of dequeue's push function
              res.first = true;

            } else { // else at capacity or disabled.
              res.second = std::move(data);
              this->size.fetch_sub(1, std::memory_order_relaxed); // failed enqueue, decrement size.
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

          size_t v = 0UL;
          bool res = false;

          do {  // loop forever, unless queue is disabled, or insert succeeded.

            v = this->size.load(std::memory_order_relaxed);

            if (v >= Base::DISABLED) { // disabled so return false
              res = false;
              break;
            } else if (v < this->capacity) { // else between 0 and capacity - 1.  enqueue.
              while(spinlock.test_and_set(std::memory_order_acq_rel));
              res = this->size.compare_exchange_strong(v, v+1, std::memory_order_relaxed);
              if (res) {
                q.push_back(data);
                spinlock.clear(std::memory_order_release);
                break;
              } // failed reservation.  wait
              spinlock.clear(std::memory_order_release);
            } // else over capacity.  wait.

            _mm_pause();
          } while (true);

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

          size_t v = 0UL;

          do {  // loop forever, unless queue is disabled, or insert succeeded.

          	v = this->size.load(std::memory_order_relaxed);

            if (v >= Base::DISABLED) { // disabled so return false
              res.first = false;
              res.second = std::move(data);
              break;
            } else if (v < this->capacity) { // else between 0 and capacity - 1.  enqueue.
              while(spinlock.test_and_set(std::memory_order_acq_rel));
              res.first = this->size.compare_exchange_strong(v, v+1, std::memory_order_relaxed);
              if (res.first) {
                q.emplace_back(std::forward<T>(data));
                spinlock.clear(std::memory_order_release);
                break;
              } // failed reservation.  wait
              spinlock.clear(std::memory_order_release);
            } // else over capacity.  wait.

            _mm_pause();
          } while (true);

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

          size_t v = this->size.load(std::memory_order_relaxed);   // dont' fetch_sub incase we are at 0.

          if ((v & Base::MAX_SIZE) == 0UL) return res;  // nothing to pop

          while (spinlock.test_and_set(std::memory_order_acq_rel));   // establish critical section

          res.first = this->size.compare_exchange_strong(v, v-1, std::memory_order_relaxed);

          if (res.first) {
            if (!q.empty()) {
              res.second = std::move(q.front());  // convert to movable reference and move-assign.
              q.pop_front();
            } else { // size is not zero but q is...
              spinlock.clear(std::memory_order_release);
              throw std::logic_error("size is non-zero, but queue is empty.");
            }
          } // else reservation to pop failed.
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

          size_t v = 0UL;

          do {
            v = this->size.load(std::memory_order_relaxed);

            if (v == Base::DISABLED) {  // disabled and empty. return.
              res.first = false;
              break;
            } else if ((v & Base::MAX_SIZE) > 0) {  // has entries
              while(spinlock.test_and_set(std::memory_order_acq_rel));  // lock access to q, and make size change and q access bundled.
              if (!q.empty()) {
                res.first = this->size.compare_exchange_strong(v, v-1, std::memory_order_relaxed);
                if (res.first)
                {
                  res.second = std::move(q.front());  // convert to movable reference and move-assign.
                  q.pop_front();
                  spinlock.clear(std::memory_order_release);
                  break;
                } // unable to reserve one to pop.

              } // size is not zero but q is.  okay because other threads may get there first.
              spinlock.clear(std::memory_order_release);
            }  // else empty, not disabled.  wait
          
            _mm_pause();
          } while (true);

          return res;
        }



    };

    template<typename T>
    using SpinLockQueue = ThreadSafeQueue<T, bliss::concurrent::LockType::SPINLOCK>;



  } /* namespace concurrent */
} /* namespace bliss */

#endif /* SPINLOCK_QUEUE_HPP_ */
