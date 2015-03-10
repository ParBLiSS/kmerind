/**
 * @file		mutexlock_queue.hpp
 * @ingroup wip
 * @author	Tony Pan <tpan7@gatech.edu>
 * @brief   Thread Safe Queue Implementation
 * @details this header file contains the templated implementation of a thread safe queue.  this class is used for MPI buffer management.
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef MUTEXLOCK_QUEUE_HPP_
#define MUTEXLOCK_QUEUE_HPP_

#include <deque>

#include "concurrent/threadsafe_queue.hpp"

// TODO: every operation that's modifying the queue is using unique lock. can this be made better with just careful atomic operations
// 		 e.g. with memory fence?
//       see http://en.wikipedia.org/wiki/Non-blocking_algorithm#cite_note-lf-queue-13, (CAS)
//       and http://www.cs.technion.ac.il/~mad/publications/ppopp2013-x86queues.pdf,	(F&A)
//

namespace bliss
{
  namespace concurrent
  {
    /**
     * @class     bliss::concurrent::ThreadSafeQueue
     * @brief     a multi-producer, multi-consumer thread safe queue with an optional capacity
     * @details   adapted from http://www.justsoftwaresolutions.co.uk/threading/implementing-a-thread-safe-queue-using-condition-variables.html
     *            with SIGNIFICANT modifications
     *            1. use c++11 std::thread and atomic constructs
     *            2. incorporating move semantics.
     *            3. support a capacity limit.
     *            4. supports multiple producer, multiple consumer.
     *            5. allows blocking of the enqueue function, (useful when draining queue or when finishing the use of the queue).
     *
     *          note that this is NOT truly concurrent.  the class serializes parallel access.  move semantic minimizes the copy operations needed.
     *
     *          DUE TO LOCKS, this is NOT fast.  but may be fast enough for MPI buffer management.
     */
    template <typename T>
    class ThreadSafeQueue<T, bliss::concurrent::LockType::MUTEX> :
      public detail::ThreadSafeQueueBase<ThreadSafeQueue<T, bliss::concurrent::LockType::MUTEX> >
    {
    	//static_assert(false, "Mutex Locking Thread Safe Queue encounters deadlock from time to time.  also slow.  Please do not use.");

      protected:
        using Derived = ThreadSafeQueue<T, bliss::concurrent::LockType::MUTEX>;
        using Base = detail::ThreadSafeQueueBase< Derived >;

        /// mutex for locking access to the queue during construction/assignment
        mutable std::mutex mutex;


        /// condition variable for event notification.  specifically, for unblocking the waitAndPop calls (queue transitions from empty to not-empty)
        std::condition_variable canPopCV;

        /// condition variable for event notification.  specifically, for unblocking the waitAndPush calls (queue transitions from full to not-full)
        std::condition_variable canPushCV;

        /// underlying queue that is not thread safe.
        std::deque<T> q;

      private:
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


        /**
         * private move constructor that requires a lock during the call, so that the source of the move is locked.   content of other is moved back.
         * @param other   the soruce ThreadSafeQueue object from which data will be moved.
         * @param l       a lock that uses the mutex of the source ThreadSafeQueue.
         */
        ThreadSafeQueue(Derived&& other, const std::lock_guard<std::mutex>& l) {

        	std::swap(q, other.q);

          this->capacity = other.capacity; other.capacity = 0;
          VAR(this->size) = VAR(other.size);   VAR(other.size) = Base::DISABLED;

        };


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
        explicit ThreadSafeQueue(Derived&& other) :
            Derived(std::move(other), std::lock_guard<std::mutex>(other.mutex)) {};

        /**
         * move assignment operator.  locks both the src and destination ThreadSafeQueue before performing the move.
         * @param other   the source ThreadSafeQueue from which to move.
         * @return
         */
        Derived& operator=(Derived&& other) {
          std::unique_lock<std::mutex> mylock(this->mutex, std::defer_lock),
              otherlock(other.mutex, std::defer_lock);
          std::lock(mylock, otherlock);

          std::swap(q, other.q);
          this->capacity = other.capacity; other.capacity = 0;
          VAR(this->size) = VAR(other.size);  VAR(other.size) = Base::DISABLED;

          return *this;
        }

        /**
         * clears the queue of all contents (discarding).
         */
        void clearImpl()
        {
          std::unique_lock<std::mutex> lock(this->mutex);
          q.clear();
          // clear size
          VAR(this->size) &= Base::DISABLED;  // keep the push bit, and set size to 0
          lock.unlock();

          CV_NOTIFY_ALL(canPushCV);
        }

        /**
         * set the queue to accept new elements
         */
        inline void enablePushImpl() { // overrides the base class non-virtual version
          // before enablePush, queue is probably empty or no one is writing to it.  so prior side effects don't need to be visible right away.
          // size itself is atomic
          std::unique_lock<std::mutex> lock(this->mutex);
          VAR(this->size) &= Base::MAX_SIZE;  // clear the push bit, and leave size as is
          lock.unlock();
          CV_NOTIFY_ALL(canPushCV);   // notify, so waitAndPush can check if it can push now.
          // not notifying canPopCV, used in waitAndPop, as that function exits the loop ONLY when queue is not empty.
        }

        /**
         * set the queue to disallow insertion of new elements.
         */
        inline void disablePushImpl() {  // overrides the base class non-virtual version
            // before disable push, should make sure that all writes are visible to all other threads, so release here.
          std::unique_lock<std::mutex> lock(this->mutex);
          VAR(this->size) |= Base::DISABLED;   // set the push bit, and leave size as is.
          lock.unlock();
          CV_NOTIFY_ALL(canPushCV);   // notify, so waitAndPush can check if it can push now.  (only happens with a full buffer, allows waitToPush to return)
                                    // this allows waitAndPush to fail if push is disabled before the queue becomes not full.
          CV_NOTIFY_ALL(canPopCV);   // notify, so waitAndPop can exit because there is no more data coming in. (only happens with an empty buffer, allows waitToPop to return)
        }

        /**
         * Non-blocking - pushes an element by constant reference (copy).
         * Returns true only if push was successful.
         * If queue is full, or if queue is not accepting new element inserts, return false.
         *    Does not modify the element if cannot push onto queue.
         *
         * @param data    data element to be pushed onto the thread safe queue
         * @return        whether push was successful.
         */
        bool tryPushImpl (T const& data) {
          std::unique_lock<std::mutex> lock(this->mutex);
          int64_t s = VAR(this->size);
          if (reinterpret_cast<uint64_t&>(s) < static_cast<uint64_t>(this->capacity)) {
        	  q.push_back(data);  // insert using predefined copy version of dequeue's push function
        	  ++VAR(this->size);
        	  lock.unlock();
              CV_NOTIFY_ONE(canPopCV);
              return true;
          }  // else at capacity.


          lock.unlock();

          return false;
        }

        /**
         * Non-blocking - pushes an element by constant reference (move).
         * Returns true only if push was successful.
         * If queue is full, or if queue is not accepting new element inserts, return false.
         *    Does not modify the element if cannot push onto queue.
         *
         * @param data    data element to be pushed onto the thread safe queue
         * @return        whether push was successful.
         */
        std::pair<bool, T> tryPushImpl (T && data) {
          std::pair<bool, T> res(false, T());
          std::unique_lock<std::mutex> lock(this->mutex);

          int64_t s = VAR(this->size);

          if (reinterpret_cast<uint64_t&>(s) < static_cast<uint64_t>(this->capacity)) {
        	  q.push_back(std::forward<T>(data));  // insert using predefined copy version of dequeue's push function
        	  ++VAR(this->size);
        	  lock.unlock();
              CV_NOTIFY_ONE(canPopCV);

            res.first = true;
            return std::move(res);
          }  // else at capacity.

          lock.unlock();
          res.second = std::move(data);
          return std::move(res);;

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


          if (!this->canPush()) return false;   // if finished, then no more insertion.  return.
          std::unique_lock<std::mutex> lock(this->mutex);
          int64_t s = VAR(this->size);

          while (reinterpret_cast<uint64_t&>(s) >= reinterpret_cast<uint64_t&>(this->capacity)) {
            // full q.  wait for someone to signal (not full && canPush, or !canPush).
            CV_WAIT(canPushCV, lock);
            // to get here, have to have one of these conditions changed:  pushEnabled, !full
            if (VAR(this->size) < 0) {  // blocked.
              lock.unlock();
              return false;  // if finished, then no more insertion.  return.
            }
            s = VAR(this->size);
          }
          ++VAR(this->size);
          q.push_back(data);   // insert using predefined copy version of deque's push function

          lock.unlock();
          CV_NOTIFY_ONE(canPopCV);
          return true;

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

            if (!this->canPush()) {
              res.second = std::move(data);
          	  return std::move(res);  // if finished, then no more insertion.  return.
            }

            std::unique_lock<std::mutex> lock(this->mutex);
            int64_t s = VAR(this->size);
            while (reinterpret_cast<uint64_t&>(s) >= reinterpret_cast<uint64_t&>(this->capacity)) {
              // full q.  wait for someone to signal (not full && canPush, or !canPush).
              CV_WAIT(canPushCV, lock);

              // to get here, have to have one of these conditions changed:  pushEnabled, !full
              if (VAR(this->size) < 0) {
                lock.unlock();
                res.second = std::move(data);
                return std::move(res);  // if finished, then no more insertion.  return.
              }
              s = VAR(this->size);
            }
            ++VAR(this->size);
            q.push_back(std::forward<T>(data));   // insert using predefined copy version of deque's push function

            lock.unlock();
            CV_NOTIFY_ONE(canPopCV);
            res.first = true;
            return std::move(res);

        }





        /**
         * Non-blocking - remove the first element in the queue and return it to the calling thread.
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

          std::unique_lock<std::mutex> lock(this->mutex);
          if (!this->isEmpty()) {
        	  res.second = std::move(q.front());  // convert to movable reference and move-assign.
        	  q.pop_front();
        	    --VAR(this->size);
              res.first = true;
              CV_NOTIFY_ONE(canPushCV);
          }

          lock.unlock();
          return res;

        }

        /**
         * Semi-blocking - remove the first element in the queue and return it to the calling thread.
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


          // if !canPush and queue is empty, then return false.
          if (!this->canPop()) return res;
          // else either canPush (queue empty or not), or !canPush and queue is not empty.

          // while the queue is not disabled nor empty
          std::unique_lock<std::mutex> lock(this->mutex);

          while (this->isEmpty()) {
            // empty q.  wait for someone to signal (when !isEmpty, or !canPop)
            CV_WAIT(canPopCV, lock);

            // if !canPush and queue is empty, then return false.
            if (!this->canPop()) {
              lock.unlock();
              return res;
            }
          }

          --VAR(this->size);

          res.second = std::move(q.front());  // convert to movable reference and move-assign.
          q.pop_front();
          res.first = true;

          lock.unlock();

          CV_NOTIFY_ONE(canPushCV);

          return res;
        }



    };


    template<typename T>
    using MutexLockQueue = bliss::concurrent::ThreadSafeQueue<T, bliss::concurrent::LockType::MUTEX>;

  } /* namespace concurrent */
} /* namespace bliss */

#endif /* MUTEXLOCK_QUEUE_HPP_ */
