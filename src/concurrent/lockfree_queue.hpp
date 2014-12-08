/**
 * @file		lockfree_queue.hpp
 * @ingroup bliss::concurrent
 * @author	tpan
 * @brief   Thread Safe Queue Implementation
 * @details this header file contains the templated implementation of a thread safe queue.  this class is used for MPI buffer management.
 *
 *    // can't use boost's lockfree queue:  expect T to have copy constuctor, tribial assignment operator, and trivial destuctor.
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef THREADSAFE_QUEUE_HPP_
#define THREADSAFE_QUEUE_HPP_

#include <thread>
#include <mutex>
#include <limits>
#include <atomic>
#include <stdexcept>
#include <xmmintrin.h>
#include <concurrentqueue/concurrentqueue.h>


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
    class ThreadSafeQueue
    {
      protected:

        /// mutex for locking access to the queue
        mutable std::mutex mutex;

        /// underlying queue that is not thread safe.
        boost::lockfree::queue<T> q;

        /// capacity of the queue.  if set to std::numeric_limits<size_t>::max() indicates unlimited size queue
        size_t capacity;

        /// atomic boolean variable to indicate whether a calling thread can push into this queue.  use when suspending or terminating a queue.
        std::atomic<bool> pushEnabled;

      private:
        /**
         * private move constructor that requires a lock during the call, so that the source of the move is locked.   content of other is moved back.
         * @param other   the soruce ThreadSafeQueue object from which data will be moved.
         * @param l       a lock that uses the mutex of the source ThreadSafeQueue.
         */
        ThreadSafeQueue(ThreadSafeQueue<T>&& other, const std::lock_guard<std::mutex>& l) :
            capacity(other.capacity) {
          other.capacity = 0;
          pushEnabled.exchange(other.pushEnabled.exchange(false));

          q.unsynchronized_push(other.q.)

        };


        /**
         * copy constructor, DISABLED
         * @param other   the source ThreadSafeQueue from which to copy.
         */
        explicit ThreadSafeQueue(const ThreadSafeQueue<T>& other) = delete;

        /**
         * copy assignment operator
         * @param other   the source ThreadSafeQueue from which to copy.
         * @return
         */
        ThreadSafeQueue<T>& operator=(const ThreadSafeQueue<T>& other) = delete;

      public:
        /// maximum possible size of a thread safe queue.  initialized to maximum size_t value.
        static constexpr size_t MAX_SIZE  = std::numeric_limits<size_t>::max();

        /**
         * normal constructor allowing the caller to specify an optional capacity parameter.
         * @param _capacity   The maximum capacity for the thread safe queue.
         */
        explicit ThreadSafeQueue(const size_t &_capacity = MAX_SIZE) :
               capacity(_capacity), qsize(0), pushEnabled(true)
        {
          if (capacity == 0)
            throw std::invalid_argument("ThreadSafeQueue constructor parameter capacity is given as 0");
        };

        /**
         * move constructor.  mutex locks the src ThreadSafeQueue first before delegating to the private constructor.
         * @param other   the source ThreadSafeQueue from which to move.
         */
        explicit ThreadSafeQueue(ThreadSafeQueue<T>&& other) :
            ThreadSafeQueue<T>(std::move(other), std::lock_guard<std::mutex>(other.mutex)) {};

        /**
         * move assignment operator.  locks both the src and destination ThreadSafeQueue before performing the move.
         * @param other   the source ThreadSafeQueue from which to move.
         * @return
         */
        ThreadSafeQueue<T>& operator=(ThreadSafeQueue<T>&& other) {
          std::unique_lock<std::mutex> mylock(mutex, std::defer_lock), otherlock(other.mutex, std::defer_lock);
          std::lock(mylock, otherlock);
          q = std::move(other.q);
          capacity = other.capacity; other.capacity = 0;
          qsize.exchange(other.qsize.exchange(0));             // relaxed, since we have a lock.
          pushEnabled.exchange(other.pushEnabled.exchange(false));
           return *this;
        }

        /**
         * get the capacity of the thread safe queue
         * @return    capacity of the queue
         */
        const size_t& getCapacity() const {
          return capacity;
        }

        /**
         * check if the thread safe queue is full.
         * @return    boolean - whether the queue is full.
         */
        bool isFull() const {
        	if (capacity >= MAX_SIZE) return false;
        
          return qsize.load(std::memory_order_seq_cst) >= capacity;
        }

        /**
         * check if the thread safe queue is empty
         * @return    boolean - whether the queue is empty.
         */
        bool isEmpty() const
        {
          return qsize.load(std::memory_order_seq_cst) == 0;
        }

        /**
         * get the current size of the queue
         * @return    the current size of the queue
         */
        size_t getSize() const
        {
          return qsize.load(std::memory_order_seq_cst);
        }

        /**
         * clears the queue of all contents (discarding).
         */
        void clear()
        {
          while(spinlock.test_and_set());
          q.clear();
          qsize.store(0, std::memory_order_seq_cst);        // have lock.  relaxed.
          spinlock.clear();

        }

        /**
         * set the queue to accept new elements
         */
        void enablePush() {
//          while(spinlock.test_and_set());  simple atomic upate to pushEnabled.  no need for lock
          pushEnabled.store(true, std::memory_order_seq_cst);
          // not notifying canPopCV, used in waitAndPop, as that function exits the loop ONLY when queue is not empty.
        }

        /**
         * set the queue to disallow insertion of new elements.
         */
        void disablePush() {
//          while(spinlock.test_and_set());   simple atomic update to pushEnabled.  no need for lock
          pushEnabled.store(false, std::memory_order_seq_cst);
          //while(spinlock.test_and_set());
          //full_cv = false;   // notify, so waitAndPush can check if it can push now.  (only happens with a full buffer, allows waitToPush to return)
                                    // this allows waitAndPush to fail if push is disabled before the queue becomes not full.
          //empty_cv = false;   // notify, so waitAndPop can exit because there is no more data coming in. (only happens with an empty buffer, allows waitToPop to return)
          //spinlock.clear();
        }

        /**
         * check if the thread safe queue is accepting new elements. (full or not)
         * @return    boolean - queue insertion allowed or not.
         */
        bool canPush() {
          return pushEnabled.load(std::memory_order_seq_cst);
        }

        /**
         * check if the thread safe queue can produce an element now or in near future.  ThreadSafeQueue can pop only if it has elements in the base queue. or
         * if additional items can be pushed in.
         * @return    boolean - queue pop is allowed or not.
         */
        bool canPop() {
          return canPush() || !isEmpty();
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
        bool tryPush (T const& data) {
          while(spinlock.test_and_set());

          if (!canPush() || isFull()) {
            spinlock.clear();
            return false;
          }

          //while(spinlock.test_and_set());
          qsize.fetch_add(1, std::memory_order_seq_cst);
          q.push_back(data);   // insert using predefined copy version of dequeue's push function

          spinlock.clear();
          return true;
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
        std::pair<bool,T> tryPush (T && data) {
          while(spinlock.test_and_set());

          if (!canPush() || isFull()) {
            spinlock.clear();
            return std::move(std::make_pair(false, std::move(data)));;
          }

//          while(spinlock.test_and_set());
          qsize.fetch_add(1, std::memory_order_seq_cst);
          q.push_back(std::move(data));    // insert using predefined move version of deque's push function
          spinlock.clear();

          return std::move(std::make_pair(true, std::move(T())));
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
        bool waitAndPush (T const& data) {

          if (!canPush()) return false;   // if finished, then no more insertion.  return.

          while(spinlock.test_and_set());
          while (!canPush() || isFull()) {
            spinlock.clear();
            _mm_pause();
            // full q.  wait for someone to signal (not full && canPush, or !canPush).

            // to get here, have to have one of these conditions changed:  pushEnabled, !full
            if (!canPush()) {
              return false;  // if finished, then no more insertion.  return.
            }

            while(spinlock.test_and_set());
          }
          qsize.fetch_add(1, std::memory_order_seq_cst);
          q.push_back(data);   // insert using predefined copy version of deque's push function

          spinlock.clear();

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
        std::pair<bool,T> waitAndPush (T && data) {
          if (!canPush())
        	  return std::move(std::make_pair(false, std::move(data)));  // if finished, then no more insertion.  return.

          while(spinlock.test_and_set());
          while (!canPush() || isFull()) {
            // full q.  wait for someone to signal  (not full && canPush, or !canPush).
            spinlock.clear();
            _mm_pause();

            // to get here, have to have one of these conditions changed:  pushEnabled, !full
            if (!canPush()) {
              return std::move(std::make_pair(false, std::move(data)));  // if finished, then no more insertion.  return.
            }
            while(spinlock.test_and_set());
          }

          qsize.fetch_add(1, std::memory_order_seq_cst);
          q.push_back(std::move(data));    // insert using predefined move version of deque's push function

          spinlock.clear();
          return std::move(std::make_pair(true, std::move(T())));
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
        bool tryPushFront (T const& data) {
          while(spinlock.test_and_set());

          if (!canPush() || isFull()) {
            spinlock.clear();
            return false;
          }

          //while(spinlock.test_and_set());
          qsize.fetch_add(1, std::memory_order_seq_cst);
          q.push_front(data);   // insert using predefined copy version of dequeue's push function
          spinlock.clear();

          return true;
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
        std::pair<bool, T> tryPushFront (T && data) {
          while(spinlock.test_and_set());

          if (!canPush() || isFull()) {
            spinlock.clear();
            return std::move(std::make_pair(false, std::move(data)));
          }

//          while(spinlock.test_and_set());
          qsize.fetch_add(1, std::memory_order_seq_cst);
          q.push_front(std::move(data));    // insert using predefined move version of deque's push function
          spinlock.clear();

          return std::move(std::make_pair(true, std::move(T())));
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
        bool waitAndPushFront (T const& data) {

          if (!canPush()) return false;   // if finished, then no more insertion.  return.

          while(spinlock.test_and_set());
          while (!canPush() || isFull()) {
            // full q.  wait for someone to signal (not full && canPush, or !canPush).
            spinlock.clear();
            _mm_pause();

            // to get here, have to have one of these conditions changed:  pushEnabled, !full
            if (!canPush()) {
              return false;  // if finished, then no more insertion.  return.
            }
            while(spinlock.test_and_set());

          }
          qsize.fetch_add(1, std::memory_order_seq_cst);
          q.push_front(data);   // insert using predefined copy version of deque's push function

          spinlock.clear();
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
        std::pair<bool, T> waitAndPushFront (T && data) {
          if (!canPush()) return std::move(std::make_pair(false, std::move(data)));  // if finished, then no more insertion.  return.

          while(spinlock.test_and_set());
          while (!canPush() || isFull()) {
            // full q.  wait for someone to signal  (not full && canPush, or !canPush).
            spinlock.clear();
            _mm_pause();

            // to get here, have to have one of these conditions changed:  pushEnabled, !full
            if (!canPush()) {
              return std::move(std::make_pair(false, std::move(data)));  // if finished, then no more insertion.  return.
            }
            while(spinlock.test_and_set());

          }

          qsize.fetch_add(1, std::memory_order_seq_cst);
          q.push_front(std::move(data));    // insert using predefined move version of deque's push function

          spinlock.clear();

          return std::move(std::make_pair(true, std::move(T())));
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
        std::pair<bool, T> tryPop() {
          std::pair<bool, T> output;
          output.first = false;

          while(spinlock.test_and_set());
          if (isEmpty()) {
            spinlock.clear();
            return output;
          }

          //while(spinlock.test_and_set());
          output.second = std::move(q.front());  // convert to movable reference and move-assign.
          q.pop_front();
          qsize.fetch_sub(1, std::memory_order_seq_cst);
          spinlock.clear();

          output.first = true;

          return output;

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
        std::pair<bool, T> waitAndPop() {
          std::pair<bool, T> output;
          output.first = false;

          // if !canPush and queue is empty, then return false.
          if (!canPop()) return output;

          // else either canPush (queue empty or not), or !canPush and queue is not empty.

          while(spinlock.test_and_set());
          while (isEmpty()) {
            // empty q.  wait for someone to signal (when !isEmpty, or !canPop)
            spinlock.clear();
            _mm_pause();

            // if !canPush and queue is empty, then return false.
            if (!canPop()) {
              return output;
            }
            while(spinlock.test_and_set());
          }

          qsize.fetch_sub(1, std::memory_order_seq_cst);
          output.second = std::move(q.front());  // convert to movable reference and move-assign.
          q.pop_front();

          spinlock.clear();
          output.first = true;


          return output;
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
        std::pair<bool, T> tryPopBack() {
          std::pair<bool, T> output;
          output.first = false;

          while(spinlock.test_and_set());
          if (isEmpty()) {
            spinlock.clear();
            return output;
          }

          //while(spinlock.test_and_set());
          output.second = std::move(q.back());  // convert to movable reference and move-assign.
          q.pop_back();
          qsize.fetch_sub(1, std::memory_order_seq_cst);
          spinlock.clear();

          output.first = true;

          return output;

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
        std::pair<bool, T> waitAndPopBack() {
          std::pair<bool, T> output;
          output.first = false;

          // if !canPush and queue is empty, then return false.
          if (!canPop()) return output;

          // else either canPush (queue empty or not), or !canPush and queue is not empty.

          while(spinlock.test_and_set());

          while (isEmpty()) {
            // empty q.  wait for someone to signal (when !isEmpty, or !canPop)
            spinlock.clear();
            _mm_pause();

            // if !canPush and queue is empty, then return false.
            if (!canPop()) {
              return output;
            }
            while(spinlock.test_and_set());
          }

          qsize.fetch_sub(1, std::memory_order_seq_cst);
          output.second = std::move(q.back());  // convert to movable reference and move-assign.
          q.pop_back();

          spinlock.clear();
          output.first = true;


          return output;
        }


    };

    /**
     * static templated MAX_SIZE definition.
     */
    template<typename T> constexpr size_t ThreadSafeQueue<T>::MAX_SIZE;



  } /* namespace concurrent */
} /* namespace bliss */

#endif /* THREADSAFE_QUEUE_HPP_ */
