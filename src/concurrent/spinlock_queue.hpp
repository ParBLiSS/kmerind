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
#ifndef THREADSAFE_QUEUE_HPP_
#define THREADSAFE_QUEUE_HPP_

#include <cassert>
#include <thread>
#include <mutex>
//#include <deque>
#include <list>
#include <limits>
#include <atomic>
#include <stdexcept>
#include <xmmintrin.h>

#include "concurrent/concurrent.hpp"
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
    class ThreadSafeQueue
    {
      protected:

    	/// spinlock for locking during enqueue, dequeue.
        std::atomic_flag spinlock = ATOMIC_FLAG_INIT;

        /// internal queue, protected by spinlock.
        //std::deque<T> q;
        std::list<T> q;

        /// capacity of the queue.  if set to std::numeric_limits<size_t>::max() indicates unlimited size queue
        mutable int64_t capacity;

        /// the current size of the underlying queue.  using atomic data type avoids having to lock the queue prior to checking its size.
        volatile std::atomic<int64_t> size;

      private:

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
        static constexpr int64_t MAX_SIZE  = std::numeric_limits<int64_t>::max();
        static constexpr int64_t DISABLE_FLAG = std::numeric_limits<int64_t>::lowest();

        /**
         * normal constructor allowing the caller to specify an optional capacity parameter.
         * @param _capacity   The maximum capacity for the thread safe queue.
         */
        explicit ThreadSafeQueue(const size_t &_capacity = static_cast<size_t>(MAX_SIZE)) :
               capacity(static_cast<int64_t>(_capacity)), size(0)
        {
           assert(_capacity <= static_cast<size_t>(MAX_SIZE));

           if (capacity == 0)
            throw std::invalid_argument("ThreadSafeQueue constructor parameter capacity is given as 0");
        };

        /**
         * move constructor.  mutex locks the src ThreadSafeQueue first before delegating to the private constructor.
         * @param other   the source ThreadSafeQueue from which to move.
         */
        explicit ThreadSafeQueue(ThreadSafeQueue<T>&& other) {
        	while (other.spinlock.test_and_set(bliss::concurrent::MO_ACQ_REL));

    			std::swap(q, other.q);
		    	capacity = other.capacity; other.capacity = 0;
			    // ordering may seem strange, but see http://en.cppreference.com/w/cpp/atomic/memory_order re. operator ordering within same thread.
			    size.exchange(other.size.exchange(DISABLE_FLAG, bliss::concurrent::MO_ACQUIRE), bliss::concurrent::MO_RELEASE);
    			other.spinlock.clear(bliss::concurrent::MO_RELEASE);
        };

        /**
         * move assignment operator.  locks both the src and destination ThreadSafeQueue before performing the move.
         * @param other   the source ThreadSafeQueue from which to move.
         * @return
         */
        ThreadSafeQueue<T>& operator=(ThreadSafeQueue<T>&& other) {

        	while (spinlock.test_and_set(bliss::concurrent::MO_ACQ_REL));
        	while (other.spinlock.test_and_set(bliss::concurrent::MO_ACQ_REL));

        	std::swap(q, other.q);
        	capacity = other.capacity; other.capacity = 0;
        	size.exchange(other.size.exchange(DISABLE_FLAG, bliss::concurrent::MO_ACQUIRE), bliss::concurrent::MO_RELEASE);

        	other.spinlock.clear(bliss::concurrent::MO_RELEASE);
        	spinlock.clear(bliss::concurrent::MO_RELEASE);

           return *this;
        }

        /**
         * get the capacity of the thread safe queue
         * @return    capacity of the queue
         */
        inline size_t getCapacity() const {
          return capacity;
        }

        /// check if queue is fixed size or growable.
        inline bool isFixedSize() const {
        	return getCapacity() < MAX_SIZE;
        }

        /**
         * check if the thread safe queue is full.
         * @return    boolean - whether the queue is full.
         */
        inline bool isFull() const {
          return (isFixedSize()) && (getSize() >= getCapacity());
        }

        /**
         * check if the thread safe queue is empty
         * @return    boolean - whether the queue is empty.
         */
        inline bool isEmpty() const
        {
          return getSize() == 0;
        }

        /**
         * get the current size of the queue
         * @return    the current size of the queue
         */
        inline size_t getSize() const
        {
          return static_cast<size_t>(size.load(bliss::concurrent::MO_RELAXED) & MAX_SIZE);   // size is atomic, so don't need strong memory ordering itself.
        }

        /**
         * clears the queue of all contents (discarding).
         */
        void clear()
        {
          while(spinlock.test_and_set(bliss::concurrent::MO_ACQ_REL));
          q.clear();
          // clear size
          size.fetch_and(DISABLE_FLAG, bliss::concurrent::MO_ACQ_REL);  // keep the push bit, and set size to 0

          spinlock.clear(bliss::concurrent::MO_RELEASE);
        }

        /**
         * set the queue to accept new elements
         */
        inline void enablePush() {
          // before enablePush, queue is probably empty or no one is writing to it.  so prior side effects don't need to be visible right away.
          // size itself is atomic
          while(spinlock.test_and_set(bliss::concurrent::MO_ACQ_REL));
          size.fetch_and(MAX_SIZE, bliss::concurrent::MO_ACQ_REL);  // clear the push bit, and leave size as is
          spinlock.clear(bliss::concurrent::MO_RELEASE);
        }

        /**
         * set the queue to disallow insertion of new elements.
         */
        inline void disablePush() {
          // before disable push, should make sure that all writes are visible to all other threads, so release here.
          while(spinlock.test_and_set(bliss::concurrent::MO_ACQ_REL));
          size.fetch_or(DISABLE_FLAG, bliss::concurrent::MO_ACQ_REL);   // set the push bit, and leave size as is.
          spinlock.clear(bliss::concurrent::MO_RELEASE);
        }

        /**
         * check if the thread safe queue can accept new elements. (full or not)
         * @return    boolean - queue insertion allowed or not.
         */
        inline bool canPush() {
          // pushing thread does not immediately care what other threads did before this
          // size is >= 0 (not disabled), and less than capacity.
          // note that if we reinterpret_cast size to size_t, then disabled (< 0) will have MSB set to 1, so > max(int64_t).
          return size.load(bliss::concurrent::MO_RELAXED) >= 0;   // int highest bit set means negative, and means cannot push
        }

        /**
         * @brief    check if the thread safe queue can produce an element now or in near future.
         * @details  ThreadSafeQueue can pop only if it has elements in the base queue. or
         *            if additional items can be pushed in.
         * @return    boolean - queue pop is allowed or not.
         */
        inline bool canPop() {
          // canPop == first bit is 0, OR has some elements (not 0 for remaining bits).  so basically, not 1000000000...
          // popping thread should have visibility of all changes to the queue
          return size.load(bliss::concurrent::MO_RELAXED) != DISABLE_FLAG;
        }

      protected:
        /**
         * check if the thread safe queue can accept new elements. (full or not)
         * @return    boolean - queue insertion allowed or not.
         */
        inline bool canPushAndHasRoom() {
          // if we reinterpret this number as a uint64_t, then we only need to check less than capacity, since now highest bits are all way higher.
          int64_t v = size.load(bliss::concurrent::MO_RELAXED);  // since uint64_t&, need a variable
          return reinterpret_cast<uint64_t&>(v) < capacity;
        }

      public:

        /**
         * @brief Non-blocking - pushes an element by constant reference (copy).
         * Returns true only if push was successful.
         * If queue is full, or if queue is not accepting new element inserts, return false.
         *    Does not modify the element if cannot push onto queue.
         *
         * @param data    data element to be pushed onto the thread safe queue
         * @return        whether push was successful.
         */
        bool tryPush (T const& data) {
          bool res = false;
          while(spinlock.test_and_set(bliss::concurrent::MO_ACQ_REL));    // critical block instructions confined here via acquire/release

          int64_t v = size.load(bliss::concurrent::MO_ACQUIRE);        // load and RMW are ordered via acq/release at the fence.
          if (reinterpret_cast<uint64_t&>(v) < static_cast<uint64_t>(capacity)) {
            // following 2 instructions cannot be reordered to above this line,
            // and the if check cannot be reordered to below this line
            std::atomic_thread_fence(bliss::concurrent::MO_ACQ_REL);            

      			q.push_back(data);  // insert using predefined copy version of dequeue's push function
            res = true;
            size.fetch_add(1, bliss::concurrent::MO_ACQ_REL);
          } // else at capacity.
          spinlock.clear(bliss::concurrent::MO_RELEASE);
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
        std::pair<bool,T> tryPush (T && data) {
            std::pair<bool, T> res;
            res.first = false;

            while(spinlock.test_and_set(bliss::concurrent::MO_ACQ_REL));   // see const ref version for memory order info.

            int64_t v = size.load(bliss::concurrent::MO_ACQUIRE);
            if (reinterpret_cast<uint64_t&>(v) < static_cast<uint64_t>(capacity)) {
              std::atomic_thread_fence(bliss::concurrent::MO_ACQ_REL);

          	  q.emplace_back(std::forward<T>(data));  // insert using predefined copy version of dequeue's push function
              res.first = true;
              size.fetch_add(1, bliss::concurrent::MO_ACQ_REL);
            } else { // else at capacity.
              res.second = std::move(data);
            }

            spinlock.clear(bliss::concurrent::MO_RELEASE);
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
        bool waitAndPush (T const& data) {

          volatile int64_t v = 0L;
          bool res = false;

          do {  // loop forever, unless queue is disabled, or insert succeeded.
        	  while(spinlock.test_and_set(bliss::concurrent::MO_ACQ_REL));

            v = size.load(bliss::concurrent::MO_ACQUIRE);

            if (v < 0L) { // disabled so return false
             	spinlock.clear(bliss::concurrent::MO_RELEASE);
              break;
            }

            if (v >= capacity) {   // over capacity.  wait.
             	spinlock.clear(bliss::concurrent::MO_RELEASE);
            	_mm_pause();
              continue;
            }

            // else between 0 and capacity - 1.  enqueue.

          	std::atomic_thread_fence(bliss::concurrent::MO_ACQ_REL);    // separate the conditional from the actual enqueue

           	q.push_back(data);  // successfully enqueued.
            res = true;
            size.fetch_add(1, bliss::concurrent::MO_ACQ_REL);
  
          	spinlock.clear(bliss::concurrent::MO_RELEASE);
            break;                        
            
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
        std::pair<bool,T> waitAndPush (T && data) {

          std::pair<bool, T> res;
          res.first = false;

          volatile int64_t v = 0L;

          do {  // loop forever, unless queue is disabled, or insert succeeded.
          	while(spinlock.test_and_set(bliss::concurrent::MO_ACQ_REL));

            v = size.load(bliss::concurrent::MO_ACQUIRE);

            if (v < 0L) { // disabled so return false
              res.second = std::move(data);
            	spinlock.clear(bliss::concurrent::MO_RELEASE);
              break;
            }

            if (v >= capacity) {   // over capacity.  wait.
             	spinlock.clear(bliss::concurrent::MO_RELEASE);
            	_mm_pause();
              continue;
            }

            // else between 0 and capacity - 1.  enqueue.

          	std::atomic_thread_fence(bliss::concurrent::MO_ACQ_REL);

           	q.emplace_back(std::forward<T>(data));  // successfully enqueued.
            res.first = true;
            size.fetch_add(1, bliss::concurrent::MO_ACQ_REL);
            	
          	spinlock.clear(bliss::concurrent::MO_RELEASE);
            break;

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
        std::pair<bool, T> tryPop() {
          std::pair<bool, T> output;
          output.first = false;

          if (isEmpty()) return output;

          while (spinlock.test_and_set(bliss::concurrent::MO_ACQ_REL));   // establish critical section

          if (!isEmpty()) {
        	  std::atomic_thread_fence(bliss::concurrent::MO_ACQ_REL);    // fence separate conditional from dequeue.
        	  output.second = std::move(q.front());  // convert to movable reference and move-assign.
            std::atomic_thread_fence(bliss::concurrent::MO_ACQ_REL);    // fence enforce q.front and q.pop_front
        	  q.pop_front();
            output.first = true;
            size.fetch_sub(1, bliss::concurrent::MO_ACQ_REL);
          }

          spinlock.clear(bliss::concurrent::MO_RELEASE);

          return output;

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
        std::pair<bool, T> waitAndPop() {
          std::pair<bool, T> output;
          output.first = false;


          volatile int64_t v = 0L;
          bool breaking = false;

          do {

            while (spinlock.test_and_set(bliss::concurrent::MO_ACQ_REL));

//            printf("checking..."); fflush(stdout);
            v = size.load(bliss::concurrent::MO_ACQUIRE);


//            if (v < 0) printf("-"); fflush(stdout);
            if (v == DISABLE_FLAG) {  // disabled and empty. return.
//              printf("-"); fflush(stdout);
              breaking = true;
            } else if (v == 0L) {  // empty but not disabled. so wait
//              printf("%ld", v); fflush(stdout);
            } else {

//            printf("popping..."); fflush(stdout);

              //  else has some entries
            
              std::atomic_thread_fence(bliss::concurrent::MO_ACQ_REL);   // fence to separate conditional and dequeue
              output.second = std::move(q.front());  // convert to movable reference and move-assign.
              std::atomic_thread_fence(bliss::concurrent::MO_ACQ_REL);
              q.pop_front();
              
              output.first = true;
              size.fetch_sub(1, bliss::concurrent::MO_ACQ_REL);
              breaking = true;
            }
            spinlock.clear(bliss::concurrent::MO_RELEASE);

            if (breaking) break;
          
            _mm_pause();
          } while (!output.first);

//          printf("done..."); fflush(stdout);

          return output;
        }



    };

    /**
     * static templated MAX_SIZE definition.
     */
    template<typename T> constexpr int64_t ThreadSafeQueue<T>::MAX_SIZE;



  } /* namespace concurrent */
} /* namespace bliss */

#endif /* THREADSAFE_QUEUE_HPP_ */
