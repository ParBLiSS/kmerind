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
#ifndef THREADSAFE_QUEUE_HPP_
#define THREADSAFE_QUEUE_HPP_

#include <thread>
#include <mutex>
#include <condition_variable>
#include <deque>
#include <limits>
#include <atomic>
#include <stdexcept>
#include <cassert>

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
    	//static_assert(false, "Mutex Locking Thread Safe Queue encounters deadlock from time to time.  also slow.  Please do not use.");

      protected:

        /// mutex for locking access to the queue
        mutable std::mutex mutex;

        /// condition variable for event notification.  specifically, for unblocking the waitAndPop calls (queue transitions from empty to not-empty)
        std::condition_variable canPopCV;

        /// condition variable for event notification.  specifically, for unblocking the waitAndPush calls (queue transitions from full to not-full)
        std::condition_variable canPushCV;


        /// underlying queue that is not thread safe.
        std::deque<T> q;

        /// capacity of the queue.  if set to std::numeric_limits<size_t>::max() indicates unlimited size queue
        mutable int64_t capacity;

        /// the current size of the underlying queue.  using atomic data type avoids having to lock the queue prior to checking its size.
        std::atomic<int64_t> size;

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


        /**
         * private move constructor that requires a lock during the call, so that the source of the move is locked.   content of other is moved back.
         * @param other   the soruce ThreadSafeQueue object from which data will be moved.
         * @param l       a lock that uses the mutex of the source ThreadSafeQueue.
         */
        ThreadSafeQueue(ThreadSafeQueue<T>&& other, const std::lock_guard<std::mutex>& l) {

        	std::atomic_thread_fence(std::memory_order_acquire);

			std::swap(q, other.q);
			capacity = other.capacity;  other.capacity = 0;
			size.exchange(other.size.exchange(std::numeric_limits<int64_t>::lowest(), std::memory_order_relaxed), std::memory_order_relaxed);
			// relaxed, since we have a lock.
        	std::atomic_thread_fence(std::memory_order_release);

        };


      public:
        /// maximum possible size of a thread safe queue.  initialized to maximum size_t value.
        static constexpr size_t MAX_SIZE  = std::numeric_limits<int64_t>::max();

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
        explicit ThreadSafeQueue(ThreadSafeQueue<T>&& other) : ThreadSafeQueue<T>(std::move(other), std::lock_guard<std::mutex>(other.mutex)) {};

        /**
         * move assignment operator.  locks both the src and destination ThreadSafeQueue before performing the move.
         * @param other   the source ThreadSafeQueue from which to move.
         * @return
         */
        ThreadSafeQueue<T>& operator=(ThreadSafeQueue<T>&& other) {
          std::unique_lock<std::mutex> mylock(mutex, std::defer_lock), otherlock(other.mutex, std::defer_lock);
          std::lock(mylock, otherlock);

      	std::atomic_thread_fence(std::memory_order_acquire);

      	std::swap(q, other.q);
      	capacity = other.capacity; other.capacity = 0;
      	size.exchange(other.size.exchange(std::numeric_limits<int64_t>::lowest(), std::memory_order_relaxed), std::memory_order_relaxed);

      	std::atomic_thread_fence(std::memory_order_release);

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
          return static_cast<size_t>(size.load(std::memory_order_relaxed) & MAX_SIZE);   // size is atomic, so don't need strong memory ordering itself.
        }


        /**
         * clears the queue of all contents (discarding).
         */
        void clear()
        {
          std::unique_lock<std::mutex> lock(mutex);
          q.clear();
          // clear size
          size.fetch_and(std::numeric_limits<int64_t>::lowest(), std::memory_order_release);  // keep the push bit, and set size to 0
          lock.unlock();

          canPushCV.notify_all();
        }

        /**
         * set the queue to accept new elements
         */
        inline void enablePush() {
			// before enablePush, queue is probably empty or no one is writing to it.  so prior side effects don't need to be visible right away.
			// size itself is atomic
			size.fetch_and(MAX_SIZE, std::memory_order_relaxed);  // clear the push bit, and leave size as is
          canPushCV.notify_all();   // notify, so waitAndPush can check if it can push now.
          // not notifying canPopCV, used in waitAndPop, as that function exits the loop ONLY when queue is not empty.
        }

        /**
         * set the queue to disallow insertion of new elements.
         */
        inline void disablePush() {
            // before disable push, should make sure that all writes are visible to all other threads, so release here.
            size.fetch_or(std::numeric_limits<int64_t>::lowest(), std::memory_order_relaxed);   // set the push bit, and leave size as is.
          canPushCV.notify_all();   // notify, so waitAndPush can check if it can push now.  (only happens with a full buffer, allows waitToPush to return)
                                    // this allows waitAndPush to fail if push is disabled before the queue becomes not full.
          canPopCV.notify_all();   // notify, so waitAndPop can exit because there is no more data coming in. (only happens with an empty buffer, allows waitToPop to return)
        }

        /**
         * check if the thread safe queue can accept new elements. (full or not)
         * @return    boolean - queue insertion allowed or not.
         */
        inline bool canPush() {
          // pushing thread does not immediately care what other threads did before this
          // size is >= 0 (not disabled), and less than capacity.
          // note that if we reinterpret_cast size to size_t, then disabled (< 0) will have MSB set to 1, so > max(int64_t).
          return size.load(std::memory_order_relaxed) >= 0;   // int highest bit set means negative, and means cannot push
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
          return size.load(std::memory_order_relaxed) != std::numeric_limits<int64_t>::lowest();
        }

      protected:
        /**
         * check if the thread safe queue can accept new elements. (full or not)
         * @return    boolean - queue insertion allowed or not.
         */
        inline bool canPushAndHasRoom() {
          // if we reinterpret this number as a uint64_t, then we only need to check less than capacity, since now highest bits are all way higher.
          int64_t v = size.load(std::memory_order_relaxed);
          return reinterpret_cast<uint64_t&>(v) < capacity;
        }

      public:

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
          std::unique_lock<std::mutex> lock(mutex);

          int64_t v = size.fetch_add(1, std::memory_order_acquire);
          if (reinterpret_cast<uint64_t&>(v) < static_cast<uint64_t>(capacity)) {
        	  q.push_back(data);  // insert using predefined copy version of dequeue's push function
        	  std::atomic_thread_fence(std::memory_order_release);
        	  lock.unlock();
              canPopCV.notify_one();
              return true;
          }  // else at capacity.
          size.fetch_sub(1, std::memory_order_relaxed); // failed enqueue, decrement size.

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
        std::pair<bool,T> tryPush (T && data) {
          std::unique_lock<std::mutex> lock(mutex);


          int64_t v = size.fetch_add(1, std::memory_order_acquire);
          if (reinterpret_cast<uint64_t&>(v) < static_cast<uint64_t>(capacity)) {
        	  q.push_back(std::forward<T>(data));  // insert using predefined copy version of dequeue's push function
        	  std::atomic_thread_fence(std::memory_order_release);
        	  lock.unlock();
              canPopCV.notify_one();
            return std::move(std::make_pair(true, std::forward<T>(data)));
          }  // else at capacity.
          size.fetch_sub(1, std::memory_order_relaxed); // failed enqueue, decrement size.

          lock.unlock();
          return std::move(std::make_pair(false, std::forward<T>(data)));;

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
          std::unique_lock<std::mutex> lock(mutex);
          int64_t v = size.load(std::memory_order_relaxed);
          while (reinterpret_cast<uint64_t&>(v) >= reinterpret_cast<uint64_t&>(capacity)) {
            // full q.  wait for someone to signal (not full && canPush, or !canPush).
            canPushCV.wait(lock);
            v = size.load(std::memory_order_relaxed);
            // to get here, have to have one of these conditions changed:  pushEnabled, !full
            if (v < 0) {  // blocked.
              lock.unlock();
              return false;  // if finished, then no more insertion.  return.
            }
          }
          size.fetch_add(1, std::memory_order_acquire);
          q.push_back(data);   // insert using predefined copy version of deque's push function
          std::atomic_thread_fence(std::memory_order_release);

          lock.unlock();
          canPopCV.notify_one();
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
          	  return std::move(std::make_pair(false, std::forward<T>(data)));  // if finished, then no more insertion.  return.

            std::unique_lock<std::mutex> lock(mutex);
            int64_t v = size.load(std::memory_order_relaxed);
            while (reinterpret_cast<uint64_t&>(v) >= reinterpret_cast<uint64_t&>(capacity)) {
              // full q.  wait for someone to signal (not full && canPush, or !canPush).
              canPushCV.wait(lock);
              v = size.load(std::memory_order_relaxed);
              // to get here, have to have one of these conditions changed:  pushEnabled, !full
              if (v < 0) {
                lock.unlock();
                return std::move(std::make_pair(false, std::forward<T>(data)));  // if finished, then no more insertion.  return.
              }
            }
            size.fetch_add(1, std::memory_order_acquire);
            q.push_back(std::forward<T>(data));   // insert using predefined copy version of deque's push function
            std::atomic_thread_fence(std::memory_order_release);

            lock.unlock();
            canPopCV.notify_one();
            return std::move(std::make_pair(true, std::forward<T>(data)));

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

          if (isEmpty()) return output;

          std::unique_lock<std::mutex> lock(mutex);
          if (!isEmpty()) {
        	  std::atomic_thread_fence(std::memory_order_acquire);
        	  output.second = std::move(q.front());  // convert to movable reference and move-assign.
        	  q.pop_front();
              size.fetch_sub(1, std::memory_order_release);
              output.first = true;
              canPushCV.notify_one();
          }

          lock.unlock();
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

          // while the queue is not disabled nor empty
          std::unique_lock<std::mutex> lock(mutex);

          while (isEmpty()) {
            // empty q.  wait for someone to signal (when !isEmpty, or !canPop)
            canPopCV.wait(lock);

            // if !canPush and queue is empty, then return false.
            if (!canPop()) {
              lock.unlock();
              return output;
            }
          }

          size.fetch_sub(1, std::memory_order_acquire);

          output.second = std::move(q.front());  // convert to movable reference and move-assign.
          q.pop_front();
          std::atomic_thread_fence(std::memory_order_release);
          output.first = true;

          lock.unlock();

          canPushCV.notify_one();

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
