/**
 * @file		concurrent_queue.hpp
 * @ingroup
 * @author	tpan
 * @brief
 * @details adapted from http://www.justsoftwaresolutions.co.uk/threading/implementing-a-thread-safe-queue-using-condition-variables.html
 *          with modifications
 *            1. use c++11 std::thread constructs
 *            2. incorporating move semantics.
 *
 *          note that this is NOT truly concurrent.  the class serializes parallel access.  move semantic should minimize the copy operations needed.
 *
 *          DUE TO LOCKS, this is NOT fast.  but may be fast enough for MPI buffer management.
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


namespace bliss
{
  namespace concurrent
  {
    /**
     * @class     bliss::concurrent::ThreadSafeQueue
     * @brief
     * @details
     *
     */
    template <typename T>
    class ThreadSafeQueue
    {
      private:
        mutable std::mutex mutex;
        std::condition_variable canPopCV;
        std::condition_variable canPushCV;

        std::deque<T> q;
        size_t capacity;
        std::atomic<size_t> qsize;
        std::atomic<bool> pushEnabled;
        std::atomic<bool> popEnabled;

        ThreadSafeQueue(ThreadSafeQueue<T>&& other, const std::lock_guard<std::mutex>&) : q(std::move(other.q)),
            capacity(other.capacity), qsize(0), pushEnabled(false), popEnabled(false) {
          other.capacity = 0;
          qsize.exchange(other.qsize.load(std::memory_order_relaxed), std::memory_order_relaxed);              // relaxed, since we have a lock.
          pushEnabled.exchange(other.pushEnabled.load(std::memory_order_relaxed), std::memory_order_relaxed);
          popEnabled.exchange(other.popEnabled.load(std::memory_order_relaxed), std::memory_order_relaxed);
        };


      public:
        static constexpr size_t MAX_SIZE  = std::numeric_limits<size_t>::max();
        explicit ThreadSafeQueue(const size_t &_capacity = MAX_SIZE) : capacity(_capacity), qsize(0), pushEnabled(true), popEnabled(true) {};

        explicit ThreadSafeQueue(const ThreadSafeQueue<T>& other) = delete;
        explicit ThreadSafeQueue(ThreadSafeQueue<T>&& other) : ThreadSafeQueue<T>(std::move(other), std::lock_guard<std::mutex>(other.mutex)) {};
        ThreadSafeQueue<T>& operator=(ThreadSafeQueue<T>&& other) {
          std::unique_lock<std::mutex> mylock(mutex, std::defer_lock), otherlock(other.mutex, std::defer_lock);
          std::lock(mylock, otherlock);
          q = std::move(other.q);
          capacity = other.capacity; other.capacity = 0;
          qsize.exchange(other.qsize.load(std::memory_order_relaxed), std::memory_order_relaxed);             // relaxed, since we have a lock.
          pushEnabled.exchange(other.pushEnabled.load(std::memory_order_relaxed), std::memory_order_relaxed);
          popEnabled.exchange(other.popEnabled.load(std::memory_order_relaxed), std::memory_order_relaxed);
          return *this;
        }
        ThreadSafeQueue<T>& operator=(const ThreadSafeQueue<T>& other) = delete;


        const size_t& getCapacity() const {
          return capacity;
        }

        bool full() const {
        	if (capacity == MAX_SIZE) return false;
        
//          std::unique_lock<std::mutex> lock(mutex);
          return qsize.load(std::memory_order_consume) >= capacity;
        }
        bool empty() const
        {
          //std::unique_lock<std::mutex> lock(mutex);
          return qsize.load(std::memory_order_consume) == 0;
        }

        size_t size() const
        {
          return qsize.load(std::memory_order_consume);
        }

        void clear()
        {
          std::unique_lock<std::mutex> lock(mutex);
          q.clear();
          qsize.store(0, std::memory_order_relaxed);        // have lock.  relaxed.
        }

        void enablePush() {
          pushEnabled.store(true, std::memory_order_release);
//          printf("push enabled\n");
          canPushCV.notify_all();   // notify, so waitAndPush can check if it can push now.
        }
        void disablePush() {
          pushEnabled.store(false, std::memory_order_release);
//          printf("push disabled\n");
          canPopCV.notify_all();   // notify, so waitAndPop can exit because there is no more data coming in.
        }
        bool canPush() {
          return pushEnabled.load(std::memory_order_acquire);
        }


//        void enablePop() {
//          popEnabled.store(true, std::memory_order_release);
//          printf("pop enabled\n");
//          canPopCV.notify_all();
//        }
//        void disablePop() {
//          popEnabled.store(false, std::memory_order_release);
//          printf("pop disabled\n");
//          canPopCV.notify_all();
//        }
        bool canPop() {
          //return popEnabled.load(std::memory_order_acquire);
          return canPush() || !empty();
        }




        bool tryPush (T const& data) {
          std::unique_lock<std::mutex> lock(mutex);

          if (!canPush() || full()) return false;

          qsize.fetch_add(1, std::memory_order_acq_rel);
          q.push_back(data);   // move using predefined copy refernece version of push

          lock.unlock();
          canPopCV.notify_one();
          return true;
        }

        bool tryPush (T && data) {
          std::unique_lock<std::mutex> lock(mutex);
          if (!canPush() || full()) return false;

          qsize.fetch_add(1, std::memory_order_acq_rel);
          q.push_back(std::move(data));    // move using predefined move reference version of push

          lock.unlock();
          canPopCV.notify_one();
          return true;
        }

        bool waitAndPush (T const& data) {

          if (!canPush()) return false;   // if finished, then no more pushing.  return.

          std::unique_lock<std::mutex> lock(mutex);
          while (!canPush() || full()) {
            // full q.  wait for someone to signal.
            canPushCV.wait(lock);

            // to get here, have to have one of these conditions changed:  !finished, pushEnabled, !full
            if (!canPush()) return false;  // if finished, then no more pushing.  return.
          }
          qsize.fetch_add(1, std::memory_order_acq_rel);
          q.push_back(data);   // move using predefined copy refernece version of push

          lock.unlock();
          canPopCV.notify_one();
          return true;
        }

        bool waitAndPush (T && data) {
          if (!canPush()) return false;  // if finished, then no more pushing.  return.

          std::unique_lock<std::mutex> lock(mutex);


          while (!canPush() || full()) {
            // full q.  wait for someone to signal.
            canPushCV.wait(lock);

            // to get here, have to have one of these conditions changed:  !finished, pushEnabled, !full
            if (!canPush()) return false;  // if finished, then no more pushing.  return.
          }

          qsize.fetch_add(1, std::memory_order_acq_rel);
          q.push_back(std::move(data));    // move using predefined move reference version of push

          lock.unlock();
          canPopCV.notify_one();
          return true;
        }



        std::pair<bool, T> tryPop() {
          std::pair<bool, T> output;
          output.first = false;

          std::unique_lock<std::mutex> lock(mutex);
          if (empty()) return output;

          qsize.fetch_sub(1, std::memory_order_acq_rel);
          output.second = std::move(q.front());  // convert to movable reference and move-assign.
          q.pop_front();

          lock.unlock();
          output.first = true;

          canPushCV.notify_one();
          return output;

        }

        std::pair<bool, T> waitAndPop() {
          std::pair<bool, T> output;
          output.first = false;

          // if finished and queue is empty, then return false.
          if (!canPop()) return output;  // if finished, then no more pushing.  return.

          // else either not finished (queue empty or not), or finished and queue is not empty.

          std::unique_lock<std::mutex> lock(mutex);
          while (empty()) {
            // empty q.  wait for someone to signal.
            canPopCV.wait(lock);

            if (!canPop()) return output;  // if finished, then no more pushing.  return.
          }
          // when cond_var is notified, then lock will be acquired and q will be examined.

          qsize.fetch_sub(1, std::memory_order_acq_rel);
          output.second = std::move(q.front());  // convert to movable reference and move-assign.
          q.pop_front();

          lock.unlock();
          output.first = true;

          canPushCV.notify_one();

          return output;
        }

    };

    template<typename T> constexpr size_t ThreadSafeQueue<T>::MAX_SIZE;



  } /* namespace concurrent */
} /* namespace bliss */

#endif /* THREADSAFE_QUEUE_HPP_ */
