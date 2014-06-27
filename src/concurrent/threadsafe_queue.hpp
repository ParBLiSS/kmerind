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
        std::deque<T> q;
        mutable std::mutex mutex;
        std::condition_variable empty_cv;
        std::condition_variable full_cv;
        size_t capacity;

        ThreadSafeQueue(ThreadSafeQueue<T>&& other, const std::lock_guard<std::mutex>&) : q(std::move(other.q)),
            capacity(other.capacity) {
          other.capacity = 0;
        };


      public:
        static constexpr size_t MAX_SIZE = std::numeric_limits<size_t>::max();
        explicit ThreadSafeQueue(const size_t &_capacity = MAX_SIZE) : capacity(_capacity) {};

        explicit ThreadSafeQueue(const ThreadSafeQueue<T>& other) = delete;
        explicit ThreadSafeQueue(ThreadSafeQueue<T>&& other) : ThreadSafeQueue<T>(std::move(other), std::lock_guard<std::mutex>(other.mutex)) {};
        ThreadSafeQueue<T>& operator=(ThreadSafeQueue<T>&& other) {
          capacity = other.capacity; other.capacity = 0;
          std::unique_lock<std::mutex> mylock(mutex, std::defer_lock), otherlock(other.mutex, std::defer_lock);
          std::lock(mylock, otherlock);
          q = std::move(other.q);
          return *this;
        }
        ThreadSafeQueue<T>& operator=(const ThreadSafeQueue<T>& other) = delete;


        const size_t& getCapacity() const {
          return capacity;
        }

        bool full() const {
        	if (capacity == MAX_SIZE) return false;
        
          std::unique_lock<std::mutex> lock(mutex);
          return q.size() >= capacity;
        }


        bool tryPush (T const& data) {
          std::unique_lock<std::mutex> lock(mutex);
          if (q.size() >= capacity) {
            // empty q.  wait for someone to signal.
            return false;
          }

          q.push_back(data);   // move using predefined copy refernece version of push
          lock.unlock();
          empty_cv.notify_one();
          return true;
        }

        bool tryPush (T && data) {
          std::unique_lock<std::mutex> lock(mutex);
          if (q.size() >= capacity) {
            // empty q.  wait for someone to signal.
            return false;
          }

          q.push_back(std::move(data));    // move using predefined move reference version of push
          lock.unlock();
          empty_cv.notify_one();
          return true;
        }

        void waitAndPush (T const& data) {
          std::unique_lock<std::mutex> lock(mutex);
          while (q.size() >= capacity) {
            // full q.  wait for someone to signal.
            full_cv.wait(lock);
          }

          q.push_back(data);   // move using predefined copy refernece version of push
          lock.unlock();
          empty_cv.notify_one();
        }

        void waitAndPush (T && data) {
          std::unique_lock<std::mutex> lock(mutex);
          while (q.size() >= capacity) {
            // full q.  wait for someone to signal.
            full_cv.wait(lock);
          }

          q.push_back(std::move(data));    // move using predefined move reference version of push
          lock.unlock();
          empty_cv.notify_one();
        }

        size_t size() const
        {
          std::unique_lock<std::mutex> lock(mutex);
          return q.size();
        }

        void clear()
        {
          std::unique_lock<std::mutex> lock(mutex);
          q.clear();
        }

        bool empty() const
        {
          std::unique_lock<std::mutex> lock(mutex);
          return q.empty();
        }

        std::pair<bool, T> tryPop() {
          std::pair<bool, T> output;
          std::unique_lock<std::mutex> lock(mutex);
          if (q.empty()) {
            lock.unlock();
            // empty q, nothing.
            output.first = false;
            return output;
          } else {
            output.second = std::move(q.front());  // convert to movable reference and move-assign.
            q.pop_front();
            lock.unlock();
            output.first = true;
            full_cv.notify_one();
            return output;
          }
        }

        T waitAndPop() {
          T output;
          std::unique_lock<std::mutex> lock(mutex);
          while (q.empty()) {
            // empty q.  wait for someone to signal.
//            printf("empty wait\n");
            empty_cv.wait(lock);
//            printf("empty wait over\n");
          }
          // when cond_var is notified, then lock will be acquired and q will be examined.

          output = std::move(q.front());  // convert to movable reference and move-assign.
          q.pop_front();
          lock.unlock();
          full_cv.notify_one();

          return output;
        }

    };





  } /* namespace concurrent */
} /* namespace bliss */

#endif /* THREADSAFE_QUEUE_HPP_ */
