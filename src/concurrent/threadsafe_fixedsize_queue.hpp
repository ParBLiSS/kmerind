/**
 * @file		threadsafe_fixed_size_queue.hpp
 * @ingroup
 * @author	tpan
 * @brief
 * @details adapted from http://www.justsoftwaresolutions.co.uk/threading/implementing-a-thread-safe-queue-using-condition-variables.html
 *          with modifications
 *            1. use c++11 std::thread constructs
 *            2. incorporating move semantics.
 *            3. support fixed size queue - blocking if queue is full.
 *
 *            note that this is NOT truly concurrent.  the class serializes parallel access.  move semantic should minimize the copy operations needed.
 *
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef THREADSAFE_FIXEDSIZE_QUEUE_HPP_
#define THREADSAFE_FIXEDSIZE_QUEUE_HPP_

#include <thread>
#include <mutex>
#include <condition_variable>
#include <queue>

namespace bliss
{
  namespace concurrent
  {

    /**
     * @class			bliss::concurrent::ThreadSafeFixedSizeQueue
     * @brief
     * @details
     *
     */
    template <typename T>
    class ThreadSafeFixedSizeQueue
    {
      private:
        std::queue<T> q;
        mutable std::mutex mutex;
        std::condition_variable empty_cv;
        std::condition_variable full_cv;
        const size_t maxSize;

      public:

        ThreadSafeFixedSizeQueue(size_t max_size) : maxSize(max_size) {
        }

        const size_t& getMaxSize() const {
          return maxSize;
        }

        bool full() const
        {
          std::unique_lock<std::mutex> lock(mutex);
          return q.size() == maxSize;
        }


        bool tryPush (T const& data) {
          std::unique_lock<std::mutex> lock(mutex);
          if (q.size() == maxSize) {
            // empty q.  wait for someone to signal.
            return false;
          }

          q.push(data);   // move using predefined copy refernece version of push
          lock.unlock();
          empty_cv.notify_one();
          return true;
        }

        bool tryPush (T && data) {
          std::unique_lock<std::mutex> lock(mutex);
          if (q.size() == maxSize) {
            // empty q.  wait for someone to signal.
            return false;
          }

          q.push(data);    // move using predefined move reference version of push
          lock.unlock();
          empty_cv.notify_one();
          return true;
        }

        void waitAndPush (T const& data) {
          std::unique_lock<std::mutex> lock(mutex);
          while (q.size() == maxSize) {
            // full q.  wait for someone to signal.
            full_cv.wait(lock);
          }

          q.push(data);   // move using predefined copy refernece version of push
          lock.unlock();
          empty_cv.notify_one();
        }

        void waitAndPush (T && data) {
          std::unique_lock<std::mutex> lock(mutex);
          while (q.size() == maxSize) {
            // full q.  wait for someone to signal.
            full_cv.wait(lock);
//            printf("full wait over\n");
            //printf(".");
          }
          //printf("\n");

          q.push(data);    // move using predefined move reference version of push
          lock.unlock();
          empty_cv.notify_one();
        }

        size_t size() const
        {
          std::unique_lock<std::mutex> lock(mutex);
          return q.size();
        }

        bool empty() const
        {
          std::unique_lock<std::mutex> lock(mutex);
          return q.empty();
        }

        bool tryPop(T& output) {
          std::unique_lock<std::mutex> lock(mutex);
          if (q.empty()) {
            // empty q, nothing.
            return false;
          }

          output = std::move(q.front());  // convert to movable reference and move-assign.
          q.pop();
          lock.unlock();
          full_cv.notify_one();
          return true;
        }

        void waitAndPop(T& output) {
          std::unique_lock<std::mutex> lock(mutex);
          while (q.empty()) {
            // empty q.  wait for someone to signal.
//            printf("empty wait\n");
            empty_cv.wait(lock);
//            printf("empty wait over\n");
          }
          // when cond_var is notified, then lock will be acquired and q will be examined.

          output = std::move(q.front());  // convert to movable reference and move-assign.
          q.pop();
          lock.unlock();
          full_cv.notify_one();
        }

    };

  } /* namespace concurrent */
} /* namespace bliss */

#endif /* THREADSAFE_FIXEDSIZE_QUEUE_HPP_ */
