/**
 * @file		test_threadsafe_queue.cpp
 * @ingroup
 * @author	Tony Pan <tpan7@gatech.edu>
 * @brief
 * @details
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */


#include <iostream>
#include <cstdio>
#include <unistd.h>
#include <chrono>

#if defined(BLISS_NONE)
#include <queue>
#elif defined( BLISS_MUTEX)
#include "wip/mutexlock_queue.hpp"
#elif defined( BLISS_SPINLOCK )
#include "concurrent/spinlock_queue.hpp"
#else   //if defined( BLISS_LOCKFREE )
#include "concurrent/lockfree_queue.hpp"
#endif


#include "omp.h"



#if defined(BLISS_NONE)

template<typename T>
void timeSTDQueueSerial(const std::string &message, std::queue<T>&& q, const int entries) {
  printf("=== %s, capacity unlimited, single thread enqueue dequeue test: ", message.c_str());

  int count = 0;

  T output;
  T result = 0;
  for (int i = 0; i < entries; ++i) {

    if (i % 3 == 0)
      q.push(T(i));
    else {
      if (!q.empty()) {
        output = q.front();
        q.pop();
        result ^= output;
        ++count;
      }
    }
  }
  if (count != entries / 3) printf("FAIL: entries %d, numPopped: %d, expected: %d", entries, count, entries/3);
  else printf("PASS");

  printf(" result = %d\n", result);

};

template<typename T>
void timeSTDQueueThreaded(const std::string &message, std::queue<T>&& q, const int entries, const int nProducer, const int nConsumer) {
  printf("=== %s, capacity unlimited, multithread tests:", message.c_str());

  int count = 0, count2 = 0;
  volatile T result = 0;

  omp_lock_t lock;
  omp_init_lock(&lock);

  volatile bool done = false;

#pragma omp parallel sections num_threads(2) shared(q, done, result)
  {
#pragma omp section
    {
      volatile int localCount = entries-2;
      // insert
#pragma omp parallel num_threads(nProducer) shared(q, localCount) reduction(+:count)
      {
        int id;

        while (true) {
#pragma omp atomic capture
          {
            id = localCount;
            --localCount;
          }

          if (id < 0) break;

          omp_set_lock(&lock);
          q.push(std::move(T(id)));
          omp_unset_lock(&lock);
          ++count;
        }
      }

#pragma omp critical
      done = true;
    }


#pragma omp section
    {
      // insert
#pragma omp parallel num_threads(nConsumer) shared(q, done) reduction(+:count2) reduction(^:result)
      {
        while (!done || !q.empty()) {

          omp_set_lock(&lock);
          if (!q.empty()) {
            result ^= q.front();
            q.pop();
            ++count2;
          }
          omp_unset_lock(&lock);

        }
      }


    }
  }

  if ((count != entries-1) || (count2 != entries-1)) printf("FAIL: entries %d, numPushed = %d, numPopped= %d, expected= %d", entries-1, count, count2, entries-1);
  else printf("PASS");

  printf(" result = %d\n", result);

  omp_destroy_lock(&lock);
};


#else

using namespace bliss::concurrent;


template<typename T>
void timeTSQueueSingleThread(const std::string &message, bliss::concurrent::ThreadSafeQueue<T>&& q, const int entries) {
  printf("=== %s, capacity %lu, single thread enqueue dequeue test: ", message.c_str(), q.getCapacity());

  int count = 0;
  T result = 0;

  for (int i = 0; i < entries; ++i) {

    if (i % 3 == 0)
      q.tryPush(std::move(T(i)));
    else
    {
      auto r2 = q.tryPop();
      if (r2.first) {
        result ^= r2.second;
        ++count;
      }
    }
  }
  if (count != entries / 3) printf("FAIL: entries %d, numPopped: %d, expected: %d", entries, count, entries/3);
  else printf("PASS");

  printf(" result = %d\n", result);


};


template<typename T>
void timeTSQueueThreaded(const std::string &message, bliss::concurrent::ThreadSafeQueue<T>&& q, const int entries, const int nProducer, const int nConsumer) {
  printf("=== %s, capacity %lu, multithread tests:", message.c_str(), q.getCapacity());

  int count = 0, count2 = 0;
  volatile T result = 0;


#pragma omp parallel sections num_threads(2) shared(q, count, count2, result)
  {
#pragma omp section
    {
      // push
#pragma omp parallel for default(none) num_threads(nProducer) shared(q) reduction(+:count)
      for (int i = 0; i < entries-1; ++i) {
        q.tryPush(std::move(T(i)));
        ++count;
      }

      q.disablePush();


    }

#pragma omp section
    {
      T results[nConsumer];
      // pop
#pragma omp parallel num_threads(nConsumer) default(none) shared(q, results) reduction(+: count2)
      {
        results[omp_get_thread_num()] = 0;

        while (q.canPop()) {
          auto r2 = q.tryPop();

          if (r2.first) {
            results[omp_get_thread_num()] ^= r2.second;
            ++count2;
          }
        }
      }


    #pragma omp parallel for default(none) num_threads(nConsumer) shared(results) reduction(^:result)
      for (int i = 0; i < nConsumer; ++i) {
          result = results[omp_get_thread_num()];
      }

    }
  }

  if ((count != entries-1) || (count2 != entries-1)) printf("FAIL: entries %d, numPushed = %d, numPopped= %d, expected= %d", entries-1, count, count2, entries-1);
  else printf("PASS");

  printf(" result = %d\n", result);

};


#endif

int main(int argc, char** argv) {

  std::chrono::high_resolution_clock::time_point t1, t2;
  std::chrono::duration<double> time_span;

  omp_set_nested(1);
  omp_set_dynamic(0);

  int count = 10000000;

#if defined(BLISS_NONE)
//  ========= STD VERSION


    t1 = std::chrono::high_resolution_clock::now();
    timeSTDQueueSerial("TIME STD 1 thread", std::move(std::queue<int>()), count);
    t2 = std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    printf("  time: entries %d nthread %d time %f\n", count, 1, time_span.count());

    t1 = std::chrono::high_resolution_clock::now();
    timeSTDQueueThreaded("TIME STD 4 thread", std::move(std::queue<int>()), count, 4, 1);
    t2 = std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    printf("  time: entries %d nthread %d,%d time %f\n", count, 4, 1, time_span.count());

    t1 = std::chrono::high_resolution_clock::now();
    timeSTDQueueThreaded("TIME STD 4 thread", std::move(std::queue<int>()), count, 1, 4);
    t2 = std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    printf("  time: entries %d nthread %d,%d time %f\n", count, 1, 4, time_span.count());

    t1 = std::chrono::high_resolution_clock::now();
    timeSTDQueueThreaded("TIME STD 4 thread", std::move(std::queue<int>()), count, 4, 4);
    t2 = std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    printf("  time: entries %d nthread %d,%d time %f\n", count, 4, 4, time_span.count());

#else

    typedef bliss::concurrent::ThreadSafeQueue<int> QueueType;

  t1 = std::chrono::high_resolution_clock::now();
  timeTSQueueSingleThread("TIME TSQ 1 thread ", std::move(QueueType()), count);
  t2 = std::chrono::high_resolution_clock::now();
  time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
  printf("  time: entries %d nthread %d time %f\n", count, 1, time_span.count());


  t1 = std::chrono::high_resolution_clock::now();
  timeTSQueueThreaded("TIME TSQ 4 threads ", std::move(QueueType()), count, 4, 1);
  t2 = std::chrono::high_resolution_clock::now();
  time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
  printf("  time: entries %d nthread %d,%d time %f\n", count, 4, 1, time_span.count());


  t1 = std::chrono::high_resolution_clock::now();
  timeTSQueueThreaded("TIME TSQ 4 threads ", std::move(QueueType()), count, 1, 4);
  t2 = std::chrono::high_resolution_clock::now();
  time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
  printf("  time: entries %d nthread %d,%d time %f\n", count, 1, 4, time_span.count());

  t1 = std::chrono::high_resolution_clock::now();
  timeTSQueueThreaded("TIME TSQ 4 threads ", std::move(QueueType()), count, 4, 4);
  t2 = std::chrono::high_resolution_clock::now();
  time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
  printf("  time: entries %d nthread %d,%d time %f\n", count, 4, 4, time_span.count());

#endif
};
