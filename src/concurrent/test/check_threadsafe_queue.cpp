/**
 * @file		test_threadsafe_queue.cpp
 * @ingroup
 * @author	tpan
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
#include <queue>

#include "concurrent/threadsafe_queue.hpp"

#include "omp.h"

using namespace bliss::concurrent;





template<typename T>
void testWaitAndPush(bliss::concurrent::ThreadSafeQueue<T> &queue, const int entries, const int nProducer) {
  usleep(1000);
  int count = 0;
#pragma omp parallel for default(none) num_threads(nProducer) shared(queue) reduction(+:count)
  for (int i = 0; i < entries; ++i) {
    if (queue.waitAndPush(T(i)))
//    if (i % 1000 == 0) usleep(5);
      ++count;
    usleep(5);
  }
  queue.disablePush();

  if (count != entries) printf("TSQueue capacity %lu, finished waitAndPush with %d entries.\n", queue.getCapacity(), count);
};

template<typename T>
void testTryPush(bliss::concurrent::ThreadSafeQueue<T> &queue, const int entries, const int nProducer) {
  usleep(1000);
  int count = 0, count2 = 0;
#pragma omp parallel for default(none) num_threads(nProducer) shared(queue) reduction(+: count, count2)
  for (int i = 0; i < entries; ++i) {
    if (queue.tryPush(T(i)))
      ++count;
    else
      ++count2;
//    if (i % 1000 == 0) usleep(20);
    usleep(5);
  }
  queue.disablePush();

  if (count + count2 != entries || count == 0)
    printf("TSQueue capacity %lu, finished tryPush. %d successful, %d failed\n", queue.getCapacity(), count, count2);

};

template<typename T>
void testWaitAndPop1(bliss::concurrent::ThreadSafeQueue<T> &queue) {
  int count = 0, count3 = 0;

  while (queue.canPop()) {
    if (queue.waitAndPop().first)
      ++count;

    //if (count %1000) usleep(20);
  }
  // now empty it
  while (queue.tryPop().first) {
    ++count3;
  }
  if (!queue.isEmpty() || count3 != 0) printf("TSQueue capacity %lu, finished waitAndPop 1 thread. %d successful, %d flushed.  empty? %s\n", queue.getCapacity(), count, count3, (queue.isEmpty() ? "yes" : "no"));
};

template<typename T>
void testTryPop1(bliss::concurrent::ThreadSafeQueue<T> &queue) {
  int count = 0, count2 = 0, count3 = 0;

  while (queue.canPop()) {
    if (queue.tryPop().first)
      ++count;
    else
      ++count2;

    //if (count %1000) usleep(20);
  }
  // now empty it
  while (queue.tryPop().first) {
    ++count3;
  }
  if (!queue.isEmpty() || (count == 0 && count3 == 0)) printf("TSQueue capacity %lu, finished tryPop 1 thread. %d successful, %d failed, %d flushed.  empty? %s\n", queue.getCapacity(), count, count2, count3, (queue.isEmpty() ? "yes" : "no"));

};

template<typename T>
void testTryPop(bliss::concurrent::ThreadSafeQueue<T> &queue, const int nConsumer) {


  int counts[nConsumer];
  int counts2[nConsumer];
  int counts3[nConsumer];

  // while producer is producing, don't know the size.  spawn tasks.
#pragma omp parallel num_threads(nConsumer) default(none) shared(queue, counts, counts2, counts3)
  {
    counts[omp_get_thread_num()] = 0;
    counts2[omp_get_thread_num()] = 0;
    counts3[omp_get_thread_num()] = 0;

#pragma omp single nowait
    {
      while(queue.canPop()) {

        // deferred task to pop.
#pragma omp task default(none) shared(queue, counts, counts2)
        {
          if (queue.tryPop().first)
            ++counts[omp_get_thread_num()];
          else
            ++counts2[omp_get_thread_num()];

          //if ((counts[omp_get_thread_num()] % 1000) == 0) usleep(20);
        }

      };
    }

  }  // implicit barrier.  won't have issue with reporting too early.
  // NO WASTED TIME here
  // tasks create and complete fast, then some wasted time in flush task queue.  data queue may be empty by now.
  // tasks slow but task create fast, then task queue long, take extra time to flush task queue, data queue may be empty here
  // or tasks slow and task create slow, then no wasted task or time in this phase.  data queue not empty
  // or tasks fast and task create slow, then task queue short, and data queue not empty

  // now use parallel for to flush the data queue
  size_t size = queue.getSize();
#pragma omp parallel for num_threads(nConsumer) default(none) shared(queue, size, counts2, counts3)
  for (int i = 0; i < size; ++i) {
    if (queue.tryPop().first)
      ++counts3[omp_get_thread_num()];
    else
      ++counts2[omp_get_thread_num()];
  }


  // summarize
  int count = 0, count2 = 0, count3 = 0;
#pragma omp parallel num_threads(nConsumer) default(none) shared(counts, counts2, counts3) reduction(+:count, count2, count3)
  {
    count = counts[omp_get_thread_num()];
    count2 = counts2[omp_get_thread_num()];
    count3 = counts3[omp_get_thread_num()];
  }
  if (!queue.isEmpty() || (count == 0 && count3 == 0)) printf("TSQueue capacity %lu, finished tryPop. %d successful, %d failed, %d flushed.  empty? %s\n", queue.getCapacity(), count, count2, count3, (queue.isEmpty() ? "yes" : "no"));

};

// no waitAndPop.  if task creation is too fast, then we may have more waitAndPop calls than there are entries in queue, then deadlock.



template<typename T>
void testTSQueue(const std::string &message, bliss::concurrent::ThreadSafeQueue<T>&& queue, const int nProducer, const int nConsumer) {

  printf("=== TEST %s: %d producers, %d consumers, capacity %lu\n", message.c_str(), nProducer, nConsumer, queue.getCapacity());

  int entries = (queue.getCapacity() == bliss::concurrent::ThreadSafeQueue<T>::MAX_SIZE) ? 128 : queue.getCapacity();

  int i = 0;
  int count = 0, count2 = 0;


  printf("  CHECK tryPop on empty\n");
  queue.clear();
#pragma omp parallel for num_threads(nConsumer) private(i) shared(queue, entries) default(none) reduction(+: count, count2)
  for (i = 0; i < entries; ++i) {
    if ( queue.tryPop().first )
      ++count;
    else
      ++count2;
  }
  if (count2 != entries) printf("TSQueue capacity %lu, finished tryPop on empty at iteration %d, success %d, fail %d\n", queue.getCapacity(), i, count, count2);

  printf("  CHECK tryPush too much\n");
  count = 0;
  count2 = 0;
  queue.clear();
#pragma omp parallel for num_threads(nProducer) private(i) shared(queue, entries) default(none) reduction(+: count, count2)
  for (i = 0; i < (entries + 2); ++i) {
    if (queue.tryPush(T(i)))
      ++count;
    else
      ++count2;
  }
  if (count != std::min(queue.getCapacity(), 2UL + entries) || (count + count2 != (entries+2)))  printf("TSQueue capacity %lu, finished tryPush until full at iteration %d, success %d, fail %d. \n", queue.getCapacity(), i, count, count2);

  printf("  CHECK tryPop too much\n");
  count = 0;
  count2 = 0;
#pragma omp parallel for num_threads(nConsumer) private(i) shared(queue, entries) default(none) reduction(+: count, count2)
  for (i = 0; i < (entries + 2); ++i) {
    if ( queue.tryPop().first)
      ++count;
    else
      ++count2;
  }
  if (count != std::min(queue.getCapacity(), 2UL + entries) || (count + count2 != (entries+2))) printf("TSQueue capacity %lu, finished tryPop from full at iteration %d, success %d, fail %d\n",queue.getCapacity(),  i, count, count2);



  printf("  CHECK  waitAndPush then do waitAndPop in parallel\n");
  queue.clear();
  queue.enablePush();
  count = 0;
#pragma omp parallel for num_threads(nProducer) private(i) shared(queue, entries) default(none) reduction(+: count, count2)
  for (i = 0; i < entries; ++i) {
    if (queue.waitAndPush(T(i)))
      ++count;
  }
  if (count != entries) printf("TSQueue capacity %lu, finished waitAndPush with %d entries\n", queue.getCapacity(), count);

  count2 = 0;
#pragma omp parallel for num_threads(nConsumer) private(i) shared(queue, count) default(none) reduction(+: count2)
  for (i = 0; i < count; ++i) {
    if (queue.waitAndPop().first)
      ++count2;
  }
  if (count2 != count) printf("TSQueue capacity %lu, finished waitAndPop with, %d entries.  should be %d\n", queue.getCapacity(), count2, entries);




  printf("  CHECK tryPush, and tryPop.\n");
  queue.clear();
  queue.enablePush();
#pragma omp parallel sections num_threads(2) shared(queue, entries) default(none)
  {
#pragma omp section
    {
      testTryPush(queue, entries, nProducer);
      //printf("tryPush done correctly.  final queue size (with consumer): %lu\n", queue.size());
    }

#pragma omp section
    {
      if (nConsumer > 1)
        testTryPop(queue, nConsumer);
      else
        testTryPop1(queue);
    }
  }


  printf("  CHECK waitAndPush, and tryPop.\n");
  queue.clear();
  queue.enablePush();
#pragma omp parallel sections num_threads(2) shared(queue, entries) default(none)
  {
#pragma omp section
    {
      testWaitAndPush(queue, entries, nProducer);
    }

#pragma omp section
    {
      if (nConsumer > 1)
        testTryPop(queue, nConsumer);
      else
        testTryPop1(queue);
    }
  }

  //TODO: can have !done, waitAndPop, and done=true in other thread, so pop thread never gets to check "done" ->  deadlock.

    // not testing this with more than 1 consumer thread.  there could be a lot more consumer tasks generated than
    //   there are elements because the produce thread is generating slower than consume threads.  the excess threads result in wait forever.
    // this could still happen if we slow down producer and speed up consumers



    printf("  CHECK tryPush, and waitAndPop.\n");
    queue.clear();
    queue.enablePush();
  #pragma omp parallel sections num_threads(2) shared(queue, entries) default(none)
    {
  #pragma omp section
      {
        testTryPush(queue, entries, nProducer);
      }

  #pragma omp section
      {
        testWaitAndPop1(queue);
      }
    }



    printf("  CHECK waitAndPush, and waitAndPop.\n");
    queue.clear();
    queue.enablePush();
  #pragma omp parallel sections num_threads(2) shared(queue, entries) default(none)
    {
  #pragma omp section
      {
        testWaitAndPush(queue, entries, nProducer);
      }

  #pragma omp section
      {
        testWaitAndPop1(queue);
      }
    }
};



template<typename T>
void timeSTDQueueSerial(const std::string &message, std::queue<T>&& q, const int entries) {
  printf("=== %s, capacity unlimited, entries %d, results ", message.c_str(), entries);

  T output;
  T result = 0;
  for (int i = 0; i < entries-1; ++i) {

    if (i % 2 == 0)
      q.push(T(i));
    else {
      output = q.front();
      q.pop();
      result ^= output;
    }
  }
  printf("%d\n", result);

};

template<typename T>
void timeSTDQueueThreaded(const std::string &message, std::queue<T>&& q, const int entries, const int nProducer, const int nConsumer) {
  printf("=== %s, capacity unlimited, entries %d, results ", message.c_str(), entries);

  omp_lock_t lock;

  omp_init_lock(&lock);

#pragma omp parallel sections num_threads(2) shared(q)
  {
#pragma omp section
    {
      int localCount = entries;
      // insert
#pragma omp parallel num_threads(nProducer) shared(q, localCount)
      {
        int count;
        while (true) {
#pragma omp atomic capture
          {
            count = localCount;
            --localCount;
          }

          if (count <= 0) break;

          omp_set_lock(&lock);
          q.push(std::move(T(count)));
          omp_unset_lock(&lock);
        }
      }
    }


#pragma omp section
    {
      int localCount = entries - 1;
      T result;
      // insert
#pragma omp parallel num_threads(nConsumer) shared(q, localCount, result)
      {
        int count;
        T output;
        while (true) {
#pragma omp atomic capture
          {
            count = localCount;
            --localCount;
          }

          if (count <= 0) break;

          omp_set_lock(&lock);
          if (!q.empty()) {
            output = q.front();
            q.pop();
          }
          omp_unset_lock(&lock);

#pragma omp atomic
          result ^= output;
        }
      }
      printf("%d\n", result);


    }
  }

  omp_destroy_lock(&lock);
};




template<typename T>
void timeTSQueueSingleThread(const std::string &message, bliss::concurrent::ThreadSafeQueue<T>&& q, const int entries) {
  printf("=== %s, capacity %lu, ", message.c_str(), q.getCapacity());

  T result = 0;

  for (int i = 0; i < entries - 1; ++i) {

    if (i % 2 == 0)
      q.tryPush(std::move(T(i)));
    else
    {
      auto r2 = q.tryPop();
      if (r2.first) {
        result ^= r2.second;

        //printf("%d %d\n", r2.second, result);
      }
    }
  }
  printf("entries %d, xor result: %d\n", entries, result);
};


template<typename T>
void timeTSQueueThreaded(const std::string &message, bliss::concurrent::ThreadSafeQueue<T>&& q, const int entries, const int nProducer, const int nConsumer) {
  printf("=== %s, capacity %lu, ", message.c_str(), q.getCapacity());

  int count = 0, count2 = 0;

#pragma omp parallel sections num_threads(2) shared(q, count, count2)
  {
#pragma omp section
    {
#pragma omp parallel for default(none) num_threads(nProducer) shared(q) reduction(+:count)
      for (int i = 0; i < entries - 1; ++i) {
        q.tryPush(std::move(T(i)));
        ++count;
      }
      q.disablePush();
    }

#pragma omp section
    {
      T results[nConsumer];
      // insert
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


      T result = 0;
    #pragma omp parallel for default(none) num_threads(nConsumer) shared(results) reduction(^:result)
      for (int i = 0; i < nConsumer; ++i) {
          result = results[omp_get_thread_num()];
      }
      printf("entries %d, pushed %d, popped %d, xor result: %d\n", entries, count, count2, result);

    }
  }



};




int main(int argc, char** argv) {

  std::chrono::high_resolution_clock::time_point t1, t2;
  std::chrono::duration<double> time_span;

  omp_set_nested(1);
  omp_set_dynamic(0);

  typedef bliss::concurrent::ThreadSafeQueue<int> QueueType;


  timeTSQueueSingleThread("TSQ 1 thread growable", std::move(QueueType()),1000);
  timeTSQueueSingleThread("TSQ 1 thread 10 element", std::move(QueueType(10)), 1000);

  testTSQueue("TSQ nthread growable", std::move(QueueType()), 1, 1);
  testTSQueue("TSQ nthread growable", std::move(QueueType()), 2, 1);
  testTSQueue("TSQ nthread growable", std::move(QueueType()), 3, 1);
  testTSQueue("TSQ nthread growable", std::move(QueueType()), 4, 1);
  testTSQueue("TSQ nthread growable", std::move(QueueType()), 1, 2);
  testTSQueue("TSQ nthread growable", std::move(QueueType()), 1, 3);
  testTSQueue("TSQ nthread growable", std::move(QueueType()), 1, 4);
  testTSQueue("TSQ nthread growable", std::move(QueueType()), 2, 2);
  testTSQueue("TSQ nthread growable", std::move(QueueType()), 2, 3);
  testTSQueue("TSQ nthread growable", std::move(QueueType()), 3, 2);

  testTSQueue("TSQ nthread 100 elements", std::move(QueueType(100)), 1, 1);
  testTSQueue("TSQ nthread 100 elements", std::move(QueueType(100)), 2, 1);
  testTSQueue("TSQ nthread 100 elements", std::move(QueueType(100)), 3, 1);
  testTSQueue("TSQ nthread 100 elements", std::move(QueueType(100)), 4, 1);
  testTSQueue("TSQ nthread 100 elements", std::move(QueueType(100)), 1, 2);
  testTSQueue("TSQ nthread 100 elements", std::move(QueueType(100)), 1, 3);
  testTSQueue("TSQ nthread 100 elements", std::move(QueueType(100)), 1, 4);
  testTSQueue("TSQ nthread 100 elements", std::move(QueueType(100)), 2, 2);
  testTSQueue("TSQ nthread 100 elements", std::move(QueueType(100)), 2, 3);
  testTSQueue("TSQ nthread 100 elements", std::move(QueueType(100)), 3, 2);

  int count = 10000000;

  t1 = std::chrono::high_resolution_clock::now();
  timeTSQueueSingleThread("TIME TSQ 1 thread ", std::move(QueueType()), count);
  t2 = std::chrono::high_resolution_clock::now();
  time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
  printf("  time: entries %d nthread %d time %f\n", 128, 1, time_span.count());

  t1 = std::chrono::high_resolution_clock::now();
  timeSTDQueueSerial("TIME STD 1 thread", std::move(std::queue<int>()), count);
  t2 = std::chrono::high_resolution_clock::now();
  time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
  printf("  time: entries %d nthread %d time %f\n", 128, 1, time_span.count());


  t1 = std::chrono::high_resolution_clock::now();
  timeTSQueueThreaded("TIME TSQ 4 threads ", std::move(QueueType()), count, 4, 1);
  t2 = std::chrono::high_resolution_clock::now();
  time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
  printf("  time: entries %d nthread %d time %f\n", 128, 4, time_span.count());

  t1 = std::chrono::high_resolution_clock::now();
  timeSTDQueueThreaded("TIME STD 4 thread", std::move(std::queue<int>()), count, 4, 1);
  t2 = std::chrono::high_resolution_clock::now();
  time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
  printf("  time: entries %d nthread %d time %f\n", 128, 4, time_span.count());




};
