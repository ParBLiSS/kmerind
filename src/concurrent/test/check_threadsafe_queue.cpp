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
void testSTDQueueSerial(const std::string &message, std::queue<T>&& q, const int entries) {
  T output;
  T result = 0;
  for (int i = 0; i < entries; ++i) {

    if (i % 2 == 0)
      q.push(T(i));
    else {
      output = q.front();
      q.pop();
      result ^= output;
    }
  }
  printf("result: %d\n", result);

};

template<typename T>
void testSTDQueueThreaded(const std::string &message, std::queue<T>&& q, const int entries, const int nProducer, const int nConsumer) {
  T output;
  T result = 0;

  omp_lock_t lock;

  omp_init_lock(&lock);

#pragma omp parallel sections num_threads(2) shared(queue, entries)
  {
#pragma omp section
    {
      int localCount = entries;
      // insert
#pragma omp parallel num_threads(nProducer) shared(queue, localCount)
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
          q.push(T(count));
          omp_unset_lock(&lock);
        }
      }
    }


#pragma omp section
    {
      int localCount = entries;
      T result;
      // insert
#pragma omp parallel num_threads(nProducer) shared(queue, localCount, result)
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
          output = q.front();
          q.pop();
          omp_unset_lock(&lock);

#pragma omp atomic
          result ^= output;
        }
      }
      printf("result: %d\n", result);


    }
  }

  omp_destroy_lock(&lock);
};





template<typename T>
void testWaitAndPush(bliss::concurrent::ThreadSafeQueue<T> &queue, const int entries, const int nProducer, bool& done) {
  usleep(1000);
#pragma omp parallel for default(none) num_threads(nProducer) shared(queue, entries)
  for (int i = 0; i < entries; ++i) {
    queue.waitAndPush(T(i));
//    if (i % 1000 == 0) usleep(5);
    usleep(5);
  }
  done = true;
#pragma omp flush(done)

  printf("TSQueue finished waitAndPush.\n");
};

template<typename T>
void testTryPush(bliss::concurrent::ThreadSafeQueue<T> &queue, const int entries, const int nProducer, bool& done) {
  usleep(1000);
  int count = 0, count2 = 0;
#pragma omp parallel for default(none) num_threads(nProducer) shared(queue, entries) reduction(+: count)
  for (int i = 0; i < entries; ++i) {
    if (queue.tryPush(T(i)))
      ++count;
    else
      ++count2;
//    if (i % 1000 == 0) usleep(20);
    usleep(5)
  }
  done = true;
#pragma omp flush(done)

  printf("TSQueue finished tryPush. %d successful, %d failed\n", count, count2);

};

template<typename T>
void testWaitAndPop1(bliss::concurrent::ThreadSafeQueue<T> &queue, bool &done) {
  int count = 0, count3 = 0;

#pragma omp flush(done)
  while (!done) {
    queue.waitAndPop();
    ++count;

    //if (count %1000) usleep(20);
#pragma omp flush(done)
  }
  // now empty it
  while (queue.tryPop().first) {
    ++count3;
  }
  printf("TSQueue finished waitAndPop 1 thread. %d successful, %d flushed.  empty? %s\n", count, count3, (queue.empty() ? "yes" : "no"));
};

template<typename T>
void testTryPop1(bliss::concurrent::ThreadSafeQueue<T> &queue, bool &done) {
  int count = 0, count2 = 0, count3 = 0;

#pragma omp flush(done)
  while (!done) {
    if (queue.tryPop().first)
      ++count;
    else
      ++count2;

    //if (count %1000) usleep(20);
#pragma omp flush(done)
  }
  // now empty it
  while (queue.tryPop().first) {
    ++count3;
  }
  printf("TSQueue finished tryPop 1 thread. %d successful, %d failed, %d flushed.  empty? %s\n", count, count2, count3, (queue.empty() ? "yes" : "no"));

};

template<typename T>
void testTryPop(bliss::concurrent::ThreadSafeQueue<T> &queue, const int nConsumer, bool &done ) {


  int counts[nConsumer];
  int counts2[nConsumer];
  int counts3[nConsumer];

  // while producer is producing, don't know the size.  spawn tasks.
#pragma omp parallel num_threads(nConsumer) default(none) shared(queue, done, counts, counts2)
  {
    counts[omp_get_thread_num()] = 0;
    counts2[omp_get_thread_num()] = 0;
    counts3[omp_get_thread_num()] = 0;

#pragma omp single nowait
    {
      do {

        // deferred task to pop.
#pragma omp task default(none) shared(queue, counts, counts2)
        {
          if (queue.tryPop().first)
            ++counts[omp_get_thread_num()];
          else
            ++counts2[omp_get_thread_num()];

          //if ((counts[omp_get_thread_num()] % 1000) == 0) usleep(20);
        }


#pragma omp flush(done)
      } while (!done);
    }

  }  // implicit barrier.  won't have issue with reporting too early.
  // NO WASTED TIME here
  // tasks create and complete fast, then some wasted time in flush task queue.  data queue may be empty by now.
  // tasks slow but task create fast, then task queue long, take extra time to flush task queue, data queue may be empty here
  // or tasks slow and task create slow, then no wasted task or time in this phase.  data queue not empty
  // or tasks fast and task create slow, then task queue short, and data queue not empty

  // now use parallel for to flush the data queue
  size_t size = queue.size();
#pragma omp parallel for num_threads(nConsumer) default(none) shared(queue, counts2, counts3)
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
  printf("TSQueue finished tryPop. %d successful, %d failed, %d flushed.  empty? %s\n", count, count2, count3, (queue.empty() ? "yes" : "no"));

};

// no waitAndPop.  if task creation is too fast, then we may have more waitAndPop calls than there are entries in queue, then deadlock.



template<typename T>
void testTSQueue(const std::string &message, bliss::concurrent::ThreadSafeQueue<T>&& queue, const int nProducer, const int nConsumer) {

  printf("TEST %s: %d producers, %d consumers\n", message.c_str(), nProducer, nConsumer);

  int entries = (queue.getCapacity() == bliss::concurrent::ThreadSafeQueue<T>::MAX_SIZE) ? 128 : queue.getCapacity();

  int i = 0;
  int count = 0, count2 = 0;


  // check tryPop on empty
  queue.clear();
#pragma omp parallel for num_threads(nConsumer) private(i) shared(queue, entries) default(none) reduction(+: count, count2)
  for (i = 0; i < entries; ++i) {
    if ( queue.tryPop().first )
      ++count;
    else
      ++count2;
  }
  printf("TSFixedQueue finished tryPop on empty at iteration %d, success %d, fail %d\n", i, count, count2);

  // check tryPush too much
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
  printf("TSFixedQueue finished tryPush until full at iteration %d, success %d, fail %d\n", i, count, count2);

  // check tryPop too much
  count = 0;
  count2 = 0;
#pragma omp parallel for num_threads(nConsumer) private(i) shared(queue, entries) default(none) reduction(+: count, count2)
  for (i = 0; i < (entries + 2); ++i) {
    if ( queue.tryPop().first)
      ++count;
    else
      ++count2;
  }
  printf("TSFixedQueue finished tryPop from full at iteration %d, success %d, fail %d\n", i, count, count2);



  // check  waitAndPush then do waitAndPop in parallel
  queue.clear();
  count = 0;
#pragma omp parallel for num_threads(nProducer) private(i) shared(queue, entries) default(none) reduction(+: count, count2)
  for (i = 0; i < entries; ++i) {
    queue.waitAndPush(T(i));
    ++count;
  }
  printf("TSFixedQueue finished waitAndPush with %d entries\n", count);

  count2 = 0;
#pragma omp parallel for num_threads(nthreads) private(i) shared(queue, entries) default(none) reduction(+: count)
  for (i = 0; i < count; ++i) {
    T output = queue.waitAndPop();
    ++count2;
  }
  printf("TSFixedQueue finished waitAndPop with, %d entries\n", count2);




  // check tryPush, and tryPop.
  bool done = false;
#pragma omp parallel sections num_threads(2) shared(queue, entries, done, nProducer, nConsumer) default(none)
  {
#pragma omp section
    {
      testTryPush(queue, entries, nProducer, done);
    }

#pragma omp section
    {
      if (nConsumer > 1)
        testTryPop(queue, nConsumer, done);
      else
        testTryPop1(queue, done);
    }
  }



  // check waitAndPush, and tryPop.
  done = false;
#pragma omp parallel sections num_threads(2) shared(queue, entries, done, nProducer, nConsumer) default(none)
  {
#pragma omp section
    {
      testWaitAndPush(queue, entries, nProducer, done);
    }

#pragma omp section
    {
      if (nConsumer > 1)
        testTryPop(queue, nConsumer, done);
      else
        testTryPop1(queue, done);
    }
  }

  if (nConsumer == 1) {

    // not testing this with more than 1 consumer thread.  there could be a lot more consumer tasks generated than
    //   there are elements because the produce thread is generating slower than consume threads.  the excess threads result in wait forever.
    // this could still happen if we slow down producer and speed up consumers



    // check tryPush, and waitAndPop.
    done = false;
  #pragma omp parallel sections num_threads(2) shared(queue, entries, done, nProducer, nConsumer) default(none)
    {
  #pragma omp section
      {
        testTryPush(queue, entries, nProducer, done);
      }

  #pragma omp section
      {
        testWaitAndPop1(queue, done);
      }
    }



    // check waitAndPush, and waitAndPop.
    done = false;
  #pragma omp parallel sections num_threads(2) shared(queue, entries, done, nProducer, nConsumer) default(none)
    {
  #pragma omp section
      {
        testWaitAndPush(queue, entries, nProducer, done);
      }

  #pragma omp section
      {
        testWaitAndPop1(queue, done);
      }
    }
  }
};


template<typename T>
void testTSQueueSingleThread(const std::string &message, bliss::concurrent::ThreadSafeQueue<T>&& q) {
  T output;
  T result = 0;
  int entries = (q.getCapacity() == bliss::concurrent::ThreadSafeQueue<T>::MAX_SIZE) ? 128 : q.getCapacity();

  for (int i = 0; i < entries; ++i) {

    if (i % 2 == 0)
      q.tryPush(std::move(T(i)));
    else
    {
      q.tryPop(output);
      result ^= output;
    }
  }
  printf("result: %d\n", result);
};






int main(int argc, char** argv) {

  int elements = 100000;
  std::chrono::high_resolution_clock::time_point t1, t2;
  std::chrono::duration<double> time_span;

  omp_set_nested(1);
  omp_set_dynamic(0);




  typedef QueueTest<int> QTestType;
  for (int i = 1; i <=4; ++i) {
    QTestType qt(i);
    qt.run(elements, qt.getQueue());

    t1 = std::chrono::high_resolution_clock::now();
    qt.run2(elements);
    t2 = std::chrono::high_resolution_clock::now();
    time_span =
        std::chrono::duration_cast<std::chrono::duration<double>>(
            t2 - t1);

    printf("queueTest nthread %d time %f\n", i, time_span.count());

  }

  typedef FixedSizeQueueTest<int, 1000> FSQTestType;
  for (int i = 1; i <=4; ++i) {
    FSQTestType qt(i);
    qt.run(elements,  qt.getQueue());

    t1 = std::chrono::high_resolution_clock::now();
    qt.run2(elements);
    t2 = std::chrono::high_resolution_clock::now();
    time_span =
        std::chrono::duration_cast<std::chrono::duration<double>>(
            t2 - t1);

    printf("fixed queueTest size %d nthread %d time %f\n", 1000, i, time_span.count());

  }


  typedef STDQueueTest<int> SQTestType;
  SQTestType qt;

  t1 = std::chrono::high_resolution_clock::now();
  qt.run2(elements);
  t2 = std::chrono::high_resolution_clock::now();
  time_span =
      std::chrono::duration_cast<std::chrono::duration<double>>(
          t2 - t1);

  printf("std queueTest time %f\n", time_span.count());


};
