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

#include "bliss-config.hpp"

#include <iostream>
#include <unistd.h>

#include "utils/logging.h"


#if defined( BLISS_MUTEX)
#include "concurrent/mutexlock_queue.hpp"
#elif defined( BLISS_SPINLOCK )
#include "concurrent/spinlock_queue.hpp"
#else   //if defined( BLISS_LOCKFREE )
#include "concurrent/lockfree_queue.hpp"
#endif

#include "omp.h"

using namespace bliss::concurrent;





template<typename T, bliss::concurrent::LockType LT>
void testWaitAndPush(bliss::concurrent::ThreadSafeQueue<T, LT> &queue, const int entries, const int nProducer) {
  //usleep(1000);
  int count = 0;
#pragma omp parallel for OMP_SHARE_DEFAULT num_threads(nProducer) shared(queue) reduction(+:count)
  for (int i = 0; i < entries; ++i) {
    if (queue.waitAndPush(T(i)).first)
//    if (i % 1000 == 0) usleep(5);
      ++count;
    //usleep(5);
  }
  queue.disablePush();

  if (count != entries) FATALF("FAIL: TSQueue capacity %lu, finished waitAndPush with %d entries, expected %d.", queue.getCapacity(), count, entries);
  else INFOF("PASS,");
};


template<typename T, bliss::concurrent::LockType LT>
void testWaitAndPushSome(bliss::concurrent::ThreadSafeQueue<T, LT> &queue, const int entries, const int nProducer) {
  //usleep(1000);
  int count = 0;
#pragma omp parallel for OMP_SHARE_DEFAULT num_threads(nProducer) shared(queue) reduction(+:count)
  for (int i = 0; i < entries; ++i) {
    if (queue.waitAndPush(T(i)).first)
//    if (i % 1000 == 0) usleep(5);
      ++count;
    //usleep(5);
  }
  queue.disablePush();


  INFOF("PASS, pushed %d of %d", count, entries);
};

template<typename T, bliss::concurrent::LockType LT>
void testTryPush(bliss::concurrent::ThreadSafeQueue<T, LT> &queue, const int entries, const int nProducer) {
  //usleep(1000);
  int count = 0, count2 = 0;
#pragma omp parallel for OMP_SHARE_DEFAULT num_threads(nProducer) shared(queue) reduction(+: count, count2)
  for (int i = 0; i < entries; ++i) {
    if (queue.tryPush(T(i)).first)
      ++count;
    else
      ++count2;
//    if (i % 1000 == 0) usleep(20);
    //usleep(5);
  }
  queue.disablePush();

  if (count + count2 != entries || count == 0)
    FATALF("FAIL: TSQueue capacity %lu, finished tryPush. %d successful, %d failed", queue.getCapacity(), count, count2);
  else INFOF("PASS,");

};

template<typename T, bliss::concurrent::LockType LT>
void testWaitAndPop1(bliss::concurrent::ThreadSafeQueue<T, LT> &queue) {
  int count = 0, count3 = 0;

  while (queue.canPop()) {
//    INFOF("queue waitAndPopping. "); fflush(stdout);
    if (queue.waitAndPop().first)
      ++count;
//    INFOF("queue size = %lu", queue.getSize());  fflush(stdout);

    //if (count %1000) usleep(20);
  }
      INFOF("done with waitAndPop()");  fflush(stdout);


  // now empty it
  while (queue.tryPop().first) {
    ++count3;

    INFOF("emptying. queue size = %lu", queue.getSize());

  }
  if (!queue.isEmpty() || count3 != 0) FATALF("FAIL: TSQueue capacity %lu, finished waitAndPop 1 thread. %d successful, %d flushed.  empty? %s", queue.getCapacity(), count, count3, (queue.isEmpty() ? "yes" : "no"));
  else INFOF("PASS,");

};

template<typename T, bliss::concurrent::LockType LT>
void testTryPop1(bliss::concurrent::ThreadSafeQueue<T, LT> &queue) {
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
  if (!queue.isEmpty() || (count == 0 && count3 == 0)) FATALF("FAIL: TSQueue capacity %lu, finished tryPop 1 thread. %d successful, %d failed, %d flushed.  empty? %s", queue.getCapacity(), count, count2, count3, (queue.isEmpty() ? "yes" : "no"));
  else INFOF("PASS,");

};

template<typename T, bliss::concurrent::LockType LT>
void testTryPop(bliss::concurrent::ThreadSafeQueue<T, LT> &queue, const int nConsumer) {


  int counts[nConsumer];
  int counts2[nConsumer];
  int counts3[nConsumer];

  // while producer is producing, don't know the size.  spawn tasks.
#pragma omp parallel num_threads(nConsumer) OMP_SHARE_DEFAULT shared(queue, counts, counts2, counts3)
  {
    counts[omp_get_thread_num()] = 0;
    counts2[omp_get_thread_num()] = 0;
    counts3[omp_get_thread_num()] = 0;

//#pragma omp single nowait
//    {
      while(queue.canPop()) {

        // deferred task to pop.
//#pragma omp task OMP_SHARE_DEFAULT shared(queue, counts, counts2)
//        {
          if (queue.tryPop().first)
            ++counts[omp_get_thread_num()];
          else
            ++counts2[omp_get_thread_num()];

          //if ((counts[omp_get_thread_num()] % 1000) == 0) usleep(20);
//        }

      };
//    }

  }  // implicit barrier.  won't have issue with reporting too early.
  // NO WASTED TIME here
  // tasks create and complete fast, then some wasted time in flush task queue.  data queue may be empty by now.
  // tasks slow but task create fast, then task queue long, take extra time to flush task queue, data queue may be empty here
  // or tasks slow and task create slow, then no wasted task or time in this phase.  data queue not empty
  // or tasks fast and task create slow, then task queue short, and data queue not empty

  // now use parallel for to flush the data queue
  size_t size = queue.getSize();
#pragma omp parallel for num_threads(nConsumer) OMP_SHARE_DEFAULT shared(queue, size, counts2, counts3)
  for (size_t i = 0; i < size; ++i) {
    if (queue.tryPop().first)
      ++counts3[omp_get_thread_num()];
    else
      ++counts2[omp_get_thread_num()];
  }


  // summarize
  int count = 0, count2 = 0, count3 = 0;
#pragma omp parallel num_threads(nConsumer) OMP_SHARE_DEFAULT shared(counts, counts2, counts3) reduction(+:count, count2, count3)
  {
    count = counts[omp_get_thread_num()];
    count2 = counts2[omp_get_thread_num()];
    count3 = counts3[omp_get_thread_num()];
  }
  if (!queue.isEmpty() || (count == 0 && count3 == 0)) FATALF("FAIL: TSQueue capacity %lu, finished tryPop. %d successful, %d failed, %d flushed.  empty? %s", queue.getCapacity(), count, count2, count3, (queue.isEmpty() ? "yes" : "no"));
  else INFOF("PASS,");

};

// no waitAndPop.  if task creation is too fast, then we may have more waitAndPop calls than there are entries in queue, then deadlock.



template<typename T, bliss::concurrent::LockType LT>
void testTSQueue(const std::string &message, bliss::concurrent::ThreadSafeQueue<T, LT>&& queue, const int nProducer, const int nConsumer) {

  size_t entries = (!queue.isFixedSize()) ? 10000 : queue.getCapacity();

  INFOF("=== TEST %s: %d producers, %d consumers, capacity %lu, entries %d", message.c_str(), nProducer, nConsumer, queue.getCapacity(), entries);


  size_t i = 0;
  size_t count = 0, count2 = 0;


  INFOF("  CHECK tryPop on empty: ");  fflush(stdout);
  queue.clear();
#pragma omp parallel for num_threads(nConsumer) private(i) shared(queue, entries) OMP_SHARE_DEFAULT reduction(+: count, count2)
  for (i = 0; i < entries; ++i) {
    if ( queue.tryPop().first )
      ++count;
    else
      ++count2;
  }
  if (count2 != entries) FATALF("FAIL: TSQueue capacity %lu, finished tryPop on empty at iteration %lu, success %lu, fail %lu", queue.getCapacity(), i, count, count2);
  else INFOF("PASS");

  INFOF("  CHECK tryPush too much: ");  fflush(stdout);
  count = 0;
  count2 = 0;
  queue.clear();
#pragma omp parallel for num_threads(nProducer) private(i) shared(queue, entries) OMP_SHARE_DEFAULT reduction(+: count, count2)
  for (i = 0; i < (entries + 2); ++i) {
    if (queue.tryPush(T(i)).first)
      ++count;
    else
      ++count2;
  }
  auto expected = (!queue.isFixedSize()) ? entries+2 : entries;
  if (count != std::min(queue.getCapacity(), 2UL + entries) || (count + count2 != (entries+2)))  FATALF("FAIL: TSQueue capacity %lu, finished tryPush until full, expected %lu, success %lu, fail %lu. ", queue.getCapacity(), expected, count, count2);
  else INFOF("PASS");

  INFOF("  CHECK tryPop too much: ");  fflush(stdout);
  count = 0;
  count2 = 0;
#pragma omp parallel for num_threads(nConsumer) private(i) shared(queue, entries) OMP_SHARE_DEFAULT reduction(+: count, count2)
  for (i = 0; i < (entries + 2); ++i) {
    if ( queue.tryPop().first)
      ++count;
    else
      ++count2;
  }
  expected = (!queue.isFixedSize()) ? entries+2 : entries;
  if ((count + count2 != (entries+2))) FATALF("FAIL: TSQueue capacity %lu, finished tryPop from full, expected %lu, success %lu, fail %lu",queue.getCapacity(),  expected, count, count2);
  else INFOF("PASS");



  INFOF("  CHECK waitAndPush then do waitAndPop in parallel: ");   fflush(stdout);
  queue.clear();
  queue.enablePush();
  count = 0;
#pragma omp parallel for num_threads(nProducer) private(i) shared(queue, entries) OMP_SHARE_DEFAULT reduction(+: count, count2)
  for (i = 0; i < entries; ++i) {
    if (queue.waitAndPush(T(i)).first)
      ++count;
  }

  count2 = 0;
#pragma omp parallel for num_threads(nConsumer) private(i) shared(queue, count) OMP_SHARE_DEFAULT reduction(+: count2)
  for (i = 0; i < count; ++i) {
    if (queue.waitAndPop().first)
      ++count2;
  }
  if (count != entries ||  count2 != count) FATALF("FAIL: TSQueue capacity %lu, finished waitAndPush with %lu entries, finished waitAndPop with %lu entries.  should be %lu", queue.getCapacity(), count, count2, entries);
  else INFOF("PASS");




  INFOF("  CHECK tryPush, and tryPop: ");  fflush(stdout);
  queue.clear();
  queue.enablePush();
#pragma omp parallel sections num_threads(2) shared(queue, entries) OMP_SHARE_DEFAULT
  {
#pragma omp section
    {
      testTryPush(queue, entries, nProducer);
      //INFOF("tryPush done correctly.  final queue size (with consumer): %lu", queue.size());
    }

#pragma omp section
    {
      if (nConsumer > 1)
        testTryPop(queue, nConsumer);
      else
        testTryPop1(queue);
    }
  }


  INFOF("  CHECK waitAndPush, and tryPop: ");   fflush(stdout);
  queue.clear();
  queue.enablePush();
#pragma omp parallel sections num_threads(2) shared(queue, entries) OMP_SHARE_DEFAULT
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


  INFOF("  CHECK waitAndPush, and disablePush: ");   fflush(stdout);
  queue.clear();
  queue.enablePush();
#pragma omp parallel sections num_threads(2) shared(queue, entries) OMP_SHARE_DEFAULT
  {
#pragma omp section
    {
      testWaitAndPushSome(queue, entries, nProducer);
    }

#pragma omp section
    {
	usleep(1000);
	queue.disablePush();
//#pragma omp flush(queue)
    }
  }

  //TODO: can have !done, waitAndPop, and done=true in other thread, so pop thread never gets to check "done" ->  deadlock.

    // not testing this with more than 1 consumer thread.  there could be a lot more consumer tasks generated than
    //   there are elements because the produce thread is generating slower than consume threads.  the excess threads result in wait forever.
    // this could still happen if we slow down producer and speed up consumers



    INFOF("  CHECK tryPush, and waitAndPop: ");  fflush(stdout);
    queue.clear();
    queue.enablePush();
  #pragma omp parallel sections num_threads(2) shared(queue, entries) OMP_SHARE_DEFAULT
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




    INFOF("  CHECK waitAndPush, and waitAndPop: ");  fflush(stdout);
    queue.clear();
    queue.enablePush();
  #pragma omp parallel sections num_threads(2) shared(queue, entries) OMP_SHARE_DEFAULT
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



int main(int argc, char** argv) {


  omp_set_nested(1);
  omp_set_dynamic(0);

#if defined(BLISS_MUTEX)
  constexpr bliss::concurrent::LockType lt = bliss::concurrent::LockType::MUTEX;
#elif defined(BLISS_SPINLOCK)
  constexpr bliss::concurrent::LockType lt = bliss::concurrent::LockType::SPINLOCK;
#else
  constexpr bliss::concurrent::LockType lt = bliss::concurrent::LockType::LOCKFREE;
#endif

  typedef bliss::concurrent::ThreadSafeQueue<int, lt> QueueType;

  testTSQueue("TSQ nthread 100 elements", QueueType(100), 1, 1);
  testTSQueue("TSQ nthread 100 elements", QueueType(100), 2, 1);
  testTSQueue("TSQ nthread 100 elements", QueueType(100), 3, 1);
  testTSQueue("TSQ nthread 100 elements", QueueType(100), 4, 1);
  testTSQueue("TSQ nthread 100 elements", QueueType(100), 1, 2);
  testTSQueue("TSQ nthread 100 elements", QueueType(100), 1, 3);
  testTSQueue("TSQ nthread 100 elements", QueueType(100), 1, 4);
  testTSQueue("TSQ nthread 100 elements", QueueType(100), 2, 2);
  testTSQueue("TSQ nthread 100 elements", QueueType(100), 2, 3);
  testTSQueue("TSQ nthread 100 elements", QueueType(100), 3, 2);

  testTSQueue("TSQ nthread growable", QueueType(), 1, 1);
  testTSQueue("TSQ nthread growable", QueueType(), 2, 1);
  testTSQueue("TSQ nthread growable", QueueType(), 3, 1);
  testTSQueue("TSQ nthread growable", QueueType(), 4, 1);
  testTSQueue("TSQ nthread growable", QueueType(), 1, 2);
  testTSQueue("TSQ nthread growable", QueueType(), 1, 3);
  testTSQueue("TSQ nthread growable", QueueType(), 1, 4);
  testTSQueue("TSQ nthread growable", QueueType(), 2, 2);
  testTSQueue("TSQ nthread growable", QueueType(), 2, 3);
  testTSQueue("TSQ nthread growable", QueueType(), 3, 2);



};
