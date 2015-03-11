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
#include <unistd.h>

#include <limits>
#include <algorithm>

#include "concurrent/mutexlock_queue.hpp"
#include "concurrent/spinlock_queue.hpp"
#include "concurrent/lockfree_queue.hpp"

#include "omp.h"

using namespace bliss::concurrent;

// defined here so concurrent queue does not deadlock.
#if defined(BLISS_LOCKFREE)
#define THREAD_EXIT_NOTIFIER_START RelacyThreadExitNotifier::notify_relacy_thread_start()
#define THREAD_EXIT_NOTIFIER_END  RelacyThreadExitNotifier::notify_relacy_thread_exit()
#else
#define THREAD_EXIT_NOTIFIER_START
#define THREAD_EXIT_NOTIFIER_END
#endif



// note that relacy simulate threads.
template<typename T, bliss::concurrent::LockType LT, int nProducers, int nConsumers, int64_t capacity, int64_t entries>
struct testWaitAndPush :
    rl::test_suite< testWaitAndPush<T, LT, nProducers, nConsumers, capacity, entries>,
    nProducers + nConsumers >
{

      bliss::concurrent::ThreadSafeQueue<T, LT> queue;
      std::array<int64_t, nProducers + nConsumers> lsuccess;
      std::array<int64_t, nProducers + nConsumers> lfail;
      static constexpr int64_t expected = (capacity > ((entries/nProducers) * nProducers) ? ((entries/nProducers) * nProducers) : capacity);


      testWaitAndPush() : queue(capacity) {};



      void before() {
        for (int i = 0; i < nProducers + nConsumers; ++i) {
          lsuccess[i] = 0;
          lfail[i] = 0;
        }
      }


      void thread(unsigned thread_index)
      {
        THREAD_EXIT_NOTIFIER_START;

        if (thread_index < nProducers) {
          for (int64_t i = 0; i < entries / nProducers; ++i) {
            if (queue.waitAndPush(T(i * nProducers + thread_index)).first)
              ++lsuccess[thread_index];
            else
              ++lfail[thread_index];
          }
        } else { //if (thread_index == (nProducers + nConsumers - 1)){
          while (queue.getSize() < expected);
          queue.disablePush();
        }

        THREAD_EXIT_NOTIFIER_END;
      }

      void after() {
        int64_t success = 0;
        int64_t fail = 0;
        queue.disablePush();
        for (int i = 0; i < nProducers + nConsumers; ++i) {
          success += lsuccess[i];
          fail += lfail[i];
        }

        //INFOF("Test WaitAndPush: LT %d nProd %d nCons %d cap %ld elems %ld.  capacity %ld, expected %lu  actual %ld, failed %ld", LT, nProducers, nConsumers, queue.getCapacity(), entries, capacity, expected, success, fail);
        assert(success == expected);

      }

      void invariant() {

      }

};

template<typename T, bliss::concurrent::LockType LT, int nProducers, int nConsumers, int64_t capacity, int64_t entries>
struct testTryPush :
    rl::test_suite< testTryPush<T, LT, nProducers, nConsumers, capacity, entries>,
    nProducers + nConsumers >
{

      bliss::concurrent::ThreadSafeQueue<T, LT> queue;
      std::array<int64_t, nProducers + nConsumers> lsuccess;
      std::array<int64_t, nProducers + nConsumers> lfail;
      static constexpr int64_t expected = (capacity > ((entries/nProducers) * nProducers) ? ((entries/nProducers) * nProducers) : capacity);


      testTryPush() : queue(capacity) {};



      void before() {
        for (int i = 0; i < nProducers + nConsumers; ++i) {
          lsuccess[i] = 0;
          lfail[i] = 0;
        }
      }


      void thread(unsigned thread_index)
      {
        THREAD_EXIT_NOTIFIER_START;

        if (thread_index < nProducers) {
          for (int64_t i = 0; i < entries / nProducers; ++i) {
            if (queue.tryPush(T(i * nProducers + thread_index)).first)
              ++lsuccess[thread_index];
            else
              ++lfail[thread_index];
          }
        } else { //if (thread_index == (nProducers + nConsumers - 1)){
          while (queue.getSize() < nProducers * 2);
          queue.disablePush();
        }

        THREAD_EXIT_NOTIFIER_END;
      }

      void after() {
        int64_t count = 0;
        int64_t fail = 0;
        queue.disablePush();
        for (int i = 0; i < nProducers + nConsumers; ++i) {
          count += lsuccess[i];
          fail += lfail[i];
        }

        //INFOF("Test WaitAndPush: LT %d nProd %d nCons %d cap %ld elems %ld.  capacity %ld, expected %lu  actual %ld, failed %ld", LT, nProducers, nConsumers, queue.getCapacity(), entries, capacity, expected, count, fail);
        assert(count == expected);

      }

      void invariant() {

      }

};


    template<typename T, bliss::concurrent::LockType LT, int nProducers, int nConsumers, int64_t capacity, int64_t entries>
    struct testWaitAndPop :
        rl::test_suite< testWaitAndPop<T, LT, nProducers, nConsumers, capacity, entries>,
        nProducers + nConsumers >
    {

          bliss::concurrent::ThreadSafeQueue<T, LT> queue;
          std::array<int64_t, nProducers + nConsumers> lsuccess;
          std::array<int64_t, nProducers + nConsumers> lfail;
          static constexpr int64_t expected = entries;


          testWaitAndPop() : queue(capacity) {};



          void before() {
            for (int i = 0; i < nProducers + nConsumers; ++i) {
              lsuccess[i] = 0;
              lfail[i] = 0;
            }

            for (int64_t i = 0; i < entries; ++i) {
              queue.waitAndPush(T(i));
            }
            queue.disablePush();
          }


          void thread(unsigned thread_index)
          {
            THREAD_EXIT_NOTIFIER_START;

            if (thread_index >= nProducers) {
              while (queue.canPop()) {
                if (queue.waitAndPop().first)
                  ++lsuccess[thread_index];
                else
                  ++lfail[thread_index];
              }
            }

            THREAD_EXIT_NOTIFIER_END;
          }

          void after() {
            int64_t count = 0;
            int64_t fail = 0;
            queue.disablePush();
            for (int i = 0; i < nProducers + nConsumers; ++i) {
              count += lsuccess[i];
              fail += lfail[i];
            }

            //INFOF("Test WaitAndPush: LT %d nProd %d nCons %d cap %ld elems %ld.  capacity %ld, expected %lu  actual %ld, failed %ld", LT, nProducers, nConsumers, queue.getCapacity(), entries, capacity, expected, count, fail);
            assert(count == expected);

          }

          void invariant() {

          }

    };



        template<typename T, bliss::concurrent::LockType LT, int nProducers, int nConsumers, int64_t capacity, int64_t entries>
        struct testTryPop :
            rl::test_suite< testTryPop<T, LT, nProducers, nConsumers, capacity, entries>,
            nProducers + nConsumers >
        {

              bliss::concurrent::ThreadSafeQueue<T, LT> queue;
              std::array<int64_t, nProducers + nConsumers> lsuccess;
              std::array<int64_t, nProducers + nConsumers> lfail;
              static constexpr int64_t expected = entries;


              testTryPop() : queue(capacity) {};



              void before() {
                for (int i = 0; i < nProducers + nConsumers; ++i) {
                  lsuccess[i] = 0;
                  lfail[i] = 0;
                }

                for (int64_t i = 0; i < entries; ++i) {
                  queue.waitAndPush(T(i));
                }
                queue.disablePush();
              }


              void thread(unsigned thread_index)
              {
                THREAD_EXIT_NOTIFIER_START;

                if (thread_index >= nProducers) {
                  while (queue.canPop()) {
                    if (queue.tryPop().first)
                      ++lsuccess[thread_index];
                    else
                      ++lfail[thread_index];
                  }
                }

                THREAD_EXIT_NOTIFIER_END;
              }

              void after() {
                int64_t count = 0;
                int64_t fail = 0;
                queue.disablePush();
                for (int i = 0; i < nProducers + nConsumers; ++i) {
                  count += lsuccess[i];
                  fail += lfail[i];
                }

                //INFOF("Test WaitAndPush: LT %d nProd %d nCons %d cap %ld elems %ld.  capacity %ld, expected %lu  actual %ld, failed %ld", LT, nProducers, nConsumers, queue.getCapacity(), entries, capacity, expected, count, fail);
                assert(count == expected);

              }

              void invariant() {

              }

        };

/*



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
  if (!queue.isEmpty() || (count == 0 && count3 == 0)) FATALF("FAIL: TSQueue capacity %lu, finished tryPop 1 thread. %d successful, %d failed, %d flushed.  empty? %s\n", queue.getCapacity(), count, count2, count3, (queue.isEmpty() ? "yes" : "no"));
  else INFOF("PASS,");

};

template<typename T, bliss::concurrent::LockType LT>
void testTryPop(bliss::concurrent::ThreadSafeQueue<T, LT> &queue, const int nConsumer) {


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
  if (!queue.isEmpty() || (count == 0 && count3 == 0)) FATALF("FAIL: TSQueue capacity %lu, finished tryPop. %d successful, %d failed, %d flushed.  empty? %s\n", queue.getCapacity(), count, count2, count3, (queue.isEmpty() ? "yes" : "no"));
  else INFOF("PASS,");

};

// no waitAndPop.  if task creation is too fast, then we may have more waitAndPop calls than there are entries in queue, then deadlock.



template<typename T, bliss::concurrent::LockType LT>
void testTSQueue(const std::string &message, bliss::concurrent::ThreadSafeQueue<T, LT>&& queue, const int nProducer, const int nConsumer) {

  int entries = (!queue.isFixedSize()) ? 10000 : queue.getCapacity();

  INFOF("=== TEST %s: %d producers, %d consumers, capacity %lu, entries %d\n", message.c_str(), nProducer, nConsumer, queue.getCapacity(), entries);


  int i = 0;
  int count = 0, count2 = 0;


  INFOF("  CHECK tryPop on empty: ");  fflush(stdout);
  queue.clear();
#pragma omp parallel for num_threads(nConsumer) private(i) shared(queue, entries) default(none) reduction(+: count, count2)
  for (i = 0; i < entries; ++i) {
    if ( queue.tryPop().first )
      ++count;
    else
      ++count2;
  }
  if (count2 != entries) FATALF("FAIL: TSQueue capacity %lu, finished tryPop on empty at iteration %d, success %d, fail %d\n", queue.getCapacity(), i, count, count2);
  else INFOF("PASS\n");

  INFOF("  CHECK tryPush too much: ");  fflush(stdout);
  count = 0;
  count2 = 0;
  queue.clear();
#pragma omp parallel for num_threads(nProducer) private(i) shared(queue, entries) default(none) reduction(+: count, count2)
  for (i = 0; i < (entries + 2); ++i) {
    if (queue.tryPush(T(i)).first)
      ++count;
    else
      ++count2;
  }
  int expected = (!queue.isFixedSize()) ? entries+2 : entries;
  if (count != std::min(queue.getCapacity(), 2UL + entries) || (count + count2 != (entries+2)))  FATALF("FAIL: TSQueue capacity %lu, finished tryPush until full, expected %d, success %d, fail %d. \n", queue.getCapacity(), expected, count, count2);
  else INFOF("PASS\n");

  INFOF("  CHECK tryPop too much: ");  fflush(stdout);
  count = 0;
  count2 = 0;
#pragma omp parallel for num_threads(nConsumer) private(i) shared(queue, entries) default(none) reduction(+: count, count2)
  for (i = 0; i < (entries + 2); ++i) {
    if ( queue.tryPop().first)
      ++count;
    else
      ++count2;
  }
  expected = (!queue.isFixedSize()) ? entries+2 : entries;
  if (count != std::min(queue.getCapacity(), 2UL + entries) || (count + count2 != (entries+2))) FATALF("FAIL: TSQueue capacity %lu, finished tryPop from full, expected %d, success %d, fail %d\n",queue.getCapacity(),  expected, count, count2);
  else INFOF("PASS\n");



  INFOF("  CHECK waitAndPush then do waitAndPop in parallel: ");   fflush(stdout);
  queue.clear();
  queue.enablePush();
  count = 0;
#pragma omp parallel for num_threads(nProducer) private(i) shared(queue, entries) default(none) reduction(+: count, count2)
  for (i = 0; i < entries; ++i) {
    if (queue.waitAndPush(T(i)).first)
      ++count;
  }

  count2 = 0;
#pragma omp parallel for num_threads(nConsumer) private(i) shared(queue, count) default(none) reduction(+: count2)
  for (i = 0; i < count; ++i) {
    if (queue.waitAndPop().first)
      ++count2;
  }
  if (count != entries ||  count2 != count) FATALF("FAIL: TSQueue capacity %lu, finished waitAndPush with %d entries, finished waitAndPop with %d entries.  should be %d\n", queue.getCapacity(), count, count2, entries);
  else INFOF("PASS\n");




  INFOF("  CHECK tryPush, and tryPop: ");  fflush(stdout);
  queue.clear();
  queue.enablePush();
#pragma omp parallel sections num_threads(2) shared(queue, entries) default(none)
  {
#pragma omp section
    {
      testTryPush(queue, entries, nProducer);
      //INFOF("tryPush done correctly.  final queue size (with consumer): %lu\n", queue.size());
    }

#pragma omp section
    {
      if (nConsumer > 1)
        testTryPop(queue, nConsumer);
      else
        testTryPop1(queue);
    }
  }
  INFOF("\n");

  INFOF("  CHECK waitAndPush, and tryPop: ");   fflush(stdout);
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
  INFOF("\n");

  INFOF("  CHECK waitAndPush, and disablePush: ");   fflush(stdout);
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
	usleep(10000);
	queue.disablePush();
//#pragma omp flush(queue)
    }
  }
  INFOF("\n");
  //TODO: can have !done, waitAndPop, and done=true in other thread, so pop thread never gets to check "done" ->  deadlock.

    // not testing this with more than 1 consumer thread.  there could be a lot more consumer tasks generated than
    //   there are elements because the produce thread is generating slower than consume threads.  the excess threads result in wait forever.
    // this could still happen if we slow down producer and speed up consumers



    INFOF("  CHECK tryPush, and waitAndPop: ");  fflush(stdout);
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
    INFOF("\n");



    INFOF("  CHECK waitAndPush, and waitAndPop: ");  fflush(stdout);
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
    INFOF("\n");
};


*/
int main(int argc, char** argv) {

  int iterations = 1000;
  if (argc > 1) {
    iterations = atoi(argv[1]);
  }

  rl::test_params p;
  p.search_type = rl::sched_random;
  p.iteration_count = iterations;
  p.output_history = true;
  p.execution_depth_limit = 4000;

#if defined(BLISS_MUTEX)
  constexpr bliss::concurrent::LockType lt = bliss::concurrent::LockType::MUTEX;

#elif defined(BLISS_SPINLOCK)
  constexpr bliss::concurrent::LockType lt = bliss::concurrent::LockType::SPINLOCK;

#else
  constexpr bliss::concurrent::LockType lt = bliss::concurrent::LockType::LOCKFREE;

#endif




  rl::simulate<testWaitAndPush<int, lt, 1, 1, 20, 29> >(p);
  rl::simulate<testWaitAndPush<int, lt, 2, 1, 20, 29> >(p);
  rl::simulate<testWaitAndPush<int, lt, 3, 1, 20, 29> >(p);
  rl::simulate<testWaitAndPush<int, lt, 4, 1, 20, 29> >(p);
  rl::simulate<testWaitAndPush<int, lt, 1, 2, 20, 29> >(p);
  rl::simulate<testWaitAndPush<int, lt, 1, 3, 20, 29> >(p);
  rl::simulate<testWaitAndPush<int, lt, 1, 4, 20, 29> >(p);
  rl::simulate<testWaitAndPush<int, lt, 2, 2, 20, 29> >(p);
  rl::simulate<testWaitAndPush<int, lt, 2, 3, 20, 29> >(p);
  rl::simulate<testWaitAndPush<int, lt, 3, 2, 20, 29> >(p);

  rl::simulate<testTryPush<int, lt, 1, 1, 20, 29> >(p);
  rl::simulate<testTryPush<int, lt, 2, 1, 20, 29> >(p);
  rl::simulate<testTryPush<int, lt, 3, 1, 20, 29> >(p);
  rl::simulate<testTryPush<int, lt, 4, 1, 20, 29> >(p);
  rl::simulate<testTryPush<int, lt, 1, 2, 20, 29> >(p);
  rl::simulate<testTryPush<int, lt, 1, 3, 20, 29> >(p);
  rl::simulate<testTryPush<int, lt, 1, 4, 20, 29> >(p);
  rl::simulate<testTryPush<int, lt, 2, 2, 20, 29> >(p);
  rl::simulate<testTryPush<int, lt, 2, 3, 20, 29> >(p);
  rl::simulate<testTryPush<int, lt, 3, 2, 20, 29> >(p);

  rl::simulate<testWaitAndPop<int, lt, 1, 1, 20, 29> >(p);
  rl::simulate<testWaitAndPop<int, lt, 2, 1, 20, 29> >(p);
  rl::simulate<testWaitAndPop<int, lt, 3, 1, 20, 29> >(p);
  rl::simulate<testWaitAndPop<int, lt, 4, 1, 20, 29> >(p);
  rl::simulate<testWaitAndPop<int, lt, 1, 2, 20, 29> >(p);
  rl::simulate<testWaitAndPop<int, lt, 1, 3, 20, 29> >(p);
  rl::simulate<testWaitAndPop<int, lt, 1, 4, 20, 29> >(p);
  rl::simulate<testWaitAndPop<int, lt, 2, 2, 20, 29> >(p);
  rl::simulate<testWaitAndPop<int, lt, 2, 3, 20, 29> >(p);
  rl::simulate<testWaitAndPop<int, lt, 3, 2, 20, 29> >(p);

  rl::simulate<testTryPop<int, lt, 1, 1, 20, 29> >(p);
  rl::simulate<testTryPop<int, lt, 2, 1, 20, 29> >(p);
  rl::simulate<testTryPop<int, lt, 3, 1, 20, 29> >(p);
  rl::simulate<testTryPop<int, lt, 4, 1, 20, 29> >(p);
  rl::simulate<testTryPop<int, lt, 1, 2, 20, 29> >(p);
  rl::simulate<testTryPop<int, lt, 1, 3, 20, 29> >(p);
  rl::simulate<testTryPop<int, lt, 1, 4, 20, 29> >(p);
  rl::simulate<testTryPop<int, lt, 2, 2, 20, 29> >(p);
  rl::simulate<testTryPop<int, lt, 2, 3, 20, 29> >(p);
  rl::simulate<testTryPop<int, lt, 3, 2, 20, 29> >(p);

  rl::simulate<testWaitAndPush<int, lt, 1, 1, std::numeric_limits<int64_t>::max(), 29> >(p);
  rl::simulate<testWaitAndPush<int, lt, 2, 1, std::numeric_limits<int64_t>::max(), 29> >(p);
  rl::simulate<testWaitAndPush<int, lt, 3, 1, std::numeric_limits<int64_t>::max(), 29> >(p);
  rl::simulate<testWaitAndPush<int, lt, 4, 1, std::numeric_limits<int64_t>::max(), 29> >(p);
  rl::simulate<testWaitAndPush<int, lt, 1, 2, std::numeric_limits<int64_t>::max(), 29> >(p);
  rl::simulate<testWaitAndPush<int, lt, 1, 3, std::numeric_limits<int64_t>::max(), 29> >(p);
  rl::simulate<testWaitAndPush<int, lt, 1, 4, std::numeric_limits<int64_t>::max(), 29> >(p);
  rl::simulate<testWaitAndPush<int, lt, 2, 2, std::numeric_limits<int64_t>::max(), 29> >(p);
  rl::simulate<testWaitAndPush<int, lt, 2, 3, std::numeric_limits<int64_t>::max(), 29> >(p);
  rl::simulate<testWaitAndPush<int, lt, 3, 2, std::numeric_limits<int64_t>::max(), 29> >(p);

  rl::simulate<testTryPush<int, lt, 1, 1, std::numeric_limits<int64_t>::max(), 29> >(p);
  rl::simulate<testTryPush<int, lt, 2, 1, std::numeric_limits<int64_t>::max(), 29> >(p);
  rl::simulate<testTryPush<int, lt, 3, 1, std::numeric_limits<int64_t>::max(), 29> >(p);
  rl::simulate<testTryPush<int, lt, 4, 1, std::numeric_limits<int64_t>::max(), 29> >(p);
  rl::simulate<testTryPush<int, lt, 1, 2, std::numeric_limits<int64_t>::max(), 29> >(p);
  rl::simulate<testTryPush<int, lt, 1, 3, std::numeric_limits<int64_t>::max(), 29> >(p);
  rl::simulate<testTryPush<int, lt, 1, 4, std::numeric_limits<int64_t>::max(), 29> >(p);
  rl::simulate<testTryPush<int, lt, 2, 2, std::numeric_limits<int64_t>::max(), 29> >(p);
  rl::simulate<testTryPush<int, lt, 2, 3, std::numeric_limits<int64_t>::max(), 29> >(p);
  rl::simulate<testTryPush<int, lt, 3, 2, std::numeric_limits<int64_t>::max(), 29> >(p);


  rl::simulate<testWaitAndPop<int, lt, 1, 1, std::numeric_limits<int64_t>::max(), 29> >(p);
  rl::simulate<testWaitAndPop<int, lt, 2, 1, std::numeric_limits<int64_t>::max(), 29> >(p);
  rl::simulate<testWaitAndPop<int, lt, 3, 1, std::numeric_limits<int64_t>::max(), 29> >(p);
  rl::simulate<testWaitAndPop<int, lt, 4, 1, std::numeric_limits<int64_t>::max(), 29> >(p);
  rl::simulate<testWaitAndPop<int, lt, 1, 2, std::numeric_limits<int64_t>::max(), 29> >(p);
  rl::simulate<testWaitAndPop<int, lt, 1, 3, std::numeric_limits<int64_t>::max(), 29> >(p);
  rl::simulate<testWaitAndPop<int, lt, 1, 4, std::numeric_limits<int64_t>::max(), 29> >(p);
  rl::simulate<testWaitAndPop<int, lt, 2, 2, std::numeric_limits<int64_t>::max(), 29> >(p);
  rl::simulate<testWaitAndPop<int, lt, 2, 3, std::numeric_limits<int64_t>::max(), 29> >(p);
  rl::simulate<testWaitAndPop<int, lt, 3, 2, std::numeric_limits<int64_t>::max(), 29> >(p);


  rl::simulate<testTryPop<int, lt, 1, 1, std::numeric_limits<int64_t>::max(), 29> >(p);
  rl::simulate<testTryPop<int, lt, 2, 1, std::numeric_limits<int64_t>::max(), 29> >(p);
  rl::simulate<testTryPop<int, lt, 3, 1, std::numeric_limits<int64_t>::max(), 29> >(p);
  rl::simulate<testTryPop<int, lt, 4, 1, std::numeric_limits<int64_t>::max(), 29> >(p);
  rl::simulate<testTryPop<int, lt, 1, 2, std::numeric_limits<int64_t>::max(), 29> >(p);
  rl::simulate<testTryPop<int, lt, 1, 3, std::numeric_limits<int64_t>::max(), 29> >(p);
  rl::simulate<testTryPop<int, lt, 1, 4, std::numeric_limits<int64_t>::max(), 29> >(p);
  rl::simulate<testTryPop<int, lt, 2, 2, std::numeric_limits<int64_t>::max(), 29> >(p);
  rl::simulate<testTryPop<int, lt, 2, 3, std::numeric_limits<int64_t>::max(), 29> >(p);
  rl::simulate<testTryPop<int, lt, 3, 2, std::numeric_limits<int64_t>::max(), 29> >(p);


};
