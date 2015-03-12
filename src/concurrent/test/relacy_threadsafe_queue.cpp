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
template<typename T, bliss::concurrent::LockType LT, int nProducers,
    int nConsumers, int64_t capacity, int64_t entries>
struct testWaitAndPush : rl::test_suite<
    testWaitAndPush<T, LT, nProducers, nConsumers, capacity, entries>,
    nProducers + nConsumers>
{

    bliss::concurrent::ThreadSafeQueue<T, LT> queue;
    std::array<int64_t, nProducers + nConsumers> lsuccess;
    std::array<int64_t, nProducers + nConsumers> lfail;
    static constexpr int64_t expected = (
        capacity > ((entries / nProducers) * nProducers) ?
            ((entries / nProducers) * nProducers) : capacity);

    testWaitAndPush()
        : queue(capacity)
    {
    }
    ;

    void before()
    {
      for (int i = 0; i < nProducers + nConsumers; ++i)
      {
        lsuccess[i] = 0;
        lfail[i] = 0;
      }
    }

    void thread(unsigned thread_index)
    {
      THREAD_EXIT_NOTIFIER_START;

      if (thread_index < nProducers)
      {
        for (int64_t i = 0; i < entries / nProducers; ++i)
        {
          if (queue.waitAndPush(T(i * nProducers + thread_index)).first)
            ++lsuccess[thread_index];
          else
            ++lfail[thread_index];
        }
      }
      else
      { //if (thread_index == (nProducers + nConsumers - 1)){
        while (queue.getSize() < expected)
          ;
        queue.disablePush();
      }

      THREAD_EXIT_NOTIFIER_END;
    }

    void after()
    {
      int64_t success = 0;
      int64_t fail = 0;
      for (int i = 0; i < nProducers + nConsumers; ++i)
      {
        success += lsuccess[i];
        fail += lfail[i];
      }

      if (success != expected)
        INFOF(
            "Test WaitAndPush: LT %d nProd %d nCons %d cap %ld elems %ld.  capacity %ld, expected %lu  actual %ld, failed %ld",
            LT, nProducers, nConsumers, queue.getCapacity(), entries, capacity,
            expected, success, fail);
      assert(success == expected);

    }

    void invariant()
    {

    }

};

template<typename T, bliss::concurrent::LockType LT, int nProducers,
    int nConsumers, int64_t capacity, int64_t entries>
struct testTryPush : rl::test_suite<
    testTryPush<T, LT, nProducers, nConsumers, capacity, entries>,
    nProducers + nConsumers>
{

    bliss::concurrent::ThreadSafeQueue<T, LT> queue;
    std::array<int64_t, nProducers + nConsumers> lsuccess;
    std::array<int64_t, nProducers + nConsumers> lfail;
    static constexpr int64_t expected = (
        capacity > ((entries / nProducers) * nProducers) ?
            ((entries / nProducers) * nProducers) : capacity);

    testTryPush()
        : queue(capacity)
    {
    }
    ;

    void before()
    {
      for (int i = 0; i < nProducers + nConsumers; ++i)
      {
        lsuccess[i] = 0;
        lfail[i] = 0;
      }
    }

    void thread(unsigned thread_index)
    {
      THREAD_EXIT_NOTIFIER_START;

      if (thread_index < nProducers)
      {
        for (int64_t i = 0; i < entries / nProducers; ++i)
        {
          if (queue.tryPush(T(i * nProducers + thread_index)).first)
            ++lsuccess[thread_index];
          else
          {
            --i; // retry.
            ++lfail[thread_index];
          }
        }
      }
      else
      { //if (thread_index == (nProducers + nConsumers - 1)){
        while (queue.getSize() < expected)
          ;
        queue.disablePush();
      }

      THREAD_EXIT_NOTIFIER_END;
    }

    void after()
    {
      int64_t success = 0;
      int64_t fail = 0;
      for (int i = 0; i < nProducers + nConsumers; ++i)
      {
        success += lsuccess[i];
        fail += lfail[i];
      }

      if (success != expected)
        INFOF(
            "Test TryPush: LT %d nProd %d nCons %d cap %ld elems %ld.  capacity %ld, expected %lu  actual %ld, failed %ld",
            LT, nProducers, nConsumers, queue.getCapacity(), entries, capacity,
            expected, success, fail);
      assert(success == expected);

    }

    void invariant()
    {

    }

};

template<typename T, bliss::concurrent::LockType LT, int nProducers,
    int nConsumers, int64_t capacity, int64_t entries>
struct testWaitAndPop : rl::test_suite<
    testWaitAndPop<T, LT, nProducers, nConsumers, capacity, entries>,
    nProducers + nConsumers>
{

    bliss::concurrent::ThreadSafeQueue<T, LT> queue;
    std::array<int64_t, nProducers + nConsumers> lsuccess;
    std::array<int64_t, nProducers + nConsumers> lfail;
    static constexpr int64_t expected =
        (capacity > entries ? entries : capacity);

    testWaitAndPop()
        : queue(capacity)
    {
    }
    ;

    void before()
    {
      for (int i = 0; i < nProducers + nConsumers; ++i)
      {
        lsuccess[i] = 0;
        lfail[i] = 0;
      }

      for (int64_t i = 0; i < expected; ++i)
      {
        queue.waitAndPush(T(i));
      }
      queue.disablePush();
    }

    void thread(unsigned thread_index)
    {
      THREAD_EXIT_NOTIFIER_START;

      if (thread_index >= nProducers)
      {
        while (queue.canPop())
        {
          if (queue.waitAndPop().first)
            ++lsuccess[thread_index];
          else
            ++lfail[thread_index];
        }
      }

      THREAD_EXIT_NOTIFIER_END;
    }

    void after()
    {
      int64_t success = 0;
      int64_t fail = 0;
      for (int i = 0; i < nProducers + nConsumers; ++i)
      {
        success += lsuccess[i];
        fail += lfail[i];
      }

      if (success != expected)
        INFOF(
            "Test WaitAndPop: LT %d nProd %d nCons %d cap %ld elems %ld.  capacity %ld, expected %lu  actual %ld, failed %ld",
            LT, nProducers, nConsumers, queue.getCapacity(), entries, capacity,
            expected, success, fail);
      assert(success == expected);

    }

    void invariant()
    {

    }

};

template<typename T, bliss::concurrent::LockType LT, int nProducers,
    int nConsumers, int64_t capacity, int64_t entries>
struct testTryPop : rl::test_suite<
    testTryPop<T, LT, nProducers, nConsumers, capacity, entries>,
    nProducers + nConsumers>
{

    bliss::concurrent::ThreadSafeQueue<T, LT> queue;
    std::array<int64_t, nProducers + nConsumers> lsuccess;
    std::array<int64_t, nProducers + nConsumers> lfail;
    static constexpr int64_t expected =
        (capacity > entries ? entries : capacity);

    testTryPop()
        : queue(capacity)
    {
    }
    ;

    void before()
    {
      for (int i = 0; i < nProducers + nConsumers; ++i)
      {
        lsuccess[i] = 0;
        lfail[i] = 0;
      }

      for (int64_t i = 0; i < capacity; ++i)
      {
        queue.waitAndPush(T(i));
      }
      queue.disablePush();
    }

    void thread(unsigned thread_index)
    {
      THREAD_EXIT_NOTIFIER_START;

      if (thread_index >= nProducers)
      {
        while (queue.canPop())
        {
          if (queue.tryPop().first)
            ++lsuccess[thread_index];
          else
          {
            ++lfail[thread_index];
          }
        }
      }

      THREAD_EXIT_NOTIFIER_END;
    }

    void after()
    {
      int64_t success = 0;
      int64_t fail = 0;
      for (int i = 0; i < nProducers + nConsumers; ++i)
      {
        success += lsuccess[i];
        fail += lfail[i];
      }

      if (success != expected)
        INFOF(
            "Test TryPop: LT %d nProd %d nCons %d cap %ld elems %ld.  capacity %ld, expected %lu  actual %ld, failed %ld",
            LT, nProducers, nConsumers, queue.getCapacity(), entries, capacity,
            expected, success, fail);
      assert(success == expected);

    }

    void invariant()
    {

    }

};

template<typename T, bliss::concurrent::LockType LT, int nProducers,
    int nConsumers, int64_t capacity, int64_t entries>
struct testTryPushAndTryPop : rl::test_suite<
    testTryPushAndTryPop<T, LT, nProducers, nConsumers, capacity, entries>,
    nProducers + nConsumers>
{

    bliss::concurrent::ThreadSafeQueue<T, LT> queue;
    std::array<int64_t, nProducers + nConsumers> lsuccess;
    std::array<int64_t, nProducers + nConsumers> lfail;
    static constexpr int64_t expectedMin = (
        capacity > ((entries / nProducers) * nProducers) ?
            ((entries / nProducers) * nProducers) : capacity);
    static constexpr int64_t expectedMax = entries;

    std::atomic<int> doneProducers;

    testTryPushAndTryPop()
        : queue(capacity)
    {
    }
    ;

    void before()
    {
      for (int i = 0; i < nProducers + nConsumers; ++i)
      {
        lsuccess[i] = 0;
        lfail[i] = 0;
      }
      doneProducers.store(0, std::memory_order_relaxed);
    }

    void thread(unsigned thread_index)
    {
      THREAD_EXIT_NOTIFIER_START;

      int x;

      if (thread_index < nProducers)
      {
        for (int i = 0; i < entries / nProducers; ++i)
        {
          if (queue.tryPush(T(i * nProducers + thread_index)).first)
            ++lsuccess[thread_index];
          else
          {
            --i; // retry.
            ++lfail[thread_index];
          }
        }
        x = doneProducers.fetch_add(1, std::memory_order_relaxed);
        if (x + 1 == nProducers)
        {
          queue.disablePush();
        }

      }
      else
      {
        while (queue.canPop())
        {
          if (queue.tryPop().first)
            ++lsuccess[thread_index];
          else
          {
            ++lfail[thread_index];
          }
        }
      }

      THREAD_EXIT_NOTIFIER_END;
    }

    void after()
    {
      int64_t pushsuccess = 0;
      int64_t pushfail = 0;
      for (int i = 0; i < nProducers; ++i)
      {
        pushsuccess += lsuccess[i];
        pushfail += lfail[i];
      }

      int64_t popsuccess = 0;
      int64_t popfail = 0;
      for (int i = 0; i < nConsumers; ++i)
      {
        popsuccess += lsuccess[i + nProducers];
        popfail += lfail[i + nProducers];
      }

      if (pushsuccess < expectedMin || pushsuccess > expectedMax
          || pushsuccess != popsuccess)
        INFOF(
            "Test TryPushAndTryPop: LT %d nProd %d nCons %d cap %ld elems %ld.  capacity %ld, expected [%lu - %lu]  push(Y/N) %ld/%ld, pop(Y/N) %ld/%ld",
            LT, nProducers, nConsumers, queue.getCapacity(), entries, capacity,
            expectedMin, expectedMax, pushsuccess, pushfail, popsuccess,
            popfail);
      assert(pushsuccess == popsuccess);
      assert(pushsuccess >= expectedMin && pushsuccess <= expectedMax);

    }

    void invariant()
    {

    }

};

template<typename T, bliss::concurrent::LockType LT, int nProducers,
    int nConsumers, int64_t capacity, int64_t entries>
struct testWaitPushAndTryPop : rl::test_suite<
    testWaitPushAndTryPop<T, LT, nProducers, nConsumers, capacity, entries>,
    nProducers + nConsumers>
{

    bliss::concurrent::ThreadSafeQueue<T, LT> queue;
    std::array<int64_t, nProducers + nConsumers> lsuccess;
    std::array<int64_t, nProducers + nConsumers> lfail;
    static constexpr int64_t expectedMin = (
        capacity > ((entries / nProducers) * nProducers) ?
            ((entries / nProducers) * nProducers) : capacity);
    static constexpr int64_t expectedMax = entries;

    std::atomic<int> doneProducers;

    testWaitPushAndTryPop()
        : queue(capacity)
    {
    }
    ;

    void before()
    {
      for (int i = 0; i < nProducers + nConsumers; ++i)
      {
        lsuccess[i] = 0;
        lfail[i] = 0;
      }
      doneProducers.store(0, std::memory_order_relaxed);
    }

    void thread(unsigned thread_index)
    {
      THREAD_EXIT_NOTIFIER_START;

      int x;

      if (thread_index < nProducers)
      {
        for (int i = 0; i < entries / nProducers; ++i)
        {
          if (queue.waitAndPush(T(i * nProducers + thread_index)).first)
            ++lsuccess[thread_index];
          else
          {
            ++lfail[thread_index];
          }
        }
        x = doneProducers.fetch_add(1, std::memory_order_relaxed);
        if (x + 1 == nProducers)
        {
          queue.disablePush();
        }

      }
      else
      {
        while (queue.canPop())
        {
          if (queue.tryPop().first)
            ++lsuccess[thread_index];
          else
          {
            ++lfail[thread_index];
          }
        }
      }

      THREAD_EXIT_NOTIFIER_END;
    }

    void after()
    {
      int64_t pushsuccess = 0;
      int64_t pushfail = 0;
      for (int i = 0; i < nProducers; ++i)
      {
        pushsuccess += lsuccess[i];
        pushfail += lfail[i];
      }

      int64_t popsuccess = 0;
      int64_t popfail = 0;
      for (int i = 0; i < nConsumers; ++i)
      {
        popsuccess += lsuccess[i + nProducers];
        popfail += lfail[i + nProducers];
      }

      if (pushsuccess < expectedMin || pushsuccess > expectedMax
          || pushsuccess != popsuccess || pushfail > nProducers)
        INFOF(
            "Test WaitPushAndTryPop: LT %d nProd %d nCons %d cap %ld elems %ld.  capacity %ld, expected [%lu - %lu]  push(Y/N) %ld/%ld, pop(Y/N) %ld/%ld",
            LT, nProducers, nConsumers, queue.getCapacity(), entries, capacity,
            expectedMin, expectedMax, pushsuccess, pushfail, popsuccess,
            popfail);
      assert(pushsuccess == popsuccess);
      //assert(pushfail <= nProducers);
      assert(pushsuccess >= expectedMin && pushsuccess <= expectedMax);

    }

    void invariant()
    {

    }

};

template<typename T, bliss::concurrent::LockType LT, int nProducers,
    int nConsumers, int64_t capacity, int64_t entries>
struct testTryPushAndWaitPop : rl::test_suite<
    testTryPushAndWaitPop<T, LT, nProducers, nConsumers, capacity, entries>,
    nProducers + nConsumers>
{

    bliss::concurrent::ThreadSafeQueue<T, LT> queue;
    std::array<int64_t, nProducers + nConsumers> lsuccess;
    std::array<int64_t, nProducers + nConsumers> lfail;
    static constexpr int64_t expectedMin = (
        capacity > ((entries / nProducers) * nProducers) ?
            ((entries / nProducers) * nProducers) : capacity);
    static constexpr int64_t expectedMax = entries;

    std::atomic<int> doneProducers;

    testTryPushAndWaitPop()
        : queue(capacity)
    {
    }
    ;

    void before()
    {
      for (int i = 0; i < nProducers + nConsumers; ++i)
      {
        lsuccess[i] = 0;
        lfail[i] = 0;
      }
      doneProducers.store(0, std::memory_order_relaxed);
    }

    void thread(unsigned thread_index)
    {
      THREAD_EXIT_NOTIFIER_START;

      int x;

      if (thread_index < nProducers)
      {
        for (int i = 0; i < entries / nProducers; ++i)
        {
          if (queue.tryPush(T(i * nProducers + thread_index)).first)
            ++lsuccess[thread_index];
          else
          {
            --i; // retry.
            ++lfail[thread_index];
          }
        }
        x = doneProducers.fetch_add(1, std::memory_order_relaxed);
        if (x + 1 == nProducers)
        {
          queue.disablePush();
        }

      }
      else
      {
        while (queue.canPop())
        {
          if (queue.waitAndPop().first)
            ++lsuccess[thread_index];
          else
          {
            ++lfail[thread_index];
          }
        }
      }

      THREAD_EXIT_NOTIFIER_END;
    }

    void after()
    {
      int64_t pushsuccess = 0;
      int64_t pushfail = 0;
      for (int i = 0; i < nProducers; ++i)
      {
        pushsuccess += lsuccess[i];
        pushfail += lfail[i];
      }

      int64_t popsuccess = 0;
      int64_t popfail = 0;
      for (int i = 0; i < nConsumers; ++i)
      {
        popsuccess += lsuccess[i + nProducers];
        popfail += lfail[i + nProducers];
      }

      if (pushsuccess < expectedMin || pushsuccess > expectedMax
          || pushsuccess != popsuccess || popfail > nConsumers)

        INFOF(
            "Test TryPushAndWaitPop: LT %d nProd %d nCons %d cap %ld elems %ld.  capacity %ld, expected [%lu - %lu]  push(Y/N) %ld/%ld, pop(Y/N) %ld/%ld",
            LT, nProducers, nConsumers, queue.getCapacity(), entries, capacity,
            expectedMin, expectedMax, pushsuccess, pushfail, popsuccess,
            popfail);
      assert(pushsuccess == popsuccess);
      //assert(popfail <= nConsumers);
      assert(pushsuccess >= expectedMin && pushsuccess <= expectedMax);

    }

    void invariant()
    {

    }

};

template<typename T, bliss::concurrent::LockType LT, int nProducers,
    int nConsumers, int64_t capacity, int64_t entries>
struct testWaitPushAndWaitPop : rl::test_suite<
    testWaitPushAndWaitPop<T, LT, nProducers, nConsumers, capacity, entries>,
    nProducers + nConsumers>
{

    bliss::concurrent::ThreadSafeQueue<T, LT> queue;
    std::array<int64_t, nProducers + nConsumers> lsuccess;
    std::array<int64_t, nProducers + nConsumers> lfail;
    static constexpr int64_t expectedMin = (
        capacity > ((entries / nProducers) * nProducers) ?
            ((entries / nProducers) * nProducers) : capacity);
    static constexpr int64_t expectedMax = entries;

    std::atomic<int> doneProducers;

    testWaitPushAndWaitPop()
        : queue(capacity)
    {
    }
    ;

    void before()
    {
      for (int i = 0; i < nProducers + nConsumers; ++i)
      {
        lsuccess[i] = 0;
        lfail[i] = 0;
      }
      doneProducers.store(0, std::memory_order_relaxed);
    }

    void thread(unsigned thread_index)
    {
      THREAD_EXIT_NOTIFIER_START;

      int x;

      if (thread_index < nProducers)
      {
        for (int i = 0; i < entries / nProducers; ++i)
        {
          if (queue.waitAndPush(T(i * nProducers + thread_index)).first)
            ++lsuccess[thread_index];
          else
          {
            ++lfail[thread_index];
          }
        }
        x = doneProducers.fetch_add(1, std::memory_order_relaxed);
        if (x + 1 == nProducers)
        {
          queue.disablePush();
        }

      }
      else
      {
        while (queue.canPop())
        {
          if (queue.waitAndPop().first)
            ++lsuccess[thread_index];
          else
          {
            ++lfail[thread_index];
          }
        }
      }

      THREAD_EXIT_NOTIFIER_END;
    }

    void after()
    {
      int64_t pushsuccess = 0;
      int64_t pushfail = 0;
      for (int i = 0; i < nProducers; ++i)
      {
        pushsuccess += lsuccess[i];
        pushfail += lfail[i];
      }

      int64_t popsuccess = 0;
      int64_t popfail = 0;
      for (int i = 0; i < nConsumers; ++i)
      {
        popsuccess += lsuccess[i + nProducers];
        popfail += lfail[i + nProducers];
      }

      if (pushsuccess < expectedMin || pushsuccess > expectedMax
          || pushsuccess != popsuccess || pushfail > nProducers
          || popfail > nConsumers)
        INFOF(
            "Test WaitPushAndWaitPop: LT %d nProd %d nCons %d cap %ld elems %ld.  capacity %ld, expected [%lu - %lu]  push(Y/N) %ld/%ld, pop(Y/N) %ld/%ld",
            LT, nProducers, nConsumers, queue.getCapacity(), entries, capacity,
            expectedMin, expectedMax, pushsuccess, pushfail, popsuccess,
            popfail);
      assert(pushsuccess == popsuccess);
      //assert(pushfail <= nProducers);
      //assert(popfail <= nConsumers);
      assert(pushsuccess >= expectedMin && pushsuccess <= expectedMax);

    }

    void invariant()
    {

    }

};

template<bliss::concurrent::LockType LT, int nProducers, int nConsumers>
void simulate(rl::test_params p)
{

  rl::simulate<testWaitAndPush<int, LT, nProducers, nConsumers, 20, 29> >(p);
  rl::simulate<testTryPush<int, LT, nProducers, nConsumers, 20, 29> >(p);
  rl::simulate<testWaitAndPop<int, LT, nProducers, nConsumers, 20, 29> >(p);
  rl::simulate<testTryPop<int, LT, nProducers, nConsumers, 20, 29> >(p);

  rl::simulate<testTryPushAndTryPop<int, LT, nProducers, nConsumers, 20, 29> >(
      p);
  rl::simulate<testWaitPushAndTryPop<int, LT, nProducers, nConsumers, 20, 29> >(
      p);
  rl::simulate<testTryPushAndWaitPop<int, LT, nProducers, nConsumers, 20, 29> >(
      p);
  rl::simulate<testWaitPushAndWaitPop<int, LT, nProducers, nConsumers, 20, 29> >(
      p);

  rl::simulate<
      testWaitAndPush<int, LT, nProducers, nConsumers,
          std::numeric_limits<int64_t>::max(), 29> >(p);
  rl::simulate<
      testTryPush<int, LT, nProducers, nConsumers,
          std::numeric_limits<int64_t>::max(), 29> >(p);
  rl::simulate<
      testWaitAndPop<int, LT, nProducers, nConsumers,
          std::numeric_limits<int64_t>::max(), 29> >(p);
  rl::simulate<
      testTryPop<int, LT, nProducers, nConsumers,
          std::numeric_limits<int64_t>::max(), 29> >(p);

  rl::simulate<
      testTryPushAndTryPop<int, LT, nProducers, nConsumers,
          std::numeric_limits<int64_t>::max(), 29> >(p);
  rl::simulate<
      testWaitPushAndTryPop<int, LT, nProducers, nConsumers,
          std::numeric_limits<int64_t>::max(), 29> >(p);
  rl::simulate<
      testTryPushAndWaitPop<int, LT, nProducers, nConsumers,
          std::numeric_limits<int64_t>::max(), 29> >(p);
  rl::simulate<
      testWaitPushAndWaitPop<int, LT, nProducers, nConsumers,
          std::numeric_limits<int64_t>::max(), 29> >(p);
}

int main(int argc, char** argv)
{

  int iterations = 1000;
  if (argc > 1)
  {
    iterations = atoi(argv[1]);
  }

  rl::test_params p;
  p.search_type = rl::sched_random;
  p.iteration_count = iterations;
  p.execution_depth_limit = 4000;

#if defined(BLISS_MUTEX)
  constexpr bliss::concurrent::LockType lt = bliss::concurrent::LockType::MUTEX;

#elif defined(BLISS_SPINLOCK)
  constexpr bliss::concurrent::LockType lt = bliss::concurrent::LockType::SPINLOCK;

#else
  constexpr bliss::concurrent::LockType lt =
      bliss::concurrent::LockType::LOCKFREE;

#endif

  simulate<lt, 1, 1>(p);
  simulate<lt, 2, 1>(p);
  simulate<lt, 3, 1>(p);
  simulate<lt, 1, 2>(p);
  simulate<lt, 1, 3>(p);
  simulate<lt, 2, 2>(p);
  simulate<lt, 2, 3>(p);
  simulate<lt, 3, 2>(p);

}
;
