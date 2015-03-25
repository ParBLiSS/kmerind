/**
 * @file    testconcurrentqueue.cpp
 * @ingroup
 * @author  Tony Pan <tpan7@gatech.edu>
 * @brief
 * @details
 *
 * Copyright (c) 2015 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */



#include "concurrentqueue/concurrentqueue.h"
#include <cassert>
#include <omp.h>

#include <thread>
#include "utils/logging.h"
#include "utils/iterator_test_utils.hpp"


using namespace moodycamel;

void test1() {

  // Implicit
  const int MAX_THREADS = 48;
  ConcurrentQueue<int> q(4096 * (MAX_THREADS + 1));

#pragma omp parallel num_threads(4) shared (q)
  {
    q.enqueue(omp_get_thread_num());
  }

#pragma omp parallel num_threads(4) shared (q)
  {
    int v = -1;
    bool r = q.try_dequeue(v);
    assert(r);
    INFOF("tid %d dequeued %d.  %s", omp_get_thread_num(), v, (r? "success" : "failure"));
  }

#pragma omp parallel num_threads(4) shared (q)
  {
    q.enqueue(omp_get_thread_num());

#pragma omp barrier
    int v = -1;
    bool r = q.try_dequeue(v);
    assert(r);
    INFOF("tid %d dequeued %d. %s", omp_get_thread_num(), v, (r? "success" : "failure"));
  }




#pragma omp parallel num_threads(MAX_THREADS) shared(q)
    {
      for (volatile int i = 0; i != 4096; ++i) {
        continue;
      }
      q.enqueue(omp_get_thread_num());

      INFOF("tid %d enqueue.", omp_get_thread_num());
      for (volatile int i = 0; i != 4096; ++i) {
        continue;
      }

#pragma omp barrier
    }


  std::vector<bool> seenIds(MAX_THREADS, false);
  int v = -1;
  bool r = true;
  for (std::size_t i = 0; i != MAX_THREADS; ++i) {
	  r &= q.try_dequeue(v);
    assert(r);
    if (seenIds[v]) INFOF("already seen %d", v);
    else INFOF("haven't seen %d", v);
    seenIds[v] = true;
  }
  if (!r) INFOF("there was a failed dequeue.");
  for (std::size_t i = 0; i != MAX_THREADS; ++i) {
    assert(seenIds[i]);
  }
}


void test2() {

  // Test many threads and implicit queues being created and destroyed concurrently
  int nThreads = 32;
  std::vector<bool> success(nThreads, true);
#pragma omp parallel num_threads(nThreads) shared(success)
  {
    for (int i = 0; i != 5; ++i) {
      ConcurrentQueue<int> q(1);
      q.enqueue(i);
    }

    ConcurrentQueue<int> q(15);
    for (int i = 0; i != 100; ++i) {
      q.enqueue(i);
    }
    int item = -1;
    for (int i = 0; i != 100; ++i) {
      if (!q.try_dequeue(item) || item != i) {
        success[omp_get_thread_num()] = false;
      }
    }
    if (q.size_approx() != 0) {
      success[omp_get_thread_num()] = false;
    }
#pragma omp barrier
  }
  for (int tid = 0; tid != nThreads; ++tid) {
    assert(success[tid]);
  }
}

// this one is not openmp based.
void test3() {

    // Implicit
    const int MAX_THREADS = 48;
    ConcurrentQueue<int> q(4096 * (MAX_THREADS + 1));

    std::thread t0([&]() { q.enqueue(0); });
    t0.join();

    std::thread t1([&]() { q.enqueue(1); });
    t1.join();

    std::thread t2([&]() { q.enqueue(2); });
    t2.join();

    q.enqueue(3);

    int item = -1;
    int i = 0;
    std::vector<int> vals;
    while (q.try_dequeue(item)) {
      vals.push_back(item);
    }
    assert(vals.size() == 4);
    std::sort(vals.begin(), vals.end());
    int s = vals.size();
    for (; i < s; ++i) {
      assert(vals[i] == i);
    }

    std::vector<std::thread> threads(MAX_THREADS);
    for (int rep = 0; rep != 2; ++rep) {
      for (std::size_t tid = 0; tid != threads.size(); ++tid) {
        threads[tid] = std::thread([&](std::size_t tid) {
          for (volatile int i = 0; i != 4096; ++i) {
            continue;
          }
          q.enqueue((int)tid);
          for (volatile int i = 0; i != 4096; ++i) {
            continue;
          }
        }, tid);
      }
      for (std::size_t tid = 0; tid != threads.size(); ++tid) {
        threads[tid].join();
      }
      std::vector<bool> seenIds(threads.size());
      bool r = true;
      for (std::size_t i = 0; i != threads.size(); ++i) {
    	  r &= q.try_dequeue(item);
        assert(!seenIds[item]);
        seenIds[item] = true;
      }
      if (!r) INFOF("there was a failed dequeue.");
      for (std::size_t i = 0; i != seenIds.size(); ++i) {
        assert(seenIds[i]);
      }
    }
}

void test4 () {

    // Test many threads and implicit queues being created and destroyed concurrently
    std::vector<std::thread> threads(32);
    std::vector<bool> success(threads.size(), true);
    for (std::size_t tid = 0; tid != threads.size(); ++tid) {
      threads[tid] = std::thread([&](std::size_t tid) {
        for (int i = 0; i != 5; ++i) {
          ConcurrentQueue<int> q(1);
          q.enqueue(i);
        }

        ConcurrentQueue<int> q(15);
        for (int i = 0; i != 100; ++i) {
          q.enqueue(i);
        }
        int item = -1;
        for (int i = 0; i != 100; ++i) {
          if (!q.try_dequeue(item) || item != i) {
            success[tid] = false;
          }
        }
        if (q.size_approx() != 0) {
          success[tid] = false;
        }
      }, tid);
    }
    for (std::size_t tid = 0; tid != threads.size(); ++tid) {
      threads[tid].join();
      assert(success[tid]);
    }


}


void test5() {

  int elements = 100000;
  // Test many threads and implicit queues being created and destroyed concurrently
  int nThreads = 4;
  std::vector< int > input(elements);

  for (int i = 0; i < elements; ++i) {
    input[i] = i;
  }

  std::vector<std::vector<int> > outputs(nThreads);
  for (int i = 0; i < nThreads; ++i ) {
    outputs[i].clear();
  }

  std::atomic<bool> finished(false);
  ConcurrentQueue<int> q;

#pragma omp parallel sections default(none) shared(q, input, outputs, nThreads, elements, finished) num_threads(2)
  {
#pragma omp section
    {
#pragma omp parallel default(none) shared(q, input, elements, finished) num_threads(nThreads)
      {
         int tid = omp_get_thread_num();
         int nt = omp_get_num_threads();

         for (int i = tid; i < elements; i += nt) {
           q.enqueue(input[i]);
         }
      }

      finished.store(true, std::memory_order_release);
    }

#pragma omp section
    {
#pragma omp parallel default(none) shared(q, outputs, finished) num_threads(nThreads)
      {
        int tid = omp_get_thread_num();
        int v = 0;
        bool r = false;

        while (!finished.load(std::memory_order_acquire) || (r = q.try_dequeue(v))) {
          if (r) outputs[tid].push_back(v);
        }
      }
    }
  }


  std::vector<int> output;
  for (int i = 0; i < nThreads; ++i) {
    output.insert(output.end(), outputs[i].begin(), outputs[i].end());
  }

  bool result = compareUnorderedSequences(input.begin(), output.begin(), elements);
  INFOF("enqueue/dequeue result is same? %s", (result ? "y" : "n"));
}

void testMove5() {

  int elements = 100000;
  // Test many threads and implicit queues being created and destroyed concurrently
  int nThreads = 4;
  std::vector< int > input(elements);

  for (int i = 0; i < elements; ++i) {
    input[i] = i;
  }

  std::vector<std::vector<int> > outputs(nThreads);
  for (int i = 0; i < nThreads; ++i ) {
    outputs[i].clear();
  }

  std::atomic<bool> finished(false);
  ConcurrentQueue<int> q;

#pragma omp parallel sections default(none) shared(q, input, outputs, nThreads, elements, finished) num_threads(2)
  {
#pragma omp section
    {
#pragma omp parallel default(none) shared(q, input, elements, finished) num_threads(nThreads)
      {
         int tid = omp_get_thread_num();
         int nt = omp_get_num_threads();
         bool res = false;

         for (int i = tid; i < elements; i += nt) {
           int v = input[i];
           res = q.enqueue(std::move(v));
         }
      }

      finished.store(true, std::memory_order_release);
    }

#pragma omp section
    {
#pragma omp parallel default(none) shared(q, outputs, finished) num_threads(nThreads)
      {
        int tid = omp_get_thread_num();
        int v = 0;
        bool r = false;

        while (!finished.load(std::memory_order_acquire) || (r = q.try_dequeue(v))) {
          if (r) outputs[tid].push_back(v);
        }
      }
    }
  }


  std::vector<int> output;
  for (int i = 0; i < nThreads; ++i) {
    output.insert(output.end(), outputs[i].begin(), outputs[i].end());
  }

  bool result = compareUnorderedSequences(input.begin(), output.begin(), elements);
  INFOF("enqueue/dequeue result is same? %s", (result ? "y" : "n"));
}





int main (int argc, char** argv) {


//   test3();
//   test4();
   test1();
   test2();

   test5();
   testMove5();

}
