/**
 * @file    testconcurrentqueue.cpp
 * @ingroup
 * @author  tpan
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
    int v;
    assert(q.try_dequeue(v));
    printf("tid %d dequeued %d\n", omp_get_thread_num(), v);
  }

#pragma omp parallel num_threads(4) shared (q)
  {
    q.enqueue(omp_get_thread_num());

#pragma omp barrier
    int v;
    assert(q.try_dequeue(v));
    printf("tid %d dequeued %d\n", omp_get_thread_num(), v);
  }




#pragma omp parallel num_threads(MAX_THREADS) shared(q)
    {
      for (volatile int i = 0; i != 4096; ++i) {
        continue;
      }
      q.enqueue(omp_get_thread_num());

      printf("tid %d enqueue.\n", omp_get_thread_num());
      for (volatile int i = 0; i != 4096; ++i) {
        continue;
      }

#pragma omp barrier
    }


  std::vector<bool> seenIds(MAX_THREADS, false);
  int v;
  for (std::size_t i = 0; i != MAX_THREADS; ++i) {
    assert(q.try_dequeue(v));
    if (seenIds[v]) printf("already seen %d\n", v);
    else printf("haven't seen %d\n", v);
    seenIds[v] = true;
  }
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
    int item;
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

    int item;
    int i = 0;
    while (q.try_dequeue(item)) {
      assert(item == i);
      ++i;
    }
    assert(i == 4);

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
      for (std::size_t i = 0; i != threads.size(); ++i) {
        assert(q.try_dequeue(item));
        assert(!seenIds[item]);
        seenIds[item] = true;
      }
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
        int item;
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



int main (int argc, char** argv) {

   test3();
   test4();
   test1();
   test2();
}
