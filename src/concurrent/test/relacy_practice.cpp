/**
 * @file		test_copyable_atomics.cpp
 * @ingroup
 * @author	Tony Pan <tpan7@gatech.edu>
 * @brief
 * @details
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */

//#include <cstdio>  // this causes a long jun to uninitialized stack frame.
//#include "utils/logging.h"

#include "config/relacy_config.hpp"


// template parameter '2' is number of threads
struct race_test : rl::test_suite<race_test, 2>
{
    std::atomic<int> a;
    VAR_T(int) x;
    TLS_T(double) y;

    // executed in single thread before main thread function
    void before()
    {
        VAR(a) = 0;
        VAR(x) = 0;
    }

    // main thread function
    void thread(unsigned thread_index)
    {
        VAR(y) = 0.0;
        if (0 == thread_index)
        {
          VAR(y) = 1.0;
          VAR(x) = 1;
          a.store(1, std::memory_order_relaxed);   // see relacy_config.hpp for important information about variable usage.

        }
        else
        {
          VAR(y) = 2.0;
            if (1 == a.load(std::memory_order_relaxed))
              VAR(x) = 2;
        }
    }

    // executed in single thread after main thread function
    void after()
    {
    }

    // executed in single thread after every 'visible' action in main threads
    // disallowed to modify any state
    void invariant()
    {
    }
};


// template parameter '2' is number of threads
struct norace_test : rl::test_suite<norace_test, 2>
{
    std::atomic<int> a;
    VAR_T(int) x;
    TLS_T(double) y;

    // executed in single thread before main thread function
    void before()
    {
        VAR(a) = 0;
        VAR(x) = 0;
    }

    // main thread function
    void thread(unsigned thread_index)
    {
        VAR(y) = 0.0;
        if (0 == thread_index)
        {
          VAR(y) = 1.0;
          VAR(x) = 1;
          a.store(1, std::memory_order_release);   // see relacy_config.hpp for important information about variable usage.

        }
        else
        {
          VAR(y) = 2.0;
            if (1 == a.load(std::memory_order_acquire)) {
              VAR(x) = 2;
              std::atomic_thread_fence(std::memory_order_release);
            }
        }
    }

    // executed in single thread after main thread function
    void after()
    {
    }

    // executed in single thread after every 'visible' action in main threads
    // disallowed to modify any state
    void invariant()
    {
    }
};


// template parameter '2' is number of threads
template<int nthreads, int elems>
struct condvar_test : rl::test_suite<condvar_test<nthreads, elems>, nthreads>
{
    std::mutex m;
    std::condition_variable cv;
    VAR_T(int) x;
    std::array<int, nthreads> results;

    // executed in single thread before main thread function
    void before()
    {
      VAR(x) = 0;

      for (int i = 0; i < nthreads; ++i) {
        results[i] = 0;
      }
    }

    // main thread function
    void thread(unsigned thread_index)
    {
      std::unique_lock<std::mutex> lock(m, std::defer_lock);
      if (thread_index == 0) {
        for (int i = 0; i < elems; ++i) {
          lock.lock();
          ++VAR(x);
          lock.unlock();
        }
        CV_NOTIFY_ALL(cv);
        for (int i = 0; i < elems; ++i) {
          lock.lock();
          ++VAR(x);
          lock.unlock();
        }
        results[0] = VAR(x);
        CV_NOTIFY_ALL(cv);  // what happens if notify_all is called again?

      } else {
        int v = 0;
        lock.lock();
        while ((v = VAR(x)) < elems) {
          CV_WAIT(cv, lock);
        }
        lock.unlock();
        results[thread_index] = v;
      }
    }

    // executed in single thread after main thread function
    void after()
    {
      assert(VAR(x) == 2*elems);
      assert(results[0] == 2*elems);
      for (int i = 1; i < nthreads; ++i) {
        assert(results[i] >= elems);
      }
    }

    // executed in single thread after every 'visible' action in main threads
    // disallowed to modify any state
    void invariant()
    {
    }
};



int main()
{

    rl::simulate<race_test>();
    rl::simulate<norace_test>();
    rl::simulate<condvar_test<2, 100> >();
    rl::simulate<condvar_test<3, 100> >();
    rl::simulate<condvar_test<4, 100> >();
    rl::simulate<condvar_test<2, 200> >();
    rl::simulate<condvar_test<3, 200> >();
    rl::simulate<condvar_test<4, 200> >();
    rl::simulate<condvar_test<2, 400> >();
    rl::simulate<condvar_test<3, 400> >();
    rl::simulate<condvar_test<4, 400> >();
    rl::simulate<condvar_test<2, 800> >();
    rl::simulate<condvar_test<3, 800> >();
    rl::simulate<condvar_test<4, 800> >();

    rl::test_params p;
    p.search_type = rl::sched_random;
    p.execution_depth_limit = 4000;

    rl::simulate<condvar_test<2, 800> >(p);
    rl::simulate<condvar_test<3, 800> >(p);
    rl::simulate<condvar_test<4, 800> >(p);

}

