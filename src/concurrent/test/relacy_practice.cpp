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

int main()
{
    rl::simulate<race_test>();
}

