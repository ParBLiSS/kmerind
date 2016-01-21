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

#include <vector>
#include <unordered_map>
#include "concurrent/copyable_atomic.hpp"

#include "utils/logging.h"

using namespace bliss::concurrent;



// template parameter '2' is number of threads
template<typename T, int nThreads, int nelems>
struct copyableAtomicsThreadLocalTest : rl::test_suite<copyableAtomicsThreadLocalTest<T, nThreads, nelems>, nThreads>
{
  std::unordered_map<int, copyable_atomic<int> > atoms;


    // executed in single thread before main thread function
    void before()
    {
      for (int i = 0; i < nelems; ++i) {
          //atoms.emplace(std::make_pair(i, std::move(copyable_atomic<T>(i))));
        atoms.emplace(std::make_pair(i, copyable_atomic<T>(i)));
          atoms.emplace(i + nelems, copyable_atomic<T>(1));
          atoms[i + 2 * nelems] = copyable_atomic<T>(1);
          atoms[i + 3 * nelems] = copyable_atomic<T>(1);
      }
    }

    // main thread function
    void thread(unsigned thread_index)
    {
      for (int i = thread_index; i < nelems; i+= nThreads) {
          atoms.at(i + 3*nelems).store(i, std::memory_order_relaxed);
          atoms.at(i + nelems) = copyable_atomic<T>(i);  // copy assignment
          //atoms.at(i + 2*nelems) = std::move(copyable_atomic<T>(i));  // move assignment
          atoms.at(i + 2*nelems) = copyable_atomic<T>(i);  // move assignment

          atoms.at(i).fetch_xor(atoms.at(i+ nelems).load(std::memory_order_relaxed), std::memory_order_relaxed);
      }
    }

    // executed in single thread after main thread function
    void after()
    {
      bool success1 = true, success2 = true, success3 = true, success4 = true;

      for (int i = 0; i < nelems; ++i) {
        success1 &= atoms.at(i).load(std::memory_order_relaxed) == 0;
        success2 &= atoms.at(i + nelems).load(std::memory_order_relaxed) == i;
        success3 &= atoms.at(i + 2 * nelems).load(std::memory_order_relaxed) == i;
        success4 &= atoms.at(i + 3 * nelems).load(std::memory_order_relaxed) == i;
      }

      assert(success1);
      assert(success2);
      assert(success3);
      assert(success4);


    }

    // executed in single thread after every 'visible' action in main threads
    // disallowed to modify any state
    void invariant()
    {
    }
};



// template parameter '2' is number of threads
template<typename T, int nThreads, int nelems>
struct copyableAtomicsTest : rl::test_suite<copyableAtomicsTest<T, nThreads, nelems>, nThreads>
{
  std::unordered_map<int, copyable_atomic<int> > atoms;  // unordered map requires copy constructible elems

    // executed in single thread before main thread function
    void before()
    {
      for (int i = 0; i < nelems; ++i) {
          //atoms.emplace(std::make_pair(i, std::move(copyable_atomic<T>(0))));
        atoms.emplace(std::make_pair(i, copyable_atomic<T>(0)));
      }
    }

    // main thread function
    void thread(unsigned thread_index)
    {
      for (int i = nelems-1; i >= thread_index; --i) {
          atoms.at(i).fetch_add(1, std::memory_order_relaxed);
      }
    }

    // executed in single thread after main thread function
    void after()
    {
      bool success = true;

      for (int i = 0; i < nelems; ++i) {
        success &= atoms.at(i).load(std::memory_order_relaxed) == i;
      }

      assert(success);
    }

    // executed in single thread after every 'visible' action in main threads
    // disallowed to modify any state
    void invariant()
    {
    }
};


template<int nThreads>
void simulate(rl::test_params& p) {
  rl::simulate<copyableAtomicsThreadLocalTest<int, nThreads, 100> >(p);
  rl::simulate<copyableAtomicsTest<int, nThreads, 100> >(p);
}

int main(int argc, char** argv) {

  int iterations = 1000;
  if (argc > 1)
  {
    iterations = atoi(argv[1]);
  }

  rl::test_params p;
  p.search_type = rl::sched_random;
  p.iteration_count = iterations;
  p.execution_depth_limit = 4000;


  simulate<1>(p);
  simulate<2>(p);
  simulate<3>(p);
  simulate<4>(p);
  simulate<8>(p);
  simulate<16>(p);

}
