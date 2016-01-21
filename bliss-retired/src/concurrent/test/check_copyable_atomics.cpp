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
#include "bliss-config.hpp"
#include <vector>
#include <unordered_map>
#include "concurrent/copyable_atomic.hpp"

#include <omp.h>

#include "utils/logging.h"

using namespace bliss::concurrent;


template<int nThreads>
struct testTSAN
{
    int i;

    void before() {
      i = 0;
    }

    void thread(unsigned thread_index) {
      ++i;
    }


    void after() {
      if (i != nThreads) DEBUGF("i =  %d", i);
    }

};



// template parameter '2' is number of threads
template<typename T, int nThreads, int nelems>
struct copyableAtomicsThreadLocalTest
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

};



// template parameter '2' is number of threads
template<typename T, int nThreads, int nelems>
struct copyableAtomicsTest
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
      int t = thread_index;
      for (int i = 0; i < nelems; ++i) {
          if (i % nThreads >= t) atoms.at(i).fetch_add(1, std::memory_order_relaxed);
      }
    }

    // executed in single thread after main thread function
    void after()
    {
      bool success = true;
      bool elsucc = true;
      for (int i = 0; i < nelems; ++i) {
        elsucc = atoms.at(i).load(std::memory_order_relaxed) == (i % nThreads) + 1;
        success &= elsucc;
        if (!elsucc) DEBUGF("atoms[%d] = %d, expect %d", i, atoms.at(i).load(std::memory_order_relaxed), (i%nThreads) + 1  );
      }

      if (!success) ERRORF("incorrect result for copyableAtomicTest");

    }

};


template<int nThreads>
void simulate(int iterations) {

  for (int i = 0; i < iterations; ++i) {

//    testTSAN<nThreads> test;
//    test.before();
//#pragma omp parallel OMP_SHARE_DEFAULT num_threads(nThreads) shared(test)
//    {
//      test.thread(omp_get_thread_num());
//    }
//    test.after();

    copyableAtomicsThreadLocalTest<int, nThreads, 100>  test1;
    test1.before();

#pragma omp parallel OMP_SHARE_DEFAULT num_threads(nThreads) shared(test1)
    {
      test1.thread(omp_get_thread_num());
    }

    test1.after();


    copyableAtomicsTest<int, nThreads, 100>  test2;
    test2.before();

#pragma omp parallel OMP_SHARE_DEFAULT num_threads(nThreads) shared(test2)
    {
      test2.thread(omp_get_thread_num());
    }

    test2.after();


  }

}

int main(int argc, char** argv) {

  int p = 1000;
  if (argc > 1)
  {
    p = atoi(argv[1]);
  }



  simulate<1>(p);
  simulate<2>(p);
  simulate<3>(p);
  simulate<4>(p);
  simulate<8>(p);
  simulate<16>(p);

}


