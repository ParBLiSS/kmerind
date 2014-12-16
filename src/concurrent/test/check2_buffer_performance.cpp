/**
 * @file		check_buffer.cpp
 * @ingroup
 * @author	tpan
 * @brief
 * @details
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */


#include "concurrent/buffer.hpp"


#include "omp.h"
#include <atomic>
#include <cassert>
#include <cstdio>
#include <chrono>
#include <iterator>  // for ostream_iterator
#include <iostream>   // for cout
#include <sstream>
#include <xmmintrin.h>

#include <utils/test_utils.hpp>



template<bliss::concurrent::LockType TS, int64_t CAP, int NumThreads = 1>
void appendTimed(const int iterations) {

  printf("PROFILING: %d threads, locktype %d buffer append with %ld elements and %d iterations\n", NumThreads, static_cast<int>(TS), CAP, iterations);
  bliss::io::Buffer<TS, CAP> buf;

  std::chrono::high_resolution_clock::time_point t1, t2;
  std::chrono::duration<double> time_span;

  int success;
  int failure;
  int i;
  int count = CAP / sizeof(int);

  for (int k = 0; k < iterations; ++k) {

      buf.clear_and_unblock_writes();
      t1 = std::chrono::high_resolution_clock::now();

      i = 0;
      success = 0;
      failure = 0;


#pragma omp parallel for num_threads(NumThreads) default(none) shared(buf, count) private(i) reduction(+:success, failure)
      for (i = 0; i < count + NumThreads; ++i) {
        int data = static_cast<int>(i);
        int result = buf.append(&data, sizeof(int));

        if (result & 0x1) ++success;
        else ++failure;

        // not swapping.  this is run with 1000000 entries.
      }

      assert(success == count && failure == NumThreads);

      t2 = std::chrono::high_resolution_clock::now();

      time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
      printf("Append %d (success) and %d (failure) elements duration = %f\n", success, failure, time_span.count());

  }

}


int main(int argc, char** argv) {


#if defined( BLISS_MUTEX)
  constexpr bliss::concurrent::LockType lt = bliss::concurrent::LockType::MUTEX;
#elif defined(BLISS_SPINLOCK)
  constexpr bliss::concurrent::LockType lt = bliss::concurrent::LockType::SPINLOCK;
#elif defined(BLISS_LOCKFREE)
  constexpr bliss::concurrent::LockType lt = bliss::concurrent::LockType::LOCKFREE;
#else //if defined(BLISS_LOCKFREE2)
  constexpr bliss::concurrent::LockType lt = bliss::concurrent::LockType::LOCKFREE2;
#endif

    // timing tests

    ////////////// timing.  the insert before this is to warm up.
    appendTimed<bliss::concurrent::LockType::NONE, 10000000, 1>( 10);

    appendTimed<lt, 10000000, 1>(10);
    appendTimed<lt, 10000000, 2>(10);
    appendTimed<lt, 10000000, 3>(10);
    appendTimed<lt, 10000000, 4>(10);
    appendTimed<lt, 10000000, 5>(10);
    appendTimed<lt, 10000000, 6>(10);
    appendTimed<lt, 10000000, 7>(10);
    appendTimed<lt, 10000000, 8>(10);
}
