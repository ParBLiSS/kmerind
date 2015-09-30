/**
 * @file		check_buffer.cpp
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
#include "concurrent/buffer.hpp"


#include "omp.h"
#include <atomic>
#include <cassert>
#include <chrono>
#include <iterator>  // for ostream_iterator
#include <iostream>   // for cout
#include <sstream>
#include <xmmintrin.h>

#include "utils/iterator_test_utils.hpp"

#include "utils/logging.h"


template<bliss::concurrent::LockType TS, size_t CAP, int NumThreads = 1>
void appendTimed(const int iterations) {

  INFOF("PROFILING: %d threads, locktype %d buffer append with %ld elements and %d iterations", NumThreads, static_cast<int>(TS), CAP, iterations);
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


#pragma omp parallel for num_threads(NumThreads) OMP_SHARE_DEFAULT shared(buf, count) private(i) reduction(+:success, failure)
      for (i = 0; i < count + NumThreads; ++i) {
        int data = static_cast<int>(i);
        int result = buf.append(&data, 1);

        if (result & 0x1) ++success;
        else ++failure;

        // not swapping.  this is run with 1000000 entries.
      }

      assert(success == count && failure == NumThreads);

      t2 = std::chrono::high_resolution_clock::now();

      time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
      INFOF("Append %d (success) and %d (failure) elements duration = %f", success, failure, time_span.count());

  }

}


int main(int argc, char** argv) {


#if defined( BLISS_MUTEX)
  constexpr bliss::concurrent::LockType lt = bliss::concurrent::LockType::MUTEX;
#elif defined(BLISS_SPINLOCK)
  constexpr bliss::concurrent::LockType lt = bliss::concurrent::LockType::SPINLOCK;
//#elif defined(BLISS_LOCKFREE)
//  constexpr bliss::concurrent::LockType lt = bliss::concurrent::LockType::LOCKFREE;
#else //if defined(BLISS_LOCKFREE)
  constexpr bliss::concurrent::LockType lt = bliss::concurrent::LockType::LOCKFREE;
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
