/**
 * @file		check_bufferpool.cpp
 * @ingroup
 * @author	tpan
 * @brief
 * @details
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */

#include <unistd.h>  // for usleep

#include "omp.h"
#include <cassert>
#include <chrono>
#include <vector>
#include <cstdlib>   // for rand
#include <atomic>
#include <memory>

#include "utils/test_utils.hpp"

#include "concurrent/buffer.hpp"
#include "concurrent/object_pool.hpp"



template<typename PoolType>
void timeAppendMultipleBuffers(const int NumThreads, const int total_count, bliss::concurrent::LockType poollt, bliss::concurrent::LockType bufferlt, const int64_t buffer_cap) {


//  printf("TESTING: %d threads, pool lock %d buffer lock %d append with %ld bufferSize and %d total counts from unlimited pool\n",
//         NumThreads, poollt, bufferlt, buffer_cap, total_count);


  PoolType pool;

  std::vector<int> gold;
  std::vector<int> stored;

  int data = 12;
  unsigned int result = 0;

  int success = 0;
  int failure = 0;
  int swap = 0;
  int i = 0;

  auto buf_ptr = pool.tryAcquireObject();
  buf_ptr->clear_and_unblock_writes();

#pragma omp parallel for num_threads(NumThreads) default(none) shared(buf_ptr, gold, stored, pool) private(i, data, result) reduction(+:success, failure, swap)
  for (i = 0; i < total_count; ++i) {

    auto sptr = buf_ptr;

    if (sptr) {  // valid ptr

      data = static_cast<int>(i);
      result = sptr->append(&data, sizeof(int));
    } else {  // expired ptr
      result = 0x0;
    }

    if (result & 0x1) {
      ++success;
    }
    else ++failure;

    if (result & 0x2) {
      ++swap;

      // swap in a new one.
      auto new_buf_ptr = pool.tryAcquireObject();
      if (new_buf_ptr) new_buf_ptr->clear_and_unblock_writes();

      sptr = buf_ptr;
      buf_ptr = new_buf_ptr;
#pragma omp flush(buf_ptr)

      // and release - few threads doing this, and full.

      pool.releaseObject(sptr);

    }

  }

  if (buf_ptr) buf_ptr->block_and_flush();

  // compare unordered buffer content.
  pool.releaseObject(buf_ptr);



}



int main(int argc, char** argv) {

  // construct, acquire, access, release
#if defined( BLISS_MUTEX )
  constexpr bliss::concurrent::LockType lt = bliss::concurrent::LockType::MUTEX;
#elif defined( BLISS_SPINLOCK )
  constexpr bliss::concurrent::LockType lt = bliss::concurrent::LockType::SPINLOCK;
#elif defined( BLISS_NONE )
  constexpr bliss::concurrent::LockType lt = bliss::concurrent::LockType::NONE;
#endif

  std::chrono::high_resolution_clock::time_point t1, t2;
  std::chrono::duration<double> time_span;


  //////////////  unbounded version

  /// thread unsafe.  test in single thread way.
  int count = 10000000;

#if defined(BLISS_NONE)

  for (int i = 1; i <= 8; ++i) {
    t1 = std::chrono::high_resolution_clock::now();
    timeAppendMultipleBuffers<bliss::io::ObjectPool<lt, bliss::io::Buffer<bliss::concurrent::LockType::LOCKFREE, 8192> > >(1, count, lt, bliss::concurrent::LockType::LOCKFREE, 8192);
    t2 = std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
//    printf("  time: entries %d nthread %d time %f\n", count, 1, time_span.count());
  }

#else
  for (int i = 1; i <= 8; ++i) {
    t1 = std::chrono::high_resolution_clock::now();
    timeAppendMultipleBuffers<bliss::io::ObjectPool<lt, bliss::io::Buffer<bliss::concurrent::LockType::LOCKFREE, 8192> > >(i, count, lt, bliss::concurrent::LockType::LOCKFREE, 8192);
    t2 = std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
//    printf("  time: entries %d nthread %d time %f\n", count, 1, time_span.count());
  }
#endif


}
