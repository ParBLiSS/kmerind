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


#include "io/buffer.hpp"
#include "omp.h"
#include <cassert>
#include <chrono>
#include <iterator>  // for ostream_iterator
#include <iostream>   // for cout
#include <sstream>

#include <utils/test_utils.hpp>



template<bliss::concurrent::ThreadSafety TS1, bliss::concurrent::ThreadSafety TS2>
void moveTest( const int count = 1000, const int capacity = 8192) {
  printf("TESTING move from thread %s to thread %s\n", TS1 ? "SAFE" : "UNSAFE", TS2 ? "SAFE" : "UNSAFE");


  bliss::io::Buffer<TS1> b1(capacity);
  assignSequence(b1.operator int*(), count);
  b1.block();

  // verify that what went in is correct
  if (!checkSequence(b1.operator int*(), count)) {
    printf("FAIL: append is not working correctly\n");
  }


  size_t bufferCapBefore = b1.getCapacity();
  size_t bufferCapAfter = 0;

  int bufferSizeBefore = b1.getFinalSize();
  int bufferSizeAfter = 0;

  int* bufferPtrBefore = b1.operator int*();
  int* bufferPtrAfter = nullptr;


  printf("TEST move ctor: ");
  bliss::io::Buffer<TS2> b2(std::move(b1));
  if (b1.getCapacity()  !=  bufferCapAfter  ||
      b1.getApproximateSize()      != bufferSizeAfter ||
      b1.operator int*()      !=  bufferPtrAfter  ||
      b2.getCapacity() !=  bufferCapBefore  ||
      b2.getApproximateSize()     != bufferSizeBefore ||
      b2.operator int*()     !=  bufferPtrBefore)
    printf(" FAIL: size %u->%u, capacity %u->%u, pointer = %p->%p\n",
         b1.getApproximateSize(), b2.getApproximateSize(),
         b1.getCapacity(),        b2.getCapacity(),
         b1.operator int*(),            b2.operator int*());
  else {

    if (checkSequence(b2.operator int*(), count)) {
      printf("PASS\n");
    } else {
      printf("FAIL: content mismatch\n");
    }
  }

  bliss::io::Buffer<TS1> b3(8192);

  assignSequence(b3.operator int*(), count);
  b3.block();

  bufferCapBefore = b3.getCapacity();
  bufferCapAfter = 0;

  bufferSizeBefore = b3.getFinalSize();
  bufferSizeAfter = 0;

  bufferPtrBefore = b3.operator int*();
  bufferPtrAfter = nullptr;


  printf("TEST move = operator\n");
  b2 = std::move(b3);
  if (b3.getCapacity()  !=  bufferCapAfter  ||
      b3.getApproximateSize()      != bufferSizeAfter ||
      b3.operator int*()      !=  bufferPtrAfter  ||
      b2.getCapacity() !=  bufferCapBefore  ||
      b2.getApproximateSize()     != bufferSizeBefore ||
      b2.operator int*()     !=  bufferPtrBefore)
    printf(" FAIL: size %u->%u, capacity %u->%u, pointer = %p->%p\n",
         b3.getApproximateSize(), b2.getApproximateSize(),
         b3.getCapacity(),        b2.getCapacity(),
         b3.operator int*(),            b2.operator int*());
  else {

    if (checkSequence(b2.operator int*(), count)) {
      printf("PASS\n");
    } else {
      printf("FAIL: content mismatch\n");
    }
  }
}

template<bliss::concurrent::ThreadSafety TS>
void append(const int nthreads, bliss::io::Buffer<TS>& buf, const unsigned int start, const unsigned int end, int& success, int& failure, int& swap, std::vector<int>& gold) {

  unsigned int i;
  int lsuccess = 0;
  int lfailure = 0;
  int lswap = 0;
#pragma omp parallel for num_threads(nthreads) default(none) private(i) shared(buf, gold) reduction(+: lsuccess, lfailure, lswap)
  for (i = start; i < end; ++i) {
    int data = static_cast<int>(i);
    std::pair<bool, bool> result = buf.append(&data, sizeof(int));

    if (result.first) {
      ++lsuccess;

#pragma omp critical
      gold.push_back(data);

    } else {
      ++lfailure;
    }

    if (result.second) {
      ++lswap;
    }
  }

  success += lsuccess;
  failure += lfailure;
  swap += lswap;

}


template<bliss::concurrent::ThreadSafety TS, int NumThreads = 1>
void appendTest(const int capacity = 8192) {
  static_assert(NumThreads > 0, "instantiated with NumThreads < 1");
  static_assert(TS || NumThreads == 1, "instantiated with Thread Unsafe version and NumThreads != 1");

  printf("TESTING operations on thread %s buffer\n", TS ? "SAFE" : "UNSAFE");


  // create a buffer.
  bliss::io::Buffer<TS> b1(capacity);
  int nelems = capacity / sizeof(int);

  int success = 0;
  int failure = 0;
  int swap = 0;
  std::vector<int> gold;

  printf("TEST insert under capacity: ");
  append(NumThreads, b1, 0, nelems/2, success, failure, swap, gold);
  if ((success != nelems/2) || failure != 0 || swap != 0) printf("FAIL: (actual/expected) success (%d/%d), failure (%d/%d), swap(%d/%d)\n", success, nelems/2, failure, 0, swap, 0);
  else {
    // compare unordered buffer content.
    if (compareUnorderedSequences(b1.operator int*(), gold.begin(), success)) {
      printf("PASS\n");
    } else {
      printf("FAIL: content not matching\n");
    }
  }


  int success2 = 0;
  int failure2 = 0;
  int swap2 = 0;
  printf("TEST insert over capacity: ");

  append(NumThreads, b1, nelems/2, nelems * 2, success2, failure2, swap2, gold);

  success += success2;
  failure += failure2;
  swap += swap2;

  if ((success != nelems) || failure != nelems || swap != 1) printf("FAIL: (actual/expected) success (%d/%d), failure (%d/%d), swap(%d/%d)\n", success, nelems, failure, nelems, swap, 1);
  else {
    // compare unordered buffer content.
    if (compareUnorderedSequences(b1.operator int*(), gold.begin(), success)) {
      printf("PASS\n");
    } else {
      printf("FAIL: content not matching\n");
    }
  }


  printf("TEST clear: ");
  b1.clear();
  if (b1.getFinalSize() != -1 || b1.getApproximateSize() != 0) printf("\tFAIL: NOT empty:  finalSize: %d, approx Size: %d\n", b1.getFinalSize(), b1.getApproximateSize());
  else printf("PASS\n");



  success = 0;
  failure = 0;
  swap = 0;
  gold.clear();
  b1.unblock();

  printf("TEST insert AT capacity: ");

  append(NumThreads, b1, 0, nelems, success, failure, swap, gold);

  if ((success != nelems) || failure != 0 || swap != 0) printf("FAIL: (actual/expected) success (%d/%d), failure (%d/%d), swap(%d/%d)\n", success, nelems, failure, 0, swap, 0);
  else {
    // compare unordered buffer content.
    if (compareUnorderedSequences(b1.operator int*(), gold.begin(), success)) {
      printf("PASS\n");
    } else {
      printf("FAIL: content not matching\n");
    }
  }

  b1.clear();
  b1.unblock();

  success = 0;
  failure = 0;
  swap = 0;
  gold.clear();

  printf("TEST insert JUST OVER capacity: ");

  append(NumThreads, b1, 0, nelems + NumThreads, success, failure, swap, gold);

  if ((success != nelems) || failure != NumThreads || swap != 1) printf("FAIL: (actual/expected) success (%d/%d), failure (%d/%d), swap(%d/%d)\n", success, nelems, failure, NumThreads, swap, 1);
  else {
    // compare unordered buffer content.
    if (compareUnorderedSequences(b1.operator int*(), gold.begin(), success)) {
      printf("PASS\n");
    } else {
      printf("FAIL: content not matching\n");
    }
  }


  printf("TEST blocked buffer: ");
  b1.clear();
  b1.block();

  success = 0;
  failure = 0;
  swap = 0;
  gold.clear();

  append(NumThreads, b1, 0, nelems, success, failure, swap, gold);

  if ((success != 0) || failure != nelems || swap != 0) printf("FAIL: (actual/expected) success (%d/%d), failure (%d/%d), swap(%d/%d)\n", success, 0, failure, nelems, swap, 0);
  else {
    // compare unordered buffer content.
    if (compareUnorderedSequences(b1.operator int*(), gold.begin(), success)) {
      printf("PASS\n");
    } else {
      printf("FAIL: content not matching\n");
    }
  }

  printf("TEST unblock buffer: ");

  b1.clear();
  b1.unblock();

  success = 0;
  failure = 0;
  swap = 0;
  gold.clear();

  append(NumThreads, b1, 0, nelems, success, failure, swap, gold);

  if ((success != nelems) || failure != 0 || swap != 0) printf("FAIL: (actual/expected) success (%d/%d), failure (%d/%d), swap(%d/%d)\n", success, nelems, failure, 0, swap, 0);
  else {
    // compare unordered buffer content.
    if (compareUnorderedSequences(b1.operator int*(), gold.begin(), success)) {
      printf("PASS\n");
    } else {
      printf("FAIL: content not matching\n");
    }
  }

}


template<bliss::concurrent::ThreadSafety TS, int NumThreads>
void testAppendMultipleBuffers(const int buffer_capacity, const int total_count) {

  printf("TESTING: %d threads, thread %s append with %d bufferSize and %d total counts\n", NumThreads, TS ? "SAFE" : "UNSAFE", buffer_capacity, total_count);


  printf("TEST: save full buffers and process at end: ");
  std::vector<std::unique_ptr<bliss::io::Buffer<TS> > > full;

  std::vector<int> gold;
  std::vector<int> stored;

  int data = 0;
  std::pair<bool, bool> result { false, false };

  int success = 0;
  int failure = 0;
  int swap = 0;
  int i = 0;

  std::unique_ptr<bliss::io::Buffer<TS> > buf_ptr(new bliss::io::Buffer<TS>(buffer_capacity));

#pragma omp parallel for num_threads(NumThreads) default(none) shared(buf_ptr, full, gold) private(i, data, result) reduction(+:success, failure, swap)
  for (i = 0; i < total_count; ++i) {

    bliss::io::Buffer<TS>* ptr = buf_ptr.get();

    data = static_cast<int>(i);
    result = ptr->append(&data, sizeof(int));

    if (result.first) {
      ++success;
#pragma omp critical
      gold.push_back(data);
    }
    else ++failure;

    if (result.second) {
      ++swap;

      // swap in a new one.
      std::unique_ptr<bliss::io::Buffer<TS> > new_buf_ptr(new bliss::io::Buffer<TS>(buffer_capacity));
      buf_ptr.swap(new_buf_ptr);
      // save the old buffer
      full.push_back(std::move(new_buf_ptr));

    }

  }

  buf_ptr->block();

  for (int i = 0; i < full.size(); ++i) {
		full.at(i)->waitForAllUpdates();
    stored.insert(stored.end(), full.at(i)->operator int*(), full.at(i)->operator int*() + full.at(i)->getFinalSize() / sizeof(int));
  }

buf_ptr->waitForAllUpdates();
  // compare unordered buffer content.
  stored.insert(stored.end(), buf_ptr->operator int*(), buf_ptr->operator int*() + buf_ptr->getFinalSize() / sizeof(int));
  int stored_count = stored.size();


  if (swap != full.size()  || swap != success / (buffer_capacity / sizeof(int)) || success != stored_count)
      printf("FAIL: (actual/expected)  success (%d/%d), failure (%d/?), swap(%ld,%lu/%d).\n", success, stored_count, failure, success / (buffer_capacity / sizeof(int)), full.size(), swap);
  else {
    printf("INFO: success %d, failure %d, swap %d, total %d\n", success, failure, swap, total_count);

    if (compareUnorderedSequences(stored.begin(), gold.begin(), stored_count)) {
      printf("PASS\n");
    } else {
      printf("FAIL: content not matching\n");
    }
  }



  printf("TEST: process full buffers along the way: ");

  full.clear();

  gold.clear();
  stored.clear();

  data = 0;

  success = 0;
  failure = 0;
  swap = 0;
  i = 0;

  buf_ptr.reset(new bliss::io::Buffer<TS>(buffer_capacity));

#pragma omp parallel for num_threads(NumThreads) default(none) shared(buf_ptr, gold, stored) private(i, data, result) reduction(+:success, failure, swap)
  for (i = 0; i < total_count; ++i) {

    std::atomic_thread_fence(std::memory_order_seq_cst);

    bliss::io::Buffer<TS>* ptr = buf_ptr.get();



    data = static_cast<int>(i);
    result = ptr->append(&data, sizeof(int));

    if (result.first) {
      ++success;
#pragma omp critical
      gold.push_back(data);
    }
    else ++failure;

    if (result.second) {
      ++swap;

			buf_ptr->waitForAllUpdates();
      // swap in a new one.
      std::unique_ptr<bliss::io::Buffer<TS> > new_buf_ptr(new bliss::io::Buffer<TS>(buffer_capacity));
      buf_ptr.swap(new_buf_ptr);
      // save the old buffer

#pragma omp critical
      {
        stored.insert(stored.end(), new_buf_ptr->operator int*(), new_buf_ptr->operator int*() + new_buf_ptr->getFinalSize() / sizeof(int));
      }

    }

  }

  buf_ptr->block();


  // compare unordered buffer content.
  stored.insert(stored.end(), buf_ptr->operator int*(), buf_ptr->operator int*() + buf_ptr->getFinalSize() / sizeof(int));
  stored_count = stored.size();


  if ( swap != success / (buffer_capacity / sizeof(int)) || success != stored_count)
    printf("FAIL: (actual/expected)  success (%d/%d), failure (%d/?), swap(%ld/%d).\n", success, stored_count, failure, success / (buffer_capacity / sizeof(int)), swap);
  else {
    printf("INFO: success %d, failure %d, swap %d, total %d\n", success, failure, swap, total_count);

    if (compareUnorderedSequences(stored.begin(), gold.begin(), stored_count)) {
      printf("PASS\n");
    } else {
      printf("FAIL: content not matching\n");
    }
  }


}




template<bliss::concurrent::ThreadSafety TS, int NumThreads>
void appendTimed(const int capacity, const int iterations) {

  printf("PROFILING: %d threads, thread %s buffer append with %d elements and %d iterations\n", NumThreads, TS ? "SAFE" : "UNSAFE", capacity, iterations);
  bliss::io::Buffer<TS> buf(capacity);

  std::chrono::high_resolution_clock::time_point t1, t2;
  std::chrono::duration<double> time_span;

  int success;
  int failure;
  int i;
  int count = capacity / sizeof(int);
  int data;
  std::pair<bool, bool> result;

  for (int k = 0; k < iterations; ++k) {

      buf.clear();
      buf.unblock();
      t1 = std::chrono::high_resolution_clock::now();

      i = 0;
      success = 0;
      failure = 0;


#pragma omp parallel for num_threads(NumThreads) default(none) shared(buf, count) private(i, data, result) reduction(+:success, failure)
      for (i = 0; i < count + NumThreads; ++i) {
        data = static_cast<int>(i);
        result = buf.append(&data, sizeof(int));

        if (result.first) ++success;
        else ++failure;
      }

      assert(success == count && failure == NumThreads);

      t2 = std::chrono::high_resolution_clock::now();

      time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
      printf("Append %d (success) and %d (failure) elements duration = %f\n", success, failure, time_span.count());

  }

}


int main(int argc, char** argv) {

  // test move
  moveTest<bliss::concurrent::THREAD_UNSAFE, bliss::concurrent::THREAD_UNSAFE>();
  moveTest<bliss::concurrent::THREAD_UNSAFE, bliss::concurrent::THREAD_SAFE>();
  moveTest<bliss::concurrent::THREAD_SAFE, bliss::concurrent::THREAD_SAFE>();
  moveTest<bliss::concurrent::THREAD_SAFE, bliss::concurrent::THREAD_UNSAFE>();


  // test append
  appendTest<bliss::concurrent::THREAD_UNSAFE, 1>(8192);

  appendTest<bliss::concurrent::THREAD_SAFE, 1>(8192);
  appendTest<bliss::concurrent::THREAD_SAFE, 2>(8192);
  appendTest<bliss::concurrent::THREAD_SAFE, 3>(8192);
  appendTest<bliss::concurrent::THREAD_SAFE, 4>(8192);
  appendTest<bliss::concurrent::THREAD_SAFE, 5>(8192);
  appendTest<bliss::concurrent::THREAD_SAFE, 6>(8192);
  appendTest<bliss::concurrent::THREAD_SAFE, 7>(8192);
  appendTest<bliss::concurrent::THREAD_SAFE, 8>(8192);

  // test append with buffer that is not multple of element size.
  appendTest<bliss::concurrent::THREAD_UNSAFE, 1>(8191);

  appendTest<bliss::concurrent::THREAD_SAFE, 1>(8191);
  appendTest<bliss::concurrent::THREAD_SAFE, 2>(8191);
  appendTest<bliss::concurrent::THREAD_SAFE, 3>(8191);
  appendTest<bliss::concurrent::THREAD_SAFE, 4>(8191);
  appendTest<bliss::concurrent::THREAD_SAFE, 5>(8191);
  appendTest<bliss::concurrent::THREAD_SAFE, 6>(8191);
  appendTest<bliss::concurrent::THREAD_SAFE, 7>(8191);
  appendTest<bliss::concurrent::THREAD_SAFE, 8>(8191);


  // timing tests

  ////////////// timing.  the insert before this is to warm up.
  appendTimed<bliss::concurrent::THREAD_UNSAFE, 1>(1000000, 3);

  appendTimed<bliss::concurrent::THREAD_SAFE, 1>(1000000, 3);
  appendTimed<bliss::concurrent::THREAD_SAFE, 2>(1000000, 3);
  appendTimed<bliss::concurrent::THREAD_SAFE, 3>(1000000, 3);
  appendTimed<bliss::concurrent::THREAD_SAFE, 4>(1000000, 3);
  appendTimed<bliss::concurrent::THREAD_SAFE, 5>(1000000, 3);
  appendTimed<bliss::concurrent::THREAD_SAFE, 6>(1000000, 3);
  appendTimed<bliss::concurrent::THREAD_SAFE, 7>(1000000, 3);
  appendTimed<bliss::concurrent::THREAD_SAFE, 8>(1000000, 3);


  // multiple buffer swap test.

  ////////////// timing.  the insert before this is to warm up.
  testAppendMultipleBuffers<bliss::concurrent::THREAD_UNSAFE, 1>(8191, 1000000);

  testAppendMultipleBuffers<bliss::concurrent::THREAD_SAFE, 1>(8191, 1000000);
  testAppendMultipleBuffers<bliss::concurrent::THREAD_SAFE, 2>(8191, 1000000);
  testAppendMultipleBuffers<bliss::concurrent::THREAD_SAFE, 3>(8191, 1000000);
  testAppendMultipleBuffers<bliss::concurrent::THREAD_SAFE, 4>(8191, 1000000);
  testAppendMultipleBuffers<bliss::concurrent::THREAD_SAFE, 5>(8191, 1000000);
  testAppendMultipleBuffers<bliss::concurrent::THREAD_SAFE, 6>(8191, 1000000);
  testAppendMultipleBuffers<bliss::concurrent::THREAD_SAFE, 7>(8191, 1000000);
  testAppendMultipleBuffers<bliss::concurrent::THREAD_SAFE, 8>(8191, 1000000);


}
