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


template<bliss::concurrent::LockType TS, int64_t CAP>
void append(const int nthreads, bliss::io::Buffer<TS, CAP>& buf, const unsigned int start, const unsigned int end, int& success, int& failure, int& swap, std::vector<int>& gold) {


  unsigned int i;
  int lsuccess = 0;
  int lfailure = 0;
  int lswap = 0;
#pragma omp parallel for num_threads(nthreads) default(none) private(i) shared(buf, gold, stdout) reduction(+: lsuccess, lfailure, lswap)
  for (i = start; i < end; ++i) {
    int data = static_cast<int>(i);
    int result = buf.append(&data, sizeof(int));

    if ((result & 0x1) > 0) {
      ++lsuccess;

#pragma omp critical
      gold.push_back(data);

    } else {
      ++lfailure;
    }

    if ((result & 0x2) > 0) {
      if (!buf.is_read_only()) {
        fprintf(stdout, "FAIL append: at this point the buffer should be in read state.\n");

      }

      ++lswap;
    }
  }

  success += lsuccess;
  failure += lfailure;
  swap += lswap;
  //printf("DEBUG: threads %d (actual,added/expected) success (%d/%d), failure (%d/%d), swap(%d)\n", nthreads, success, end - start, failure, end-start, swap);

}


template<bliss::concurrent::LockType TS, int64_t CAP, int NumThreads = 1>
void appendTest() {
  static_assert(NumThreads > 0, "instantiated with NumThreads < 1");
  static_assert(TS != bliss::concurrent::LockType::NONE || NumThreads == 1, "instantiated with Thread Unsafe version and NumThreads != 1");

  printf("TESTING operations on locktype %d buffer\n", static_cast<int>(TS) );


  // create a buffer.
  bliss::io::Buffer<TS, CAP> b1;
  b1.clear_and_unblock_writes();

  int nelems = CAP / sizeof(int);

  int success = 0;
  int failure = 0;
  int swap = 0;
  std::vector<int> gold;

  printf("TEST insert under capacity: ");
  append(NumThreads, b1, 0, nelems/2, success, failure, swap, gold);
  if (success == 0 || (success != nelems/2) || failure != 0 || swap != 0) printf("FAIL: (actual,added/expected) success (%d,%d/%d), failure (%d,%d/%d), swap(%d,%d/%d)\n", success, success, nelems/2, failure, failure, 0, swap, swap, 0);
  else {
    // compare unordered buffer content.
    if (compareUnorderedSequences(b1.operator int*(), gold.begin(), success)) {
      printf("PASS success %d failure %d swap %d\n", success, failure, swap);
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

  if (success == 0 || (success != nelems) || failure != nelems || swap != 1) printf("FAIL: (actual,added/expected) success (%d,%d/%d), failure (%d,%d/%d), swap(%d,%d/%d)\n", success, success2, nelems, failure, failure2, nelems, swap, swap2, 1);
  else {
    // compare unordered buffer content.
    if (compareUnorderedSequences(b1.operator int*(), gold.begin(), success)) {
      printf("PASS success %d failure %d swap %d\n", success, failure, swap);
    } else {
      printf("FAIL: content not matching\n");
    }
  }


  printf("TEST clear: ");
  b1.clear_and_block_writes();
  if (b1.getSize() != 0) printf("\tFAIL: NOT empty:  Size: %ld\n", b1.getSize());
  else printf("PASS\n");



  success = 0;
  failure = 0;
  swap = 0;
  gold.clear();
  b1.unblock_writes();

  printf("TEST insert AT capacity: ");

  append(NumThreads, b1, 0, nelems, success, failure, swap, gold);

  if (success == 0 || (success != nelems) || failure != 0 || swap != 0) printf("FAIL: (actual/expected) success (%d/%d), failure (%d/%d), swap(%d/%d)\n", success, nelems, failure, 0, swap, 0);
  else {
    // compare unordered buffer content.
    if (compareUnorderedSequences(b1.operator int*(), gold.begin(), success)) {
      printf("PASS success %d failure %d swap %d\n", success, failure, swap);
    } else {
      printf("FAIL: content not matching\n");
    }
  }

  b1.clear_and_unblock_writes();

  success = 0;
  failure = 0;
  swap = 0;
  gold.clear();

  printf("TEST insert JUST OVER capacity: ");

  append(NumThreads, b1, 0, nelems + NumThreads, success, failure, swap, gold);

  if (success == 0 || (success != nelems) || failure != NumThreads || swap != 1) printf("FAIL: (actual/expected) success (%d/%d), failure (%d/%d), swap(%d/%d)\n", success, nelems, failure, NumThreads, swap, 1);
  else {
    // compare unordered buffer content.
    if (compareUnorderedSequences(b1.operator int*(), gold.begin(), success)) {
      printf("PASS success %d failure %d swap %d\n", success, failure, swap);
    } else {
      printf("FAIL: content not matching\n");
    }
  }


  printf("TEST blocked buffer: ");
  b1.clear_and_block_writes();
  b1.block_and_flush();

  success = 0;
  failure = 0;
  swap = 0;
  gold.clear();

  append(NumThreads, b1, 0, nelems, success, failure, swap, gold);

  if ((success != 0) || failure != nelems || swap != 0) printf("FAIL: (actual/expected) success (%d/%d), failure (%d/%d), swap(%d/%d)\n", success, 0, failure, nelems, swap, 0);
  else {
    // compare unordered buffer content.
    if (compareUnorderedSequences(b1.operator int*(), gold.begin(), success)) {
      printf("PASS success %d failure %d swap %d\n", success, failure, swap);
    } else {
      printf("FAIL: content not matching\n");
    }
  }

  printf("TEST unblock buffer: ");

  b1.clear_and_unblock_writes();

  success = 0;
  failure = 0;
  swap = 0;
  gold.clear();

  append(NumThreads, b1, 0, nelems, success, failure, swap, gold);

  if (success == 0 || (success != nelems) || failure != 0 || swap != 0) printf("FAIL: (actual/expected) success (%d/%d), failure (%d/%d), swap(%d/%d)\n", success, nelems, failure, 0, swap, 0);
  else {
    // compare unordered buffer content.
    if (compareUnorderedSequences(b1.operator int*(), gold.begin(), success)) {
      printf("PASS success %d failure %d swap %d\n", success, failure, swap);
    } else {
      printf("FAIL: content not matching\n");
    }
  }

}



template<bliss::concurrent::LockType TS, int64_t CAP, int NumThreads>
void testAppendMultipleBuffersAtomicPtrs(const int total_count) {

  printf("TESTING atomic_ptrs: %d threads, locktype %d append with %ld bufferSize and %d total counts\n", NumThreads, static_cast<int>(TS), CAP, total_count);
  omp_lock_t writelock;
  omp_init_lock(&writelock);
  omp_lock_t writelock2;
  omp_init_lock(&writelock2);
  omp_lock_t writelock3;
  omp_init_lock(&writelock3);


  printf("TEST: save full buffers and process at end: ");
  std::vector<std::unique_ptr<bliss::io::Buffer<TS, CAP> > > full;

  std::vector<int> gold;
  std::vector<int> stored;

  int success = 0;
  int failure = 0;
  int swap = 0;
  int i = 0;

  std::atomic<bliss::io::Buffer<TS, CAP>* > ptr(new bliss::io::Buffer<TS, CAP>());                            // ensure atomicity
  ptr.load()->unblock_writes();

#pragma omp parallel for num_threads(NumThreads) default(none) shared(ptr, full, gold, stderr, stdout, std::cout, writelock, writelock2, writelock3) private(i) reduction(+:success, failure, swap)
  for (i = 0; i < total_count; ++i) {

    int data = static_cast<int>(i);
    //auto buf = ptr.load();

    std::atomic_thread_fence(std::memory_order_seq_cst);

    unsigned int result = ptr.load()->append(&data, sizeof(int));

    if (result & 0x1) {
      ++success;

      omp_set_lock(&writelock2);
      gold.push_back(data);
      omp_unset_lock(&writelock2);

    } else {
      ++failure;
      _mm_pause();  // slow it down a little.
    }

    if (result & 0x2 || result & 0x4) {

//      if (result & 0x4) std::cout << "SWAPPING: " << *(buf) << std::endl << std::flush;
//
//      if (!buf->is_read_only()) {
//        fprintf(stdout, "FAIL atomic batched proc: at this point the buffer should be in read state.\n");
//        fflush(stdout);
//        std::cout << "buffer: " << *(buf) << std::endl << std::flush;
//
//      }

      // swap in a new one.
      bliss::io::Buffer<TS, CAP>* new_ptr = new bliss::io::Buffer<TS, CAP>();  // manage new buffer
	    new_ptr->unblock_writes();

      bliss::io::Buffer<TS, CAP>* old_ptr = nullptr;

      old_ptr = ptr.exchange(new_ptr);
#pragma omp flush(ptr)
//
      // save the old buffer



      // this is showing a possible spurious wakeup...
      int oldsize = old_ptr->getSize() / sizeof(int);
      if (oldsize != CAP / sizeof(int)) {
        fprintf(stdout, "FAIL atomic DID NOT GET 2047 elements 1. local swap = %d, i = %d\n", swap, i);
        std::cout << "   atomic old buf: " << *(old_ptr) << std::endl
            << "   atomic new buf: " << *(ptr.load()) << std::endl << std::flush;
      }


      omp_set_lock(&writelock);
      full.push_back(std::move(std::unique_ptr<bliss::io::Buffer<TS, CAP> >(old_ptr)));
      omp_unset_lock(&writelock);

      ++swap;
    }

  }
  //printf("LAST BUFFER 1\n");

  ptr.load()->block_and_flush();
  int last = ptr.load()->getSize();
  if (last == (CAP / sizeof(int))) {
    ++swap;
  }

  auto b = ptr.exchange(nullptr);
  full.push_back(std::move(std::unique_ptr<bliss::io::Buffer<TS, CAP> >(b)));

//  printf("DEBUG: atomic 1 success %d, failure %d, swap %d, total %d, full count %ld\n", success, failure, swap, total_count, full.size());
//  std::cout << " buffer: " << *(ptr.load()) << std::endl << std::flush;

  for (int i = 0; i < full.size(); ++i) {

    stored.insert(stored.end(), full.at(i)->operator int*(), full.at(i)->operator int*() + full.at(i)->getSize() / sizeof(int));
  }
  int stored_count = stored.size();


  if (success == 0 || swap != full.size() - 1  || swap != success / (CAP / sizeof(int)) || success != stored_count)
    printf("FAIL atomic: (actual/expected)  success (%d/%d), failure (%d/?), last %d, swap(%d,%ld/%ld), last buf size %d, content match? %s.\n", stored_count, success, failure, last, swap, full.size(), success / (CAP / sizeof(int)), last, compareUnorderedSequences(stored.begin(), gold.begin(), stored_count) ? "same" : "diff");

  else {
    if (compareUnorderedSequences(stored.begin(), gold.begin(), stored_count)) {
      printf("PASS: atomic success %d, failure %d, swap %d, total %d\n", success, failure, swap, total_count);
    } else {
      printf("FAIL: atomic success %d, failure %d, swap %d, total %d, content not matching\n", success, failure, swap, total_count);
    }
  }



  printf("TEST: process full buffers along the way (SAVE IN VECTOR): ");

  omp_set_lock(&writelock);
  full.clear();   // deletes all the buffers in it.
  omp_unset_lock(&writelock);

  gold.clear();
  stored.clear();

  success = 0;
  failure = 0;
  swap = 0;
  i = 0;
  int success2 = 0;

  b = ptr.exchange(new bliss::io::Buffer<TS, CAP>());  // old pointer was managed by unique ptr.
  ptr.load()->unblock_writes();

#pragma omp parallel for num_threads(NumThreads) default(none) shared(ptr, gold, stored, stderr, stdout, std::cout,writelock, writelock2, writelock3, full) private(i) reduction(+:success, failure, swap, success2)
  for (i = 0; i < total_count; ++i) {

    std::atomic_thread_fence(std::memory_order_seq_cst);


    int data = static_cast<int>(i);
    //auto buf = ptr.load();
    std::atomic_thread_fence(std::memory_order_seq_cst);


    int res = ptr.load()->append(&data, sizeof(int));

    if (res & 0x1) {
      ++success;

      omp_set_lock(&writelock2);
      gold.push_back(data);
      omp_unset_lock(&writelock2);

    } else {
       ++failure;
//       _mm_pause();  // slow it down a little.
    }

    if (res & 0x2 || res & 0x4) {

//      if (res & 0x4) std::cout << "SWAPPING: " << *(buf) << std::endl << std::flush;
//    	// TODO: issue here:  if a large number of threads call append, and most of them are rescheduled, so that we reach calc
//    	// of pointer for a large number of threads in progress.  Then we could have the "just overflowing" thread executing and returning
//    	// 0x2 before all the memcpy are completed.  thus we could get is_read_only() failed while result is 0x2, and also observe a large
//    	// number of writes after result is set to 0x2 (and before that the flush bit is set)
//    	// this is a theory.
//
//      if (!(buf->is_read_only())) {
//        fprintf(stdout, "FAIL atomic incremental proc: at this point the buffer should be in read state.  res= %d\n", res);
//        fflush(stdout);
//        std::cout << "buffer: " << *(buf) << std::endl << std::flush;
//      }

      bliss::io::Buffer<TS, CAP>* new_ptr = new bliss::io::Buffer<TS, CAP>();  // manage new buffer
      //std::cout << "   new buf before assing: " << *(new_ptr) << std::endl <<  std::flush;

      bliss::io::Buffer<TS, CAP>* old_ptr = nullptr;


      new_ptr->unblock_writes();

      old_ptr = ptr.exchange(new_ptr);             //
#pragma omp flush(ptr)
      // save the old buffer

        if (old_ptr != nullptr) {
          ++swap;
          int oldsize = old_ptr->getSize() / sizeof(int);
  //        int newsize = buf_ptr.load()->getSize() / sizeof(int);
          if (oldsize != CAP / sizeof(int) || !(old_ptr->is_read_only())) {
            fprintf(stdout, "FAIL: atomic DID NOT GET 2047 elements 2. local swap = %d, i = %d\n", swap, i);
            std::cout << "   old buf: " << *(old_ptr) << std::endl <<  std::flush;
          }
          success2 += oldsize;

          omp_set_lock(&writelock);
            stored.insert(stored.end(), old_ptr->operator int*(), old_ptr->operator int*() + oldsize);
            full.push_back(std::move(std::unique_ptr<bliss::io::Buffer<TS, CAP> >(old_ptr)));
          omp_unset_lock(&writelock);

        }

    }

  }



  //printf("LAST BUFFER 2\n");
  ptr.load()->block_and_flush();
  last = ptr.load()->getSize();
  if (last == (CAP / sizeof(int))) {
    ++swap;
  }

  stored_count = stored.size();

  //printf("DEBUG: atomic before last buffer (actual/expected)  success (%d,%d/%d), failure (%d/?), swap(%d/%ld). content match? %s\n", stored_count, success2, success, failure, swap, success / (CAP / sizeof(int)), compareUnorderedSequences(stored.begin(), gold.begin(), stored_count) ? "same" : "diff");


  // compare unordered buffer content.
    stored.insert(stored.end(), ptr.load()->operator int*(), ptr.load()->operator int*() + ptr.load()->getSize() / sizeof(int));
    full.push_back(std::move(std::unique_ptr<bliss::io::Buffer<TS, CAP> >(ptr.load())));

  stored_count = stored.size();
  success2 += ptr.load()->getSize() / sizeof(int);


  //printf("DEBUG: atomic after last buffer (actual/expected)  success (%d,%d/%d), failure (%d/?), swap(%d/%ld), final buf size %d, content match? %s\n", stored_count, success2, success, failure, swap, success / (CAP / sizeof(int)), last, compareUnorderedSequences(stored.begin(), gold.begin(), stored_count) ? "same" : "diff");

  if ( success == 0 || swap != success / (CAP / sizeof(int)) || success != stored_count)
    printf("FAIL atomic: (actual/expected)  success (%d,%d/%d), failure (%d/?), last %d, swap(%d/%ld). content match? %s\n", stored_count, success2, success, failure, last, swap, success / (CAP / sizeof(int)), compareUnorderedSequences(stored.begin(), gold.begin(), stored_count) ? "same" : "diff");
  else {

    if (compareUnorderedSequences(stored.begin(), gold.begin(), stored_count)) {
      printf("PASS: atomic success %d, failure %d, swap %d, total %d\n", success, failure, swap, total_count);
    } else {
      printf("FAIL: atomic success %d, failure %d, swap %d, total %d, content not matching\n", success, failure, swap, total_count);
    }
  }


  ptr.exchange(nullptr);
  full.clear();

  omp_destroy_lock(&writelock);
  omp_destroy_lock(&writelock2);
  omp_destroy_lock(&writelock3);
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




  // test append
  appendTest<bliss::concurrent::LockType::NONE, 8192, 1>();

  appendTest<lt, 8192, 1>();
  appendTest<lt, 8192, 2>();
  appendTest<lt, 8192, 3>();
  appendTest<lt, 8192, 4>();
  appendTest<lt, 8192, 5>();
  appendTest<lt, 8192, 6>();
  appendTest<lt, 8192, 7>();
  appendTest<lt, 8192, 8>();

  // test append with buffer that is not multple of element size.
  appendTest<bliss::concurrent::LockType::NONE, 8191, 1>();

  appendTest<lt, 8191, 1>();
  appendTest<lt, 8191, 2>();
  appendTest<lt, 8191, 3>();
  appendTest<lt, 8191, 4>();
  appendTest<lt, 8191, 5>();
  appendTest<lt, 8191, 6>();
  appendTest<lt, 8191, 7>();
  appendTest<lt, 8191, 8>();



  // multiple buffer swap test.

  ////////////// timing.  the insert before this is to warm up.
  testAppendMultipleBuffersAtomicPtrs<bliss::concurrent::LockType::NONE, 8191, 1>(1000000);

  testAppendMultipleBuffersAtomicPtrs<lt, 8191, 1>(1000000);
  testAppendMultipleBuffersAtomicPtrs<lt, 8191, 2>(1000000);
  testAppendMultipleBuffersAtomicPtrs<lt, 8191, 3>(1000000);
  testAppendMultipleBuffersAtomicPtrs<lt, 8191, 4>(1000000);
  testAppendMultipleBuffersAtomicPtrs<lt, 8191, 5>(1000000);
  testAppendMultipleBuffersAtomicPtrs<lt, 8191, 6>(1000000);
  testAppendMultipleBuffersAtomicPtrs<lt, 8191, 7>(1000000);
  testAppendMultipleBuffersAtomicPtrs<lt, 8191, 8>(1000000);


}
