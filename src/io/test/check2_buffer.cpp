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


#ifdef LOCKING
#include "io/locking_buffer.hpp"
#endif
#ifdef SPINLOCKING
#include "io/spinlocking_buffer.hpp"
#endif
#ifdef LOCKFREE
#include "io/lockfree_buffer.hpp"
#endif
#ifdef LOCKFREE_PTR
#include "io/buffer.hpp"
#endif



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



template<bliss::concurrent::ThreadSafety TS1, bliss::concurrent::ThreadSafety TS2>
void moveTest( const int count = 1000, const int capacity = 8192) {
  printf("TESTING move from thread %s to thread %s\n", TS1 ? "SAFE" : "UNSAFE", TS2 ? "SAFE" : "UNSAFE");


  bliss::io::Buffer<TS1> b1(capacity);
  b1.unblock();

//  std::cout << " empty buffer: " << b1 << std::endl << std::flush;
  for (size_t i = 0; i < count; ++i) {
	  int data = i;
    b1.append(&data, sizeof(int));
  }
//  std::cout << " initialized buffer: " << b1 << std::endl << std::flush;
//  b1.block_and_flush();

  // verify that what went in is correct
  if (!checkSequence(b1.operator int*(), count)) {
    printf("FAIL: append is not working correctly\n");
  }

  b1.block();

  int64_t bufferCapBefore = b1.getCapacity();
  int64_t bufferCapAfter = 0;

  int64_t bufferSizeBefore = b1.getSize();
  int64_t bufferSizeAfter = 0;

  int* bufferPtrBefore = b1.operator int*();
  int* bufferPtrAfter = nullptr;


  printf("TEST move ctor: ");
  bliss::io::Buffer<TS2> b2(std::move(b1));


  if (b1.getCapacity()  !=  bufferCapAfter  ||
      b1.getSize()      != bufferSizeAfter ||
      b1.operator int*()      !=  bufferPtrAfter  ||
      b2.getCapacity() !=  bufferCapBefore  ||
      b2.getSize()     != bufferSizeBefore ||
      b2.operator int*()     !=  bufferPtrBefore)
    printf(" FAIL: size %ld (%ld->%ld), capacity %ld (%ld->%ld), pointer = %p (%p->%p)\n",
         b2.getSize(), bufferSizeBefore, b1.getSize(),
         b2.getCapacity(),        bufferCapBefore,  b1.getCapacity(),
         b2.operator int*(),      bufferPtrBefore,   b1.operator int*()
         );
  else {

    if (checkSequence(b2.operator int*(), count)) {
      printf("PASS\n");
    } else {
      printf("FAIL: content mismatch\n");
    }
  }

  bliss::io::Buffer<TS1> b3(8192);
  b3.unblock();

  for (size_t i = 0; i < count; ++i) {
	  int data = i;
    b3.append(&data, sizeof(int));
  }
//  assignSequence(b3.operator int*(), count);
  b3.block_and_flush();

  bufferCapBefore = b3.getCapacity();
  bufferSizeBefore = b3.getSize();
  bufferPtrBefore = b3.operator int*();

  bufferCapAfter  = 0;
  bufferSizeAfter = 0;
  bufferPtrAfter  = nullptr;


  printf("TEST move = operator: ");
  b2 = std::move(b3);
  if (b3.getCapacity()  !=  bufferCapAfter  ||
      b3.getSize()      != bufferSizeAfter ||
      b3.operator int*()      !=  bufferPtrAfter  ||
      b2.getCapacity() !=  bufferCapBefore  ||
      b2.getSize()     != bufferSizeBefore ||
      b2.operator int*()     !=  bufferPtrBefore)
    printf(" FAIL: size %ld (%ld->%ld), capacity %ld (%ld->%ld), pointer = %p (%p->%p)\n",
         b2.getSize(), bufferSizeBefore, b3.getSize(),
         b2.getCapacity(),        bufferCapBefore,  b3.getCapacity(),
         b2.operator int*(),      bufferPtrBefore,   b3.operator int*()
         );
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
      if (!buf.is_reading()) {
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


template<bliss::concurrent::ThreadSafety TS, int NumThreads = 1>
void appendTest(const int capacity = 8192) {
  static_assert(NumThreads > 0, "instantiated with NumThreads < 1");
  static_assert(TS || NumThreads == 1, "instantiated with Thread Unsafe version and NumThreads != 1");

  printf("TESTING operations on thread %s buffer\n", TS ? "SAFE" : "UNSAFE");


  // create a buffer.
  bliss::io::Buffer<TS> b1(capacity);
  b1.clear();
  b1.unblock();

  int nelems = capacity / sizeof(int);

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
  b1.clear();
  b1.block();
  if (b1.getSize() != 0) printf("\tFAIL: NOT empty:  Size: %ld\n", b1.getSize());
  else printf("PASS\n");



  success = 0;
  failure = 0;
  swap = 0;
  gold.clear();
  b1.unblock();

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

  b1.clear();
  b1.unblock();

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
  b1.clear();
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

  b1.clear();
  b1.unblock();

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

//
//template<bliss::concurrent::ThreadSafety TS, int NumThreads>
//void testAppendMultipleBuffers(const int buffer_capacity, const int total_count) {
//
//  printf("TESTING: %d threads, thread %s append with %d bufferSize and %d total counts\n", NumThreads, TS ? "SAFE" : "UNSAFE", buffer_capacity, total_count);
//
//
//  printf("TEST: save full buffers and process at end: ");
//  std::vector<std::unique_ptr<bliss::io::Buffer<TS> > > full;
//
//  std::vector<int> gold;
//  std::vector<int> stored;
//
//  int data = 0;
//  unsigned int result = 0;
//
//  int success = 0;
//  int failure = 0;
//  int swap = 0;
//  int i = 0;
//
//  std::unique_ptr<bliss::io::Buffer<TS> > buf_ptr(new bliss::io::Buffer<TS>(buffer_capacity));
//  buf_ptr->unblock();
//
//#pragma omp parallel for num_threads(NumThreads) default(none) shared(buf_ptr, full, gold, stderr,stdout, std::cout) private(i, data, result) reduction(+:success, failure, swap)
//  for (i = 0; i < total_count; ++i) {
//
//    bliss::io::Buffer<TS>* ptr = buf_ptr.get();
//
//    data = static_cast<int>(i);
//    result = ptr->append(&data, sizeof(int));
//
//    if ((result & 0x1) > 0) {
//      ++success;
//#pragma omp critical
//      gold.push_back(data);
//    }
//    else {
//      ++failure;
//      _mm_pause();  // slow it down a little.
//    }
//
//    if ((result & 0x2) > 0) {
//      if (!ptr->is_reading()) {
//        fprintf(stdout, "FAIL batched proc: at this point the buffer should be in read state.\n");
//      }
//
//      // swap in a new one.
//      std::unique_ptr<bliss::io::Buffer<TS> > new_buf_ptr(new bliss::io::Buffer<TS>(buffer_capacity));
//      new_buf_ptr->block();
//
//      buf_ptr.swap(new_buf_ptr);
//      buf_ptr->unblock();
//      // save the old buffer
//
//      // this is showing a possible spurious wakeup...
//      int oldsize = new_buf_ptr->getSize() / sizeof(int);
//      if (oldsize != buffer_capacity / sizeof(int)) {
//        fprintf(stdout, "DID NOT GET 2047 elements 1. local swap = %d, i = %d\n", swap, i);
//        std::cout << "old buf: " << *(new_buf_ptr.get()) << std::endl
//            << "new buf: " << *(buf_ptr) << std::endl << std::flush;
//      }
//#pragma omp critical
//      {
//        full.push_back(std::move(new_buf_ptr));
//      }
//
//      ++swap;
//      buf_ptr->unblock();
//    }
//
//  }
//
//  buf_ptr->block_and_flush();
//
//  for (int i = 0; i < full.size(); ++i) {
//    stored.insert(stored.end(), full.at(i)->operator int*(), full.at(i)->operator int*() + full.at(i)->getSize() / sizeof(int));
//  }
//
//  buf_ptr->block_and_flush();
//  // compare unordered buffer content.
//  stored.insert(stored.end(), buf_ptr->operator int*(), buf_ptr->operator int*() + buf_ptr->getSize() / sizeof(int));
//  int stored_count = stored.size();
//
//
//  if (success == 0 || swap != full.size()  || swap != success / (buffer_capacity / sizeof(int)) || success != stored_count)
//    printf("FAIL: (actual/expected)  success (%d/%d), failure (%d/?), swap(%d,%ld/%ld), content match? %s.\n", stored_count, success, failure, swap, full.size(), success / (buffer_capacity / sizeof(int)), compareUnorderedSequences(stored.begin(), gold.begin(), stored_count) ? "same" : "diff");
//
//  else {
//    printf("INFO: success %d, failure %d, swap %d, total %d\n", success, failure, swap, total_count);
//
//    if (compareUnorderedSequences(stored.begin(), gold.begin(), stored_count)) {
//      printf("PASS\n");
//    } else {
//      printf("FAIL: content not matching\n");
//    }
//  }
//
//
//
//  printf("TEST: process full buffers along the way: ");
//
//  full.clear();
//
//  gold.clear();
//  stored.clear();
//
//  data = 0;
//
//  success = 0;
//  failure = 0;
//  swap = 0;
//  i = 0;
//  int success2 = 0;
//
//  buf_ptr.reset(new bliss::io::Buffer<TS>(buffer_capacity));
//  buf_ptr->unblock();
//
//#pragma omp parallel for num_threads(NumThreads) default(none) shared(buf_ptr, gold, stored, stderr, stdout, std::cout) private(i, data, result) reduction(+:success, failure, swap, success2)
//  for (i = 0; i < total_count; ++i) {
//
//    std::atomic_thread_fence(std::memory_order_seq_cst);
//
//    bliss::io::Buffer<TS>* ptr = buf_ptr.get();
//
//
//
//    data = static_cast<int>(i);
//    result = ptr->append(&data, sizeof(int));
//
//    if ((result & 0x1) > 0) {
//      ++success;
//#pragma omp critical
//      gold.push_back(data);
//    }
//    else {
//       ++failure;
//       _mm_pause();  // slow it down a little.
//    }
//
//    if ((result & 0x2) > 0) {
//
//      if (!ptr->is_reading()) {
//        fprintf(stdout, "FAIL incremental proc: at this point the buffer should be in read state.\n");
//      }
//
//
//      std::unique_ptr<bliss::io::Buffer<TS> > new_buf_ptr(new bliss::io::Buffer<TS>(buffer_capacity));
//
//        // try locking this part...  -
//      new_buf_ptr->block();
//
//      // swap in a new one.
//      buf_ptr.swap(new_buf_ptr);
//      buf_ptr->unblock();
//      // save the old buffer
//
//        // this part shows that some threads end up
//
//      // this is showing a possible spurious wake up.
//        int oldsize = new_buf_ptr->getSize() / sizeof(int);
////        int newsize = buf_ptr->getSize() / sizeof(int);
//        if (oldsize != buffer_capacity / sizeof(int)) {
//          fprintf(stdout, "DID NOT GET 2047 elements 2. local swap = %d, i = %d\n", swap, i);
//          std::cout << "old buf: " << *(new_buf_ptr.get()) << std::endl
//              << "new buf: " << *(buf_ptr) << std::endl << std::flush;
//        }
//
//#pragma omp atomic
//        success2 += oldsize;
//
//#pragma omp critical
//      {
//        stored.insert(stored.end(), new_buf_ptr->operator int*(), new_buf_ptr->operator int*() + oldsize);
//      }
//
//      ++swap;
//      buf_ptr->unblock();
//
//    }
//
//  }
//
//  buf_ptr->block_and_flush();
//
//
//  // compare unordered buffer content.
//  stored.insert(stored.end(), buf_ptr->operator int*(), buf_ptr->operator int*() + buf_ptr->getSize() / sizeof(int));
//  stored_count = stored.size();
//  success2 += buf_ptr->getSize() / sizeof(int);
//
//  if ( success == 0 || swap != success / (buffer_capacity / sizeof(int)) || success != stored_count)
//    printf("FAIL: (actual/expected)  success (%d,%d/%d), failure (%d/?), swap(%d/%ld). content match? %s\n", stored_count, success2, success, failure, swap, success / (buffer_capacity / sizeof(int)), compareUnorderedSequences(stored.begin(), gold.begin(), stored_count) ? "same" : "diff");
//  else {
//    printf("INFO: success %d, failure %d, swap %d, total %d\n", success, failure, swap, total_count);
//
//    if (compareUnorderedSequences(stored.begin(), gold.begin(), stored_count)) {
//      printf("PASS\n");
//    } else {
//      printf("FAIL: content not matching\n");
//    }
//  }
//
//
//}


template<bliss::concurrent::ThreadSafety TS, int NumThreads>
void testAppendMultipleBuffersAtomicPtrs(const int buffer_capacity, const int total_count) {

  printf("TESTING atomic_ptrs: %d threads, thread %s append with %d bufferSize and %d total counts\n", NumThreads, TS ? "SAFE" : "UNSAFE", buffer_capacity, total_count);
  omp_lock_t writelock;
  omp_init_lock(&writelock);
  omp_lock_t writelock2;
  omp_init_lock(&writelock2);
  omp_lock_t writelock3;
  omp_init_lock(&writelock3);


  printf("TEST: save full buffers and process at end: ");
  std::vector<std::unique_ptr<bliss::io::Buffer<TS> > > full;

  std::vector<int> gold;
  std::vector<int> stored;

  int success = 0;
  int failure = 0;
  int swap = 0;
  int i = 0;

  std::atomic<bliss::io::Buffer<TS>* > ptr(new bliss::io::Buffer<TS>(buffer_capacity));                            // ensure atomicity
  ptr.load()->unblock();

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
//      if (!buf->is_reading()) {
//        fprintf(stdout, "FAIL atomic batched proc: at this point the buffer should be in read state.\n");
//        fflush(stdout);
//        std::cout << "buffer: " << *(buf) << std::endl << std::flush;
//
//      }

      // swap in a new one.
      bliss::io::Buffer<TS>* new_ptr = new bliss::io::Buffer<TS>(buffer_capacity);  // manage new buffer
	    new_ptr->unblock();

      bliss::io::Buffer<TS>* old_ptr = nullptr;

      old_ptr = ptr.exchange(new_ptr);
#pragma omp flush(ptr)
//
      // save the old buffer



      // this is showing a possible spurious wakeup...
      int oldsize = old_ptr->getSize() / sizeof(int);
      if (oldsize != buffer_capacity / sizeof(int)) {
        fprintf(stdout, "FAIL atomic DID NOT GET 2047 elements 1. local swap = %d, i = %d\n", swap, i);
        std::cout << "   atomic old buf: " << *(old_ptr) << std::endl
            << "   atomic new buf: " << *(ptr.load()) << std::endl << std::flush;
      }


      omp_set_lock(&writelock);
      full.push_back(std::move(std::unique_ptr<bliss::io::Buffer<TS> >(old_ptr)));
      omp_unset_lock(&writelock);

      ++swap;
    }

  }
  //printf("LAST BUFFER 1\n");

  ptr.load()->block_and_flush();
  int last = ptr.load()->getSize();
  if (last == (buffer_capacity / sizeof(int))) {
    ++swap;
  }

  auto b = ptr.exchange(nullptr);
  full.push_back(std::move(std::unique_ptr<bliss::io::Buffer<TS> >(b)));

//  printf("DEBUG: atomic 1 success %d, failure %d, swap %d, total %d, full count %ld\n", success, failure, swap, total_count, full.size());
//  std::cout << " buffer: " << *(ptr.load()) << std::endl << std::flush;

  for (int i = 0; i < full.size(); ++i) {
    stored.insert(stored.end(), full.at(i)->operator int*(), full.at(i)->operator int*() + full.at(i)->getSize() / sizeof(int));
  }
  int stored_count = stored.size();


  if (success == 0 || swap != full.size() - 1  || swap != success / (buffer_capacity / sizeof(int)) || success != stored_count)
    printf("FAIL atomic: (actual/expected)  success (%d/%d), failure (%d/?), swap(%d,%ld/%ld), last buf size %d, content match? %s.\n", stored_count, success, failure, swap, full.size(), success / (buffer_capacity / sizeof(int)), last, compareUnorderedSequences(stored.begin(), gold.begin(), stored_count) ? "same" : "diff");

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

  b = ptr.exchange(new bliss::io::Buffer<TS>(buffer_capacity));  // old pointer was managed by unique ptr.
  ptr.load()->unblock();

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
//    	// 0x2 before all the memcpy are completed.  thus we could get is_reading() failed while result is 0x2, and also observe a large
//    	// number of writes after result is set to 0x2 (and before that the flush bit is set)
//    	// this is a theory.
//
//      if (!(buf->is_reading())) {
//        fprintf(stdout, "FAIL atomic incremental proc: at this point the buffer should be in read state.  res= %d\n", res);
//        fflush(stdout);
//        std::cout << "buffer: " << *(buf) << std::endl << std::flush;
//      }

      bliss::io::Buffer<TS>* new_ptr = new bliss::io::Buffer<TS>(buffer_capacity);  // manage new buffer
      //std::cout << "   new buf before assing: " << *(new_ptr) << std::endl <<  std::flush;

      bliss::io::Buffer<TS>* old_ptr = nullptr;


      new_ptr->unblock();


      old_ptr = ptr.exchange(new_ptr);             //
#pragma omp flush(ptr)
      // save the old buffer
      ++swap;


      if (old_ptr != nullptr) {
        int oldsize = old_ptr->getSize() / sizeof(int);
//        int newsize = buf_ptr.load()->getSize() / sizeof(int);
        if (oldsize != buffer_capacity / sizeof(int) || !(old_ptr->is_reading())) {
          fprintf(stdout, "FAIL: atomic DID NOT GET 2047 elements 2. local swap = %d, i = %d\n", swap, i);
          std::cout << "   old buf: " << *(old_ptr) << std::endl <<  std::flush;
        }
        success2 += oldsize;

        omp_set_lock(&writelock);
        stored.insert(stored.end(), old_ptr->operator int*(), old_ptr->operator int*() + oldsize);
        full.push_back(std::move(std::unique_ptr<bliss::io::Buffer<TS> >(old_ptr)));
        omp_unset_lock(&writelock);

      }


    }

  }



  //printf("LAST BUFFER 2\n");
  ptr.load()->block_and_flush();
  last = ptr.load()->getSize();
  if (last == (buffer_capacity / sizeof(int))) {
    ++swap;
  }

  stored_count = stored.size();

  //printf("DEBUG: atomic before last buffer (actual/expected)  success (%d,%d/%d), failure (%d/?), swap(%d/%ld). content match? %s\n", stored_count, success2, success, failure, swap, success / (buffer_capacity / sizeof(int)), compareUnorderedSequences(stored.begin(), gold.begin(), stored_count) ? "same" : "diff");


  // compare unordered buffer content.
  stored.insert(stored.end(), ptr.load()->operator int*(), ptr.load()->operator int*() + ptr.load()->getSize() / sizeof(int));
  full.push_back(std::move(std::unique_ptr<bliss::io::Buffer<TS> >(ptr.load())));

  stored_count = stored.size();
  success2 += ptr.load()->getSize() / sizeof(int);


  //printf("DEBUG: atomic after last buffer (actual/expected)  success (%d,%d/%d), failure (%d/?), swap(%d/%ld), final buf size %d, content match? %s\n", stored_count, success2, success, failure, swap, success / (buffer_capacity / sizeof(int)), last, compareUnorderedSequences(stored.begin(), gold.begin(), stored_count) ? "same" : "diff");

  if ( success == 0 || swap != success / (buffer_capacity / sizeof(int)) || success != stored_count)
    printf("FAIL atomic: (actual/expected)  success (%d,%d/%d), failure (%d/?), last %d, swap(%d/%ld). content match? %s\n", stored_count, success2, success, failure, last, swap, success / (buffer_capacity / sizeof(int)), compareUnorderedSequences(stored.begin(), gold.begin(), stored_count) ? "same" : "diff");
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

  for (int k = 0; k < iterations; ++k) {

      buf.clear();
      buf.unblock();
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



  // multiple buffer swap test.

  ////////////// timing.  the insert before this is to warm up.
  testAppendMultipleBuffersAtomicPtrs<bliss::concurrent::THREAD_UNSAFE, 1>(8191, 1000000);

  testAppendMultipleBuffersAtomicPtrs<bliss::concurrent::THREAD_SAFE, 1>(8191, 1000000);
  testAppendMultipleBuffersAtomicPtrs<bliss::concurrent::THREAD_SAFE, 2>(8191, 1000000);
  testAppendMultipleBuffersAtomicPtrs<bliss::concurrent::THREAD_SAFE, 3>(8191, 1000000);
  testAppendMultipleBuffersAtomicPtrs<bliss::concurrent::THREAD_SAFE, 4>(8191, 1000000);
  testAppendMultipleBuffersAtomicPtrs<bliss::concurrent::THREAD_SAFE, 5>(8191, 1000000);
  testAppendMultipleBuffersAtomicPtrs<bliss::concurrent::THREAD_SAFE, 6>(8191, 1000000);
  testAppendMultipleBuffersAtomicPtrs<bliss::concurrent::THREAD_SAFE, 7>(8191, 1000000);
  testAppendMultipleBuffersAtomicPtrs<bliss::concurrent::THREAD_SAFE, 8>(8191, 1000000);

//  testAppendMultipleBuffers<bliss::concurrent::THREAD_UNSAFE, 1>(8191, 1000000);
//
//  testAppendMultipleBuffers<bliss::concurrent::THREAD_SAFE, 1>(8191, 1000000);
//  testAppendMultipleBuffers<bliss::concurrent::THREAD_SAFE, 2>(8191, 1000000);
//  testAppendMultipleBuffers<bliss::concurrent::THREAD_SAFE, 3>(8191, 1000000);
//  testAppendMultipleBuffers<bliss::concurrent::THREAD_SAFE, 4>(8191, 1000000);
//  testAppendMultipleBuffers<bliss::concurrent::THREAD_SAFE, 5>(8191, 1000000);
//  testAppendMultipleBuffers<bliss::concurrent::THREAD_SAFE, 6>(8191, 1000000);
//  testAppendMultipleBuffers<bliss::concurrent::THREAD_SAFE, 7>(8191, 1000000);
//  testAppendMultipleBuffers<bliss::concurrent::THREAD_SAFE, 8>(8191, 1000000);


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
}
