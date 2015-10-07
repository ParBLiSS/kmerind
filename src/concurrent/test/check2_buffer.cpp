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
#include <cstdio>
#include <chrono>
#include <iterator>  // for ostream_iterator
#include <iostream>   // for cout
#include <sstream>
#include <xmmintrin.h>

#include "utils/iterator_test_utils.hpp"
#include "utils/logging.h"

#include <deque>
#include <list>

template<bliss::concurrent::LockType TS, size_t CAP, size_t MDSize>
void append(const int nthreads, bliss::io::Buffer<TS, CAP, MDSize>& buf,
            const unsigned int start, const unsigned int end, size_t& success,
            size_t& failure, size_t& swap, std::vector<std::vector<int> >& gold)
{

  unsigned int i;
  size_t lsuccess = 0;
  size_t lfailure = 0;
  size_t lswap = 0;
  // TODO: OMP default(none) should be compatible with shared, but clang complains.
#pragma omp parallel for OMP_SHARE_DEFAULT num_threads(nthreads) shared(buf, gold, stdout) private(i) reduction(+: lsuccess, lfailure, lswap)
  for (i = start; i < end; ++i)
  {
    int data = static_cast<int>(i);
    int result = buf.append(&data, 1);

    if ((result & 0x1) > 0)
    {
      ++lsuccess;

      gold[omp_get_thread_num()].push_back(data);

    }
    else
    {
      ++lfailure;
    }

    if ((result & 0x2) > 0)
    {
      if (!buf.is_read_only())
      {
        ERRORF(
            "FAIL append: at this point the buffer should be in read state.");

      }

      ++lswap;
    }
  }

  success += lsuccess;
  failure += lfailure;
  swap += lswap;
  //INFOF("DEBUG: threads %d (actual,added/expected) success (%d/%d), failure (%d/%d), swap(%d)", nthreads, success, end - start, failure, end-start, swap);

}

template<bliss::concurrent::LockType TS, size_t CAP, size_t MDSize,
    int NumThreads = 1>
void appendTest()
{
  static_assert(NumThreads > 0, "instantiated with NumThreads < 1");
  static_assert(TS != bliss::concurrent::LockType::NONE || NumThreads == 1, "instantiated with Thread Unsafe version and NumThreads != 1");

  INFOF("TESTING operations on locktype %d buffer", static_cast<int>(TS));

  // create a buffer.
  bliss::io::Buffer<TS, CAP, MDSize> b1;
  b1.clear_and_unblock_writes();

  size_t nelems = CAP / sizeof(int);
  size_t remainder = CAP % sizeof(int);

  size_t success = 0;
  size_t failure = 0;
  size_t swap = 0;
  std::vector<std::vector<int> > gold(NumThreads);
  std::vector<int> ggold;

  INFOF("TEST insert under capacity: ");
  append(NumThreads, b1, 0, nelems / 2, success, failure, swap, gold);

  for (int i = 0; i < NumThreads; ++i)
  {
    ggold.insert(ggold.end(), gold[i].begin(), gold[i].end());
  }
  assert(ggold.size() == nelems/2);

  if (success == 0 || (success != nelems / 2) || failure != 0 || swap != 0)
  {
    ERRORF(
        "FAIL: (actual,added/expected) success (%lu,%lu/%lu), failure (%lu,%lu/%d), swap(%lu,%lu/%d)",
        success, success, nelems/2, failure, failure, 0, swap, swap, 0);
  }
  else
  {
    // compare unordered buffer content.
    INFOF("success : %lu, gold size %lu", success, ggold.size());

    if (compareUnorderedSequences(b1.template begin<int>(), ggold.begin(), success))
    {
      INFOF("PASS success %lu failure %lu swap %lu", success, failure, swap);
    }
    else
    {
      FATALF("FAIL: content not matching");
    }
  }

  size_t success2 = 0;
  size_t failure2 = 0;
  size_t swap2 = 0;
  INFOF("TEST insert over capacity: ");

  append(NumThreads, b1, nelems / 2, nelems * 2, success2, failure2, swap2,
         gold);

  assert(success2 > 0);
  success += success2;
  failure += failure2;
  swap += swap2;

  ggold.clear();
  for (int i = 0; i < NumThreads; ++i)
  {
    ggold.insert(ggold.end(), gold[i].begin(), gold[i].end());
  }
  assert(ggold.size() == nelems);

  if (success == 0 || (success != nelems) || failure != nelems || swap != 1)
  {
    ERRORF(
        "FAIL: (actual,added/expected) success (%lu,%lu/%lu), failure (%lu,%lu/%lu), swap(%lu,%lu/%d)",
        success, success2, nelems, failure, failure2, nelems, swap, swap2, 1);
  }
  else
  {
    INFOF("success : %lu, gold size %lu", success, ggold.size());


    // compare unordered buffer content.
    if (compareUnorderedSequences(b1.template begin<int>(), ggold.begin(), success))
    {
      INFOF("PASS success %lu failure %lu swap %lu", success, failure, swap);
    }
    else
    {
      FATALF("FAIL: content not matching");
    }
  }

  INFOF("TEST clear: ");
  b1.clear_and_block_writes();
  if (b1.getSize() != 0)
    ERRORF("FAIL: NOT empty:  Size: %ld", b1.getSize());
  else
    INFOF("PASS");

  success = 0;
  failure = 0;
  swap = 0;
  for (int i = 0; i < NumThreads; ++i)
  {
    gold[i].clear();
  }
  ggold.clear();
  b1.unblock_writes();

  INFOF("TEST insert AT capacity: ");

  append(NumThreads, b1, 0, nelems, success, failure, swap, gold);

  size_t swap_exp = (remainder > 0 ? 0 : 1);

  for (int i = 0; i < NumThreads; ++i)
  {
    ggold.insert(ggold.end(), gold[i].begin(), gold[i].end());
  }


  if (success == 0 || (success != nelems) || failure != 0 || swap != swap_exp)
  {
    ERRORF(
        "FAIL: (actual/expected) success (%lu/%lu), failure (%lu/%d), swap(%lu/%lu)",
        success, nelems, failure, 0, swap, swap_exp);
  }
  else
  {

    // compare unordered buffer content.
    if (compareUnorderedSequences(b1.template begin<int>(), ggold.begin(), success))
    {
      INFOF("PASS success %lu failure %lu swap %lu", success, failure, swap);
    }
    else
    {
      FATALF("FAIL: content not matching");
    }
  }

  b1.clear_and_unblock_writes();

  success = 0;
  failure = 0;
  swap = 0;
  for (int i = 0; i < NumThreads; ++i)
  {
    gold[i].clear();
  }
  ggold.clear();

  INFOF("TEST insert JUST OVER capacity: ");

  append(NumThreads, b1, 0, nelems + NumThreads, success, failure, swap, gold);

  for (int i = 0; i < NumThreads; ++i)
  {
    ggold.insert(ggold.end(), gold[i].begin(), gold[i].end());
  }

  if (success == 0 || (success != nelems) || failure != NumThreads || swap != 1)
  {
    ERRORF(
        "FAIL: (actual/expected) success (%lu/%lu), failure (%lu/%d), swap(%lu/%d)",
        success, nelems, failure, NumThreads, swap, 1);
  }
  else
  {

    // compare unordered buffer content.
    if (compareUnorderedSequences(b1.template begin<int>(), ggold.begin(), success))
    {
      INFOF("PASS success %lu failure %lu swap %lu", success, failure, swap);
    }
    else
    {
      FATALF("FAIL: content not matching");
    }
  }

  INFOF("TEST blocked buffer: ");
  b1.clear_and_block_writes();
  b1.block_and_flush();

  success = 0;
  failure = 0;
  swap = 0;
  for (int i = 0; i < NumThreads; ++i)
  {
    gold[i].clear();
  }
  ggold.clear();

  append(NumThreads, b1, 0, nelems, success, failure, swap, gold);

  for (int i = 0; i < NumThreads; ++i)
  {
    ggold.insert(ggold.end(), gold[i].begin(), gold[i].end());
  }


  if ((success != 0) || failure != nelems || swap != 0)
  {
    ERRORF(
        "FAIL: (actual/expected) success (%lu/%d), failure (%lu/%lu), swap(%lu/%d)",
        success, 0, failure, nelems, swap, 0);
  }
  else
  {
    // compare unordered buffer content.
    if (compareUnorderedSequences(b1.template begin<int>(), ggold.begin(), success))
    {
      INFOF("PASS success %lu failure %lu swap %lu", success, failure, swap);
    }
    else
    {
      FATALF("FAIL: content not matching");
    }
  }

  INFOF("TEST unblock buffer: ");

  b1.clear_and_unblock_writes();

  success = 0;
  failure = 0;
  swap = 0;
  for (int i = 0; i < NumThreads; ++i)
  {
    gold[i].clear();
  }
  ggold.clear();

  append(NumThreads, b1, 0, nelems, success, failure, swap, gold);

  for (int i = 0; i < NumThreads; ++i)
  {
    ggold.insert(ggold.end(), gold[i].begin(), gold[i].end());
  }


  if (success == 0 || (success != nelems) || failure != 0 || swap != swap_exp)
  {
    ERRORF(
        "FAIL: (actual/expected) success (%lu/%lu), failure (%lu/%d), swap(%lu/%lu)",
        success, nelems, failure, 0, swap, swap_exp);
  }
  else
  {
    // compare unordered buffer content.
    if (compareUnorderedSequences(b1.template begin<int>(), ggold.begin(), success))
    {
      INFOF("PASS success %lu failure %lu swap %lu", success, failure, swap);
    }
    else
    {
      FATALF("FAIL: content not matching");
    }
  }

}

template<bliss::concurrent::LockType TS, size_t CAP, size_t MDSize,
    int NumThreads>
void testAppendMultipleBuffersAtomicPtrs(const size_t total_count)
{

  INFOF(
      "TESTING atomic_ptrs: %d threads, locktype %d append with %ld bufferSize and %lu total counts",
      NumThreads, static_cast<int>(TS), CAP, total_count);

  constexpr int elSize = sizeof(int);
  constexpr size_t capInEl = CAP / sizeof(int);

  INFOF("TEST: save full buffers and process at end: ");
  std::vector<std::vector<std::unique_ptr<bliss::io::Buffer<TS, CAP, MDSize> > > > full(NumThreads);

  std::vector<std::vector<int> > gold(NumThreads);
  std::vector<std::vector<int> > lstored(NumThreads);
  std::vector<int> stored;
  std::vector<int> ggold;

  size_t success = 0;
  size_t failure = 0;
  size_t swap = 0;

  std::atomic<bliss::io::Buffer<TS, CAP, MDSize>*> ptr(
      new bliss::io::Buffer<TS, CAP, MDSize>());             // ensure atomicity
  ptr.load(std::memory_order_relaxed)->unblock_writes();

#pragma omp parallel for num_threads(NumThreads) OMP_SHARE_DEFAULT shared(ptr, full, gold, stderr, stdout, std::cout) reduction(+:success, failure, swap)
  for (size_t i = 0; i < total_count; ++i)
  {

    int data = static_cast<int>(i);
    //auto buf = ptr.load();

    std::atomic_thread_fence(std::memory_order_seq_cst);

    unsigned int result = ptr.load(std::memory_order_relaxed)->append(&data, 1);

    if (result & 0x1)
    {
      ++success;

      gold[omp_get_thread_num()].push_back(data);

    }
    else
    {
      ++failure;
      _mm_pause();  // slow it down a little.
    }

    if (result & 0x2)
    {

//      if (result & 0x4) INFO( "SWAPPING: " << *(buf) );
//
//      if (!buf->is_read_only()) {
//        FATALF("FAIL atomic batched proc: at this point the buffer should be in read state.");
//        fflush(stdout);
//        INFO( "buffer: " << *(buf) );
//
//      }

      // swap in a new one.
      bliss::io::Buffer<TS, CAP, MDSize>* new_ptr = new bliss::io::Buffer<TS,
          CAP, MDSize>();  // manage new buffer
      new_ptr->unblock_writes();

      bliss::io::Buffer<TS, CAP, MDSize>* old_ptr = nullptr;

      old_ptr = ptr.exchange(new_ptr, std::memory_order_acq_rel);
#pragma omp flush(ptr)
//
      // save the old buffer

      // this is showing a possible spurious wakeup...
      size_t oldsize = old_ptr->getSize() / elSize;
      if (oldsize != capInEl)
      {
        ERRORF(
            "FAIL 1 atomic DID NOT GET %lu elements in cap %lu bytes. got %lu in %lu bytes. local swap = %lu, i = %lu",
            capInEl, old_ptr->getCapacity(), oldsize, old_ptr->getSize(), swap,
            i);
        INFO(
            "   atomic old buf: " << *(old_ptr) << std::endl << "   atomic new buf: " << *(ptr.load(std::memory_order_relaxed)));
      }

      full[omp_get_thread_num()].push_back(
          std::move(std::unique_ptr<bliss::io::Buffer<TS, CAP, MDSize> >(old_ptr)));

      ++swap;
    }

  }
  //INFOF("LAST BUFFER 1");

  ptr.load(std::memory_order_relaxed)->block_and_flush();
  auto last = ptr.load(std::memory_order_relaxed)->getSize();
  if (last == (capInEl))
  {
    ++swap;
  }

  auto b = ptr.exchange(nullptr, std::memory_order_acq_rel);
  full[0].push_back(
      std::move(std::unique_ptr<bliss::io::Buffer<TS, CAP, MDSize> >(b)));

//  INFOF("DEBUG: atomic 1 success %lu, failure %lu, swap %lu, total %lu, full count %ld", success, failure, swap, total_count, full.size());
//  INFO( " buffer: " << *(ptr.load()) );
  size_t fullsize = 0;
  for (size_t i = 0; i < full.size(); ++i)
  {
    fullsize += full[i].size();
    for (size_t j = 0; j < full[i].size(); ++j)
    {
      stored.insert(stored.end(), full[i][j]->template begin<int>(),
                    full[i][j]->template end<int>());
    }
  }
  auto stored_count = stored.size();

  for (int i = 0; i < NumThreads; ++i)
  {
    ggold.insert(ggold.end(), gold[i].begin(), gold[i].end());
  }

  if (success == 0 || swap != fullsize - 1 || swap != success / (capInEl)
      || success != stored_count)
  {
    ERRORF(
        "FAIL atomic: (actual/expected)  success (%lu/%lu), failure (%lu/?), last %lu, swap(%lu,%lu/%lu), last buf size %lu, content match? %s.",
        stored_count,
        success,
        failure,
        last,
        swap,
        fullsize,
        success / (capInEl),
        last,
        compareUnorderedSequences(stored.begin(), ggold.begin(), stored_count) ? "same" : "diff");
  }
  else
  {

    if (compareUnorderedSequences(stored.begin(), ggold.begin(), stored_count))
    {
      INFOF("PASS: atomic success %lu, failure %lu, swap %lu, total %lu", success,
            failure, swap, total_count);
    }
    else
    {
      FATALF(
          "FAIL: atomic success %lu, failure %lu, swap %lu, total %lu, content not matching",
          success, failure, swap, total_count);
    }
  }

  INFOF("TEST: process full buffers along the way (SAVE IN VECTOR): ");

  for (int i = 0; i < NumThreads; ++i)
  {
    gold[i].clear();
    full[i].clear();
  }
  ggold.clear();
  stored.clear();

  success = 0;
  failure = 0;
  swap = 0;

  size_t success2 = 0;

  b = ptr.exchange(new bliss::io::Buffer<TS, CAP, MDSize>(), std::memory_order_acq_rel); // old pointer was managed by unique ptr.
  ptr.load(std::memory_order_relaxed)->unblock_writes();

#pragma omp parallel for num_threads(NumThreads) OMP_SHARE_DEFAULT shared(ptr, gold, lstored, stderr, stdout, std::cout) reduction(+:success, failure, swap, success2)
  for (size_t i = 0; i < total_count; ++i)
  {

    std::atomic_thread_fence(std::memory_order_seq_cst);

    int data = static_cast<int>(i);
    //auto buf = ptr.load();
    std::atomic_thread_fence(std::memory_order_seq_cst);

    int res = ptr.load(std::memory_order_relaxed)->append(&data, 1);

    if (res & 0x1)
    {
      ++success;

      gold[omp_get_thread_num()].push_back(data);

    }
    else
    {
      ++failure;
//       _mm_pause();  // slow it down a little.
    }

    if (res & 0x2)
    {

//      if (res & 0x4) INFO( "SWAPPING: " << *(buf) );
//    	// TODO: issue here:  if a large number of threads call append, and most of them are rescheduled, so that we reach calc
//    	// of pointer for a large number of threads in progress.  Then we could have the "just overflowing" thread executing and returning
//    	// 0x2 before all the memcpy are completed.  thus we could get is_read_only() failed while result is 0x2, and also observe a large
//    	// number of writes after result is set to 0x2 (and before that the flush bit is set)
//    	// this is a theory.
//
//      if (!(buf->is_read_only())) {
//        FATALF("FAIL atomic incremental proc: at this point the buffer should be in read state.  res= %lu", res);
//        fflush(stdout);
//        INFO( "buffer: " << *(buf) );
//      }

      bliss::io::Buffer<TS, CAP, MDSize>* new_ptr = new bliss::io::Buffer<TS,
          CAP, MDSize>();  // manage new buffer
      //INFO( "   new buf before assing: " << *(new_ptr) );

      bliss::io::Buffer<TS, CAP, MDSize>* old_ptr = nullptr;

      new_ptr->unblock_writes();

      old_ptr = ptr.exchange(new_ptr, std::memory_order_acq_rel);             //
#pragma omp flush(ptr)
      // save the old buffer

      if (old_ptr != nullptr)
      {
        ++swap;
        size_t oldsize = old_ptr->getSize() / elSize;
        //        int newsize = buf_ptr.load()->getSize() / elSize;
        if (oldsize != capInEl || !(old_ptr->is_read_only()))
        {
          ERRORF(
              "FAIL 2 atomic DID NOT GET %lu elements. actual %lu (%lu). local swap = %lu, i = %lu",
              capInEl, oldsize, old_ptr->getSize(), swap, i);
          INFO("   old buf: " << *(old_ptr));
        }
        success2 += oldsize;

        int* d = old_ptr->template begin<int>();

        lstored[omp_get_thread_num()].insert(
            lstored[omp_get_thread_num()].end(), d,
            (d + old_ptr->getSize() / sizeof(int)));

      }

    }

  }

  //INFOF("LAST BUFFER 2");
  ptr.load(std::memory_order_relaxed)->block_and_flush();
  last = ptr.load(std::memory_order_relaxed)->getSize();
  if (last == (capInEl))
  {
    ++swap;
  }

  //INFOF("DEBUG: atomic before last buffer (actual/expected)  success (%lu,%lu/%lu), failure (%lu/?), swap(%lu/%ld). content match? %s", stored_count, success2, success, failure, swap, success / (capInEl), compareUnorderedSequences(stored.begin(), gold.begin(), stored_count) ? "same" : "diff");

  for (int i = 0; i < NumThreads; ++i)
  {
    stored.insert(stored.end(), lstored[i].begin(), lstored[i].end());
    ggold.insert(ggold.end(), gold[i].begin(), gold[i].end());
  }

  // compare unordered buffer content.
  stored.insert(stored.end(), ptr.load(std::memory_order_relaxed)->template begin<int>(),
                ptr.load(std::memory_order_relaxed)->template end<int>());

  stored_count = stored.size();
  success2 += ptr.load(std::memory_order_relaxed)->getSize() / elSize;

  //INFOF("DEBUG: atomic after last buffer (actual/expected)  success (%lu,%lu/%lu), failure (%lu/?), swap(%lu/%ld), final buf size %lu, content match? %s", stored_count, success2, success, failure, swap, success / (capInEl), last, compareUnorderedSequences(stored.begin(), gold.begin(), stored_count) ? "same" : "diff");

  if (success == 0 || swap != success / (capInEl) || success != stored_count)
  {
    ERRORF(
        "FAIL atomic: (actual/expected)  success (%lu,%lu/%lu), failure (%lu/?), last %lu, swap(%lu/%lu). content match? %s",
        stored_count,
        success2,
        success,
        failure,
        last,
        swap,
        success / (capInEl),
        compareUnorderedSequences(stored.begin(), ggold.begin(), stored_count) ? "same" : "diff");
  }
  else
  {
    if (compareUnorderedSequences(stored.begin(), ggold.begin(), stored_count))
    {
      INFOF("PASS: atomic success %lu, failure %lu, swap %lu, total %lu", success,
            failure, swap, total_count);
    }
    else
    {
      FATALF(
          "FAIL: atomic success %lu, failure %lu, swap %lu, total %lu, content not matching",
          success, failure, swap, total_count);
    }
  }

  ptr.exchange(nullptr, std::memory_order_acq_rel);
  stored.clear();

}

template<bliss::concurrent::LockType TS, size_t CAP, size_t MDSize,
    int NumThreads>
void stressTestAppendMultipleBuffersAtomicPtrs(const size_t total_count)
{

  INFOF(
      "TESTING atomic_ptrs: stress %d threads, locktype %d append with %ld bufferSize and %lu total counts",
      NumThreads, static_cast<int>(TS), CAP, total_count);

  constexpr size_t elSize = sizeof(size_t);
  constexpr size_t capInEl = CAP / elSize;

  size_t success = 0;
  size_t failure = 0;
  size_t failure2 = 0;
  size_t failure3 = 0;
  size_t swap = 0;
  size_t i = 0;

  std::atomic<bliss::io::Buffer<TS, CAP, MDSize>*> ptr(
      new bliss::io::Buffer<TS, CAP, MDSize>());             // ensure atomicity
  ptr.load(std::memory_order_relaxed)->clear_and_unblock_writes();

  std::deque<bliss::io::Buffer<TS, CAP, MDSize>*> full;

#pragma omp parallel for num_threads(NumThreads) OMP_SHARE_DEFAULT shared(ptr, stdout, full) private(i) reduction(+:success, failure, swap, failure2, failure3)
  for (i = 0; i < total_count; ++i)
  {

    size_t data = i;
    void* out = nullptr;
    auto localptr = ptr.load(std::memory_order_consume);
    unsigned int result = localptr->append(&data, 1, &out);

    if ((result & 0x1) > 0)
    {
      ++success;

      if (out == nullptr)
      {
        FATALF("ERROR: successful append but no pointer returned.");
        fflush(stdout);
        ++failure2;
      }
      else
      {
        size_t od = *(reinterpret_cast<size_t*>(out));
        if (od != data)
        {
          FATALF(
              "ERROR: thread %d successful append but value is not correctly stored: expected %lu, actual %lu. insert buf %p, curr buffer %p, insert dataptr %p, data ptr %p, curr data ptr %p, returned %p, offset %ld",
              omp_get_thread_num(), data, od, localptr, ptr.load(std::memory_order_relaxed), localptr->template begin<char>(),
              localptr->template begin<char>(), ptr.load(std::memory_order_relaxed)->template begin<char>(),
              (char*)out, (char*)out - (localptr->template begin<char>()));
          fflush(stdout);
          ++failure3;
        }
      }

    }
    else
    {
      ++failure;
      _mm_pause();  // slow it down a little.
    }

    if ((result & 0x2) > 0)
    {

      // swap in a new one.
      bliss::io::Buffer<TS, CAP, MDSize>* new_ptr = new bliss::io::Buffer<TS,
          CAP, MDSize>();  // manage new buffer
      new_ptr->clear_and_unblock_writes();

      bliss::io::Buffer<TS, CAP, MDSize>* old_ptr = localptr;
      bool exchanged = ptr.compare_exchange_strong(localptr, new_ptr,
                                                   std::memory_order_acq_rel);
//#pragma omp flush(ptr)
      //INFOF("SWAP: old buf %p, new buf %p", old_ptr, ptr.load());

      // save the old buffer

      if (exchanged)
      {
        //if (omp_get_num_threads() > 1) INFOF("INFO: exchanged. thread %d/%d,  old %p, new %p, ptr %p", omp_get_thread_num(), omp_get_num_threads(), old_ptr, new_ptr, ptr.load(std::memory_order_relaxed));
        // this is showing a possible spurious wakeup...
        size_t oldsize = old_ptr ? old_ptr->getSize() / elSize : 0;
        if (oldsize != capInEl)
        {
          ERRORF(
              "FAIL 3 thread %d/%d atomic DID NOT GET %lu elements, actual %lu. local swap = %lu, i = %lu. oldbuf %p, newbuf %p",
              omp_get_thread_num(), omp_get_num_threads(), capInEl, oldsize,
              swap, i, old_ptr, ptr.load(std::memory_order_relaxed));
        }

//        delete old_ptr;
        if (full.size() > NumThreads * NumThreads)
        {   // picked t^2 arbitrarily
          if (full.front())
            delete full.front();
          full.pop_front();
        }
        full.push_back(old_ptr);
        ++swap;
      }
      else
      {
        FATALF(
            "FAIL: thread %d/%d atomic buffer ptr swap failed, orig %p, new %p, curr %p",
            omp_get_thread_num(), omp_get_num_threads(), old_ptr, new_ptr,
            ptr.load(std::memory_order_relaxed));
        delete new_ptr;
      }
    }

  }
  //INFOF("LAST BUFFER 1");

  ptr.load(std::memory_order_relaxed)->block_and_flush();
  size_t last = ptr.load(std::memory_order_relaxed)->getSize() / elSize;
  if (last == (capInEl))
  {
    ++swap;
  }

  auto b = ptr.exchange(nullptr, std::memory_order_acq_rel);
  delete b;

  if (failure2 > 0 || failure3 > 0)
  {
    ERRORF(
        "FAIL: bad inserts present: count of nullptr returned %lu, count of bad value %lu",
        failure2, failure3);
  }

  if (success == 0 || swap != success / (capInEl))
  {
    ERRORF(
        "FAIL atomic: success (%lu), failure (%lu/%lu/%lu), swap(%lu/%ld), last buf size %lu.",
        success, failure, failure2, failure3, swap, success / (capInEl), last);
  }
  else
    INFOF("PASS: atomic success %lu, failure %lu/%lu/%lu, swap %lu, total %lu",
          success, failure, failure2, failure3, swap, total_count);

}
//
//template<bliss::concurrent::LockType TS, size_t CAP, size_t MDSize>
//struct AllocatingBufferPool {
//    using BufferType = bliss::io::Buffer<TS, CAP, MDSize>;
//    using BufferPtrType = std::shared_ptr< BufferType >;
//
//    BufferPtrType buf;  // need to be replaced.
//
//    BufferPtrType acquireBuffer() {
//      buf = BufferPtrType(new BufferType());  // new ptr assigned to buf.  buf's previous content is lost in this scope but exists in scope of other threads.
//      buf->clear_and_unblock_writes();
//      return BufferPtrType(buf);  // copy and return;
//    }
//
//    BufferPtrType getBuffer() {
//      return BufferPtrType(buf);  // copy and return
//    }
//};
//

//
//template<bliss::concurrent::LockType TS, size_t CAP, size_t MDSize, int NumThreads>
//void stressTestAppendMultipleBuffersSharedPtrs(const size_t total_count) {
//
//  INFOF("TESTING shared_ptr: stress %d threads, locktype %d append with %ld bufferSize and %lu total counts", NumThreads, static_cast<int>(TS), CAP, total_count);
//
//  constexpr size_t elSize = sizeof(size_t);
//  constexpr size_t capInEl = CAP / elSize;
//
//  int success = 0;
//  int failure = 0;
//  int failure2 = 0;
//  int failure3 = 0;
//  int swap = 0;
//  size_t i = 0;
//
//  AllocatingBufferPool<TS, CAP, MDSize> pool;
//  {
//    pool.acquireBuffer();
//  }
//
//#pragma omp parallel for num_threads(NumThreads) OMP_SHARE_DEFAULT shared(pool, stdout) private(i) reduction(+:success, failure, swap, failure2, failure3)
//  for (i = 0; i < total_count; ++i) {
//
//    std::shared_ptr<bliss::io::Buffer<TS, CAP, MDSize>> ptr = pool.getBuffer();
//
//    size_t data = i;
//    void* out = nullptr;
//
//    unsigned int result = ptr->append(&data, 1, &out);
//
//    if (result & 0x1) {
//      ++success;
//
//      if (out == nullptr) {
//        FATALF("ERROR: successful append but no pointer returned.");
//        fflush(stdout);
//        ++failure2;
//      } else {
//        size_t od = *((size_t*)out);
//        if (od != data) {
//          FATALF("ERROR: thread %d successful append but value is not correctly stored: expected %lu, actual %lu. buffer %p data ptr %p, result ptr %p, offset %ld",
//                 omp_get_thread_num(), data, od, ptr.get(), ptr->template begin<char>(), (char*) out, (char*)out - ptr->template begin<char>());
//          fflush(stdout);
//          ++failure3;
//        }
//      }
//
//    } else {
//      ++failure;
//      _mm_pause();  // slow it down a little.
//    }
//
//    if (result & 0x2) {
//
//      // swap in a new one.
//
//      auto new_ptr = pool.acquireBuffer();
//
////      bliss::io::Buffer<TS, CAP, MDSize>* new_ptr = new bliss::io::Buffer<TS, CAP, MDSize>();  // manage new buffer
////      new_ptr->clear_and_unblock_writes();
////
////      bliss::io::Buffer<TS, CAP, MDSize>* old_ptr = gptr.get();
////
////      gptr.reset(new_ptr);
//#pragma omp flush(pool)
//      //INFOF("SWAP: old buf %p, new buf %p", old_ptr, ptr.load());
//
//      // save the old buffer
//
//      // this is showing a possible spurious wakeup...
//      int oldsize = ptr ? ptr->getSize() / elSize : 0;
//      if (oldsize != capInEl) {
//        FATALF("FAIL shared DID NOT GET 2047 elements 1. local swap = %d, i = %lu. oldbuf %p, newbuf %p", swap, i, ptr.get(), new_ptr.get());
//      }
//
//
////      delete old_ptr;
//      ++swap;
//    }
//
//  }
//  //INFOF("LAST BUFFER 1");
//  auto gptr = pool.getBuffer();
//  int last = 0;
//  if (gptr) {
//    gptr->block_and_flush();
//    last = gptr->getSize();
//    if (last == (capInEl)) {
//      ++swap;
//    }
//  }
//
////  auto b = gptr.get();
////  delete b;
//
//  if (failure2 > 0 || failure3 > 0 ) {
//    FATALF("FAIL: bad inserts present: count of nullptr returned %d, count of bad value %d", failure2, failure3);
//  }
//
//  if (success == 0 || swap != success / (capInEl))
//    FATALF("FAIL shared: success (%d), failure (%d/%d/%d), swap(%d/%ld), last buf size %d.", success, failure, failure2, failure3, swap, success / (capInEl), last);
//
//  else
//    INFOF("PASS: shared success %d, failure %d/%d/%d, swap %d, total %lu", success, failure, failure2, failure3, swap, total_count);
//
//
//}

int main(int argc, char** argv)
{

#if defined( BLISS_MUTEX)
  constexpr bliss::concurrent::LockType lt = bliss::concurrent::LockType::MUTEX;
//#elif defined(BLISS_SPINLOCK)
//  constexpr bliss::concurrent::LockType lt = bliss::concurrent::LockType::SPINLOCK;
//#elif defined(BLISS_LOCKFREE)
//  constexpr bliss::concurrent::LockType lt = bliss::concurrent::LockType::LOCKFREE;
#else //if defined(BLISS_LOCKFREE)
  constexpr bliss::concurrent::LockType lt =
      bliss::concurrent::LockType::LOCKFREE;
#endif

  // test append
  appendTest<bliss::concurrent::LockType::NONE, 8192, 0, 1>();

  appendTest<lt, 8192, 0, 1>();
  appendTest<lt, 8192, 0, 2>();
  appendTest<lt, 8192, 0, 3>();
  appendTest<lt, 8192, 0, 4>();
  appendTest<lt, 8192, 0, 5>();
  appendTest<lt, 8192, 0, 6>();
  appendTest<lt, 8192, 0, 7>();
  appendTest<lt, 8192, 0, 8>();

  // test append with buffer that is not multple of element size.
  appendTest<bliss::concurrent::LockType::NONE, 8191, 0, 1>();

  appendTest<lt, 8191, 0, 1>();
  appendTest<lt, 8191, 0, 2>();
  appendTest<lt, 8191, 0, 3>();
  appendTest<lt, 8191, 0, 4>();
  appendTest<lt, 8191, 0, 5>();
  appendTest<lt, 8191, 0, 6>();
  appendTest<lt, 8191, 0, 7>();
  appendTest<lt, 8191, 0, 8>();

  // multiple buffer swap test.

  ////////////// timing.  the insert before this is to warm up.
  testAppendMultipleBuffersAtomicPtrs<bliss::concurrent::LockType::NONE, 8191,
      0, 1>(1000000);

  testAppendMultipleBuffersAtomicPtrs<lt, 8191, 0, 1>(1000000);
  testAppendMultipleBuffersAtomicPtrs<lt, 8191, 0, 2>(1000000);
  testAppendMultipleBuffersAtomicPtrs<lt, 8191, 0, 3>(1000000);
  testAppendMultipleBuffersAtomicPtrs<lt, 8191, 0, 4>(1000000);
  testAppendMultipleBuffersAtomicPtrs<lt, 8191, 0, 5>(1000000);
  testAppendMultipleBuffersAtomicPtrs<lt, 8191, 0, 6>(1000000);
  testAppendMultipleBuffersAtomicPtrs<lt, 8191, 0, 7>(1000000);
  testAppendMultipleBuffersAtomicPtrs<lt, 8191, 0, 8>(1000000);

  // multiple buffer swap test.

  ////////////// timing.  the insert before this is to warm up.
  testAppendMultipleBuffersAtomicPtrs<bliss::concurrent::LockType::NONE, 8192,
      0, 1>(1000000);

  testAppendMultipleBuffersAtomicPtrs<lt, 8192, 0, 1>(1000000);
  testAppendMultipleBuffersAtomicPtrs<lt, 8192, 0, 2>(1000000);
  testAppendMultipleBuffersAtomicPtrs<lt, 8192, 0, 3>(1000000);
  testAppendMultipleBuffersAtomicPtrs<lt, 8192, 0, 4>(1000000);
  testAppendMultipleBuffersAtomicPtrs<lt, 8192, 0, 5>(1000000);
  testAppendMultipleBuffersAtomicPtrs<lt, 8192, 0, 6>(1000000);
  testAppendMultipleBuffersAtomicPtrs<lt, 8192, 0, 7>(1000000);
  testAppendMultipleBuffersAtomicPtrs<lt, 8192, 0, 8>(1000000);

  // no swapping.  - insert 10M elements into buffer of 100MB.
  stressTestAppendMultipleBuffersAtomicPtrs<bliss::concurrent::LockType::NONE,
      100000000, 0, 1>(10000000);

  stressTestAppendMultipleBuffersAtomicPtrs<lt, 100000000, 0, 1>(10000000);
  stressTestAppendMultipleBuffersAtomicPtrs<lt, 100000000, 0, 2>(10000000);
  stressTestAppendMultipleBuffersAtomicPtrs<lt, 100000000, 0, 3>(10000000);
  stressTestAppendMultipleBuffersAtomicPtrs<lt, 100000000, 0, 4>(10000000);
  stressTestAppendMultipleBuffersAtomicPtrs<lt, 100000000, 0, 5>(10000000);
  stressTestAppendMultipleBuffersAtomicPtrs<lt, 100000000, 0, 6>(10000000);
  stressTestAppendMultipleBuffersAtomicPtrs<lt, 100000000, 0, 7>(10000000);
  stressTestAppendMultipleBuffersAtomicPtrs<lt, 100000000, 0, 8>(10000000);

//  // swapping a lot.  - DATA RACE That should NOT be resolved by mutex.
//  stressTestAppendMultipleBuffersAtomicPtrs<bliss::concurrent::LockType::NONE, 2048, 0, 1>(1000000000);
//
//  stressTestAppendMultipleBuffersAtomicPtrs<lt, 2048, 0, 1>(1000000000);
//  stressTestAppendMultipleBuffersAtomicPtrs<lt, 2048, 0, 2>(1000000000);
//  stressTestAppendMultipleBuffersAtomicPtrs<lt, 2048, 0, 3>(1000000000);
//  stressTestAppendMultipleBuffersAtomicPtrs<lt, 2048, 0, 4>(1000000000);
//  stressTestAppendMultipleBuffersAtomicPtrs<lt, 2048, 0, 5>(1000000000);
//  stressTestAppendMultipleBuffersAtomicPtrs<lt, 2048, 0, 6>(1000000000);
//  stressTestAppendMultipleBuffersAtomicPtrs<lt, 2048, 0, 7>(1000000000);
//  stressTestAppendMultipleBuffersAtomicPtrs<lt, 2048, 0, 8>(1000000000);

}
