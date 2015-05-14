/**
 * @file		check_messagebuffers.cpp
 * @ingroup
 * @author	Tony Pan <tpan7@gatech.edu>
 * @brief
 * @details
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */


//#include <unistd.h>  // for usleep

#include "utils/logging.h"

#include "io/message_buffers.hpp"
#include "concurrent/referenced_object_pool.hpp"
#include "concurrent/mutexlock_queue.hpp"
//#include "concurrent/lockfree_queue.hpp"

#include "omp.h"
#include <cassert>
#include <chrono>
#include <cstdlib>   // for rand

int nelems = 100;
int bufferSize = 2048;
std::string data("this is a test.  this a test of the emergency broadcast system.  this is only a test. ");


template<typename BuffersType>
void testBuffers(BuffersType && buffers, bliss::concurrent::LockType poollt, bliss::concurrent::LockType bufferlt, int nthreads) {

  INFOF("*** TESTING Buffers lock %d buffer lock %d: ntargets = %lu, pool threads %d", poollt, bufferlt, buffers.getNumDests(), nthreads);


  INFOF("TEST append until full: ");
  typedef typename BuffersType::BufferPtrType BufferPtrType;
  bool op_suc = false;
  BufferPtrType ptr = nullptr;
  bliss::concurrent::ThreadSafeQueue<typename BuffersType::BufferPtrType, bliss::concurrent::LockType::MUTEX> fullBuffers;

  //INFOF("test string is \"%s\", length %lu", data.c_str(), data.length());
  int id = 0;

  int i;
  int success = 0;
  int failure = 0;
  int sswap = 0;
  int fswap = 0;


#pragma omp parallel for num_threads(nthreads) default(none) private(i, op_suc, ptr) firstprivate(id) shared(buffers, data, fullBuffers, nelems, bufferSize) reduction(+ : success, failure, sswap, fswap)
  for (i = 0; i < nelems; ++i) {
    //INFOF("insert %lu chars into %d", data.length(), id);
    uint32_t count = data.length();
    char* data_remain = nullptr;
    uint32_t count_remain = 0;

    std::tie(op_suc, ptr) = buffers.append(data.c_str(), count, data_remain, count_remain, id);

    if (op_suc) {
      ++success; // success
      if (ptr) {
        ++sswap;     // full buffer
        fullBuffers.waitAndPush(std::move(ptr));  // full buffer
      }
    } else {
      ++failure; // failure
      if (ptr) {
        ++fswap;     // full buffer
        fullBuffers.waitAndPush(std::move(ptr));  // full buffer
      }
    }

  }
//  if ((count + count2) != nelems) ERRORF("FAIL: number of successful inserts should be %d.  actual %d", nelems, count);
//  else if (count2 != 0) ERRORF("FAIL: number of failed insert overall should be 0. actual %d", count2);
//  else
//    if (!(count5 <= count && count <= (count5 + count3))) ERRORF("FAIL: number of successful inserts should be close to successful inserts without full buffers.");

  // compute 2 expected since append returns full buffer only on failed insert and if there are n inserts that brings it to just before full, then it depends on timing
  // as to when the buffer becomes full.  (other threads will fail on append until swap happens, but not return a full buffer.)
  int expectedFullMax = success / (bufferSize/data.length());
  int expectedFullMin = success / (bufferSize/data.length()) - nthreads;

  if (fullBuffers.getSize() != sswap + fswap) ERRORF("FAIL: number of full Buffers do not match: fullbuffer size %ld  full count %d + %d", fullBuffers.getSize(), sswap, fswap);
  // buffer at 23 entries (86 bytes each, 2048 bytes per buffer) will not show as full until the next iterator.
  else if (((fswap + sswap) > expectedFullMax) || ((sswap + fswap) < expectedFullMin)) ERRORF("FAIL: number of full Buffers is not right: %d+%d should be between %d and %d", sswap, fswap, expectedFullMin, expectedFullMax);
  //else if (count4 != 0) ERRORF("FAIL: number of failed insert due to no buffer should be 0. actual %d", count4);
  else INFOF("PASS");
  INFOF("Number of failed attempt to append to buffer is %d, success %d. full buffers size: %lu.  success swapped = %d, fail swapped = %d.", failure, success, fullBuffers.getSize(), sswap, fswap);



  INFOF("TEST release: ");
  int sswap2 = 0, fswap2 = 0, error = 0, over = 0;
//  INFOF("buffer ids 0 to %lu initially in use.", buffers.getSize() - 1);
//  INFOF("releasing: ");

  int iterations = (sswap + fswap) + 10;
#pragma omp parallel for num_threads(nthreads) default(none) private(id, op_suc, ptr) shared(buffers, data, fullBuffers, bufferSize, iterations) reduction(+ : sswap2, fswap2, error, over)
  for (id = 0; id < iterations; ++id) {
    try {
      std::tie(op_suc, ptr) = fullBuffers.tryPop();
      if (op_suc) {
//        INFOF("%d ", result.second);
        if (ptr) {
          buffers.releaseBuffer(std::move(ptr));
          ++sswap2;    // successful pop
        } else {
          ++fswap2;   // successful pop but no actual buffer to release
        }
      } else {

        ++over;     // failed pop.
      }
    } catch(const std::invalid_argument & e)
    {
      ERRORF("FAIL with %s", e.what());
      ++error;       // error during pop
    }
  }


  if (error != 0) ERRORF("FAIL: invalid argument exception during pop.  count = %d", error);
  else if (over != 10) ERRORF("FAIL: failed on pop %d times, expected 10", over);
  else if (fswap2 != 0) ERRORF("FAIL: succeeded in pop but not full buffer. %d", fswap2);
  else if (sswap2 != (sswap + fswap)) ERRORF("FAIL: successful pops. expected %d.  actual %d", (sswap+fswap), sswap2);
  else
    INFOF("PASS");



  buffers.reset();

  INFOF("TEST all operations together: ");
  int success3 = 0, failure3 = 0, sswap3 = 0, fswap3 = 0, bytes3 = 0, chars3 = 0;
  int success4 = 0, bytes4 = 0;
  int gbytes = 0;
  //INFOF("full buffer: ");
  id = 0;
#pragma omp parallel for num_threads(nthreads) default(none) private(i, op_suc, ptr) firstprivate(id) shared(buffers, data, nelems, bufferSize) reduction(+ : success3, success4, bytes4, failure3, sswap3, fswap3, bytes3, chars3)
  for (i = 0; i < nelems; ++i) {
    uint32_t count = data.length();
    char* data_remain = nullptr;
    uint32_t count_remain = 0;


    std::tie(op_suc, ptr) = buffers.append(data.c_str(), count, data_remain, count_remain, id);

    if (op_suc) {
      ++success3;

      if (ptr) {
        ++sswap3;
//        count7 = count1;  // save the number of successful inserts so far.
        bool updating = ptr->is_writing();
        if (updating) INFOF("  FULLBUFFER1: size %ld updating? %s, blocked? %s", ptr->getSize(), (updating ? "Y" : "N"), (ptr->is_read_only() ? "Y" : "N"));

        bytes3 += ptr->getSize();
        chars3 += strlen(ptr->operator char*());
        buffers.releaseBuffer(std::move(ptr));
      }
    } else {
      if (count_remain == count) {
        ++failure3;
      } else {
        ++success4;
        bytes4 += (count - count_remain);
      }

      if (ptr) {
        ++fswap3;
//        count7 = count1;  // save the number of successful inserts so far.
        bool updating = ptr->is_writing();
        if (updating) INFOF("  FULLBUFFER: size %ld blocked? %s", ptr->getSize(), (ptr->is_read_only() ? "Y" : "N"));

        bytes3 += ptr->getSize();
        chars3 += strlen(ptr->operator char*());

        buffers.releaseBuffer(std::move(ptr));
      }
    }
  }

  buffers.at(id)->block_and_flush();
  bool updating = buffers.at(id)->is_writing();
  if (updating) INFOF("  PreFLUSH: size %ld updating? %s, blocked? %s", buffers.at(id)->getSize(), (updating ? "Y" : "N"), (buffers.at(id)->is_read_only() ? "Y" : "N"));

  std::vector<BufferPtrType> finals = buffers.flushBufferForRank(id);
  if (bufferlt == bliss::concurrent::LockType::NONE && finals.size() != nthreads) ERRORF("FAIL: expected %d threads have %lu actual.", nthreads, finals.size());
  for (auto final : finals) {
    gbytes += (final == nullptr) ? 0 : final->getSize();
    buffers.releaseBuffer(std::move(final));
  }
  finals.clear();

  if ((bytes3 + gbytes) != success3 * data.length() + bytes4) {
    ERRORF("FAIL: total bytes %d (%d + %d).  expected %ld bytes (%ld + %d).",
	 (bytes3 + gbytes), bytes3, gbytes, success3 * data.length() + bytes4, success3 * data.length(), bytes4);
  }

  //if (count7 != count/data.length()) ERRORF("FAIL: append count = %d, actual data inserted is %ld", count7, count/data.length() );
  finals = buffers.flushBufferForRank(id);
  int gbytes2 = 0;
  for (auto final : finals) {
    gbytes2 += (final == nullptr) ? 0 : final->getSize();
    buffers.releaseBuffer(std::move(final));
  }
  finals.clear();

  if (gbytes2 != 0) {
    ERRORF("FAIL: number of bytes STILL in message buffers is %d",gbytes2);
  }


//  expectedFull = success3 / (bufferSize/data.length()) - (success3 % (bufferSize/data.length()) == 0 ? 1 : 0);
//  expectedFull2 = success3 / (bufferSize/data.length());
//
////  if (count1 != nelems) ERRORF("FAIL: number of successful inserts should be %d.  actual %d", nelems, count1);
////  else if (count2 != 0) ERRORF("FAIL: number of failed insert overall should be 0. actual %d", count2);
////  else
//  if ((sswap3 + fswap3) != expectedFull && (fswap3 + sswap3) != expectedFull2) {
//    ERRORF("FAIL: number of full Buffers from failed insert is not right: %d+%d should be %d or %d", fswap3, sswap3, expectedFull, expectedFull2);
//    fflush(stdout);
//  }
//  else
    INFOF("PASS");


  INFOF("Number of failed attempt to append to buffer is %d, success %d. full buffers size: %lu.  success swapped = %d, fail swapped = %d. total bytes %d + %d", failure3, success3, fullBuffers.getSize(), sswap3, fswap3, bytes3, gbytes);


};



template<typename BuffersType>
void testBuffersWaitForInsert(BuffersType && buffers, bliss::concurrent::LockType poollt, bliss::concurrent::LockType bufferlt, int nthreads) {

  INFOF("*** TESTING Buffers WaitForInsert  lock %d buffer lock %d: ntargets = %lu, pool threads %d", poollt, bufferlt, buffers.getSize(), nthreads);

  typedef typename BuffersType::BufferPtrType BufferPtrType;
  bool op_suc = false;
  BufferPtrType ptr = nullptr;
  bliss::concurrent::ThreadSafeQueue<typename BuffersType::BufferPtrType, bliss::concurrent::LockType::MUTEX> fullBuffers;

  //INFOF("test string is \"%s\", length %lu", data.c_str(), data.length());
  int id = 0;

  int i;
  int swap = 0;

  buffers.reset();

  INFOF("TEST all operations together: ");
  int bytes= 0;
  int gbytes = 0, gbytes2 = 0;
  //INFOF("full buffer: ");
  id = 0;
  int attempts = 0;
#pragma omp parallel for num_threads(nthreads) default(none) private(i, op_suc, ptr) firstprivate(id) shared(buffers, data, nelems, bufferSize) reduction(+ : swap, bytes, attempts)
  for (i = 0; i < nelems; ++i) {

    do {
      uint32_t count = data.length();
      char* data_remain = nullptr;
      uint32_t count_remain = 0;

      std::tie(op_suc, ptr) = buffers.append(data.c_str(), count, data_remain, count_remain, id);
      ++attempts;

      if (ptr) {
        ++swap;
//        count7 = count1;  // save the number of successful inserts so far.
        bool updating = ptr->is_writing();
        if (updating) INFOF("  FULLBUFFER1: size %ld updating? %s, blocked? %s", ptr->getSize(), (updating ? "Y" : "N"), (ptr->is_read_only() ? "Y" : "N"));

        bytes += ptr->getSize();
        buffers.releaseBuffer(std::move(ptr));
      }


    } while (!op_suc);

  }

  ptr = buffers.at(id);
  ptr->block_and_flush();
  bool updating = ptr->is_writing();
  if (updating) INFOF("  PreFLUSH: size %ld updating? %s, blocked? %s", buffers.at(id)->getSize(), (updating ? "Y" : "N"), (buffers.at(id)->is_read_only() ? "Y" : "N"));

  std::vector<BufferPtrType> finals = buffers.flushBufferForRank(id);
  if (bufferlt == bliss::concurrent::LockType::NONE && finals.size() != nthreads) ERRORF("FAIL: expected %d threads have %lu actual.", nthreads, finals.size());
  for (auto final : finals) {
    gbytes += (final == nullptr) ? 0 : final->getSize();
    buffers.releaseBuffer(std::move(final));
  }
  finals.clear();

  if ((bytes + gbytes) != nelems * data.length()) {
    ERRORF("FAIL: total bytes %d (%d + %d) for %ld entries.  expected %d entries.", (bytes + gbytes), bytes, gbytes, (bytes + gbytes)/data.length(), nelems);
  }

  //if (count7 != count/data.length()) ERRORF("FAIL: append count = %d, actual data inserted is %ld", count7, count/data.length() );
  finals = buffers.flushBufferForRank(id);
  for (auto final : finals) {
    gbytes2 += (final == nullptr) ? 0 : final->getSize();
    buffers.releaseBuffer(std::move(final));
  }
  finals.clear();

  if (gbytes2 != 0) {
    ERRORF("FAIL: number of bytes STILL in message buffers is %d",gbytes2);
  }

  INFOF("PASS");

  INFOF("Number appended to buffer is %d, total attempts is %d. full buffers size: %lu.  success swapped = %d. total bytes %d + %d = %d", nelems, attempts, fullBuffers.getSize(), swap, bytes, gbytes, (bytes + gbytes));


};

int main(int argc, char** argv) {

  // construct, acquire, access, release
  if (argc > 1) {
    nelems = atoi(argv[1]);
  }


  constexpr bliss::concurrent::LockType lt = bliss::concurrent::LockType::THREADLOCAL;
  constexpr bliss::concurrent::LockType lt2 = bliss::concurrent::LockType::NONE;

#if defined( BLISS_MUTEX )
  constexpr bliss::concurrent::LockType lt1 = bliss::concurrent::LockType::MUTEX;

//#elif defined( BLISS_THREADLOCAL_SPINLOCK_NONE )
//  constexpr bliss::concurrent::LockType lt = bliss::concurrent::LockType::THREADLOCAL;
//  constexpr bliss::concurrent::LockType lt1 = bliss::concurrent::LockType::SPINLOCK;
//  constexpr bliss::concurrent::LockType lt2 = bliss::concurrent::LockType::NONE;


/// DISABLED BECAUSE SWAPPING BUFFER PTRS IN MUTLITHREADED ENVIRONMENT IS NOT SAFE, because threads hold on to ptrs to perform tasks.
//#elif defined( BLISS_MUTEX_LOCKFREE )
//  constexpr bliss::concurrent::LockType lt = bliss::concurrent::LockType::MUTEX;
//  constexpr bliss::concurrent::LockType lt2 = bliss::concurrent::LockType::LOCKFREE;
//#elif defined( BLISS_SPINLOCK_LOCKFREE )
//  constexpr bliss::concurrent::LockType lt = bliss::concurrent::LockType::SPINLOCK;
//  constexpr bliss::concurrent::LockType lt2 = bliss::concurrent::LockType::LOCKFREE;
#endif



  /// thread unsafe.  test in single thread way.
//while(true) {

  bliss::concurrent::ObjectPool<bliss::io::Buffer<lt2, 2047, 0>, lt1 > pool;

  for (int i = 1; i <= 8; ++i) {  // num targets
    //testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::LockType::NONE, bliss::concurrent::LockType::NONE, 2047>(i,1)), bliss::concurrent::LockType::NONE, bliss::concurrent::LockType::NONE, 1);

    for (int j = 1; j <= 8; ++j) {  // num threads
      testBuffers(std::move(bliss::io::SendMessageBuffers<lt, bliss::concurrent::ObjectPool<bliss::io::Buffer<lt2, 2047, 0>, lt1  > >(pool, i, j)), lt, lt2, j);

      //testBuffersWaitForInsert(std::move(bliss::io::SendMessageBuffers<lt, bliss::concurrent::ObjectPool<bliss::io::Buffer<lt2, 2047, 0> , lt1> >(pool, i, j)), lt, lt2, j);

    }
  }

  // this one is not defined because it's not logical.  not compilable.
  // bliss::io::BufferPool<bliss::concurrent::THREAD_UNSAFE, bliss::concurrent::THREAD_SAFE> tsusPool(8192, 8);


//}



}
