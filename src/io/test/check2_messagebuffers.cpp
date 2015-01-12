/**
 * @file		check_messagebuffers.cpp
 * @ingroup
 * @author	tpan
 * @brief
 * @details
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */

// TODO: replace content here to use messagebuffers.


#include <unistd.h>  // for usleep

#include "io/message_buffers.hpp"
#include "concurrent/lockfree_queue.hpp"
#include "omp.h"
#include <cassert>
#include <chrono>
#include <cstdlib>   // for rand

int nelems = 11;
int bufferSize = 2048;
std::string data("this is a test.  this a test of the emergency broadcast system.  this is only a test. ");


template<typename BuffersType>
void testPool(BuffersType && buffers, bliss::concurrent::LockType poollt, bliss::concurrent::LockType bufferlt, int nthreads) {

  printf("*** TESTING pool lock %d buffer lock %d: ntargets = %lu, pool threads %d\n", poollt, bufferlt, buffers.getSize(), nthreads);


  printf("TEST append until full: ");
  typedef typename BuffersType::BufferPtrType BufferPtrType;
  bool op_suc = false;
  BufferPtrType ptr = nullptr;
  bliss::concurrent::ThreadSafeQueue<typename BuffersType::BufferPtrType> fullBuffers;

  //printf("test string is \"%s\", length %lu\n", data.c_str(), data.length());
  int id = 0;

  int i;
  int success = 0;
  int failure = 0;
  int sswap = 0;
  int fswap = 0;


#pragma omp parallel for num_threads(nthreads) default(none) private(i, op_suc, ptr) firstprivate(id) shared(buffers, data, fullBuffers, nelems, bufferSize) reduction(+ : success, failure, sswap, fswap)
  for (i = 0; i < nelems; ++i) {
    //printf("insert %lu chars into %d\n", data.length(), id);

    std::tie(op_suc, ptr) = buffers.append(data.c_str(), data.length(), id);

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
//  if ((count + count2) != nelems) printf("\nFAIL: number of successful inserts should be %d.  actual %d", nelems, count);
//  else if (count2 != 0) printf("\nFAIL: number of failed insert overall should be 0. actual %d", count2);
//  else
//    if (!(count5 <= count && count <= (count5 + count3))) printf("\nFAIL: number of successful inserts should be close to successful inserts without full buffers.");

  // compute 2 expected since append returns full buffer only on failed insert and if there are n inserts that brings it to just before full, then it depends on timing
  // as to when the buffer becomes full.  (other threads will fail on append until swap happens, but not return a full buffer.)
  int expectedFull = success / (bufferSize/data.length()) - (success % (bufferSize/data.length()) == 0 ? 1 : 0);
  int expectedFull2 = success / (bufferSize/data.length());

  if (fullBuffers.getSize() != sswap + fswap) printf("\nFAIL: number of full Buffers do not match: fullbuffer size %ld  full count %d + %d", fullBuffers.getSize(), sswap, fswap);
  // buffer at 23 entries (86 bytes each, 2048 bytes per buffer) will not show as full until the next iterator.
  else if (fswap + sswap != expectedFull && sswap + fswap != expectedFull2) printf("\nFAIL: number of full Buffers is not right: %d+%d should be %d or %d", sswap, fswap, expectedFull, expectedFull2);
  //else if (count4 != 0) printf("\nFAIL: number of failed insert due to no buffer should be 0. actual %d", count4);
  else printf("PASS");
  printf("\n");
  printf("Number of failed attempt to append to buffer is %d, success %d. full buffers size: %lu.  success swapped = %d, fail swapped = %d.\n", failure, success, fullBuffers.getSize(), sswap, fswap);



  printf("TEST release: ");
  int sswap2 = 0, fswap2 = 0, error = 0, over = 0;
//  printf("buffer ids 0 to %lu initially in use.\n", buffers.getSize() - 1);
//  printf("releasing: ");

  int iterations = (sswap + fswap) + 10;
#pragma omp parallel for num_threads(nthreads) default(none) private(id, op_suc, ptr) shared(buffers, data, fullBuffers, bufferSize, iterations) reduction(+ : sswap2, fswap2, error, over)
  for (id = 0; id < iterations; ++id) {
    try {
      std::tie(op_suc, ptr) = fullBuffers.tryPop();
      if (op_suc) {
//        printf("%d ", result.second);
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
      printf("\nFAIL with %s", e.what());
      ++error;       // error during pop
    }
  }

  expectedFull = success / (bufferSize/data.length()) - (success % (bufferSize/data.length()) == 0 ? 1 : 0);
  expectedFull2 = success / (bufferSize/data.length());

  if (error != 0) printf("\nFAIL: invalid argument exception during pop.  count = %d", error);
  else if (over != 10) printf("\nFAIL: failed on pop %d times, expected 10", over);
  else if (fswap2 != 0) printf("\nFAIL: succeeded in pop but not full buffer. %d", fswap2);
  else if (sswap2 != expectedFull && sswap2 != expectedFull2) printf("FAIL: expected %d or %d full buffers, but received %d", expectedFull, expectedFull2, sswap2);
  else if (sswap2 != (sswap + fswap)) printf("\nFAIL: successful pops. expected %d.  actual %d", (sswap+fswap), sswap2);
  else
    printf("PASS");
  printf("\n");


  buffers.reset();

  printf("TEST all operations together: ");
  int success3 = 0, failure3 = 0, sswap3 = 0, fswap3 = 0, bytes3 = 0, chars3 = 0;
  int gbytes = 0, gchars = 0;
  //printf("full buffer: ");
  id = 0;
#pragma omp parallel for num_threads(nthreads) default(none) private(i, op_suc, ptr) firstprivate(id) shared(buffers, data, nelems, bufferSize) reduction(+ : success3, failure3, sswap3, fswap3, bytes3, chars3)
  for (i = 0; i < nelems; ++i) {
    std::tie(op_suc, ptr) = buffers.append(data.c_str(), data.length(), id);

    if (op_suc) {
      ++success3;

      if (ptr) {
        ++sswap3;
//        count7 = count1;  // save the number of successful inserts so far.
        bool updating = ptr->is_writing();
        if (updating) printf("  FULLBUFFER1: size %ld updating? %s, blocked? %s\n", ptr->getSize(), (updating ? "Y" : "N"), (ptr->is_read_only() ? "Y" : "N"));

        bytes3 += ptr->getSize();
        chars3 += strlen(ptr->operator char*());
        buffers.releaseBuffer(std::move(ptr));
      }
    } else {
      ++failure3;

      if (ptr) {
        ++fswap3;
//        count7 = count1;  // save the number of successful inserts so far.
        bool updating = ptr->is_writing();
        if (updating) printf("  FULLBUFFER: size %ld blocked? %s\n", ptr->getSize(), (ptr->is_read_only() ? "Y" : "N"));

        bytes3 += ptr->getSize();
        chars3 += strlen(ptr->operator char*());

        buffers.releaseBuffer(std::move(ptr));
      }
    }
  }

  buffers.at(id)->block_and_flush();
  bool updating = buffers.at(id)->is_writing();
  if (updating) printf("  PreFLUSH: size %ld updating? %s, blocked? %s\n", buffers.at(id)->getSize(), (updating ? "Y" : "N"), (buffers.at(id)->is_read_only() ? "Y" : "N"));

  std::vector<BufferPtrType> finals = buffers.flushBufferForRank(id);
  if (finals.size() != nthreads) printf("\nFAIL: expected %d threads have %lu actual.\n", nthreads, finals.size());
  for (auto final : finals) {
    gbytes += (final == nullptr) ? 0 : final->getSize();
    buffers.releaseBuffer(std::move(final));
  }
  finals.clear();

  if ((bytes3 + gbytes) != success3 * data.length()) {
    printf("\nFAIL: total bytes %d (%d + %d) for %ld entries.  expected %d entries.\n", (bytes3 + gbytes), bytes3, gbytes, (bytes3 + gbytes)/data.length(), success3);
  }

  //if (count7 != count/data.length()) printf("\nFAIL: append count = %d, actual data inserted is %ld", count7, count/data.length() );
  finals = buffers.flushBufferForRank(id);
  int gbytes2 = 0;
  for (auto final : finals) {
    gbytes2 += (final == nullptr) ? 0 : final->getSize();
    buffers.releaseBuffer(std::move(final));
  }
  finals.clear();

  if (gbytes2 != 0) {
    printf("\nFAIL: number of bytes STILL in message buffers is %d\n",gbytes2);
  }


//  expectedFull = success3 / (bufferSize/data.length()) - (success3 % (bufferSize/data.length()) == 0 ? 1 : 0);
//  expectedFull2 = success3 / (bufferSize/data.length());
//
////  if (count1 != nelems) printf("\nFAIL: number of successful inserts should be %d.  actual %d", nelems, count1);
////  else if (count2 != 0) printf("\nFAIL: number of failed insert overall should be 0. actual %d", count2);
////  else
//  if ((sswap3 + fswap3) != expectedFull && (fswap3 + sswap3) != expectedFull2) {
//    printf("\nFAIL: number of full Buffers from failed insert is not right: %d+%d should be %d or %d", fswap3, sswap3, expectedFull, expectedFull2);
//    fflush(stdout);
//  }
//  else
    printf("PASS");
  printf("\n");

  //printf("\n");
  printf("Number of failed attempt to append to buffer is %d, success %d. full buffers size: %lu.  success swapped = %d, fail swapped = %d. total bytes %d + %d\n", failure3, success3, fullBuffers.getSize(), sswap3, fswap3, bytes3, gbytes);


};


int main(int argc, char** argv) {

  // construct, acquire, access, release
  if (argc > 1) {
    nelems = atoi(argv[1]);
  }


#if defined( BLISS_MUTEX_LOCKFREE )
  constexpr bliss::concurrent::LockType lt = bliss::concurrent::LockType::MUTEX;
  constexpr bliss::concurrent::LockType lt2 = bliss::concurrent::LockType::LOCKFREE;
#elif defined( BLISS_SPINLOCK_LOCKFREE )
  constexpr bliss::concurrent::LockType lt = bliss::concurrent::LockType::SPINLOCK;
  constexpr bliss::concurrent::LockType lt2 = bliss::concurrent::LockType::LOCKFREE;
#elif defined( BLISS_MUTEX_NONE )
  constexpr bliss::concurrent::LockType lt = bliss::concurrent::LockType::MUTEX;
  constexpr bliss::concurrent::LockType lt2 = bliss::concurrent::LockType::NONE;
#elif defined( BLISS_SPINLOCK_NONE )
  constexpr bliss::concurrent::LockType lt = bliss::concurrent::LockType::SPINLOCK;
  constexpr bliss::concurrent::LockType lt2 = bliss::concurrent::LockType::NONE;
#endif



  /// thread unsafe.  test in single thread way.
//while(true) {

  for (int i = 1; i <= 8; ++i) {  // num targets
    //testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::LockType::NONE, bliss::concurrent::LockType::NONE, 2047>(i)), bliss::concurrent::LockType::NONE, bliss::concurrent::LockType::NONE, 1);

    for (int j = 1; j <= 8; ++j) {  // num threads
      testPool(std::move(bliss::io::SendMessageBuffers<lt, lt2, 2047>(i)), lt, lt2, j);
    }
  }

  // this one is not defined because it's not logical.  not compilable.
  // bliss::io::BufferPool<bliss::concurrent::THREAD_UNSAFE, bliss::concurrent::THREAD_SAFE> tsusPool(8192, 8);


//}



}
