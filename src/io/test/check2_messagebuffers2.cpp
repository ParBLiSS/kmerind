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

// TODO: replace content here to use messagebuffers.


#include <unistd.h>  // for usleep

#include "utils/logging.h"

#include "io/message_buffers.hpp"
#include "concurrent/referenced_object_pool.hpp"
#include "concurrent/lockfree_queue.hpp"
#include "omp.h"
#include <cassert>
#include <chrono>
#include <cstdlib>   // for rand
#include <algorithm>

#include "io/message_types.hpp"

int nelems = 11;
int bufferSize = 2047;



template<typename BuffersType>
void testBuffers(BuffersType && buffers, bliss::concurrent::LockType poollt, bliss::concurrent::LockType bufferlt, int nthreads) {

  INFOF("*** TESTING Buffers lock %d buffer lock %d: ntargets = %lu, pool threads %d", poollt, bufferlt, buffers.getSize(), nthreads);


  INFOF("TEST append until full: ");
  typedef typename BuffersType::BufferPtrType BufferPtrType;
  bool op_suc = false;
  BufferPtrType ptr = nullptr;
  bliss::concurrent::ThreadSafeQueue<typename BuffersType::BufferPtrType, bliss::concurrent::LockType::LOCKFREE> fullBuffers;

  //INFOF("test string is \"%s\", length %lu", data.c_str(), sizeof(int));
  int id = 0;

  int i;
  int success = 0;
  int failure = 0;
  int sswap = 0;
  int fswap = 0;


#pragma omp parallel for num_threads(nthreads) default(none) private(i, op_suc, ptr) firstprivate(id) shared(buffers, fullBuffers, nelems, bufferSize) reduction(+ : success, failure, sswap, fswap)
  for (i = 0; i < nelems; ++i) {
    //INFOF("insert %lu chars into %d", data.length(), id);
    //INFOF("insert %lu chars into %d", sizeof(int), id);
    int data = i;
    std::tie(op_suc, ptr) = buffers.append(&data, sizeof(int), id);

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
  int expectedFullMin = success / (bufferSize/sizeof(int)) - nthreads;
  int expectedFullMax = success / (bufferSize/sizeof(int));

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
#pragma omp parallel for num_threads(nthreads) default(none) private(id, op_suc, ptr) shared(buffers, fullBuffers, bufferSize, iterations) reduction(+ : sswap2, fswap2, error, over)
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
  int success3 = 0, failure3 = 0, sswap3 = 0, fswap3 = 0, bytes3 = 0;
  int gbytes = 0;
  //INFOF("full buffer: ");

  std::vector< std::vector<int> > stored(nthreads);
  std::vector< std::vector<int> > appended(nthreads);


  id = 0;
#pragma omp parallel for num_threads(nthreads) default(none) private(i, op_suc, ptr) firstprivate(id) shared(stored, appended, buffers, nelems, bufferSize) reduction(+ : success3, failure3, sswap3, fswap3, bytes3)
  for (i = 0; i < nelems; ++i) {
    int data = i;
    std::tie(op_suc, ptr) = buffers.append(&data, sizeof(int), id);

    if (op_suc) {
      ++success3;
      appended[omp_get_thread_num()].push_back(data);


      if (ptr) {
      	++sswap3;
        bool updating = ptr->is_writing();
        if (updating) INFOF("  FULLBUFFER1: size %ld updating? %s, blocked? %s", ptr->getSize(), (updating ? "Y" : "N"), (ptr->is_read_only() ? "Y" : "N"));

        bytes3 += ptr->getSize();
        stored[omp_get_thread_num()].insert(stored[omp_get_thread_num()].end(), ptr->operator int*(), ptr->operator int*() + ptr->getSize() / sizeof(int));
        buffers.releaseBuffer(std::move(ptr));
      }
    } else {
      ++failure3;

      if (ptr) {
        ++fswap3;
       // count7 = count1;  // save the number of successful inserts so far.
        bool updating = ptr->is_writing();
        if (updating) INFOF("  FULLBUFFER: size %ld blocked? %s", ptr->getSize(), (ptr->is_read_only() ? "Y" : "N"));

        // result.second->lock_read();

        bytes3 += ptr->getSize();

        stored[omp_get_thread_num()].insert(stored[omp_get_thread_num()].end(), ptr->operator int*(), ptr->operator int*() + ptr->getSize() / sizeof(int));
        buffers.releaseBuffer(std::move(ptr));
      }
    }
  }

  // concatenate the vectors
  std::vector<int> allstored;
  std::vector<int> allappended;
  for (int k = 0; k < nthreads; ++k) {
    allstored.insert(allstored.end(), stored[k].begin(), stored[k].end());
    allappended.insert(allappended.end(), appended[k].begin(), appended[k].end());
  }

  buffers.at(id)->block_and_flush();
  bool updating = buffers.at(id)->is_writing();
  if (updating) INFOF("  PreFLUSH: size %ld updating? %s, blocked? %s", buffers.at(id)->getSize(), (updating ? "Y" : "N"), (buffers.at(id)->is_read_only() ? "Y" : "N"));

  std::vector<BufferPtrType> finals = buffers.flushBufferForRank(id);
  if (bufferlt == bliss::concurrent::LockType::NONE && finals.size() != nthreads) ERRORF("FAIL: expected %d threads have %lu actual.", nthreads, finals.size());
  for (auto final : finals) {
    gbytes += (final == nullptr) ? 0 : final->getSize();
    if (final) {
      allstored.insert(allstored.end(), final->operator int*(),  final->operator int*() + final->getSize() / sizeof(int));
      buffers.releaseBuffer(std::move(final));
    }
    //if (count7 != count/sizeof(int)) ERRORF("FAIL: append count = %d, actual data inserted is %ld", count7, count/sizeof(int) );
  }
  finals.clear();

  if ((bytes3 + gbytes) != success3 * sizeof(int)) {
    ERRORF("FAIL: total bytes (%d + %d) bytes.  expected %ld bytes.", bytes3, gbytes, success3 * sizeof(int));
  }

  //if (count7 != count/data.length()) ERRORF("FAIL: append count = %d, actual data inserted is %ld", count7, count/data.length() );
  finals = buffers.flushBufferForRank(id);
  int gbytes2 = 0;
  for (auto final : finals) {
    gbytes2 += (final == nullptr) ? 0 : final->getSize();
      allstored.insert(allstored.end(), final->operator int*(),  final->operator int*() + final->getSize() / sizeof(int));
    buffers.releaseBuffer(std::move(final));
  }
  finals.clear();

  if (gbytes2 != 0) {
    ERRORF("FAIL: number of bytes STILL in message buffers is %d",gbytes2);
  }


  // now sort it an check to see if we are missing anything
  std::sort(allstored.begin(), allstored.end());
  std::sort(allappended.begin(), allappended.end());
  std::vector<int> appendedButNotStored; std::vector<int> storedButNotAppended;
  std::set_difference(allappended.begin(), allappended.end(), allstored.begin(), allstored.end(), std::inserter(appendedButNotStored, appendedButNotStored.begin()));
  std::set_difference(allstored.begin(), allstored.end(), allappended.begin(), allappended.end(), std::inserter(storedButNotAppended, storedButNotAppended.begin()));



  if (appendedButNotStored.size() > 0) {
    ERRORF("FAIL!: %ld elements appended but not stored: ", appendedButNotStored.size());
    for (i = 0; i < appendedButNotStored.size(); ++i) {
      INFOF("%d, ", appendedButNotStored[i]);
    }
  }
  if (storedButNotAppended.size() > 0) {
    ERRORF("FAIL?: %ld elements stored but not appended: ", storedButNotAppended.size());
    for (i = 0; i < storedButNotAppended.size(); ++i) {
      INFOF("%d, ", storedButNotAppended[i]);
    }
  }


  if ((bytes3 + gbytes) != success3 * sizeof(int)) {
    ERRORF("FAIL: total bytes (%d + %d) for expected %ld bytes.", bytes3 , gbytes, success*sizeof(int));

    throw std::logic_error("missing entries.");
  }




//  expectedFull = (count1 - 1 + (bufferSize/sizeof(int))) / (bufferSize/sizeof(int)) - (count1 % (bufferSize/sizeof(int)) == 0 ? 1 : 0);
//  expectedFull2 = count1 / (bufferSize/sizeof(int));
//
////  if (count1 != nelems) ERRORF("FAIL: number of successful inserts should be %d.  actual %d", nelems, count1);
////  else if (count2 != 0) ERRORF("FAIL: number of failed insert overall should be 0. actual %d", count2);
////  else
//  if (count3 != 0) ERRORF("FAIL: number of full Buffers from successful insert is not right: %d should be 0", count3);
//  else if (count4 != expectedFull && count4 != expectedFull2) {
//    ERRORF("FAIL: number of full Buffers from failed insert is not right: %d should be %d or %d", count4, expectedFull, expectedFull2);
//    fflush(stdout);
//  }
//  else INFOF("PASS");

  INFOF("Number of failed attempt to append to buffer is %d, success %d. full buffers size: %lu.  success swapped = %d, fail swapped = %d. total bytes %d + %d", failure3, success3, fullBuffers.getSize(), sswap3, fswap3, bytes3, gbytes);


};



template<typename BuffersType>
void testBuffersWaitForInsert(BuffersType && buffers, bliss::concurrent::LockType poollt, bliss::concurrent::LockType bufferlt, int nthreads) {

  INFOF("*** TESTING Buffers WaitForInsert lock %d buffer lock %d: ntargets = %lu, pool threads %d", poollt, bufferlt, buffers.getSize(), nthreads);

  typedef typename BuffersType::BufferPtrType BufferPtrType;
  bool op_suc = false;
  BufferPtrType ptr = nullptr;
  bliss::concurrent::ThreadSafeQueue<typename BuffersType::BufferPtrType, bliss::concurrent::LockType::LOCKFREE> fullBuffers;

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
#pragma omp parallel for num_threads(nthreads) default(none) private(i, op_suc, ptr) firstprivate(id) shared(buffers, nelems, bufferSize) reduction(+ : swap, bytes, attempts)
  for (i = 0; i < nelems; ++i) {

    do {
      int data = i;
      std::tie(op_suc, ptr) = buffers.append(&data, sizeof(int), id);
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

  if ((bytes + gbytes) != nelems * sizeof(int)) {
    ERRORF("FAIL: total bytes %d (%d + %d) for %ld entries.  expected %d entries.", (bytes + gbytes), bytes, gbytes, (bytes + gbytes)/sizeof(int), nelems);
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


#if defined( BLISS_THREADLOCAL_MUTEX_NONE )
  constexpr bliss::concurrent::LockType lt = bliss::concurrent::LockType::THREADLOCAL;
  constexpr bliss::concurrent::LockType lt1 = bliss::concurrent::LockType::MUTEX;
  constexpr bliss::concurrent::LockType lt2 = bliss::concurrent::LockType::NONE;
#elif defined( BLISS_THREADLOCAL_SPINLOCK_NONE )
  constexpr bliss::concurrent::LockType lt = bliss::concurrent::LockType::THREADLOCAL;
  constexpr bliss::concurrent::LockType lt1 = bliss::concurrent::LockType::SPINLOCK;
  constexpr bliss::concurrent::LockType lt2 = bliss::concurrent::LockType::NONE;

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
      testBuffers(std::move(bliss::io::SendMessageBuffers<lt, bliss::concurrent::ObjectPool<bliss::io::Buffer<lt2, 2047, 0>, lt1 >>(pool, i, j)), lt, lt2, j);

      testBuffersWaitForInsert(std::move(bliss::io::SendMessageBuffers<lt, bliss::concurrent::ObjectPool<bliss::io::Buffer<lt2, 2047, 0>, lt1 >>(pool, i, j)), lt, lt2, j);

    }
  }


  // this one is not defined because it's not logical.  not compilable.
  // bliss::io::BufferPool<bliss::concurrent::THREAD_UNSAFE, bliss::concurrent::THREAD_SAFE> tsusPool(8192, 8);


//}



}
