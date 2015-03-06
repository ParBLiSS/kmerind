/**
 * @file		check_bufferpool.cpp
 * @ingroup
 * @author	Tony Pan <tpan7@gatech.edu>
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

#include "utils/iterator_test_utils.hpp"

#include "concurrent/buffer.hpp"
#include "utils/logging.h"

#include "concurrent/unreferenced_object_pool.hpp"


template<typename PoolType>
void testAppendMultipleBuffers(const int NumThreads, const int total_count, bliss::concurrent::LockType poollt, bliss::concurrent::LockType bufferlt, const int64_t buffer_cap) {
  omp_lock_t writelock;
  omp_init_lock(&writelock);
  omp_lock_t writelock2;
  omp_init_lock(&writelock2);
  omp_lock_t writelock3;
  omp_init_lock(&writelock3);
//  std::atomic_flag writelock = ATOMIC_FLAG_INIT;
//  std::atomic_flag writelock2 = ATOMIC_FLAG_INIT;
//  std::atomic_flag writelock3 = ATOMIC_FLAG_INIT;


  INFOF("TESTING: %d threads, pool lock %d buffer lock %d append with %ld bufferSize and %d total counts from unlimited pool\n",
         NumThreads, poollt, bufferlt, buffer_cap, total_count);


  PoolType pool;

  std::vector<int> gold;
  std::vector<int> stored;

  int data = 0;
  unsigned int result = 0;

  int success = 0;
  int failure = 0;
  int swap = 0;
  int i = 0;



  auto buf_ptr = pool.tryAcquireObject();
  if (buf_ptr) buf_ptr->clear_and_unblock_writes();


#pragma omp parallel for num_threads(NumThreads) default(none) shared(buf_ptr, gold, stored, writelock, writelock2, writelock3, pool) private(i, data, result) reduction(+:success, failure, swap)
  for (i = 0; i < total_count; ++i) {

//    while (writelock2.test_and_set());
//    omp_set_lock(&writelock2);
    auto sptr = buf_ptr;
//    omp_unset_lock(&writelock2);
//    writelock2.clear();

    if (sptr) {  // valid ptr

      data = static_cast<int>(i);
      result = sptr->append(&data, sizeof(int));
    } else {  // expired ptr
      result = 0x0;
    }

    if (result & 0x1) {
      ++success;
//      while (writelock.test_and_set());
      omp_set_lock(&writelock);
      gold.push_back(data);
      omp_unset_lock(&writelock);
//      writelock.clear();
    }
    else ++failure;

    if (result & 0x2) {
      ++swap;

      // swap in a new one.
      auto new_buf_ptr = pool.tryAcquireObject();
      if (new_buf_ptr) new_buf_ptr->clear_and_unblock_writes();

//      while (writelock2.test_and_set());
//      omp_set_lock(&writelock2);
      sptr = buf_ptr;
      buf_ptr = new_buf_ptr;
#pragma omp flush(buf_ptr)
      //new_buf_ptr = tmp;
//      omp_unset_lock(&writelock2);
//      writelock2.clear();

      // process the old buffer
//      while (writelock3.test_and_set());
      //sptr = new_buf_ptr;
      omp_set_lock(&writelock3);
      if (sptr)
      stored.insert(stored.end(), sptr->operator int*(), sptr->operator int*() + sptr->getSize() / sizeof(int));
      omp_unset_lock(&writelock3);
//      writelock3.clear();

      // and release - few threads doing this, and full.

      pool.releaseObject(sptr);

    }

  }

  auto sptr = buf_ptr;
  if (sptr) {sptr->block_and_flush();

    // compare unordered buffer content.
    stored.insert(stored.end(), sptr->operator int*(), sptr->operator int*() + sptr->getSize() / sizeof(int));
  }
  pool.releaseObject(buf_ptr);
  int stored_count = stored.size();


  if ( swap != success / (buffer_cap / sizeof(int)) || success != stored_count)
      FATALF("FAIL: (actual/expected)  success (%d/%d), failure (%d/?), swap(%ld/%d).\n", success, stored_count, failure, success / (buffer_cap / sizeof(int)), swap);
  else {
    INFOF("INFO: success %d, failure %d, swap %d, total %d\n", success, failure, swap, total_count);

    if (compareUnorderedSequences(stored.begin(), gold.begin(), stored_count)) {
      INFOF("PASS\n");
    } else {
      FATALF("FAIL: content not matching\n");
    }
  }
  omp_destroy_lock(&writelock);
  omp_destroy_lock(&writelock2);
  omp_destroy_lock(&writelock3);
}


template<typename PoolType>
void testPool(PoolType && pool, bliss::concurrent::LockType poollt, bliss::concurrent::LockType bufferlt, int pool_threads, int buffer_threads) {


  INFOF("TESTING pool lock %d buffer lock %d %s: pool threads %d, buffer threads %d\n", poollt, bufferlt, (pool.isUnlimited() ? "GROW" : "FIXED"),  pool_threads, buffer_threads);

  INFOF("TEST acquire: ");
  int expected;
  int i = 0;
  int count = 0;
  int mx = pool.isUnlimited() ? 100 : pool.getCapacity();
#pragma omp parallel for num_threads(pool_threads) default(none) private(i) shared(pool, mx) reduction(+ : count)
  for (i = 0; i < mx; ++i) {
	  auto ptr = pool.tryAcquireObject();
    if (! ptr) {
      ++count;
      delete ptr;  // clean up since we are not tracking what was acquired.
    }
  }
  expected = 0;
  if (count != expected) FATALF("FAIL: number of failed attempt to acquire buffer should be %d, actual %d.  pool capacity %lu, remaining: %lu \n", expected, count, pool.getCapacity(), pool.getAvailableCount());
  else INFOF("PASSED.\n");
  pool.reset();

  INFOF("TEST acquire with growth: ");
  i = 0;
  count = 0;
  mx = pool.isUnlimited() ? 100 : pool.getCapacity();
#pragma omp parallel for num_threads(pool_threads) default(none) private(i) shared(pool, mx) reduction(+ : count)
  for (i = 0; i <= mx; ++i) {  // <= so we get 1 extra
		auto ptr = pool.tryAcquireObject();
	    if (! ptr) {
	      ++count;
	      delete ptr;  // clean up since we are not tracking what was acquired.
	    }
  }
  expected = pool.isUnlimited() ? 0 : 1;
  if (count != expected) FATALF("FAIL: number of failed attempt to acquire buffer should be %d, actual %d.  pool remaining: %lu \n", expected, count, pool.getAvailableCount());
  else INFOF("PASSED.\n");
  pool.reset();

  INFOF("TEST release: ");
  count = 0;
  mx = pool.isUnlimited() ? 100 : pool.getCapacity();
  // and create some dummy buffers to insert
  std::vector<typename PoolType::ObjectPtrType> temp;
  // first drain the pool
  for (i = 0; i < mx; ++i) {
    typename PoolType::ObjectPtrType ptr = pool.tryAcquireObject();
    ptr->block_and_flush();
    temp.push_back(ptr);
  }
  pool.reset();  // reset so can get another 100 or up to capacity.
  for (i = 0; i < mx; ++i) {
    typename PoolType::ObjectPtrType ptr = pool.tryAcquireObject();
    ptr->block_and_flush();
    temp.push_back(ptr);
  }
#pragma omp parallel for num_threads(pool_threads) default(none) shared(pool, mx, temp) private(i) reduction(+ : count)
  for (i = 0; i < mx * 2; ++i) {

    typename PoolType::ObjectPtrType ptr = temp[i];
    if (ptr) {
      if (! pool.releaseObject(ptr)) {
        ++count; // failed release
        delete ptr;
      }
    }
  }
  expected = mx;  // unlimited or not, can only push back in as much as taken out.
  if (count != expected) FATALF("FAIL: number of failed attempt to release buffer should be %d, actual %d. pool remaining: %lu \n", expected, count, pool.getAvailableCount());
  else INFOF("PASSED.\n");
  pool.reset();
  temp.clear();


  // NOTE THAT WE CAN'T RELEASE THE SAME OBJECT more than 1x

  INFOF("TEST access by multiple threads, each a separate buffer: ");

  count = 0;
  int count1 = 0;
  int count2 = 0;
#pragma omp parallel num_threads(pool_threads) default(none) shared(pool, std::cout) reduction(+ : count, count1, count2)
  {
    int v = omp_get_thread_num() + 5;
    auto ptr = pool.tryAcquireObject();

    if (! ptr) ++count2;
    else {
      ptr->clear_and_unblock_writes();

      int res = ptr->append(&v, sizeof(int));


      if (! (res & 0x1)) {
        ++count;
      }

      int u = ptr->operator int*()[0];
      if (v != u) {
        ++count1;
      }

      ptr->block_and_flush();
      pool.releaseObject(ptr);
    }
  }
  if (count2 != 0) FATALF("FAIL: acquire failed\n");
  else if (count != 0) FATALF("FAIL: append failed\n");
  else if (count1 != 0) FATALF("FAIL: inserted and got back wrong values\n");
  else INFOF("PASSED.\n");
  pool.reset();

  INFOF("TEST access by multiple threads, all to same buffer: ");
  auto ptr = pool.tryAcquireObject();
  ptr->clear_and_unblock_writes();

#pragma omp parallel num_threads(buffer_threads) default(none) shared(pool, ptr)
  {
    int v = 7;
    ptr->append(&v, sizeof(int));
  }

  bool same = true;
  for (int i = 0; i < buffer_threads ; ++i) {
    same &= ptr->operator int*()[i] == 7;
  }
  if (!same) FATALF("FAIL: inserted not same\n");
  else INFOF("PASSED.\n");
  pool.releaseObject(ptr);
  pool.reset();


  omp_set_nested(1);

  INFOF("TEST all operations together: ");
#pragma omp parallel num_threads(pool_threads) default(none) shared(pool, pool_threads, buffer_threads, std::cout)
  {
    // Id range is 0 to 100
    int iter;
    int j = 0;
    for (int i = 0; i < 100; ++i) {
      // acquire
      auto buf = pool.tryAcquireObject();
//      INFOF("acquiring ");
      while (!buf) {
        _mm_pause();
//        INFOF(".");
        buf = pool.tryAcquireObject();
      }
//      INFOF(" done.\n");
      buf->clear_and_unblock_writes();

      // access
      iter = rand() % 100;
      int count = 0;
#pragma omp parallel for num_threads(buffer_threads) default(none) shared(buf, iter) private(j) reduction(+:count)
      for (j = 0; j < iter; ++j) {
        bool res = buf->append(&j, sizeof(int));
        if (! (res & 0x1)) {
          count++;
        }
      }

      // random sleep
      for (int i = 0; i < rand() % 1000; ++i) {
        _mm_pause();
      }
      // clear buffer
//      INFO( "before block and flush: " << *buf );
      buf->block_and_flush();
//      INFO( "after block and flush: " << *buf );

//      INFOF("count = %d\n", count);

      if (buf->getSize() != sizeof(int) * iter  || count != 0)
        FATALF("FAIL: thread %d/%d buffer size is %ld, expected %lu\n", omp_get_thread_num() + 1, pool_threads, buf->getSize(), sizeof(int) * iter);
 // else INFOF("PASSED.\n");

      //release
      pool.releaseObject(buf);
      //if (i % 25 == 0)
//      INFOF("thread %d released buffer %d\n", omp_get_thread_num(), id);

    }
  }




};


int main(int argc, char** argv) {

  // construct, acquire, access, release
#if defined( BLISS_NONE )
  constexpr bliss::concurrent::LockType lt = bliss::concurrent::LockType::NONE;
#else // #ifdef BLISS_LOCKFREE
  constexpr bliss::concurrent::LockType lt = bliss::concurrent::LockType::LOCKFREE;
#endif



  //////////////  unbounded version

  /// thread unsafe.  test in single thread way.


  testPool(std::move(bliss::concurrent::ObjectPool< bliss::io::Buffer<bliss::concurrent::LockType::NONE, 8192>, bliss::concurrent::LockType::NONE, false >()), bliss::concurrent::LockType::NONE,bliss::concurrent::LockType::NONE, 1, 1);
  testPool(std::move(bliss::concurrent::ObjectPool< bliss::io::Buffer<bliss::concurrent::LockType::NONE, 8192>, bliss::concurrent::LockType::NONE, false >(16)), bliss::concurrent::LockType::NONE,bliss::concurrent::LockType::NONE, 1, 1);

  for (int i = 1; i <= 8; ++i) {
	if (i == 5 || i == 6 || i == 7) continue;

    // okay to test.  in real life, pools would not be single threaded.
    testPool(std::move(bliss::concurrent::ObjectPool<bliss::io::Buffer<bliss::concurrent::LockType::MUTEX, 8192>   , bliss::concurrent::LockType::NONE, false >()),    bliss::concurrent::LockType::NONE, bliss::concurrent::LockType::MUTEX, 1, i);
    testPool(std::move(bliss::concurrent::ObjectPool<bliss::io::Buffer<bliss::concurrent::LockType::SPINLOCK, 8192>, bliss::concurrent::LockType::NONE, false >()), bliss::concurrent::LockType::NONE, bliss::concurrent::LockType::SPINLOCK, 1, i);
    testPool(std::move(bliss::concurrent::ObjectPool<bliss::io::Buffer<bliss::concurrent::LockType::LOCKFREE, 8192>, bliss::concurrent::LockType::NONE, false >()), bliss::concurrent::LockType::NONE, bliss::concurrent::LockType::LOCKFREE, 1, i);

    // okay to test.  in real life, pools would not be single threaded.
    testPool(std::move(bliss::concurrent::ObjectPool<bliss::io::Buffer<bliss::concurrent::LockType::MUTEX, 8192>   , bliss::concurrent::LockType::NONE, false >(16)),    bliss::concurrent::LockType::NONE, bliss::concurrent::LockType::MUTEX, 1, i);
    testPool(std::move(bliss::concurrent::ObjectPool<bliss::io::Buffer<bliss::concurrent::LockType::SPINLOCK, 8192>, bliss::concurrent::LockType::NONE, false >(16)), bliss::concurrent::LockType::NONE, bliss::concurrent::LockType::SPINLOCK, 1, i);
    testPool(std::move(bliss::concurrent::ObjectPool<bliss::io::Buffer<bliss::concurrent::LockType::LOCKFREE, 8192>, bliss::concurrent::LockType::NONE, false >(16)), bliss::concurrent::LockType::NONE, bliss::concurrent::LockType::LOCKFREE, 1, i);


    for (int j = 1; j <= 4; ++j) {
	if (i * j > 16) continue;
      testPool(std::move(bliss::concurrent::ObjectPool<bliss::io::Buffer<bliss::concurrent::LockType::MUTEX, 8192>   , lt, false >()),    lt, bliss::concurrent::LockType::MUTEX, j, i);
      testPool(std::move(bliss::concurrent::ObjectPool<bliss::io::Buffer<bliss::concurrent::LockType::SPINLOCK, 8192>, lt, false >()), lt, bliss::concurrent::LockType::SPINLOCK, j, i);
      testPool(std::move(bliss::concurrent::ObjectPool<bliss::io::Buffer<bliss::concurrent::LockType::LOCKFREE, 8192>, lt, false >()), lt, bliss::concurrent::LockType::LOCKFREE, j, i);

      testPool(std::move(bliss::concurrent::ObjectPool<bliss::io::Buffer<bliss::concurrent::LockType::MUTEX, 8192>   , lt, false >(16)),    lt, bliss::concurrent::LockType::MUTEX, j, i);
      testPool(std::move(bliss::concurrent::ObjectPool<bliss::io::Buffer<bliss::concurrent::LockType::SPINLOCK, 8192>, lt, false >(16)), lt, bliss::concurrent::LockType::SPINLOCK, j, i);
      testPool(std::move(bliss::concurrent::ObjectPool<bliss::io::Buffer<bliss::concurrent::LockType::LOCKFREE, 8192>, lt, false >(16)), lt, bliss::concurrent::LockType::LOCKFREE, j, i);

    }

    testAppendMultipleBuffers<bliss::concurrent::ObjectPool<bliss::io::Buffer<bliss::concurrent::LockType::LOCKFREE, 8192>, lt, false > >(i, 1000000, lt, bliss::concurrent::LockType::LOCKFREE, 8192);


    // no multithread pool single thread buffer test right now.
    //testPool(std::move(bliss::concurrent::ObjectPool<lt, bliss::io::Buffer<bliss::concurrent::LockType::NONE, 8192> >()), lt, bliss::concurrent::LockType::NONE, i, 1);
    //testPool(std::move(bliss::concurrent::ObjectPool<lt, bliss::io::Buffer<bliss::concurrent::LockType::NONE, 8192> >(16)), lt, bliss::concurrent::LockType::NONE, i, 1);
  }


}
