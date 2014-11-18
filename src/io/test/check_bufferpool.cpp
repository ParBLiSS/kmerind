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

#include "io/buffer_pool.hpp"
#include "omp.h"
#include <cassert>
#include <chrono>
#include <vector>
#include <cstdlib>   // for rand
#include <atomic>


#include <utils/test_utils.hpp>


template<typename PoolType, int NumThreads>
void testAppendMultipleBuffers(const int buffer_capacity, const int total_count) {
  omp_lock_t writelock;
  omp_init_lock(&writelock);
  omp_lock_t writelock2;
  omp_init_lock(&writelock2);
  omp_lock_t writelock3;
  omp_init_lock(&writelock3);
//  std::atomic_flag writelock = ATOMIC_FLAG_INIT;
//  std::atomic_flag writelock2 = ATOMIC_FLAG_INIT;
//  std::atomic_flag writelock3 = ATOMIC_FLAG_INIT;


  printf("TESTING: %d threads, pool thread %s buffer thread %s append with %d bufferSize and %d total counts from unlimited pool\n",
         NumThreads, PoolType::poolTS ? "SAFE" : "UNSAFE", PoolType::bufferTS ? "SAFE" : "UNSAFE", buffer_capacity, total_count);


  PoolType pool(buffer_capacity);

  std::vector<int> gold;
  std::vector<int> stored;

  int data = 0;
  unsigned int result = 0;

  int success = 0;
  int failure = 0;
  int swap = 0;
  int i = 0;

  auto buf_ptr = std::move(pool.tryAcquireBuffer());

#pragma omp parallel for num_threads(NumThreads) default(none) shared(buf_ptr, gold, stored, writelock, writelock2, writelock3) private(i, data, result) reduction(+:success, failure, swap)
  for (i = 0; i < total_count; ++i) {

//    while (writelock2.test_and_set());
    omp_set_lock(&writelock2);
    auto ptr = buf_ptr.get();
    omp_unset_lock(&writelock2);
//    writelock2.clear();

    data = static_cast<int>(i);
    result = ptr->append(&data, sizeof(int));

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
      auto new_buf_ptr = std::move(pool.tryAcquireBuffer());

//      while (writelock2.test_and_set());
      omp_set_lock(&writelock2);
      buf_ptr.swap(new_buf_ptr);
#pragma omp flush(new_buf_ptr)
      omp_unset_lock(&writelock2);
//      writelock2.clear();

      // process the old buffer
//      while (writelock3.test_and_set());
      omp_set_lock(&writelock3);
      stored.insert(stored.end(), new_buf_ptr->operator int*(), new_buf_ptr->operator int*() + new_buf_ptr->getSize() / sizeof(int));
      omp_unset_lock(&writelock3);
//      writelock3.clear();

      // and release - few threads doing this, and full.
      pool.releaseBuffer(std::move(new_buf_ptr));

    }

  }

  buf_ptr->force_lock_read();

  // compare unordered buffer content.
  stored.insert(stored.end(), buf_ptr->operator int*(), buf_ptr->operator int*() + buf_ptr->getSize() / sizeof(int));
  pool.releaseBuffer(std::move(buf_ptr));
  int stored_count = stored.size();


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
  omp_destroy_lock(&writelock);
  omp_destroy_lock(&writelock2);
  omp_destroy_lock(&writelock3);
}


template<typename PoolType>
void testPool(PoolType && pool, const std::string &name, int pool_threads, int buffer_threads) {

  using BufferType = typename PoolType::BufferType;


  printf("TESTING %s %s: pool threads %d, buffer threads %d\n", name.c_str(), (pool.isUnlimited() ? "GROW" : "FIXED"),  pool_threads, buffer_threads);

  printf("TEST acquire\n");
  int expected;
  int i = 0;
  int count = 0;
  int mx = pool.isUnlimited() ? 100 : pool.getCapacity();
#pragma omp parallel for num_threads(pool_threads) default(none) private(i) shared(pool, mx) reduction(+ : count)
  for (i = 0; i < mx; ++i) {
	  auto ptr = std::move(pool.tryAcquireBuffer());
    if (! ptr) {
      ++count;
    }
  }
  expected = 0;
  if (count != expected) printf("ERROR: number of failed attempt to acquire buffer should be %d, actual %d.  pool capacity %lu, remaining: %lu \n", expected, count, pool.getCapacity(), pool.getAvailableCount());
  pool.reset();

  printf("TEST acquire with growth\n");
  i = 0;
  count = 0;
  mx = pool.isUnlimited() ? 100 : pool.getCapacity();
#pragma omp parallel for num_threads(pool_threads) default(none) private(i) shared(pool, mx) reduction(+ : count)
  for (i = 0; i <= mx; ++i) {  // <= so we get 1 extra
		auto ptr = std::move(pool.tryAcquireBuffer());
	    if (! ptr) {
	      ++count;
	    }
  }
  expected = pool.isUnlimited() ? 0 : 1;
  if (count != expected) printf("ERROR: number of failed attempt to acquire buffer should be %d, actual %d.  pool remaining: %lu \n", expected, count, pool.getAvailableCount());

  pool.reset();

  printf("TEST release\n");
  count = 0;
  mx = pool.isUnlimited() ? 100 : pool.getCapacity();
  // first drain the pool
  for (i = 0; i < mx; ++i) {
    auto ptr = std::move(pool.tryAcquireBuffer());
  }
  // and create some dummy buffers to insert
  std::vector<typename PoolType::BufferPtrType> temp;
  for (i = 0; i < mx; ++i) {
    temp.push_back(std::move(std::unique_ptr<BufferType>(new BufferType(pool.getBufferCapacity()))));
    temp.push_back(std::move(std::unique_ptr<BufferType>(new BufferType(pool.getBufferCapacity()))));
  }
#pragma omp parallel for num_threads(pool_threads) default(none) shared(pool, mx, temp) private(i) reduction(+ : count)
  for (i = 0; i < mx * 2; ++i) {

    auto ptr = std::move(temp[i]);
    if (ptr) {
      ptr->block_and_flush();
      if (! pool.releaseBuffer(std::move(ptr))) {
        ++count; // failed release
      }
    }
  }
  expected = pool.isUnlimited() ? 0 : mx;  // unlimited or not, can only push back in as much as taken out.
  if (count != expected) printf("ERROR: number of failed attempt to release buffer should be %d, actual %d. pool remaining: %lu \n", expected, count, pool.getAvailableCount());
  pool.reset();
  temp.clear();

  printf("TEST access by multiple threads, each a separate buffer.\n");


  count = 0;
  int count1 = 0;
#pragma omp parallel num_threads(pool_threads) default(none) shared(pool) reduction(+ : count, count1)
  {
    int v = omp_get_thread_num() + 5;
    auto ptr = std::move(pool.tryAcquireBuffer());
    if (! (ptr->append(&v, sizeof(int)) & 0x1)) {
      ++count;
    }

    int u = ptr->operator int*()[0];
    if (v != u) {
      ++count1;
    }

    ptr->block_and_flush();
    pool.releaseBuffer(std::move(ptr));
  }
  if (count != 0) printf("ERROR: append failed\n");
  else if (count1 != 0) printf("ERROR: inserted and got back\n");
  pool.reset();

  printf("TEST access by multiple threads, all to same buffer.\n");


  auto ptr = std::move(pool.tryAcquireBuffer());
#pragma omp parallel num_threads(buffer_threads) default(none) shared(pool, ptr)
  {
    int v = 7;
    ptr->append(&v, sizeof(int));
  }

  bool same = true;
  for (int i = 0; i < buffer_threads ; ++i) {
    same &= ptr->operator int*()[i] == 7;
  }
  if (!same) printf("ERROR: inserted not same\n");
  pool.reset();


  omp_set_nested(1);

  printf("TEST all operations together\n");
#pragma omp parallel num_threads(pool_threads) default(none) shared(pool, pool_threads, buffer_threads)
  {
    // Id range is 0 to 100
    int iter;
    int j = 0;
    for (int i = 0; i < 100; ++i) {
      // acquire
      auto buf = pool.tryAcquireBuffer();
      while (!buf) {
        usleep(50);
        buf = pool.tryAcquireBuffer();
      }

      // access
      iter = rand() % 100;
#pragma omp parallel for num_threads(buffer_threads) default(none) shared(pool, buf, iter) private(j)
      for (j = 0; j < iter; ++j) {
        buf->append(&j, sizeof(int));
      }

      // random sleep
      usleep(rand() % 1000);
      // clear buffer
      buf->block_and_flush();

      if (buf->getSize() != sizeof(int) * iter)
        printf("ERROR: thread %d/%d buffer size is %ld, expected %lu\n", omp_get_thread_num(), pool_threads, buf->getSize(), sizeof(int) * iter);

      //release
      pool.releaseBuffer(std::move(buf));
      //if (i % 25 == 0)
//      printf("thread %d released buffer %d\n", omp_get_thread_num(), id);

    }
  }




};


int main(int argc, char** argv) {

  // construct, acquire, access, release

  //////////////  unbounded version

  /// thread unsafe.  test in single thread way.
//  bliss::io::BufferPool<bliss::concurrent::THREAD_UNSAFE> usPool(8192);


  testPool(std::move(bliss::io::BufferPool<bliss::concurrent::THREAD_UNSAFE>(8192)), "thread unsafe pool", 1, 1);




  testPool(std::move(bliss::io::BufferPool<bliss::concurrent::THREAD_SAFE>(8192)), "thread safe pool, thread safe buffer", 1, 1);
  testPool(std::move(bliss::io::BufferPool<bliss::concurrent::THREAD_SAFE>(8192)), "thread safe pool, thread safe buffer", 1, 2);
  testPool(std::move(bliss::io::BufferPool<bliss::concurrent::THREAD_SAFE>(8192)), "thread safe pool, thread safe buffer", 1, 3);
  testPool(std::move(bliss::io::BufferPool<bliss::concurrent::THREAD_SAFE>(8192)), "thread safe pool, thread safe buffer", 1, 4);
  testPool(std::move(bliss::io::BufferPool<bliss::concurrent::THREAD_SAFE>(8192)), "thread safe pool, thread safe buffer", 1, 8);

  testPool(std::move(bliss::io::BufferPool<bliss::concurrent::THREAD_SAFE>(8192)), "thread safe pool, thread safe buffer", 2, 1);
  testPool(std::move(bliss::io::BufferPool<bliss::concurrent::THREAD_SAFE>(8192)), "thread safe pool, thread safe buffer", 2, 2);
  testPool(std::move(bliss::io::BufferPool<bliss::concurrent::THREAD_SAFE>(8192)), "thread safe pool, thread safe buffer", 2, 3);
  testPool(std::move(bliss::io::BufferPool<bliss::concurrent::THREAD_SAFE>(8192)), "thread safe pool, thread safe buffer", 2, 4);

  testPool(std::move(bliss::io::BufferPool<bliss::concurrent::THREAD_SAFE>(8192)), "thread safe pool, thread safe buffer", 3, 1);
  testPool(std::move(bliss::io::BufferPool<bliss::concurrent::THREAD_SAFE>(8192)), "thread safe pool, thread safe buffer", 3, 2);
  testPool(std::move(bliss::io::BufferPool<bliss::concurrent::THREAD_SAFE>(8192)), "thread safe pool, thread safe buffer", 3, 3);






  testPool(std::move(bliss::io::BufferPool<bliss::concurrent::THREAD_SAFE, bliss::concurrent::THREAD_UNSAFE>(8192)), "thread safe pool, thread unsafe buffer", 1, 1);
  testPool(std::move(bliss::io::BufferPool<bliss::concurrent::THREAD_SAFE, bliss::concurrent::THREAD_UNSAFE>(8192)), "thread safe pool, thread unsafe buffer", 2, 1);
  testPool(std::move(bliss::io::BufferPool<bliss::concurrent::THREAD_SAFE, bliss::concurrent::THREAD_UNSAFE>(8192)), "thread safe pool, thread unsafe buffer", 3, 1);
  testPool(std::move(bliss::io::BufferPool<bliss::concurrent::THREAD_SAFE, bliss::concurrent::THREAD_UNSAFE>(8192)), "thread safe pool, thread unsafe buffer", 4, 1);




  // this one is not defined because it's not logical.  not compilable.
  // bliss::io::BufferPool<bliss::concurrent::THREAD_UNSAFE, bliss::concurrent::THREAD_SAFE> tsusPool(8192, 8);




  /////////////  fixed size version.



  /// testing buffer_pool. thread unsafe


  testPool(std::move(bliss::io::BufferPool<bliss::concurrent::THREAD_UNSAFE>(8, 8192)), "thread unsafe pool", 1, 1);




  /// testing buffer_pool. thread safe

  testPool(std::move(bliss::io::BufferPool<bliss::concurrent::THREAD_SAFE>(8, 8192)), "thread safe pool, thread safe buffer", 1, 1);
  testPool(std::move(bliss::io::BufferPool<bliss::concurrent::THREAD_SAFE>(8, 8192)), "thread safe pool, thread safe buffer", 1, 2);
  testPool(std::move(bliss::io::BufferPool<bliss::concurrent::THREAD_SAFE>(8, 8192)), "thread safe pool, thread safe buffer", 1, 3);
  testPool(std::move(bliss::io::BufferPool<bliss::concurrent::THREAD_SAFE>(8, 8192)), "thread safe pool, thread safe buffer", 1, 4);
  testPool(std::move(bliss::io::BufferPool<bliss::concurrent::THREAD_SAFE>(8, 8192)), "thread safe pool, thread safe buffer", 1, 8);

  testPool(std::move(bliss::io::BufferPool<bliss::concurrent::THREAD_SAFE>(8, 8192)), "thread safe pool, thread safe buffer", 2, 1);
  testPool(std::move(bliss::io::BufferPool<bliss::concurrent::THREAD_SAFE>(8, 8192)), "thread safe pool, thread safe buffer", 2, 2);
  testPool(std::move(bliss::io::BufferPool<bliss::concurrent::THREAD_SAFE>(8, 8192)), "thread safe pool, thread safe buffer", 2, 3);
  testPool(std::move(bliss::io::BufferPool<bliss::concurrent::THREAD_SAFE>(8, 8192)), "thread safe pool, thread safe buffer", 2, 4);

  testPool(std::move(bliss::io::BufferPool<bliss::concurrent::THREAD_SAFE>(8, 8192)), "thread safe pool, thread safe buffer", 3, 1);
  testPool(std::move(bliss::io::BufferPool<bliss::concurrent::THREAD_SAFE>(8, 8192)), "thread safe pool, thread safe buffer", 3, 2);
  testPool(std::move(bliss::io::BufferPool<bliss::concurrent::THREAD_SAFE>(8, 8192)), "thread safe pool, thread safe buffer", 3, 3);



  /// testing buffer_pool, thread safe pool, unsafe buffers.

  testPool(std::move(bliss::io::BufferPool<bliss::concurrent::THREAD_SAFE, bliss::concurrent::THREAD_UNSAFE>(8, 8192)), "thread safe pool, thread unsafe buffer", 1, 1);
  testPool(std::move(bliss::io::BufferPool<bliss::concurrent::THREAD_SAFE, bliss::concurrent::THREAD_UNSAFE>(8, 8192)), "thread safe pool, thread unsafe buffer", 2, 1);
  testPool(std::move(bliss::io::BufferPool<bliss::concurrent::THREAD_SAFE, bliss::concurrent::THREAD_UNSAFE>(8, 8192)), "thread safe pool, thread unsafe buffer", 3, 1);
  testPool(std::move(bliss::io::BufferPool<bliss::concurrent::THREAD_SAFE, bliss::concurrent::THREAD_UNSAFE>(8, 8192)), "thread safe pool, thread unsafe buffer", 4, 1);


  // this one is not defined because it's not logical.  not compilable.
  // bliss::io::BufferPool<bliss::concurrent::THREAD_UNSAFE, bliss::concurrent::THREAD_SAFE> tsusPool(8, 8192);





}
