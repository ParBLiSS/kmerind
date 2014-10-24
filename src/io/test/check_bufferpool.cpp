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
#include <cstdlib>   // for rand

template<typename PoolType>
void testPool(PoolType && pool, const std::string &name, int pool_threads, int buffer_threads) {

  printf("TESTING %s %s: pool threads %d, buffer threads %d\n", name.c_str(), (pool.isFixedSize() ? "FIXED" : "GROW"),  pool_threads, buffer_threads);

  printf("TEST acquire\n");
  int expected;
  int i = 0;
  typename PoolType::IdType id = 0;
  int count = 0;
#pragma omp parallel for num_threads(pool_threads) default(none) private(i, id) shared(pool) reduction(+ : count)
  for (i = 0; i < 100; ++i) {
    if (!pool.tryAcquireBuffer().first) {
      ++count;
    }
  }
  expected = pool.getCapacity();
  expected = std::max(0, 100 - expected);
  if (count != expected) printf("ERROR: number of failed attempt to acquire buffer should be %d, actual %d.  pool capacity %d, size: %d \n", expected, count, pool.getCapacity(), pool.getSize());
  pool.reset();

  printf("TEST acquire with growth\n");
  id = 0;
  i = 0;
  count = 0;
#pragma omp parallel for num_threads(pool_threads) default(none) private(i, id) shared(pool) reduction(+ : count)
  for (i = 0; i < 150; ++i) {
    if (!pool.tryAcquireBuffer().first) {
      ++count;
    }
  }
  expected = pool.getCapacity();
  expected = std::max(0, 150 - expected);
  if (count != expected) printf("ERROR: number of failed attempt to acquire buffer should be %d, actual %d.  pool size: %d \n", expected, count, pool.getSize());



  printf("TEST release\n");
  count = 0;
  id = 0;
#pragma omp parallel for num_threads(pool_threads) default(none) private(id) shared(pool) reduction(+ : count)
  for (id = 0; id < 350; ++id) {
    try {
      pool[id].block();
      pool.releaseBuffer(id);
    } catch (const std::out_of_range & e)
    {
      ++count;
    }
  }
  expected = 350 - pool.getSize();
  if (count != expected) printf("ERROR: number of failed attempt to release buffer should be %d, actual %d. pool size: %d \n", expected, count, pool.getSize());

#pragma omp parallel for num_threads(pool_threads) default(none) private(id) shared(pool) reduction(+ : count)
  for (id = 0; id < pool.getSize(); ++id) {
    pool[id].unblock();
  }

  printf("TEST access by multiple threads, each a separate buffer.\n");
#pragma omp parallel num_threads(pool_threads) default(none) shared(pool)
  {
    int v = omp_get_thread_num() + 5;
    pool[omp_get_thread_num()+ 1].append(&v, sizeof(int));
  }
#pragma omp parallel num_threads(pool_threads) default(none) shared(pool)
  {
    int u = reinterpret_cast<const int*>(pool[omp_get_thread_num() + 1].getData())[0];
    if (omp_get_thread_num() + 5 != u) printf("ERROR: %d value inserted was %d, getting %d\n", omp_get_thread_num(), omp_get_thread_num() + 5, u);
  }
  pool.reset();

  printf("TEST access by multiple threads, all to same buffer.\n");
#pragma omp parallel num_threads(buffer_threads) default(none) shared(pool)
  {
    int v = omp_get_thread_num() + 7;
    pool[1].append(&v, sizeof(int));
  }
  printf("values inserted were ");
  for (int i = 0; i < buffer_threads ; ++i) {
    int u = reinterpret_cast<const int*>(pool[1].getData())[i];
    printf(", %d", u);
  }
  printf("\n");

  pool.reset();
  omp_set_nested(1);

  printf("TEST all operations together\n");
#pragma omp parallel num_threads(pool_threads) default(none) shared(pool, pool_threads, buffer_threads)
  {
    // Id range is 0 to 100
    typename PoolType::IdType id;
    int iter;
    int j = 0;
    for (int i = 0; i < 100; ++i) {
      // acquire
      auto newBuf = pool.tryAcquireBuffer();
      while (!newBuf.first) usleep(50);
      //if (i % 25 == 0)
//      printf("thread %d acquired buffer %d\n", omp_get_thread_num(), id);
      id = newBuf.second;
      // access
      iter = rand() % 100;
#pragma omp parallel for num_threads(buffer_threads) default(none) shared(pool, id, iter) private(j)
      for (j = 0; j < iter; ++j) {
        pool[id].append(&j, sizeof(int));
      }

      // random sleep
      usleep(rand() % 1000);
      if (pool[id].getSize() != sizeof(int) * iter)
        printf("ERROR: thread %d/%d buffer %d size is %lu, expected %lu\n", omp_get_thread_num(), pool_threads, id, pool[id].getSize(), sizeof(int) * iter);

      // clear buffer
      pool[id].block();
      pool[id].clear();

      //release
      pool.releaseBuffer(id);
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
