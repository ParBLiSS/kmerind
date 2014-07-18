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
#include "concurrent/threadsafe_queue.hpp"
#include "omp.h"
#include <cassert>
#include <chrono>
#include <cstdlib>   // for rand

template<typename BuffersType>
void testPool(BuffersType && buffers, const std::string &name, int nthreads) {

  printf("TESTING %s: ntargets = %lu, pool threads %d\n", name.c_str(), buffers.getSize(), nthreads);


  printf("TEST append until full\n");
  typedef typename BuffersType::BufferIdType BufferIdType;
  std::pair<bool, BufferIdType > result(false, -1);
  bliss::concurrent::ThreadSafeQueue<typename BuffersType::BufferIdType> fullBuffers;

  std::string data("this is a test.  this a test of the emergency broadcast system.  this is only a test. ");
  printf("test string is \"%s\", length %lu\n", data.c_str(), data.length());
  int id = 0;

  int i;
  int count = 0;
  int count2 = 0;
  int count3 = 0;
  int count4 = 0;



#pragma omp parallel for num_threads(nthreads) default(none) private(i, id, result) shared(buffers, data, fullBuffers) reduction(+ : count, count2, count3, count4)
  for (i = 0; i < 1000; ++i) {
    id = 0;
    //printf("insert %lu chars into %d\n", data.length(), id);

    result = buffers.append(data.c_str(), data.length(), id);

    if (result.first) {
      ++count;
    } else {
      ++count2;
    }

    if (result.second != -1) {
      ++count3;
    } else {
      if (!result.first) {
        ++count4;
      }
    }

    if (result.second != -1)
      fullBuffers.waitAndPush(result.second);
  }
  printf("Number of failed attempt to append to buffer is %d, success %d. full buffers size: %lu.  numFullBuffers = %d.  num failed append due to no buffer = %d\n", count2, count, fullBuffers.getSize(), count3, count4);

  printf("TEST release\n");
  count = 0;
//  printf("buffer ids 0 to %lu initially in use.\n", buffers.getSize() - 1);
//  printf("releasing: ");
#pragma omp parallel for num_threads(nthreads) default(none) private(id, result) shared(buffers, data, fullBuffers) reduction(+ : count)
  for (id = 0; id < 350; ++id) {
    try {
      result = fullBuffers.tryPop();
      if (result.first) {
//        printf("%d ", result.second);
        if (result.second != -1)
          buffers.releaseBuffer(result.second);
      }
    } catch(const std::invalid_argument & e)
    {
      ++count;
    }
  }
//;  printf("\n");
  int expected = 0;
  if (count != expected) printf("ERROR: number of failed attempt to release buffer should be %d, actual %d. buffers size: %lu \n", expected, count, fullBuffers.getSize());

  buffers.reset();

  printf("TEST all operations together\n");
  count3 = 0;
  count = 0;
  count2 = 0;
  //printf("full buffer: ");
#pragma omp parallel for num_threads(nthreads) default(none) private(id, result) shared(buffers, data) reduction(+ : count, count2, count3)
  for (i = 0; i < 1000; ++i) {
    id = 0;
    result = buffers.append(data.c_str(), data.length(), id);

    if (!result.first) {
      ++count;
    } else {
      ++count2;
    }

    if (result.second != -1) {
      usleep(300);
      ++count3;
      //printf("%d ", result.second);
      buffers.releaseBuffer(result.second);
    }
  }
  //printf("\n");
  printf("Number of failed attempt to append to buffer is %d, success %d. full buffers size: %lu, released %d\n", count, count2, fullBuffers.getSize(), count3);

};


int main(int argc, char** argv) {

  // construct, acquire, access, release


  /// thread unsafe.  test in single thread way.


  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_UNSAFE>(1, 2048)), "thread unsafe buffers", 1);
  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_UNSAFE>(2, 2048)), "thread unsafe buffers", 1);
  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_UNSAFE>(3, 2048)), "thread unsafe buffers", 1);
  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_UNSAFE>(4, 2048)), "thread unsafe buffers", 1);
  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_UNSAFE>(5, 2048)), "thread unsafe buffers", 1);
  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_UNSAFE>(6, 2048)), "thread unsafe buffers", 1);
  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_UNSAFE>(7, 2048)), "thread unsafe buffers", 1);
  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_UNSAFE>(8, 2048)), "thread unsafe buffers", 1);



  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_SAFE>(1, 2048)), "thread safe buffers", 1);
  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_SAFE>(1, 2048)), "thread safe buffers", 2);
  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_SAFE>(1, 2048)), "thread safe buffers", 3);
  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_SAFE>(1, 2048)), "thread safe buffers", 4);
  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_SAFE>(1, 2048)), "thread safe buffers", 8);

  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_SAFE>(2, 2048)), "thread safe buffers", 1);
  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_SAFE>(2, 2048)), "thread safe buffers", 2);
  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_SAFE>(2, 2048)), "thread safe buffers", 3);
  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_SAFE>(2, 2048)), "thread safe buffers", 4);

  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_SAFE>(3, 2048)), "thread safe buffers", 1);
  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_SAFE>(3, 2048)), "thread safe buffers", 2);
  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_SAFE>(3, 2048)), "thread safe buffers", 3);

  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_SAFE>(4, 2048)), "thread safe buffers", 1);
  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_SAFE>(4, 2048)), "thread safe buffers", 2);


  // this one is not defined because it's not logical.  not compilable.
  // bliss::io::BufferPool<bliss::concurrent::THREAD_UNSAFE, bliss::concurrent::THREAD_SAFE> tsusPool(8192, 8);






}
