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


#include "io/buffer.hpp"
#include "omp.h"
#include <cassert>

int main(int argc, char** argv) {


  printf(" check thread local buffer.\n");

  bliss::io::Buffer<bliss::concurrent::THREAD_UNSAFE> tlBuffer(8192);

  // check insertion.
  int start = -1, end = -1;
  int i = 0;
  for (i = 0; i < 9000; i += 4) {
    if (! tlBuffer.append(&i, sizeof(int))) {
       if (start == -1) start = i;
    }
  }
  end = i;
  printf("can't insert (full) from %d to %d!\n", start, end);

  // check clear
  tlBuffer.clear();
  if (tlBuffer.getSize() == 0) printf("empty\n");
  else printf("not empty\n");


  // check isFull
  start = -1;
  end = -1;
  i = 0;
  for (i = 0; i < 9000; i += 4) {
    if (tlBuffer.isFull()) {
      if (start == -1) start = i;
    }
    tlBuffer.append(&i, sizeof(int));
  }
  end = i;
  printf("isFULL from %d to %d!\n", start, end);

  const int* temp = reinterpret_cast<const int*>(tlBuffer.getData());
  for (unsigned int j = 0; j < (8192/sizeof(int)); ++j) {
    assert(temp[j] == static_cast<int>(j*4));
  }
  printf("CONTENT Passes check after insert\n");


  printf("check move ctor from thread unsafe to thread unsafe.\n");
  bliss::io::Buffer<bliss::concurrent::THREAD_UNSAFE> tlBuffer2(std::move(tlBuffer));
  printf("before  size %lu, capacity = %lu, pointer = %p; after   size %lu, capacity = %lu, pointer = %p\n",
         tlBuffer.getSize(), tlBuffer.getCapacity(), tlBuffer.getData(),
         tlBuffer2.getSize(), tlBuffer2.getCapacity(), tlBuffer2.getData());

  printf(" check move = from thread unsafe to thread unsafe.\n");
  tlBuffer = std::move(tlBuffer2);
  printf("before  size %lu, capacity = %lu, pointer = %p; after   size %lu, capacity = %lu, pointer = %p\n",
         tlBuffer2.getSize(), tlBuffer2.getCapacity(), tlBuffer2.getData(),
         tlBuffer.getSize(), tlBuffer.getCapacity(), tlBuffer.getData());




  printf(" check move ctor from thread unsafe to thread safe.\n");
  bliss::io::Buffer<bliss::concurrent::THREAD_SAFE> tsBuffer(std::move(tlBuffer));
  printf("before  size %lu, capacity = %lu, pointer = %p; after   size %lu, capacity = %lu, pointer = %p\n",
         tlBuffer.getSize(), tlBuffer.getCapacity(), tlBuffer.getData(),
         tsBuffer.getSize(), tsBuffer.getCapacity(), tsBuffer.getData());


  printf(" check move ctor from thread safe to thread safe.\n");
  bliss::io::Buffer<bliss::concurrent::THREAD_SAFE> tsBuffer2(std::move(tsBuffer));
  printf("before  size %lu, capacity = %lu, pointer = %p; after   size %lu, capacity = %lu, pointer = %p\n",
         tsBuffer.getSize(), tsBuffer.getCapacity(), tsBuffer.getData(),
         tsBuffer2.getSize(), tsBuffer2.getCapacity(), tsBuffer2.getData());

  printf(" check move = from thread safe to thread safe.\n");
  tsBuffer = std::move(tsBuffer2);
  printf("before  size %lu, capacity = %lu, pointer = %p; after   size %lu, capacity = %lu, pointer = %p\n",
         tsBuffer2.getSize(), tsBuffer2.getCapacity(), tsBuffer2.getData(),
         tsBuffer.getSize(), tsBuffer.getCapacity(), tsBuffer.getData());


  printf(" check move ctor from thread safe to thread unsafe.\n");
  bliss::io::Buffer<bliss::concurrent::THREAD_UNSAFE> tlBuffer3(std::move(tsBuffer));
  printf("before  size %lu, capacity = %lu, pointer = %p; after   size %lu, capacity = %lu, pointer = %p\n",
         tsBuffer.getSize(), tsBuffer.getCapacity(), tsBuffer.getData(),
         tlBuffer3.getSize(), tlBuffer3.getCapacity(), tlBuffer3.getData());


  printf(" check move = from thread unsafe to thread safe.\n");
  tsBuffer = std::move(tlBuffer3);
  printf("before  size %lu, capacity = %lu, pointer = %p; after   size %lu, capacity = %lu, pointer = %p\n",
         tlBuffer3.getSize(), tlBuffer3.getCapacity(), tlBuffer3.getData(),
         tsBuffer.getSize(), tsBuffer.getCapacity(), tsBuffer.getData());

  printf(" check move = from thread unsafe to thread safe.\n");
  tlBuffer = std::move(tsBuffer);
  printf("before  size %lu, capacity = %lu, pointer = %p; after   size %lu, capacity = %lu, pointer = %p\n",
         tsBuffer.getSize(), tsBuffer.getCapacity(), tsBuffer.getData(),
         tlBuffer.getSize(), tlBuffer.getCapacity(), tlBuffer.getData());



  temp = reinterpret_cast<const int*>(tlBuffer.getData());
  for (unsigned int j = 0; j < (8192/sizeof(int)); ++j) {
    assert(temp[j] == static_cast<int>(j*4));
  }
  printf("CONTENT Passes check after move\n");




  // check concurrent insert
  printf(" check thread local buffer.\n");


  // check clear
  bliss::io::Buffer<bliss::concurrent::THREAD_SAFE> tsBuffer4(8192);


  // check insertion.
  end = -1;
  std::atomic<int> start2(-1);
  i = 0;
#pragma omp parallel for num_threads(4) default(none) shared(tsBuffer4, start2, end) private(i)
  for (i = 0; i < 9000; i += 4) {
    if (! tsBuffer4.append(&i, sizeof(int))) {
       start2.compare_exchange_strong(end, i);
    }
  }
  end = i;
  printf("can't insert (full) from %d to %d!\n", start2.load(), end);
//
//  temp = reinterpret_cast<const int*>(tsBuffer4.getData());
//  for (unsigned int j = 0; j < (8192/sizeof(int)); ++j) {
//    printf("%d %d\n", j, temp[j]);
//  }
//  printf("\n\n");
//

  // check clear
  tsBuffer4.clear();
  if (tsBuffer4.getSize() == 0) printf("empty\n");
  else printf("not empty\n");


  // check isFull
  start2 = -1;
  end = -1;
  i = 0;

#pragma omp parallel for num_threads(4) default(none) shared(tsBuffer4, start2, end) private(i)
  for (i = 0; i < 9000; i += 4) {
    if (tsBuffer4.isFull()) {
      start2.compare_exchange_strong(end, i);
    }
    tsBuffer4.append(&i, sizeof(int));
  }
  end = i;
  printf("isFULL from %d to %d!\n", start2.load(), end);


  temp = reinterpret_cast<const int*>(tsBuffer4.getData());
  for (unsigned int j = 0; j < (8192/sizeof(int)); ++j) {
    printf("[%d %d], ", j, temp[j]);
    if ((j % 8) == 7) printf("\n");
  }
  printf("\n\n");






}
