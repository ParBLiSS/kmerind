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
#include <chrono>
#include <iterator>  // for ostream_iterator
#include <iostream>   // for cout
#include <algorithm>  // for sort

int main(int argc, char** argv) {

  typedef int valType;
  int typeSize = sizeof(valType);
  valType val;


  printf(" check thread local buffer.\n");

  bliss::io::Buffer<bliss::concurrent::THREAD_UNSAFE> tlBuffer(8192);

  // check insertion.
  int start = -1, end = -1;
  int i = 0;
  for (i = 0; i < 9000; ++i) {
    val = static_cast<valType>(i);
    if (! tlBuffer.append(&val, sizeof(valType))) {
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
  for (i = 0; i < 9000; ++i) {
    if (tlBuffer.isFull()) {
      if (start == -1) start = i;
    }
    val = static_cast<valType>(i);
    tlBuffer.append(&val, sizeof(valType));
  }
  end = i;
  printf("isFULL from %d to %d!\n", start, end);

  const valType* temp = reinterpret_cast<const valType*>(tlBuffer.getData());
  for (int j = 0; j < (8192/typeSize); ++j) {
    assert(temp[j] == j);
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



  temp = reinterpret_cast<const valType*>(tlBuffer.getData());
  for (int j = 0; j < (8192/typeSize); ++j) {
    assert(temp[j] == j);
  }
  printf("CONTENT Passes check after move\n");




  // check concurrent insert
  printf(" check thread safe buffer.\n");


  // check clear
  int bufferSize = 32768;
  bliss::io::Buffer<bliss::concurrent::THREAD_SAFE> tsBuffer4(bufferSize);


  // check insertion.
  std::vector<int> failed;

  i = 0;
  int fail = 0;
#pragma omp parallel for num_threads(4) default(none) shared(tsBuffer4, failed) private(i, val) reduction(+: fail)
  for (i = 0; i < 9000; ++i) {
    val = static_cast<valType>(i);

    if (! tsBuffer4.append(&val, sizeof(valType))) {
      ++fail;
#pragma omp critical
      failed.push_back(i);
    }
  }


//  temp = reinterpret_cast<const valType*>(tsBuffer4.getData());
//
//  std::vector<valType> content;
//
//  for (int j = 0; j < (bufferSize/typeSize); ++j) {
//    content.push_back(temp[j]);
//    printf("[%d %d], ", j, temp[j]);
//    if ((j % 8) == 7) printf("\n");
//  }
//  printf("\n\n");
//  std::sort(content.begin(), content.end());
//  for (int j = 0; j < (bufferSize/typeSize); ++j) {
//    printf("%d, ", content[j]);
//    if ((j % 8 ) == 7) printf("\n");
//  }
//  printf("\n\n");

  std::sort(failed.begin(), failed.end());

  printf("FAILED insertion ids:");

  std::ostream_iterator<int> oit(std::cout, ",");
  std::copy(failed.begin(), failed.end(), oit);
  std::cout << std::endl;

  printf("concurrent insert failed %d times \n", fail);



  // check clear
  tsBuffer4.clear();
  if (tsBuffer4.getSize() == 0) printf("empty\n");
  else printf("not empty\n");


  // check isFull
  i = 0;
  fail = 0;
#pragma omp parallel for num_threads(4) default(none) shared(tsBuffer4) private(i, val) reduction(+ : fail)
  for (i = 0; i < 9000; ++i) {
    if (tsBuffer4.isFull()) {
      ++fail;
    } else {
      val = static_cast<valType>(i);

      if (!tsBuffer4.append(&val, sizeof(valType))) {
        ++fail;
      }
    }
  }

  printf("isFULL for %d iterations\n", fail);





  ////////////// timing.  the insert before this is to warm up.
  int cap = 10000;

  bliss::io::Buffer<bliss::concurrent::THREAD_SAFE> tsBuffer5(cap);



    std::chrono::high_resolution_clock::time_point t1, t2;
    std::chrono::duration<double> time_span;

    for (int k = 0; k < 3; ++k) {

      tsBuffer5.clear();
      t1 = std::chrono::high_resolution_clock::now();

      // check insertion.
      i = 0;
      fail = 0;
    #pragma omp parallel for num_threads(4) default(none) shared(tsBuffer5, cap, typeSize) private(i, val) reduction(+:fail)
      for (i = 0; i < cap/typeSize; ++i) {
        val = static_cast<valType>(i);

        if (! tsBuffer5.append(&val, sizeof(valType))) {
          ++fail;
        }
      }
      t2 = std::chrono::high_resolution_clock::now();

      time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
      printf("locking append can't insert (full) %d iterations!  duration = %f\n", fail, time_span.count());


      tsBuffer5.clear();
      t1 = std::chrono::high_resolution_clock::now();

      // check insertion.
      i = 0;
      fail = 0;
    #pragma omp parallel for num_threads(4) default(none) shared(tsBuffer5, cap, typeSize) private(i, val) reduction(+:fail)
      for (i = 0; i < cap/typeSize; ++i) {
        val = static_cast<valType>(i);

        if (! tsBuffer5.append_lockfree(&val, sizeof(valType))) {
          ++fail;
        }
      }
      t2 = std::chrono::high_resolution_clock::now();

      time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
      printf("lockfree append can't insert (full) %d iterations!  duration = %f\n", fail, time_span.count());

    }
}
