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
  valType val;


  printf("THREAD LOCAL BUFFERS\n");

  bliss::io::Buffer<bliss::concurrent::THREAD_UNSAFE> tlBuffer(8192);

  printf("TEST insert\n");
  int start = -1, end = -1;
  unsigned int i;
  for (i = 0; i < 9000; ++i) {
    val = static_cast<valType>(i);
    if (! (tlBuffer.append(&val, sizeof(valType)).first) ) {
       if (start == -1) start = i;
    }
  }
  end = i;
  if (start != 8192/sizeof(valType)) printf("\tFAILED can't insert (full) from %d to %d!\n", start, end);
  const valType* temp = reinterpret_cast<const valType*>(tlBuffer.getData());
  for (int j = 0; j < int(8192/sizeof(valType)); ++j) {
    if (temp[j] != j) {
      printf("ERROR: result is not expected.  expected: %d actual %d\n", j, temp[j]);
      //assert(false);
    }
  }


  printf("TEST clear\n");
  tlBuffer.clear();
  if (tlBuffer.getFinalSize() != 0) printf("\tERROR: NOT empty\n");


  printf("TEST isFull\n");
  start = -1;
  end = -1;
  for (i = 0; i < 9000; ++i) {
    if (tlBuffer.isFull()) {
      if (start == -1) start = i;
    }
    val = static_cast<valType>(i);
    tlBuffer.append(&val, sizeof(valType));
  }
  end = i;
  if (start != 8192/sizeof(valType)) printf("\tERROR: isFULL from %d to %d!\n", start, end);


  temp = reinterpret_cast<const valType*>(tlBuffer.getData());
  for (int j = 0; j < int(8192/sizeof(valType)); ++j) {
    if (temp[j] != j) {
      //printf("ERROR: result is not expected.  expected: %d actual %d\n", j, temp[j]);
      assert(false);
    }
  }
  printf("CONTENT Passes check after insert\n");


  size_t bufferCapBefore = tlBuffer.getCapacity();
  size_t bufferCapAfter = 0;

  size_t bufferSizeBefore = tlBuffer.getApproximateSize();
  size_t bufferSizeAfter = 0;

  const void* bufferPtrBefore = tlBuffer.getData();
  const void* bufferPtrAfter = nullptr;

  printf("TEST move ctor from thread unsafe to thread unsafe.\n");
  bliss::io::Buffer<bliss::concurrent::THREAD_UNSAFE> tlBuffer2(std::move(tlBuffer));
  if (tlBuffer.getCapacity()  !=  bufferCapAfter  ||
      tlBuffer.getApproximateSize()      != bufferSizeAfter ||
      tlBuffer.getData()      !=  bufferPtrAfter  ||
      tlBuffer2.getCapacity() !=  bufferCapBefore  ||
      tlBuffer2.getApproximateSize()     != bufferSizeBefore ||
      tlBuffer2.getData()     !=  bufferPtrBefore)
  printf("before  size %u, capacity = %u, pointer = %p; after   size %u, capacity = %u, pointer = %p\n",
         tlBuffer.getApproximateSize(), tlBuffer.getCapacity(), tlBuffer.getData(),
         tlBuffer2.getApproximateSize(), tlBuffer2.getCapacity(), tlBuffer2.getData());

  printf("TEST move = from thread unsafe to thread unsafe.\n");
  tlBuffer = std::move(tlBuffer2);
  if (tlBuffer2.getCapacity()  !=  bufferCapAfter  ||
      tlBuffer2.getApproximateSize()      != bufferSizeAfter ||
      tlBuffer2.getData()      !=  bufferPtrAfter  ||
      tlBuffer.getCapacity() !=  bufferCapBefore  ||
      tlBuffer.getApproximateSize()     != bufferSizeBefore ||
      tlBuffer.getData()     !=  bufferPtrBefore)
  printf("before  size %u, capacity = %u, pointer = %p; after   size %u, capacity = %u, pointer = %p\n",
         tlBuffer2.getApproximateSize(), tlBuffer2.getCapacity(), tlBuffer2.getData(),
         tlBuffer.getApproximateSize(), tlBuffer.getCapacity(), tlBuffer.getData());




  printf("TEST move ctor from thread unsafe to thread safe.\n");
  bliss::io::Buffer<bliss::concurrent::THREAD_SAFE> tsBuffer(std::move(tlBuffer));
  if (tlBuffer.getCapacity()  !=  bufferCapAfter  ||
      tlBuffer.getApproximateSize()      != bufferSizeAfter ||
      tlBuffer.getData()      !=  bufferPtrAfter  ||
      tsBuffer.getCapacity() !=  bufferCapBefore  ||
      tsBuffer.getApproximateSize()     != bufferSizeBefore ||
      tsBuffer.getData()     !=  bufferPtrBefore)
  printf("before  size %u, capacity = %u, pointer = %p; after   size %u, capacity = %u, pointer = %p\n",
         tlBuffer.getApproximateSize(), tlBuffer.getCapacity(), tlBuffer.getData(),
         tsBuffer.getApproximateSize(), tsBuffer.getCapacity(), tsBuffer.getData());


  printf("TEST move ctor from thread safe to thread safe.\n");
  bliss::io::Buffer<bliss::concurrent::THREAD_SAFE> tsBuffer2(std::move(tsBuffer));
  if (tsBuffer.getCapacity()  !=  bufferCapAfter  ||
      tsBuffer.getApproximateSize()      != bufferSizeAfter ||
      tsBuffer.getData()      !=  bufferPtrAfter  ||
      tsBuffer2.getCapacity() !=  bufferCapBefore  ||
      tsBuffer2.getApproximateSize()     != bufferSizeBefore ||
      tsBuffer2.getData()     !=  bufferPtrBefore)
  printf("before  size %u, capacity = %u, pointer = %p; after   size %u, capacity = %u, pointer = %p\n",
         tsBuffer.getApproximateSize(), tsBuffer.getCapacity(), tsBuffer.getData(),
         tsBuffer2.getApproximateSize(), tsBuffer2.getCapacity(), tsBuffer2.getData());

  printf("TEST move = from thread safe to thread safe.\n");
  tsBuffer = std::move(tsBuffer2);
  if (tsBuffer2.getCapacity()  !=  bufferCapAfter  ||
      tsBuffer2.getApproximateSize()      != bufferSizeAfter ||
      tsBuffer2.getData()      !=  bufferPtrAfter  ||
      tsBuffer.getCapacity() !=  bufferCapBefore  ||
      tsBuffer.getApproximateSize()     != bufferSizeBefore ||
      tsBuffer.getData()     !=  bufferPtrBefore)
  printf("before  size %u, capacity = %u, pointer = %p; after   size %u, capacity = %u, pointer = %p\n",
         tsBuffer2.getApproximateSize(), tsBuffer2.getCapacity(), tsBuffer2.getData(),
         tsBuffer.getApproximateSize(), tsBuffer.getCapacity(), tsBuffer.getData());


  printf("TEST move ctor from thread safe to thread unsafe.\n");
  bliss::io::Buffer<bliss::concurrent::THREAD_UNSAFE> tlBuffer3(std::move(tsBuffer));
  if (tsBuffer.getCapacity()  !=  bufferCapAfter  ||
      tsBuffer.getApproximateSize()      != bufferSizeAfter ||
      tsBuffer.getData()      !=  bufferPtrAfter  ||
      tlBuffer3.getCapacity() !=  bufferCapBefore  ||
      tlBuffer3.getApproximateSize()     != bufferSizeBefore ||
      tlBuffer3.getData()     !=  bufferPtrBefore)
  printf("before  size %u, capacity = %u, pointer = %p; after   size %u, capacity = %u, pointer = %p\n",
         tsBuffer.getApproximateSize(), tsBuffer.getCapacity(), tsBuffer.getData(),
         tlBuffer3.getApproximateSize(), tlBuffer3.getCapacity(), tlBuffer3.getData());


  printf("TEST move = from thread unsafe to thread safe.\n");
  tsBuffer = std::move(tlBuffer3);
  if (tlBuffer3.getCapacity()  !=  bufferCapAfter  ||
      tlBuffer3.getApproximateSize()      != bufferSizeAfter ||
      tlBuffer3.getData()      !=  bufferPtrAfter  ||
      tsBuffer.getCapacity() !=  bufferCapBefore  ||
      tsBuffer.getApproximateSize()     != bufferSizeBefore ||
      tsBuffer.getData()     !=  bufferPtrBefore)
  printf("before  size %u, capacity = %u, pointer = %p; after   size %u, capacity = %u, pointer = %p\n",
         tlBuffer3.getApproximateSize(), tlBuffer3.getCapacity(), tlBuffer3.getData(),
         tsBuffer.getApproximateSize(), tsBuffer.getCapacity(), tsBuffer.getData());

  printf("TEST move = from thread unsafe to thread safe.\n");
  tlBuffer = std::move(tsBuffer);
  if (tsBuffer.getCapacity()  !=  bufferCapAfter  ||
      tsBuffer.getApproximateSize()      != bufferSizeAfter ||
      tsBuffer.getData()      !=  bufferPtrAfter  ||
      tlBuffer.getCapacity() !=  bufferCapBefore  ||
      tlBuffer.getApproximateSize()     != bufferSizeBefore ||
      tlBuffer.getData()     !=  bufferPtrBefore)
  printf("before  size %u, capacity = %u, pointer = %p; after   size %u, capacity = %u, pointer = %p\n",
         tsBuffer.getApproximateSize(), tsBuffer.getCapacity(), tsBuffer.getData(),
         tlBuffer.getApproximateSize(), tlBuffer.getCapacity(), tlBuffer.getData());



  temp = reinterpret_cast<const valType*>(tlBuffer.getData());
  for (int j = 0; j < int(8192/sizeof(valType)); ++j) {
    assert(temp[j] == j);
  }
  printf("CONTENT Passes check after move\n");




  // check concurrent insert
  printf("THREAD SAFE BUFFERS\n");


  // check clear
  int bufferSize = 32768;
  bliss::io::Buffer<bliss::concurrent::THREAD_SAFE> tsBuffer4(bufferSize);


  // check insertion.
  std::vector<int> failed;

  i = 0;
  unsigned int fail = 0;
#pragma omp parallel for num_threads(4) default(none) shared(tsBuffer4, failed) private(i, val) reduction(+: fail)
  for (i = 0; i < 9000; ++i) {
    val = static_cast<valType>(i);

    if (! (tsBuffer4.append(&val, sizeof(valType))).first) {
      ++fail;
#pragma omp critical
      failed.push_back(i);
    }
  }

  if (fail != 9000-(bufferSize/sizeof(valType))) {

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

  printf("ERROR: concurrent insert failed %d times \n", fail);
  }


  printf("TEST: clear\n");
  tsBuffer4.clear();
  if (tsBuffer4.getFinalSize() != 0) printf("\tERROR: NOT empty\n");



  printf("TEST: par blocked\n");
  i = 0;
  unsigned int failedFull = 0;
  fail = 0;
  tsBuffer4.block();
#pragma omp parallel for num_threads(4) default(none) shared(tsBuffer4) private(i, val) reduction(+ : fail, failedFull)
  for (i = 0; i < 100; ++i) {

    if (tsBuffer4.isFull()) {
      ++failedFull;
    } else {
      val = static_cast<valType>(i);

      if (!tsBuffer4.append(&val, sizeof(valType)).first) {
        ++fail;
      }
    }
  }

  if (failedFull != 0 || fail != 100) printf("\tERROR: BLOCKED for %d iterations, with fail due to full %d iterations \n", fail, failedFull);


  printf("TEST: par Full\n");
  tsBuffer4.clear();
  tsBuffer4.unblock();
  i = 0;
  fail = 0;
  failedFull = 0;
#pragma omp parallel for num_threads(4) default(none) shared(tsBuffer4) private(i, val) reduction(+ : fail, failedFull)
  for (i = 0; i < 9000; ++i) {
    if (tsBuffer4.isFull()) {
      ++failedFull;
    } else {
      val = static_cast<valType>(i);

      if (!tsBuffer4.append(&val, sizeof(valType)).first) {
        ++fail;
      }
    }
  }

  if ((failedFull + fail) != (9000 - (bufferSize/sizeof(valType)))) printf("\tERROR: isFULL for %d iterations, failed append for %d iterations\n", failedFull, fail);





  ////////////// timing.  the insert before this is to warm up.
  int cap = 1000000;

  bliss::io::Buffer<bliss::concurrent::THREAD_SAFE> tsBuffer5(cap);



    std::chrono::high_resolution_clock::time_point t1, t2;
    std::chrono::duration<double> time_span;

    for (int k = 0; k < 3; ++k) {

      tsBuffer5.clear();
      t1 = std::chrono::high_resolution_clock::now();

      // check insertion.
      i = 0;
      fail = 0;
    #pragma omp parallel for num_threads(4) default(none) shared(tsBuffer5, cap) private(i, val) reduction(+:fail)
      for (i = 0; i < cap/sizeof(valType) + 1; ++i) {
        val = static_cast<valType>(i);

        if (! tsBuffer5.append(&val, sizeof(valType)).first) {
          ++fail;
        }
      }
      t2 = std::chrono::high_resolution_clock::now();

      time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
      int expected = 1;
      printf("locking append can't insert (full) %d iterations, expected %d iterations.  duration = %f\n", fail, expected, time_span.count());


//      tsBuffer5.clear();
//      t1 = std::chrono::high_resolution_clock::now();
//
//      // check insertion.
//      i = 0;
//      fail = 0;
//    #pragma omp parallel for num_threads(4) default(none) shared(tsBuffer5, cap, typeSize) private(i, val) reduction(+:fail)
//      for (i = 0; i < cap/typeSize + 1; ++i) {
//        val = static_cast<valType>(i);
//
//        if (! tsBuffer5.append_lockfree(&val, sizeof(valType))) {
//          ++fail;
//        }
//      }
//      t2 = std::chrono::high_resolution_clock::now();
//
//      time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
//      expected = 1;
//      printf("lockfree append can't insert (full) %d iterations, expected %d iterations.  duration = %f\n", fail, expected, time_span.count());

    }
}
