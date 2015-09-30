/**
 * @file    check_tsan.cpp
 * @ingroup
 * @author  tpan
 * @brief
 * @details  testing some of the common race conditions that Thread Sanitizer points out.
 *           see https://code.google.com/p/thread-sanitizer/wiki/PopularDataRaces
 *
 * Copyright (c) 2015 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#include "bliss-config.hpp"
#include <omp.h>

#include <vector>
#include <unordered_map>
#include <sstream>
#include <iterator>
#include <algorithm>
#include <xmmintrin.h>
#include <atomic>
#include <mutex>


void testVectorBad() {
  printf("testVectorBad\n");

  std::vector<int> data(4);

#pragma omp parallel for OMP_SHARE_DEFAULT num_threads(4) shared(data)
  for (int i = 0; i < 1000000; ++i) {
    ++data[omp_get_thread_num()];
  }

  std::stringstream ss;
  std::ostream_iterator<int> osi(ss, ",");
  std::copy(data.begin(), data.end(), osi);

  printf("test vector for race result [%s]\n", ss.str().c_str());

}


void testMapBad() {
  printf("testMapBad\n");
  std::unordered_map<int, int> data;

  for (int i = 0; i < 4; ++i) {
    data[i] = 0;
  }


#pragma omp parallel for OMP_SHARE_DEFAULT num_threads(4) shared(data)
  for (int i = 0; i < 1000000; ++i) {
    ++data[omp_get_thread_num()];
  }

  std::stringstream ss;
  for (int i = 0; i < 4; ++i) {
    ss << "(" << i << "=" << data[i] << "), ";
  }

  printf("test vector for race result [%s]\n", ss.str().c_str());
}


//
//void testNotificationBad() {
//  printf("testNotificationBad\n");
//
//  bool done = false;
//  int count1, count2;
//
//
//#pragma omp parallel sections OMP_SHARE_DEFAULT num_threads(4) shared(done, count1, count2)
//  {
//#pragma omp section
//    {
//      while (!done || count1 < 10000000) {   // compiler optimization could cause this to become an infinite loop.
//        count1++;
//      }
//
//      printf("count1 = %d\n", count1);
//    }
//#pragma omp section
//    {
//      count2++;
//      done = true;
//      printf("count2 = %d\n", count2);
//    }
//  }
//
//  printf("count1 = %d\n", count1);
//  printf("count2 = %d\n", count2);
//}


void testNotificationRace() {
  printf("testNotificationRace\n");

  bool done = false;
  int count1 = 0, count2 = 0;


#pragma omp parallel sections OMP_SHARE_DEFAULT num_threads(2) shared(done, count1, count2)
  {
#pragma omp section
    {
      std::atomic_thread_fence(std::memory_order_acquire);
      while (!done) {
        count1++;
        std::atomic_thread_fence(std::memory_order_acquire);
      }

      printf("count1 = %d\n", count1);
    }
#pragma omp section
    {
      count2++;
      done = true;
      std::atomic_thread_fence(std::memory_order_release);
      printf("count2 = %d\n", count2);
    }
  }

  printf("count1 = %d\n", count1);
  printf("count2 = %d\n", count2);
}

void testNotificationGood() {
  printf("testNotificationGood\n");

  std::atomic<bool> done(false);
  int count1 = 0, count2 = 0;


#pragma omp parallel sections OMP_SHARE_DEFAULT num_threads(2) shared(done, count1, count2)
  {
#pragma omp section
    {
      while (!done.load(std::memory_order_relaxed)) {
        count1++;
      }

      printf("count1 = %d\n", count1);
    }
#pragma omp section
    {
      count2++;
      done.store(true, std::memory_order_relaxed);
      printf("count2 = %d\n", count2);
    }
  }

  printf("count1 = %d\n", count1);
  printf("count2 = %d\n", count2);
}



void testPublishBad() {
  printf("testPublishBad\n");

  std::pair<int, int>* obj = nullptr;


#pragma omp parallel sections OMP_SHARE_DEFAULT num_threads(2) shared(obj)
  {
#pragma omp section
    {
      obj = new std::pair<int, int>(1, 2);
    }
#pragma omp section
    {
      while (obj == nullptr) _mm_pause();
      printf("object is (%d,%d)\n", obj->first, obj->second);
    }
  }
}

void testPublishGood() {
  printf("testPublishGood\n");

  std::atomic<std::pair<int, int>* > obj(nullptr);


#pragma omp parallel sections OMP_SHARE_DEFAULT num_threads(2) shared(obj)
  {
#pragma omp section
    {
      obj.store(new std::pair<int, int>(1, 2), std::memory_order_release);
    }
#pragma omp section
    {
      while (obj.load(std::memory_order_acquire) == nullptr) _mm_pause();
      auto o = obj.load(std::memory_order_acquire);
      printf("object is (%d,%d)\n", o->first, o->second);
    }
  }
}


void testBitFieldBad() {
  printf("testBitFieldBad\n");
  struct bits {
      int a:4, b:4;
  } x;

#pragma omp parallel sections OMP_SHARE_DEFAULT num_threads(2) shared(x)
  {
#pragma omp section
    {
      x.a++;
    }
#pragma omp section
    {
      x.b++;
    }
  }
}



void testDoubleCheckedLockingBad() {
  printf("testDoubleCheckedLockingBad\n");

  bool inited = false;
  int x = 0;
  std::mutex mu;

#pragma omp parallel OMP_SHARE_DEFAULT num_threads(2) shared(inited, mu, x)
  {
    // May be called by multiple threads.
    if (!inited) {
      std::lock_guard<std::mutex> lock(mu);
      if (!inited) {
        x++;
      }
      inited = true;
    }
  }
  printf("x = %d\n", x);
}


void testSpinlockBad() {
  printf("testSpinLockBad\n");
// occasional race.   lock.clear does not prevent internals from being optimized into another's spin
  std::atomic_flag lock = ATOMIC_FLAG_INIT;
  int x = 0;

#pragma omp parallel OMP_SHARE_DEFAULT num_threads(128) shared(lock, x)
  {
    _mm_pause();
    while (lock.test_and_set(std::memory_order_acquire));
    x += 2;
    lock.clear(std::memory_order_release);
  }
  printf("x = %d\n", x);
}




int main(int argc, char** argv) {
  testVectorBad();

  testMapBad();

//  testNotificationBad();  infinite loop.

  testNotificationRace();
  testNotificationGood();

  testPublishBad();

  testBitFieldBad();

  testDoubleCheckedLockingBad();

  testSpinlockBad();
}


