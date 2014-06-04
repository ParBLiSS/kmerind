/**
 * @file		test_threadsafe_queue.cpp
 * @ingroup
 * @author	tpan
 * @brief
 * @details
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */



/**
 * range_test.cpp
 * Test range class
 *  Created on: Feb 18, 2014
 *      Author: tpan
 */

#include <iostream>
#include <cstdio>

#include "concurrent/threadsafe_queue.hpp"
#include "concurrent/threadsafe_fixedsize_queue.hpp"

#include "omp.h"

using namespace bliss::concurrent;



template<typename T>
class STDQueueTest
{
  protected:
    std::queue<T> q;
    int nthreads;

  public:

    void run2(int entries) {
      T output;
      for (int i = 0; i < entries; ++i) {

        if (i % 2 == 0)
          q.push(T(i));
        else {
          T v = q.front();
          q.pop();
        }

      }
    }
};




template<typename T>
class QueueTest
{
  protected:
    bliss::concurrent::ThreadSafeQueue<T> q;
    int nthreads;

  public:
    QueueTest(int nThreads) : q(), nthreads(nThreads) {};

    void run(int entries) {

#pragma omp parallel for num_threads(nthreads) shared(std::cout, q, entries) default(none)
      for (int i = 0; i < entries; ++i) {
        int tid = omp_get_thread_num();

        if (i < 100) {
          T output;
          bool result = q.tryPop(output);
          std::cout << "thread " << tid << " try pop: " << result << std::endl;
        }

        q.push(T(i));

        if (i % 3 == 0) {
          T output;
          q.waitAndPop(output);
          std::cout << "thread " << tid << " pop result " << output << std::endl;
        }

      }

      std::cout << "size = " << q.size() << std::endl;
    }
    void run2(int entries) {
      T output;
      for (int i = 0; i < entries; ++i) {

        if (i % 2 == 0)
          q.push(T(i));
        else
            q.tryPop(output);
      }
    }
};

template<typename T>
class FixedSizeQueueTest
{
  protected:
    bliss::concurrent::ThreadSafeFixedSizeQueue<T> q;
    int nthreads;
    size_t size;

  public:
    FixedSizeQueueTest(int nThreads, int _size) : q(_size), nthreads(nThreads), size(_size) {};

    void run(int entries) {

#pragma omp parallel for default(none) shared(std::cout, q, entries) num_threads(nthreads)
      for (int i = 0; i < entries; ++i) {
        int tid = omp_get_thread_num();

        if (i < 100) {
          T output;
          bool result = q.tryPop(output);
          std::cout << "thread " << tid << " try pop: " << result << std::endl;
        }

        if (i >= 100 && i <= (size + 200)) {
          bool result = q.tryPush(T(i));
          std::cout << "thread " << tid << " try push: " << result << std::endl;
        }


        if (i % 2 == 0) {
          T output;
          q.waitAndPop(output);
          std::cout << "thread " << tid << " pop result " << output << std::endl;
        } else {
          q.waitAndPush(T(i));
        }
      }

      std::cout << "size = " << q.size() << std::endl;
    }

    void run2(int entries) {
      T output;
      for (int i = 0; i < entries; ++i) {

        if (i % 2 == 0)
          q.tryPush(T(i));
        else
            q.tryPop(output);
      }
    }

};



int main(int argc, char** argv) {

  int elements = 10000;
  std::chrono::high_resolution_clock::time_point t1, t2;
  std::chrono::duration<double> time_span;

  typedef QueueTest<int> QTestType;
  for (int i = 1; i <=4; ++i) {
    QTestType qt(i);
    qt.run(elements);

    t1 = std::chrono::high_resolution_clock::now();
    qt.run2(elements);
    t2 = std::chrono::high_resolution_clock::now();
    time_span =
        std::chrono::duration_cast<std::chrono::duration<double>>(
            t2 - t1);

    printf("queueTest nthread %d time %f\n", i, time_span.count());

  }

  int size = 1000;
  typedef FixedSizeQueueTest<int> FSQTestType;
  for (int i = 1; i <=4; ++i) {
    FSQTestType qt(i, size);
    qt.run(elements);

    t1 = std::chrono::high_resolution_clock::now();
    qt.run2(elements);
    t2 = std::chrono::high_resolution_clock::now();
    time_span =
        std::chrono::duration_cast<std::chrono::duration<double>>(
            t2 - t1);

    printf("fixed queueTest size %d nthread %d time %f\n", size, i, time_span.count());

  }


  typedef STDQueueTest<int> SQTestType;
  SQTestType qt;

  t1 = std::chrono::high_resolution_clock::now();
  qt.run2(elements);
  t2 = std::chrono::high_resolution_clock::now();
  time_span =
      std::chrono::duration_cast<std::chrono::duration<double>>(
          t2 - t1);

  printf("std queueTest time %f\n", time_span.count());


};
