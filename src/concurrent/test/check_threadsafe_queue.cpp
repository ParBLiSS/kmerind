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


#include <iostream>
#include <cstdio>
#include <unistd.h>
#include <queue>

#include "concurrent/threadsafe_queue.hpp"

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
      T result = 0;
      for (int i = 0; i < entries; ++i) {

        if (i % 2 == 0)
          q.push(std::move(T(i)));
        else {
          output = q.front();
          q.pop();
          result ^= output;
        }

      }
      printf("result: %d\n", result);
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

    bliss::concurrent::ThreadSafeQueue<T>& getQueue() {
      return q;
    }

    void run(const int &entries, bliss::concurrent::ThreadSafeQueue<T>& queue) {

      // check tryPop
      bool result = false;
      int i = 0;
      T output;
      for (i = 0; (i < entries) && !result; ++i) {
        result = queue.tryPop(output);
      }
      printf("TSQueue finished tryPop on empty at iteration %d\n", i);

      // now push and pop
      for (i = 0; (i < entries); ++i) {
        queue.tryPush(T(i));
      }
      result = true;
      for (i = 0; (i < entries); ++i) {
        result &= queue.tryPop(output);
      }
      printf("TSQueue finished tryPop on full at iteration %d\n", i);

      // now push in parallel
#pragma omp parallel for num_threads(nthreads) private(i) shared(queue, entries) default(none)
      for (i = 0; (i < entries); ++i) {
        queue.tryPush(T(i));
      }
      printf("TSQueue finished push in parallel with %lu entries\n", queue.size());

      // now pop in parallel
      int count = 0;
#pragma omp parallel for num_threads(nthreads) private(i) shared(queue, entries) default(none) reduction(+: count)
      for (i = 0; (i < entries); ++i) {
        T output;
        if (queue.tryPop(output)) ++count;
      }
      printf("TSQueue finished tryPop in parallel with %d entries\n", count);


      // now push then do waitAndPop in parallel
      for (i = 0; (i < entries); ++i) {
        queue.tryPush(T(i));
      }
      count = 0;
#pragma omp parallel for num_threads(nthreads) private(i) shared(queue, entries) default(none) reduction(+: count)
      for (i = 0; (i < entries); ++i) {
        T output;
        queue.waitAndPop(output);
        ++count;
      }
      printf("TSQueue finished waitAndPop in parallel with %d entries\n", count);

      // 1 thread, finished.
      if (nthreads == 1) return;


      // now have 1 thread wait and push, and 1 thread tryPop
      count = 0;
      int count2 = 0;
      int count3 = 0;

      bool done = false;

#pragma omp parallel sections num_threads(2) shared(queue, entries, done) default(none) reduction(+: count, count2, count3)
      {
#pragma omp section
        {
          usleep(1000);
          for (int i = 0; (i < entries); ++i) {
            queue.tryPush(T(i));
          }
          done = true;
#pragma omp flush(done)
        }

#pragma omp section
        {

          T output;
          while (!done) {
            if (queue.tryPop(output)) ++count;
            else ++count2;
#pragma omp flush(done)
          }

          // now empty it
          bool test = queue.tryPop(output);
          while (test) {
            test = queue.tryPop(output);
            ++count3;
          }

        }
      }
      printf("TSQueue finished 1 producer 1 consumer with %d success %d failed tryPops, and %d final tryPops\n", count, count2, count3);


      printf("TSQueue is empty? %s\n", (queue.empty() ? "yes" : "no"));

      // now have 1 thread wait and push, and 1 thread waitAndPop
      count = 0;
      count3 = 0;
      done = false;
#pragma omp parallel sections num_threads(2) shared(queue, entries, done) default(none) reduction(+: count, count3)
      {
#pragma omp section
        {
          usleep(1000);
          for (int i = 0; (i < entries); ++i) {
            queue.tryPush(T(i));
          }
          done = true;
#pragma omp flush(done)
        }

#pragma omp section
        {
          T output;
          while (!done) {
            queue.waitAndPop(output);
            ++count;
#pragma omp flush(done)
          }
          // now empty it
          bool test = queue.tryPop(output);
          while (test) {
            test = queue.tryPop(output);
            ++count3;
          }
        }
      }
      printf("TSQueue finished 1 producer 1 consumer with %d waitAndPops, and %d final tryPops\n", count, count3);

      printf("TSQueue is empty? %s\n", (queue.empty() ? "yes" : "no"));



      if (nthreads == 2) return;

      // now have n-1 thread wait and push, and 1 thread tryPop
      count = 0;
      count2 = 0;
      count3 = 0;
      done = false;
#pragma omp parallel sections num_threads(2) shared(queue, entries, done) default(none) reduction(+: count, count2,count3)
      {
#pragma omp section
        {
          usleep(1000);
          int nt = nthreads - 1;
#pragma omp parallel for default(none) num_threads(nt) shared(queue, entries)
          for (int i = 0; (i < entries); ++i) {
            queue.tryPush(T(i));
          }
          done = true;
#pragma omp flush(done)
        }

#pragma omp section
        {
          T output;
          while (!done) {
            if (queue.tryPop(output)) ++count;
            else ++count2;
#pragma omp flush(done)
          }
          // now empty it
          bool test = queue.tryPop(output);
          while (test) {
            test = queue.tryPop(output);
            ++count3;
          }
        }
      }
      printf("TSQueue finished %d producer 1 consumer with %d successful tryPops and %d failed, %d final\n", nthreads - 1, count, count2, count3);

      printf("TSQueue is empty? %s\n", (queue.empty() ? "yes" : "no"));


      // now have 1 thread wait and tryPush, and n-1 thread tryPop.
      count = 0;
      count2 = 0;
      count3 = 0;
      done = false;
#pragma omp parallel sections num_threads(2) shared(queue, entries, done) default(none) reduction(+: count, count2, count3)
      {
#pragma omp section
        {
          usleep(1000);
          for (int i = 0; (i < entries); ++i) {
            queue.tryPush(T(i));
          }
          done = true;
#pragma omp flush(done)
        }

#pragma omp section
        {
          int nt = nthreads -1;

          int counts[nt];
          int counts2[nt];
          int counts3[nt];


#pragma omp parallel num_threads(nt) default(none) shared(queue, entries, done, counts, counts2, counts3) reduction(+:count, count2, count3)
          {
            counts[omp_get_thread_num()] = 0;
            counts2[omp_get_thread_num()] = 0;
            counts3[omp_get_thread_num()] = 0;

#pragma omp single nowait
            {
              while (!done) {
  #pragma omp task default(none) shared(queue, counts, counts2)
                {
                  T output;
                  if (queue.tryPop(output))
                    ++counts[omp_get_thread_num()];
                  else
                    ++counts2[omp_get_thread_num()];
                }
  #pragma omp flush(done)
              }


              while (!queue.empty()) {
  #pragma omp task default(none) shared(queue, counts3)
                {
                  T output;
                  if (queue.tryPop(output))
                    ++counts3[omp_get_thread_num()];
                }
              }
            }

#pragma omp barrier
// require barrier else we accumulate too early.

//            printf("tid %d success %d fail %d final %d \n", omp_get_thread_num(), counts[omp_get_thread_num()], counts2[omp_get_thread_num()], counts3[omp_get_thread_num()]);
          count = counts[omp_get_thread_num()];
          count2 = counts2[omp_get_thread_num()];
          count3 = counts3[omp_get_thread_num()];

          }
        }
      }
      printf("TSQueue finished 1 producer %d consumer with %d successful tryPops and failed %d, final %d\n", nthreads - 1, count, count2, count3);
      printf("TSQueue is empty? %s\n", (queue.empty() ? "yes" : "no"));



      // now have n-1 thread wait and push, and 1 thread tryPop.  waiting is blocking. then the
      count = 0;
      count3 = 0;
      done = false;
#pragma omp parallel sections num_threads(2) shared(queue, entries, done) default(none) reduction(+: count, count3)
      {
#pragma omp section
        {
          usleep(1000);
          int nt = nthreads - 1;
#pragma omp parallel for default(none) num_threads(nt) shared(queue, entries)
          for (int i = 0; (i < entries); ++i) {
            queue.tryPush(T(i));
          }
          done = true;
#pragma omp flush(done)

        }

#pragma omp section
        {
          T output;
          while (!done) {
            queue.waitAndPop(output);
            ++count;
#pragma omp flush(done)
          }
          // now empty it
          bool test = queue.tryPop(output);
          while (test) {
            test = queue.tryPop(output);
            ++count3;
          }

        }
      }
      printf("TSQueue finished %d producer 1 consumer with %d waitAndPops, %d final \n", nthreads - 1, count, count3);
      printf("TSQueue is empty? %s\n", (queue.empty() ? "yes" : "no"));


      // now have 1 thread wait and push, and n-1 thread waitAndPop  - can't do this.  omp task (pop) creation finishes "done."  there could be a lot more tasks generated than
      //   there are elements because the produce thread is generating slower than consume threads.  the excess threads result in wait forever.
      // removed.

      if (nthreads == 3) return;

      // now have n/2 thread wait and push, and n/2 thread tryPop.
      count = 0;
      count2 = 0;
      count3 = 0;
      done = false;
#pragma omp parallel sections num_threads(2) shared(queue, entries, done) default(none) reduction(+: count, count2, count3)
      {
#pragma omp section
        {
          usleep(1000);
          int nt = nthreads / 2;
#pragma omp parallel for default(none) num_threads(nt) shared(queue, entries)
          for (int i = 0; (i < entries); ++i) {
            queue.tryPush(T(i));
          }
          done = true;
#pragma omp flush(done)
        }

#pragma omp section
        {
          int nt = (nthreads + 1) / 2;

          int counts[nt];
          int counts2[nt];
          int counts3[nt];

#pragma omp parallel num_threads(nt) default(none) shared(queue, entries, done, counts, counts2, counts3) reduction(+:count, count2, count3)
          {
            counts[omp_get_thread_num()] = 0;
            counts2[omp_get_thread_num()] = 0;
            counts3[omp_get_thread_num()] = 0;

#pragma omp single nowait
            {
              while (!done) {
  #pragma omp task default(none) shared(queue, counts, counts2)
                {
                  T output;
                  if (queue.tryPop(output))
                    ++counts[omp_get_thread_num()];
                  else
                    ++counts2[omp_get_thread_num()];
                }
  #pragma omp flush(done)
              }


              while (!queue.empty()) {
  #pragma omp task default(none) shared(queue, counts3)
                {
                  T output;
                  if (queue.tryPop(output))
                    ++counts3[omp_get_thread_num()];
                }
              }
            }

#pragma omp barrier
// require barrier else we accumulate too early.

//            printf("tid %d success %d fail %d final %d \n", omp_get_thread_num(), counts[omp_get_thread_num()], counts2[omp_get_thread_num()], counts3[omp_get_thread_num()]);
          count = counts[omp_get_thread_num()];
          count2 = counts2[omp_get_thread_num()];
          count3 = counts3[omp_get_thread_num()];

          }
        }
      }
      printf("TSQueue finished %d producer %d consumer with %d successful tryPops and failed %d, final %d\n", nthreads/2, (nthreads + 1)/2, count, count2, count3);
      printf("TSQueue is empty? %s\n", (queue.empty() ? "yes" : "no"));


    }

    void run2(int entries) {
      T output = 0;
      T result = 0;
      for (int i = 0; i < entries; ++i) {

        if (i % 2 == 0)
          q.tryPush(std::move(T(i)));
        else {
          q.tryPop(output);
          result ^= output;
        }
      }
      printf("result: %d\n", result);
    }
};

template<typename T, size_t N>
class FixedSizeQueueTest
{
  protected:
    bliss::concurrent::ThreadSafeQueue<T> q;
    int nthreads;

  public:
    FixedSizeQueueTest(int nThreads) : q(N), nthreads(nThreads) {};

    bliss::concurrent::ThreadSafeQueue<T>& getQueue() {
      return q;
    }

    void run(int entries, bliss::concurrent::ThreadSafeQueue<T>& queue) {
      // check tryPop

      // check tryPop
      bool result = false;
      int i = 0;
      T output;
      for (i = 0; (i < entries) && !result; ++i) {
        result = queue.tryPop(output);
      }
      printf("TSFixedQueue finished tryPop on empty at iteration %d\n", i);

      // now push and pop
      result = true;
      for (i = 0; (i < entries) && result; ++i) {
        result = queue.tryPush(std::move(T(i)));
      }
      printf("TSFixedQueue finished tryPush until full at iteration %d\n", i-1);
      result = true;
      for (i = 0; (i < entries) && result; ++i) {
        result = queue.tryPop(output);
      }
      printf("TSFixedQueue finished tryPop from full at iteration %d\n", i-1);

      // now push in parallel
#pragma omp parallel for num_threads(nthreads) private(i) shared(queue, entries) default(none)
      for (i = 0; (i < entries); ++i) {
        queue.tryPush(T(i));
      }
      printf("TSFixedQueue finished tryPush in parallel with %lu entries\n", queue.size());

      // now pop in parallel
      int count = 0;
#pragma omp parallel for num_threads(nthreads) private(i) shared(queue, entries) default(none) reduction(+: count)
      for (i = 0; (i < entries); ++i) {
        T output;
        if (queue.tryPop(output)) ++count;
      }
      printf("TSFixedQueue finished tryPop in parallel with %d entries\n", count);


      // now push then do waitAndPop in parallel
      for (i = 0; (i < entries); ++i) {
        queue.tryPush(T(i));
      }
      count = 0;
#pragma omp parallel for num_threads(nthreads) private(i) shared(queue, entries) default(none) reduction(+: count)
      for (i = 0; (i < queue.getMaxSize()); ++i) {
        T output;
        queue.waitAndPop(output);
        ++count;
      }
      printf("TSFixedQueue finished waitAndPop in parallel with %d entries\n", count);

      // 1 thread, finished.
      if (nthreads == 1) return;


      // now have 1 thread tryPush, and 1 thread tryPop
      count = 0;
      int count2 = 0;
      int count3 = 0;

      bool done = false;

#pragma omp parallel sections num_threads(2) shared(queue, entries, done) default(none) reduction(+: count, count2, count3)
      {
#pragma omp section
        {
          usleep(1000);
          for (int i = 0; (i < entries); ++i) {
            queue.tryPush(T(i));
          }
          done = true;
#pragma omp flush(done)
        }

#pragma omp section
        {

          T output;
          while (!done) {
            if (queue.tryPop(output)) ++count;
            else ++count2;

            if (count %1000) usleep(20);

#pragma omp flush(done)
          }


          // now empty it
          bool test = queue.tryPop(output);
          while (test) {
            test = queue.tryPop(output);
            ++count3;
          }
        }
      }
      printf("TSFixedQueue finished 1 producer 1 consumer with %d success %d failed tryPops, and %d final tryPops\n", count, count2, count3);


      printf("TSFixedQueue is empty? %s\n", (queue.empty() ? "yes" : "no"));

      // now have 1 thread tryPush, and 1 thread waitAndPop
      count = 0;
      count3 = 0;
      done = false;
#pragma omp parallel sections num_threads(2) shared(queue, entries, done) default(none) reduction(+: count, count3)
      {
#pragma omp section
        {
          usleep(1000);
          for (int i = 0; (i < entries); ++i) {
            queue.tryPush(T(i));
          }
          done = true;
#pragma omp flush(done)
        }

#pragma omp section
        {
          T output;
          while (!done) {
            queue.waitAndPop(output);
            ++count;

            if (count %1000) usleep(20);

#pragma omp flush(done)
          }
          // now empty it
          bool test = queue.tryPop(output);
          while (test) {
            test = queue.tryPop(output);
            ++count3;
          }
        }
      }
      printf("TSFixedQueue finished 1 producer 1 consumer with %d waitAndPops, and %d final tryPops\n", count, count3);

      printf("TSFixedQueue is empty? %s\n", (queue.empty() ? "yes" : "no"));




      // now have 1 thread waitAndPush, and 1 thread tryPop
      count = 0;
      count2 = 0;
      count3 = 0;

      done = false;

#pragma omp parallel sections num_threads(2) shared(queue, entries, done) default(none) reduction(+: count, count2, count3)
      {
#pragma omp section
        {
          usleep(1000);
          for (int i = 0; (i < entries); ++i) {
            queue.waitAndPush(T(i));
          }
          done = true;
#pragma omp flush(done)
        }

#pragma omp section
        {

          T output;
          while (!done) {
            if (queue.tryPop(output)) ++count;
            else ++count2;

            if (count %1000) usleep(20);

#pragma omp flush(done)
          }


          // now empty it
          bool test = queue.tryPop(output);
          while (test) {
            test = queue.tryPop(output);
            ++count3;
          }
        }
      }
      printf("TSFixedQueue finished 1 producer waitAndPush 1 consumer with %d success %d failed tryPops, and %d final tryPops\n", count, count2, count3);


      printf("TSFixedQueue is empty? %s\n", (queue.empty() ? "yes" : "no"));

      // now have 1 thread waitAndPush, and 1 thread waitAndPop
      count = 0;
      count3 = 0;
      done = false;
#pragma omp parallel sections num_threads(2) shared(queue, entries, done) default(none) reduction(+: count, count3)
      {
#pragma omp section
        {
          usleep(1000);
          for (int i = 0; (i < entries); ++i) {
            queue.waitAndPush(T(i));
          }
          done = true;
#pragma omp flush(done)
        }

#pragma omp section
        {
          T output;
          while (!done) {
            queue.waitAndPop(output);
            ++count;

            if (count %1000) usleep(20);

#pragma omp flush(done)
          }
          // now empty it
          bool test = queue.tryPop(output);
          while (test) {
            test = queue.tryPop(output);
            ++count3;
          }
        }
      }
      printf("TSFixedQueue finished 1 producer waitAndPush 1 consumer with %d waitAndPops, and %d final tryPops\n", count, count3);

      printf("TSFixedQueue is empty? %s\n", (queue.empty() ? "yes" : "no"));

      if (nthreads == 2) return;

      // now have n-1 thread wait and push, and 1 thread tryPop
      count = 0;
      count2 = 0;
      count3 = 0;
      done = false;
#pragma omp parallel sections num_threads(2) shared(queue, entries, done) default(none) reduction(+: count, count2,count3)
      {
#pragma omp section
        {
          usleep(1000);
          int nt = nthreads - 1;
#pragma omp parallel for default(none) num_threads(nt) shared(queue, entries)
          for (int i = 0; (i < entries); ++i) {
            queue.tryPush(T(i));
          }
          done = true;
#pragma omp flush(done)
        }

#pragma omp section
        {
          T output;
          while (!done) {
            if (queue.tryPop(output)) ++count;
            else ++count2;

            if (count %1000) usleep(20);

#pragma omp flush(done)
          }
          // now empty it
          bool test = queue.tryPop(output);
          while (test) {
            test = queue.tryPop(output);
            ++count3;
          }
        }
      }
      printf("TSFixedQueue finished %d producer 1 consumer with %d successful tryPops and %d failed, %d final\n", nthreads - 1, count, count2, count3);

      printf("TSFixedQueue is empty? %s\n", (queue.empty() ? "yes" : "no"));


      // now have 1 thread wait and push, and n-1 thread tryPop.
      count = 0;
      count2 = 0;
      count3 = 0;
      done = false;
#pragma omp parallel sections num_threads(2) shared(queue, entries, done) default(none) reduction(+: count, count2, count3)
      {
#pragma omp section
        {
          usleep(1000);
          for (int i = 0; (i < entries); ++i) {
            queue.tryPush(T(i));
          }
          done = true;
#pragma omp flush(done)
        }

#pragma omp section
        {
          int nt = nthreads -1;

          int counts[nt];
          int counts2[nt];
          int counts3[nt];


#pragma omp parallel num_threads(nt) default(none) shared(queue, entries, done, counts, counts2, counts3) reduction(+:count, count2, count3)
          {
            counts[omp_get_thread_num()] = 0;
            counts2[omp_get_thread_num()] = 0;
            counts3[omp_get_thread_num()] = 0;

#pragma omp single nowait
            {
              while (!done) {
  #pragma omp task default(none) shared(queue, counts, counts2)
                {
                  T output;
                  if (queue.tryPop(output))
                    ++counts[omp_get_thread_num()];
                  else
                    ++counts2[omp_get_thread_num()];

                  if (counts[omp_get_thread_num()] %1000) usleep(20);
                }
  #pragma omp flush(done)
              }


              while (!queue.empty()) {
  #pragma omp task default(none) shared(queue, counts3)
                {
                  T output;
                  if (queue.tryPop(output))
                    ++counts3[omp_get_thread_num()];
                }
              }
            }

#pragma omp barrier
// require barrier else we accumulate too early.

//            printf("tid %d success %d fail %d final %d \n", omp_get_thread_num(), counts[omp_get_thread_num()], counts2[omp_get_thread_num()], counts3[omp_get_thread_num()]);
          count = counts[omp_get_thread_num()];
          count2 = counts2[omp_get_thread_num()];
          count3 = counts3[omp_get_thread_num()];

          }
        }
      }
      printf("TSFixedQueue finished 1 producer %d consumer with %d successful tryPops and failed %d, final %d\n", nthreads - 1, count, count2, count3);
      printf("TSFixedQueue is empty? %s\n", (queue.empty() ? "yes" : "no"));



      // now have n-1 thread wait and push, and 1 thread tryPop.  waiting is blocking. then the
      count = 0;
      count3 = 0;
      done = false;
#pragma omp parallel sections num_threads(2) shared(queue, entries, done) default(none) reduction(+: count, count3)
      {
#pragma omp section
        {
          usleep(1000);
          int nt = nthreads - 1;
#pragma omp parallel for default(none) num_threads(nt) shared(queue, entries)
          for (int i = 0; (i < entries); ++i) {
            queue.tryPush(T(i));
          }
          done = true;
#pragma omp flush(done)

        }

#pragma omp section
        {
          T output;
          while (!done) {
            queue.waitAndPop(output);
            ++count;

            if (count %1000) usleep(20);

#pragma omp flush(done)
          }
          // now empty it
          bool test = queue.tryPop(output);
          while (test) {
            test = queue.tryPop(output);
            ++count3;
          }

        }
      }
      printf("TSFixedQueue finished %d producer 1 consumer with %d waitAndPops, %d final \n", nthreads - 1, count, count3);
      printf("TSFixedQueue is empty? %s\n", (queue.empty() ? "yes" : "no"));


      // now have 1 thread wait and push, and n-1 thread waitAndPop  - can't do this.  omp task (pop) creation finishes "done."  there could be a lot more tasks generated than
      //   there are elements because the produce thread is generating slower than consume threads.  the excess threads result in wait forever.
      // removed.

      // now have n-1 thread wait and push, and 1 thread tryPop
      count = 0;
      count2 = 0;
      count3 = 0;
      done = false;
#pragma omp parallel sections num_threads(2) shared(queue, entries, done) default(none) reduction(+: count, count2,count3)
      {
#pragma omp section
        {
          usleep(1000);
          int nt = nthreads - 1;
#pragma omp parallel for default(none) num_threads(nt) shared(queue, entries)
          for (int i = 0; (i < entries); ++i) {
            queue.waitAndPush(T(i));
          }
          done = true;
#pragma omp flush(done)
        }

#pragma omp section
        {
          T output;
          while (!done) {
            if (queue.tryPop(output)) ++count;
            else ++count2;

            if (count %1000) usleep(20);

#pragma omp flush(done)
          }
          // now empty it
          bool test = queue.tryPop(output);
          while (test) {
            test = queue.tryPop(output);
            ++count3;
          }
        }
      }
      printf("TSFixedQueue finished %d wait producer 1 consumer with %d successful tryPops and %d failed, %d final\n", nthreads - 1, count, count2, count3);

      printf("TSFixedQueue is empty? %s\n", (queue.empty() ? "yes" : "no"));


      // now have 1 thread wait and push, and n-1 thread tryPop.
      count = 0;
      count2 = 0;
      count3 = 0;
      done = false;
#pragma omp parallel sections num_threads(2) shared(queue, entries, done) default(none) reduction(+: count, count2, count3)
      {
#pragma omp section
        {
          usleep(1000);
          for (int i = 0; (i < entries); ++i) {
            queue.waitAndPush(T(i));
          }
          done = true;
#pragma omp flush(done)
        }

#pragma omp section
        {
          int nt = nthreads -1;

          int counts[nt];
          int counts2[nt];
          int counts3[nt];


#pragma omp parallel num_threads(nt) default(none) shared(queue, entries, done, counts, counts2, counts3) reduction(+:count, count2, count3)
          {
            counts[omp_get_thread_num()] = 0;
            counts2[omp_get_thread_num()] = 0;
            counts3[omp_get_thread_num()] = 0;

#pragma omp single nowait
            {
              while (!done) {
  #pragma omp task default(none) shared(queue, counts, counts2)
                {
                  T output;
                  if (queue.tryPop(output))
                    ++counts[omp_get_thread_num()];
                  else
                    ++counts2[omp_get_thread_num()];

                  if (counts[omp_get_thread_num()] %1000) usleep(20);
                }
  #pragma omp flush(done)
              }


              while (!queue.empty()) {
  #pragma omp task default(none) shared(queue, counts3)
                {
                  T output;
                  if (queue.tryPop(output))
                    ++counts3[omp_get_thread_num()];
                }
              }
            }

#pragma omp barrier
// require barrier else we accumulate too early.

//            printf("tid %d success %d fail %d final %d \n", omp_get_thread_num(), counts[omp_get_thread_num()], counts2[omp_get_thread_num()], counts3[omp_get_thread_num()]);
          count = counts[omp_get_thread_num()];
          count2 = counts2[omp_get_thread_num()];
          count3 = counts3[omp_get_thread_num()];

          }
        }
      }
      printf("TSFixedQueue finished 1 wait producer %d consumer with %d successful tryPops and failed %d, final %d\n", nthreads - 1, count, count2, count3);
      printf("TSFixedQueue is empty? %s\n", (queue.empty() ? "yes" : "no"));



      // now have n-1 thread wait and push, and 1 thread tryPop.  waiting is blocking. then the
      count = 0;
      count3 = 0;
      done = false;
#pragma omp parallel sections num_threads(2) shared(queue, entries, done) default(none) reduction(+: count, count3)
      {
#pragma omp section
        {
          usleep(1000);
          int nt = nthreads - 1;
#pragma omp parallel for default(none) num_threads(nt) shared(queue, entries)
          for (int i = 0; (i < entries); ++i) {
            queue.waitAndPush(T(i));
          }
          done = true;
#pragma omp flush(done)

        }

#pragma omp section
        {
          T output;
          while (!done) {
            queue.waitAndPop(output);
            ++count;

            if (count %1000) usleep(20);

#pragma omp flush(done)
          }
          // now empty it
          bool test = queue.tryPop(output);
          while (test) {
            test = queue.tryPop(output);
            ++count3;
          }

        }
      }
      printf("TSFixedQueue finished %d wait producer 1 consumer with %d waitAndPops, %d final \n", nthreads - 1, count, count3);
      printf("TSFixedQueue is empty? %s\n", (queue.empty() ? "yes" : "no"));



      if (nthreads == 3) return;

      // now have n/2 thread wait and push, and n/2 thread tryPop.
      count = 0;
      count2 = 0;
      count3 = 0;
      done = false;
#pragma omp parallel sections num_threads(2) shared(queue, entries, done) default(none) reduction(+: count, count2, count3)
      {
#pragma omp section
        {
          usleep(1000);
          int nt = nthreads / 2;
#pragma omp parallel for default(none) num_threads(nt) shared(queue, entries)
          for (int i = 0; (i < entries); ++i) {
            queue.tryPush(T(i));
          }
          done = true;
#pragma omp flush(done)
        }

#pragma omp section
        {
          int nt = (nthreads + 1) / 2;

          int counts[nt];
          int counts2[nt];
          int counts3[nt];

#pragma omp parallel num_threads(nt) default(none) shared(queue, entries, done, counts, counts2, counts3) reduction(+:count, count2, count3)
          {
            counts[omp_get_thread_num()] = 0;
            counts2[omp_get_thread_num()] = 0;
            counts3[omp_get_thread_num()] = 0;

#pragma omp single nowait
            {
              while (!done) {
  #pragma omp task default(none) shared(queue, counts, counts2)
                {
                  T output;
                  if (queue.tryPop(output))
                    ++counts[omp_get_thread_num()];
                  else
                    ++counts2[omp_get_thread_num()];

                  if ((counts[omp_get_thread_num()] % 1000) == 0) usleep(20);
                }
  #pragma omp flush(done)
              }


              while (!queue.empty()) {
  #pragma omp task default(none) shared(queue, counts3)
                {
                  T output;
                  if (queue.tryPop(output))
                    ++counts3[omp_get_thread_num()];
                }
              }
            }

#pragma omp barrier
// require barrier else we accumulate too early.

//            printf("tid %d success %d fail %d final %d \n", omp_get_thread_num(), counts[omp_get_thread_num()], counts2[omp_get_thread_num()], counts3[omp_get_thread_num()]);
          count = counts[omp_get_thread_num()];
          count2 = counts2[omp_get_thread_num()];
          count3 = counts3[omp_get_thread_num()];

          }
        }
      }
      printf("TSFixedQueue finished %d producer %d consumer with %d successful tryPops and failed %d, final %d\n", nthreads/2, (nthreads + 1)/2, count, count2, count3);
      printf("TSFixedQueue is empty? %s\n", (queue.empty() ? "yes" : "no"));



      // now have n/2 thread wait and push, and n/2 thread tryPop.
      count = 0;
      count2 = 0;
      count3 = 0;
      done = false;
#pragma omp parallel sections num_threads(2) shared(queue, entries, done) default(none) reduction(+: count, count2, count3)
      {
#pragma omp section
        {
          usleep(1000);
          int nt = nthreads / 2;
#pragma omp parallel for default(none) num_threads(nt) shared(queue, entries)
          for (int i = 0; (i < entries); ++i) {
            queue.waitAndPush(T(i));
          }
          done = true;
#pragma omp flush(done)
        }

#pragma omp section
        {
          int nt = (nthreads + 1) / 2;

          int counts[nt];
          int counts2[nt];
          int counts3[nt];

#pragma omp parallel num_threads(nt) default(none) shared(queue, entries, done, counts, counts2, counts3) reduction(+:count, count2, count3)
          {
            counts[omp_get_thread_num()] = 0;
            counts2[omp_get_thread_num()] = 0;
            counts3[omp_get_thread_num()] = 0;

#pragma omp single nowait
            {
              while (!done) {
  #pragma omp task default(none) shared(queue, counts, counts2)
                {
                  T output;
                  if (queue.tryPop(output))
                    ++counts[omp_get_thread_num()];
                  else
                    ++counts2[omp_get_thread_num()];

                  if ((counts[omp_get_thread_num()] % 1000) == 0) usleep(20);
                }
  #pragma omp flush(done)
              }


              while (!queue.empty()) {
  #pragma omp task default(none) shared(queue, counts3)
                {
                  T output;
                  if (queue.tryPop(output))
                    ++counts3[omp_get_thread_num()];
                }
              }
            }

#pragma omp barrier
// require barrier else we accumulate too early.

//            printf("tid %d success %d fail %d final %d \n", omp_get_thread_num(), counts[omp_get_thread_num()], counts2[omp_get_thread_num()], counts3[omp_get_thread_num()]);
          count = counts[omp_get_thread_num()];
          count2 = counts2[omp_get_thread_num()];
          count3 = counts3[omp_get_thread_num()];

          }
        }
      }
      printf("TSFixedQueue finished %d wait producer %d consumer with %d successful tryPops and failed %d, final %d\n", nthreads/2, (nthreads + 1)/2, count, count2, count3);
      printf("TSFixedQueue is empty? %s\n", (queue.empty() ? "yes" : "no"));



    }

    void run2(int entries) {
      T output;
      T result = 0;
      for (int i = 0; i < entries; ++i) {

        if (i % 2 == 0)
          q.tryPush(std::move(T(i)));
        else
        {
          q.tryPop(output);
          result ^= output;

        }
      }
      printf("result: %d\n", result);
    }

};



int main(int argc, char** argv) {

  int elements = 100000;
  std::chrono::high_resolution_clock::time_point t1, t2;
  std::chrono::duration<double> time_span;

  omp_set_nested(1);
  omp_set_dynamic(0);


  typedef QueueTest<int> QTestType;
  for (int i = 1; i <=4; ++i) {
    QTestType qt(i);
    qt.run(elements, qt.getQueue());

    t1 = std::chrono::high_resolution_clock::now();
    qt.run2(elements);
    t2 = std::chrono::high_resolution_clock::now();
    time_span =
        std::chrono::duration_cast<std::chrono::duration<double>>(
            t2 - t1);

    printf("queueTest nthread %d time %f\n", i, time_span.count());

  }

  typedef FixedSizeQueueTest<int, 1000> FSQTestType;
  for (int i = 1; i <=4; ++i) {
    FSQTestType qt(i);
    qt.run(elements,  qt.getQueue());

    t1 = std::chrono::high_resolution_clock::now();
    qt.run2(elements);
    t2 = std::chrono::high_resolution_clock::now();
    time_span =
        std::chrono::duration_cast<std::chrono::duration<double>>(
            t2 - t1);

    printf("fixed queueTest size %d nthread %d time %f\n", 1000, i, time_span.count());

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
