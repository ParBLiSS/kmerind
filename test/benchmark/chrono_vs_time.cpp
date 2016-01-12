/**
 * @file    chrono_vs_time.cpp
 * @ingroup
 * @author  tpan
 * @brief
 * @details
 *
 * Copyright (c) 2016 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */



#include <ctime>
#include <iostream>
#include <cstring>
#include "utils/timer.hpp"

timespec operator-(timespec end, timespec start)
{
  timespec temp;
  if ((end.tv_nsec-start.tv_nsec)<0) {
    temp.tv_sec = end.tv_sec-start.tv_sec-1;
    temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
  } else {
    temp.tv_sec = end.tv_sec-start.tv_sec;
    temp.tv_nsec = end.tv_nsec-start.tv_nsec;
  }
  return temp;
}

timespec operator+(timespec lhs, timespec rhs)
{
  timespec temp;
  temp.tv_sec = lhs.tv_sec + rhs.tv_sec;
  temp.tv_nsec = lhs.tv_nsec + rhs.tv_nsec;
  if (temp.tv_nsec >= 1000000000) {
    temp.tv_sec += 1;
    temp.tv_nsec -= 1000000000;
  }
  return temp;
}



int main(int argc, char** argv) {

  int iters = 10000000;

  if (argc > 2) {
    iters = atoi(argv[2]);
  }

  bool chron = true;
  if (argc > 1) {
    if (strncmp(argv[1], "-t", 2) == 0) chron = false;
  }
  double dummy = 0;
  if (chron) {


    TIMER_INIT(test);

    TIMER_LOOP_START(test);
    for (int i = 0; i < iters; ++i) {
      TIMER_LOOP_RESUME(test);

      dummy += 0.00001;

      TIMER_LOOP_PAUSE(test);
    }

    TIMER_LOOP_END(test, "chrono in loop", iters);
    TIMER_REPORT(test, 0);

  } else {

    timespec start, end, duration;

    duration.tv_sec = 0; duration.tv_nsec = 0;

    for (int i = 0; i < iters; ++i) {
      clock_gettime(CLOCK_MONOTONIC, &start);

      dummy += 0.00001;

      clock_gettime(CLOCK_MONOTONIC, &end);
      duration = duration + (end - start);
    }

    std::cout << "elapsed: " << duration.tv_sec << "s " << duration.tv_nsec << "ns" << std::endl;


  }
  printf("dummy %f\n", dummy);




  return 0;
}
