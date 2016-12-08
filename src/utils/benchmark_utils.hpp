/*
 * Copyright 2015 Georgia Institute of Technology
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/**
 * @file    system_utils.hpp
 * @ingroup
 * @author  tpan
 * @brief
 * @details
 *
 */
#ifndef SRC_WIP_SYSTEM_UTILS_HPP_
#define SRC_WIP_SYSTEM_UTILS_HPP_

#include "bliss-logger_config.hpp"

#include "utils/timer.hpp"
#include "utils/memory_usage.hpp"

#if BL_BENCHMARK == 1

  #define BL_BENCH_INIT(title)                            BL_TIMER_INIT(title);  BL_MEMUSE_INIT(title); do { BL_MEMUSE_MARK(title, "begin");  } while (0)
  #define BL_BENCH_RESET(title)                           do { BL_TIMER_RESET(title); BL_MEMUSE_RESET(title); } while (0)
  #define BL_BENCH_LOOP_START(title, id)                      do { BL_TIMER_LOOP_START(title, id); } while (0)
  #define BL_BENCH_LOOP_RESUME(title, id)                     do { BL_TIMER_LOOP_RESUME(title, id); } while (0)
  #define BL_BENCH_LOOP_PAUSE(title, id)                      do { BL_TIMER_LOOP_PAUSE(title, id); } while (0)
  #define BL_BENCH_LOOP_END(title, id, name, n_elem)          do { BL_TIMER_LOOP_END(title, id, name, n_elem); BL_MEMUSE_MARK(title, name); } while (0)
  #define BL_BENCH_START(title)                           do { BL_TIMER_START(title); } while (0)
  #define BL_BENCH_COLLECTIVE_START(title, name, comm)    do { BL_TIMER_COLLECTIVE_START(title, name, comm); } while (0)
  #define BL_BENCH_COLLECTIVE_END(title, name, n_elem, comm)    do { BL_TIMER_COLLECTIVE_END(title, name, n_elem, comm); BL_MEMUSE_MARK(title, name); } while (0)
  #define BL_BENCH_END(title, name, n_elem)               do { BL_TIMER_END(title, name, n_elem); BL_MEMUSE_MARK(title, name); } while (0)
  #define BL_BENCH_REPORT(title, rank)                    do { BL_TIMER_REPORT(title); BL_MEMUSE_REPORT(title); } while (0)
  #define BL_BENCH_REPORT_MPI(title, rank, comm)          do { BL_TIMER_REPORT_MPI(title, comm); BL_MEMUSE_REPORT_MPI(title, comm); } while (0)
  #define BL_BENCH_REPORT_NAMED(title, name)                    do { BL_TIMER_REPORT_NAMED(title, name); BL_MEMUSE_REPORT_NAMED(title, name); } while (0)
  #define BL_BENCH_REPORT_MPI_NAMED(title, name, comm)          do { BL_TIMER_REPORT_MPI_NAMED(title, name, comm); BL_MEMUSE_REPORT_MPI_NAMED(title, name, comm); } while (0)

#else

  #define BL_BENCH_INIT(title)
  #define BL_BENCH_RESET(title)
  #define BL_BENCH_LOOP_START(title, id)
  #define BL_BENCH_LOOP_RESUME(title, id)
  #define BL_BENCH_LOOP_PAUSE(title, id)
  #define BL_BENCH_LOOP_END(title, id, name, n_elem)
  #define BL_BENCH_START(title)
  #define BL_BENCH_COLLECTIVE_START(title, name, comm)
  #define BL_BENCH_COLLECTIVE_END(title, name, n_elem, comm)
  #define BL_BENCH_END(title, name, n_elem)
  #define BL_BENCH_REPORT(title, rank)
  #define BL_BENCH_REPORT_MPI(title, rank, comm)
  #define BL_BENCH_REPORT_NAMED(title, name)
  #define BL_BENCH_REPORT_MPI_NAMED(title, name, comm)
#endif

#endif /* SRC_WIP_SYSTEM_UTILS_HPP_ */

