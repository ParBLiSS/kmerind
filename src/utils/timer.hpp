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
 * @file    timer.hpp
 * @ingroup
 * @author  tpan
 * @brief
 * @details
 *
 */
#ifndef SRC_UTILS_TIMER_HPP_
#define SRC_UTILS_TIMER_HPP_

#include "bliss-logger_config.hpp"

#if BENCHMARK == 1

#include <chrono>   // clock
//#include <cstdio>   // printf
#include <vector>
#include <string>
#include <algorithm>  // std::min

#include "io/mxx_support.hpp"

//TODO: make this class friendly.

struct Timer {
    std::chrono::steady_clock::time_point t1, t2;
    std::vector<std::string> names;
    std::vector<double> durations; 
    std::vector<double> counts;
    std::chrono::duration<double> time_span;
};


#define TIMER_INIT(timing)      Timer timing##_timer;


#define TIMER_RESET(timing)     do {  timing##_timer.names.clear(); \
                                       timing##_timer.durations.clear(); \
                                       timing##_timer.count.clear(); \
                                 } while (0)

#define TIMER_LOOP_START(timing)     do { timing##_timer.time_span = std::chrono::duration<double>::zero(); } while (0)
#define TIMER_LOOP_RESUME(timing)    do { timing##_timer.t1 = std::chrono::steady_clock::now(); } while (0)
#define TIMER_LOOP_PAUSE(timing)     do { timing##_timer.t2 = std::chrono::steady_clock::now(); \
                                      timing##_timer.time_span += (std::chrono::duration_cast<std::chrono::duration<double> >(timing##_timer.t2 - timing##_timer.t1)); } while (0)
#define TIMER_LOOP_END(timing, name, n_elem) do { timing##_timer.names.push_back(name); \
                                             timing##_timer.durations.push_back(timing##_timer.time_span.count()); \
                                             timing##_timer.counts.push_back(n_elem); \
                                        } while (0)

#define TIMER_START(timing)     do { timing##_timer.t1 = std::chrono::steady_clock::now(); } while (0)
#define TIMER_COLLECTIVE_START(timing, name, comm)     do { \
  timing##_timer.t1 = std::chrono::steady_clock::now(); \
  MPI_Barrier(comm); \
  timing##_timer.t2 = std::chrono::steady_clock::now(); \
  timing##_timer.time_span = (std::chrono::duration_cast<std::chrono::duration<double> >(timing##_timer.t2 - timing##_timer.t1)); \
  timing##_timer.names.push_back("bar_" name); \
  timing##_timer.durations.push_back(timing##_timer.time_span.count()); \
  timing##_timer.counts.push_back(0); \
  timing##_timer.t1 = std::chrono::steady_clock::now(); } while (0)
#define TIMER_END(timing, name, n_elem) do { timing##_timer.t2 = std::chrono::steady_clock::now(); \
                                             timing##_timer.time_span = (std::chrono::duration_cast<std::chrono::duration<double> >(timing##_timer.t2 - timing##_timer.t1)); \
                                             timing##_timer.names.push_back(name); \
                                             timing##_timer.durations.push_back(timing##_timer.time_span.count()); \
                                             timing##_timer.counts.push_back(n_elem); \
                                        } while (0)


#define TIMER_REPORT(timing, rank) \
        do { \
          std::stringstream output; \
          output << std::fixed; \
          output << "R " << rank << " " << #timing << "\theader\t[,"; \
          std::ostream_iterator<std::string> nit(output, ","); \
          std::copy(timing##_timer.names.begin(), timing##_timer.names.end(), nit); \
          output << "]\n" << #timing << "\tdur\t[,"; \
          output.precision(9); \
          std::ostream_iterator<double> dit(output, ","); \
          std::copy(timing##_timer.durations.begin(), timing##_timer.durations.end(), dit); \
          output << "]\n" << #timing << "\tcount\t[,"; \
          output.precision(0); \
          std::ostream_iterator<double> cit(output, ","); \
          std::copy(timing##_timer.counts.begin(), timing##_timer.counts.end(), cit); \
          output << "]"; \
          fflush(stdout); \
          printf("[BENCH] %s\n", output.str().c_str()); \
          fflush(stdout); \
        } while (0)

// do not use this function for now.  mxx::min_element and max_element has invalid read problem, reported by valgrind.
#define TIMER_REPORT_MPI_LOC(timing, rank, comm) \
        do { \
          int p; \
          MPI_Comm_size(comm, &p); \
          \
          printf("durations size = %lu\n", timing##_timer.durations.size());\
          auto dur_mins = mxx::min_element(timing##_timer.durations, comm); \
          auto dur_maxs = mxx::max_element(timing##_timer.durations, comm); \
          auto dur_means = ::mxx::reduce(timing##_timer.durations, 0, ::std::plus<double>(),comm); \
          ::std::for_each(timing##_timer.durations.begin(), timing##_timer.durations.end(), [](double &x) { x = x*x; }); \
          auto dur_stdevs = ::mxx::reduce(timing##_timer.durations, 0, ::std::plus<double>(), comm); \
          \
          printf("count size = %lu\n", timing##_timer.durations.size());\
          auto cnt_mins = ::mxx::min_element(timing##_timer.counts, comm); \
          auto cnt_maxs = ::mxx::max_element(timing##_timer.counts, comm); \
          auto cnt_means = ::mxx::reduce(timing##_timer.counts, 0, ::std::plus<double>(), comm); \
          ::std::for_each(timing##_timer.counts.begin(), timing##_timer.counts.end(), [](double &x) { x = x*x; }); \
          auto cnt_stdevs = ::mxx::reduce(timing##_timer.counts, 0, ::std::plus<double>(), comm); \
          \
          if (rank == 0) { \
            ::std::for_each(dur_means.begin(), dur_means.end(), [p](double & x) { x /= p; }); \
            ::std::transform(dur_stdevs.begin(), dur_stdevs.end(), dur_means.begin(), dur_stdevs.begin(), \
                             [p](double const & x, double const & y) { return ::std::sqrt(x / p - y * y); }); \
            ::std::for_each(cnt_means.begin(), cnt_means.end(), [p](double & x) { x /= p; }); \
            ::std::transform(cnt_stdevs.begin(), cnt_stdevs.end(), cnt_means.begin(), cnt_stdevs.begin(), \
                             [p](double const & x, double const & y) { return ::std::sqrt(x / p - y * y); }); \
            \
            std::stringstream output; \
            std::ostream_iterator<double> dit(output, ","); \
            std::ostream_iterator<int> iit(output, ","); \
            output << std::fixed; \
            output << "R " << rank << "/" << p << " " << #timing << "\theader\t[,"; \
            std::ostream_iterator<std::string> nit(output, ","); \
            std::copy(timing##_timer.names.begin(), timing##_timer.names.end(), nit); \
            output.precision(9); \
            auto get_first = [](const std::pair<double, int>& x){return x.first;};\
            auto get_second = [](const std::pair<double, int>& x){return x.second;};\
            output << "]\n" << #timing << "\tdur_min\t[,"; \
            std::transform(dur_mins.begin(), dur_mins.end(), dit, get_first); \
            output << "]\n" << #timing << "\tdur_min_idx\t[,"; \
            std::transform(dur_mins.begin(), dur_mins.end(), iit, get_second); \
            output << "]\n" << #timing << "\tdur_max\t[,"; \
            std::transform(dur_maxs.begin(), dur_maxs.end(), dit, get_first); \
            output << "]\n" << #timing << "\tdur_max_idx\t[,"; \
            std::transform(dur_maxs.begin(), dur_maxs.end(), iit, get_second); \
            output << "]\n" << #timing << "\tdur_mean\t[,"; \
            std::copy(dur_means.begin(), dur_means.end(), dit); \
            output << "]\n" << #timing << "\tdur_stdev\t[,"; \
            std::copy(dur_stdevs.begin(), dur_stdevs.end(), dit); \
            output.precision(0); \
            output << "]\n" << #timing << "\tcnt_min\t[,"; \
            std::transform(cnt_mins.begin(), cnt_mins.end(), dit, get_first); \
            output << "]\n" << #timing << "\tcnt_min_idx\t[,"; \
            std::transform(cnt_mins.begin(), cnt_mins.end(), iit, get_second); \
            output << "]\n" << #timing << "\tcnt_max\t[,"; \
            std::transform(cnt_maxs.begin(), cnt_maxs.end(), dit, get_first); \
            output << "]\n" << #timing << "\tcnt_max_idx\t[,"; \
            std::transform(cnt_maxs.begin(), cnt_maxs.end(), iit, get_second); \
            output.precision(2); \
            output << "]\n" << #timing << "\tcnt_mean\t[,"; \
            std::copy(cnt_means.begin(), cnt_means.end(), dit); \
            output << "]\n" << #timing << "\tcnt_stdev\t[,"; \
            std::copy(cnt_stdevs.begin(), cnt_stdevs.end(), dit); \
            output << "]"; \
            fflush(stdout); \
            printf("[BENCH] %s\n", output.str().c_str()); \
            fflush(stdout); \
          } \
          MPI_Barrier(comm); \
        } while (0)

#define TIMER_REPORT_MPI(timing, rank, comm) \
        do { \
          int p; \
          MPI_Comm_size(comm, &p); \
          \
          auto dur_mins = mxx::reduce(timing##_timer.durations, 0, [](double const & x, double const & y) { return ::std::min(x, y); }, comm); \
          auto dur_maxs = mxx::reduce(timing##_timer.durations, 0, [](double const & x, double const & y) { return ::std::max(x, y); }, comm); \
          auto dur_means = ::mxx::reduce(timing##_timer.durations, 0, ::std::plus<double>(),comm); \
          ::std::for_each(timing##_timer.durations.begin(), timing##_timer.durations.end(), [](double &x) { x = x*x; }); \
          auto dur_stdevs = ::mxx::reduce(timing##_timer.durations, 0, ::std::plus<double>(), comm); \
          \
          printf("count size = %lu\n", timing##_timer.durations.size());\
          auto cnt_mins = ::mxx::reduce(timing##_timer.counts, 0, [](double const & x, double const & y) { return ::std::min(x, y); }, comm); \
          auto cnt_maxs = ::mxx::reduce(timing##_timer.counts, 0, [](double const & x, double const & y) { return ::std::max(x, y); }, comm); \
          auto cnt_means = ::mxx::reduce(timing##_timer.counts, 0, ::std::plus<double>(), comm); \
          ::std::for_each(timing##_timer.counts.begin(), timing##_timer.counts.end(), [](double &x) { x = x*x; }); \
          auto cnt_stdevs = ::mxx::reduce(timing##_timer.counts, 0, ::std::plus<double>(), comm); \
          \
          if (rank == 0) { \
            ::std::for_each(dur_means.begin(), dur_means.end(), [p](double & x) { x /= p; }); \
            ::std::transform(dur_stdevs.begin(), dur_stdevs.end(), dur_means.begin(), dur_stdevs.begin(), \
                             [p](double const & x, double const & y) { return ::std::sqrt(x / p - y * y); }); \
            ::std::for_each(cnt_means.begin(), cnt_means.end(), [p](double & x) { x /= p; }); \
            ::std::transform(cnt_stdevs.begin(), cnt_stdevs.end(), cnt_means.begin(), cnt_stdevs.begin(), \
                             [p](double const & x, double const & y) { return ::std::sqrt(x / p - y * y); }); \
            \
            std::stringstream output; \
            std::ostream_iterator<double> dit(output, ","); \
            std::ostream_iterator<int> iit(output, ","); \
            output << std::fixed; \
            output << "R " << rank << "/" << p << " " << #timing << "\theader\t[,"; \
            std::ostream_iterator<std::string> nit(output, ","); \
            std::copy(timing##_timer.names.begin(), timing##_timer.names.end(), nit); \
            output.precision(9); \
            output << "]\n" << #timing << "\tdur_min\t[,"; \
            std::copy(dur_mins.begin(), dur_mins.end(), dit); \
            output << "]\n" << #timing << "\tdur_max\t[,"; \
            std::copy(dur_maxs.begin(), dur_maxs.end(), dit); \
            output << "]\n" << #timing << "\tdur_mean\t[,"; \
            std::copy(dur_means.begin(), dur_means.end(), dit); \
            output << "]\n" << #timing << "\tdur_stdev\t[,"; \
            std::copy(dur_stdevs.begin(), dur_stdevs.end(), dit); \
            output.precision(0); \
            output << "]\n" << #timing << "\tcnt_min\t[,"; \
            std::copy(cnt_mins.begin(), cnt_mins.end(), dit); \
            output << "]\n" << #timing << "\tcnt_max\t[,"; \
            std::copy(cnt_maxs.begin(), cnt_maxs.end(), dit); \
            output.precision(2); \
            output << "]\n" << #timing << "\tcnt_mean\t[,"; \
            std::copy(cnt_means.begin(), cnt_means.end(), dit); \
            output << "]\n" << #timing << "\tcnt_stdev\t[,"; \
            std::copy(cnt_stdevs.begin(), cnt_stdevs.end(), dit); \
            output << "]"; \
            fflush(stdout); \
            printf("[BENCH] %s\n", output.str().c_str()); \
            fflush(stdout); \
          } \
          MPI_Barrier(comm); \
        } while (0)


#else

#define TIMER_INIT(timing)
#define TIMER_RESET(timing)
#define TIMER_LOOP_START(timing)
#define TIMER_LOOP_RESUME(timing)
#define TIMER_LOOP_PAUSE(timing)
#define TIMER_LOOP_END(timing, name, n_elem)
#define TIMER_START(timing)
#define TIMER_COLLECTIVE_START(timing, name, comm)
#define TIMER_END(timing, name, n_elem)
#define TIMER_REPORT(timing, rank)
#define TIMER_REPORT_MPI(timing, rank, comm)

#endif


#endif /* SRC_UTILS_TIMER_HPP_ */
