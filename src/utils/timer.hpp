/**
 * @file    timer.hpp
 * @ingroup
 * @author  tpan
 * @brief
 * @details
 *
 * Copyright (c) 2015 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
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



#define TIMER_INIT(timing)      std::chrono::steady_clock::time_point timing##_t1, timing##_t2; \
                                 std::vector<std::string> timing##_names; \
                                 std::vector<double> timing##_durations; \
                                 std::vector<double> timing##_counts; \
                                 std::chrono::duration<double> timing##_time_span;


#define TIMER_RESET(timing)     do {  timing##_names.clear(); \
                                       timing##_durations.clear(); \
                                       timing##_count.clear(); \
                                 } while (0)

#define TIMER_LOOP(timing)     do { timing##_time_span = std::chrono::duration<double>::zero(); } while (0)
#define TIMER_LOOP_RESUME(timing)    do { timing##_t1 = std::chrono::steady_clock::now(); } while (0)
#define TIMER_LOOP_PAUSE(timing)     do { timing##_t2 = std::chrono::steady_clock::now(); \
                                      timing##_time_span += (std::chrono::duration_cast<std::chrono::duration<double> >(timing##_t2 - timing##_t1)); } while (0)
#define TIMER_LOOP_END(timing, name, n_elem) do { timing##_names.push_back(name); \
                                             timing##_durations.push_back(timing##_time_span.count()); \
                                             timing##_counts.push_back(n_elem); \
                                        } while (0)

#define TIMER_START(timing)     do { timing##_t1 = std::chrono::steady_clock::now(); } while (0)
#define TIMER_COLLECTIVE_START(timing, name, comm)     do { \
  timing##_t1 = std::chrono::steady_clock::now(); \
  MPI_Barrier(comm); \
  timing##_t2 = std::chrono::steady_clock::now(); \
  timing##_time_span = (std::chrono::duration_cast<std::chrono::duration<double> >(timing##_t2 - timing##_t1)); \
  timing##_names.push_back("bar_" name); \
  timing##_durations.push_back(timing##_time_span.count()); \
  timing##_counts.push_back(0); \
  timing##_t1 = std::chrono::steady_clock::now(); } while (0)
#define TIMER_END(timing, name, n_elem) do { timing##_t2 = std::chrono::steady_clock::now(); \
                                             timing##_time_span = (std::chrono::duration_cast<std::chrono::duration<double> >(timing##_t2 - timing##_t1)); \
                                             timing##_names.push_back(name); \
                                             timing##_durations.push_back(timing##_time_span.count()); \
                                             timing##_counts.push_back(n_elem); \
                                        } while (0)


#define TIMER_REPORT(timing, rank) \
        do { \
          std::stringstream output; \
          output << std::fixed; \
          output << "R " << rank << " " << #timing << " header\t["; \
          std::ostream_iterator<std::string> nit(output, ","); \
          std::copy(timing##_names.begin(), timing##_names.end(), nit); \
          output << "]\n" << #timing << "\tdur\t\t["; \
          output.precision(9); \
          std::ostream_iterator<double> dit(output, ","); \
          std::copy(timing##_durations.begin(), timing##_durations.end(), dit); \
          output << "]\n" << #timing << "\tcount\t\t["; \
          output.precision(0); \
          std::ostream_iterator<double> cit(output, ","); \
          std::copy(timing##_counts.begin(), timing##_counts.end(), cit); \
          output << "]"; \
          fflush(stdout); \
          printf("[BENCH]\t%s\n", output.str().c_str()); \
          fflush(stdout); \
        } while (0)

// not used.
//template <typename T>
//::std::vector<::std::tuple<T, T, T, T> > compute_elem_stats(::std::vector<T> const & v, MPI_Comm comm) {
//  int p;
//  MPI_Comm_size(comm, &p);
//
//  ::std::vector<::std::tuple<T, T, T, T> > u(v.size());
//  ::std::transform(v.begin(), v.end(), u.begin(), [](T const & x) { return ::std::make_tuple(x, x, x, x*x); });
//  auto stats = ::mxx2::reduce(u, [](::std::tuple<T, T, T, T> const &x, ::std::tuple<T, T, T, T> const &y) {
//    return ::std::make_tuple(::std::min(::std::get<0>(x), ::std::get<0>(y)),
//                             ::std::max(::std::get<1>(x), ::std::get<1>(y)),
//                             ::std::plus<T>(::std::get<2>(x), ::std::get<2>(y)),
//                             ::std::plus<T>(::std::get<3>(x), ::std::get<3>(y)));
//    }, comm, 0);
//  ::std::transform(stats.begin(), stats.end(), stats.begin(), [p](::std::tuple<T, T, T, T> x) {
//    ::std::get<2>(x) /= p;
//    ::std::get<3>(x) = ::std::sqrt(::std::get<3>(x) / p - ::std::get<2>(x) * ::std::get<2>(x));
//    return x;
//  });
//
//  return stats;
//}

#define TIMER_REPORT_MPI(timing, rank, comm) \
        do { \
          int p; \
          MPI_Comm_size(comm, &p); \
          \
          auto dur_mins = ::mxx2::reduce_loc(timing##_durations, MPI_MINLOC, comm, 0); \
          auto dur_maxs = ::mxx2::reduce_loc(timing##_durations, MPI_MAXLOC, comm, 0); \
          auto dur_means = ::mxx2::reduce(timing##_durations, ::std::plus<double>(), comm, 0); \
          ::std::for_each(timing##_durations.begin(), timing##_durations.end(), [](double &x) { x = x*x; }); \
          auto dur_stdevs = ::mxx2::reduce(timing##_durations, ::std::plus<double>(), comm, 0); \
          \
          auto cnt_mins = ::mxx2::reduce_loc(timing##_counts, MPI_MINLOC, comm, 0); \
          auto cnt_maxs = ::mxx2::reduce_loc(timing##_counts, MPI_MAXLOC, comm, 0); \
          auto cnt_means = ::mxx2::reduce(timing##_counts, ::std::plus<double>(), comm, 0); \
          ::std::for_each(timing##_counts.begin(), timing##_counts.end(), [](double &x) { x = x*x; }); \
          auto cnt_stdevs = ::mxx2::reduce(timing##_counts, ::std::plus<double>(), comm, 0); \
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
            output << "R " << rank << "/" << p << " " << #timing << " header\t["; \
            std::ostream_iterator<std::string> nit(output, ","); \
            std::copy(timing##_names.begin(), timing##_names.end(), nit); \
            output.precision(9); \
            output << "]\n" << #timing << "\tdur_min\t\t["; \
            std::copy(dur_mins.first.begin(), dur_mins.first.end(), dit); \
            output << "]\n" << #timing << "\tdur_min_idx\t\t["; \
            std::copy(dur_mins.second.begin(), dur_mins.second.end(), iit); \
            output << "]\n" << #timing << "\tdur_max\t\t["; \
            std::copy(dur_maxs.first.begin(), dur_maxs.first.end(), dit); \
            output << "]\n" << #timing << "\tdur_max_idx\t\t["; \
            std::copy(dur_maxs.second.begin(), dur_maxs.second.end(), iit); \
            output << "]\n" << #timing << "\tdur_mean\t["; \
            std::copy(dur_means.begin(), dur_means.end(), dit); \
            output << "]\n" << #timing << "\tdur_stdev\t["; \
            std::copy(dur_stdevs.begin(), dur_stdevs.end(), dit); \
            output.precision(0); \
            output << "]\n" << #timing << "\tcnt_min\t\t["; \
            std::copy(cnt_mins.first.begin(), cnt_mins.first.end(), dit); \
            output << "]\n" << #timing << "\tcnt_min_idx\t\t["; \
            std::copy(cnt_mins.second.begin(), cnt_mins.second.end(), iit); \
            output << "]\n" << #timing << "\tcnt_max\t\t["; \
            std::copy(cnt_maxs.first.begin(), cnt_maxs.first.end(), dit); \
            output << "]\n" << #timing << "\tcnt_max_idx\t\t["; \
            std::copy(cnt_maxs.second.begin(), cnt_maxs.second.end(), iit); \
            output.precision(2); \
            output << "]\n" << #timing << "\tcnt_mean\t["; \
            std::copy(cnt_means.begin(), cnt_means.end(), dit); \
            output << "]\n" << #timing << "\tcnt_stdev\t["; \
            std::copy(cnt_stdevs.begin(), cnt_stdevs.end(), dit); \
            output << "]"; \
            fflush(stdout); \
            printf("[BENCH]\t%s\n", output.str().c_str()); \
            fflush(stdout); \
          } \
          MPI_Barrier(comm); \
        } while (0)


#else

#define TIMER_INIT(timing)
#define TIMER_RESET(timing)
#define TIMER_LOOP(timing)
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
