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

#define BENCHMARK 1

#if defined(BENCHMARK)

#include <chrono>   // clock
#include <cstdio>   // printf
#include <vector>
#include <string>
#include <algorithm>  // std::min

#include "io/mxx_support.hpp"


#define TIMER_INIT(session)      std::chrono::steady_clock::time_point session##_t1, session##_t2; \
                                 std::vector<std::string> session##_names; \
                                 std::vector<double> session##_durations; \
                                 std::vector<double> session##_counts;

#define TIMER_RESET(session)     do {  session##_names.clear(); \
                                       session##_durations.clear(); \
                                       session##_count.clear(); \
                                 } while (0)

#define TIMER_START(session)     do { session##_t1 = std::chrono::steady_clock::now(); } while (0)

#define TIMER_END(session, name, n_elem) do { session##_t2 = std::chrono::steady_clock::now(); \
                                             std::chrono::duration<double> time_span = (std::chrono::duration_cast<std::chrono::duration<double> >(session##_t2 - session##_t1)); \
                                             session##_names.push_back(name); \
                                             session##_durations.push_back(time_span.count()); \
                                             session##_counts.push_back(n_elem); \
                                        } while (0)

#define TIMER_REPORT(session, rank) \
        do { \
          std::stringstream output; \
          output << std::fixed; \
          output << "R " << rank << " " << #session << " | ["; \
          std::ostream_iterator<std::string> nit(output, ","); \
          std::copy(session##_names.begin(), session##_names.end(), nit); \
          output << "] | ["; \
          output.precision(9); \
          std::ostream_iterator<double> dit(output, ","); \
          std::copy(session##_durations.begin(), session##_durations.end(), dit); \
          output << "] | ["; \
          output.precision(0); \
          std::ostream_iterator<double> cit(output, ","); \
          std::copy(session##_counts.begin(), session##_counts.end(), cit); \
          output << "]"; \
          fflush(stdout); \
          printf("%s\n", output.str().c_str()); \
          fflush(stdout); \
          MPI_Barrier(comm); \
        } while (0)


template <typename T>
::std::vector<::std::tuple<T, T, T, T> > compute_elem_stats(::std::vector<T> const & v, MPI_Comm comm) {
  int p;
  MPI_Comm_size(comm, &p);

  ::std::vector<::std::tuple<T, T, T, T> > u(v.size());
  ::std::transform(v.begin(), v.end(), u.begin(), [](T const & x) { return ::std::make_tuple(x, x, x, x*x); });
  auto stats = ::mxx2::reduce(u, [](::std::tuple<T, T, T, T> const &x, ::std::tuple<T, T, T, T> const &y) {
    return ::std::make_tuple(::std::min(::std::get<0>(x), ::std::get<0>(y)),
                             ::std::max(::std::get<1>(x), ::std::get<1>(y)),
                             ::std::plus<T>(::std::get<2>(x), ::std::get<2>(y)),
                             ::std::plus<T>(::std::get<3>(x), ::std::get<3>(y)));
    }, comm, 0);
  ::std::transform(stats.begin(), stats.end(), stats.begin(), [p](::std::tuple<T, T, T, T> x) {
    ::std::get<2>(x) /= p;
    ::std::get<3>(x) = ::std::sqrt(::std::get<3>(x) / p - ::std::get<2>(x) * ::std::get<2>(x));
    return x;
  });

  return stats;
}

#define TIMER_REPORT_MPI(session, rank, comm) \
        do { \
          int p; \
          MPI_Comm_size(comm, &p); \
          \
          auto dur_mins = ::mxx2::reduce(session##_durations, [](double const &x, double const &y) { return ::std::min(x,y); }, comm, 0); \
          auto dur_maxs = ::mxx2::reduce(session##_durations, [](double const &x, double const &y) { return ::std::max(x,y); }, comm, 0); \
          auto dur_means = ::mxx2::reduce(session##_durations, ::std::plus<double>(), comm, 0); \
          ::std::for_each(session##_durations.begin(), session##_durations.end(), [](double &x) { x = x*x; }); \
          auto dur_stdevs = ::mxx2::reduce(session##_durations, ::std::plus<double>(), comm, 0); \
          \
          auto cnt_mins = ::mxx2::reduce(session##_counts, [](double const &x, double const &y) { return ::std::min(x,y); }, comm, 0); \
          auto cnt_maxs = ::mxx2::reduce(session##_counts, [](double const &x, double const &y) { return ::std::max(x,y); }, comm, 0); \
          auto cnt_means = ::mxx2::reduce(session##_counts, ::std::plus<double>(), comm, 0); \
          ::std::for_each(session##_counts.begin(), session##_counts.end(), [](double &x) { x = x*x; }); \
          auto cnt_stdevs = ::mxx2::reduce(session##_counts, ::std::plus<double>(), comm, 0); \
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
            output << std::fixed; \
            output << "R " << rank << "/" << p << " " << #session << " header\t["; \
            std::ostream_iterator<std::string> nit(output, ","); \
            std::copy(session##_names.begin(), session##_names.end(), nit); \
            output.precision(9); \
            output << "]\n" << #session << "\tdur_min\t\t["; \
            std::ostream_iterator<double> dit(output, ","); \
            std::copy(dur_mins.begin(), dur_mins.end(), dit); \
            output << "]\n" << #session << "\tdur_max\t\t["; \
            std::copy(dur_maxs.begin(), dur_maxs.end(), dit); \
            output << "]\n" << #session << "\tdur_mean\t["; \
            std::copy(dur_means.begin(), dur_means.end(), dit); \
            output << "]\n" << #session << "\tdur_stdev\t["; \
            std::copy(dur_stdevs.begin(), dur_stdevs.end(), dit); \
            output.precision(0); \
            output << "]\n" << #session << "\tcnt_min\t\t["; \
            std::copy(cnt_mins.begin(), cnt_mins.end(), dit); \
            output << "]\n" << #session << "\tcnt_max\t\t["; \
            std::copy(cnt_maxs.begin(), cnt_maxs.end(), dit); \
            output.precision(2); \
            output << "]\n" << #session << "\tcnt_mean\t["; \
            std::copy(cnt_means.begin(), cnt_means.end(), dit); \
            output << "]\n" << #session << "\tcnt_stdev\t["; \
            std::copy(cnt_stdevs.begin(), cnt_stdevs.end(), dit); \
            output << "]"; \
            fflush(stdout); \
            printf("%s\n", output.str().c_str()); \
            fflush(stdout); \
          } \
          MPI_Barrier(comm); \
        } while (0)


#else

#define TIMER_INIT(session)
#define TIMER_RESET(session)
#define TIMER_START(session)
#define TIMER_END(session, name, count)
#define TIMER_REPORT(session, rank)
#define TIMER_REPORT_MPI(session, rank, comm)

#endif


#endif /* SRC_UTILS_TIMER_HPP_ */
