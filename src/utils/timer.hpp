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

#define COLLECT_TIMES 1

#if defined(COLLECT_TIMES)

#include <chrono>   // clock
#include <cstdio>   // printf
#include <vector>
#include <string>

#define TIMER_INIT(session)      std::chrono::steady_clock::time_point session##_t1, session##_t2; \
                                 std::chrono::duration<double> session##_time_span; \
                                 std::vector<std::string> session##_dur_names; \
                                 std::vector<double> session##_durations;

#define TIMER_START(session)     session##_t1 = std::chrono::steady_clock::now()
#define TIMER_END(session, name) session##_t2 = std::chrono::steady_clock::now(); \
                                 session##_time_span = std::chrono::duration_cast<std::chrono::duration<double>>(session##_t2 - session##_t1); \
                                 session##_dur_names.push_back(name); \
                                 session##_durations.push_back(session##_time_span.count());

#define TIMER_REPORT(session) \
        { \
          std::stringstream names; \
          std::ostream_iterator<std::string> nit(names, ", "); \
          std::copy(session##_names.begin(), session##_names.end(), nit); \
          std::stringstream durations; \
          std::ostream_iterator<double> dit(durations, ", "); \
          std::copy(session##_durations.begin(), session##_durations.end(), dit); \
          printf("%s\n\t[\s]\n\t[\s]\n", #session, names.str().c_str(), durations.std().c_str()); \
        };
#else

#define TIMER_INIT(session)
#define TIMER_START(session)
#define TIMER_END(session, name)
#define TIMER_REPORT(session)

#endif


#endif /* SRC_UTILS_TIMER_HPP_ */
