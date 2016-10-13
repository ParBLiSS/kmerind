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

#if BL_BENCHMARK_TIME == 1

#include <chrono>   // clock
//#include <cstdio>   // printf
#include <vector>
#include <string>
#include <algorithm>  // std::min
#include <sstream>
#include <cmath>

#include <mxx/reduction.hpp>


class Timer {
  protected:
    std::chrono::steady_clock::time_point t1, t2;
    std::vector<std::string> names;
    std::vector<double> durations; 
    std::vector<double> counts;
    std::chrono::duration<double> time_span;

  public:
    void reset() {
      names.clear();
      durations.clear();
      counts.clear();
    }


//=========== loop stuff.
    void loop_start() { time_span = std::chrono::duration<double>::zero(); }
    void loop_resume() { t1 = std::chrono::steady_clock::now(); }
    void loop_pause() {
    	t2 = std::chrono::steady_clock::now();
    	time_span += (std::chrono::duration_cast<std::chrono::duration<double> >(t2 - t1));
    }
    void loop_end(::std::string const & name, double const & n_elem) {
    	names.push_back(name);
    	durations.push_back(time_span.count());
		counts.push_back(n_elem);
    }

//============ timer start
    void start() { t1 = std::chrono::steady_clock::now(); }
    void collective_start(::std::string const & name, ::mxx::comm const & comm) {

      // time a barrier.
      t1 = std::chrono::steady_clock::now();
      comm.barrier();
      t2 = std::chrono::steady_clock::now();
      time_span = (std::chrono::duration_cast<std::chrono::duration<double> >(t2 - t1));

      ::std::string tmp("barrier_"); tmp.append(name);
      names.push_back(tmp);
      durations.push_back(time_span.count());
      counts.push_back(0);

      t1 = std::chrono::steady_clock::now();
    }
    void end(::std::string const & name, double const & n_elem) {
      t2 = std::chrono::steady_clock::now();
      time_span = (std::chrono::duration_cast<std::chrono::duration<double> >(t2 - t1));

      names.push_back(name);
      durations.push_back(time_span.count());
      counts.push_back(n_elem);
    }
    void collective_end(::std::string const & name, double const & n_elem, ::mxx::comm const & comm) {

		// time a barrier.
		comm.barrier();

		end(name, n_elem);
    }
    void report(::std::string const & title) {
        std::stringstream output;

        std::ostream_iterator<std::string> nit(output, ",");
        std::ostream_iterator<double> dit(output, ",");

        output << std::fixed;
        output << "[TIME] " << title << "\theader (s)\t[,";
        std::copy(names.begin(), names.end(), nit);
        output << "]" << std::endl;

        output.precision(9);

        output << "[TIME] " << title << "\tdur\t[,";
        std::copy(durations.begin(), durations.end(), dit);
        output << "]" << ::std::endl;

        output << "[TIME] " << title << "\tcum\t[,";
        std::partial_sum(durations.begin(), durations.end(), dit);
        output << "]" << ::std::endl;

        output.precision(0);
        output << "[TIME] " << title << "\tcount\t[,";
        std::copy(counts.begin(), counts.end(), dit);
        output << "]";

        // print pending stuff, then print entire string at once (minimizes multiple threads/processes mixing output )
        fflush(stdout);
        printf("%s\n", output.str().c_str());
        fflush(stdout);
    }

#if 0
    // do not use this function for now.  mxx::min_element and max_element has invalid
    //  read problem, reported by valgrind.
    void report_loc(::std::string const & title, ::mxx::comm const & comm) {

        auto dur_mins = mxx::min_element(durations, comm);
        auto dur_maxs = mxx::max_element(durations, comm);
        auto dur_means = ::mxx::reduce(durations, 0, ::std::plus<double>(),comm);
        ::std::for_each(durations.begin(), durations.end(), [](double &x) { x = x*x; });
        auto dur_stdevs = ::mxx::reduce(durations, 0, ::std::plus<double>(), comm);

        auto cnt_mins = ::mxx::min_element(counts, comm);
        auto cnt_maxs = ::mxx::max_element(counts, comm);
        auto cnt_means = ::mxx::reduce(counts, 0, ::std::plus<double>(), comm);
        ::std::for_each(counts.begin(), counts.end(), [](double &x) { x = x*x; });
        auto cnt_stdevs = ::mxx::reduce(counts, 0, ::std::plus<double>(), comm);

    	int rank = comm.rank();

        if (rank == 0) {
        	int p = comm.size();
            auto get_first = [](const std::pair<double, int>& x){return x.first;};
            auto get_second = [](const std::pair<double, int>& x){return x.second;};


          ::std::for_each(dur_means.begin(), dur_means.end(), [&p](double & x) { x /= p; });
          ::std::transform(dur_stdevs.begin(), dur_stdevs.end(), dur_means.begin(), dur_stdevs.begin(),
                           [&p](double const & x, double const & y) { return ::std::sqrt(x / p - y * y); });

          ::std::for_each(cnt_means.begin(), cnt_means.end(), [&p](double & x) { x /= p; });
          ::std::transform(cnt_stdevs.begin(), cnt_stdevs.end(), cnt_means.begin(), cnt_stdevs.begin(),
                           [&p](double const & x, double const & y) { return ::std::sqrt(x / p - y * y); });

          std::stringstream output;
          std::ostream_iterator<std::string> nit(output, ",");
          std::ostream_iterator<double> dit(output, ",");
          std::ostream_iterator<int> iit(output, ",");

          output << std::fixed;
          output << "[TIME] " << "R " << rank << "/" << p << " ";

          output << "[TIME] " << title << "\theader (s)\t[,";
          std::copy(names.begin(), names.end(), nit);
          output << "]" << std::endl;

          output.precision(9);
          output << "[TIME] " << title << "\tdur_min\t[,";
          std::transform(dur_mins.begin(), dur_mins.end(), dit, get_first);
          output << "]" << std::endl;
          output << "[TIME] " << title << "\tdur_min_idx\t[,";
          std::transform(dur_mins.begin(), dur_mins.end(), iit, get_second);
          output << "]" << std::endl;

          output << "[TIME] " << title << "\tdur_max\t[,";
          std::transform(dur_maxs.begin(), dur_maxs.end(), dit, get_first);
          output << "]" << std::endl;
          output << "[TIME] " << title << "\tdur_max_idx\t[,";
          std::transform(dur_maxs.begin(), dur_maxs.end(), iit, get_second);
          output << "]" << std::endl;

          output << "[TIME] " << title << "\tdur_mean\t[,";
          std::copy(dur_means.begin(), dur_means.end(), dit);
          output << "]" << std::endl;

          output << "[TIME] " << title << "\tdur_stdev\t[,";
          std::copy(dur_stdevs.begin(), dur_stdevs.end(), dit);
          output << "]" << std::endl;

          output.precision(2);
          output << "[TIME] " << title << "\tcnt_min\t[,";
          std::transform(cnt_mins.begin(), cnt_mins.end(), dit, get_first);
          output << "]" << std::endl;
          output << "[TIME] " << title << "\tcnt_min_idx\t[,";
          std::transform(cnt_mins.begin(), cnt_mins.end(), iit, get_second);
          output << "]" << std::endl;

          output << "[TIME] " << title << "\tcnt_max\t[,";
          std::transform(cnt_maxs.begin(), cnt_maxs.end(), dit, get_first);
          output << "]" << std::endl;
          output << "[TIME] " << title << "\tcnt_max_idx\t[,";
          std::transform(cnt_maxs.begin(), cnt_maxs.end(), dit, get_second);
          output << "]" << std::endl;

          output.precision(2);
          output << "[TIME] " << title << "\tcnt_mean\t[,";
          std::copy(cnt_means.begin(), cnt_means.end(), dit);
          output << "]" << std::endl;

          output << "[TIME] " << title << "\tcnt_stdev\t[,";
          std::copy(cnt_stdevs.begin(), cnt_stdevs.end(), dit);
          output << "]";


          fflush(stdout);
          printf("%s\n", output.str().c_str());
          fflush(stdout);
        }
        comm.barrier();
    }
#endif

    void report(::std::string const & title, ::mxx::comm const & comm) {
      std::vector<double> cumulative(durations);
      std::vector<double> dur_mins, dur_maxs, dur_means, dur_stdevs;
      std::vector<double> cum_mins, cum_maxs, cum_means, cum_stdevs;
      std::vector<double> cnt_mins, cnt_maxs, cnt_means, cnt_stdevs;
      int p = comm.size();
      int rank = comm.rank();

      if (durations.size() > 0) {

    	  dur_mins = ::mxx::reduce(durations, 0,
    			[](double const & x, double const & y) { return ::std::min(x, y); }, comm);
        dur_maxs = ::mxx::reduce(durations, 0,
        		[](double const & x, double const & y) { return ::std::max(x, y); }, comm);
        dur_means = ::mxx::reduce(durations, 0, ::std::plus<double>(), comm);
        ::std::for_each(durations.begin(), durations.end(), [](double &x) { x = x*x; });
        dur_stdevs = ::mxx::reduce(durations, 0, ::std::plus<double>(), comm);

        // compute partial sum of times.
        ::std::partial_sum(cumulative.begin(), cumulative.end(), cumulative.begin());
        cum_mins = ::mxx::reduce(cumulative, 0,
                                      [](double const & x, double const & y) { return ::std::min(x, y); }, comm);
        cum_maxs = ::mxx::reduce(cumulative, 0,
                                      [](double const & x, double const & y) { return ::std::max(x, y); }, comm);
        cum_means = ::mxx::reduce(cumulative, 0, ::std::plus<double>(), comm);
        ::std::for_each(cumulative.begin(), cumulative.end(), [](double &x) { x = x*x; });
        cum_stdevs = ::mxx::reduce(cumulative, 0, ::std::plus<double>(), comm);


        cnt_mins = ::mxx::reduce(counts, 0,
        		[](double const & x, double const & y) { return ::std::min(x, y); }, comm);
        cnt_maxs = ::mxx::reduce(counts, 0,
        		[](double const & x, double const & y) { return ::std::max(x, y); }, comm);
        cnt_means = ::mxx::reduce(counts, 0, ::std::plus<double>(), comm);
        ::std::for_each(counts.begin(), counts.end(), [](double &x) { x = x*x; });
        cnt_stdevs = ::mxx::reduce(counts, 0, ::std::plus<double>(), comm);

        if (rank == 0) {

          ::std::for_each(dur_means.begin(), dur_means.end(), [&p](double & x) { x /= p; });
          ::std::transform(dur_stdevs.begin(), dur_stdevs.end(), dur_means.begin(), dur_stdevs.begin(),
                           [&p](double const & x, double const & y) { return ::std::sqrt(x / p - y * y); });

          ::std::for_each(cum_means.begin(), cum_means.end(), [&p](double & x) { x /= p; });
          ::std::transform(cum_stdevs.begin(), cum_stdevs.end(), cum_means.begin(), cum_stdevs.begin(),
                           [&p](double const & x, double const & y) { return ::std::sqrt(x / p - y * y); });

          ::std::for_each(cnt_means.begin(), cnt_means.end(), [&p](double & x) { x /= p; });
          ::std::transform(cnt_stdevs.begin(), cnt_stdevs.end(), cnt_means.begin(), cnt_stdevs.begin(),
                           [&p](double const & x, double const & y) { return ::std::sqrt(x / p - y * y); });
        }


      }


    	if (rank == 0) {
          std::stringstream output;
          std::ostream_iterator<std::string> nit(output, ",");
          std::ostream_iterator<double> dit(output, ",");

          output << std::fixed;
          output << "[TIME] " << "R " << rank << "/" << p << std::endl;

          output << "[TIME] " << title << "\theader (s)\t[,";
          std::copy(names.begin(), names.end(), nit);
          output << "]" << std::endl;

          output.precision(9);
          output << "[TIME] " << title << "\tdur_min\t[,";
          std::copy(dur_mins.begin(), dur_mins.end(), dit);
          output << "]" << std::endl;

          output << "[TIME] " << title << "\tdur_max\t[,";
          std::copy(dur_maxs.begin(), dur_maxs.end(), dit);
          output << "]" << std::endl;

          output << "[TIME] " << title << "\tdur_mean\t[,";
          std::copy(dur_means.begin(), dur_means.end(), dit);
          output << "]" << std::endl;

          output << "[TIME] " << title << "\tdur_stdev\t[,";
          std::copy(dur_stdevs.begin(), dur_stdevs.end(), dit);
          output << "]" << std::endl;


          output.precision(9);
          output << "[TIME] " << title << "\tcum_min\t[,";
          std::copy(cum_mins.begin(), cum_mins.end(), dit);
          output << "]" << std::endl;

          output << "[TIME] " << title << "\tcum_max\t[,";
          std::copy(cum_maxs.begin(), cum_maxs.end(), dit);
          output << "]" << std::endl;

          output << "[TIME] " << title << "\tcum_mean\t[,";
          std::copy(cum_means.begin(), cum_means.end(), dit);
          output << "]" << std::endl;

          output << "[TIME] " << title << "\tcum_stdev\t[,";
          std::copy(cum_stdevs.begin(), cum_stdevs.end(), dit);
          output << "]" << std::endl;


          output.precision(2);
          output << "[TIME] " << title << "\tcnt_min\t[,";
          std::copy(cnt_mins.begin(), cnt_mins.end(), dit);
          output << "]" << std::endl;

          output << "[TIME] " << title << "\tcnt_max\t[,";
          std::copy(cnt_maxs.begin(), cnt_maxs.end(), dit);
          output << "]" << std::endl;

          output.precision(2);
          output << "[TIME] " << title << "\tcnt_mean\t[,";
          std::copy(cnt_means.begin(), cnt_means.end(), dit);
          output << "]" << std::endl;

          output << "[TIME] " << title << "\tcnt_stdev\t[,";
          std::copy(cnt_stdevs.begin(), cnt_stdevs.end(), dit);
          output << "]";

          fflush(stdout);
          printf("%s\n", output.str().c_str());
          fflush(stdout);
        }
        comm.barrier();
    }

};


#define BL_TIMER_INIT(title)      Timer title##_timer;
#define BL_TIMER_RESET(title)     do {  title##_timer.reset(); } while (0)

#define BL_TIMER_LOOP_START(title)     do { title##_timer.loop_start(); } while (0)
#define BL_TIMER_LOOP_RESUME(title)    do { title##_timer.loop_resume(); } while (0)
#define BL_TIMER_LOOP_PAUSE(title)     do { title##_timer.loop_pause(); } while (0)
#define BL_TIMER_LOOP_END(title, name, n_elem) do { title##_timer.loop_end(name, n_elem); } while (0)

#define BL_TIMER_START(title)     do { title##_timer.start(); } while (0)
#define BL_TIMER_END(title, name, n_elem) do { title##_timer.end(name, n_elem); } while (0)
#define BL_TIMER_COLLECTIVE_START(title, name, comm) do { title##_timer.collective_start(name, comm); } while (0)
#define BL_TIMER_COLLECTIVE_END(title, name, n_elem, comm) do { title##_timer.collective_end(name, n_elem, comm); } while (0)
#define BL_TIMER_REPORT(title) do { title##_timer.report(#title); } while (0)
#define BL_TIMER_REPORT_NAMED(title, name) do { title##_timer.report(name); } while (0)

#if 0
// do not use this function for now.  mxx::min_element and max_element has invalid read problem, reported by valgrind.
#define BL_TIMER_REPORT_MPI_LOC(title, comm) do { title##_timer.report_loc(#title, comm); } while (0)
#endif

#define BL_TIMER_REPORT_MPI(title, comm) do { title##_timer.report(#title, comm); } while (0)
#define BL_TIMER_REPORT_MPI_NAMED(title, name, comm) do { title##_timer.report(name, comm); } while (0)


#else

#define BL_TIMER_INIT(title)
#define BL_TIMER_RESET(title)
#define BL_TIMER_LOOP_START(title)
#define BL_TIMER_LOOP_RESUME(title)
#define BL_TIMER_LOOP_PAUSE(title)
#define BL_TIMER_LOOP_END(title, name, n_elem)
#define BL_TIMER_START(title)
#define BL_TIMER_END(title, name, n_elem)
#define BL_TIMER_COLLECTIVE_START(title, name, comm)
#define BL_TIMER_COLLECTIVE_end(title, name, n_elem, comm)
#define BL_TIMER_REPORT(title)
#define BL_TIMER_REPORT_MPI(title, comm)
#define BL_TIMER_REPORT_NAMED(title, name)
#define BL_TIMER_REPORT_MPI_NAMED(title, name, comm)

#endif


#endif /* SRC_UTILS_TIMER_HPP_ */
