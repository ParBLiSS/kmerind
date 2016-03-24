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
 * @file    memory_usage.hpp
 * @ingroup
 * @author  tpan
 * @brief   functions to track memory usage during program execution (for marked functional blocks)
 * @details each "mark" call snapshots the current memory usage and peak memory usage.
 *          relies on http://nadeausoftware.com/articles/2012/07/c_c_tip_how_get_process_resident_set_size_physical_memory_use#GetProcessMemoryInfonbspforpeakandcurrentresidentsetsize
 *
 *          also see http://www.linuxatemyram.com/play.html and
 *          http://stackoverflow.com/questions/349889/how-do-you-determine-the-amount-of-linux-system-ram-in-c
 *
 */
#ifndef SRC_UTILS_MEMORY_USAGE_HPP_
#define SRC_UTILS_MEMORY_USAGE_HPP_

#include "bliss-logger_config.hpp"


#if BL_BENCHMARK_MEM == 1

#include <sys/resource.h>  // getrusage
#include <sys/types.h>     // meminfo
#include <sys/sysinfo.h>    // sysinfo

#include <chrono>   // clock
//#include <cstdio>   // printf
#include <unistd.h>  // getpid()
#include <vector>
#include <cstring>	// strerror
#include <string>
#include <algorithm>  // std::min
#include <sstream>

#include <io/io_exception.hpp>
#include <mxx/reduction.hpp>

//http://nadeausoftware.com/articles/2012/07/c_c_tip_how_get_process_resident_set_size_physical_memory_use#GetProcessMemoryInfonbspforpeakandcurrentresidentsetsize
// note:  reports in bytes.
#include "getRSS.h"

// also see http://www.linuxatemyram.com/play.html

//#include "stdlib.h"
//#include "stdio.h"
//#include "sys/types.h"
//#include "sys/sysinfo.h"
//
//// implement from http://linux.die.net/man/2/getrusage
//#include <sys/time.h>
//#include <sys/resource.h>
//
///// get current process's max allocated memory
//void get_mem_allocated(long long &phyMemUsed, long long &swapUsed) {
//  struct rusage u;
//  getrusage(RUSAGE_SELF, &u);
//  phyMemUsed = u.ru_maxrss;
//}
//
///// get system wide physical and swap usage
//void get_mem_allocated_sysinfo(long long &phyMemUsed, long long &swapUsed) {
//  //from http://stackoverflow.com/questions/63166/how-to-determine-cpu-and-memory-consumption-from-inside-a-process
//  struct sysinfo memInfo;
//  sysinfo (&memInfo);
//
//  phyMemUsed = (memInfo.totalram - memInfo.freeram) * memInfo.mem_unit;
//  swapUsed = (memInfo.totalswap - memInfo.freeswap) * memInfo.mem_unit;
//}
//
///// get changes in physical memory and swap usage relative to baseline (from get_mem_allocated_sysinfo)
//void get_mem_consumed_sysinfo(const long long &phyMemBaseline, const long long &swapBaseline, long long &phyMemUsed, long long &swapUsed) {
//  long long phy;
//  long long swap;
//  get_mem_allocated_sysinfo(phy, swap);
//
//  phyMemUsed = phy - phyMemBaseline;
//  swapUsed = swap - swapBaseline;
//}


///// get the physical mem size in bytes.
//static size_t get_physical_mem() {
//  //from http://stackoverflow.com/questions/63166/how-to-determine-cpu-and-memory-consumption-from-inside-a-process
//  struct sysinfo memInfo;
//
//  sysinfo (&memInfo);
//  printf("total %ld, buffer %ld, shared %ld, free %ld\n", memInfo.totalram, memInfo.bufferram, memInfo.sharedram, memInfo.freeram);
//
//  return memInfo.totalram * memInfo.mem_unit;
//}



class MemUsage {
  protected:
    std::vector<std::string> names;
    std::vector<double> mem_curr;
    std::vector<double> mem_max;

  public:

    /// return the program usable ram in bytes.
    static size_t get_usable_mem() {
       FILE *meminfo = fopen("/proc/meminfo", "r");
       if(meminfo == NULL) {
         // if open failed, throw exception.
         ::std::stringstream ss;
         int myerr = errno;
         ss << "ERROR opening /proc/meminfo: " << myerr << ": " << strerror(myerr);

         throw ::bliss::io::IOException(ss.str());
       }
       char line[256];
       size_t total = 0;
       size_t free_ram = 0;
       size_t buffered = 0;
       size_t cached = 0;
       size_t active_file = 0;

       while(fgets(line, sizeof(line), meminfo))
        {
         sscanf(line, "MemTotal: %lu kB", &total);
         sscanf(line, "MemFree: %lu kB", &free_ram);
         sscanf(line, "Buffers: %lu kB", &buffered);
         sscanf(line, "Cached: %lu kB", &cached);
         sscanf(line, "Active(file): %lu kB", &active_file);
        }

        // If we got here, then we couldn't find the proper line in the meminfo file:
        // do something appropriate like return an error code, throw an exception, etc.
        fclose(meminfo);
        size_t avail = free_ram + ::std::max(active_file, cached);
        size_t used = total - avail;  // get the used part
        used += (used >> 2);  // increase used by 25%

        return (total - used) * 1024UL;
    }



    void reset() {
      names.clear();
      mem_curr.clear();
      mem_max.clear();
    }


//============ memory_usage start
    void mark(::std::string const & name) {
      names.push_back(name);
      mem_curr.push_back(::getCurrentRSS());
      mem_max.push_back(::getPeakRSS());
    }
    void collective_mark(::std::string const & name, ::mxx::comm const & comm) {

		// time a barrier.
		comm.barrier();

		mark(name);
    }

    void report(::std::string const & title) {
        auto BtoMB = [](double const & x) { return x / (1024.0 * 1024.0); };

    	std::stringstream output;

        output << std::fixed;
        std::ostream_iterator<std::string> nit(output, ",");
        std::ostream_iterator<double> dit(output, ",");

        output << "[MEM] " << title << "\theader (MB)\t[,";
        std::copy(names.begin(), names.end(), nit);
        output << "]" << std::endl;

        output.precision(3);

        output << "[MEM] " << title << "\tcurr\t[,";
        std::transform(mem_curr.begin(), mem_curr.end(), dit, BtoMB);
        output << "]" << ::std::endl;

        output << "[MEM] " << title << "\tmax\t[,";
        std::transform(mem_max.begin(), mem_max.end(), dit, BtoMB);
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

        auto curr_mins = mxx::min_element(mem_curr, comm);
        auto curr_maxs = mxx::max_element(mem_curr, comm);
        auto curr_means = ::mxx::reduce(mem_curr, 0, ::std::plus<double>(),comm);
        ::std::for_each(mem_curr.begin(), mem_curr.end(), [](double &x) { x = x*x; });
        auto curr_stdevs = ::mxx::reduce(mem_curr, 0, ::std::plus<double>(), comm);

        auto peak_mins = ::mxx::min_element(mem_max, comm);
        auto peak_maxs = ::mxx::max_element(mem_max, comm);
        auto peak_means = ::mxx::reduce(mem_max, 0, ::std::plus<double>(), comm);
        ::std::for_each(mem_max.begin(), mem_max.end(), [](double &x) { x = x*x; });
        auto peak_stdevs = ::mxx::reduce(mem_max, 0, ::std::plus<double>(), comm);

    	int rank = comm.rank();

        if (rank == 0) {
        	int p = comm.size();
            auto get_first = [](const std::pair<double, int>& x){ return x.first / (1024.0 * 1024.0); };
            auto get_second = [](const std::pair<double, int>& x){ return x.second / (1024.0 * 1024.0); };


          ::std::for_each(curr_means.begin(), curr_means.end(), [p](double & x) { x /= p; });
          ::std::transform(curr_stdevs.begin(), curr_stdevs.end(), curr_means.begin(), curr_stdevs.begin(),
                           [p](double const & x, double const & y) { return ::std::sqrt(x / p - y * y); });

          ::std::for_each(peak_means.begin(), peak_means.end(), [p](double & x) { x /= p; });
          ::std::transform(peak_stdevs.begin(), peak_stdevs.end(), peak_means.begin(), peak_stdevs.begin(),
                           [p](double const & x, double const & y) { return ::std::sqrt(x / p - y * y); });

          std::stringstream output;
          std::ostream_iterator<std::string> nit(output, ",");
          std::ostream_iterator<double> dit(output, ",");
          std::ostream_iterator<int> iit(output, ",");

          output << std::fixed;
          output << "[MEM] " << "R " << rank << "/" << p << " ";

          output << "[MEM] " << title << "\theader (MB)\t[,";
          std::copy(names.begin(), names.end(), nit);
          output << "]" << std::endl;

          output.precision(3);

//          output << "[MEM] " << title << "\tcurr_min\t[,";
//          std::transform(curr_mins.begin(), curr_mins.end(), dit, get_first);
//          output << "]" << std::endl;
//          output << "[MEM] " << title << "\tcurr_min_idx\t[,";
//          std::transform(curr_mins.begin(), curr_mins.end(), iit, get_second);
//          output << "]" << std::endl;
//
//          output << "[MEM] " << title << "\tcurr_max\t[,";
//          std::transform(curr_maxs.begin(), curr_maxs.end(), dit, get_first);
//          output << "]" << std::endl;
//          output << "[MEM] " << title << "\tcurr_max_idx\t[,";
//          std::transform(curr_maxs.begin(), curr_maxs.end(), iit, get_second);
//          output << "]" << std::endl;
//
//          output << "[MEM] " << title << "\tcurr_mean\t[,";
//          std::copy(curr_means.begin(), curr_means.end(), dit);
//          output << "]" << std::endl;
//
//          output << "[MEM] " << title << "\tcurr_stdev\t[,";
//          std::copy(curr_stdevs.begin(), curr_stdevs.end(), dit);
//          output << "]" << std::endl;

          output << "[MEM] " << title << "\tpeak_min\t[,";
          std::transform(peak_mins.begin(), peak_mins.end(), dit, get_first);
          output << "]" << std::endl;
          output << "[MEM] " << title << "\tpeak_min_idx\t[,";
          std::transform(peak_mins.begin(), peak_mins.end(), iit, get_second);
          output << "]" << std::endl;

          output << "[MEM] " << title << "\tpeak_max\t[,";
          std::transform(peak_maxs.begin(), peak_maxs.end(), dit, get_first);
          output << "]" << std::endl;
          output << "[MEM] " << title << "\tpeak_max_idx\t[,";
          std::transform(peak_maxs.begin(), peak_maxs.end(), dit, get_second);
          output << "]" << std::endl;

          output << "[MEM] " << title << "\tpeak_mean\t[,";
          std::copy(peak_means.begin(), peak_means.end(), dit);
          output << "]" << std::endl;

          output << "[MEM] " << title << "\tpeak_stdev\t[,";
          std::copy(peak_stdevs.begin(), peak_stdevs.end(), dit);
          output << "]";


          fflush(stdout);
          printf("%s\n", output.str().c_str());
          fflush(stdout);
        }
        comm.barrier();
    }
#endif

    void report(::std::string const & title, ::mxx::comm const & comm) {


    	auto curr_mins = ::mxx::reduce(mem_curr, 0,
    			[](double const & x, double const & y) { return ::std::min(x, y); }, comm);
        auto curr_maxs = ::mxx::reduce(mem_curr, 0,
        		[](double const & x, double const & y) { return ::std::max(x, y); }, comm);
        auto curr_means = ::mxx::reduce(mem_curr, 0, ::std::plus<double>(), comm);
        ::std::for_each(mem_curr.begin(), mem_curr.end(), [](double &x) { x = x*x; });
        auto curr_stdevs = ::mxx::reduce(mem_curr, 0, ::std::plus<double>(), comm);

        auto peak_mins = ::mxx::reduce(mem_max, 0,
        		[](double const & x, double const & y) { return ::std::min(x, y); }, comm);
        auto peak_maxs = ::mxx::reduce(mem_max, 0,
        		[](double const & x, double const & y) { return ::std::max(x, y); }, comm);
        auto peak_means = ::mxx::reduce(mem_max, 0, ::std::plus<double>(), comm);
        ::std::for_each(mem_max.begin(), mem_max.end(), [](double &x) { x = x*x; });
        auto peak_stdevs = ::mxx::reduce(mem_max, 0, ::std::plus<double>(), comm);

    	int rank = comm.rank();

        if (rank == 0) {
        	int p = comm.size();
            auto BtoMB = [](double const & x) { return x / (1024.0 * 1024.0); };

          ::std::for_each(curr_means.begin(), curr_means.end(), [p](double & x) { x /= p; });
          ::std::transform(curr_stdevs.begin(), curr_stdevs.end(), curr_means.begin(), curr_stdevs.begin(),
                           [p](double const & x, double const & y) { return ::std::sqrt(x / p - y * y); });

          ::std::for_each(peak_means.begin(), peak_means.end(), [p](double & x) { x /= p; });
          ::std::transform(peak_stdevs.begin(), peak_stdevs.end(), peak_means.begin(), peak_stdevs.begin(),
                           [p](double const & x, double const & y) { return ::std::sqrt(x / p - y * y); });

          std::stringstream output;
          std::ostream_iterator<std::string> nit(output, ",");
          std::ostream_iterator<double> dit(output, ",");

          output << std::fixed;
          output << "[MEM] " << "R " << rank << "/" << p << std::endl;

          output << "[MEM] " << title << "\theader (MB)\t[,";
          std::copy(names.begin(), names.end(), nit);
          output << "]" << std::endl;

          output.precision(3);
          output << "[MEM] " << title << "\tcurr_min\t[,";
          std::transform(curr_mins.begin(), curr_mins.end(), dit, BtoMB);
          output << "]" << std::endl;

          output << "[MEM] " << title << "\tcurr_max\t[,";
          std::transform(curr_maxs.begin(), curr_maxs.end(), dit, BtoMB);
          output << "]" << std::endl;

          output << "[MEM] " << title << "\tcurr_mean\t[,";
          std::transform(curr_means.begin(), curr_means.end(), dit, BtoMB);
          output << "]" << std::endl;

          output << "[MEM] " << title << "\tcurr_stdev\t[,";
          std::transform(curr_stdevs.begin(), curr_stdevs.end(), dit, BtoMB);
          output << "]" << std::endl;

          output << "[MEM] " << title << "\tpeak_min\t[,";
          std::transform(peak_mins.begin(), peak_mins.end(), dit, BtoMB);
          output << "]" << std::endl;

          output << "[MEM] " << title << "\tpeak_max\t[,";
          std::transform(peak_maxs.begin(), peak_maxs.end(), dit, BtoMB);
          output << "]" << std::endl;

          output << "[MEM] " << title << "\tpeak_mean\t[,";
          std::transform(peak_means.begin(), peak_means.end(), dit, BtoMB);
          output << "]" << std::endl;

          output << "[MEM] " << title << "\tpeak_stdev\t[,";
          std::transform(peak_stdevs.begin(), peak_stdevs.end(), dit, BtoMB);
          output << "]";


          fflush(stdout);
          printf("%s\n", output.str().c_str());
          fflush(stdout);
        }
        comm.barrier();
    }

};


#define BL_MEMUSE_INIT(title)      MemUsage title##_memusage;


#define BL_MEMUSE_RESET(title)     do {  title##_memusage.reset(); } while (0)

#define BL_MEMUSE_COLLECTIVE_MARK(title, name, comm) do { title##_memusage.collective_mark(name, comm); } while (0)
#define BL_MEMUSE_MARK(title, name) do { title##_memusage.mark(name); } while (0)
#define BL_MEMUSE_REPORT(title) do { title##_memusage.report(#title); } while (0)

#if 0
// do not use this function for now.  mxx::min_element and max_element has invalid read problem, reported by valgrind.
#define BL_MEMUSE_REPORT_MPI_LOC(title, comm) do { title##_memusage.report_loc(#title, comm); } while (0)
#endif

#define BL_MEMUSE_REPORT_MPI(title, comm) do { title##_memusage.report(#title, comm); } while (0)


#else

#define BL_MEMUSE_INIT(title)
#define BL_MEMUSE_RESET(title)
#define BL_MEMUSE_COLLECTIVE_MARK(title, name, comm)
#define BL_MEMUSE_MARK(title, name)
#define BL_MEMUSE_REPORT(title)
#define BL_MEMUSE_REPORT_MPI(title, comm)

#endif


#endif /* SRC_UTILS_MEMORY_USAGE_HPP_ */
