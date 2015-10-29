/**
 * @file    system_utils.hpp
 * @ingroup
 * @author  tpan
 * @brief
 * @details
 *
 * Copyright (c) 2015 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef SRC_WIP_SYSTEM_UTILS_HPP_
#define SRC_WIP_SYSTEM_UTILS_HPP_

#include "bliss-logger_config.hpp"

#if BENCHMARK == 1

//#include <cstdio>   // printf
#include <vector>
#include <string>
#include <algorithm>  // std::min

#include "utils/logging.h"

#include "stdlib.h"
#include "stdio.h"
#include "sys/types.h"
#include "sys/sysinfo.h"

// implement from http://linux.die.net/man/2/getrusage
#include <sys/time.h>
#include <sys/resource.h>


void get_mem_allocated(long long &phyMemUsed, long long &swapUsed) {
  struct rusage u;
  getrusage(RUSAGE_SELF, &u);
  phyMemUsed = u.ru_maxrss;
}


void get_mem_allocated_sysinfo(long long &phyMemUsed, long long &swapUsed) {
  //from http://stackoverflow.com/questions/63166/how-to-determine-cpu-and-memory-consumption-from-inside-a-process
  struct sysinfo memInfo;
  sysinfo (&memInfo);

  phyMemUsed = (memInfo.totalram - memInfo.freeram) * memInfo.mem_unit;
  swapUsed = (memInfo.totalswap - memInfo.freeswap) * memInfo.mem_unit;
}

void get_mem_consumed_sysinfo(const long long &phyMemBaseline, const long long &swapBaseline, long long &phyMemUsed, long long &swapUsed) {
  long long phy;
  long long swap;
  get_mem_allocated_sysinfo(phy, swap);

  phyMemUsed = phy - phyMemBaseline;
  swapUsed = swap - swapBaseline;
}


#define MEM_INIT(memsession)     long long memsession##_vmem1, memsession##_pmem1, memsession##_vmem2, memsession##_pmem2; \
                                 std::vector<std::string> memsession##_memnames; \
                                 std::vector<double> memsession##_vmems; \
                                 std::vector<double> memsession##_pmems; \
                                 std::vector<double> memsession##_memcounts;

#define MEM_RESET(memsession)     do {  memsession##_memnames.clear(); \
                                       memsession##_vmems.clear(); \
                                       memsession##_pmems.clear(); \
                                       memsession##_count.clear(); \
                                 } while (0)

#define MEM_START(memsession)     do { get_mem_allocated(memsession##_pmem1, memsession##_vmem1); \
                                 } while (0)

#define MEM_END(memsession, name, n_elem) do { get_mem_allocated(memsession##_pmem2, memsession##_vmem2); \
                                             memsession##_memnames.push_back(name); \
                                             memsession##_pmems.push_back(memsession##_pmem2 - memsession##_pmem1); \
                                             memsession##_vmems.push_back(memsession##_vmem2 - memsession##_vmem1); \
                                             memsession##_memcounts.push_back(n_elem); \
                                        } while (0)


#define MEM_REPORT(memsession, rank) \
        do { \
          std::stringstream output; \
          output << std::fixed; \
          output << "R " << rank << " " << #memsession << " header\t["; \
          std::ostream_iterator<std::string> nit(output, ","); \
          std::copy(memsession##_memnames.begin(), memsession##_memnames.end(), nit); \
          output << "]\n" << #memsession << "\tpmem\t\t["; \
          output.precision(9); \
          std::ostream_iterator<double> dit(output, ","); \
          std::copy(memsession##_pmems.begin(), memsession##_pmems.end(), dit); \
          output << "]\n" << #memsession << "\tvmem\t\t["; \
          output.precision(9); \
          std::copy(memsession##_vmems.begin(), memsession##_vmems.end(), dit); \
          output << "]\n" << #memsession << "\tcount\t\t["; \
          output.precision(0); \
          std::ostream_iterator<double> cit(output, ","); \
          std::copy(memsession##_memcounts.begin(), memsession##_memcounts.end(), cit); \
          output << "]"; \
          fflush(stdout); \
          INFOF("%s\n", output.str().c_str()); \
          fflush(stdout); \
        } while (0)



#define MEM_REPORT_MPI(memsession, reduce_op, rank, comm) \
        do { \
          int p; \
          MPI_Comm_size(comm, &p); \
          \
          auto pmem_mins = reduce_op(memsession##_pmems, [](double const &x, double const &y) { return ::std::min(x,y); }, comm, 0); \
          auto pmem_maxs = reduce_op(memsession##_pmems, [](double const &x, double const &y) { return ::std::max(x,y); }, comm, 0); \
          auto pmem_means = reduce_op(memsession##_pmems, ::std::plus<double>(), comm, 0); \
          ::std::for_each(memsession##_pmems.begin(), memsession##_pmems.end(), [](double &x) { x = x*x; }); \
          auto pmem_stdevs = reduce_op(memsession##_pmems, ::std::plus<double>(), comm, 0); \
          \
          auto vmem_mins = reduce_op(memsession##_vmems, [](double const &x, double const &y) { return ::std::min(x,y); }, comm, 0); \
          auto vmem_maxs = reduce_op(memsession##_vmems, [](double const &x, double const &y) { return ::std::max(x,y); }, comm, 0); \
          auto vmem_means = reduce_op(memsession##_vmems, ::std::plus<double>(), comm, 0); \
          ::std::for_each(memsession##_vmems.begin(), memsession##_vmems.end(), [](double &x) { x = x*x; }); \
          auto vmem_stdevs = reduce_op(memsession##_vmems, ::std::plus<double>(), comm, 0); \
          \
          auto cnt_mins = reduce_op(memsession##_memcounts, [](double const &x, double const &y) { return ::std::min(x,y); }, comm, 0); \
          auto cnt_maxs = reduce_op(memsession##_memcounts, [](double const &x, double const &y) { return ::std::max(x,y); }, comm, 0); \
          auto cnt_means = reduce_op(memsession##_memcounts, ::std::plus<double>(), comm, 0); \
          ::std::for_each(memsession##_memcounts.begin(), memsession##_memcounts.end(), [](double &x) { x = x*x; }); \
          auto cnt_stdevs = reduce_op(memsession##_memcounts, ::std::plus<double>(), comm, 0); \
          \
          if (rank == 0) { \
            ::std::for_each(pmem_means.begin(), pmem_means.end(), [p](double & x) { x /= p; }); \
            ::std::transform(pmem_stdevs.begin(), pmem_stdevs.end(), pmem_means.begin(), pmem_stdevs.begin(), \
                             [p](double const & x, double const & y) { return ::std::sqrt(x / p - y * y); }); \
            ::std::for_each(vmem_means.begin(), vmem_means.end(), [p](double & x) { x /= p; }); \
            ::std::transform(vmem_stdevs.begin(), vmem_stdevs.end(), vmem_means.begin(), vmem_stdevs.begin(), \
                             [p](double const & x, double const & y) { return ::std::sqrt(x / p - y * y); }); \
            ::std::for_each(cnt_means.begin(), cnt_means.end(), [p](double & x) { x /= p; }); \
            ::std::transform(cnt_stdevs.begin(), cnt_stdevs.end(), cnt_means.begin(), cnt_stdevs.begin(), \
                             [p](double const & x, double const & y) { return ::std::sqrt(x / p - y * y); }); \
            \
            std::stringstream output; \
            output << std::fixed; \
            output << "R " << rank << "/" << p << " " << #memsession << " header\t["; \
            std::ostream_iterator<std::string> nit(output, ","); \
            std::copy(memsession##_memnames.begin(), memsession##_memnames.end(), nit); \
            output.precision(0); \
            output << "]\n" << #memsession << "\tpmem_min\t\t["; \
            std::ostream_iterator<double> dit(output, ","); \
            std::copy(pmem_mins.begin(), pmem_mins.end(), dit); \
            output << "]\n" << #memsession << "\tpmem_max\t\t["; \
            std::copy(pmem_maxs.begin(), pmem_maxs.end(), dit); \
            output.precision(2); \
            output << "]\n" << #memsession << "\tpmem_mean\t["; \
            std::copy(pmem_means.begin(), pmem_means.end(), dit); \
            output << "]\n" << #memsession << "\tpmem_stdev\t["; \
            std::copy(pmem_stdevs.begin(), pmem_stdevs.end(), dit); \
            output.precision(0); \
            output << "]\n" << #memsession << "\tvmem_min\t\t["; \
            std::copy(vmem_mins.begin(), vmem_mins.end(), dit); \
            output << "]\n" << #memsession << "\tvmem_max\t\t["; \
            std::copy(vmem_maxs.begin(), vmem_maxs.end(), dit); \
            output.precision(2); \
            output << "]\n" << #memsession << "\tvmem_mean\t["; \
            std::copy(vmem_means.begin(), vmem_means.end(), dit); \
            output << "]\n" << #memsession << "\tvmem_stdev\t["; \
            std::copy(vmem_stdevs.begin(), vmem_stdevs.end(), dit); \
            output.precision(0); \
            output << "]\n" << #memsession << "\tcnt_min\t\t["; \
            std::copy(cnt_mins.begin(), cnt_mins.end(), dit); \
            output << "]\n" << #memsession << "\tcnt_max\t\t["; \
            std::copy(cnt_maxs.begin(), cnt_maxs.end(), dit); \
            output.precision(2); \
            output << "]\n" << #memsession << "\tcnt_mean\t["; \
            std::copy(cnt_means.begin(), cnt_means.end(), dit); \
            output << "]\n" << #memsession << "\tcnt_stdev\t["; \
            std::copy(cnt_stdevs.begin(), cnt_stdevs.end(), dit); \
            output << "]"; \
            fflush(stdout); \
            INFOF("%s\n", output.str().c_str()); \
            fflush(stdout); \
          } \
          MPI_Barrier(comm); \
        } while (0)


#else

#define MEM_INIT(memsession)
#define MEM_RESET(memsession)
#define MEM_START(memsession)
#define MEM_END(memsession, name, n_elem)
#define MEM_REPORT(memsession, rank)
#define MEM_REPORT_MPI(memsession, rank, comm)

#endif

#endif /* SRC_WIP_SYSTEM_UTILS_HPP_ */
