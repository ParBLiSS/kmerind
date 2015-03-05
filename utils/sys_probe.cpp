/**
 * sys_probe.cpp
 *
 *  Created on: Jan 22, 2014
 *      Author: Tony Pan <tpan7@gatech.edu>
 */

#include <unistd.h>
#include <errno.h>  // for errno
#include <stdio.h>  // for stderr
#include <stdlib.h> // for EXIT_FAILURE
#include <string.h> // for strerror
#include <stdint.h> // for uint32_t
#include <sys/stat.h> // for stat

#include "utils/logging.h"

int checkNProcs()
{
  long nprocs = -1;
  long nprocs_max = -1;

  INFOF("Does not yet compute the number of CPUs\n");

#ifdef _SC_NPROCESSORS_ONLN
  nprocs = sysconf(_SC_NPROCESSORS_ONLN);
  if (nprocs < 1)
  {
    INFOF("Could not determine number of Cores online:\n%s\n",
            strerror(errno));
    return EXIT_FAILURE;
  }
  nprocs_max = sysconf(_SC_NPROCESSORS_CONF);
  if (nprocs_max < 1)
  {
    INFOF("Could not determine number of Cores configured:\n%s\n",
            strerror(errno));
    return EXIT_FAILURE;
  }
  INFOF("%ld of %ld processors online\n", nprocs, nprocs_max);
  return EXIT_SUCCESS;
#else
  INFOF("Could not determine number of Cores\n");
  return EXIT_FAILURE;
#endif
  return 0;
}

int checkMainMem()
{
  long pages = sysconf(_SC_PHYS_PAGES);
  long page_size = sysconf(_SC_PAGE_SIZE);

  INFOF("%ld pages at %ld bytes per page, total %ld bytes\n", pages,
          page_size, pages * page_size);
  return 0;
}

int checkMPIBuffer()
{
  INFOF("MPI tests not yet implemented\n");
  return 0;
}

int checkFileSystem()
{
  INFOF("%d buffer size\n", BUFSIZ);

  struct stat fileStat;
  if (stat(".", &fileStat) < 0)
    return EXIT_FAILURE;

  INFOF("on disk file block size: %ld\n", fileStat.st_size);

  return EXIT_SUCCESS;
}

// from http://en.wikipedia.org/wiki/CPUID
void cpuid(unsigned info, unsigned *eax, unsigned *ebx, unsigned *ecx,
           unsigned *edx)
{
  __asm__(
      "xchg %%ebx, %%edi;" /* 32bit PIC: don't clobber ebx */
      "cpuid;"
      "xchg %%ebx, %%edi;"
      :"=a" (*eax), "=D" (*ebx), "=c" (*ecx), "=d" (*edx)
      :"0" (info)
  );
}

void cpu_features()
{
  uint32_t eax, ebx, ecx, edx;

  eax = 1; // cpu features

  cpuid(1, &eax, &ebx, &ecx, &edx);  // get cpu features
  bool mmx = edx && (1 << 23);
  bool sse = edx && (1 << 25);
  bool sse2 = edx && (1 << 26);
  bool sse3 = ecx && (1);
  bool ssse3 = ecx && (1 << 9);
  bool sse41 = ecx && (1 << 19);
  bool sse42 = ecx && (1 << 20);
  bool avx = ecx && (1 << 28);
  bool ia64 = edx && (1 << 30);
  bool hyperthread = edx && (1 << 28);
  bool fma = ecx && (1 << 12);
  bool cx16 = ecx && (1 << 13);
  bool popcnt = ecx && (1 << 23);

  INFOF("mmx: %s\n", (mmx ? "true" : "false"));
  INFOF("sse: %s\n", (sse ? "true" : "false"));
  INFOF("sse2: %s\n", (sse2 ? "true" : "false"));
  INFOF("sse3: %s\n", (sse3 ? "true" : "false"));
  INFOF("ssse3: %s\n", (ssse3 ? "true" : "false"));
  INFOF("sse41: %s\n", (sse41 ? "true" : "false"));
  INFOF("sse42: %s\n", (sse42 ? "true" : "false"));
  INFOF("avx: %s\n", (avx ? "true" : "false"));
  INFOF("ia64: %s\n", (ia64 ? "true" : "false"));
  INFOF("hyperthread: %s\n", (hyperthread ? "true" : "false"));
  INFOF("fma: %s\n", (fma ? "true" : "false"));
  INFOF("cx16: %s\n", (cx16 ? "true" : "false"));
  INFOF("popcnt: %s\n", (popcnt ? "true" : "false"));

}

void i386_cpuid_caches()
{
  int i;
  for (i = 0; i < 32; i++)
  {

    // Variables to hold the contents of the 4 i386 legacy registers
    uint32_t eax, ebx, ecx, edx;

    eax = 4; // get cache info
    ecx = i; // cache id

    __asm__ (
        "cpuid" // call i386 cpuid instruction
        : "+a" (eax)// contains the cpuid command code, 4 for cache query
        , "=b" (ebx)
        , "+c" (ecx)// contains the cache id
        , "=d" (edx)
    );
    // generates output in 4 registers eax, ebx, ecx and edx

    // taken from http://download.intel.com/products/processor/manual/325462.pdf Vol. 2A 3-149
    int cache_type = eax & 0x1F;

    if (cache_type == 0) // end of valid cache identifiers
      break;

    const char * cache_type_string;
    switch (cache_type)
    {
      case 1:
        cache_type_string = "Data Cache";
        break;
      case 2:
        cache_type_string = "Instruction Cache";
        break;
      case 3:
        cache_type_string = "Unified Cache";
        break;
      default:
        cache_type_string = "Unknown Type Cache";
        break;
    }

    int cache_level = (eax >>= 5) & 0x7;

    int cache_is_self_initializing = (eax >>= 3) & 0x1; // does not need SW initialization
    int cache_is_fully_associative = (eax >>= 1) & 0x1;

    // taken from http://download.intel.com/products/processor/manual/325462.pdf 3-166 Vol. 2A
    // ebx contains 3 integers of 10, 10 and 12 bits respectively
    unsigned int cache_sets = ecx + 1;
    unsigned int cache_coherency_line_size = (ebx & 0xFFF) + 1;
    unsigned int cache_physical_line_partitions = ((ebx >>= 12) & 0x3FF) + 1;
    unsigned int cache_ways_of_associativity = ((ebx >>= 10) & 0x3FF) + 1;

    // Total cache size is the product
    size_t cache_total_size = cache_ways_of_associativity
        * cache_physical_line_partitions * cache_coherency_line_size
        * cache_sets;

    INFOF("Cache ID %d:\n"
           "- Level: %d\n"
           "- Type: %s\n"
           "- Sets: %d\n"
           "- System Coherency Line Size: %d bytes\n"
           "- Physical Line partitions: %d\n"
           "- Ways of associativity: %d\n"
           "- Total Size: %zu bytes (%zu kb)\n"
           "- Is fully associative: %s\n"
           "- Is Self Initializing: %s\n"
           "\n",
           i, cache_level, cache_type_string, cache_sets,
           cache_coherency_line_size, cache_physical_line_partitions,
           cache_ways_of_associativity, cache_total_size,
           cache_total_size >> 10,
           cache_is_fully_associative ? "true" : "false",
           cache_is_self_initializing ? "true" : "false");
  }
}

int main(int argc, char* argv[])
{

  checkNProcs();
  checkMainMem();
  checkFileSystem();

  checkMPIBuffer();
  cpu_features();
  i386_cpuid_caches();

  return 0;
}

