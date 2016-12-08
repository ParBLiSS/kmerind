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
 * @file    clear_cache.cpp
 * @ingroup
 * @author  tpan
 * @brief   clears cache on a linux system
 * @details for use before doing file io benchmark.
 *
 * Copyright (c) 2016 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#include "bliss-config.hpp"

#include "utils/memory_usage.hpp"
#include <cstring>

#ifdef USE_MPI
#include <mxx/env.hpp>
#include <mxx/comm.hpp>
#endif

#ifdef USE_OPENMP
#include <omp.h>
#endif

/**
 * @brief  clear the disk cache in linux by allocating a bunch of memory.
 * @details  in 1MB chunks to workaround memory fragmentation preventing allocation.
 *           can potentially use swap, or if there is not enough physical + swap, be killed.
 *
 */
void clear_cache() {
  size_t avail = MemUsage::get_usable_mem();
  size_t rem = avail;

  size_t minchunk = 1UL << 24;    // 16MB chunks
  size_t chunk = minchunk;

  size_t maxchunk = std::min(1UL << 36, (avail >> 4));
  for ( ;chunk < maxchunk; chunk <<= 1) ;   // keep increasing chunks until larger than 1/16 of avail, or until 1GB.

  std::vector<size_t *> dummy;

  size_t nchunks;
  size_t j = 0, lj;

#if defined(USE_OPENMP)
	int max_threads = omp_get_max_threads();

	if (max_threads > 8) max_threads -= 2;
	else if (max_threads > 4) max_threads -= 1;
#else
	int max_threads = 1;
#endif


  printf("begin clearing %lu bytes using %d threads\n", avail, max_threads);
  size_t iter_cleared = 0;

  while ((chunk >= minchunk) && (rem > minchunk)) {
    nchunks = rem / chunk;

    iter_cleared = 0;
    lj = 0;
#if defined(USE_OPENMP)
#pragma omp parallel for num_threads(max_threads) shared(nchunks, chunk, dummy) reduction(+:lj, iter_cleared)
#endif
      for (size_t i = 0; i < nchunks; ++i) {
      // (c|m)alloc/free seems to be optimized out.  using new works.
      size_t * ptr = new size_t[(chunk / sizeof(size_t))];

      iter_cleared += chunk;
      memset(ptr, 0, chunk);
      ptr[0] = i;

#if defined(USE_OPENMP)
#pragma omp critical
#endif
      {
        dummy.push_back(ptr);
      }

      ++lj;
    }

    j += lj;

    //    // check for available memory, every 64 MBs
    //    if ((i % 64) == 0) {
    //      avail = MemUsage::get_usable_mem();
    //      if (avail < (64 << 20)) break;   // less than 64MB available.
    //    }
    rem -= iter_cleared;

    printf("cleared %lu bytes using %lu chunk %lu bytes. total cleared %lu bytes, rem %lu bytes \n", iter_cleared, nchunks, chunk, avail - rem, rem);
    fflush(stdout);


    // reduce the size of the chunk by 4
    chunk >>= 4;

  }
  printf("finished clearing %lu/%lu bytes with %lu remaining\n", avail - rem, avail, rem);
  fflush(stdout);

  size_t sum = 0;
  size_t ii = 0;
  size_t *ptr;
  for (; ii < dummy.size(); ++ii) {
    ptr = dummy[ii];

    if (ptr != nullptr) {
//      if ((ii % 16) == 0) {
//        printf("%lu ", ii);  // 16 outputs.
//        fflush(stdout);
//      }
      sum += ptr[ii >> 10];
      delete [] ptr;
      //free(ptr);
      dummy[ii] = nullptr;
    } else {
      break;
    }
  }
  printf("\n");
  printf("disk cache cleared (dummy %lu). %lu blocks %lu bytes\n", sum, j, avail - rem);
}


int main(int argc, char** argv) {

#ifdef USE_MPI
  ::mxx::env e(argc, argv);
  ::mxx::comm world;


  if (world.size() > 1) {
    ::mxx::comm node = world.split_shared();

    if (node.rank() == 0) {
	printf("rank %d clearing\n", world.rank());
	clear_cache();
    } else {
	//printf("rank %d waiting\n", world.rank());
    }

    node.barrier();
    if (node.rank() == 0) 
	printf("rank %d cache cleared\n", world.rank());
  } else {
//	std::cout << "world size is " << world.size() << std::endl;

//#else
//	std::cout << "non-mpi." << std::endl;
#endif
	
  clear_cache();

#ifdef USE_MPI
}
#endif

  return 0;
}
