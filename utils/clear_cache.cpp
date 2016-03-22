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

#include "utils/memory_usage.hpp"

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

  printf("begin clearing %lu bytes\n", avail);
  size_t iter_cleared = 0;

  while ((chunk >= minchunk) && (rem > minchunk)) {
    nchunks = rem / chunk;

    iter_cleared = 0;
    lj = 0;
#pragma omp parallel for shared(nchunks, chunk, dummy) reduction(+:lj, iter_cleared)
    for (size_t i = 0; i < nchunks; ++i) {

      // (c|m)alloc/free seems to be optimized out.  using new works.
      size_t * ptr = new size_t[(chunk / sizeof(size_t))];

      iter_cleared += chunk;
      memset(ptr, 0, chunk);
      ptr[0] = i;

#pragma omp critical
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

  clear_cache();

  return 0;
}
