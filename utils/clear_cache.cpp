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

/**
 * @brief  clear the disk cache in linux by allocating a bunch of memory.
 * @details  in 1MB chunks to workaround memory fragmentation preventing allocation.
 *           can potentially use swap, or if there is not enough physical + swap, be killed.
 *
 */
void clear_cache() {
  size_t avail = MemUsage::get_usable_mem();
  size_t rem = avail;

  size_t minchunk = 1 << 23;
  size_t chunk = minchunk;

  while (chunk < (avail >> 4)) chunk <<= 1;   // keep increasing chunks until larger than 1/16 of avail.

  std::vector<size_t *> dummy;

  size_t nchunks;
  size_t i;
  size_t *ptr;
  size_t j = 0;

  printf("begin clearing %lu bytes\n", avail);

  while ((chunk >= minchunk) && (rem > chunk)) {
    nchunks = rem / chunk;

    for (i = 0; i < nchunks; ++i, ++j) {
      if ((j % 16) == 0) {
        printf("%lu ", j);  // 16 outputs.
        fflush(stdout);
      }

      // (c|m)alloc/free seems to be optimized out.  using new works.
      ptr = new size_t[(chunk / sizeof(size_t))];

      rem -= chunk;
      memset(ptr, 0, chunk);
      ptr[0] = i;

      dummy.push_back(ptr);
    }

    //    // check for available memory, every 64 MBs
    //    if ((i % 64) == 0) {
    //      avail = MemUsage::get_usable_mem();
    //      if (avail < (64 << 20)) break;   // less than 64MB available.
    //    }

    // reduce the size of the chunk by 4
    chunk >>= 2;
  }
  printf("\n");
  printf("attempting to allocate %lu bytes\n", avail - rem);
  fflush(stdout);

  size_t sum = 0;
  size_t ii = 0;
  for (; ii < dummy.size(); ++ii) {
    ptr = dummy[ii];

    if (ptr != nullptr) {
      if ((ii % 16) == 0) {
        printf("%lu ", ii);  // 16 outputs.
        fflush(stdout);
      }
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
