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

  size_t blocks = avail >> 20;   // in 2^20 (1MB) blocks);

  std::vector<size_t *> dummy(blocks, nullptr);

  printf("attempting to allocate %lu bytes\n", avail);

  // (c|m)alloc/free seems to be optimized out.  using new works.

  size_t i = 0;
  size_t *ptr;
  for (; i < blocks; ++i) {

    if ((i % (blocks >> 4)) == 0) {
      printf("%lu ", i);  // 16 outputs.
      fflush(stdout);
    }

    // check for available memory, every 64 MBs
    if ((i % 64) == 0) {
      avail = MemUsage::get_usable_mem();
      if (avail < (64 << 20)) break;   // less than 64MB available.
    }

//    // try malloc.
//    ptr = (size_t*)calloc((1 << 20) / sizeof(size_t), sizeof(size_t));
//    if (ptr == nullptr) {  // out of mem
//      break;
//    }
    ptr = new size_t[((1 << 20) / sizeof(size_t))];

    memset(ptr, 0, 1 << 20);
    ptr[i >> 10] = i;

    dummy[i] = ptr;
  }
  printf("\n");

  printf("mem allocated. %lu blocks %lu bytes\n", i, i * (1 << 20));
  fflush(stdout);

  size_t sum = 0;
  size_t ii = 0;
  for (; ii < dummy.size(); ++ii) {
    ptr = dummy[ii];

    if (ptr != nullptr) {
      if ((ii % (blocks >> 4)) == 0) {
        printf("%lu ", ptr[ii >> 10]);  // 16 outputs.
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
  printf("disk cache cleared (dummy %lu). %lu blocks %lu bytes\n", sum, ii, ii * (1 << 20));
}


int main(int argc, char** argv) {

  clear_cache();

  return 0;
}
