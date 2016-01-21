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
 * @file    kmer_reverse.cpp
 * @ingroup
 * @author  tpan
 * @brief
 * @details
 *
 */


#include <random>
#include <iostream>
#include <cstring>
#include <cassert>

/// DEREFERENCING POINTERS IS SLIGHTLY SLOWER.
uint8_t use_pointer(uint8_t * const & out, uint8_t const * const & in, size_t len, uint8_t blah) {
  assert(len == 8);
  *(reinterpret_cast<uint64_t * const &>(out)) = *(reinterpret_cast<uint64_t const * const &>(in));
  return rand() % 255;
}

uint8_t use_reference(uint64_t & out, uint64_t const & in, uint8_t blah) {
  out = in;
  return rand() % 255;
}

constexpr size_t iters = 10000000;


int main(int argc, char** argv) {

  bool bitgroup = true;
  if (argc > 1) {
    bitgroup = (strcmp(argv[1], "-ref") == 0);
  }


  uint64_t *data = new uint64_t[iters];
  for (size_t i = 0; i < iters; ++i) {
    data[i] = rand() % 100;
  }
  uint64_t *output = new uint64_t[iters];
  memset(output, 0, sizeof(uint64_t) * iters);

  uint8_t off = 0;
  for (size_t i = 0; i < iters; ++i) {

    if (bitgroup)
      off = use_reference(output[i], data[i], off);
    else {

      off = use_pointer(reinterpret_cast<uint8_t*>(&(output[i])), reinterpret_cast<uint8_t const *>(&(data[i])), 8, off);
    }

    if (i % 100000 == 0) std::cout << i / 100000 << "%" << std::endl;
  }

  bool same = true;
  for (size_t i = 0; i < iters; ++i) {
    same &= (output[i] == data[i]);
  }
  if (!same ) std::cerr << "not matched" << std::endl;

  delete [] data;
  delete [] output;

}
