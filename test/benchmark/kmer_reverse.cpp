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



#include "common/kmer.hpp"
#include "common/test/kmer_reverse_helper.hpp"

#include <random>
#include <iostream>
#include <cstring>

int main(int argc, char** argv) {

  bool bitgroup = true;
  if (argc > 1) {
    bitgroup = (strcmp(argv[1], "-bg") == 0);
  }

  using KmerType = ::bliss::common::Kmer< 64, bliss::common::DNA,   uint64_t>;

  KmerType km, tmp, rev;

  for (size_t i = 0; i < KmerType::size; ++i) {
    km.nextFromChar(rand() % KmerType::KmerAlphabet::SIZE);
  }




  bliss::common::test::KmerReverseHelper<KmerType> helper;

  size_t iters = 10000000;
  for (size_t i = 0; i < iters; ++i) {

    if (bitgroup)
      tmp = km.reverse();
    else
      tmp = helper.reverse_bswap(km);


    rev ^= tmp;

    km.nextFromChar(rand() % KmerType::KmerAlphabet::SIZE);

    if (i % 100000 == 0) std::cout << i / 100000 << "%" << std::endl;

  }

  ::std::cout << rev << std::endl;


}
