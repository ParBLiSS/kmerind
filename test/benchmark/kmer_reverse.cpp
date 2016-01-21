/**
 * @file    kmer_reverse.cpp
 * @ingroup
 * @author  tpan
 * @brief
 * @details
 *
 * Copyright (c) 2016 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
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
