/**
 * @file		kmer_functors.cpp
 * @ingroup
 * @author	tpan
 * @brief
 * @details
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */

#include "index/kmer_functors.hpp"

namespace bliss
{
  namespace index
  {
    // can't use auto keyword.  declare and initialize in class declaration
    // then "define" but not initialize outside class declaration, again.
    template<typename ENCODING, typename Iterator, typename TO, int K>
    constexpr std::array<typename ENCODING::value_type, ENCODING::size> bliss::index::generate_qual<
        ENCODING, Iterator, TO, K>::lut;
  }
}

