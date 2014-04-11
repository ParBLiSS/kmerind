/**
 * @file		kmer_index_functors.cpp
 * @ingroup
 * @author	tpan
 * @brief
 * @details
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */

#include "index/kmer_index_functors.hpp"

namespace bliss
{
  namespace index
  {
    // can't use auto keyword.  declare and initialize in class declaration
    // then "define" but not initialize outside class declaration, again.
    template<typename Sequence, typename KmerSize, typename QualityType, typename Encoding>
    constexpr std::array<typename Encoding::value_type, Encoding::size>
      bliss::index::generate_qual<Sequence, KmerSize, QualityType, Encoding>::lut;
  }
}

