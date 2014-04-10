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
    template<typename Sequence, typename KmerIndex, typename Encoding,
    typename HasQual = typename std::enable_if<std::is_same<KmerIndex, KmerIndexElementWithIdAndQuality>::value>::type >
    constexpr std::array<typename Encoding::value_type, Encoding::size>
      bliss::index::generate_qual<Sequence, KmerIndex, Encoding, HasQual >::lut;
  }
}

