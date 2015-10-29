/**
 * @file    mxx_support.hpp
 * @ingroup
 * @author  tpan
 * @brief
 * @details
 *
 * Copyright (c) 2015 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef SRC_IO_MXX_SUPPORT_HPP_
#define SRC_IO_MXX_SUPPORT_HPP_

#include <mxx/datatypes.hpp>

#include "common/kmer.hpp"

//#include "utils/system_utils.hpp"

#include <algorithm>  // for std::min

#include <unistd.h>   // for gethostname

#include "partition/range.hpp"
#include "common/sequence.hpp"
#include "debruijn/de_bruijn_node_trait.hpp"

namespace mxx {

  template<unsigned int size, typename A, typename WT>
  class datatype<typename bliss::common::Kmer<size, A, WT> > :
    public datatype_contiguous<typename bliss::common::Kmer<size, A, WT>::KmerWordType,
      bliss::common::Kmer<size, A, WT>::nWords> {};

  template<unsigned int size, typename A, typename WT>
  class datatype<const typename bliss::common::Kmer<size, A, WT> > :
    public datatype_contiguous<typename bliss::common::Kmer<size, A, WT>::KmerWordType,
      bliss::common::Kmer<size, A, WT>::nWords> {};

  template<>
  class datatype<bliss::common::LongSequenceKmerId > :
    public datatype<decltype(bliss::common::LongSequenceKmerId::id)> {};
  template<>
  class datatype<const bliss::common::LongSequenceKmerId > :
    public datatype<decltype(bliss::common::LongSequenceKmerId::id)> {};

  template<>
  class datatype<bliss::common::ShortSequenceKmerId > :
    public datatype<decltype(bliss::common::ShortSequenceKmerId::id)> {};

  template<>
  class datatype<const bliss::common::ShortSequenceKmerId > :
    public datatype<decltype(bliss::common::ShortSequenceKmerId::id)> {};


  template <typename T>
  class datatype<bliss::partition::range<T> > :
    public datatype_contiguous<T, sizeof(bliss::partition::range<T>) / sizeof(T)> {};

  template <typename T>
  class datatype<const bliss::partition::range<T> > :
    public datatype_contiguous<T, sizeof(bliss::partition::range<T>) / sizeof(T)> {};


  template <typename A, typename T>
  class datatype<bliss::de_bruijn::node::edge_counts<A, T> > :
    public datatype<decltype(bliss::de_bruijn::node::edge_counts<A, T>::counts)> {};

  template <typename A, typename T>
  class datatype<const bliss::de_bruijn::node::edge_counts<A, T> > :
    public datatype<decltype(bliss::de_bruijn::node::edge_counts<A, T>::counts)> {};

  template <typename A>
  class datatype<bliss::de_bruijn::node::edge_exists<A> > :
    public datatype<decltype(bliss::de_bruijn::node::edge_exists<A>::counts)> {};

  template <typename A>
  class datatype<const bliss::de_bruijn::node::edge_exists<A> > :
    public datatype<decltype(bliss::de_bruijn::node::edge_exists<A>::counts)> {};

}

std::ostream &operator<<(std::ostream &os, uint8_t const &t) {
  return os << static_cast<uint32_t>(t);
}
std::ostream &operator<<(std::ostream &os, int8_t const &t) {
  return os << static_cast<int32_t>(t);
}

template <typename T1, typename T2>
std::ostream &operator<<(std::ostream &os, std::pair<T1, T2> const &t) {
    return os << t.first << ":" << t.second;
}

#endif /* SRC_IO_MXX_SUPPORT_HPP_ */
