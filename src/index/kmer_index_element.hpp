/**
 * @file		kmer_index_element.hpp
 * @ingroup
 * @author	tpan
 * @brief
 * @details
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef KMER_INDEX_ELEMENT_HPP_
#define KMER_INDEX_ELEMENT_HPP_

namespace bliss
{
  namespace index
  {

    /**
     * @class			bliss::io::KMerIndexElement
     * @brief
     * @details
     *
     */

    template<typename T1, typename T2 = double>
    struct kmer_index_element
    {
        typedef T1 kmer_type;
        typedef T2 qual_type;

        bliss::iterator::read_id id;
        T1 kmer;
        T2 qual;
    };

  } /* namespace io */
} /* namespace bliss */

#endif /* KMER_INDEX_ELEMENT_HPP_ */
