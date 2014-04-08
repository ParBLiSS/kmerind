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
     * @class     bliss::io::kmer_size
     * @brief     represents the size of the kmer.
     * @details
     */
    template<int K>
    struct kmer_size
    {
        constexpr int size = K;
    };

    /**
     * @class			bliss::io::kmer_index_element
     * @brief     basic struct for storing position and quality of a kmer.
     * @details
     *  K is a typename so that we do not use memory.
     */
    template<typename K, typename KmerT, typename IdT, typename QualT>
    struct kmer_index_element
    {
        typedef K size;
        typedef KmerT kmer_type;
        typedef IdT position_type;
        typedef QualT qual_type;

        KmerT kmer;
        IdT id;
        QualT qual;
    };

    /**
     * @class     bliss::io::kmer_index_element
     * @brief     basic struct for storing position and quality of a kmer.
     * @details
     *  K is a typename so that we do not use memory.
     */
    template<typename K, typename KmerT, typename IdT>
    struct kmer_index_element<K, KmerT, IdT, void>
    {
        typedef K size;
        typedef KmerT kmer_type;
        typedef IdT position_type;

        KmerT kmer;
        IdT id;
    };

    /**
     * @class     bliss::io::kmer_index_element
     * @brief     basic struct for storing position and quality of a kmer.
     * @details
     *  K is a typename so that we do not use memory.
     */
    template<typename K, typename KmerT>
    struct kmer_index_element<K, KmerT, void, void>
    {
        typedef K size;
        typedef KmerT kmer_type;

        KmerT kmer;
    };

  } /* namespace io */
} /* namespace bliss */

#endif /* KMER_INDEX_ELEMENT_HPP_ */
