/**
 * @file    kmer_index_element.hpp
 * @ingroup
 * @author  tpan
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
    struct KmerSize
    {
        /**
         * size of kmer
         */
        static constexpr int size = K;
    };


    /**
     * @class     bliss::io::KmerIndexElement
     * @brief     basic struct for storing a kmer (no position or quality).
     * @details   K is a typename so that we do not use memory.
     */
    template<typename K, typename Kmer>
    struct KmerIndexElement
    {
        /**
         * size of kmer
         */
        typedef K SizeType;
        /**
         * kmer type
         */
        typedef Kmer KmerType;

        /**
         * kmer value
         */
        Kmer kmer;
    };


    /**
     * @class     bliss::io::KmerIndexElement
     * @brief     basic struct for storing position of a kmer (no quality score).
     * @details   K is a typename so that we do not use memory.
     *            inherits from KmerIndexElement, which has a kmer member variable.
     */
    template<typename K, typename Kmer, typename Id>
    struct KmerIndexElementWithId : public KmerIndexElement<K, Kmer>
    {
        typedef K SizeType;
        typedef Kmer KmerType;
        typedef Id PositionType;

        /**
         * kmer position/id.
         */
        Id id;
    };


    /**
     * @class     bliss::io::KmerIndexElement
     * @brief     basic struct for storing position and quality of a kmer.
     * @details   K is a typename so that we do not use memory.
     *            inherits from KmerIndexElementWithId, which has kmer and id member variables.
     */
    template<typename K, typename Kmer, typename Id, typename Qual>
    struct KmerIndexElementWithIdAndQuality : public KmerIndexElementWithId<K, Kmer, Id>
    {
        typedef K SizeType;
        typedef Kmer KmerType;
        typedef Id PositionType;
        typedef Qual QualityType;

        /**
         * quality score for the entire kmer
         */
        Qual qual;
    };


  } /* namespace index */
} /* namespace bliss */

#endif /* KMER_INDEX_ELEMENT_HPP_ */
