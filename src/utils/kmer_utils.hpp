/**
 * @file    kmer_utils.hpp
 * @ingroup utils
 * @author  Tony Pan <tpan7@gatech.edu>
 * @brief   Utility functions and functors for kmers
 * @details
 *
 * Copyright (c) 2015 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef KMERUTILS_HPP_
#define KMERUTILS_HPP_

#include "common/kmer.hpp"
#include "common/bit_ops.hpp"

namespace bliss
{
  namespace utils
  {
    /**
     * @class  bliss::utils::KmoleculeToCanonicalKmerFunctor
     * @brief  choose the smaller kmer from a Kmolecule (std::pair of Kmer and its reverse complement)
     *
     */
    template<typename KmerType>
    struct KmoleculeToCanonicalKmerFunctor {
    	/// DEFINE kmolecule type
        typedef std::pair<KmerType, KmerType> KmoleculeType;

        /// operator for choosing the lexigraphically smaller kmer/reverse complement.
        KmerType operator()(const KmoleculeType& input) {
          return (input.first < input.second ? input.first : input.second);
        }
    };

    /**
     * @class    bliss::io::KmerUtils
     * @brief	 collection of kmer related utility functions.
     * @details
     *
     */
    class KmerUtils
    {

      public:
    	/**
    	 * @brief convenience function to convert from kmer to ascii representation.  minimizes shifting.
    	 */
        template<typename Kmer>
        static std::string toASCIIString(const Kmer &kmer)
        {
          using Alphabet = typename Kmer::KmerAlphabet;

          static_assert(Kmer::bitsPerChar == bliss::common::AlphabetTraits<Alphabet>::getBitsPerChar(), "Kmer's bits Per Char is different than Alphabet's bits per char.");

          /* return the char representation of the data array values */
          std::string result;
          result.resize(Kmer::size);
          Kmer cpy(kmer);
          size_t forBitMask = (1 << Kmer::bitsPerChar) - 1;
          for (unsigned int i = 0; i < Kmer::size; ++i)
          {
            result[Kmer::size-i-1] = Alphabet::TO_ASCII[static_cast<size_t>(forBitMask & cpy.data[0])];
            cpy.do_right_shift(Kmer::bitsPerChar);
          }

          return result;

        }

    };

  } /* namespace io */
} /* namespace bliss */

#endif /* KMERUTILS_HPP_ */
