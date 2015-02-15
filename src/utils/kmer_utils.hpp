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

          static_assert(Kmer::bitsPerChar == bliss::AlphabetTraits<Alphabet>::getBitsPerChar(), "Kmer's bits Per Char is different than Alphabet's bits per char.");
          std::stringstream ss;

          using word_type = typename Kmer::KmerWordType;
          word_type w = 0;
          int bit_pos = 0;
          int word_pos = 0;
          int bit_pos_in_word = 0;
          int offset = 0;
          for (int i = Kmer::size - 1; i >= 0; --i)
          {
            // start pos
            bit_pos = i * Kmer::bitsPerChar;
            word_pos = bit_pos / (sizeof(WordType) * 8);
            bit_pos_in_word = bit_pos % (sizeof(WordType) * 8);

            w = (kmer.data[word_pos] >> bit_pos_in_word);

            // now check if char crosses bit boundary.
            if (bit_pos_in_word + Kmer::bitsPerChar > (sizeof(WordType) * 8)) {  // == means that there were just enough bits in the previous word.
              offset = (sizeof(WordType) * 8) - bit_pos_in_word;
              ++word_pos;

              // now need to get the next part.
              if (word_pos < Kmer::nWords) w = w | (kmer.data[word_pos] << offset);

            }
            w = w & getLeastSignificantBitsMask<word_type>(Kmer::bitsPerChar);

            ss << Alphabet::TO_ASCII[w];

          }

          return ss.str();
        }

    };

  } /* namespace io */
} /* namespace bliss */

#endif /* KMERUTILS_HPP_ */
