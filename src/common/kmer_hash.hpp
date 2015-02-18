/**
 * @file    kmer_hash.hpp
 * @ingroup bliss::hash
 * @author  tpan
 * @brief   collections of hash functions defined for kmers.
 * @details support the following:  raw bits directly extracted; std::hash version; murmurhash;
 *
 * Copyright (c) 2015 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef KMER_HASH_HPP_
#define KMER_HASH_HPP_

#include <tuple>  // for hash - std::pair
#include <exception>  // for hash - std::system_error
#include <algorithm>

// includ the murmurhash code.
#include "smhasher/MurmurHash3.cpp"
#include "common/kmer.hpp"


// Kmer specialization for std::hash
namespace std {
  /**
   * @brief std::hash specialization for kmer.
   * @note  note that WORD_TYPE may be much smaller than uint64_t.  std::hash does a static_cast internally,
   *        so high bits are more likely to be 0.
   *
   *        solution is to reinterpret the data as uint64_t.
   */
  template<unsigned int KMER_SIZE, typename ALPHABET, typename WORD_TYPE>
  class hash<::bliss::common::Kmer<KMER_SIZE, ALPHABET, WORD_TYPE> > {
    protected:
      using KmerType = ::bliss::common::Kmer<KMER_SIZE, ALPHABET, WORD_TYPE>;

    public:
      /**
       * @brief function to compute the kmer's hash
       * @details  note that if the kmer uses a word type that is different than uint64_t,
       *        we need to make sure that we use all 64 bits of the hash value..
       *        if the number of bits in kmer is smaller than 64, then we also need to
       *        add an offset.
       *
       * @param kmer
       * @return
       */
      uint64_t operator()(KmerType const& kmer) const
      {
        constexpr int tuplesize = sizeof(uint64_t) / sizeof(WORD_TYPE);
        if (tuplesize == 0) throw std::logic_error("WORD TYPE is larger than 64 bit");

        constexpr int tuples = KmerType::nWords / tuplesize;
        constexpr int leftover = (tuplesize <= 1) ? 0 : KmerType::nWords % tuplesize;

        // first compute the hash, from 64 bits of input at a time
        uint64_t const * data = reinterpret_cast<uint64_t const*>(kmer.getData());
        uint64_t h = std::hash<uint64_t>()(data[0]);
        uint64_t hp;
        for (int i = 1; i < tuples; ++i) {
          hp = std::hash<uint64_t>()(data[i]);
          h ^= (hp << 1);
        }

        // the remainder bits, if any, needs to be dealt with now.
        if (leftover > 0) {
          hp = 0;
          memcpy(&hp, kmer.getData() + tuples * tuplesize, leftover * sizeof(WORD_TYPE));
          hp = std::hash<uint64_t>()(hp);
          h ^= (hp << 1);
        }

        return h;
      }
  };

}  // namespace std

namespace bliss {


  namespace hash {


    namespace std {
      /**
       * @brief  Kmer hash, returns the most significant NumBits from std hash.
       * @note   since the number of buckets is not known ahead of time, can't have nbit be a type
       */
      template <typename KMER>
      class KmerPrefixHash {
        protected:
          const int shift;
          ::std::hash<KMER> hashf;
        public:
          /// constructor
          KmerPrefixHash(const unsigned int nBits = ::std::min(KMER::nBits, 64U)) : shift(::std::min(KMER::nBits, 64U) - nBits) {
            if ((nBits == 0) || (nBits > ::std::min(KMER::nBits, 64U))) throw ::std::invalid_argument("ERROR: does not support more than 64 bits for hash");

          };
          /// operator to compute the hash
          uint64_t operator()(const KMER & kmer) const {
            return hashf(kmer) >> shift;
          };
      };

      /**
       * @brief  Kmer hash, returns the middle NumBits, offset from MSB, directly as identity hash.
       * @note   since the number of buckets is not known ahead of time, can't have nbit be a type
       */
      template<typename KMER>
      class KmerInfixHash {
        protected:
          const int shift;
          const uint64_t mask;
          ::std::hash<KMER> hashf;
        public:
          /// constructor
          KmerInfixHash(const unsigned int nBits = ::std::min(KMER::nBits, 64U), const unsigned int _offset = 0) :
              shift(::std::min(KMER::nBits, 64U) - nBits - _offset), mask(::std::numeric_limits<uint64_t>::max() >> (::std::min(KMER::nBits, 64U) - nBits)) {
            if ((nBits == 0) || (nBits > ::std::min(KMER::nBits, 64U))) throw ::std::invalid_argument("ERROR: does not support more than 64 bits for hash");
            if ((nBits + _offset) > ::std::min(KMER::nBits, 64U) ) throw ::std::invalid_argument("ERROR: nBits and _offset together cannot exceed 64 bits");

          };
          /// operator to compute the hash
          uint64_t operator()(const KMER & kmer) const {
            return (hashf(kmer) >> shift) & mask;
          }
      };

      /**
       * @brief  Kmer hash, returns the least significant NumBits directly as identity hash.
       * @note   since the number of buckets is not known ahead of time, can't have nbit be a type
       */
      template <typename KMER>
      class KmerSuffixHash {
        protected:
          const uint64_t mask;
          ::std::hash<KMER> hashf;
        public:
          /// constructor
          KmerSuffixHash(const unsigned int nBits = ::std::min(KMER::nBits, 64U)) :
              mask(::std::numeric_limits<uint64_t>::max() >> (::std::min(KMER::nBits, 64U) - nBits)) {
            if ((nBits == 0) || (nBits > ::std::min(KMER::nBits, 64U))) throw ::std::invalid_argument("ERROR: does not support more than 64 bits for hash");
          };
          /// operator to compute the hash
          uint64_t operator()(const KMER & kmer) const {
            return hashf(kmer) & mask;
          }
      };

    }



    namespace identity {
      /**
       * @brief  Kmer hash, returns the most significant NumBits directly as identity hash.
       * @note   since the number of buckets is not known ahead of time, can't have nbit be a type
       */
      template <typename KMER>
      class KmerPrefixHash {
        protected:
          const unsigned int NumBits;
        public:
          /// constructor
          KmerPrefixHash(const unsigned int nBits = ::std::min(KMER::nBits, 64U)) : NumBits(nBits) {
            if ((nBits == 0) || (nBits > ::std::min(KMER::nBits, 64U))) throw ::std::invalid_argument("ERROR: does not support more than 64 bits for hash");

          };
          /// operator to compute hash value
          uint64_t operator()(const KMER & kmer) const {
            return kmer.getPrefix(NumBits);
          };
      };

      /**
       * @brief  Kmer hash, returns the middle NumBits, offset from MSB, directly as identity hash.
       * @note   since the number of buckets is not known ahead of time, can't have nbit be a type
       */
      template<typename KMER>
      struct KmerInfixHash {
        protected:
          const unsigned int NumBits;
          const unsigned int offset;
        public:
          /// constructor
          KmerInfixHash(const unsigned int nBits = ::std::min(KMER::nBits, 64U), const unsigned int _offset = 0) : NumBits(nBits), offset(_offset) {
            if ((nBits == 0) || (nBits > ::std::min(KMER::nBits, 64U))) throw ::std::invalid_argument("ERROR: does not support more than 64 bits for hash");
            if ((nBits + _offset) > KMER::nBits ) throw ::std::invalid_argument("ERROR: nBits and _offset together cannot exceed total number of bits for the kmer bits");
          };
          /// operator to compute hash value
          uint64_t operator()(const KMER & kmer) const {
            return kmer.getInfix(NumBits, offset);
          }
      };

      /**
       * @brief  Kmer hash, returns the least significant NumBits directly as identity hash.
       * @note   since the number of buckets is not known ahead of time, can't have nbit be a type
       */
      template <typename KMER>
      struct KmerSuffixHash {
        protected:
          const unsigned int NumBits;
        public:
          /// constructor
          KmerSuffixHash(const unsigned int nBits = ::std::min(KMER::nBits, 64U)) : NumBits(nBits) {
            if ((nBits == 0) || (nBits > ::std::min(KMER::nBits, 64U))) throw ::std::invalid_argument("ERROR: does not support more than 64 bits for hash");

          };
          /// operator to compute hash value
          uint64_t operator()(const KMER & kmer) const {
            return kmer.getSuffix(NumBits);
          }
      };


    }


    namespace murmur {

      template<typename T>
      class hash;

      /**
       * @brief Kmer specialization for MurmurHash.
       */
      template<unsigned int KMER_SIZE, typename ALPHABET, typename WORD_TYPE>
      class hash<bliss::common::Kmer<KMER_SIZE, ALPHABET, WORD_TYPE> > {
        protected:
          using KmerType = bliss::common::Kmer<KMER_SIZE, ALPHABET, WORD_TYPE>;

        public:
          ::std::pair<uint64_t, uint64_t> operator()(KmerType const& kmer) const
          {
            uint64_t h[2];
            // let compiler optimize out all except one of these.
            if (sizeof (void*) == 8)
              MurmurHash3_x64_128(kmer.getData(), sizeof(WORD_TYPE) * KmerType::nWords, 42, h);
            else if (sizeof(void*) == 4)
              MurmurHash3_x86_128(kmer.getData(), sizeof(WORD_TYPE) * KmerType::nWords, 42, h);
            else
              throw ::std::logic_error("ERROR: neither 32 bit nor 64 bit system");

            return ::std::make_pair(h[0], h[1]);
          }
      };


      /**
       * @brief  Kmer hash, returns the most significant NumBits from murmur hash.
       * @note   since the number of buckets is not known ahead of time, can't have nbit be a type
       */
      template <typename KMER>
      class KmerPrefixHash {
        protected:
          const int shift;
          ::bliss::hash::murmur::hash<KMER> hashf;
        public:
          /// constructor
          KmerPrefixHash(const unsigned int nBits = 64) : shift(sizeof(uint64_t) * 8 - nBits) {
            if ((nBits == 0) || (nBits > sizeof(uint64_t) * 8)) throw ::std::invalid_argument("ERROR: does not support more than 64 bits for hash");
          };
          /// operator to compute actual hash
          uint64_t operator()(const KMER & kmer) const {
            uint64_t hl, hh;
            ::std::tie(hl, hh) = hashf(kmer);
            return hh >> shift;
          };
      };

      /**
       * @brief  Kmer hash, returns the middle NumBits, offset from MSB, from murmur hash.
       * @note   since the number of buckets is not known ahead of time, can't have nbit be a type
       */
      template<typename KMER>
      class KmerInfixHash {
        protected:
          const int shift1;
          const int shift2;
          const uint64_t mask;

          /// function pointer so we can choose the logic based on nBits and offset.
          uint64_t (KmerInfixHash::*extract) (uint64_t& hl, uint64_t& hh) const;
          ::bliss::hash::murmur::hash<KMER> hashf;
        public:
          // constructor
          KmerInfixHash(const unsigned int nBits = 64, const unsigned int _offset = 0) :
              shift1(sizeof(uint64_t) * 8 * 2 - nBits - _offset), shift2(sizeof(uint64_t) * 8 - nBits - _offset),
              mask(::std::numeric_limits<uint64_t>::max() >> (sizeof(uint64_t) * 8 - nBits)){

            if ((nBits == 0) || (nBits > sizeof(uint64_t) * 8)) throw ::std::invalid_argument("ERROR: does not support more than 64 bits for hash");
            if ((nBits + _offset) > (sizeof(uint64_t) * 8 * 2) ) throw ::std::invalid_argument("ERROR: nBits and _offset together cannot exceed 128 bits");

            // make some algo choices.
            if (_offset >= sizeof(uint64_t) * 8) { // in lower hash
              extract = &KmerInfixHash::extract1;
            } else if (_offset + nBits <  sizeof(uint64_t) * 8) { // in high hash
              extract = &KmerInfixHash::extract2;
            } else { // straddling
              extract = &KmerInfixHash::extract3;
            }

          };

          /// operator to compute hash
          uint64_t operator()(const KMER & kmer) const {
            uint64_t hl, hh;
            ::std::tie(hl, hh) = hashf(kmer);
            return (this->*extract)(hl, hh) & mask;
          }
        protected:
          //=== specializations based on whether the bits to extract spans the 64 bit boundary or not.
          uint64_t extract1(uint64_t& hl, uint64_t& hh) const {
            return (hl >> shift1);
          }
          uint64_t extract2(uint64_t& hl, uint64_t& hh) const {
            return (hh >> shift2);
          }
          uint64_t extract3(uint64_t& hl, uint64_t& hh) const {
            return (hl >> shift1) | (hh << (-shift2));
          }

      };

      /**
       * @brief  Kmer hash, returns the least significant NumBits from murmur hash.
       * @note   since the number of buckets is not known ahead of time, can't have nbit be a type
       */
      template <typename KMER>
      class KmerSuffixHash {
        protected:

          const uint64_t mask;
          ::bliss::hash::murmur::hash<KMER> hashf;
        public:
          /// constructor
          KmerSuffixHash(const unsigned int nBits = 64) :
              mask(::std::numeric_limits<uint64_t>::max() >> (sizeof(uint64_t) * 8 - nBits)){
            if ((nBits == 0) || (nBits > sizeof(uint64_t) * 8)) throw ::std::invalid_argument("ERROR: does not support more than 64 bits for hash");
          };

          /// operator to compute hash
          uint64_t operator()(const KMER & kmer) const {
            uint64_t hl, hh;
            ::std::tie(hl, hh) = hashf(kmer);
            return hl & mask;
          }
      };


    } // namespace murmur
  } // namespace hash
} // namespace bliss



#endif /* KMER_HASH_HPP_ */
