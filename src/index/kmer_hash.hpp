/**
 * @file    kmer_hash.hpp
 * @ingroup bliss::hash
 * @author  tpan
 * @brief   collections of hash functions defined for kmers.
 * @details support the following:  raw bits directly extracted; std::hash version; murmurhash; and farm hash
 *
 *          assuming the use is in a distributed hash table with total buckets N,  N = p * t * l,
 *          where p is number of processes, t is number of threads, and l is number of local buckets.
 *
 *          then the key is assigned to a bucket via hash(key) % N.
 *            the process assignment is: hash(key) / (t * l)
 *            the thread assignment is:  (hash(key) / l) % t
 *            the local bucket assignment is: hash(key) % l;
 *
 *          there are unlikely to be more than 2^64 local buckets, so we can limit the hash(key) % l to be the lower 64bit of hash(key).
 *          this also means that if the hash key is 64 bits, then no shifting or bit masking is needed, which improves performance for local hashtable lookup.
 *
 *          l is a variable that is decided by the local hash table based on number of entries.
 *          we should instead look at the first ceil(log (p*t)) bits of the hash(key).  let's call this "prefix".
 *
 *          process assignment is then:  prefix / t
 *          thread assignment is then:   prefix % t.
 *
 *          prefix for process assignment can be changed to use ceil(log(p)) bits, let's call this "pre-prefix".
 *
 *          process assignment is then:  pre-prefix % p.
 *
 *          2 functions are sufficient then:  prefix_hash(), and suffix_hash().
 *            we restrict our hash() functions to return 64 bits, and make suffix hash() to be the same as hash().
 *
 *
 *          namespace bliss::hash::kmer has a generic hash function that can work with kmer, kmer xor rev comp, computed or provided.
 *          the generic hash function also allows customization via bliss::hash::kmer::detail::{std,identity,murmur,farm} specializations
 *
 *          as stated above, 2 versions for each specialization: hash() and hash_prefix().  the specialization is especially for identity and murmur hashes,
 *          as murmur hash produces 128 bit value, and identity hash uses the original kmer.
 *
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
#include <type_traits>  // enable_if

#include "common/kmer.hpp"

// includ the murmurhash code.
#ifndef _MURMURHASH3_H_
#include <smhasher/MurmurHash3.cpp>
#endif

// and farm hash
#ifndef FARM_HASH_H_
#include <farmhash/src/farmhash.cc>
#endif

//// Kmer specialization for std::hash
//namespace std {
//  /**
//   * @brief std::hash specialization for kmer.
//   * @note  note that WORD_TYPE may be much smaller than uint64_t.  std::hash does a static_cast internally,
//   *        so high bits are more likely to be 0.
//   *
//   *        solution is to reinterpret the data as uint64_t.
//   */
//  template<unsigned int KMER_SIZE, typename ALPHABET, typename WORD_TYPE>
//  struct hash<::bliss::common::Kmer<KMER_SIZE, ALPHABET, WORD_TYPE> > {
//      using KmerType = ::bliss::common::Kmer<KMER_SIZE, ALPHABET, WORD_TYPE>;
//
//      static constexpr size_t tuples = KmerType::nBytes / sizeof(uint64_t);
//      static constexpr size_t leftover = KmerType::nBytes % sizeof(uint64_t);
//
//      using value_type = KmerType;
//      static const unsigned int default_init_value = KmerType::nBits < 64U ? KmerType::nBits : 64U;
//
//
//      /**
//       * @brief function to compute the kmer's hash
//       * @details  note that if the kmer uses a word type that is different than uint64_t,
//       *        we need to make sure that we use all 64 bits of the hash value..
//       *        if the number of bits in kmer is smaller than 64, then we also need to
//       *        add an offset.
//       *
//       *        note that on 32 bit system, this operator does not conform to std::hash declaration.
//       *
//       * @param kmer
//       * @return
//       */
//      inline uint64_t operator()(value_type const& kmer) const
//      {
//
//        // first compute the hash, from 64 bits of input at a time
//        uint64_t const * data = reinterpret_cast<uint64_t const*>(kmer.getConstData());
//        uint64_t h = std::hash<uint64_t>()(data[0]);
//        uint64_t hp;
//        for (size_t i = 1; i < tuples; ++i) {
//          hp = std::hash<uint64_t>()(data[i]);
//          h ^= (hp << 1);
//        }
//
//        // the remainder bits, if any, needs to be dealt with now.
//        if (leftover > 0) {
//          hp = 0;
//          memcpy(&hp, kmer.getConstData() + tuples * sizeof(uint64_t), leftover);
//          hp = std::hash<uint64_t>()(hp);
//          h ^= (hp << 1);
//        }
//
//        return h;
//      }
//  };
//
//
//}  // namespace std
//



namespace bliss {

  namespace kmer
  {



    namespace hash
    {


      /**
       * @brief  Kmer hash, returns the least significant NumBits directly as identity hash.
       * @note   since the number of buckets is not known ahead of time, can't have nbit be a type.  max is 64bits.
       */
      template<typename KMER, bool Prefix = false>
      class cpp_std {
        protected:
          static constexpr size_t tuples = KMER::nWords * sizeof(typename KMER::KmerWordType) / sizeof(size_t);
          static constexpr size_t leftover = KMER::nWords * sizeof(typename KMER::KmerWordType) % sizeof(size_t);

          int64_t shift;
          ::std::hash<size_t> op;

        public:
          static constexpr unsigned int default_init_value = (KMER::nBits < 64U) ? KMER::nBits : 64U;

          cpp_std(const unsigned int prefix_bits = default_init_value) : shift(default_init_value - prefix_bits) {
            if (prefix_bits > default_init_value)  throw ::std::invalid_argument("ERROR: prefix larger than available hash size or kmer size.");
          };

          /// operator to compute hash
          inline size_t operator()(const KMER & kmer) const {

              size_t const * data = reinterpret_cast<size_t const*>(kmer.getConstData());
              size_t h = 0;

			  // first deal with the remainder bits, if any.
			  if (leftover > 0) {
				memcpy(&h, kmer.getConstData() + tuples * sizeof(size_t), leftover);
				h = op(h);
			  }

              size_t hp;
            // then compute the hash, from 64 bits of input at a time
            for (size_t i = 0; i < tuples; ++i) {
              hp = op(data[i]);
              h ^= (hp << 1);
            }


            if (Prefix)
              // if multiword, nBits > 64, then the returned hash will have meaningful msb.  min(nBits, 64) == 64.
              // if single word, nBits <= 64, then min(nBits, 64) = nBits, and shift is the right thing.
              return h >> shift;
            else
              // if multiword, nBits > 64, then the returned hash will have meaningful msb.  min(nBits, 64) == 64.
              // if single word, nBits <= 64, then min(nBits, 64) = nBits, and shift is the right thing.
              return h;
          }

      };


      /**
       * @brief  Kmer hash, returns the least significant NumBits directly as identity hash.
       * @note   since the number of buckets is not known ahead of time, can't have nbit be a type
       */
      template <typename KMER, bool Prefix = false>
      class identity {

        protected:
          unsigned int bits;
        public:
          static const unsigned int default_init_value = KMER::nBits < 64U ? KMER::nBits : 64U;

          /// constructor
          identity(const unsigned int prefix_bits = default_init_value) : bits(prefix_bits) {
            if ((prefix_bits == 0) || (prefix_bits > default_init_value)) throw ::std::invalid_argument("ERROR: prefix size outside of supported range");
          };

          /// operator to compute hash value
          inline uint64_t operator()(const KMER & kmer) const {
            if (Prefix)
              return kmer.getPrefix(bits);  // get the first bits number of bits from kmer.
            else
              // if nBits < 64U, then it's padded with 0.  just get it.  else we want the entire 64bit.
              return kmer.getSuffix(default_init_value);
          }
      };

      /**
       * @brief Kmer specialization for MurmurHash.  generated hash is 128 bit.
       */
      template <typename KMER, bool Prefix = false>
      class murmur {

        protected:
          static constexpr unsigned int nBytes = (KMER::nBits + 7) / 8;

        public:
          static const unsigned int default_init_value = 64U;

          murmur(const unsigned int prefix_bits = default_init_value) {};

          inline uint64_t operator()(const KMER & kmer) const
          {
            // produces 128 bit hash.
            uint64_t h[2];
            // let compiler optimize out all except one of these.
            if (sizeof(void*) == 8)
              MurmurHash3_x64_128(kmer.getConstData(), nBytes, 42, h);
            else if (sizeof(void*) == 4)
              MurmurHash3_x86_128(kmer.getConstData(), nBytes, 42, h);
            else
              throw ::std::logic_error("ERROR: neither 32 bit nor 64 bit system");

            // use the upper 64 bits.
            if (Prefix)
              return h[1];
            else
              return h[0];
          }

      };

      /**
       * @brief  Kmer hash, returns the least significant NumBits from murmur hash.
       * @note   since the number of buckets is not known ahead of time, can't have nbit be a type
       *
       * @tparam Prefix:  to get prefix,
       */
      template <typename KMER, bool Prefix = false>
      class farm {

        protected:
          static constexpr unsigned int nBytes = (KMER::nBits + 7) / 8;
          size_t shift;

        public:
          static const unsigned int default_init_value = 64U;

          farm(const unsigned int prefix_bits = default_init_value) : shift(64U - prefix_bits) {
            if (prefix_bits > default_init_value)  throw ::std::invalid_argument("ERROR: prefix larger than 64 bit.");
          };

          /// operator to compute hash.  64 bit again.
          inline uint64_t operator()(const KMER & kmer) const {
            if (Prefix)
              return ::util::Hash(reinterpret_cast<const char*>(kmer.getConstData()), nBytes) >> shift;
            else
              return ::util::Hash(reinterpret_cast<const char*>(kmer.getConstData()), nBytes);
          }

      };

      // usage  apply transform first, then apply hash.
    } // namespace hash
  } // namespace kmer
} // namespace bliss



#endif /* KMER_HASH_HPP_ */
