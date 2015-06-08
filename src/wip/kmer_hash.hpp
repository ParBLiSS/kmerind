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
#include "smhasher/MurmurHash3.cpp"
// and farm hash
#include "farmhash.cc"


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

      static constexpr size_t tuplesize = sizeof(uint64_t) / sizeof(WORD_TYPE);
      static constexpr size_t tuples = KmerType::nWords / tuplesize;
      static constexpr size_t leftover = (tuplesize <= 1) ? 0 : (KmerType::nWords % tuplesize) * sizeof(WORD_TYPE);

    public:

      using value_type = KmerType;
      static const unsigned int default_init_value = ::std::min(KmerType::nBits, 64U);
      hash(const unsigned int) {};  // ignored.

      /**
       * @brief function to compute the kmer's hash
       * @details  note that if the kmer uses a word type that is different than uint64_t,
       *        we need to make sure that we use all 64 bits of the hash value..
       *        if the number of bits in kmer is smaller than 64, then we also need to
       *        add an offset.
       *
       *        note that on 32 bit system, this operator does not conform to std::hash declaration.
       *
       * @param kmer
       * @return
       */
      inline uint64_t operator()(value_type const& kmer) const
      {

        // first compute the hash, from 64 bits of input at a time
        uint64_t const * data = reinterpret_cast<uint64_t const*>(kmer.getData());
        uint64_t h = std::hash<uint64_t>()(data[0]);
        uint64_t hp;
        for (size_t i = 1; i < tuples; ++i) {
          hp = std::hash<uint64_t>()(data[i]);
          h ^= (hp << 1);
        }

        // the remainder bits, if any, needs to be dealt with now.
        if (leftover > 0) {
          hp = 0;
          memcpy(&hp, kmer.getData() + tuples * tuplesize, leftover);
          hp = std::hash<uint64_t>()(hp);
          h ^= (hp << 1);
        }

        return h;
      }
  };


}  // namespace std




namespace bliss {

  namespace kmer
  {

    namespace transform {

      // QUESTION:  xor of hash, or hash of xor?.  second is faster.  Also if there is GC-AT imbalance, xor of raw sequence kind of flattens the distribution, so hash input is now more even.

      template <typename KMER>
      struct identity {
          inline KMER operator()(KMER const & x) const {
            return x;
          }
      };

      template <typename KMER>
      struct xor_rev_comp {
          inline KMER operator()(KMER const & x) const {
            return x ^ x.reverse_complement();
          }
          inline KMER operator()(::std::pair<KMER, KMER> const & x) const  {
            return x.first ^ x.second;
          }
      };

      template <typename KMER>
      struct lex_less {
          inline KMER operator()(KMER const & x) const  {
            auto y = x.reverse_complement();
            return (x < y) ? x : y;
          }
          inline KMER operator()(::std::pair<KMER, KMER> const & x) const  {
            return (x.first < x.second) ? x.first : x.second;
          }
      };

      template <typename KMER>
      struct lex_greater {
          inline KMER operator()(KMER const & x) const  {
            auto y = x.reverse_complement();
            return (x > y) ? x : y;
          }
          inline KMER operator()(::std::pair<KMER, KMER> const & x) const  {
            return (x.first > x.second) ? x.first : x.second;
          }
      };


    } // namespace transform


    namespace hash
    {


      /**
       * @brief  Kmer hash, returns the least significant NumBits directly as identity hash.
       * @note   since the number of buckets is not known ahead of time, can't have nbit be a type.  max is 64bits.
       */
      template<typename KMER, bool Prefix = false>
      class cpp_std {
        protected:
          const int64_t shift;
        public:
          static constexpr unsigned int default_init_value = ::std::min(KMER::nBits, 64U);

          cpp_std(const unsigned int prefix_bits = default_init_value) : shift(default_init_value - prefix_bits) {
            if (prefix_bits > default_init_value)  throw ::std::invalid_argument("ERROR: prefix larger than available hash size or kmer size.");
          };

          /// operator to compute hash
          inline size_t operator()(const KMER & kmer) const {
            if (Prefix)
              // if multiword, nBits > 64, then the returned hash will have meaningful msb.  min(nBits, 64) == 64.
              // if single word, nBits <= 64, then min(nBits, 64) = nBits, and shift is the right thing.
              return ::std::hash<KMER>(kmer) >> shift;
            else
              // if multiword, nBits > 64, then the returned hash will have meaningful msb.  min(nBits, 64) == 64.
              // if single word, nBits <= 64, then min(nBits, 64) = nBits, and shift is the right thing.
              return ::std::hash<KMER>(kmer);
          }

      };


      /**
       * @brief  Kmer hash, returns the least significant NumBits directly as identity hash.
       * @note   since the number of buckets is not known ahead of time, can't have nbit be a type
       */
      template <typename KMER, bool Prefix = false>
      class identity {

        protected:
          const unsigned int bits;
        public:
          static const unsigned int default_init_value = ::std::min(KMER::nBits, 64U);

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
          murmur(const unsigned int) {};

          inline uint64_t operator()(KMER const& kmer) const
          {
            // produces 128 bit hash.
            uint64_t h[2];
            // let compiler optimize out all except one of these.
            if (sizeof(void*) == 8)
              MurmurHash3_x64_128(kmer.getData(), nBytes, 42, h);
            else if (sizeof(void*) == 4)
              MurmurHash3_x86_128(kmer.getData(), nBytes, 42, h);
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
          const size_t shift;

        public:
          static const unsigned int default_init_value = 64U;

          farm(const unsigned int prefix_bits = default_init_value) : shift(64U - prefix_bits) {
            if (prefix_bits > default_init_value)  throw ::std::invalid_argument("ERROR: prefix larger than 64 bit.");
          };

          /// operator to compute hash.  64 bit again.
          inline uint64_t operator()(const KMER & kmer) const {
            if (Prefix)
              return ::util::Hash(reinterpret_cast<const char*>(kmer.getData()), nBytes) >> shift;
            else
              return ::util::Hash(reinterpret_cast<const char*>(kmer.getData()), nBytes);
          }

      };

      // usage  apply transform first, then apply hash.
    } // namespace hash
  } // namespace kmer
} // namespace bliss



#endif /* KMER_HASH_HPP_ */
