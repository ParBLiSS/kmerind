/*
 * Copyright 2015 Georgia Institute of Technology
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

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
 */
#ifndef KMER_HASH_HPP_
#define KMER_HASH_HPP_

#include <tuple>  // for hash - std::pair
#include <exception>  // for hash - std::system_error
#include <algorithm>
#include <type_traits>  // enable_if

#include "common/alphabets.hpp"
#include "common/kmer.hpp"

#include "utils/transform_utils.hpp"

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
//        uint64_t const * data = reinterpret_cast<uint64_t const*>(kmer.getData());
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
//          memcpy(&hp, reinterpret_cast<uint64_t *>(kmer.getData()) + tuples, leftover);
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
          static constexpr uint8_t batch_size = 1;

          // if nBits is more than 32, then default prefix should be 32.
          static constexpr unsigned int default_init_value = 32U;

          cpp_std(const unsigned int prefix_bits = default_init_value) : shift(std::min(KMER::nBits, 64U) - std::min(prefix_bits, KMER::nBits)) {};

          /// operator to compute hash
          inline size_t operator()(const KMER & kmer) const {

              size_t const * data = reinterpret_cast<size_t const*>(kmer.getData());
              size_t h = 0;

              // first deal with the remainder bits, if any.
              if (leftover > 0) {
                memcpy(&h, data + tuples, leftover);
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
              return h;  // suffix.  just return the whole thing.
          }

      };
      template<typename KMER, bool Prefix>
      constexpr uint8_t cpp_std<KMER, Prefix>::batch_size;

      /**
       * @brief  Kmer hash, returns the least significant NumBits directly as identity hash.
       * @note   since the number of buckets is not known ahead of time, can't have nbit be a type
       */
      template <typename KMER, bool Prefix = false>
      class identity {

        protected:
          unsigned int bits;
        public:
          static constexpr uint8_t batch_size = 1;

          static constexpr unsigned int default_init_value = 24U;
          static constexpr unsigned int suffix_bits = (KMER::nBits > 64U) ? 64U : KMER::nBits;

          /// constructor
          identity(const unsigned int prefix_bits = default_init_value) : bits(std::min(KMER::nBits, prefix_bits)) {
          };

          /// operator to compute hash value
          inline uint64_t operator()(const KMER & kmer) const {
            if (Prefix)
              return kmer.getPrefix(bits);  // get the first 32 bits from kmer, or maybe the whole thing if nBits is less than 32.
            else
              // get the whole thing
              return kmer.getSuffix(suffix_bits);
          }
      };
      template<typename KMER, bool Prefix>
      constexpr uint8_t identity<KMER, Prefix>::batch_size;

      /**
       * @brief Kmer specialization for MurmurHash.  generated hash is 128 bit.
       *
       * TODO: move KMER type template param to operator.
       * TODO: change h to member variable.
       */
      template <typename KMER, bool Prefix = false>
      class murmur {


        protected:
          static constexpr unsigned int nBytes = (KMER::nBits + 7) / 8;
          uint32_t seed;

        public:
          static constexpr uint8_t batch_size = 1;

          static const unsigned int default_init_value = 24U;  // allow 16M processors.  but it's ignored here.

          murmur(const unsigned int prefix_bits = default_init_value, uint32_t const & _seed = 42 ) : seed(_seed) {};

          inline uint64_t operator()(const KMER & kmer) const
          {
            // produces 128 bit hash.
            uint64_t h[2];
            // let compiler optimize out all except one of these.
            if (sizeof(void*) == 8)
              MurmurHash3_x64_128(kmer.getData(), nBytes, seed, h);
            else if (sizeof(void*) == 4)
              MurmurHash3_x86_128(kmer.getData(), nBytes, seed, h);
            else
              throw ::std::logic_error("ERROR: neither 32 bit nor 64 bit system");

            // use the upper 64 bits.
            if (Prefix)
              return h[1];
            else
              return h[0];
          }

      };
      template<typename KMER, bool Prefix>
      constexpr uint8_t murmur<KMER, Prefix>::batch_size;

      /**
       * @brief  Kmer hash, returns the least significant NumBits from murmur hash.
       * @note   since the number of buckets is not known ahead of time, can't have nbit be a type
       *		IMPORTANT: upper and lower halves of the hash values are not independent, when farm hash is used with a given seed.
       *		therefore, different seed must be applied.  to adhere to existing api, we generate a different seed for "prefix" version.
       * @tparam Prefix:
       */
      template <typename KMER, bool Prefix = false>
      class farm {

        protected:
          static constexpr unsigned int nBytes = (KMER::nBits + 7) / 8;
          size_t shift;
          uint32_t seed;

        public:
          static constexpr uint8_t batch_size = 1;

          static const unsigned int default_init_value = 24U;   // this allows 16M processors.

          farm(const unsigned int prefix_bits = default_init_value, uint32_t const & _seed = 42 ) : shift(64U - std::min(prefix_bits, 64U)), seed(_seed)  {
          };

          /// operator to compute hash.  64 bit again.
          inline uint64_t operator()(const KMER & kmer) const {
            if (Prefix)
              return ::util::Hash64WithSeed(reinterpret_cast<const char*>(kmer.getData()), nBytes, (seed << 1) - 1);
            else
              return ::util::Hash64WithSeed(reinterpret_cast<const char*>(kmer.getData()), nBytes, seed);
          }

      };
      template<typename KMER, bool Prefix>
      constexpr uint8_t farm<KMER, Prefix>::batch_size;


      namespace sparsehash {
      	  //  ===============
      	  //  Sparse hash specific, kmer related stuff



      	  // sparsehash empty keys and delete keys.  these are kmers that NEVER appears in input data.
      	  // NOTE:  single strand and bimolecule have kmers that can be on either strand.  canonical has a more restricted set.
      	  // so for DNA and DNA16, full k-mer (no padding bits), use a key that's compatible with canonical.
      	  // for the other strand cases, use the bitwise complement of that key (in implementation directly.)
      	  // 			not-full			full, Single stranded or bimolecule	     full, Canonical
      	  //   DNA		use high bits		split into high and low hashtables		 no TTTTTTTTTT, or TTTTTTTTG - chosen because revcomp is rev + neg
      	  //   DNA16	use high bits		split into high and low hashtables       no NNNNNNNNN(~A), or NNNNNNNN(~A~C) - chosen because revcomp is bit reverse
      	  //   DNA5		use unused values : 2 (010) and 5 (101).  set at low bit
      	  // note that full word kmers (even kmers) should be infrequent.


      	  /// base class to store some convenience functions
      	  template <typename Kmer>
      	  class special_keys_base {

      	  public:
      		  inline Kmer invert(Kmer const & x) {
      			  Kmer out(false);
      		      for (size_t i = 0; i < Kmer::nWords; ++i) {
      		        out.getDataRef()[i] = ~(x.getData()[i]);
      		      }
      		      return out;
      		  }

      	  protected:
      		  inline Kmer get_min() {  // 00000000000000000
      			  return Kmer();
      		  }

      		  inline Kmer get_max() {   // 11111111111111111
      			  Kmer out(false);
      			  std::memset(out.getDataRef(), 0xFF, Kmer::nWords * sizeof(typename Kmer::KmerWordType));
      			  return out;
      		  }

      		  inline Kmer get_middle() {  // 10000000000000000
      			  Kmer out;
      			  out.getDataRef()[Kmer::nWords - 1] = static_cast<typename Kmer::KmerWordType>(~(::std::numeric_limits<typename Kmer::KmerWordType>::max() >> 1));
      			  return out;
      		  }

      	  };

      	  template <typename Kmer, bool canonical = false>
      	  class special_keys;

      	  template <unsigned int K, typename WordType, bool canonical>   // full specialization, using a templated class as parameter.
      	  class special_keys< ::bliss::common::Kmer<K, ::bliss::common::DNA6, WordType>, canonical > :
		  	  public special_keys_base<::bliss::common::Kmer<K, ::bliss::common::DNA6, WordType> > {

		  	  protected:
      		  	  using Kmer = ::bliss::common::Kmer<K, ::bliss::common::DNA6, WordType>;
      		  	  using Base = ::bliss::kmer::hash::sparsehash::special_keys_base<Kmer >;

		  	  public:

				  /// kmer empty key for DNA5  000000000010  or 0000000000101  - for use as empty and deleted.
				  inline Kmer generate(uint8_t id = 0) {
					  Kmer em(true); // create an empty one
					  em.getDataRef()[0] = (id == 0) ? 0x010 : 0x101;
					  return em;
				  }

				  using Base::invert;  // normal invert

				  inline Kmer get_splitter() {   // 11111111111111
					  return this->get_max();
				  }

				  static constexpr bool need_to_split = false;

      	  };

      	  template <unsigned int K, typename WordType, bool canonical>   // full specialization, using a templated class as parameter.
      	  class special_keys <::bliss::common::Kmer<K, ::bliss::common::DNA, WordType>, canonical > :
		  	  public special_keys_base<::bliss::common::Kmer<K, ::bliss::common::DNA, WordType> > {

		  	  protected:
      		  	  using Kmer = ::bliss::common::Kmer<K, ::bliss::common::DNA, WordType>;
      		  	  using Base = ::bliss::kmer::hash::sparsehash::special_keys_base<Kmer >;

		  	  public:
      		  	  static constexpr bool full = ((Kmer::nWords * sizeof(typename Kmer::KmerWordType) * 8 - Kmer::nBits) <= 1);
				  static constexpr bool need_to_split = !canonical && full;

          		  /// kmer empty key for DNA, unfull kmer 10000000000000, or 110000000000
          		  template <typename KM = Kmer, bool s = full, typename std::enable_if<!s, int>::type = 0>
          		  inline KM generate(uint8_t id = 0) {
          			  assert(id < 2);
          			  KM em(true); // create an empty one
          			  em.getDataRef()[KM::nWords - 1] = static_cast<typename KM::KmerWordType>(~(::std::numeric_limits<typename KM::KmerWordType>::max() >> (id + 1)));
          			  return em;
          		  }

          		  /// kmer empty key for DNA, full kmer 111111111111111, or 111111111111110
          		  template <typename KM = Kmer, bool s = full, typename std::enable_if<s, int>::type = 0>
          		  inline KM generate(uint8_t id = 0) {
          			  assert(id < 2);
          			  KM em; // create an empty one
          			  em.getDataRef()[0] = (id == 0) ? ::std::numeric_limits<typename KM::KmerWordType>::max() :
          					  ::std::numeric_limits<typename KM::KmerWordType>::max() ^ 0x1;
          			  for (size_t i = 1; i < KM::nWords; ++i) {
          				  em.getDataRef()[i] = ::std::numeric_limits<typename KM::KmerWordType>::max();
          			  }
          			  return em;
          		  }

				  using Base::invert;  // normal invert

				  inline Kmer get_splitter() {   // 10000000000000 for full kmer, 1111111111 for unfull
					  return (need_to_split) ? this->get_middle() : this->get_max();
				  }

      	  };


      	  template <unsigned int K, typename WordType, bool canonical>   // full specialization, using a templated class as parameter.
      	  class special_keys< ::bliss::common::Kmer<K, ::bliss::common::DNA16, WordType>, canonical > :
		  	  public special_keys_base<::bliss::common::Kmer<K, ::bliss::common::DNA16, WordType> > {

		  	  protected:
      		  	  using Kmer = ::bliss::common::Kmer<K, ::bliss::common::DNA16, WordType>;
      		  	  using Base = ::bliss::kmer::hash::sparsehash::special_keys_base<Kmer >;

		  	  public:
      		  	  static constexpr bool full = ((Kmer::nWords * sizeof(typename Kmer::KmerWordType) * 8 - Kmer::nBits) <= 1);
				  static constexpr bool need_to_split = !canonical && full;

          		  /// kmer empty key for DNA, unfull kmer 10000000000000, or 110000000000
          		  template <typename KM = Kmer, bool s = full, typename std::enable_if<!s, int>::type = 0>
          		  inline KM generate(uint8_t id = 0) {
          			  assert(id < 2);
          			  KM em(true); // create an empty one
          			  em.getDataRef()[KM::nWords - 1] = static_cast<typename KM::KmerWordType>(~(::std::numeric_limits<typename KM::KmerWordType>::max() >> (id + 1)));
          			  return em;
          		  }

          		  /// kmer empty key for DNA, full kmer 111111111111110, or 1111111111111100
          		  template <typename KM = Kmer, bool s = full, typename std::enable_if<s, int>::type = 0>
          		  inline KM generate(uint8_t id = 0) {
          			  assert(id < 2);
          			  KM em; // create an empty one
          			  em.getDataRef()[0] = ::std::numeric_limits<typename KM::KmerWordType>::max() ^ (id == 0 ? 0x1 : 0x11);
          			  for (size_t i = 1; i < KM::nWords; ++i) {
          				  em.getDataRef()[i] = ::std::numeric_limits<typename KM::KmerWordType>::max();
          			  }
          			  return em;
          		  }


				  using Base::invert;  // normal invert

				  inline Kmer get_splitter() {   // 10000000000000 if full, 11111111111 if not full
					  return (need_to_split) ? this->get_middle() : this->get_max();
				  }

      	  };

//
//      	  template <typename KMER>
//      	  struct empty_key {
//
//      		  /// kmer empty key for DNA5  000000000010
//      		  template <typename KM = KMER,
//      				  typename std::enable_if<
//					  ::std::is_same<typename KM::KmerAlphabet, ::bliss::common::DNA5>::value, int>::type = 0>
//      		  static inline KMER generate() {
//      			  KMER em(true); // create an empty one
//      			  em.getDataRef()[0] = 0x010;
//      			  return em;
//      		  }
//
//      		  /// kmer empty key for DNA, unfull kmer 10000000
//      		  template <typename KM = KMER,
//      				  typename std::enable_if<
//					  (::std::is_same<typename KM::KmerAlphabet, ::bliss::common::DNA>::value ||
//							  ::std::is_same<typename KM::KmerAlphabet, ::bliss::common::DNA16>::value)	  &&
//					  (KM::nBits < (KM::nWords * sizeof(typename KM::KmerWordType) * 8)), int>::type = 0>
//      		  static inline KMER generate() {
//      			  KMER em(true); // create an empty one
//      			  em.getDataRef()[KM::nWords - 1] = static_cast<typename KM::KmerWordType>(~(::std::numeric_limits<typename KM::KmerWordType>::max() >> 1));
//      			  return em;
//      		  }
//
//      		  /// kmer empty key for DNA, full kmer.  11111111111111
//      		  template <typename KM = KMER,
//      				  typename std::enable_if<
//					  ::std::is_same<typename KM::KmerAlphabet, ::bliss::common::DNA>::value &&
//					  (KM::nBits == (KM::nWords * sizeof(typename KM::KmerWordType) * 8)), int>::type = 0>
//      		  static inline KMER generate() {
//      			  KMER em; // create an empty one
//      			  for (size_t i = 0; i < KM::nWords; ++i) {
//      				  em.getDataRef()[i] = ::std::numeric_limits<typename KM::KmerWordType>::max();
//      			  }
//      			  return em;
//      		  }
//
//      		  /// kmer empty key for DNA16, full kmer.  111111111111110
//      		  template <typename KM = KMER,
//      				  typename std::enable_if<
//					  ::std::is_same<typename KM::KmerAlphabet, ::bliss::common::DNA16>::value &&
//					  (KM::nBits == (KM::nWords * sizeof(typename KM::KmerWordType) * 8)), int>::type = 0>
//      		  static inline KMER generate() {
//      			  KMER em; // create an empty one
//      			  for (size_t i = 1; i < KM::nWords; ++i) {
//      				  em.getDataRef()[i] = ::std::numeric_limits<typename KM::KmerWordType>::max();
//      			  }
//      			  em.getDataRef()[0] = ::std::numeric_limits<typename KM::KmerWordType>::max() ^ 0x1;  // all 1s is its own complement and maps to N, so is a possible value
//      			  return em;
//      		  }
//
//      	  };
//
//      	  template <typename KMER>
//      	  struct deleted_key {
//
//      		  /// kmer empty key for DNA5  0000000000101
//      		  template <typename KM = KMER,
//      				  typename std::enable_if<
//					  ::std::is_same<typename KM::KmerAlphabet, ::bliss::common::DNA5>::value, int>::type = 0>
//      		  static inline KMER generate() {
//      			  KMER em(true); // create an empty one
//      			  em.getDataRef()[0] = 0x101;
//      			  return em;
//      		  }
//
//      		  /// kmer empty key for DNA or DNA16, unfull kmer 110000000
//      		  template <typename KM = KMER,
//      				  typename std::enable_if<
//					  (::std::is_same<typename KM::KmerAlphabet, ::bliss::common::DNA>::value ||
//					  	::std::is_same<typename KM::KmerAlphabet, ::bliss::common::DNA16>::value) &&
//					  ((KM::nWords * sizeof(typename KM::KmerWordType) * 8 - KM::nBits) > 1), int>::type = 0>
//      		  static inline KMER generate() {
//      			  KMER em(true); // create an empty one
//      			  em.getDataRef()[KM::nWords - 1] = static_cast<typename KM::KmerWordType>(~(::std::numeric_limits<typename KM::KmerWordType>::max() >> 2));
//      			  return em;
//      		  }
//
//      		  /// kmer empty key for DNA, full kmer.  11111111111110
//      		  template <typename KM = KMER,
//      				  typename std::enable_if<
//					  ::std::is_same<typename KM::KmerAlphabet, ::bliss::common::DNA>::value  &&
//					  (KM::nBits == (KM::nWords * sizeof(typename KM::KmerWordType) * 8)), int>::type = 0>
//      		  static inline KMER generate() {
//      			  KMER em; // create an empty one
//      			  for (size_t i = 1; i < KM::nWords; ++i) {
//      				  em.getDataRef()[i] = ::std::numeric_limits<typename KM::KmerWordType>::max();
//      			  }
//      			 em.getDataRef()[0] = ::std::numeric_limits<typename KM::KmerWordType>::max() ^ 0x1;
//      			  return em;
//      		  }
//
//
//      		  /// kmer empty key for DNA16, full kmer.  111111111111100
//      		  template <typename KM = KMER,
//      				  typename std::enable_if<
//					  ::std::is_same<typename KM::KmerAlphabet, ::bliss::common::DNA16>::value &&
//					  (KM::nBits == (KM::nWords * sizeof(typename KM::KmerWordType) * 8)), int>::type = 0>
//      		  static inline KMER generate() {
//      			  KMER em; // create an empty one
//      			  for (size_t i = 1; i < KM::nWords; ++i) {
//      				  em.getDataRef()[i] = ::std::numeric_limits<typename KM::KmerWordType>::max();
//      			  }
//      			 em.getDataRef()[0] = ::std::numeric_limits<typename KM::KmerWordType>::max() ^ 0x11;
//      			  return em;
//      		  }
//
//      	  };
//
//      	  template <typename KMER>
//      	  struct split_key {
//
//      		  static constexpr bool need_to_split =
//      				  (::std::is_same<typename KMER::KmerAlphabet, ::bliss::common::DNA>::value ||
//  						::std::is_same<typename KMER::KmerAlphabet, ::bliss::common::DNA16>::value)  &&
//  					  ((KMER::nWords * sizeof(typename KMER::KmerWordType) * 8 - KMER::nBits) <= 1);
//
//      		  /// kmer empty key for DNA5  11111111111111111111
//      		  template <typename KM = KMER,
//      				  typename std::enable_if<
//					  ::std::is_same<typename KM::KmerAlphabet, ::bliss::common::DNA5>::value, int>::type = 0>
//      		  static inline KMER generate() {
//      			  KMER em; // create an empty one
//      			  for (size_t i = 0; i < KM::nWords; ++i) {
//      				  em.getDataRef()[i] = ::std::numeric_limits<typename KM::KmerWordType>::max();
//      			  }
//      			  return em;
//      		  }
//
//      		  /// kmer empty key for DNA or DNA16, unfull kmer 11111111111111111
//      		  template <typename KM = KMER,
//      				  typename std::enable_if<
//					  (::std::is_same<typename KM::KmerAlphabet, ::bliss::common::DNA>::value ||
//					  	::std::is_same<typename KM::KmerAlphabet, ::bliss::common::DNA16>::value) &&
//					  ((KM::nWords * sizeof(typename KM::KmerWordType) * 8 - KM::nBits) > 1), int>::type = 0>
//      		  static inline KMER generate() {
//      			  KMER em; // create an empty one
//      			  for (size_t i = 0; i < KM::nWords; ++i) {
//      				  em.getDataRef()[i] = ::std::numeric_limits<typename KM::KmerWordType>::max();
//      			  }
//      			  return em;
//      		  }
//
//      		  /// kmer empty key for DNA, full kmer.  10000000000000
//      		  template <typename KM = KMER,
//      				  typename std::enable_if<
//					  (::std::is_same<typename KM::KmerAlphabet, ::bliss::common::DNA>::value ||
//						::std::is_same<typename KM::KmerAlphabet, ::bliss::common::DNA16>::value)  &&
//					  ((KM::nWords * sizeof(typename KM::KmerWordType) * 8 - KM::nBits) <= 1), int>::type = 0>
//      		  static inline KMER generate() {
//      			  KMER em(true); // create an empty one
//      			  em.getDataRef()[KM::nWords - 1] = static_cast<typename KM::KmerWordType>(~(::std::numeric_limits<typename KM::KmerWordType>::max() >> 1));
//      			  return em;
//      		  }
//
//      	  };


      	  // this covers the Canonical case when kmer is full and DNA is 2bit or 4 bit, which allows lex_less
      	  // for less_greater, we can use the negation of these keys for DNA and DNA16, full kmer.

//
//
//      	  /// special equal_to just for sparsehash.  do not transform keys (MSB not 00) for partial DNA and DNA16 kmers
//      	  template <typename Kmer, template <typename> class Transform, bool forLower>
//      	  struct transformed_equal_to {
//      	      std::equal_to<Kmer> comp;
//      	      Transform<Kmer> trans;
//
//      	      // 110000000
//      	      static constexpr typename Kmer::KmerWordType highmask =
//      	          static_cast<typename Kmer::KmerWordType>(~(::std::numeric_limits<typename Kmer::KmerWordType>::max() >> 2));
//
//      	      transformed_equal_to(std::equal_to<Kmer> const & _cmp = std::equal_to<Kmer>(),
//      	          Transform<Kmer> const &_trans = Transform<Kmer>()) : comp(_cmp), trans(_trans) {};
//
//      	      // TODO: short circuit more ...
//
//      	      /// single or canonical (no transform), then compare raw values.  all k-mers are compared.
//      	      template <typename KM = Kmer,
//      	          typename std::enable_if<::std::is_same<Transform<KM>, ::bliss::transform::identity<KM> >::value, int>::type = 0>
//      	      inline bool operator()(Kmer const & x, Kmer const & y) const {
//      	        return comp(x, y);
//      	      }
//
//      	      /// DNA or DNA16, partial k-mer (using unused bits), bimolecule (not identity transform).  keys are outside of input k-mer space
//      	      template <typename KM = Kmer,
//                  typename std::enable_if<
//                  (!::std::is_same<Transform<Kmer>, ::bliss::transform::identity<Kmer> >::value) &&
//                      (::std::is_same<typename KM::KmerAlphabet, ::bliss::common::DNA>::value ||
//                      ::std::is_same<typename KM::KmerAlphabet, ::bliss::common::DNA16>::value) &&
//                       ((KM::nWords * sizeof(typename KM::KmerWordType) * 8 - KM::nBits) > 1), int>::type = 0>
//      	      inline bool operator()(Kmer const & x, Kmer const & y) const {
//      	        typename KM::KmerWordType x_key_bits = (x.getData()[KM::nWords - 1] & highmask);
//      	        typename KM::KmerWordType y_key_bits = (y.getData()[KM::nWords - 1] & highmask);
//      	        return (x_key_bits != y_key_bits) ?
//      	            false :   // one key, one value.
//      	            (x_key_bits > 0) ? comp(x, y) :   // compare keys
//      	                comp(trans(x), trans(y));     // compare values
//      	      }
//
//              /// comparator for DNA6, bimolecule (not identity transform).  keys are outside of input k-mer space
//              template <typename KM = Kmer,
//                  typename std::enable_if<
//                  (!::std::is_same<Transform<Kmer>, ::bliss::transform::identity<Kmer> >::value) &&
//                  ::std::is_same<typename KM::KmerAlphabet, ::bliss::common::DNA6>::value, int>::type = 0>
//              inline bool operator()(Kmer const & x, Kmer const & y) const {
//                typename KM::KmerWordType x_key_bits = (x.getData()[0] & 0x111);
//                typename KM::KmerWordType y_key_bits = (y.getData()[0] & 0x111);
//                bool x_is_key = (x_key_bits == 0x101) || (x_key_bits == 0x010);
//                bool y_is_key = (y_key_bits == 0x101) || (y_key_bits == 0x010);
//                return (x_is_key != y_is_key) ? false :   // one key, one value
//                  x_is_key ? comp(x, y) : comp(trans(x), trans(y));
//              }
//
//
//      	      /// comparator for full DNA/DNA16 k-mers.  keys are inside input k-mer space.  compare explitictly.
//              template <typename KM = Kmer,
//                  typename std::enable_if<
//                  (!::std::is_same<Transform<Kmer>, ::bliss::transform::identity<Kmer> >::value) &&
//                      (::std::is_same<typename KM::KmerAlphabet, ::bliss::common::DNA>::value ||
//                      ::std::is_same<typename KM::KmerAlphabet, ::bliss::common::DNA16>::value) &&
//                       ((KM::nWords * sizeof(typename KM::KmerWordType) * 8 - KM::nBits) <= 1), int>::type = 0>
//              inline bool operator()(Kmer const & x, Kmer const & y) const {
//
//                return comp(trans(x), trans(y));
//              }
//
//
//      	      template<typename V>
//      	      inline bool operator()(::std::pair<Kmer, V> const & x, Kmer const & y) const {
//      	        return this->operator()(x.first, y);
//      	      }
//      	      template<typename V>
//      	      inline bool operator()(::std::pair<const Kmer, V> const & x, Kmer const & y) const {
//      	        return this->operator()(x.first, y);
//      	      }
//      	      template<typename V>
//      	      inline bool operator()(Kmer const & x, ::std::pair<Kmer, V> const & y) const {
//      	        return this->operator()(x, y.first);
//      	      }
//      	      template<typename V>
//      	      inline bool operator()(Kmer const & x, ::std::pair<const Kmer, V> const & y) const {
//      	        return this->operator()(x, y.first);
//      	      }
//      	      template<typename V>
//      	      inline bool operator()(::std::pair<Kmer, V> const & x, ::std::pair<Kmer, V> const & y) const {
//      	        return this->operator()(x.first, y.first);
//      	      }
//      	      template<typename V>
//      	      inline bool operator()(::std::pair<const Kmer, V> const & x, ::std::pair<const Kmer, V> const & y) const {
//      	        return this->operator()(x.first, y.first);
//      	      }
//      	  };




      }  // namespace sparsehash

      // usage  apply transform first, then apply hash.
    } // namespace hash
  } // namespace kmer
} // namespace bliss



#endif /* KMER_HASH_HPP_ */
