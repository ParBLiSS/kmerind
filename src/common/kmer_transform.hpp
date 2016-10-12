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
#ifndef KMER_TRANSFORM_HPP_
#define KMER_TRANSFORM_HPP_

#include <tuple>  // for hash - std::pair
#include <exception>  // for hash - std::system_error
#include <algorithm>
#include <type_traits>  // enable_if

#include "common/kmer.hpp"

namespace bliss {

  namespace kmer
  {

    namespace transform {

      // QUESTION:  xor of hash, or hash of xor?.  second is faster.  Also if there is GC-AT imbalance, xor of raw sequence kind of flattens the distribution, so hash input is now more even.
//
//      template <typename KMER>
//      struct identity {
//          inline KMER operator()(KMER const & x) const {
//            return std::forward<const KMER>(x);
//          }
//          template <typename VAL>
//          inline ::std::pair<KMER, VAL> operator()(std::pair<KMER, VAL> const & x) const {
//            return std::forward<const ::std::pair<KMER, VAL> >(x);
//          }
//          template <typename VAL>
//          inline ::std::pair<const KMER, VAL> operator()(std::pair<const KMER, VAL> const & x) const {
//            return std::forward<const ::std::pair<const KMER, VAL> >(x);
//          }
//      };

      template <typename KMER>
      struct xor_rev_comp {
          inline KMER operator()(KMER const & x) const {
            return x ^ x.reverse_complement();
          }
          inline KMER operator()(KMER const & x, KMER const & rc) const  {
            return x ^ rc;
          }
          template <typename VAL>
          inline ::std::pair<KMER, VAL> operator()(std::pair<KMER, VAL> const & x) const {
              return std::pair<KMER, VAL>(operator()(x.first), x.second);
          }
          template <typename VAL>
          inline ::std::pair<const KMER, VAL> operator()(std::pair<const KMER, VAL> const & x) const {
              return std::pair<const KMER, VAL>(operator()(x.first), x.second);
          }
      };

      template <typename KMER>
      struct lex_less {
          inline KMER operator()(KMER const & x) const  {
            auto y = x.reverse_complement();
            return (x < y) ? x : y;
          }
          inline KMER operator()(KMER const & x, KMER const & rc) const  {
            return (x < rc) ? x : rc;
          }
          template <typename VAL>
          inline ::std::pair<KMER, VAL> operator()(std::pair<KMER, VAL> const & x) const {
              return std::pair<KMER, VAL>(operator()(x.first), x.second);
          }
          template <typename VAL>
          inline ::std::pair<const KMER, VAL> operator()(std::pair<const KMER, VAL> const & x) const {
              return std::pair<const KMER, VAL>(operator()(x.first), x.second);
          }
      };

      template <typename KMER>
      struct lex_greater {
          inline KMER operator()(KMER const & x) const  {
            auto y = x.reverse_complement();
            return (x > y) ? x : y;
          }
          inline KMER operator()(KMER const & x, KMER const & rc) const  {
            return (x > rc) ? x : rc;
          }
          template <typename VAL>
          inline ::std::pair<KMER, VAL> operator()(std::pair<KMER, VAL> const & x) const {
            return std::pair<KMER, VAL>(operator()(x.first), x.second);
          }
          template <typename VAL>
          inline ::std::pair<const KMER, VAL> operator()(std::pair<const KMER, VAL> const & x) const {
            return std::pair<const KMER, VAL>(operator()(x.first), x.second);
          }

      };


//      template <typename KMER, template <typename> class TRANS>
//      struct tuple_transform {
//          TRANS<KMER> transform;
//
//          inline KMER operator()(KMER & x) {
//            x = transform(x);
//            return x;
//          }
//
//          template <typename VAL>
//          inline ::std::pair<KMER, VAL> operator()(std::pair<KMER, VAL> & x) {
//            x.first = transform(x.first);
//            return x;
//          }
//
//      };

    } // namespace transform


  } // namespace kmer
} // namespace bliss



#endif /* KMER_TRANSFORM_HPP_ */
