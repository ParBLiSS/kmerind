/**
 * @file		kmer_functors.hpp
 * @ingroup
 * @author	tpan
 * @brief
 * @details
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef KMER_FUNCTORS_HPP_
#define KMER_FUNCTORS_HPP_

#include <cmath>

#include "common/AlphabetTraits.hpp"
#include "io/fastq_iterator.hpp"
#include "index/kmer_index_element.hpp"

namespace bliss
{
  namespace index
  {
    /**
     * template Param:  T is type of data
     *                  P is number of bins to hash to.
     */
    template<typename T, int P>
    struct XorModulus
    {


    };


    template<typename Sequence, typename KmerIndex, typename HashFunction = bliss::index::XorModulus>
    struct generate_kmer
    {
        typedef Sequence SequenceType;
        typedef typename Sequence::IteratorType BaseIterType;
        typedef typename Sequence::AlphabetType AlphabetType;
        typedef KmerIndex KmerIndexType;
        typedef typename KmerIndex::KmerType KmerValueType;

        const int K = KmerIndex::KmerSize::size;

        KmerIndex kmer;

        KmerValueType revcomp;
        uint16_t pos;
        HashFunction h;

        static constexpr BitSizeType nBits =
            bliss::AlphabetTraits<AlphabetType>::getBitsPerChar();
        static constexpr BitSizeType shift =
            bliss::AlphabetTraits<AlphabetType>::getBitsPerChar() * (K - 1);
        static constexpr AlphabetSizeType max =
            bliss::AlphabetTraits<AlphabetType>::getSize() - 1;

        static constexpr int word_size = sizeof(KmerValueType) * 8;
        static constexpr KmerValueType mask_reverse = ~(static_cast<KmerValueType>(0))
            >> (word_size - shift - nBits);
        //    static constexpr KmerValueType mask_lower_half = ~(static_cast<KmerValueType>(0))
        //        >> (word_size - nBits * (K + 1) / 2);

        generate_kmer(const bliss::io::fastq_sequence_id &_rid)
            : kmer(), revcomp(0), pos(0)
        {
          kmer.id = _rid;
        }

        size_t operator()(BaseIterType &iter)
        {
          // store the kmer information.
          char val = AlphabetType::FROM_ASCII[static_cast<size_t>(*iter)];
          kmer.kmer >>= nBits;
          kmer.kmer |= (static_cast<KmerValueType>(val) << shift);
          kmer.id.components.pos = pos;

          // generate the rev complement
          char complement = max - val;
          revcomp <<= nBits;
          revcomp |= static_cast<KmerValueType>(complement);
          revcomp &= mask_reverse;

          ++iter;
          ++pos;
          return 1;
        }

        /**
         *
         * @return  pair, with Index Value type and Index structure.  indexValue may be an xor or some other computed value.
         */
        std::pair<KmerValueType, KmerIndex> operator()()
        {
          // compute the reverse complement
          KmerValueType xored = kmer.kmer ^ revcomp;

          //      KmerValueType xored_recoverable = (xored & ~mask_lower_half) | (_kmer.forward & mask_lower_half);

          return std::pair<KmerValueType, KmerIndex>(xored, kmer);
        }

    };


    /**
     * Phred Scores.  see http://en.wikipedia.org/wiki/FASTQ_format.
     *
     * creating a look up table instead of computing (by calling std::log, exp, etc.).
     * saves about 3 seconds on a 35MB read file.
     */
    template<typename T>
    struct SangerToLogProbCorrect
    {
        static_assert(std::is_floating_point<T>::value, "generate_qual output needs to be floating point type");
        static constexpr size_t size = 94;

        typedef T ValueType;

        static constexpr T offset = 33;

        static constexpr size_t min = 0;
        static constexpr size_t max = 93;

        static constexpr T log2_10DivNeg10 = std::log2(10.0) / -10.0;
        constexpr ValueType operator()(const size_t v)
        {
          // some limits: v / -10 has to be negative as this becomes probability, so v > 0
          //
          return
              v < min ? std::numeric_limits<ValueType>::lowest() :
              v > max ? std::numeric_limits<ValueType>::lowest() :
              v == min ?
                  std::numeric_limits<ValueType>::lowest() :
                  std::log2(1.0 - std::exp2(static_cast<ValueType>(v) * log2_10DivNeg10));
        }
    };

    /**
     * compute kmer quality based on phred quality score.
     *
     * Q = -10 log_10 P
     *
     * where P is probability of error.
     *
     * probability of kmer being correct is: (1-p_1)(1-p_2)...(1-p_k)
     * probability of kmer being incorrect then is 1 - (p(correct...))
     *
     * Phrad adds the Phred score, which amounts to p_1 * p_2 * ... * p_k because of the logarithm.
     *  (probability that all bases in the sequence are incorrect.)
     * Quake is computes a kmer score based on probability in the phred score.  not clear exactly how.
     *
     * value range is in ascii from !(33) to ~(126), using sanger encoding.
     *
     * // now computing k-correct with approximation gives
     * //    1 - sum_1..k(p_i) + 1/2 sum_1..k sum_1..k (p_i)(p_j) - 3rd order term + 4th order term ...
     * //
     * //   p(incorrect) is sum_1..k(p_i) - 1/2 sum_1..k sum_1..k (p_i)(p_j) + 3rd order term - 4th order term ...
     * //
     * // =====>>>  max_1..k(p_i) < p(incorrect) < sum_1..k(p_i).
     * //
     * // approximate with linear terms.  (i.e. probability of 2 or more bases being called incorrectly is low).  bound is not that tight.
     * //
     * // =====>>>  max_1..k(p_i) < p(
     * //
     * // max_1..k(p_i) ~= min_1..k(q_i)
     *
     *
     *  each term in (1-p_i) is (1- 10^(q/-10)) = (1 - e ^ ((log 10)(q/-10)))
     *  accumulate using sum(log(term))., expand using e^(sum(log(term))), calculate -10 log_10(1-cumulative).
     *
     *  This was a bottleneck.  switched to using a look up table instead of computing log all the time.
     *
     *
     *  This functor is only available if KmerIndex has Quality Score.
     */
    template<typename Sequence, typename KmerIndex, typename Encoding,
      typename HasQual = typename std::enable_if<std::is_same<KmerIndex, KmerIndexElementWithIdAndQuality>::value>::type >
    struct generate_qual
    {
        typedef typename Sequence::IteratorType BaseIterType;

        typedef KmerIndex KmerType;
        typedef typename KmerIndex::QualityType QualityType;

        constexpr int K = KmerIndex::KmerSize::size;

        int kmer_pos;
        QualityType internal;
        QualityType terms[K];


        int pos;
        int zeroCount;

        // can't use auto keyword.  declare and initialize in class declaration
        // then "define" but not initialize outside class declaration, again.
        static constexpr std::array<typename Encoding::value_type, Encoding::size> lut =
            make_array<Encoding::size>(Encoding());

        generate_qual()
            : kmer_pos(0), pos(0)
        {
          for (int i = 0; i < K; ++i)
          {
            terms[i] = std::numeric_limits<QualityType>::lowest();
          }
          zeroCount = K;

    //      for (int i = 0; i < Encoding::size; ++i)
    //      {
    //        printf("lut: %d=%lf\n", i, lut[i]);
    //      }

        }

        size_t operator()(BaseIterType &iter)
        {
          int oldpos = kmer_pos;

          // drop the old value
          QualityType oldval = terms[pos];

          // add the new value,       // update the position  - circular queue
          QualityType newval = lut[*iter - Encoding::offset]; // this is for Sanger encoding.
          terms[pos] = newval;
          pos = (pos + 1) % K;

          // save the old zero count.
          int oldZeroCount = zeroCount;

          // update the zero count
          if (newval == std::numeric_limits<QualityType>::lowest())
          {
            //        printf("ZERO!\n");
            ++zeroCount;
          }
          if (oldval == std::numeric_limits<QualityType>::lowest())
          {
            --zeroCount;
          }

          // if any is zero, then return 0.
          if (zeroCount > 0)
          {
            internal = 0.0;
          }
          else
          {
            // else there is no zero, so valid values.

            if (oldZeroCount == 1)
            {
              //printf("HAD ZEROS!\n");
              // removed a zero.  so recalculate.
              internal = 0.0;
              for (int i = 0; i < K; ++i)
              {
                internal += terms[i];
              }
            }
            else
            {
              // there was not a zero.  so update.
              internal = internal + newval - oldval;
            }

          }

          //      printf("%d qual %d oldval = %f, newval = %f, internal = %f, output val = %f\n", kmer_pos, *iter, oldval, newval, internal, value);
          //      std::fflush(stdout);

          // move the iterator.
          ++iter;
          ++kmer_pos;

          return kmer_pos - oldpos;
        }

        QualityType operator()()
        {
          // compute prob of kmer being incorrect.
          if (fabs(internal) < std::numeric_limits<QualityType>::epsilon())
          {
            //printf("confident kmer!\n");
            return Encoding::max;
          }
          else
          {
            return -10.0 * std::log10(1.0 - std::exp2(internal));
          }
        }

    };

  } /* namespace io */
} /* namespace bliss */



#endif /* KMER_FUNCTORS_HPP_ */
