/**
 * @file    kmer_index_functors.hpp
 * @ingroup
 * @author  tpan
 * @brief   contains functors for generating kmers and its quality scores.
 * @details uses a containment strategy for generating kmers w/ or w/o id and quality score.
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef KMER_INDEX_FUNCTORS_HPP_
#define KMER_INDEX_FUNCTORS_HPP_

#include <cmath>  // for quality score computation

#include <array>  // for quality score lookup table
#include <type_traits>
#include <utility>

#include <utils/constexpr_array.hpp>      // for precomputing the quality scores
#include <common/AlphabetTraits.hpp>      // for mapping char to alphabet value
#include <io/fastq_iterator.hpp>          // for Id.
#include <index/kmer_index_element.hpp>   // for the basic kmer index element structures.
#include <common/sequence.hpp>

namespace bliss
{
  namespace index
  {


    template<typename Sequence, typename Alphabet, typename KmerType>
    class generate_kmer_simple {
      public:
        typedef Sequence                            SequenceType;
        typedef typename Sequence::IteratorType     BaseIterType;
        typedef KmerType    KmerValueType;
        typedef std::pair<KmerValueType, KmerValueType> OutputType;       // key-value pair.

      protected:
        typedef Alphabet               AlphabetType;
        /// Current Kmer buffer (k-mer + index + quality)
        KmerValueType kmer;
        /// Current reverse complement (just k-mer)
        KmerValueType revcomp;

      protected:


        static constexpr BitSizeType nBits =
            bliss::AlphabetTraits<AlphabetType>::getBitsPerChar();
        static constexpr BitSizeType shift =
            bliss::AlphabetTraits<AlphabetType>::getBitsPerChar() * (KmerValueType::size - 1);
        static constexpr AlphabetSizeType max =
            bliss::AlphabetTraits<AlphabetType>::getSize() - 1;

        static constexpr int word_size = sizeof(KmerValueType) * 8;
//        static constexpr KmerValueType mask_reverse = ~(KmerValueType());
//            >> (word_size - shift - nBits);

        //.constructor taking the fastq sequence id for a read
      public:

        /// generates one k-mer per call from the underlying read
        size_t operator()(BaseIterType &iter)
        {
          // store the kmer information.
          char val = AlphabetType::FROM_ASCII[static_cast<size_t>(*iter)];
          kmer.nextFromChar(val);
//          kmer >>= nBits;
//          kmer |= (static_cast<KmerValueType>(val) << shift);

          // generate the rev complement
          char complement = AlphabetType::TO_COMPLEMENT[val];
          kmer.nextFromChar(complement);
//          revcomp <<= nBits;
//          revcomp |= static_cast<KmerValueType>(complement);
//          //revcomp &= mask_reverse;

          ++iter;

          return 1;
        }

        /**
         *
         */
        OutputType operator()()
        {
          return OutputType(kmer, revcomp);   // constructing a new result.
        }

    };


    // TODO: this is all deprecated, but some parts of it might be useful

    // KmerIndex = KmerIndexElement or KmerIndexElementWithQuality, ...
    template<typename Sequence, typename KmerIndex, typename has_id = std::false_type>
    class generate_kmer {
      public:
        typedef Sequence                            SequenceType;
        typedef typename Sequence::IteratorType     BaseIterType;
        typedef KmerIndex                           KmerIndexType;
        typedef typename KmerIndexType::KmerType    KmerValueType;
        typedef std::pair<KmerValueType, KmerIndexType> OutputType;       // key-value pair.

      protected:
        typedef typename Sequence::AlphabetType     AlphabetType;
        /// Current Kmer buffer (k-mer + index + quality)
        KmerIndexType kmer;
        /// Current reverse complement (just k-mer)
        KmerValueType revcomp;

      protected:


        static constexpr BitSizeType nBits =
            bliss::AlphabetTraits<AlphabetType>::getBitsPerChar();
        static constexpr BitSizeType shift =
            bliss::AlphabetTraits<AlphabetType>::getBitsPerChar() * (KmerIndexType::SizeType::size - 1);
        static constexpr AlphabetSizeType max =
            bliss::AlphabetTraits<AlphabetType>::getSize() - 1;

        static constexpr int word_size = sizeof(KmerValueType) * 8;
        static constexpr KmerValueType mask_reverse = ~(static_cast<KmerValueType>(0))
            >> (word_size - shift - nBits);
        //    static constexpr KmerValueType mask_lower_half = ~(static_cast<KmerValueType>(0))
        //        >> (word_size - nBits * (K + 1) / 2);

        //.constructor taking the fastq sequence id for a read
      public:
        generate_kmer(const bliss::common::FASTQSequenceId &_rid)
          : kmer(), revcomp(0)
        {}  // not setting kmer.id here because this class does not have kmer id.

        /// generates one k-mer per call from the underlying read
        size_t operator()(BaseIterType &iter)
        {
          // store the kmer information.
          char val = AlphabetType::FROM_ASCII[static_cast<size_t>(*iter)];
          kmer.kmer >>= nBits;
          kmer.kmer |= (static_cast<KmerValueType>(val) << shift);

          // generate the rev complement
          char complement = max - val;
          revcomp <<= nBits;
          revcomp |= static_cast<KmerValueType>(complement);
          revcomp &= mask_reverse;

          ++iter;

          return 1;
        }

        /**
         *
         * @return  pair, with Index Value type and Index structure.  indexValue is the revcomp.  hash function can then use this.
         */
        OutputType operator()()
        {
          //      KmerValueType xored_recoverable = (xored & ~mask_lower_half) | (_kmer.forward & mask_lower_half);

          // XOR k-mer key with it's reverse complement and combine with
          // full index structure
          return OutputType(revcomp ^ kmer.kmer, kmer);   // constructing a new result.
        }

    };

    /**
     * kmer generator operator type
     * The third template parameter acts as a member checker to make sure there is an id field and it's of type bliss::common::FASTQSequenceId
     */
    template<typename Sequence, typename KmerIndex>
    class generate_kmer<Sequence, KmerIndex,
      typename std::enable_if<!std::is_void<decltype(std::declval<KmerIndex>().id)>::value, std::true_type >::type >
      : generate_kmer<Sequence, KmerIndex, std::false_type>
    {
      public:
        typedef Sequence                            SequenceType;
        typedef typename Sequence::IteratorType     BaseIterType;
        typedef KmerIndex                           KmerIndexType;
        typedef typename KmerIndexType::KmerType    KmerValueType;

        typedef std::pair<KmerValueType, KmerIndexType> OutputType;       // key-value pair.


        /// kmer and rev complement are in the base class.
        /// current position in the read sequence
        uint16_t pos;

      protected:
        typedef typename Sequence::AlphabetType     AlphabetType;
        typedef generate_kmer<Sequence, KmerIndex, std::false_type> BaseClassType;

      public:

        generate_kmer(const bliss::common::FASTQSequenceId &_rid)
            : BaseClassType(_rid), pos(0)
        {
          this->kmer.id = _rid;
        }

        size_t operator()(BaseIterType &iter)
        {
          size_t count = BaseClassType::operator()(iter);

          this->kmer.id.pos = pos;
          ++pos;

          return count;
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
        typedef std::array<T, size> LUTType;

        static constexpr T offset = 33;

        static constexpr size_t min = 0;
        static constexpr size_t max = size - 1;

        // can't use auto keyword.  declare and initialize in class declaration
        // then "define" but not initialize outside class declaration, again.
        static constexpr LUTType lut = bliss::utils::make_array<size>(SangerToLogProbCorrect<T>());


        static constexpr T log2_10DivNeg10 = std::log2(10.0) / -10.0;
        constexpr ValueType operator()(const size_t &v)
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


    // "define" but not initialize outside class declaration, again.
    // NOTE:  Do it here is fine.  do it in class file or even where it's used is NOT.  probably due to templating.
    template<typename T> constexpr typename SangerToLogProbCorrect<T>::LUTType SangerToLogProbCorrect<T>::lut;

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
    template<typename Sequence, unsigned int KmerSize, typename Encoding>
    struct generate_qual
    {
        typedef typename Sequence::IteratorType BaseIterType;
        typedef Sequence                        SequenceType;
        typedef typename Encoding::ValueType    QualityType;

        int kmer_pos;
        QualityType internal;
        QualityType terms[KmerSize];


        int pos;
        unsigned int zeroCount;

        generate_qual()
            : kmer_pos(0), pos(0)
        {
          for (unsigned int i = 0; i < KmerSize; ++i)
          {
            terms[i] = std::numeric_limits<QualityType>::lowest();
          }
          zeroCount = KmerSize;

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
          QualityType newval = Encoding::lut[*iter - Encoding::offset]; // this is for Sanger encoding.
          terms[pos] = newval;
          pos = (pos + 1) % KmerSize;

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
              for (unsigned int i = 0; i < KmerSize; ++i)
              {
                internal += terms[i];
              }
            }
            else
            {
              // there was not a zero.  so update.
              internal += newval - oldval;
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



#endif /* KMER_INDEX_FUNCTORS_HPP_ */
