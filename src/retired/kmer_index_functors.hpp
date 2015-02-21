/**
 * @file    kmer_index_functors.hpp
 * @ingroup retired
 * @author  Tony Pan <tpan7@gatech.edu>
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
#include <common/alphabets.hpp>
#include <common/alphabet_traits.hpp>      // for mapping char to alphabet value
#include <io/sequence_iterator.hpp>          // for Id.
#include <retired/kmer_index_element.hpp>   // for the basic kmer index element structures.
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
            bliss::common::AlphabetTraits<AlphabetType>::getBitsPerChar();
        static constexpr BitSizeType shift =
            bliss::common::AlphabetTraits<AlphabetType>::getBitsPerChar() * (KmerValueType::size - 1);
        static constexpr AlphabetSizeType max =
            bliss::common::AlphabetTraits<AlphabetType>::getSize() - 1;

        static constexpr int word_size = sizeof(KmerValueType) * 8;
//        static constexpr KmerValueType mask_reverse = ~(KmerValueType());
//            >> (word_size - shift - nBits);

        //.constructor taking the fastq sequence id for a read
      public:

        /// generates one k-mer per call from the underlying read
        size_t operator()(BaseIterType &iter)
        {
          // store the kmer information.
          uint8_t val = AlphabetType::FROM_ASCII[static_cast<size_t>(*iter)];
          kmer.nextFromChar(val);
//          kmer >>= nBits;
//          kmer |= (static_cast<KmerValueType>(val) << shift);

          // generate the rev complement

          revcomp.nextReverseFromChar(AlphabetType::TO_COMPLEMENT[val]);
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
        typedef typename Sequence::IdType           SeqIdType;

      protected:
        typedef typename Sequence::AlphabetType     AlphabetType;
        /// Current Kmer buffer (k-mer + index + quality)
        KmerIndexType kmer;
        /// Current reverse complement (just k-mer)
        KmerValueType revcomp;

      protected:


        static constexpr BitSizeType nBits =
            bliss::common::AlphabetTraits<AlphabetType>::getBitsPerChar();
        static constexpr BitSizeType shift =
            bliss::common::AlphabetTraits<AlphabetType>::getBitsPerChar() * (KmerIndexType::SizeType::size - 1);
        static constexpr AlphabetSizeType max =
            bliss::common::AlphabetTraits<AlphabetType>::getSize() - 1;

        static constexpr int word_size = sizeof(KmerValueType) * 8;
        static constexpr KmerValueType mask_reverse = ~(static_cast<KmerValueType>(0))
            >> (word_size - shift - nBits);
        //    static constexpr KmerValueType mask_lower_half = ~(static_cast<KmerValueType>(0))
        //        >> (word_size - nBits * (K + 1) / 2);

        //.constructor taking the fastq sequence id for a read
      public:
        generate_kmer(const SeqIdType &_rid)
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
        typedef typename Sequence::IdType           SeqIdType;

        typedef std::pair<KmerValueType, KmerIndexType> OutputType;       // key-value pair.


        /// kmer and rev complement are in the base class.
        /// current position in the read sequence
        uint16_t pos;

      protected:
        typedef typename Sequence::AlphabetType     AlphabetType;
        typedef generate_kmer<Sequence, KmerIndex, std::false_type> BaseClassType;

      public:

        generate_kmer(const SeqIdType &_rid)
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






  } /* namespace io */
} /* namespace bliss */



#endif /* KMER_INDEX_FUNCTORS_HPP_ */
