/**
 * @file    sequence.hpp
 * @ingroup common
 * @author  Tony Pan <tpan7@gatech.edu>
 * @brief   represents a sequence of characters in computer memory.
 * @details Internally, contains a start and end iterator for traversing the characters.
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef SEQUENCE_HPP_
#define SEQUENCE_HPP_

namespace bliss
{
  namespace common
  {

    /**
     * @class     bliss::io::Sequence
     * @brief     represents a biological sequence, and provides iterators for traversing the sequence.
     * @details   internally we store 2 iterators, start and end, instead of start + length.
     *            reason is that using length is best for random access iterators and we don't have guarantee of that.
     * @tparam Iterator   allows walking through the sequence data.
     * @tparam Id     the format and components of the id of a sequence. Differs from FASTA to FASTQ and based on length of reads.  defaults to FASTQ::SequenceId
     */
    template<typename Iterator, typename Id>
    struct Sequence
    {
        /// Iterator type for traversing the sequence.
        typedef Iterator IteratorType;
        /// type for the id struct/union
        typedef Id IdType;
        /// the desired quality score's type. hardcoded to void.
        typedef void QualityType;

        /// begin iterator for the sequence
        Iterator seqBegin;
        /// end iterator for the sequence.
        Iterator seqEnd;
        /// id of the sequence.
        Id id;
    };

    /**
     * @class     bliss::io::SequenceWithQuality
     * @brief     extension of bliss::io::Sequence to include quality scores for each position.
     *
     * @tparam Iterator   allows walking through the sequence data.
     * @tparam IdType     the format and components of the id of a sequence. Differs from FASTA to FASTQ and based on length of reads.  defaults to FASTQ::SequenceId
     */
    template<typename Iterator, typename Id>
    struct SequenceWithQuality : public Sequence<Iterator, Id>
    {
        /// Iterator type for traversing the sequence.
        typedef Iterator IteratorType;
        /// type for the id struct/union
        typedef Id IdType;

        /// begin iterator for the quality scores
        Iterator qualBegin;
        /// end iterator for the quality scores
        Iterator qualEnd;
    };


  } /* namespace common */
} /* namespace bliss */

#endif /* SEQUENCE_HPP_ */
