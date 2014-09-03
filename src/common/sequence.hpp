/**
 * @file    sequence.hpp
 * @ingroup
 * @author  tpan
 * @brief
 * @details
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
     * @class     bliss::io::FASTQSequenceId
     * @brief     represents a fastq sequence's id, also used for id of the FASTQ file, and for position inside a FASTQ sequence..
     * @details    this is set up as a union to allow easy serialization
     *            and parsing of the content.
     *            this keeps a 40 bit sequence ID, broken up into a 32 bit id and 8 bit significant bit id
     *                       an 8 bit file id
     *                       a 16 bit position within the sequence.
     *
     *            A separate FASTA version will have a different fields but keeps the same 64 bit total length.
     *
     *            file size is at the moment limited to 1TB (40 bits) in number of bytes.
     */
    union FASTQSequenceId
    {
        /// the concatenation of the id components as a single unsigned 64 bit field.  should use only 40 bits
        uint64_t file_pos;

        /// the id field components.  anonymous struct
        struct
        {
            /// sequence's id, lower 32 of 40 bits (potentially as offset in the containing file)
            uint32_t seq_id;
            /// sequence's id, upper 8 of 40 bits (potentially as offset in the containing file)
            uint8_t seq_id_msb;
            /// id of fastq file
            uint8_t file_id;
            /// offset within the read.  Default 0 refers to the whole sequence
            uint16_t pos;
        };
    };

    /**
     * @class     bliss::io::FASTASequenceId
     * @brief     represents a fasta sequence's id, also used for id of the FASTA file, and for position inside a FASTA sequence..
     * @details   Since we are using the in-file position of the sequence as the sequence id, we need
     *
     *            8 bits for file id
     *            40 bits for sequence id
     *            40 bits for position within the sequence.
     *
     *            this is set up as a union of to allow easy serialization
     *            and parsing of the content.
     *
     *            file size is at the moment limited to 1TB (40 bits) in number of bytes.
     */
    struct FASTASequenceId
    {
        /// the concatenation of the id components as a single unsigned 64 bit field for sequence position inside file.
        union {
            uint64_t file_pos;

            /// the id field components.  anonymous struct
            struct
            {
                /// sequence's id, lower 32 of 40 bits (potentially as offset in the containing file)
                uint32_t seq_id;
                /// sequence's id, upper 8 of 40 bits (potentially as offset in the containing file)
                uint8_t seq_id_msb;
                /// id of fastq file
                uint8_t file_id;
            };
        };
        /// offset within the read.  Default 0 refers to the whole sequence
        uint64_t pos;
    };
    /**
     * @class     bliss::io::Sequence
     * @brief     represents a biological sequence, and provides iterators for traversing the sequence.
     * @details   internally we store 2 iterators, start and end, instead of start + length.
     *            reason is that using length is best for random access iterators and we don't have guarantee of that.
     * @tparam Iterator   allows walking through the sequence data.
     * @tparam Id     the format and components of the id of a sequence. Differs from FASTA to FASTQ and based on length of reads.  defaults to FASTQSequenceId
     */
    template<typename Iterator, typename Id=FASTQSequenceId>
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
     * @tparam Quality    data type for quality scores for each base
     * @tparam IdType     the format and components of the id of a sequence. Differs from FASTA to FASTQ and based on length of reads.  defaults to FASTQSequenceId
     */
    template<typename Iterator, typename Quality, typename Id=FASTQSequenceId>
    struct SequenceWithQuality : public Sequence<Iterator, Id>
    {
        /// Iterator type for traversing the sequence.
        typedef Iterator IteratorType;
        /// type for the id struct/union
        typedef Id IdType;
        /// the desired quality score's type, e.g. double
        typedef Quality QualityType;

        /// begin iterator for the quality scores
        Iterator qualBegin;
        /// end iterator for the quality scores
        Iterator qualEnd;
    };


  } /* namespace common */
} /* namespace bliss */

#endif /* SEQUENCE_HPP_ */
