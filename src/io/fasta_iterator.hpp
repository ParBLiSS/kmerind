/**
 * @file    fasta_iterator.hpp
 * @ingroup bliss::io
 * @author  Chirag Jain <cjain@gatech.edu>
 * @date    Oct 3, 2014
 * @brief   Contains an iterator for traversing a sequence data, record by record, a FASTA (or multiFASTA) format specific parser for use by the iterator,
 *          Sequence Id datatype definition, and 1 sequence type (without quality score iterators)
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add Licence
 *
 */


#ifndef FastaIterator_HPP_
#define FastaIterator_HPP_

// C includes
#include <cassert>
#include <cstdint>

// C++ STL includes
#include <iterator>
#include <algorithm>
#include <sstream>
#include <type_traits>

// own includes
#include <io/io_exception.hpp>
#include <utils/logging.h>
#include <iterators/function_traits.hpp>
#include <partition/range.hpp>

namespace bliss
{
  namespace io
  {
    /**
     * @class     bliss::io::FASTASequenceId
     * @brief     represents a fasta sequence's id, also used for id of the FASTA file, and for position inside a FASTA sequence.
     * @detail    this is set up as a union to allow easy serialization
     *            and parsing of the content.
     *
     *            this keeps a 48 bit sequence ID, broken up into a 32 bit id and 16 bit significant bit id
     *                       a 16 bit file id
     *                       a 64 bit position within the sequence, broken up into 32 bit id and 32 significant bit id.
     *
     *            FASTQ version has different fields but keeps the same 64 bit total length.
     *            This design of representing id can use the FASTQ parser's logic 
     *
     */
    union FASTASequenceId
    {
      /// the concatenation of the id components as a single unsigned 128 bit field
      unsigned __int128 composite;

      /// the id field components
      struct
      {  
        /// sequence's id, lower 32 of 48 bits (potentially as offset in the containing file)
        uint32_t seq_id;
        /// sequence's id, upper 16 of 48 bits (potentially as offset in the containing file)
        uint16_t seq_msb;
        /// id of fasta file
        uint16_t file_id;
        /// offset within the FaSTA record, lower 32 of 64 bytes
        uint32_t pos_id;
        /// offset within the FASTA record, upper 32 of 64 bytes 
        uint32_t pos_msb;
      } components;
    };

    template<typename Iterator, typename Alphabet>
      struct Sequence
      {
        /// Type of the sequence elements
        using ValueType = typename std::iterator_traits<Iterator>::value_type;
        /// The alphabet this sequence uses
        typedef Alphabet AlphabetType;
        /// Iterator type for traversing the sequence.
        typedef Iterator IteratorType;

        /// begin iterator for the sequence
        Iterator seqBegin;
        /// end iterator for the sequence.
        Iterator seqEnd;
      };

  }
}
