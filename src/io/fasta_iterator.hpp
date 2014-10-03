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
     *            this keeps a 16 bit sequence id indicating its index
     *                       an 8 bit file id
     *                       a 40 bit as global offset in the file, broken up into 32 bit id and 8 significant bit id       
     *
     *            A separate FASTA version will have a different fields but keeps the same 64 bit total length.
     *
     *            file size is at the moment limited to 1TB (40 bits) in number of bytes.
     */
    union FASTASequenceId
    {
        /// the concatenation of the id components as a single unsigned 64 bit field
        uint64_t composite;

        /// the id field components
        struct
        {  
            /// sequence's id, the index of the sequence in the containing file
            uint16_t seq_id;
            /// id of fasta file
            uint8_t file_id;
            /// offset within the read , upper 8 of 40 bytes  (potentially as offset in the containing file)
            uint32_t pos_msb;
            /// offset within the read , lower 32 of 40 bytes (potentially as offset in the containing file)
            uint8_t pos_id;
        } components;
    };



  }
}
