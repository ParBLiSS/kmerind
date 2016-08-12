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
 * @file    sequence.hpp
 * @ingroup common
 * @author  Tony Pan <tpan7@gatech.edu>
 * @brief   represents a sequence of characters in computer memory, and the associated sequence id (mostly as position in file)
 * @details Internally, contains a start and end iterator for traversing the characters, plus a sequence id.
 *
 */
#ifndef SEQUENCE_HPP_
#define SEQUENCE_HPP_

namespace bliss
{
  namespace common
  {
    // design:  packed values can easily be copied without multiple accessors and is compact for comm.
    //          trickles to sequenceIdIterator specializations for increment. (to add sequence types, subclass sequenceIdIterator and add additional member function specializations)
    //          convertible from SequenceId to Long or Short Sequence Ids
    //          SequenceIdIterator may need specialization  for initialization.
    //          flow:  1. file->seqParser->SequenceId (in file coordinates) in Sequence
    //                 2. Sequence -> kmerParser -> Long/Short SequenceId for Kmer position
    //               During transform 1, there is info for start/end of record, file position of record, file id, start/end of sequence, and sequence uid (via global prefix).
    //               During transform 2, it's primarily incrementing the position within sequence.
    //          FASTA reading may require that the user specify whether to use ShortSequenceKmerId or LongSequenceKmerId.
    //          ShortSequenceKmerId will encounter more constructions, so make that fast - in SequenceId, encode in ShortSequenceKmerId form.


    /**
     * @brief SequenceId is used to identify a sequence in a file.  This is the basic type that is suitable for position in files.
     * @details  this class holds the data from which a KmerId can be created.
     * @note  This type is NOT sent frequently over mpi.
     *        virtual functions break this.  so no virtual destructor or other virtual functions.
     *        In addition, compilers can align the members however.  so preferably the data content should be
     *        contiguous, otherwise a mxx::datatype specialization should be used.
     *        
     *        this keeps a 64 bit position within the file.   (16 Exabyte per file)
     *                   a 64 bit seq_id within the file      (16 x 10^18 sequences in 1 file.  )
     *                   a 16 bit file id                     (64K files)
     *
     *        Future and current KmerId classes can extract the meaningful part based on realistic limits for each user applications.
     *          e.g. 1TB files, 32K sequences per file, etc.
     *
     * @note  not using inheritance here because we don't use SequenceId polymorphically.
     */
    struct SequenceId {

        /// Id of the sequence.  Typically, seq_id is same as position in the file.
        size_t pos_in_file;
        size_t seq_id;
        uint16_t file_id;

        /**
         * @brief << operator to write out SequenceId
         * @param[in/out] ost   output stream to which the content is directed.
         * @param[in]     seq_id    sequence id object to write out
         * @return              output stream object
         */
        friend std::ostream& operator<<(std::ostream& ost, const SequenceId & seq_id)
        {
          // friend keyword signals that this overrides an externally declared function
          ost << " SeqId: file=" << seq_id.file_id << " pos=" << seq_id.pos_in_file << " seq_id=" << seq_id.seq_id;
          return ost;
        }

        SequenceId() = default;

        SequenceId(size_t const & _file_pos, size_t const & _seq_id = 0, uint16_t const & _file_id = 0) :
          pos_in_file(_file_pos), seq_id(_seq_id), file_id(_file_id) {};

        SequenceId(SequenceId const & other) : pos_in_file(other.pos_in_file), seq_id(other.seq_id), file_id(other.file_id) {}

        SequenceId& operator=(SequenceId const & other) {
          file_id = other.file_id;
          pos_in_file = other.pos_in_file;
          seq_id = other.seq_id;

          return *this;
        }

        /// getter for sequence id  (here, it's a 40 bit unique id == start of seq)
        size_t get_id() const { return seq_id; }

        /// get position in file
        size_t get_pos() const { return pos_in_file ; }

        /// getter for file id
        uint8_t get_file_id() const { return file_id; }

    };


    /**
     * @class     bliss::common::ShortSequenceKmerId
     * @brief     an Id in a file with short sequences, e.g. multiFASTA or FASTQ.  Can be used to id Sequence or Kmer.
     * @details   Accessor functions are used to extract the sub fields.  bit shift and bitwise and/or are used.  they are endian-independent.
     *
     *            this keeps a 16 bit position within the sequence.   (64K length for each sequence.  identifies kmer within a sequence)
     *                       a 40 bit sequence id                     (1 TB file, == file position. identifies the sequence from which the kmer comes)
     *                       an 8 bit file id                         (256 files)
     *
     *                         |<-------pos
     *      file_id |===[=====][========*======]====================|
     *              |<---------seq_id
     *
     *
     * @note      not using union since layout of unioned struct may be architecture or compiler dependent, and we need to send this via MPI.
     *            shift operator is reliable.
     *
     *
     */
    class ShortSequenceKmerId
    {

      public:
        /// Id of the sequence.  Typically, this is the position in the file.
        size_t id;

        /**
         * @brief << operator to write out SequenceId
         * @param[in/out] ost   output stream to which the content is directed.
         * @param[in]     seq_id    sequence id object to write out
         * @return              output stream object
         */
        friend std::ostream& operator<<(std::ostream& ost, const ShortSequenceKmerId & seq_id)
        {
          ost << " ShortSeqId: file=" << static_cast<uint32_t>(seq_id.get_file_id()) << " id=" << seq_id.get_id() << " pos=" << seq_id.get_pos();

          return ost;
        }


        ShortSequenceKmerId() = default;

        ShortSequenceKmerId(size_t const & file_pos, uint8_t file_id = 0, uint16_t const & pos_in_seq = 0 ) :
          id(((file_pos & 0x000000FFFFFFFFFF) << 16) | (static_cast<size_t>(file_id) << 56) | static_cast<size_t>(pos_in_seq) ) {}
        ShortSequenceKmerId(SequenceId const & other) : ShortSequenceKmerId(other.pos_in_file, other.file_id) {}
        ShortSequenceKmerId(ShortSequenceKmerId const & other) : id(other.id) {}

        ShortSequenceKmerId& operator=(SequenceId const & other) {
          this->id = ((other.pos_in_file & 0x000000FFFFFFFFFF) << 16) | (static_cast<size_t>(other.file_id) << 56) ;
          return *this;
        }
        ShortSequenceKmerId& operator=(ShortSequenceKmerId const & other) {
          this->id = other.id;
          return *this;
        }

        bool operator==(ShortSequenceKmerId const & other) const {
          return id == other.id;
        }

        bool operator>(ShortSequenceKmerId const & other) const {
          return id > other.id;
        }

        bool operator<(ShortSequenceKmerId const & other) const {
          return id < other.id;
        }


        void operator+=(size_t dist) {
          if (((id & 0xFFFF) + dist) > 0xFFFF) throw std::invalid_argument("ShortSequenceKmerId increment overflow.");
          id += dist;
        }

        void operator-=(size_t dist) {
          if ((id & 0xFFFF) < dist)  throw std::invalid_argument("ShortSequenceKmerId decrement underflow.");
          id -= dist;
        }

        std::ptrdiff_t operator-(ShortSequenceKmerId const & other) {
          return (id >= other.id) ? static_cast<std::ptrdiff_t>(id - other.id) :
              -(static_cast<std::ptrdiff_t>(other.id - id));
        }


        /// getter for sequence id  (here, it's a 40 bit unique id == start of seq)
        size_t get_id() const { return (this->id >> 16) & 0x000000FFFFFFFFFF; }

        /// get position in file
        size_t get_pos() const { return get_id() + (this->id & 0xFFFF) ; }

        /// getter for file id
        uint8_t get_file_id() const { return id >> 56; }

    };


    /**
     * @class     bliss::common::LongSequenceKmerId
     * @brief     an Id in a file with long sequences, e.g. genomic FASTA.  Can be used to id Sequence or Kmer
     * @details   Accessor functions are used to extract the sub fields.  bit shift and bitwise and/or are used.  they are endian-independent.
     *
     *            this keeps a 40 bit position within the sequence.  (1TBase sequence, max 1TB file size.  identifies kmer position inside sequence)
     *                       a 16 bit sequence id.                   (64K sequences per file, adjacent sequences have adjacent values (2, 3, 4 ...).
     *                                                                identifies sequence (with sequence position table))
     *                       an 8 bit file id                        (256 files)
     *
     *                         |<-------pos
     *      file_id |===[=====][========*======]====================|
     *                  seq_id seq_id
     *
     * @note      not using union since layout of unioned struct may be architecture or compiler dependent, and we need to send this via MPI.
     *            shift operator is reliable.
     */
    class LongSequenceKmerId
    {
      public:
        /// Id of the sequence.  Typically, this is the position in the file.
        size_t id;

        /**
         * @brief << operator to write out SequenceId
         * @param[in/out] ost   output stream to which the content is directed.
         * @param[in]     seq_id    sequence id object to write out
         * @return              output stream object
         */
        friend std::ostream& operator<<(std::ostream& ost, const LongSequenceKmerId & seq_id)
        {
          ost << " LongSeqId: file=" << static_cast<uint32_t>(seq_id.get_file_id()) << " id=" << seq_id.get_id() << " pos=" << seq_id.get_pos();

          return ost;
        }

        LongSequenceKmerId() = default;

        LongSequenceKmerId(size_t const & file_pos, uint8_t file_id = 0, uint16_t const & seq_id = 0) :
          id((file_pos & 0x000000FFFFFFFFFF) | (static_cast<size_t>(file_id) << 56) | (static_cast<size_t>(seq_id) << 40)) {}
        LongSequenceKmerId(SequenceId const & other) : LongSequenceKmerId(other.pos_in_file, other.file_id, other.seq_id) {}
        LongSequenceKmerId(LongSequenceKmerId const & other) : id(other.id) {}

        LongSequenceKmerId& operator=(SequenceId const & other) {
          this->id = (other.pos_in_file & 0x000000FFFFFFFFFF) | (static_cast<size_t>(other.file_id) << 56) | (static_cast<size_t>(other.seq_id) << 40);
          return *this;
        }
        LongSequenceKmerId& operator=(LongSequenceKmerId const & other) {
          this->id = other.id;
          return *this;
        }

        bool operator==(LongSequenceKmerId const & other) const {
          return id == other.id;
        }

        bool operator>(LongSequenceKmerId const & other) const {
          return id > other.id;
        }

        bool operator<(LongSequenceKmerId const & other) const {
          return id < other.id;
        }


        void operator+=(size_t dist) {
          if (((id & 0xFFFFFFFFFF) + dist) > 0xFFFFFFFFFF) throw std::invalid_argument("LongSequenceKmerId increment overflow.  please check dist parameter");
          id += dist;
        }

        void operator-=(size_t dist) {
          if ((id & 0xFFFFFFFFFF) < dist)  throw std::invalid_argument("LongSequenceKmerId decrement underflow.  please check dist parameter");
          id -= dist;
        }

        std::ptrdiff_t operator-(LongSequenceKmerId const & other) {
          return (id >= other.id) ? static_cast<std::ptrdiff_t>(id - other.id) :
              -(static_cast<std::ptrdiff_t>(other.id - id));
        }

        /// getter for sequence id  (here, it is a 16 bit unique id)
        size_t get_id() const { return (this->id >> 40) & 0xFFFF; }

        /// get position within file for current kmer
        size_t get_pos() const { return (this->id & 0x000000FFFFFFFFFF); }

        /// getter for file id
        uint8_t get_file_id() const { return id >> 56; }

    };


    /**
     * @class     bliss::io::Sequence
     * @brief     represents a biological sequence, and provides iterators for traversing the sequence.
     * @details   internally we store 2 iterators, start and end, instead of start + length.
     *            reason is that using length is best for random access iterators and we don't have guarantee of that.
     *
     *
     *                         id  seq_begin   seq_end
     *              |===[=====][==={===========}===]================|
     *                         |<-- record_size -->|
     *                         |<->| seq_begin_offset
     *
     * @tparam Iterator   allows walking through the sequence data.
     */
    template<typename Iterator, typename SeqId=SequenceId>
    class Sequence
    {
      public:
        /// Iterator type for traversing the sequence.
        typedef Iterator IteratorType;
        /// type for the id struct/union
        typedef SeqId IdType;

        /// id of the sequence:  correspond to sequence_iter.  file_pos field points to beginning of record.
        mutable IdType id;

        /// length of full record.
        mutable size_t record_size;

        /// offset of seq in entire record
        mutable size_t seq_offset;

        /// offset of seq_begin from the beginning of the record.
        /// NOTE THAT IF SEQUENCE IS TRUNCATED, seq_begin may start later than the full sequence's starting point (?)
        mutable size_t seq_begin_offset;

        /// begin iterator for the sequence (may not be at record's starting position)
        mutable Iterator seq_begin;
        /// end iterator for the sequence. (may not be at record's end position)
        mutable Iterator seq_end;

        /**
         * @brief << operator to write out SequenceId
         * @param[in/out] ost   output stream to which the content is directed.
         * @param[in]     seq_id    sequence id object to write out
         * @return              output stream object
         */
        friend std::ostream& operator<<(std::ostream& ost, const Sequence & seq)
        {
          ost << " Sequence: id=[" << seq.id << "] record_size=" << seq.record_size << " seq_offset=" << seq.seq_offset <<
              " seq_begin_offset=" << seq.seq_begin_offset << " size = " << std::distance(seq.seq_begin, seq.seq_end);
          return ost;
        }


        Sequence() = default;  // needed for copy

        Sequence(IdType const & _id, size_t const & _record_size, size_t const& _seq_offset, size_t const& _seq_begin_offset, Iterator const & _begin, Iterator const & _end) :
          id(_id), record_size(_record_size), seq_offset(_seq_offset), seq_begin_offset(_seq_begin_offset), seq_begin(_begin), seq_end(_end) {}

        Sequence(Sequence const & other) :
          id(other.id), record_size(other.record_size), seq_offset(other.seq_offset), seq_begin_offset(other.seq_begin_offset), seq_begin(other.seq_begin), seq_end(other.seq_end) {}

        Sequence& operator=(Sequence const & other) {
          id = other.id;
          record_size = other.record_size;
          seq_offset = other.seq_offset;
          seq_begin_offset = other.seq_begin_offset;
          seq_begin = other.seq_begin;
          seq_end = other.seq_end;

          return *this;
        }

        virtual ~Sequence() {};

        static constexpr bool has_quality() { return false; }

        size_t seq_size() const {
          return ::std::distance(seq_begin, seq_end);
        }

        size_t seq_global_offset() const  {
          return this->id.get_pos() + seq_begin_offset;
        }

        // indicate if this record is truncated.  never happens for FASTQ
        bool is_record_truncated_at_begin() const {
        	// current sequence's data offset relative to beginning of record, is greater than the record's first sequence data position.
        	return seq_begin_offset > seq_offset;
        }
        bool is_record_truncated_at_end() const {
        	// end of sequence data in record on this node relative to the start of record,
        	// is smaller than total record size (max distance from record start in this record).
        	return (seq_begin_offset + seq_size()) < record_size;
        }

    };

  } /* namespace common */
} /* namespace bliss */

#endif /* SEQUENCE_HPP_ */
