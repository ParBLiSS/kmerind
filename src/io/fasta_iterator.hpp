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
#include <common/kmer_iterators.hpp>


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


    /**
     * @class     bliss::io::FASTAParser
     * @brief     The purpose of this class is to expose an iterator to parse 
     *            FASTA file. The output of the iterator will be k-mer and its associated
     *            sequence id
     *
     * @details   FASTA iterator is build from the available iterators in bliss 
     *            1) K-mer iterator over raw data
     *            2) id-iterator over counting iterator
     *            3) Zip and filter to generate iterator of the pairs (Id, Kmer)
     *
     *            This class requires a vector of positions containing the positions of all the 
     *            FASTA sequence headers relevant to the local block
     *
     * @tparam    Iterator    The underlying iterator to be traversed to generate a Sequence
     * @tparam    Alphabet    allows interpretation of the sequence
     * @tparam    Kmer        bliss::Kmer type
     */
    template<typename Iterator, typename Kmer>
    class FASTAParser
    {
      protected:

        //Type for the offset
        typedef uint32_t offSetType;

        //Counting iterator
        typedef boost::counting_iterator<offSetType> 
          count_iterator;

        //Tranform iterator on top of Counting iterator
        typedef boost::transform_iterator <std::function<offSetType(const offSetType&)> , count_iterator> 
          offset_transform_iterator;

        //Tuple of two iterators id and sequence
        typedef boost::tuple<offset_transform_iterator, iterator> 
          the_offset_sequence_tuple;

        //Zip iterator to zip the tuple of sequence and id
        typedef boost::zip_iterator<the_offset_sequence_tuple> 
          the_offset_sequence_zip;

        //Filter iterator over the zip iterator 
        typedef boost::filter_iterator <std::function <bool(the_zip_iterator&) >, the_zip_iterator>
          corrected_offset_sequence_zip;

        //K-mer generating iterator on Sequence data
        typedef bliss::KmerGenerationIterator<Iterator, Kmer> 
          kmer_gen_iterator;

        //Tuple of two iterators id and sequence
        typedef boost::tuple<offset_transform_iterator, kmer_gen_iterator> 
          the_offset_kmer_tuple;

        //Zip iterator to zip the tuple of kmer and id
        typedef boost::zip_iterator<the_offset_kmer_tuple> 
          the_offset_kmer_zip;

        //Vector containing positions of fasta record start
        typedef FASTALoader::vectorType vectorType;

        /// iterator for the vector 
        typedef vectorType::iterator vectorIteratorType;

        /// iterator at the current starting point in the input data from where the current FASTQ record was parsed.
        Iterator _curr;

        /// iterator pointing to the end of the input data, not to go beyond.
        Iterator _end;

        /// vector of pairs that has the positions of fasta sequence headers
        vectorType _localStartLocStore; 

        /**
         * @brief     the offset of start and end of the search range.  This is used for generating the sequenceId.
         */
        RangeType _offsetRange;

        // Iterator to iterate over the vector of tuples
        vectorIteratorType _localStartLocStoreIter;

        //Define EOL type
        static const typename std::iterator_traits<Iterator>::value_type eol = '\n';

        /// alias for this type, for use inside this class.
        using type = FASTQParser<Iterator, Kmer>;

        /**
         * @brief     Filter predicate, to discard the invalid entries in the zipped offset and sequence iterators
         * @details   sequence is not valid if its EOL, offset is invalid if it has been set 
         *            to the max of the underlining data type
         * @param e   zip iterator to filter
         * @return    true if both offset and sequence is valid, false otherwise
         */
        auto filterFunctor = [](const the_offset_sequence_zip &e) { 
          bool b1 = *(e.first) < std::numeric_limits<offSetType>::max();
          bool b2 = *(e.second) != type::eol;
          return b1 && b2;
        }

      public:

        /**
         * @brief constructor, initializes with a start and end.  this represents the start of the output iterator
         *
         * @param start                 beginning of the data to be parsed
         * @param end                   end of the data to be parsed.
         * @param offset                position number relative to the file being processed.
         * @param localStartLocStore    vector of tuples that stores fasta sequence header's begin and index offset 
         */
        FASTAParser(const Iterator& start, const Iterator& end, 
            const RangeType &offsetRange, const vectorType &localStartLocStore)
            : _curr(start), _end(end), _offsetRange(offsetRange), _localStartLocStore(localStartLocStore), _localStartLocStoreIter(_localStartLocStore.begin())
        {
        }


        /**
         * @brief   construct the begin FASTA iterator
         * @return  start FASTA iterator as zip iterator (id, kmer)
         */
        the_offset_kmer_zip begin(){

          //Initialise the begin and end of the base sequence iterator
          Iterator sequenceBegin = _curr;
          Iterator sequenceEnd = _end;

          //Initialise the begin and end of id iterator using count_iterator 
          count_iterator countBegin = count_iterator(_offsetRange.start);
          count_iterator countEnd = count_iterator(offsetRange_.end);

          //Transform counter elements to ids by using transform iterator
          offset_transform_iterator idBeginCorrect = offset_transform_iterator(countBegin, returnOffset<offSetType>());
          offset_transform_iterator idEndCorrect = offset_transform_iterator(countEnd, returnOffset<offSetType>());

          //Construct the tuple of ids and sequence iterators
          the_offset_sequence_tuple idSeqStart = the_offset_sequence_tuple(idBeginCorrect, sequenceBegin);
          the_offset_sequence_tuple idSeqEnd = the_offset_sequence_tuple(idEndCorrect, sequenceEnd);

          //Zip the tuple 
          the_offset_sequence_zip idSeqZipStart = the_offset_sequence_zip(idSeqStart);
          the_offset_sequence_zip idSeqZipEnd = the_offset_sequence_zip(idSeqEnd);

          //Filter the zip iterator using filter iterator; need both start and end iterators for the constructor
          corrected_offset_sequence_zip idSeqZipStartCorrect = corrected_offset_sequence_zip(filterFunctor, idSeqZipStart, idSeqZipEnd);

          //"Unzip the iterator" and construct the kmer iterator from sequence iterator 
          kmer_gen_iterator KmerIteratorStart = kmer_gen_iterator(idSeqZipStartCorrect.second, true);;

          //"Rezip the iterators" to pair the id and kmer iterators
          the_offset_kmer_tuple idKmerStart = the_offset_kmer_tuple(idBeginCorrect, KmerIteratorStart);
          the_offset_kmer_zip finalIdKmerZipStart = the_offset_kmer_zip(idKmerStart);

          return finalIdKmerZipStart;
        }

        /**
         * @brief   construct the end of FASTA iterator
         * @return  end FASTA iterator as zip iterator (id, kmer)
         */
        the_offset_kmer_zip end(){

          //Initialise the end of the base sequence iterator
          Iterator sequenceEnd = _end;

          //Initialise the end of id iterator using count_iterator 
          count_iterator countEnd = count_iterator(offsetRange_.end);

          //Transform counter elements to ids by using transform iterator
          offset_transform_iterator idEndCorrect = offset_transform_iterator(countEnd, returnOffset<offSetType>());

          //Construct the tuple of ids and sequence iterators
          the_offset_sequence_tuple idSeqEnd = the_offset_sequence_tuple(idEndCorrect, sequenceEnd);

          //Zip the tuple 
          the_offset_sequence_zip idSeqZipEnd = the_offset_sequence_zip(idSeqEnd);

          //Filter the zip iterator using filter iterator; need both start and end iterators for the constructor
          corrected_offset_sequence_zip idSeqZipEndCorrect = corrected_offset_sequence_zip(filterFunctor, idSeqZipEnd, idSeqZipEnd);

          //"Unzip the iterator" and construct the kmer iterator from sequence iterator 
          kmer_gen_iterator KmerIteratorEnd = kmer_gen_iterator(idSeqZipEnd.second, true);;

          //"Rezip the iterators" to pair the id and kmer iterators
          the_offset_kmer_tuple idKmerEnd = the_offset_kmer_tuple(idEndCorrect, KmerIteratorEnd);
          the_offset_kmer_zip finalIdKmerZipEnd = the_offset_kmer_zip(idKmerEnd);

          return finalIdKmerZipEnd;
        }

      protected:

        /**
         * @brief   convert raw counting iterator offset to actual offset by parsing tuple vector
         * @return  converted offset 
         * @param   rawOffset the offset value which needs to be transformed
         * @tparam  eleType   the type definition of offset 
         */
        template<typename eleType>
        eleType returnOffset(const eleType& rawOffset)
        {
          //Advance the iterator so that the fasta header overlaps with the offset or
          //it lies ahead of the current offset
          
          leftIndex = (*_localStartLocStoreIter).first;
          rightIndex = (*_localStartLocStoreIter).second;

          //if FASTA sequence header lies before the offset
          while(rightIndex < rawOffset)
            _localStartLocStoreIter ++;

          //Next FASTA sequence header lies after the offset
          if(leftIndex > rawOffset)
            return rawOffset;
          //The offset points doesn't point to sequence data but the header
          //(leftIndex <= rawOffset && rightIndex >= rawOffset)
          else
            return std::numeric_limits<eleType>::max();
        }
    }
}
