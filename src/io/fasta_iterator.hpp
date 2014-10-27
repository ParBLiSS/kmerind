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
    template<typename Iterator, typename Alphabet, typename Kmer>
    class FASTAParser
    {
      protected:
        /// alias for this type, for use inside this class.
        using type = FASTQParser<Iterator, Alphabet, Kmer>;

        //Counting iterator
        typedef boost::counting_iterator<uint32_t> 
          count_iterator;

        //Tranform iterator on top of Counting iterator
        typedef boost::transform_iterator <std::function<uint32_t(const uint32_t&)> , count_iterator> 
          offset_transform_iterator;

        //K-mer generating iterator on Sequence data
        typedef bliss::KmerGenerationIterator<Iterator, Kmer> 
          kmer_gen_iterator;

        //Tuple of two iterators id and kmers
        typedef boost::tuple<offset_transform_iterator,kmer_gen_iterator> 
          the_iterator_tuple;

        //Zip iterator to zip the tuple 
        typedef boost::zip_iterator<the_iterator_tuple> 
          the_zip_iterator;

        //Filter iterator over the zip iterator which gives the fasta iterator 
        typedef boost::filter_iterator <std::function <bool(the_zip_iterator&) >, the_zip_iterator>
          fasta_iter;

        //Vector containing positions of fasta record start
        typedef std::vector <pair<std::size_t, std::size_t>> vectorType;

        /// iterator for the vector 
        typedef vectorType::iterator vectorIteratorType;

        /// iterator at the current starting point in the input data from where the current FASTQ record was parsed.
        Iterator _curr;

        /// iterator pointing to the end of the input data, not to go beyond.
        Iterator _end;

        /// vector of pairs that has the positions of fasta sequence headers
        vectorType _localStartLocStore; 

        // Iterator to iterate over the vector of tuples
        vectorIteratorType _localStartLocStoreIter;


        /**
         * @brief     the offset from the start of the file.  This is used as sequenceId.
         * @details   for each file with file id (fid) between 0 and 255, the seq id starts at (fid << 40).
         *            if base iterator is pointer, then units are in bytes.
         */
        size_t _offset;


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
            const size_t &offset, vectorType &localStartLocStore)
            : _curr(start), _end(end), _offset(offset), _localStartLocStore(localStartLocStore), _localStartLocStoreIter(_localStartLocStore.begin())
        {
        }


        /**
         * @brief   construct the begin FASTA iterator
         * @return  begin FASTA iterator 
         */
        fasta_iter begin(){
            fasta_iter it_begin(

                //Filter predicate, to discard invalid kmers marked by some id
                [](const the_zip_iterator& e) { return *(e.first) < std::numeric_limits<uint32_t>::max();},

                //filter iterator's start range
                the_zip_iterator(
                  the_iterator_tuple(
                    offset_transform_iterator(
                      count_iterator(_offset),
                      //Transform counter to ids
                      returnOffset<uint32_t>()
                      )
                    kmer_gen_iterator(
                      _curr,
                      true
                      )
                    )
                  ),

                //filter iterator's end range
                the_zip_iterator(
                  the_iterator_tuple(
                    offset_transform_iterator(
                      //TODO : Need to correct the end of count_iterator
                      count_iterator(_offset),
                      //Transform counter to ids
                      returnOffset<uint32_t>()
                      )
                    kmer_gen_iterator(
                      _end,
                      true
                      )
                    )
                  )
                );
            return it_begin;
        }

        /**
         * @brief   construct the end of FASTA iterator
         * @return  end FASTA iterator 
         */
        fasta_iter end(){
            fasta_iter it_end(

                //Filter predicate, to discard invalid kmers marked by some id
                [](const the_zip_iterator& e) { return *(e.first) < std::numeric_limits<uint32_t>::max();},

                //filter iterator's start range
                the_zip_iterator(
                  the_iterator_tuple(
                    offset_transform_iterator(
                      count_iterator(_offset),
                      returnOffset<uint32_t>()
                      )
                    kmer_gen_iterator(
                      _end,
                      false
                      )
                    )
                  ),

                //filter iterator's end range
                the_zip_iterator(
                  the_iterator_tuple(
                    offset_transform_iterator(
                      count_iterator(_offset),
                      returnOffset<uint32_t>()
                      )
                    kmer_gen_iterator(
                      _end,
                      false
                      )
                    )
                  )
                );
            return it_end;
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

          //FASTA sequence header lies before the offset
          if(leftIndex < rawOffset && rightIndex < rawOffset)
            _localStartLocStoreIter ++;

          //Next FASTA sequence header lies after the offset
          if(leftIndex > rawOffset && rightIndex > rawOffset)
            return rawOffset;
          //The offset points doesn't point to sequence data but the header
          else if(leftIndex < rawOffset && rightIndex > rawOffset)
            return std::numeric_limits<eleType>::max();
        }
        
    }
}
