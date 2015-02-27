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
#include <partition/range.hpp>
#include <common/kmer_iterators.hpp>
#include <common/alphabets.cpp>
#include <io/fasta_loader.hpp>

// boost iterators used
#include <boost/iterator/counting_iterator.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/iterator/zip_iterator.hpp>
#include <boost/iterator/filter_iterator.hpp>


namespace bliss
{
  namespace io
  {
    //dummy struct to indicate FASTA format
    struct FASTA {};

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
     * @tparam    Iterator      The underlying iterator to be traversed to generate a Sequence
     * @tparam    Kmer          bliss::Kmer type
     * @tparam    AlphabetType  allows interpretation of the sequence (Example : DNA, DNA5)
     */
    template<typename Iterator, typename KmerType>
    class FASTAParser
    {
      protected:
        //Type for the offset
        typedef typename FASTALoader<Iterator,KmerType>::offSetType offSetType;

        //Counting iterator
        typedef typename boost::counting_iterator<offSetType> 
          count_iterator;

        //Tranform iterator on top of Counting iterator
        typedef boost::transform_iterator <std::function<offSetType(const offSetType&)> , count_iterator> 
          offset_transform_iterator;

        //Tranform iterator on top of sequence iterator
        typedef typename std::iterator_traits<Iterator>::value_type Iterator_valueType;
        typedef boost::transform_iterator <std::function<Iterator_valueType(const Iterator_valueType&)> , Iterator> 
          transformed_seq_iterator;

        //Converting raw character Iterator to Kmer Iterator (filtering later)
        typedef bliss::common::KmerGenerationIterator<transformed_seq_iterator, KmerType> KmerIncompleteIterator;

        //Tuple of two iterators id and Kmer iterator over raw data
        typedef boost::tuple<offset_transform_iterator, KmerIncompleteIterator> 
          the_offset_rawkmer_tuple;

        //Zip iterator to zip the tuple of id and raw Kmers
        typedef boost::zip_iterator<the_offset_rawkmer_tuple> 
          the_offset_rawkmer_zip;

        //Filter iterator over the zip iterator 
        typedef typename std::iterator_traits<the_offset_rawkmer_zip>::value_type the_offset_rawkmer_zip_valueType;

      public:
        //Corrected final zip iterator of Kmer and Offset value
        typedef boost::filter_iterator <std::function <bool(the_offset_rawkmer_zip_valueType) >, the_offset_rawkmer_zip>
          corrected_offset_Kmer_zip;

      protected:

        //Vector containing positions of fasta record start
        typedef typename FASTALoader<Iterator,KmerType>::vectorType vectorType;

        /// iterator for the vector 
        typedef typename vectorType::iterator vectorIteratorType;

        //Type for specifying file range
        typedef typename FASTALoader<Iterator,KmerType>::RangeType RangeType;

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

        /**
         * @brief     the offset of start and end of the complete file range.  This is used for generating the sequenceId.
         */
        RangeType _parentRange;

        // Iterator to iterate over the vector of tuples
        vectorIteratorType _localStartLocStoreBeginIter;
        vectorIteratorType _localStartLocStoreEndIter;

        //Define EOL type
        static constexpr typename std::iterator_traits<Iterator>::value_type eol = '\n';

        /**
         * @brief     Filter predicate, to discard the invalid entries in the zipped offset and sequence iterators
         * @details   sequence is not valid if its EOL, offset is invalid if it has been set 
         *            to the max of the underlining data type
         * @param e   zip iterator to filter
         * @return    true if both offset and sequence is valid, false otherwise
         */
        typedef bool (*filterFunctorType)(the_offset_rawkmer_zip_valueType);
        
        filterFunctorType filterFunctor = [](the_offset_rawkmer_zip_valueType e) -> bool { 
          bool b1 = boost::get<0>(e) < std::numeric_limits<offSetType>::max();
          //bool b2 = boost::get<1>(e) != eol;
          return b1;
        };

      public:

        /**
         * @brief constructor, initializes with a start and end.  this represents the start of the output iterator
         *
         * @param start                 beginning of the data to be parsed
         * @param end                   end of the data to be parsed.
         * @param searchRange           offset range local to this execution
         * @param parentRange           offset range for the whole file
         * @param localStartLocStore    vector of tuples that stores fasta sequence header's begin and index offset 
         */
        FASTAParser(const Iterator& start, const Iterator& end, 
            const RangeType &searchRange, const RangeType &parentRange, const vectorType &localStartLocStore)
        {
          _curr = start;
          _end = end;
          _offsetRange = searchRange;
          _parentRange = parentRange;
          _localStartLocStore = localStartLocStore;
          _localStartLocStoreBeginIter = _localStartLocStore.begin();
          _localStartLocStoreEndIter = _localStartLocStore.end();

          //For DEBUG
          //std::cout << "Printing start tuple values \n" ;
          //for(auto iter=_localStartLocStore.begin(); iter!= _localStartLocStore.end(); iter++)
            //std::cout << (*iter).first << ", " << (*iter).second << "\n";
        }

        FASTAParser() = delete; 

        /**
         * @brief   construct the begin FASTA iterator
         * @return  start FASTA iterator as zip iterator (id, kmer)
         */
        corrected_offset_Kmer_zip begin(){
          //Initialise the begin and end of the base sequence iterator
          Iterator sequenceBegin = _curr;
          Iterator sequenceEnd = _end;

          //Initialise the begin and end of id iterator using count_iterator 
          count_iterator countBegin = count_iterator(_offsetRange.start);
          count_iterator countEnd = count_iterator(_offsetRange.end);

          //Transform counter elements to ids by using transform iterator
          //Functor for the tranformation
          returnOffset<offSetType> transformFunctor(_localStartLocStoreBeginIter, _localStartLocStoreEndIter, _offsetRange, _parentRange);

          offset_transform_iterator idBeginCorrect = offset_transform_iterator(countBegin, transformFunctor);
          offset_transform_iterator idEndCorrect = offset_transform_iterator(countEnd, transformFunctor);

          //Converted sequence iterator from ascii to concise format 
          //Functor for the tranformation
          returnSeqVal seqtransformFunctor;
          transformed_seq_iterator trSequenceBegin = transformed_seq_iterator(sequenceBegin, seqtransformFunctor);
          transformed_seq_iterator trSequenceEnd = transformed_seq_iterator(sequenceEnd, seqtransformFunctor);

          //Initialise the begin and end of raw kmer iterator
          KmerIncompleteIterator kmer_raw_begin_iterator = KmerIncompleteIterator(trSequenceBegin, true);
          KmerIncompleteIterator kmer_raw_end_iterator = KmerIncompleteIterator(trSequenceEnd, true);

          //Construct the tuple of ids and kmer iterators
          the_offset_rawkmer_tuple idKmerStart = the_offset_rawkmer_tuple(idBeginCorrect, kmer_raw_begin_iterator);
          the_offset_rawkmer_tuple idKmerEnd = the_offset_rawkmer_tuple(idEndCorrect, kmer_raw_end_iterator);

          //Zip the tuple 
          the_offset_rawkmer_zip idKmerZipStart = the_offset_rawkmer_zip(idKmerStart);
          the_offset_rawkmer_zip idKmerZipEnd = the_offset_rawkmer_zip(idKmerEnd);

          //Filter the zip iterator using filter iterator; need both start and end iterators for the constructor
          corrected_offset_Kmer_zip idSeqZipStartCorrect = corrected_offset_Kmer_zip(filterFunctor, idKmerZipStart, idKmerZipEnd);

          return idSeqZipStartCorrect;
        }

        /**
         * @brief   construct the end of FASTA iterator
         * @return  end FASTA iterator as zip iterator (id, kmer)
         */
        corrected_offset_Kmer_zip end(){

          //Initialise the end of the base sequence iterator
          Iterator sequenceEnd = _end;

          //Initialise and end of id iterator using count_iterator 
          count_iterator countEnd = count_iterator(_offsetRange.end);

          //Transform counter elements to ids by using transform iterator
          //Functor for the tranformation
          returnOffset<offSetType> transformFunctor(_localStartLocStoreBeginIter, _localStartLocStoreEndIter, _offsetRange, _parentRange);
          offset_transform_iterator idEndCorrect = offset_transform_iterator(countEnd, transformFunctor);

          //Converted sequence iterator from ascii to concise format 
          //Functor for the tranformation
          returnSeqVal seqtransformFunctor;
          transformed_seq_iterator trSequenceEnd = transformed_seq_iterator(sequenceEnd, seqtransformFunctor);

          //Initialise the begin and end of raw kmer iterator
          KmerIncompleteIterator kmer_raw_end_iterator = KmerIncompleteIterator(trSequenceEnd, true);

          //Construct the tuple of ids and kmer iterators
          the_offset_rawkmer_tuple idKmerEnd = the_offset_rawkmer_tuple(idEndCorrect, kmer_raw_end_iterator);

          //Zip the tuple 
          the_offset_rawkmer_zip idKmerZipEnd = the_offset_rawkmer_zip(idKmerEnd);

          //Filter the zip iterator using filter iterator; need both start and end iterators for the constructor
          corrected_offset_Kmer_zip idSeqZipEndCorrect = corrected_offset_Kmer_zip(filterFunctor, idKmerZipEnd, idKmerZipEnd);

          return idSeqZipEndCorrect;
        }

        offSetType getOffset(const corrected_offset_Kmer_zip &e)
        {
          return boost::get<0>(*e);
        }

        KmerType getKmer(const corrected_offset_Kmer_zip &e)
        {
          return boost::get<1>(*e);
        }

      protected:

        /**
         * @brief   convert raw counting iterator offset to actual offset by parsing tuple vector
         * @return  converted offset 
         * @param   rawOffset the offset value which needs to be transformed
         * @tparam  eleType   the type definition of offset 
         */
        template<typename eleType>
        struct returnOffset
        {
          vectorIteratorType localStartLocStoreEndIter;
          vectorIteratorType localStartLocStoreIter;
          eleType leftIndex;
          eleType rightIndex;
          RangeType offsetRange; 
          RangeType parentRange;

          //Constructor (mandatory to use this)
          returnOffset(vectorIteratorType &_localStartLocStoreBeginIter, vectorIteratorType &_localStartLocStoreEndIter, const RangeType &_offsetRange, const RangeType &_parentRange)
          {
            localStartLocStoreIter = _localStartLocStoreBeginIter;
            localStartLocStoreEndIter = _localStartLocStoreEndIter;
            leftIndex = std::numeric_limits<eleType>::min();
            rightIndex = std::numeric_limits<eleType>::min();
            offsetRange = _offsetRange;
            parentRange = _parentRange;
          }

          //Transform functor
          eleType operator()(const eleType& rawOffset)
          {
            //Advance the iterator so that the fasta header overlaps with the offset or
            //it lies ahead of the current offset

            //if FASTA sequence header lies before the offset
            while(rightIndex < rawOffset && localStartLocStoreIter != localStartLocStoreEndIter)
            {
              leftIndex = (*localStartLocStoreIter).first;
              rightIndex = (*localStartLocStoreIter).second;
              localStartLocStoreIter ++;

              //A sanity check (important in multithreaded version because here we might have headers ahead that we should ignore)
              leftIndex = std::min(offsetRange.end + KmerType::size, leftIndex);
              rightIndex = std::min(offsetRange.end + KmerType::size, rightIndex);
            }

            //If there is no FASTA sequence header ahead, we assume left and right indices at the end.
            //This makes sure that we avoid reading incomplete kmers near the end of the file(or buffer)
            if(rawOffset > rightIndex)
            {
              //Assume a dummy header (This takes care of the extra characters we might need from next partition)
              //leftIndex should point to two extra offsets from the last alphabet we are interested in
              //TODO: What if we encounter a header in the overlap?
              leftIndex = std::min(offsetRange.end + KmerType::size, parentRange.end);
              rightIndex = leftIndex + 1;
            }

            //Current offset should not overlap with the header
            //Also make sure we have atleast KMER_SIZE characters ahead (also consider eol after them)
            if(rawOffset < leftIndex && leftIndex - rawOffset > KmerType::size)
            {
              return rawOffset;
            }

            //The offset points doesn't point to sequence data but the header
            else
            {
              return std::numeric_limits<eleType>::max();
            }
            return rawOffset; 
          }
        };

        /**
         * @brief   converting sequence values to a more packed representation 
         * @return  converted value 
         * @param   rawOffset the offset value which needs to be transformed
         * @tparam  eleType   the type definition of offset 
         */
        struct returnSeqVal 
        {
          Iterator_valueType operator()(const Iterator_valueType& inValue)
          {
            return KmerType::KmerAlphabet::FROM_ASCII[static_cast<size_t>(inValue)];
          }
        };
    };
  }
}
#endif
