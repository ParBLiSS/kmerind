/**
 * @file    test_threads.cpp
 * @ingroup
 * @author  Tony Pan <tpan7@gatech.edu>
 * @brief
 * @details
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#include "config.hpp"

#include <unistd.h>  // get hostname


#include <functional>
#include "utils/logging.h"

#include "common/alphabets.hpp"


#include "index/kmer_index.hpp"
#include "common/kmer.hpp"
#include "common/base_types.hpp"
#include <string>
#include <sstream>
#include "utils/kmer_utils.hpp"
#include <chrono>

#include <mxx/collective.hpp>


namespace mxx {

  template<unsigned int size, typename A, typename WT>
  class datatype<typename bliss::common::Kmer<size, A, WT> > :
    public datatype_contiguous<typename bliss::common::Kmer<size, A, WT>::KmerWordType,
      bliss::common::Kmer<size, A, WT>::nWords> {};

  template<>
  class datatype<bliss::io::FASTQ::SequenceId > :
    public datatype_contiguous<decltype(bliss::io::FASTQ::SequenceId::file_pos),1> {};

}

/*
 * TYPE DEFINITIONS
 */

void buildPosIndex(const std::string & filename, MPI_Comm comm) {
  // ==  open file
//            distributed_file df;
//            range = df.open(filename, communicator);  // memmap internally
//            data = df.data();

  using Alphabet = bliss::common::DNA;
  using FileLoaderType = bliss::io::FASTQLoader<CharType, true, false>; // raw data type :  use CharType
  using KmerType = bliss::common::Kmer<21, Alphabet, WordType>;
  using IdType = bliss::io::FASTQ::SequenceId;
  using TupleType = std::pair<KmerType, IdType>;
  using MapType = std::unordered_multimap<KmerType, IdType, bliss::hash::farm::KmerSuffixHash<KmerType> >;

  //====  now process the file, one L1 block (block partition by MPI Rank) at a time
  // from FileLoader type, get the block iter type and range type
  using FileBlockIterType = typename FileLoaderType::L1BlockType::iterator;

  using ParserType = bliss::io::FASTQParser<FileBlockIterType, void>;
  using SeqType = typename ParserType::SequenceType;
  using SeqIterType = bliss::io::SequencesIterator<ParserType>;


  /// converter from ascii to alphabet values
  using BaseCharIterator = bliss::iterator::transform_iterator<typename SeqType::IteratorType, bliss::common::ASCII2<Alphabet> >;

  /// kmer generation iterator
  using KmerIterType = bliss::common::KmerGenerationIterator<BaseCharIterator, KmerType>;
  /// kmer position iterator type
  using IdIterType = bliss::iterator::SequenceIdIterator<IdType>;

  /// combine kmer iterator and position iterator to create an index iterator type.
  using KmerIndexIterType = bliss::iterator::ZipIterator<KmerIterType, IdIterType>;

  int commSize, commRank;
  MPI_Comm_size(comm, &commSize);
  MPI_Comm_size(comm, &commRank);

  {
    std::chrono::high_resolution_clock::time_point t1, t2;
    std::chrono::duration<double> time_span;

    t1 = std::chrono::high_resolution_clock::now();
    //==== create file Loader
    FileLoaderType loader(comm, filename, 1, sysconf(_SC_PAGE_SIZE));  // this handle is alive through the entire building process.
    typename FileLoaderType::L1BlockType partition = loader.getNextL1Block();

    t2 = std::chrono::high_resolution_clock::now();
     time_span =
         std::chrono::duration_cast<std::chrono::duration<double>>(
             t2 - t1);
     INFOF("R %d file open time: %f", commRank, time_span.count());

     //== reserve
     t1 = std::chrono::high_resolution_clock::now();
     MapType map;
    // modifying the local index directly here causes a thread safety issue, since callback thread is already running.
    // index reserve internally sends a message to itself.
     // call after getting first L1Block to ensure that file is loaded.
     size_t est_size = (loader.getKmerCountEstimate(KmerType::size) + commSize - 1) / commSize;
     map.reserve(est_size);
    t2 = std::chrono::high_resolution_clock::now();
     time_span =
         std::chrono::duration_cast<std::chrono::duration<double>>(
             t2 - t1);
     INFOF("R %d reserve time: %f", commRank, time_span.count());



     // == create kmer iterator
     //            kmer_iter start(data, range);
     //            kmer_iter end(range.second,range.second);

     t1 = std::chrono::high_resolution_clock::now();
     std::vector< TupleType > temp;
     temp.reserve(est_size);

     ParserType parser;
    //=== copy into array
    while (partition.getRange().size() > 0) {
      //== process the chunk of data
      SeqType read;

      //==  and wrap the chunk inside an iterator that emits Reads.
      SeqIterType seqs_start(parser, partition.begin(), partition.end(), partition.getRange().start);
      SeqIterType seqs_end(partition.end());


      //== loop over the reads
      for (; seqs_start != seqs_end; ++seqs_start)
      {
        // first get read
        read = *seqs_start;

        // then compute and store into index (this will generate kmers and insert into index)
        if (read.seqBegin == read.seqEnd) continue;

        //== set up the kmer generating iterators.
        KmerIterType start(BaseCharIterator(read.seqBegin, bliss::common::ASCII2<Alphabet>()), true);
        KmerIterType end(BaseCharIterator(read.seqEnd, bliss::common::ASCII2<Alphabet>()), false);

        //== set up the position iterators
        IdIterType id_start(read.id);
        IdIterType id_end(read.id);

        // ==== set up the zip iterators
        KmerIndexIterType index_start(start, id_start);
        KmerIndexIterType index_end(end, id_end);


	temp.insert(temp.end(), index_start, index_end);
//        for (auto it = index_start; it != index_end; ++it) {
//          temp.push_back(*it);
//        }
        //std::copy(index_start, index_end, temp.end());
//        INFOF("R %d inserted.  new temp size = %lu", commRank, temp.size());
      }

      partition = loader.getNextL1Block();
    }
    t2 = std::chrono::high_resolution_clock::now();
     time_span =
         std::chrono::duration_cast<std::chrono::duration<double>>(
             t2 - t1);
     INFOF("R %d local kmer array time: %f. temp size = %lu", commRank, time_span.count(), temp.size());


     // distribute
     t1 = std::chrono::high_resolution_clock::now();

     bliss::hash::farm::KmerPrefixHash<KmerType> hash(2 * ceilLog2(commSize));
	if (commSize > 1)
	     mxx::msgs_all2all(temp, [&] ( TupleType const &x) {
       		return (hash(x.first) % commSize);
     		}, comm);
     t2 = std::chrono::high_resolution_clock::now();
      time_span =
          std::chrono::duration_cast<std::chrono::duration<double>>(
              t2 - t1);
      INFOF("R %d bucket and distribute: %f. temp size = %lu", commRank, time_span.count(), temp.size());

      // == insert into distributed map.
     t1 = std::chrono::high_resolution_clock::now();
      map.insert(temp.begin(), temp.end());

     t2 = std::chrono::high_resolution_clock::now();
      time_span =
          std::chrono::duration_cast<std::chrono::duration<double>>(
              t2 - t1);
      INFOF("R %d inserted: %f. final size = %lu", commRank, time_span.count(), map.size());

//            //distributed_map m(element_count);
//            m.reserve(element_count, communicator);
//            m.insert(start, end, communicator);
//            //m.local_rehash();

  //df.close();
  }

}

void buildPosQualIndex(const std::string & filename, MPI_Comm comm) {
  // ==  open file
//            distributed_file df;
//            range = df.open(filename, communicator);  // memmap internally
//            data = df.data();

  using Alphabet = bliss::common::DNA;
  using FileLoaderType = bliss::io::FASTQLoader<CharType, true, false>; // raw data type :  use CharType
  using FileBlockIterType = typename FileLoaderType::L1BlockType::iterator;

  using KmerType = bliss::common::Kmer<21, Alphabet, WordType>;

  using IdType = bliss::io::FASTQ::SequenceId;
  using QualType = float;

  using KmerInfoType = std::pair<IdType, QualType>;
  using TupleType = std::pair<KmerType, KmerInfoType>;
  using MapType = std::unordered_multimap<KmerType, KmerInfoType, bliss::hash::farm::KmerSuffixHash<KmerType> >;

  //====  now process the file, one L1 block (block partition by MPI Rank) at a time

  using ParserType = bliss::io::FASTQParser<FileBlockIterType, QualType>;
  using SeqType = typename ParserType::SequenceType;
  using SeqIterType = bliss::io::SequencesIterator<ParserType>;



  /// converter from ascii to alphabet values
  using BaseCharIterator = bliss::iterator::transform_iterator<typename SeqType::IteratorType, bliss::common::ASCII2<Alphabet> >;

  /// kmer generation iterator
  using KmerIterType = bliss::common::KmerGenerationIterator<BaseCharIterator, KmerType>;
  /// kmer position iterator type
  using IdIterType = bliss::iterator::SequenceIdIterator<IdType>;

  using QualIterType = bliss::index::QualityScoreGenerationIterator<typename SeqType::IteratorType, 21, bliss::index::Illumina18QualityScoreCodec<QualType> >;

  /// combine kmer iterator and position iterator to create an index iterator type.
  using KmerInfoIterType = bliss::iterator::ZipIterator<IdIterType, QualIterType>;
  using KmerIndexIterType = bliss::iterator::ZipIterator<KmerIterType, KmerInfoIterType>;

  int commSize, commRank;
  MPI_Comm_size(comm, &commSize);
  MPI_Comm_size(comm, &commRank);

  {
    std::chrono::high_resolution_clock::time_point t1, t2;
    std::chrono::duration<double> time_span;

    t1 = std::chrono::high_resolution_clock::now();
    //==== create file Loader
    FileLoaderType loader(comm, filename, 1, sysconf(_SC_PAGE_SIZE));  // this handle is alive through the entire building process.
    typename FileLoaderType::L1BlockType partition = loader.getNextL1Block();

    t2 = std::chrono::high_resolution_clock::now();
     time_span =
         std::chrono::duration_cast<std::chrono::duration<double>>(
             t2 - t1);
     INFOF("R %d file open time: %f", commRank, time_span.count());

     //== reserve
     t1 = std::chrono::high_resolution_clock::now();
     MapType map;
    // modifying the local index directly here causes a thread safety issue, since callback thread is already running.
    // index reserve internally sends a message to itself.
     // call after getting first L1Block to ensure that file is loaded.
     size_t est_size = (loader.getKmerCountEstimate(KmerType::size) + commSize - 1) / commSize;
     map.reserve(est_size);
    t2 = std::chrono::high_resolution_clock::now();
     time_span =
         std::chrono::duration_cast<std::chrono::duration<double>>(
             t2 - t1);
     INFOF("R %d reserve time: %f", commRank, time_span.count());



     // == create kmer iterator
     //            kmer_iter start(data, range);
     //            kmer_iter end(range.second,range.second);

     t1 = std::chrono::high_resolution_clock::now();
     std::vector< TupleType > temp;
     temp.reserve(est_size);

     ParserType parser;
    //=== copy into array
    while (partition.getRange().size() > 0) {
      //== process the chunk of data
      SeqType read;

      //==  and wrap the chunk inside an iterator that emits Reads.
      SeqIterType seqs_start(parser, partition.begin(), partition.end(), partition.getRange().start);
      SeqIterType seqs_end(partition.end());


      //== loop over the reads
      for (; seqs_start != seqs_end; ++seqs_start)
      {
        // first get read
        read = *seqs_start;

        // then compute and store into index (this will generate kmers and insert into index)
        if (read.seqBegin == read.seqEnd || read.qualBegin == read.qualEnd) continue;

        //== set up the kmer generating iterators.
        KmerIterType start(BaseCharIterator(read.seqBegin, bliss::common::ASCII2<Alphabet>()), true);
        KmerIterType end(BaseCharIterator(read.seqEnd, bliss::common::ASCII2<Alphabet>()), false);

        //== set up the position iterators
        IdIterType id_start(read.id);
        IdIterType id_end(read.id);

	QualIterType qual_start(read.qualBegin);
	QualIterType qual_end(read.qualEnd);

	KmerInfoIterType info_start(id_start, qual_start);
	KmerInfoIterType info_end(id_end, qual_end);


        // ==== set up the zip iterators
        KmerIndexIterType index_start(start, info_start);
        KmerIndexIterType index_end(end, info_end);


	temp.insert(temp.end(), index_start, index_end);
//        for (auto it = index_start; it != index_end; ++it) {
//          temp.push_back(*it);
//        }
        //std::copy(index_start, index_end, temp.end());
//        INFOF("R %d inserted.  new temp size = %lu", commRank, temp.size());
      }

      partition = loader.getNextL1Block();
    }
    t2 = std::chrono::high_resolution_clock::now();
     time_span =
         std::chrono::duration_cast<std::chrono::duration<double>>(
             t2 - t1);
     INFOF("R %d local kmer array time: %f. temp size = %lu", commRank, time_span.count(), temp.size());


     // distribute
     t1 = std::chrono::high_resolution_clock::now();

     bliss::hash::farm::KmerPrefixHash<KmerType> hash(2 * ceilLog2(commSize));
	if (commSize > 1)
	     mxx::msgs_all2all(temp, [&] ( TupleType const &x) {
       		return (hash(x.first) % commSize);
     		}, comm);
     t2 = std::chrono::high_resolution_clock::now();
      time_span =
          std::chrono::duration_cast<std::chrono::duration<double>>(
              t2 - t1);
      INFOF("R %d bucket and distribute: %f. temp size = %lu", commRank, time_span.count(), temp.size());

      // == insert into distributed map.
     t1 = std::chrono::high_resolution_clock::now();
      map.insert(temp.begin(), temp.end());

     t2 = std::chrono::high_resolution_clock::now();
      time_span =
          std::chrono::duration_cast<std::chrono::duration<double>>(
              t2 - t1);
      INFOF("R %d inserted: %f. final size = %lu", commRank, time_span.count(), map.size());

//            //distributed_map m(element_count);
//            m.reserve(element_count, communicator);
//            m.insert(start, end, communicator);
//            //m.local_rehash();

  //df.close();
  }

}


void buildCountIndex(const std::string & filename, MPI_Comm comm) {
  // ==  open file
//            distributed_file df;
//            range = df.open(filename, communicator);  // memmap internally
//            data = df.data();

  using Alphabet = bliss::common::DNA;
  using FileLoaderType = bliss::io::FASTQLoader<CharType, true, false>; // raw data type :  use CharType
  using KmerType = bliss::common::Kmer<21, Alphabet, WordType>;
  using IdType = uint32_t;
  using TupleType = std::pair<KmerType, IdType>;
  using MapType = std::unordered_map<KmerType, IdType, bliss::hash::farm::KmerSuffixHash<KmerType> >;

  //====  now process the file, one L1 block (block partition by MPI Rank) at a time
  // from FileLoader type, get the block iter type and range type
  using FileBlockIterType = typename FileLoaderType::L1BlockType::iterator;

  using ParserType = bliss::io::FASTQParser<FileBlockIterType, void>;
  using SeqType = typename ParserType::SequenceType;
  using SeqIterType = bliss::io::SequencesIterator<ParserType>;


  /// converter from ascii to alphabet values
  using BaseCharIterator = bliss::iterator::transform_iterator<typename SeqType::IteratorType, bliss::common::ASCII2<Alphabet> >;

  /// kmer generation iterator
  using KmerIterType = bliss::common::KmerGenerationIterator<BaseCharIterator, KmerType>;

  /// combine kmer iterator and position iterator to create an index iterator type.
  using KmerIndexIterType = KmerIterType;

  int commSize, commRank;
  MPI_Comm_size(comm, &commSize);
  MPI_Comm_size(comm, &commRank);

  {
    std::chrono::high_resolution_clock::time_point t1, t2;
    std::chrono::duration<double> time_span;

    t1 = std::chrono::high_resolution_clock::now();
    //==== create file Loader
    FileLoaderType loader(comm, filename, 1, sysconf(_SC_PAGE_SIZE));  // this handle is alive through the entire building process.
    typename FileLoaderType::L1BlockType partition = loader.getNextL1Block();

    t2 = std::chrono::high_resolution_clock::now();
     time_span =
         std::chrono::duration_cast<std::chrono::duration<double>>(
             t2 - t1);
     INFOF("R %d file open time: %f", commRank, time_span.count());

     //== reserve
     t1 = std::chrono::high_resolution_clock::now();
     MapType map;
    // modifying the local index directly here causes a thread safety issue, since callback thread is already running.
    // index reserve internally sends a message to itself.
     // call after getting first L1Block to ensure that file is loaded.
     size_t est_size = (loader.getKmerCountEstimate(KmerType::size) + commSize - 1) / commSize;
     map.reserve(est_size);
    t2 = std::chrono::high_resolution_clock::now();
     time_span =
         std::chrono::duration_cast<std::chrono::duration<double>>(
             t2 - t1);
     INFOF("R %d reserve time: %f", commRank, time_span.count());



     // == create kmer iterator
     //            kmer_iter start(data, range);
     //            kmer_iter end(range.second,range.second);

     t1 = std::chrono::high_resolution_clock::now();
     std::vector< KmerType > temp;
     temp.reserve(est_size);

     ParserType parser;
    //=== copy into array
    while (partition.getRange().size() > 0) {
      //== process the chunk of data
      SeqType read;

      //==  and wrap the chunk inside an iterator that emits Reads.
      SeqIterType seqs_start(parser, partition.begin(), partition.end(), partition.getRange().start);
      SeqIterType seqs_end(partition.end());


      //== loop over the reads
      for (; seqs_start != seqs_end; ++seqs_start)
      {
        // first get read
        read = *seqs_start;

        // then compute and store into index (this will generate kmers and insert into index)
        if (read.seqBegin == read.seqEnd) continue;

        //== set up the kmer generating iterators.
        KmerIterType start(BaseCharIterator(read.seqBegin, bliss::common::ASCII2<Alphabet>()), true);
        KmerIterType end(BaseCharIterator(read.seqEnd, bliss::common::ASCII2<Alphabet>()), false);


	temp.insert(temp.end(), start, end);
//        for (auto it = index_start; it != index_end; ++it) {
//          temp.push_back(*it);
//        }
        //std::copy(index_start, index_end, temp.end());
//        INFOF("R %d inserted.  new temp size = %lu", commRank, temp.size());
      }

      partition = loader.getNextL1Block();
    }
    t2 = std::chrono::high_resolution_clock::now();
     time_span =
         std::chrono::duration_cast<std::chrono::duration<double>>(
             t2 - t1);
     INFOF("R %d local kmer array time: %f. temp size = %lu", commRank, time_span.count(), temp.size());


     // distribute
     t1 = std::chrono::high_resolution_clock::now();

     bliss::hash::farm::KmerPrefixHash<KmerType> hash(2 * ceilLog2(commSize));
	if (commSize > 1)
	     mxx::msgs_all2all(temp, [&] ( KmerType const &x) {
       		return (hash(x) % commSize);
     		}, comm);
     t2 = std::chrono::high_resolution_clock::now();
      time_span =
          std::chrono::duration_cast<std::chrono::duration<double>>(
              t2 - t1);
      INFOF("R %d bucket and distribute: %f. temp size = %lu", commRank, time_span.count(), temp.size());

      // == insert into distributed map.
     t1 = std::chrono::high_resolution_clock::now();
	for (auto it = temp.begin(); it != temp.end(); ++it) {
	      map[*it] += 1;
	}
     t2 = std::chrono::high_resolution_clock::now();
      time_span =
          std::chrono::duration_cast<std::chrono::duration<double>>(
              t2 - t1);
      INFOF("R %d inserted: %f. final size = %lu", commRank, time_span.count(), map.size());

	for (auto it = map.begin(); it != map.end(); ++it) {
		assert(it->second > 0 && it->second < est_size);
	}


//            //distributed_map m(element_count);
//            m.reserve(element_count, communicator);
//            m.insert(start, end, communicator);
//            //m.local_rehash();

  //df.close();
  }

}


//void buildCountIndex(const std::string & filename, MPI_Comm comm) {
//  // ==  open file
////            distributed_file df;
////            range = df.open(filename, communicator);  // memmap internally
////            data = df.data();
//
//  using Alphabet = bliss::common::DNA;
//  using FileLoaderType = bliss::io::FASTQLoader<CharType, true, false>; // raw data type :  use CharType
//  using KmerType = bliss::common::Kmer<21, Alphabet, WordType>;
//  using IdType = size_t;
//  using TupleType = std::pair<KmerType, IdType>;
//  using MapType = std::unordered_map<KmerType, IdType, bliss::hash::farm::KmerSuffixHash<KmerType> >;
//
//  //====  now process the file, one L1 block (block partition by MPI Rank) at a time
//  // from FileLoader type, get the block iter type and range type
//  using FileBlockIterType = typename FileLoaderType::L1BlockType::iterator;
//
//  using ParserType = bliss::io::FASTQParser<FileBlockIterType, void>;
//  using SeqType = typename ParserType::SequenceType;
//  using SeqIterType = bliss::io::SequencesIterator<ParserType>;
//
//
//  /// converter from ascii to alphabet values
//  using BaseCharIterator = bliss::iterator::transform_iterator<typename SeqType::IteratorType, bliss::common::ASCII2<Alphabet> >;
//
//  /// kmer generation iterator
//  using KmerIterType = bliss::common::KmerGenerationIterator<BaseCharIterator, KmerType>;
//  /// kmer position iterator type
//  using IdIterType = bliss::iterator::SequenceIdIterator<IdType>;
//
//  /// combine kmer iterator and position iterator to create an index iterator type.
//  using KmerIndexIterType = bliss::iterator::ZipIterator<KmerIterType, IdIterType>;
//
//  int commSize, commRank;
//  MPI_Comm_size(comm, &commSize);
//  MPI_Comm_size(comm, &commRank);
//
//  {
//    std::chrono::high_resolution_clock::time_point t1, t2;
//    std::chrono::duration<double> time_span;
//
//    t1 = std::chrono::high_resolution_clock::now();
//    //==== create file Loader
//    FileLoaderType loader(comm, filename, 1, sysconf(_SC_PAGE_SIZE));  // this handle is alive through the entire building process.
//    typename FileLoaderType::L1BlockType partition = loader.getNextL1Block();
//
//    t2 = std::chrono::high_resolution_clock::now();
//     time_span =
//         std::chrono::duration_cast<std::chrono::duration<double>>(
//             t2 - t1);
//     INFOF("R %d file open time: %f", commRank, time_span.count());
//
//     //== reserve
//     t1 = std::chrono::high_resolution_clock::now();
//     MapType map;
//    // modifying the local index directly here causes a thread safety issue, since callback thread is already running.
//    // index reserve internally sends a message to itself.
//     // call after getting first L1Block to ensure that file is loaded.
//     size_t est_size = (loader.getKmerCountEstimate(KmerType::size) + commSize - 1) / commSize;
//     map.reserve(est_size);
//    t2 = std::chrono::high_resolution_clock::now();
//     time_span =
//         std::chrono::duration_cast<std::chrono::duration<double>>(
//             t2 - t1);
//     INFOF("R %d reserve time: %f", commRank, time_span.count());
//
//
//
//     // == create kmer iterator
//     //            kmer_iter start(data, range);
//     //            kmer_iter end(range.second,range.second);
//
//     t1 = std::chrono::high_resolution_clock::now();
//     std::vector< TupleType > temp;
//     temp.reserve(est_size);
//
//     ParserType parser;
//    //=== copy into array
//    while (partition.getRange().size() > 0) {
//      //== process the chunk of data
//      SeqType read;
//
//      //==  and wrap the chunk inside an iterator that emits Reads.
//      SeqIterType seqs_start(parser, partition.begin(), partition.end(), partition.getRange().start);
//      SeqIterType seqs_end(partition.end());
//
//
//      //== loop over the reads
//      for (; seqs_start != seqs_end; ++seqs_start)
//      {
//        // first get read
//        read = *seqs_start;
//
//        // then compute and store into index (this will generate kmers and insert into index)
//        if (read.seqBegin == read.seqEnd) continue;
//
//        //== set up the kmer generating iterators.
//        KmerIterType start(BaseCharIterator(read.seqBegin, bliss::common::ASCII2<Alphabet>()), true);
//        KmerIterType end(BaseCharIterator(read.seqEnd, bliss::common::ASCII2<Alphabet>()), false);
//
//        //== set up the position iterators
//        IdIterType id_start(read.id);
//        IdIterType id_end(read.id);
//
//        // ==== set up the zip iterators
//        KmerIndexIterType index_start(start, id_start);
//        KmerIndexIterType index_end(end, id_end);
//
//        std::copy(index_start, index_end, temp.end());
//
//      }
//
//      partition = loader.getNextL1Block();
//    }
//    t2 = std::chrono::high_resolution_clock::now();
//     time_span =
//         std::chrono::duration_cast<std::chrono::duration<double>>(
//             t2 - t1);
//     INFOF("R %d local kmer array time: %f", commRank, time_span.count());
//
//
//     // distribute
//     t1 = std::chrono::high_resolution_clock::now();
//
//     bliss::hash::farm::KmerPrefixHash<KmerType> hash(ceilLog2(commSize));
//     mxx::msgs_all2all(temp, [&] ( TupleType const &x) {
//       return (hash(x.first) % commSize);
//     }, comm);
//     t2 = std::chrono::high_resolution_clock::now();
//      time_span =
//          std::chrono::duration_cast<std::chrono::duration<double>>(
//              t2 - t1);
//      INFOF("R %d bucket and distribute: %f", commRank, time_span.count());
//
//      // == insert into distributed map.
//      map.insert(temp.begin(), temp.end());
//
//
////            //distributed_map m(element_count);
////            m.reserve(element_count, communicator);
////            m.insert(start, end, communicator);
////            //m.local_rehash();
//
//  //df.close();
//      INFOF("R %d inserted %lu into map.", commRank, map.size());
//  }
//
//}


void testIndex(MPI_Comm comm, const std::string & filename, std::string testname, void (*func)(const std::string & filename, MPI_Comm comm) ) {

  int nprocs = 1;
  int rank = 0;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &rank);

  DEBUGF("test nthreads is %d", nthreads);



  std::chrono::high_resolution_clock::time_point t1, t2;
  std::chrono::duration<double> time_span;
  std::vector<std::string> timespan_names;
  std::vector<double> timespans;

  size_t entries = 0;

  INFOF("RANK %d: Testing %s", rank, testname.c_str());

  t1 = std::chrono::high_resolution_clock::now();
  // initialize index
  DEBUGF("RANK %d: ***** initializing %s.", rank, testname.c_str());

  double callback_time = 0;


  // start processing.  enclosing with braces to make sure loader is destroyed before MPI finalize.
  DEBUGF("RANK %d: ***** building index first pass.  %d threads, callback_time %f ", rank, nthreads, callback_time);

  func(filename, comm);
  //kmer_index.flush();

//  MPI_Barrier(comm);
  INFO("RANK " << rank << " Index Building 1 for " << filename );

  t2 = std::chrono::high_resolution_clock::now();
  time_span =
      std::chrono::duration_cast<std::chrono::duration<double>>(
          t2 - t1);
  timespan_names.push_back("build");
  timespans.push_back(time_span.count());



  std::stringstream ss;
  std::copy(timespan_names.begin(), timespan_names.end(), std::ostream_iterator<std::string>(ss, ","));
  std::stringstream ss2;
  std::copy(timespans.begin(), timespans.end(), std::ostream_iterator<double>(ss2, ","));


  INFOF("Rank %d Test %s phases [%s] times [%s]", rank, testname.c_str(), ss.str().c_str(), ss2.str().c_str());


}




/**
 *
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char** argv) {

  //////////////// init logging
  LOG_INIT();

  //////////////// parse parameters

  std::string filename("/home/tpan/src/bliss/test/data/test.medium.fastq");
  //std::string filename("/home/tpan/src/bliss/test/data/test.fastq");
  //std::string filename("/mnt/data/1000genome/HG00096/sequence_read/SRR077487_1.filt.fastq");
  if (argc > 1)
  {
    filename.assign(argv[1]);
  }


  int rank = 0;
	int nthreads = 1;
  //////////////// initialize MPI and openMP
#ifdef USE_MPI

  if (nthreads > 1) {

    int provided;

    // one thread will be making all MPI calls.
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);

    if (provided < MPI_THREAD_FUNNELED) {
      ERRORF("The MPI Library Does not have thread support.");
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  } else {
    MPI_Init(&argc, &argv);
  }

  MPI_Comm comm = MPI_COMM_WORLD;

  MPI_Comm_rank(comm, &rank);

  {
    char hostname[256];
    memset(hostname, 0, 256);
    gethostname(hostname, 256);
    INFOF("Rank %d hostname [%s]\n", rank, hostname);
  }
  MPI_Barrier(comm);

  if (rank == 0)
    INFOF("USE_MPI is set");
#else
  static_assert(false, "MPI used although compilation is not set to use MPI");
#endif


  testIndex(comm, filename, "single thread, position index.", buildPosIndex);

  testIndex(comm, filename, "single thread, count index.", buildCountIndex);

  testIndex(comm, filename , "single thread, pos+qual index", buildPosQualIndex);



//
//
//#if defined(KMOLECULEINDEX)
//  {
//    using IndexType = bliss::index::retired::KmerPositionIndexOld<21, bliss::common::DNA, bliss::io::FASTQ>;
//    test<IndexType>(comm, filename, nthreads, chunkSize, testQueryOld<IndexType>, "Kmolecule pos index");
//
//  }
//#elif defined(KMERINDEX)
//  {
//    using IndexType = bliss::index::KmerPositionIndex<21, bliss::common::DNA, bliss::io::FASTQ>;
//    test<IndexType>(comm, filename, nthreads, chunkSize, testQuery<IndexType>, "Kmer pos index");
//
//  }
//
//
//#endif
//  MPI_Barrier(comm);
//
//
//
//#if defined(KMOLECULEINDEX)
//  {
//  using IndexType = bliss::index::retired::KmerCountIndexOld<21, bliss::common::DNA, bliss::io::FASTQ>;
//  test<IndexType>(comm, filename, nthreads, chunkSize, testQueryOld<IndexType>, "Kmolecule count index");
//
//  }
//
//#elif defined(KMERINDEX)
//
//  {
//  using IndexType = bliss::index::KmerCountIndex<21, bliss::common::DNA, bliss::io::FASTQ>;
//  test<IndexType>(comm, filename, nthreads, chunkSize, testQuery<IndexType>, "Kmer count index");
//
//  }
//
//#endif
//  MPI_Barrier(comm);
//
//
//#if defined(KMOLECULEINDEX)
//
//  // with quality score....
//  {
//  using IndexType = bliss::index::retired::KmerPositionAndQualityIndexOld<21, bliss::common::DNA, bliss::io::FASTQ>;
//
//  test<IndexType>(comm, filename, nthreads, chunkSize, testQueryOld<IndexType>, "Kmolecule pos+qual index");
//
//  }
//#elif defined(KMERINDEX)
//
//  {
//  using IndexType = bliss::index::KmerPositionAndQualityIndex<21, bliss::common::DNA, bliss::io::FASTQ>;
//
//  test<IndexType>(comm, filename, nthreads, chunkSize, testQuery<IndexType>, "Kmer pos+qual index");
//
//  }
//
//#endif
  MPI_Barrier(comm);

  //////////////  clean up MPI.
  MPI_Finalize();

  INFOF("M Rank %d called MPI_Finalize", rank);



  return 0;
}
