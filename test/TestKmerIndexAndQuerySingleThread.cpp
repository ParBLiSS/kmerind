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
#include <random>
#include <algorithm>
#include "utils/logging.h"

#include "common/alphabets.hpp"


#include "index/kmer_index.hpp"
#include "common/kmer.hpp"
#include "common/base_types.hpp"
#include <string>
#include <sstream>
#include "utils/kmer_utils.hpp"
#include <chrono>

#include "mxx/collective.hpp"
#include "utils/timer.hpp"
#include "io/mxx_support.hpp"



template <typename KmerType>
std::vector<KmerType> readForQuery(const std::string & filename, MPI_Comm comm) {
  using FileLoaderType = bliss::io::FASTQLoader<CharType, true, false>; // raw data type :  use CharType

  //====  now process the file, one L1 block (block partition by MPI Rank) at a time
  // from FileLoader type, get the block iter type and range type
  using FileBlockIterType = typename FileLoaderType::L1BlockType::iterator;

  using ParserType = bliss::io::FASTQParser<FileBlockIterType, void>;
  using SeqType = typename ParserType::SequenceType;
  using SeqIterType = bliss::io::SequencesIterator<ParserType>;

  using Alphabet = typename KmerType::KmerAlphabet;

  /// converter from ascii to alphabet values
  using BaseCharIterator = bliss::iterator::transform_iterator<typename SeqType::IteratorType, bliss::common::ASCII2<Alphabet> >;

  /// kmer generation iterator
  using KmerIterType = bliss::common::KmerGenerationIterator<BaseCharIterator, KmerType>;
  std::vector< KmerType > query;

  int commSize, rank;
  MPI_Comm_size(comm, &commSize);
  MPI_Comm_rank(comm, &rank);


  {
    //==== create file Loader
    FileLoaderType loader(comm, filename, 1, sysconf(_SC_PAGE_SIZE));  // this handle is alive through the entire building process.
    typename FileLoaderType::L1BlockType partition = loader.getNextL1Block();
    size_t est_size = (loader.getKmerCountEstimate(KmerType::size) + commSize - 1) / commSize;

    // == create kmer iterator
    //            kmer_iter start(data, range);
    //            kmer_iter end(range.second,range.second);

    query.reserve(est_size);

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


        query.insert(query.end(), start, end);

      }

      partition = loader.getNextL1Block();
    }
  }
  return query;
}


template<typename KmerType>
void sample(std::vector<KmerType> &query, size_t n, unsigned int seed) {
  std::shuffle(query.begin(), query.begin() + ::std::min(4 * n, query.size()), std::default_random_engine(seed));
  query.erase(query.begin() + n, query.end());
}



template<typename MapType, typename KmerType, typename IdType>
std::vector<std::pair<KmerType, IdType> > testQuery(const MapType & map, std::vector<KmerType> &query, MPI_Comm comm) {
      // ==  open file
    //            distributed_file df;
    //            range = df.open(filename, communicator);  // memmap internally
    //            data = df.data();

  using TupleType = std::pair<KmerType, IdType>;

      int commSize, rank;
      MPI_Comm_size(comm, &commSize);
      MPI_Comm_rank(comm, &rank);

      std::vector<TupleType> results;

      {
        TIMER_INIT(find);

         // the code below is actual query processing code.
        TIMER_START(find);
        std::vector<int> recv_counts;
        bliss::hash::kmer::hash<KmerType, bliss::hash::kmer::detail::farm::hash, bliss::hash::kmer::LexicographicLessCombiner, true> hash(ceilLog2(commSize));
        TIMER_END(find, "begin", query.size());

         // remove duplicates
        TIMER_START(find);
         std::sort(query.begin(), query.end(), [](const KmerType& x, const KmerType& y) {
           return x < y;
         });
         auto new_end = std::unique(query.begin(), query.end(), [] (const KmerType& x, const KmerType& y) {
           return x == y;
         });
         query.erase(new_end, query.end());
         TIMER_END(find, "uniq1", query.size());

         // prepare to distribute
         TIMER_START(find);

         if (commSize > 1)
           recv_counts = mxx2::msgs_all2all(query, [&] ( KmerType const &x) {
              return (hash(x) % commSize);
            }, comm);
         else
           recv_counts.insert(recv_counts.begin(), query.size());
         TIMER_END(find, "a2a1", query.size());



         // == perform the query  - memory utilization is a potential problem.
          TIMER_START(find);
         std::vector<int> send_counts(commSize, 0);
         results.reserve(query.size() * 50);                                      // TODO:  should estimate coverage.
         printf("reserving %lu\n", query.size() * 50);

         TIMER_END(find, "reserve", query.size() * 50);

         TIMER_START(find);
         int k = 0;
         size_t before = 0;
         for (int i = 0; i < commSize; ++i) {
           // work on query from process i.
           //printf("R %d working on query from proce %d\n", commRank, i);

           before = results.size();
           for (int j = 0; j < recv_counts[i]; ++j, ++k) {
              auto range = map.equal_range(query[k]);

              for (auto it2 = range.first; it2 != range.second; ++it2) {
               results.push_back(*it2);
             }
           }
           //if (rank == 0) printf("R %d added %d results for %d queries for process %d\n", rank, send_counts[i], recv_counts[i], i);
           send_counts[i] = results.size() - before;
         }
         TIMER_END(find, "local_find", results.size());

         for (int z = 0; z < send_counts.size(); ++z) {
           printf("%d %d;", z, send_counts[z]);
         }
         printf("\n");

          // ==  send back results.
          TIMER_START(find);

          mxx2::all2all(results, send_counts, comm);
          TIMER_END(find, "a2a2", results.size());

          TIMER_REPORT_MPI(find, rank, comm);

      }

      return results;

    }





template<typename MapType, typename KmerType, typename IdType>
std::vector<std::pair<KmerType, size_t> > testCount(const MapType & map, std::vector<KmerType> &query, MPI_Comm comm) {
      // ==  open file
    //            distributed_file df;
    //            range = df.open(filename, communicator);  // memmap internally
    //            data = df.data();

  using TupleType = std::pair<KmerType, size_t>;

      int commSize, rank;
      MPI_Comm_size(comm, &commSize);
      MPI_Comm_rank(comm, &rank);

      std::vector<TupleType> results;

      {
        TIMER_INIT(count);

        TIMER_START(count);

         // the code below is actual query processing code.
        std::vector<size_t> recv_counts;

        bliss::hash::kmer::hash<KmerType, bliss::hash::kmer::detail::farm::hash, bliss::hash::kmer::LexicographicLessCombiner, true> hash(ceilLog2(commSize));
        TIMER_END(count, "begin", query.size());

         // remove duplicates
        TIMER_START(count);
         std::sort(query.begin(), query.end(), [](const KmerType& x, const KmerType& y) {
           return x < y;
         });
         auto new_end = std::unique(query.begin(), query.end(), [] (const KmerType& x, const KmerType& y) {
           return x == y;
         });
         query.erase(new_end, query.end());
         TIMER_END(count, "uniq1", query.size());

         // prepare to distribute
         TIMER_START(count);

         if (commSize > 1) {
           // distribute (communication part)
           std::vector<size_t> send_counts = mxx2::bucketing(query, [&] ( KmerType const &x) {
             return (hash(x) % commSize);
           }, commSize);

           // distribute (communication part)
           recv_counts = mxx2::all2all(query, send_counts, comm);

         } else
           recv_counts.insert(recv_counts.begin(), query.size());
         TIMER_END(count, "a2a1", query.size());


         // == perform the query  - memory utilization is a potential problem.
         TIMER_START(count);
         std::vector<int> send_counts = recv_counts;
         results.reserve(query.size());  // TODO:  should estimate coverage.

         TIMER_END(count, "reserve", query.size() * 50);

         TIMER_START(count);
         int k = 0;
         for (int i = 0; i < commSize; ++i) {
           // work on query from process i.
           //printf("R %d working on query from proce %d\n", commRank, i);

           for (int j = 0; j < recv_counts[i]; ++j, ++k) {
              auto val = map.count(query[k]);

              results.push_back(std::make_pair(query[k], val));
           }
           //if (rank == 0) printf("R %d added %d results for %d queries for process %d\n", rank, send_counts[i], recv_counts[i], i);

         }
         TIMER_END(count, "local_find", results.size());


          // ==  send back results.
          TIMER_START(count);

          mxx2::all2all(results, send_counts, comm);
          TIMER_END(count, "a2a2", results.size());

          TIMER_REPORT_MPI(count, rank, comm);
      }


      return results;
    }





/*
 * TYPE DEFINITIONS
 */
template <typename MapT>
class PositionIndex {
  public:
    MapT map;

    using MapType = MapT;
    using KmerType = typename MapType::key_type;
    using IdType = typename MapType::mapped_type;
    using TupleType = std::pair<KmerType, IdType>;

    using Alphabet = typename KmerType::KmerAlphabet;


    void build(const std::string & filename, MPI_Comm comm) {
      // ==  open file
    //            distributed_file df;
    //            range = df.open(filename, communicator);  // memmap internally
    //            data = df.data();

      using FileLoaderType = bliss::io::FASTQLoader<CharType, true, false>; // raw data type :  use CharType

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


      int commSize, rank;
      MPI_Comm_size(comm, &commSize);
      MPI_Comm_rank(comm, &rank);

      {
        TIMER_INIT(file);

        TIMER_START(file);
        //==== create file Loader
        FileLoaderType loader(comm, filename, 1, sysconf(_SC_PAGE_SIZE));  // this handle is alive through the entire building process.
        typename FileLoaderType::L1BlockType partition = loader.getNextL1Block();
        TIMER_END(file, "open", partition.getRange().size());

         //== reserve
        TIMER_START(file);
        // modifying the local index directly here causes a thread safety issue, since callback thread is already running.
        // index reserve internally sends a message to itself.
         // call after getting first L1Block to ensure that file is loaded.
         size_t est_size = (loader.getKmerCountEstimate(KmerType::size) + commSize - 1) / commSize;
         std::vector< TupleType > temp;
         temp.reserve(est_size);
         TIMER_END(file, "reserve", est_size);



         // == create kmer iterator
         //            kmer_iter start(data, range);
         //            kmer_iter end(range.second,range.second);

         TIMER_START(file);

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

          }

          partition = loader.getNextL1Block();
        }
        TIMER_END(file, "read", temp.size());

        TIMER_REPORT_MPI(file, rank, comm);


        TIMER_INIT(build);

        TIMER_START(build);
        map.reserve(est_size);
        TIMER_END(build, "reserve", est_size);


         // distribute
         TIMER_START(build);
         bliss::hash::kmer::hash<KmerType, bliss::hash::kmer::detail::farm::hash, bliss::hash::kmer::LexicographicLessCombiner, true> hash(ceilLog2(commSize));
         if (commSize > 1)
           mxx2::msgs_all2all(temp, [&] ( TupleType const &x) {
              return (hash(x.first) % commSize);
            }, comm);
          TIMER_END(build, "a2a", temp.size());

          // == insert into distributed map.
          TIMER_START(build);
          map.insert(temp.begin(), temp.end());
          TIMER_END(build, "insert", map.size());

          TIMER_REPORT_MPI(build, rank, comm);

      }

    }


};





template <typename MapT>
class PositionQualityIndex {
  public:

    MapT map;
    using MapType = MapT;

    using KmerType = typename MapType::key_type;
    using KmerInfoType = typename MapType::mapped_type;
    using TupleType = std::pair<KmerType, KmerInfoType>;

    using Alphabet = typename KmerType::KmerAlphabet;
    using IdType = typename KmerInfoType::first_type;
    using QualType = typename KmerInfoType::second_type;

    void build(const std::string & filename, MPI_Comm comm) {
  // ==  open file
//            distributed_file df;
//            range = df.open(filename, communicator);  // memmap internally
//            data = df.data();


  using FileLoaderType = bliss::io::FASTQLoader<CharType, true, false>; // raw data type :  use CharType
  using FileBlockIterType = typename FileLoaderType::L1BlockType::iterator;





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

  int commSize, rank;
  MPI_Comm_size(comm, &commSize);
  MPI_Comm_rank(comm, &rank);

  {
    TIMER_INIT(file);

    TIMER_START(file);
    //==== create file Loader
    FileLoaderType loader(comm, filename, 1, sysconf(_SC_PAGE_SIZE));  // this handle is alive through the entire building process.
    typename FileLoaderType::L1BlockType partition = loader.getNextL1Block();
    TIMER_END(file, "open", partition.getRange().size());

     //== reserve
     TIMER_START(file);
    // modifying the local index directly here causes a thread safety issue, since callback thread is already running.
    // index reserve internally sends a message to itself.
     // call after getting first L1Block to ensure that file is loaded.
     size_t est_size = (loader.getKmerCountEstimate(KmerType::size) + commSize - 1) / commSize;
     std::vector< TupleType > temp;
     temp.reserve(est_size);
     TIMER_END(file, "reserve", est_size);




     // == create kmer iterator
     //            kmer_iter start(data, range);
     //            kmer_iter end(range.second,range.second);

     TIMER_START(file);

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

      }

      partition = loader.getNextL1Block();
    }
    TIMER_END(file, "read", temp.size());

    TIMER_REPORT_MPI(file, rank, comm);

    TIMER_INIT(build);

    TIMER_START(build);
    map.reserve(est_size);
    TIMER_END(build, "reserve", est_size);


     // distribute
     TIMER_START(build);


     bliss::hash::kmer::hash<KmerType, bliss::hash::kmer::detail::farm::hash, bliss::hash::kmer::LexicographicLessCombiner, true> hash(ceilLog2(commSize));
	if (commSize > 1)
	     mxx2::msgs_all2all(temp, [&] ( TupleType const &x) {
       		return (hash(x.first) % commSize);
     		}, comm);
  TIMER_END(build, "a2a", temp.size());

      // == insert into distributed map.
  TIMER_START(build);
      map.insert(temp.begin(), temp.end());
      TIMER_END(build, "insert", map.size());

      TIMER_REPORT_MPI(build, rank, comm);


  }

}
};


template <typename MapT>
class CountIndex {
  public:
    MapT map;

    using MapType = MapT;

    void build(const std::string & filename, MPI_Comm comm) {
  // ==  open file
//            distributed_file df;
//            range = df.open(filename, communicator);  // memmap internally
//            data = df.data();


  using KmerType = typename MapType::key_type;
//  using IdType = typename MapType::mapped_type;
//  using TupleType = std::pair<KmerType, IdType>;
  using Alphabet = typename KmerType::KmerAlphabet;

  using FileLoaderType = bliss::io::FASTQLoader<CharType, true, false>; // raw data type :  use CharType

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
//  using KmerIndexIterType = KmerIterType;

  int commSize, rank;
  MPI_Comm_size(comm, &commSize);
  MPI_Comm_rank(comm, &rank);

  {
    TIMER_INIT(file);

    TIMER_START(file);
    //==== create file Loader
    FileLoaderType loader(comm, filename, 1, sysconf(_SC_PAGE_SIZE));  // this handle is alive through the entire building process.
    typename FileLoaderType::L1BlockType partition = loader.getNextL1Block();
    TIMER_END(file, "open", partition.getRange().size());

     //== reserve
    TIMER_START(file);
    // modifying the local index directly here causes a thread safety issue, since callback thread is already running.
    // index reserve internally sends a message to itself.
     // call after getting first L1Block to ensure that file is loaded.
     size_t est_size = (loader.getKmerCountEstimate(KmerType::size) + commSize - 1) / commSize;
     std::vector< KmerType > temp;
     temp.reserve(est_size);
     TIMER_END(file, "reserve", est_size);



     // == create kmer iterator
     //            kmer_iter start(data, range);
     //            kmer_iter end(range.second,range.second);

     TIMER_START(file);

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
      }

      partition = loader.getNextL1Block();
    }
    TIMER_END(file, "read", temp.size());

    TIMER_REPORT_MPI(file, rank, comm);

    TIMER_INIT(build);

    TIMER_START(build);
    map.reserve(est_size);
    TIMER_END(build, "reserve", est_size);


     // distribute
     TIMER_START(build);

     bliss::hash::kmer::hash<KmerType, bliss::hash::kmer::detail::farm::hash, bliss::hash::kmer::LexicographicLessCombiner, true> hash(ceilLog2(commSize));
	if (commSize > 1)
	     mxx2::msgs_all2all(temp, [&] ( KmerType const &x) {
       		return (hash(x) % commSize);
     		}, comm);
  TIMER_END(build, "a2a", temp.size());

      // == insert into distributed map.
  TIMER_START(build);
	for (auto it = temp.begin(); it != temp.end(); ++it) {
	      map[*it]++;
	}
  TIMER_END(build, "insert", map.size());

  TIMER_REPORT_MPI(build, rank, comm);



  }

}
};

template <typename IndexType>
void testIndex(MPI_Comm comm, const std::string & filename, std::string test ) {

  int nprocs = 1;
  int rank = 0;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &rank);


  IndexType idx;

  using MapType = typename IndexType::MapType;
  using KmerType = typename MapType::key_type;
  using ValType = typename MapType::mapped_type;

  TIMER_INIT(test);

  if (rank == 0) INFOF("RANK %d / %d: Testing %s", rank, nprocs, test.c_str());

  TIMER_START(test);
  idx.build(filename, comm);
  TIMER_END(test, "build", idx.map.size());


  TIMER_START(test);
  auto query = readForQuery<KmerType>(filename, comm);
  TIMER_END(test, "read query", query.size());

  // for testing, query 1% (else could run out of memory.  if a kmer exists r times, then we may need r^2/p total storage.
  TIMER_START(test);
  unsigned seed = rank * 23;
  sample(query, query.size() / 100, seed);
  TIMER_END(test, "select 1%", query.size());

  auto query_orig = query;

  // process query
  // query
  TIMER_START(test);
  auto results = testQuery<MapType, KmerType, ValType>(idx.map, query, comm);
  TIMER_END(test, "query 1%", results.size());

  query = query_orig;

  // count
  TIMER_START(test);
  auto results2 = testCount<MapType, KmerType, ValType>(idx.map, query, comm);
  TIMER_END(test, "count 1%", results2.size());


  query = query_orig;

  // select 1
  // for testing, query 1% (else could run out of memory.  if a kmer exists r times, then we may need r^2/p total storage.
  TIMER_START(test);
  seed = rank * 57;
  sample(query, 1, seed);
  TIMER_END(test, "select 1", query.size());

  query_orig = query;

  // query 1
  TIMER_START(test);
  auto results3 = testQuery<MapType, KmerType, ValType>(idx.map, query, comm);
  TIMER_END(test, "query 1", results3.size());

  query = query_orig;

  // query 1
  TIMER_START(test);
  auto results4 = testCount<MapType, KmerType, ValType>(idx.map, query, comm);
  TIMER_END(test, "count 1", results4.size());

  TIMER_REPORT_MPI(test, rank, comm);

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

  //std::string filename("/home/tpan/src/bliss/test/data/test.medium.fastq");
  std::string filename("/home/tpan/src/bliss/test/data/test.fastq");

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
    INFOF("Rank %d hostname [%s]", rank, hostname);
  }
  MPI_Barrier(comm);

  if (rank == 0)
    INFOF("USE_MPI is set");
#else
  static_assert(false, "MPI used although compilation is not set to use MPI");
#endif


  using Alphabet = bliss::common::DNA;
  using KmerType = bliss::common::Kmer<21, Alphabet, WordType>;

  using IdType = bliss::io::FASTQ::SequenceId;
  using MapType = std::unordered_multimap<KmerType, IdType, bliss::hash::kmer::hash<KmerType, bliss::hash::kmer::detail::farm::hash, bliss::hash::kmer::LexicographicLessCombiner> >;
  testIndex<PositionIndex<MapType> >(comm, filename, "single thread, position index.");

  MPI_Barrier(comm);


  using MapType2 = std::unordered_map<KmerType, uint32_t, bliss::hash::kmer::hash<KmerType, bliss::hash::kmer::detail::farm::hash, bliss::hash::kmer::LexicographicLessCombiner> >;
  testIndex<CountIndex<MapType2> > (comm, filename, "single thread, count index.");

  MPI_Barrier(comm);

  using QualType = float;
  using KmerInfoType = std::pair<IdType, QualType>;
  using MapType3 = std::unordered_multimap<KmerType, KmerInfoType, bliss::hash::kmer::hash<KmerType, bliss::hash::kmer::detail::farm::hash, bliss::hash::kmer::LexicographicLessCombiner> >;
  testIndex<PositionQualityIndex<MapType3> >(comm, filename , "single thread, pos+qual index");

  MPI_Barrier(comm);

  //////////////  clean up MPI.
  MPI_Finalize();

  INFOF("M Rank %d called MPI_Finalize", rank);



  return 0;
}
