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

namespace mxx2 {
  // ~ 5 to 15% faster compared to standard version, but requires more memory.
  template<typename T, typename _TargetP>
  std::vector<int> msgs_all2all(std::vector<T>& msgs, _TargetP target_p_fun, MPI_Comm comm)
  {
      // get comm parameters
      int p, rank;
      MPI_Comm_size(comm, &p);
      MPI_Comm_rank(comm, &rank);


      // bucket input by their target processor
      // TODO: in-place bucketing??
      std::vector<int> send_counts(p, 0);
      std::vector<int> pids(msgs.size());
      for (int i = 0; i < msgs.size(); ++i)
      {
          pids[i] = target_p_fun(msgs[i]);
          send_counts[pids[i]]++;
      }

      // get all2all params
      std::vector<int> recv_counts = mxx::all2all(send_counts, 1, comm);
      std::vector<int> send_displs = mxx::get_displacements(send_counts);
      std::vector<int> recv_displs = mxx::get_displacements(recv_counts);

      // copy.  need to be able to track current position within each block.
      std::vector<int> offset = send_displs;
      std::vector<T> send_buffer;
      if (msgs.size() > 0)
          send_buffer.resize(msgs.size());
      for (int i = 0; i < msgs.size(); ++i)
      {
          send_buffer[offset[pids[i]]++] = msgs[i];
      }


      // resize messages to fit recv
      std::size_t recv_size = recv_displs[p-1] + recv_counts[p-1];
      msgs.clear();
      //msgs.shrink_to_fit();
      msgs.resize(recv_size);
      //msgs = std::vector<T>(recv_size);

      // get MPI type
      mxx::datatype<T> dt;
      MPI_Datatype mpi_dt = dt.type();

      // all2all
      MPI_Alltoallv(&send_buffer[0], &send_counts[0], &send_displs[0], mpi_dt,
                    &msgs[0], &recv_counts[0], &recv_displs[0], mpi_dt, comm);
      // done, result is returned in vector of input messages

      return recv_counts;
  }



}

template<typename MapType>
    void testQuery(const MapType & map, const std::string & filename, MPI_Comm comm) {
      // ==  open file
    //            distributed_file df;
    //            range = df.open(filename, communicator);  // memmap internally
    //            data = df.data();

  using KmerType = typename MapType::key_type;
  using IdType = typename MapType::mapped_type;
  using TupleType = std::pair<KmerType, IdType>;

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

      int commSize, commRank;
      MPI_Comm_size(comm, &commSize);
      MPI_Comm_rank(comm, &commRank);

      {
        std::chrono::high_resolution_clock::time_point t1, t2;
        std::chrono::duration<double> time_span;

        t1 = std::chrono::high_resolution_clock::now();
        //==== create file Loader
        FileLoaderType loader(comm, filename, 1, sysconf(_SC_PAGE_SIZE));  // this handle is alive through the entire building process.
        typename FileLoaderType::L1BlockType partition = loader.getNextL1Block();
        size_t est_size = (loader.getKmerCountEstimate(KmerType::size) + commSize - 1) / commSize;

        t2 = std::chrono::high_resolution_clock::now();
         time_span =
             std::chrono::duration_cast<std::chrono::duration<double>>(
                 t2 - t1);
         INFOF("R %d file open time: %f", commRank, time_span.count());


         // == create kmer iterator
         //            kmer_iter start(data, range);
         //            kmer_iter end(range.second,range.second);

         t1 = std::chrono::high_resolution_clock::now();
         std::vector< KmerType > query;
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
         INFOF("R %d local kmer array time: %f. query size = %lu", commRank, time_span.count(), query.size());


         // for testing, query 10% (else could run out of memory.  if a kmer exists r times, then we may need r^2/p total storage.
         t1 = std::chrono::high_resolution_clock::now();
         unsigned seed = 1;
         std::shuffle(query.begin(), query.end(), std::default_random_engine(seed));
         query.erase(query.begin() + query.size() / 100, query.end());
         t2 = std::chrono::high_resolution_clock::now();
         time_span =
             std::chrono::duration_cast<std::chrono::duration<double>>(
                 t2 - t1);
         INFOF("R %d randomly select 1 %% of kmers as query input. time %f. new size = %lu", commRank, time_span.count(), query.size());



         // the code below is actual query processing code.


         // remove duplicates
         t1 = std::chrono::high_resolution_clock::now();
         std::sort(query.begin(), query.end(), [](const KmerType& x, const KmerType& y) {
           return x < y;
         });
         auto new_end = std::unique(query.begin(), query.end(), [] (const KmerType& x, const KmerType& y) {
           return x == y;
         });
         query.erase(new_end, query.end());
         t2 = std::chrono::high_resolution_clock::now();
         time_span =
             std::chrono::duration_cast<std::chrono::duration<double>>(
                 t2 - t1);
         INFOF("R %d removed duplicate query: %f. new size = %lu", commRank, time_span.count(), query.size());

         // prepare to distribute
         t1 = std::chrono::high_resolution_clock::now();

         std::vector<int> recv_counts;

         bliss::hash::kmer::hash<KmerType, bliss::hash::kmer::detail::farm::hash, bliss::hash::kmer::LexicographicLessCombiner, true> hash(ceilLog2(commSize));
         if (commSize > 1)
           recv_counts = mxx2::msgs_all2all(query, [&] ( KmerType const &x) {
              return (hash(x) % commSize);
            }, comm);
         else
           recv_counts.insert(recv_counts.begin(), query.size());

         t2 = std::chrono::high_resolution_clock::now();
          time_span =
              std::chrono::duration_cast<std::chrono::duration<double>>(
                  t2 - t1);
          INFOF("R %d distribute query: %f. received query size = %lu", commRank, time_span.count(), query.size());


         // == perform the query  - memory utilization is a potential problem.
         t1 = std::chrono::high_resolution_clock::now();
         std::vector<TupleType> results;
         std::vector<int> send_counts(commSize, 0);
         results.reserve(query.size() * 50);                                      // TODO:  should estimate coverage.
         int k = 0;
         int s = 0;
         for (int i = 0; i < commSize; ++i) {
           // work on query from process i.
           //printf("R %d working on query from proce %d\n", commRank, i);
           send_counts[i] = 0;

           for (int j = 0; j < recv_counts[i]; ++j, ++k) {
              auto range = map.equal_range(query[k]);

              for (auto it = range.first; it != range.second; ++it) {
                if (query[k] != it->first) WARNING(" result does not match query: " << query[k] << ", " << it->first);
              }

             s = std::distance(range.first, range.second);
             if (s > 0) {
               results.insert(results.end(), range.first, range.second);
               send_counts[i] += s;
             }
           }
           if (commRank == 0) printf("R %d added %d results for %d queries for process %d\n", commRank, send_counts[i], recv_counts[i], i);

         }

         t2 = std::chrono::high_resolution_clock::now();
          time_span =
              std::chrono::duration_cast<std::chrono::duration<double> >(
                  t2 - t1);
          INFOF("R %d query result generated: %f. size = %lu", commRank, time_span.count(), results.size());


          // ==  send back results.
          t1 = std::chrono::high_resolution_clock::now();

          mxx::all2all(results, send_counts, comm);

          t2 = std::chrono::high_resolution_clock::now();
           time_span =
               std::chrono::duration_cast<std::chrono::duration<double> >(
                   t2 - t1);
           INFOF("R %d query result sent: %f. final size = %lu", commRank, time_span.count(), results.size());


    //            //distributed_map m(element_count);
    //            m.reserve(element_count, communicator);
    //            m.insert(start, end, communicator);
    //            //m.local_rehash();

      //df.close();


      }

    }


template<typename MapType>
    void test1Query(const MapType & map, const std::string & filename, MPI_Comm comm) {
      // ==  open file
    //            distributed_file df;
    //            range = df.open(filename, communicator);  // memmap internally
    //            data = df.data();

  using KmerType = typename MapType::key_type;
  using IdType = typename MapType::mapped_type;
  using TupleType = std::pair<KmerType, IdType>;

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

      int commSize, commRank;
      MPI_Comm_size(comm, &commSize);
      MPI_Comm_rank(comm, &commRank);

      {
        std::chrono::high_resolution_clock::time_point t1, t2;
        std::chrono::duration<double> time_span;

        t1 = std::chrono::high_resolution_clock::now();
        //==== create file Loader
        FileLoaderType loader(comm, filename, 1, sysconf(_SC_PAGE_SIZE));  // this handle is alive through the entire building process.
        typename FileLoaderType::L1BlockType partition = loader.getNextL1Block();
        size_t est_size = (loader.getKmerCountEstimate(KmerType::size) + commSize - 1) / commSize;

        t2 = std::chrono::high_resolution_clock::now();
         time_span =
             std::chrono::duration_cast<std::chrono::duration<double>>(
                 t2 - t1);
         INFOF("R %d file open time: %f", commRank, time_span.count());


         // == create kmer iterator
         //            kmer_iter start(data, range);
         //            kmer_iter end(range.second,range.second);

         t1 = std::chrono::high_resolution_clock::now();
         std::vector< KmerType > query;
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
         INFOF("R %d local kmer array time: %f. query size = %lu", commRank, time_span.count(), query.size());


         // for testing, query 10% (else could run out of memory.  if a kmer exists r times, then we may need r^2/p total storage.
         t1 = std::chrono::high_resolution_clock::now();
         unsigned seed = commRank * 23;
         std::shuffle(query.begin(), query.end(), std::default_random_engine(seed));
         //if (commRank == 0)
           query.erase(query.begin() + 1, query.end());
         //else query.clear();
         t2 = std::chrono::high_resolution_clock::now();
         time_span =
             std::chrono::duration_cast<std::chrono::duration<double>>(
                 t2 - t1);
         INFOF("R %d randomly select 1 of kmers as query input time %f. new size = %lu", commRank, time_span.count(), query.size());



         // the code below is actual query processing code.


         // remove duplicates
         t1 = std::chrono::high_resolution_clock::now();
         std::sort(query.begin(), query.end(), [](const KmerType& x, const KmerType& y) {
           return x < y;
         });
         auto new_end = std::unique(query.begin(), query.end(), [] (const KmerType& x, const KmerType& y) {
           return x == y;
         });
         query.erase(new_end, query.end());
         t2 = std::chrono::high_resolution_clock::now();
         time_span =
             std::chrono::duration_cast<std::chrono::duration<double>>(
                 t2 - t1);
         INFOF("R %d removed duplicate query: %f. new size = %lu", commRank, time_span.count(), query.size());

         // prepare to distribute
         t1 = std::chrono::high_resolution_clock::now();

         std::vector<int> recv_counts;

         bliss::hash::kmer::hash<KmerType, bliss::hash::kmer::detail::farm::hash, bliss::hash::kmer::LexicographicLessCombiner, true> hash(ceilLog2(commSize));
         if (commSize > 1)
           recv_counts = mxx2::msgs_all2all(query, [&] ( KmerType const &x) {
              return (hash(x) % commSize);
            }, comm);
         else
           recv_counts.insert(recv_counts.begin(), query.size());

         t2 = std::chrono::high_resolution_clock::now();
          time_span =
              std::chrono::duration_cast<std::chrono::duration<double>>(
                  t2 - t1);
          INFOF("R %d distribute query: %f. received query size = %lu", commRank, time_span.count(), query.size());


         // == perform the query  - memory utilization is a potential problem.
         t1 = std::chrono::high_resolution_clock::now();
         std::vector<TupleType> results;
         std::vector<int> send_counts(commSize, 0);
         results.reserve(query.size() * 50);                                      // TODO:  should estimate coverage.
         int k = 0;
         int s = 0;
         for (int i = 0; i < commSize; ++i) {
           // work on query from process i.
           //printf("R %d working on query from proce %d\n", commRank, i);
           send_counts[i] = 0;

           for (int j = 0; j < recv_counts[i]; ++j, ++k) {
              auto range = map.equal_range(query[k]);

              for (auto it = range.first; it != range.second; ++it) {
                if (query[k] != it->first) WARNING(" result does not match query: " << query[k] << ", " << it->first);
              }

             s = std::distance(range.first, range.second);
             if (s > 0) {
               results.insert(results.end(), range.first, range.second);
               send_counts[i] += s;
             }
           }
           if (commRank == 0) printf("R %d added %d results for %d queries for process %d\n", commRank, send_counts[i], recv_counts[i], i);

         }

         t2 = std::chrono::high_resolution_clock::now();
          time_span =
              std::chrono::duration_cast<std::chrono::duration<double> >(
                  t2 - t1);
          INFOF("R %d query result generated: %f. size = %lu", commRank, time_span.count(), results.size());


          // ==  send back results.
          t1 = std::chrono::high_resolution_clock::now();

          mxx::all2all(results, send_counts, comm);

          t2 = std::chrono::high_resolution_clock::now();
           time_span =
               std::chrono::duration_cast<std::chrono::duration<double> >(
                   t2 - t1);
           INFOF("R %d query result sent: %f. final size = %lu", commRank, time_span.count(), results.size());


    //            //distributed_map m(element_count);
    //            m.reserve(element_count, communicator);
    //            m.insert(start, end, communicator);
    //            //m.local_rehash();

      //df.close();
//           if (results.size() > 0) {
//             for (auto res : results) {
//               if (res.first != query[0])
//                 WARNING("R " << commRank << " query " << query[0] << " result " << res.first );
//             }
//           }
      }

    }



template<typename MapType>
    void testCountQuery(const MapType & map, const std::string & filename, MPI_Comm comm) {
      // ==  open file
    //            distributed_file df;
    //            range = df.open(filename, communicator);  // memmap internally
    //            data = df.data();

  using KmerType = typename MapType::key_type;
  using IdType = size_t;
  using TupleType = std::pair<KmerType, IdType>;

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

      int commSize, commRank;
      MPI_Comm_size(comm, &commSize);
      MPI_Comm_rank(comm, &commRank);

      {
        std::chrono::high_resolution_clock::time_point t1, t2;
        std::chrono::duration<double> time_span;

        t1 = std::chrono::high_resolution_clock::now();
        //==== create file Loader
        FileLoaderType loader(comm, filename, 1, sysconf(_SC_PAGE_SIZE));  // this handle is alive through the entire building process.
        typename FileLoaderType::L1BlockType partition = loader.getNextL1Block();
        size_t est_size = (loader.getKmerCountEstimate(KmerType::size) + commSize - 1) / commSize;

        t2 = std::chrono::high_resolution_clock::now();
         time_span =
             std::chrono::duration_cast<std::chrono::duration<double>>(
                 t2 - t1);
         INFOF("R %d file open time: %f", commRank, time_span.count());


         // == create kmer iterator
         //            kmer_iter start(data, range);
         //            kmer_iter end(range.second,range.second);

         t1 = std::chrono::high_resolution_clock::now();
         std::vector< KmerType > query;
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
         INFOF("R %d local kmer array time: %f. query size = %lu", commRank, time_span.count(), query.size());


         // for testing, query 10% (else could run out of memory.  if a kmer exists r times, then we may need r^2/p total storage.
         t1 = std::chrono::high_resolution_clock::now();
         unsigned seed = 1;
         std::shuffle(query.begin(), query.end(), std::default_random_engine(seed));
         query.erase(query.begin() + query.size() / 100, query.end());
         t2 = std::chrono::high_resolution_clock::now();
         time_span =
             std::chrono::duration_cast<std::chrono::duration<double>>(
                 t2 - t1);
         INFOF("R %d randomly select 1 %% of kmers as query input time %f. new size = %lu", commRank, time_span.count(), query.size());



         // the code below is actual query processing code.


         // remove duplicates
         t1 = std::chrono::high_resolution_clock::now();
         std::sort(query.begin(), query.end(), [](const KmerType& x, const KmerType& y) {
           return x < y;
         });
         auto new_end = std::unique(query.begin(), query.end(), [] (const KmerType& x, const KmerType& y) {
           return x == y;
         });
         query.erase(new_end, query.end());
         t2 = std::chrono::high_resolution_clock::now();
         time_span =
             std::chrono::duration_cast<std::chrono::duration<double>>(
                 t2 - t1);
         INFOF("R %d removed duplicate query: %f. new size = %lu", commRank, time_span.count(), query.size());

         // prepare to distribute
         t1 = std::chrono::high_resolution_clock::now();

         std::vector<int> recv_counts;

         bliss::hash::kmer::hash<KmerType, bliss::hash::kmer::detail::farm::hash, bliss::hash::kmer::LexicographicLessCombiner, true> hash(ceilLog2(commSize));
         if (commSize > 1)
           recv_counts = mxx2::msgs_all2all(query, [&] ( KmerType const &x) {
              return (hash(x) % commSize);
            }, comm);
         else
           recv_counts.insert(recv_counts.begin(), query.size());

         t2 = std::chrono::high_resolution_clock::now();
          time_span =
              std::chrono::duration_cast<std::chrono::duration<double>>(
                  t2 - t1);
          INFOF("R %d distribute query: %f. received query size = %lu", commRank, time_span.count(), query.size());


         // == perform the query  - memory utilization is a potential problem.
         t1 = std::chrono::high_resolution_clock::now();
         std::vector<TupleType> results;
         std::vector<int> send_counts = recv_counts;
         results.reserve(query.size());  // TODO:  should estimate coverage.
         int k = 0;
         for (int i = 0; i < commSize; ++i) {
           // work on query from process i.
           //printf("R %d working on query from proce %d\n", commRank, i);

           for (int j = 0; j < recv_counts[i]; ++j, ++k) {
              auto val = map.count(query[k]);

              results.push_back(std::make_pair(query[k], val));
           }
           if (commRank == 0) printf("R %d added %d results for %d queries for process %d\n", commRank, send_counts[i], recv_counts[i], i);

         }

         t2 = std::chrono::high_resolution_clock::now();
          time_span =
              std::chrono::duration_cast<std::chrono::duration<double> >(
                  t2 - t1);
          INFOF("R %d query result generated: %f. size = %lu", commRank, time_span.count(), results.size());


          // ==  send back results.
          t1 = std::chrono::high_resolution_clock::now();

          mxx::all2all(results, send_counts, comm);

          t2 = std::chrono::high_resolution_clock::now();
           time_span =
               std::chrono::duration_cast<std::chrono::duration<double> >(
                   t2 - t1);
           INFOF("R %d query result sent: %f. final size = %lu", commRank, time_span.count(), results.size());


    //            //distributed_map m(element_count);
    //            m.reserve(element_count, communicator);
    //            m.insert(start, end, communicator);
    //            //m.local_rehash();

      //df.close();
      }

    }





/*
 * TYPE DEFINITIONS
 */
template <typename MapType>
class PositionIndex {
  public:
    MapType map;


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


      int commSize, commRank;
      MPI_Comm_size(comm, &commSize);
      MPI_Comm_rank(comm, &commRank);

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
         INFOF("R %d file open time: %f, file range [%lu, %lu)", commRank, time_span.count(), partition.getRange().start, partition.getRange().end);

         //== reserve
         t1 = std::chrono::high_resolution_clock::now();
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

         bliss::hash::kmer::hash<KmerType, bliss::hash::kmer::detail::farm::hash, bliss::hash::kmer::LexicographicLessCombiner, true> hash(ceilLog2(commSize));
      if (commSize > 1)
           mxx2::msgs_all2all(temp, [&] ( TupleType const &x) {
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


};





template <typename MapType>
class PositionQualityIndex {
  public:

    MapType map;
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

  int commSize, commRank;
  MPI_Comm_size(comm, &commSize);
  MPI_Comm_rank(comm, &commRank);

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

     bliss::hash::kmer::hash<KmerType, bliss::hash::kmer::detail::farm::hash, bliss::hash::kmer::LexicographicLessCombiner, true> hash(ceilLog2(commSize));
	if (commSize > 1)
	     mxx2::msgs_all2all(temp, [&] ( TupleType const &x) {
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
};


template <typename MapType>
class CountIndex {
  public:
    MapType map;


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

  int commSize, commRank;
  MPI_Comm_size(comm, &commSize);
  MPI_Comm_rank(comm, &commRank);

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

     bliss::hash::kmer::hash<KmerType, bliss::hash::kmer::detail::farm::hash, bliss::hash::kmer::LexicographicLessCombiner, true> hash(ceilLog2(commSize));
	if (commSize > 1)
	     mxx2::msgs_all2all(temp, [&] ( KmerType const &x) {
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
	      map[*it]++;
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
};

template <typename IndexType>
void testIndex(MPI_Comm comm, const std::string & filename, std::string testname ) {

  int nprocs = 1;
  int rank = 0;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &rank);

  DEBUGF("test nthreads is %d", nthreads);


  IndexType idx;

  std::chrono::high_resolution_clock::time_point t1, t2;
  std::chrono::duration<double> time_span;
  std::vector<std::string> timespan_names;
  std::vector<double> timespans;

  INFOF("RANK %d: Testing %s", rank, testname.c_str());

  t1 = std::chrono::high_resolution_clock::now();
  // initialize index
  DEBUGF("RANK %d: ***** initializing %s.", rank, testname.c_str());

  idx.build(filename, comm);
  //kmer_index.flush();

//  MPI_Barrier(comm);
  INFO("RANK " << rank << " Index Building 1 for " << filename );

  t2 = std::chrono::high_resolution_clock::now();
  time_span =
      std::chrono::duration_cast<std::chrono::duration<double>>(
          t2 - t1);
  timespan_names.push_back("build");
  timespans.push_back(time_span.count());

  t1 = std::chrono::high_resolution_clock::now();

  test1Query(idx.map, filename, comm);
  //kmer_index.flush();

//  MPI_Barrier(comm);
  INFO("RANK " << rank << " Index query1 for " << filename );

  t2 = std::chrono::high_resolution_clock::now();
  time_span =
      std::chrono::duration_cast<std::chrono::duration<double>>(
          t2 - t1);
  timespan_names.push_back("query1");
  timespans.push_back(time_span.count());

  t1 = std::chrono::high_resolution_clock::now();

  testQuery(idx.map, filename, comm);
  //kmer_index.flush();

//  MPI_Barrier(comm);
  INFO("RANK " << rank << " Index query for " << filename );

  t2 = std::chrono::high_resolution_clock::now();
  time_span =
      std::chrono::duration_cast<std::chrono::duration<double>>(
          t2 - t1);
  timespan_names.push_back("query");
  timespans.push_back(time_span.count());




  t1 = std::chrono::high_resolution_clock::now();

  testCountQuery(idx.map, filename, comm);
  //kmer_index.flush();

//  MPI_Barrier(comm);
  INFO("RANK " << rank << " Index count query for " << filename );

  t2 = std::chrono::high_resolution_clock::now();
  time_span =
      std::chrono::duration_cast<std::chrono::duration<double>>(
          t2 - t1);
  timespan_names.push_back("count");
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
