/**
 * @file    kmer_index.hpp
 * @ingroup index
 * @author  Tony Pan <tpan7@gatech.edu>
 * @brief	High level kmer indexing API
 * @details	3 primary Kmer Index classes are currently provided:
 * 			Kmer Count Index,
 * 			Kmer Position Index, and
 * 			Kmer Position + Quality score index.
 *
 *
 *		Currently, KmerIndex classes only supports FASTQ files, but the intent is to
 *		  allow replaceable module for file reading.
 *
 *		Data distribution is via replaceable Hash functions.  currently,
 *		  PrefixHasher is used to distribute to MPI processes based on the first log_2(P) bits of the kmer
 *		  InfixHasher is used to distribute to threads within a process, based on the next log_2(p) bits of the kmer
 *		  SuffixHasher is used for threadlocal unordered_map hashing
 *
 *		All are deterministic to allow simple lookup processes.
 *
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef KMERINDEX2_HPP_
#define KMERINDEX2_HPP_

#if defined(USE_MPI)
#include "mpi.h"
#endif

//#if defined(USE_OPENMP)
//#include "omp.h"
//#endif

#include <unistd.h>     // sysconf
#include <sys/stat.h>   // block size.
#include <tuple>        // tuple and utility functions
#include <utility>      // pair and utility functions.
#include <type_traits>

#include "io/fastq_loader.hpp"
//#include "io/fasta_loader.hpp"
//#include "io/fasta_iterator.hpp"

#include "utils/logging.h"
#include "common/alphabets.hpp"
#include "common/kmer.hpp"
#include "common/base_types.hpp"
#include "common/sequence.hpp"
#include "utils/kmer_utils.hpp"

#include "io/mxx_support.hpp"
#include <containers/distributed_unordered_map.hpp>
#include <containers/distributed_sorted_map.hpp>
// way too slow.  also not updated. #include <containers/distributed_map.hpp>

#include "io/sequence_iterator.hpp"
#include "io/sequence_id_iterator.hpp"
#include "iterators/transform_iterator.hpp"
#include "common/kmer_iterators.hpp"
#include "iterators/zip_iterator.hpp"
#include "iterators/constant_iterator.hpp"
#include "index/quality_score_iterator.hpp"
#include "wip/kmer_hash.hpp"

#include "utils/timer.hpp"

namespace bliss
{
  namespace index
  {
    namespace kmer
    {

      template <typename MapType, typename Parser>
      class Index {
        protected:

          MapType map;

          MPI_Comm comm;
          int commSize;
          int commRank;


        public:
          using KmerType = typename MapType::key_type;
          using IdType   = typename MapType::mapped_type;
          using TupleType =      std::pair<KmerType, IdType>;

          using Alphabet = typename KmerType::KmerAlphabet;

          Index(MPI_Comm _comm, int _comm_size) : map(_comm, _comm_size), comm(_comm) {
            MPI_Comm_size(_comm, &commSize);
            MPI_Comm_rank(_comm, &commRank);
          }

          virtual ~Index() {};

          MapType & get_map() {
            return map;
          }



          template <typename T, typename P = Parser>
          static size_t read_file(const std::string & filename, std::vector<T> & result, MPI_Comm _comm) {

            int p, rank;
            MPI_Comm_size(_comm, &p);
            MPI_Comm_rank(_comm, &rank);


            using FileLoaderType = bliss::io::FASTQLoader<CharType, false, false>; // raw data type :  use CharType

            //====  now process the file, one L1 block (block partition by MPI Rank) at a time

            size_t before = result.size();
            ::fsc::back_emplace_iterator<std::vector< T > > emplace_iter(result);

            TIMER_INIT(file);
            {  // ensure that fileloader is closed at the end.

              TIMER_START(file);
              //==== create file Loader
              FileLoaderType loader(_comm, filename, 1, sysconf(_SC_PAGE_SIZE));  // this handle is alive through the entire building process.
              typename FileLoaderType::L1BlockType partition = loader.getNextL1Block();
              TIMER_END(file, "open", partition.getRange().size());

              P parser;


              //== reserve
              TIMER_START(file);
              // modifying the local index directly here causes a thread safety issue, since callback thread is already running.
              // index reserve internally sends a message to itself.

              // call after getting first L1Block to ensure that file is loaded.
              size_t est_size = (loader.getKmerCountEstimate(KmerType::size) + p - 1) / p;
              result.reserve(est_size);
              TIMER_END(file, "reserve", est_size);


              TIMER_START(file);
              //=== copy into array
              while (partition.getRange().size() > 0) {

                parser(partition, emplace_iter);

                partition = loader.getNextL1Block();
              }
              TIMER_END(file, "read", result.size());
            }

            TIMER_REPORT_MPI(file, rank, _comm);
            return result.size() - before;
          }



          template <typename T, typename P = Parser>
          static size_t read_file_mpi_subcomm(const std::string & filename, std::vector<T>& result, MPI_Comm _comm) {

            int p, rank;
            MPI_Comm_size(_comm, &p);
            MPI_Comm_rank(_comm, &rank);

            // split the communcator so 1 proc from each host does the read, then redistribute.
            MPI_Comm group_leaders;
            MPI_Comm group;
            ::std::tie(group_leaders, group) = ::mxx2::split_communicator_by_host(_comm);

            int g_size, g_rank;
            MPI_Comm_rank(group, &g_rank);
            MPI_Comm_size(group, &g_size);

            int gl_size = MPI_UNDEFINED;
            int gl_rank = MPI_UNDEFINED;
            if (group_leaders != MPI_COMM_NULL) {
              MPI_Comm_rank(group_leaders, &gl_rank);
              MPI_Comm_size(group_leaders, &gl_size);
            }

            // raw data type :  use CharType.   block partition at L1 and L2.  no buffering at all, since we will be copying data to the group members anyways.
            using FileLoaderType = bliss::io::FASTQLoader<CharType, false, false, bliss::partition::BlockPartitioner<bliss::partition::range<size_t> > >;

            //====  now process the file, one L1 block (block partition by MPI Rank) at a time
            ::fsc::back_emplace_iterator<std::vector< T > > emplace_iter(result);

            size_t before = result.size();

            TIMER_INIT(file);
            {
              typename FileLoaderType::L2BlockType block;
              typename FileLoaderType::RangeType range;
              ::std::vector<CharType> data;

              ::std::vector<typename FileLoaderType::RangeType> ranges;
              ::std::vector<size_t> send_counts;
              typename FileLoaderType::L1BlockType partition;

              // first load the file using the group loader's communicator.
              TIMER_START(file);
              size_t est_size = 1;
              if (group_leaders != MPI_COMM_NULL) {  // ensure file loader is closed properly.
                //==== create file Loader. this handle is alive through the entire building process.
                FileLoaderType loader(group_leaders, filename, g_size);  // for member of group_leaders, each create g_size L2blocks.

                // modifying the local index directly here causes a thread safety issue, since callback thread is already running.

                partition = loader.getNextL1Block();

                // call after getting first L1Block to ensure that file is loaded.
                est_size = (loader.getKmerCountEstimate(KmerType::size) + p - 1) / p;

                //====  now compute the send counts and ranges to be scattered.
                ranges.resize(g_size);
                send_counts.resize(g_size);
                for (size_t i = 0; i < g_size; ++i) {
                  block = loader.getNextL2Block(i);

                  ranges[i] = block.getRange();
                  send_counts[i] = ranges[i].size();
                }

                // scatter the data .  this call here relies on FileLoader still having the memory mapped.
                data = ::mxx2::scatterv(partition.begin(), partition.end(), send_counts, group, 0);

              } else {
                // replicated here because the root's copy needs to be inside the if clause - requires FileLoader to be open for the sender.

                // scatter the data .  this call here relies on FileLoader still having the memory mapped.
                data = ::mxx2::scatterv(partition.begin(), partition.end(), send_counts, group, 0);

              }


              // scatter the ranges
              range = ::mxx2::scatter(ranges, group, 0);

              // now create the L2Blocks from the data  (reuse block)
              block.assign(&(data[0]), &(data[0]) + range.size(), range);

              TIMER_END(file, "open", data.size());



              //== reserve
              TIMER_START(file);
              // broadcast the estimated size
              mxx::datatype<size_t> size_dt;
              MPI_Bcast(&est_size, 1, size_dt.type(), 0, group);
              result.reserve(est_size);
              TIMER_END(file, "reserve", est_size);


              // == parse kmer/tuples iterator
              TIMER_START(file);

              P parser;
              parser(block, emplace_iter);

              TIMER_END(file, "read", result.size());
            }
            TIMER_REPORT_MPI(file, rank, _comm);



            INFOF("freeing group communicator");
            MPI_Comm_free(&group);
            INFOF("freeing group_leader communicator");
            if (group_leaders != MPI_COMM_NULL) MPI_Comm_free(&group_leaders);
            INFOF("DONE WITH communicator release");


            return result.size() - before;
          }


          std::vector<TupleType> find(std::vector<KmerType> &query) {
            return map.find(query);
          }

          std::vector<std::pair<KmerType, size_t> > count(std::vector<KmerType> &query) {
            return map.count(query);
          }

          void erase(std::vector<KmerType> &query) {
            map.erase(query);
          }


          template <typename Predicate>
          std::vector<TupleType> find_if(std::vector<KmerType> &query, Predicate const &pred) {
            return map.find_if(query, pred);
          }
          template <typename Predicate>
          std::vector<TupleType> find_if(Predicate const &pred) {
            return map.find_if(pred);
          }

          template <typename Predicate>
          std::vector<std::pair<KmerType, size_t> > count_if(std::vector<KmerType> &query, Predicate const &pred) {
            return map.count_if(query, pred);
          }

          template <typename Predicate>
          std::vector<std::pair<KmerType, size_t> > count_if(Predicate const &pred) {
            return map.count_if(pred);
          }


          template <typename Predicate>
          void erase_if(std::vector<KmerType> &query, Predicate const &pred) {
            map.erase_if(query, pred);
          }

          template <typename Predicate>
          void erase_if(Predicate const &pred) {
            map.erase_if(pred);
          }

          size_t local_size() {
            return map.local_size();
          }

          template <typename T>
          void insert(std::vector<T> &temp) {
            TIMER_INIT(build);

            TIMER_START(build);
            this->map.reserve(this->map.size() + temp.size());
            TIMER_END(build, "reserve", temp.size());


            // distribute
            TIMER_START(build);
            this->map.insert(temp);
            TIMER_END(build, "insert", this->map.local_size());


            TIMER_START(build);
            size_t m = this->map.update_multiplicity();
            TIMER_END(build, "multiplicity", m);

            TIMER_REPORT_MPI(build, this->commRank, this->comm);

          }

          void build(const std::string & filename, MPI_Comm comm) {
            // ==  open file
            //            distributed_file df;
            //            range = df.open(filename, communicator);  // memmap internally
            //            data = df.data();
            ::std::vector<TupleType> temp;
             this->read_file<TupleType, Parser>(filename, temp, comm);
             this->insert(temp);
          }
          void build_with_mpi_subcomm(const std::string & filename, MPI_Comm comm) {
            // ==  open file
            //            distributed_file df;
            //            range = df.open(filename, communicator);  // memmap internally
            //            data = df.data();

            ::std::vector<TupleType> temp;
            this->read_file_mpi_subcomm<TupleType, Parser>(filename, temp, comm);

            this->insert(temp);
          }

      };

      struct KmerParser {

          template <typename BlockType, typename OutputIt>
          OutputIt operator()(BlockType & partition, OutputIt output_iter) {

            using KmerType = typename ::std::iterator_traits<OutputIt>::value_type;
            using Alphabet = typename KmerType::KmerAlphabet;


              // from FileLoader type, get the block iter type and range type
              using BlockIterType = typename BlockType::iterator;
              using ParserType = bliss::io::FASTQParser<BlockIterType, void>;

              using SeqType = typename ParserType::SequenceType;
              using SeqIterType = bliss::io::SequencesIterator<ParserType>;

              /// converter from ascii to alphabet values
              using BaseCharIterator = bliss::iterator::transform_iterator<typename SeqType::IteratorType, bliss::common::ASCII2<Alphabet> >;

              /// kmer generation iterator
              using KmerIterType = bliss::common::KmerGenerationIterator<BaseCharIterator, KmerType>;

              static_assert(std::is_same<typename std::iterator_traits<KmerIterType>::value_type,
                            KmerType>::value,
                            "input iterator and output iterator's value types differ");


              //== sequence parser type
              ParserType parser;

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

                ::std::copy(start, end, output_iter);
              }

              return output_iter;
            }
      };


      /*
       * TYPE DEFINITIONS
       */
      struct KmerPositionTupleParser {


          template <typename BlockType, typename OutputIt>
          OutputIt operator()(BlockType & partition, OutputIt output_iter) {

            using TupleType = typename ::std::iterator_traits<OutputIt>::value_type;
            static_assert(::std::tuple_size<TupleType>::value == 2, "kmer-pos-qual index data type should be a pair");

            using KmerType = typename std::tuple_element<0, TupleType>::type;
            using Alphabet = typename KmerType::KmerAlphabet;



            //====  now process the file, one L1 block (block partition by MPI Rank) at a time
            // from FileLoader type, get the block iter type and range type
            using FileBlockIterType = typename BlockType::iterator;

            using ParserType = bliss::io::FASTQParser<FileBlockIterType, void>;
            using SeqType = typename ParserType::SequenceType;
            using SeqIterType = bliss::io::SequencesIterator<ParserType>;
            using IdType = typename SeqType::IdType;
            static_assert(::std::is_same<typename std::tuple_element<1, TupleType>::type, IdType>::value, "position type does not match for input and output iterators" );


            /// converter from ascii to alphabet values
            using BaseCharIterator = bliss::iterator::transform_iterator<typename SeqType::IteratorType, bliss::common::ASCII2<Alphabet> >;

            /// kmer generation iterator
            using KmerIterType = bliss::common::KmerGenerationIterator<BaseCharIterator, KmerType>;
            /// kmer position iterator type
            using IdIterType = bliss::iterator::SequenceIdIterator<IdType>;

            /// combine kmer iterator and position iterator to create an index iterator type.
            using KmerIndexIterType = bliss::iterator::ZipIterator<KmerIterType, IdIterType>;

            static_assert(std::is_same<typename std::iterator_traits<KmerIndexIterType>::value_type,
                          TupleType>::value,
                          "input iterator and output iterator's value types differ");


            ParserType parser;

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

              ::std::copy(index_start, index_end, output_iter);
            }
            return output_iter;
          }


      };




      struct KmerPositionQualityTupleParser {


          template <typename BlockType, typename OutputIt>
          OutputIt operator()(BlockType & partition, OutputIt output_iter) {

            using TupleType = typename ::std::iterator_traits<OutputIt>::value_type;
            static_assert(::std::tuple_size<TupleType>::value == 2, "kmer-pos-qual index data type should be a pair");

            using KmerType = typename std::tuple_element<0, TupleType>::type;
            using Alphabet = typename KmerType::KmerAlphabet;

            using KmerInfoType = typename std::tuple_element<1, TupleType>::type;
            static_assert(::std::tuple_size<KmerInfoType>::value == 2, "pos-qual index data type should be a pair");

            using QualType = typename std::tuple_element<1, KmerInfoType>::type;

            using FileBlockIterType = typename BlockType::iterator;

            //====  now process the file, one L1 block (block partition by MPI Rank) at a time

            using ParserType = bliss::io::FASTQParser<FileBlockIterType, QualType>;
            using SeqType = typename ParserType::SequenceType;
            using SeqIterType = bliss::io::SequencesIterator<ParserType>;
            using IdType = typename SeqType::IdType;
            static_assert(::std::is_same<typename std::tuple_element<0, KmerInfoType>::type, IdType>::value, "position type does not match for input and output iterators" );


            /// converter from ascii to alphabet values
            using BaseCharIterator = bliss::iterator::transform_iterator<typename SeqType::IteratorType, bliss::common::ASCII2<Alphabet> >;

            /// kmer generation iterator
            using KmerIterType = bliss::common::KmerGenerationIterator<BaseCharIterator, KmerType>;

            /// kmer position iterator type
            using IdIterType = bliss::iterator::SequenceIdIterator<IdType>;

            using QualIterType = bliss::index::QualityScoreGenerationIterator<typename SeqType::IteratorType, 21, bliss::index::Illumina18QualityScoreCodec<QualType> >;

            /// combine kmer iterator and position iterator to create an index iterator type.
            using KmerInfoIterType = bliss::iterator::ZipIterator<IdIterType, QualIterType>;
            static_assert(std::is_same<typename std::iterator_traits<KmerInfoIterType>::value_type,
                          KmerInfoType>::value,
                          "kmer info input iterator and output iterator's value types differ");


            using KmerIndexIterType = bliss::iterator::ZipIterator<KmerIterType, KmerInfoIterType>;

            static_assert(std::is_same<typename std::iterator_traits<KmerIndexIterType>::value_type,
                          TupleType>::value,
                          "input iterator and output iterator's value types differ");


            ParserType parser;
            //=== copy into array
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


              ::std::copy(index_start, index_end, output_iter);
            }
            return output_iter;
          }
      };





      struct KmerCountTupleParser {


        template <typename BlockType, typename OutputIt>
        OutputIt operator()(BlockType & partition, OutputIt output_iter) {

          using TupleType = typename ::std::iterator_traits<OutputIt>::value_type;
          static_assert(::std::tuple_size<TupleType>::value == 2, "count data type should be a pair");

          using KmerType = typename std::tuple_element<0, TupleType>::type;
          using CountType = typename std::tuple_element<1, TupleType>::type;

          using Alphabet = typename KmerType::KmerAlphabet;

          // from FileLoader type, get the block iter type and range type
          using BlockIterType = typename BlockType::iterator;
          using ParserType = bliss::io::FASTQParser<BlockIterType, void>;

          using SeqType = typename ParserType::SequenceType;
          using SeqIterType = bliss::io::SequencesIterator<ParserType>;

          /// converter from ascii to alphabet values
          using BaseCharIterator = bliss::iterator::transform_iterator<typename SeqType::IteratorType, bliss::common::ASCII2<Alphabet> >;

          /// kmer generation iterator
          using KmerIterType = bliss::common::KmerGenerationIterator<BaseCharIterator, KmerType>;
          using CountIterType = bliss::iterator::ConstantIterator<CountType>;

          using KmerIndexIterType = bliss::iterator::ZipIterator<KmerIterType, CountIterType>;

          static_assert(std::is_same<typename std::iterator_traits<KmerIndexIterType>::value_type, std::pair<KmerType, CountType> >::value,
                        "count zip iterator not producing the right type.");

          static_assert(std::is_same<typename std::iterator_traits<KmerIndexIterType>::value_type,
                        TupleType>::value,
                        "count: input iterator and output iterator's value types differ");


          //== sequence parser type
          ParserType parser;

          //== process the chunk of data
          SeqType read;

          //==  and wrap the chunk inside an iterator that emits Reads.
          SeqIterType seqs_start(parser, partition.begin(), partition.end(), partition.getRange().start);
          SeqIterType seqs_end(partition.end());

          CountIterType count_start(1);


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

            KmerIndexIterType istart(start, count_start);
            KmerIndexIterType iend(end, count_start);

            ::std::copy(istart, iend, output_iter);
          }

          return output_iter;
        }

      };

      template <typename MapType>
      using PositionIndex = Index<MapType, KmerPositionTupleParser >;

      template <typename MapType>
      using PositionQualityIndex = Index<MapType, KmerPositionQualityTupleParser >;

      template <typename MapType>
      using CountIndex = Index<MapType, KmerCountTupleParser >;


    } /* namespace kmer */

  } /* namespace index */
} /* namespace bliss */

#endif /* KMERINDEX2_HPP_ */
