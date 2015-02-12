/**
 * @file    KmerIndex.hpp
 * @ingroup
 * @author  tpan
 * @brief
 * @details
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef KMERINDEX_HPP_
#define KMERINDEX_HPP_

#if defined(USE_MPI)
#include "mpi.h"
#endif

#if defined(USE_OPENMP)
#include "omp.h"
#endif

#include <unistd.h>     // sysconf
#include <sys/stat.h>   // block size.

#include "common/Kmer.hpp"
#include "common/base_types.hpp"
#include "common/AlphabetTraits.hpp"
#include "io/fastq_loader.hpp"
#include "io/CommunicationLayer.hpp"
#include "index/distributed_map.hpp"
#include "iterators/zip_iterator.hpp"
#include "iterators/transform_iterator.hpp"
#include "io/sequence_id_iterator.hpp"
#include "io/sequence_iterator.hpp"
#include "common/kmer_iterators.hpp"
#include "index/quality_score_iterator.hpp"

#include "retired/kmer_index_generator.hpp"
#include "retired/kmer_index_functors.hpp"
#include "retired/buffered_transform_iterator.hpp"


namespace bliss
{
  namespace index
  {

static size_t dummy = 0;
static size_t count = 0;

    template<typename KmerType>
    void receiveCountAnswer(std::pair<KmerType, bliss::index::count_t>& answer)
    {
      KmerType key;
      bliss::index::count_t val;
      std::tie(key, val) = answer;

      if (key.toString().length() > 1) dummy += val;
      ++count;
            if ((count % 100000) == 0) INFOF("count result:  %s <=> %d", key.toString().c_str(), val);
    };



    template<typename KmerType, typename IdType>
    void receivePositionAnswer(std::pair<KmerType, IdType>& answer)
    {
      KmerType key;
      IdType val;
      std::tie(key, val) = answer;

      if (key.toString().length() > 1) dummy ^= val.file_pos;

      ++count;
            if ((count % 100000) == 0) INFOF("position result:  %s <=> [%d %d %d %d]", key.toString().c_str(), val.file_id, val.seq_id_msb, val.seq_id, val.pos);

    }



    static double score = 0;

        template<typename KmerType, typename InfoType>
        void receivePositionAndQualityAnswer(std::pair<KmerType, InfoType>& answer)
        {
          KmerType key;
          InfoType val;
          std::tie(key, val) = answer;

          if (key.toString().length() > 1) {
            dummy ^= val.first.file_pos;
            score += val.second;
          }
          ++count;
          if ((count % 100000) == 0) INFOF("position + quality result: %s <=> [%d %d %d %d] %f", key.toString().c_str(), val.first.file_id, val.first.seq_id_msb, val.first.seq_id, val.first.pos, val.second);
        }


    /**
     * @class    bliss::index::KmerIndex
     * @brief    base type  for Kmer Index
     * @details  uses simple kmer generator
     *
     */
    template<typename MapType, typename FileFormat = bliss::io::FASTQ, typename Quality = void, bool ThreadLocal = true>
    class KmerIndex {
      public:
        //==================== COMMON TYPES

        /// DEFINE file loader.  this only provides the L1 and L2 blocks, not reads.
        // raw data type :  use CharType
        using FileLoaderType = bliss::io::FASTQLoader<CharType, false, true>;

        // from FileLoader type, get the block iter type and range type
        using FileBlockIterType = typename FileLoaderType::L2BlockType::iterator;
        using RangeType = typename FileLoaderType::RangeType;

        /// DEFINE the iterator parser to get fastq records.  this type determines both the file format and whether quality scoring is used or not.
        using ParserType = bliss::io::FASTQParser<FileBlockIterType, Quality>;

        // once have ParserType, then get the sequence type and sequence id type
        using SeqType = typename ParserType::SequenceType;

        /// DEFINE the transform iterator type for parsing the FASTQ file into sequence records.
        using SeqIterType = bliss::io::SequencesIterator<ParserType>;

        using KmerType = typename MapType::value_type::first_type;
        using ValueType = typename MapType::value_type::second_type;

        //==========below are to be redefined for each index type

        // TODO: to incorporate quality, use another zip iterator.

      protected:
        /// index instance to store data
        MapType index;

        /// MPI communicator used by index and fileloader.
        MPI_Comm comm;

        /// MPI rank within the communicator
        int rank;

        int commSize;

      public:
        /**
         * initializes the index
         * @param comm
         * @param comm_size
         */
        KmerIndex(MPI_Comm _comm, int _comm_size, int num_threads = 1) :
          index(_comm, _comm_size, num_threads, bliss::KmerPrefixHasher<KmerType>(ceilLog2(_comm_size))),
          comm(_comm), commSize(_comm_size) {
          MPI_Comm_rank(comm, &rank);

        };
        virtual ~KmerIndex() {
          index.finish();
        };

        void finalize() {
          index.finish();
        }

        size_t local_size() const {
          return index.local_size();
        }

        MapType& getLocalIndex() {
          return index;
        }

//        void flushInsert() {
//          index.flushInsert();
//        }

        /**
         * build index, default to num of threads = system PAGE SIZE
         * @param filename
         * @param nthreads
         */
        void build(const std::string &filename, const int &nthreads) {
           build(filename, nthreads, sysconf(_SC_PAGE_SIZE));
        }

        /**
         * @brief build index from a file.
         * @param filename
         * @param nthreads
         * @param chunkSize
         */
        void build(const std::string &filename, const int &nthreads, const int chunkSize) {

          //////////////// now partition and open the file.
          {
            //==== create file Loader
            FileLoaderType loader(comm, filename, nthreads, chunkSize);  // this handle is alive through the entire execution.

            // modifying the index here causes a thread safety issue, since callback thread is already running.
            //  index.reserve((loader.getFileRange().size() + commSize - 1) / commSize);

            //====  now process the file .
            typename FileLoaderType::L1BlockType partition = loader.getNextL1Block();
            while (partition.getRange().size() > 0) {
              buildForL1Block(loader, index, rank, nthreads, chunkSize);
              partition = loader.getNextL1Block();
            }
            index.flush();
          }  // scope to ensure file loader is destroyed.
        }

        //============== Query API. all are blocking.  iAddQuery allows overlapped IO

        /**
         * @brief Add a new query asynchronously.  when the buffer is full, then the query is executed.  result is handled via callback
         * @param kmer
         */
        void sendQuery(const KmerType & kmer) {
          index.asyncLookup(kmer);
        }

        /**
         * @brief send all the queries to remote node.
         * @return
         */
        void flushQuery() {
          index.flush();
        }

        /**
         * @brief send a set of queries in a collective/blocking manner.
         * @note  should only be used by a single thread.
         * @param kmers
         * @return
         */
        void query(const std::vector<KmerType>& kmers) {
          // TODO: for now, use the async api
          for (auto kmer : kmers) {
            index.asyncLookup(kmer);
          }
          index.flush();
          // TODO: later, use MPI collective calls.
        }

      protected:

        void buildForL1Block(FileLoaderType &loader, MapType &index, const int& rank,
                        const int& nthreads, const int& chunkSize)
        {
          // local variables for information only.
          //int nReads = 0;    // read count
          //int nChunks = 0;    // chunk count

          // get L2Block from L1Block, spread work to all threads.
          #pragma omp parallel num_threads(nthreads) shared(loader, nthreads, index, rank) OMP_SHARE_DEFAULT
          {
            //== instantiate a local parser in each thread
            ParserType parser;

            int tid = 0;
            #ifdef USE_OPENMP
                tid = omp_get_thread_num();
            #endif

            //== local variables for loop
            typename FileLoaderType::L2BlockType chunk;

            //== process L2Blocks in a loop.  loader.getNextL2Block is atomic.  the work is shared by all threads.
            chunk = loader.getNextL2Block(tid);
            while (chunk.getRange().size() > 0) {
              buildForL2Block(chunk, parser, index);

              //== get read for next loop iteration.
              chunk = loader.getNextL2Block(tid);
//              ++nChunks;
            }
          }  // compute threads parallel

//          DEBUG("buildIndex rank=" << rank << " nReads=" << nReads << " nChunks=" << nChunks << " elapsed time: " << time_span.count() << "s.");

        }

        /**
         * @brief processes a single L2Block to build the index.  results are inserted into index (atomic)
         * @details  parses a L2Block into short sequences (e.g. for FASTQ, reads)
         *           and then call buildForSequence to process each sequence.
         * @param chunk
         * @param parser
         * @param index
         */
        void buildForL2Block(typename FileLoaderType::L2BlockType &chunk, ParserType& parser, MapType& index) {
          if (chunk.getRange().size() == 0) return;

          //== process the chunk of data
          SeqType read;

          //==  and wrap the chunk inside an iterator that emits Reads.
          SeqIterType seqs_start(parser, chunk.begin(), chunk.end(), chunk.getRange().start);
          SeqIterType seqs_end(chunk.end());

          //== loop over the reads
          for (; seqs_start != seqs_end; ++seqs_start)
          {
            // first get read
            read = *seqs_start;

            // then compute and store into index (this will generate kmers and insert into index)
            buildForSequence(read, index);
          }
        }

        /**
         * @brief processes a single FASTQ read or sequence into index, and insert (atomic)
         * @param read
         * @param index
         */
        virtual void buildForSequence(SeqType &read, MapType& index) = 0;

    };


    /**
     * @class    bliss::index::KmerIndex
     * @brief    base type  for Kmer Index
     * @details  uses kmer iterator
     *
     */
    template<unsigned int Kmer_Size, typename Alphabet, typename FileFormat = bliss::io::FASTQ, bool ThreadLocal = true>
    class KmerCountIndex : public KmerIndex<bliss::index::distributed_counting_map<bliss::Kmer<Kmer_Size, Alphabet>,
                                                                                  bliss::io::CommunicationLayer<ThreadLocal>,
                                                                                  bliss::KmerSuffixHasher<bliss::Kmer<Kmer_Size, Alphabet> >,
                                                                                  bliss::KmerPrefixHasher<bliss::Kmer<Kmer_Size, Alphabet> > >,
                                            FileFormat, void, ThreadLocal>
    {
      public:
        /// DEFINE kmer index type, also used for reverse complement
        using KmerType = bliss::Kmer<Kmer_Size, Alphabet>;

        /// define the index storage type.  using an infix of kmer since prefix is used for distributing to processors.
        using MapType = bliss::index::distributed_counting_map<KmerType,
            bliss::io::CommunicationLayer<ThreadLocal>,
            bliss::KmerSuffixHasher<KmerType>,
            bliss::KmerPrefixHasher<KmerType> >;

        using ValueType = typename MapType::value_type;
        using CountType = typename ValueType::second_type;

        using BaseIndexType = KmerIndex<MapType, FileFormat, void, ThreadLocal>;
        /// DEFINE kmer index type, also used for reverse complement

        using RangeType = typename BaseIndexType::RangeType;

        // once have ParserType, then get the sequence type and sequence id type
        using SeqType = typename BaseIndexType::SeqType;

        //==========below are to be redefined for each index type

        ///////////// INDEX TYPES
        /// converter from ascii to alphabet values
        using BaseCharIterator = bliss::iterator::transform_iterator<typename SeqType::IteratorType, bliss::ASCII2<Alphabet> >;

        /// kmer generation iterator
        typedef bliss::KmerGenerationIterator<BaseCharIterator, KmerType>             KmerIterType;

        /// combine kmer iterator and position iterator to create an index iterator type.
        using KmerIndexIterType = KmerIterType;

      public:
        /**
         * initializes the index
         * @param comm        /// DEFINE kmer index type, also used for reverse complement
         *
         * @param comm_size
         */
        KmerCountIndex(MPI_Comm _comm, int _comm_size, int num_threads = 1) : BaseIndexType(_comm, _comm_size, num_threads)
        {
          this->index.setLookupAnswerCallback(std::function<void(std::pair<KmerType, bliss::index::count_t>&)>(&receiveCountAnswer<KmerType>));
          this->index.init();
        };
        virtual ~KmerCountIndex() {};


      protected:


        /**
         * @brief processes a single FASTQ read or sequence into index, and insert (atomic)
         * @param read
         * @param index
         */
        virtual void buildForSequence(SeqType &read, MapType& index) {

          if (read.seqBegin == read.seqEnd) return;

           //== transform ascii to coded value
          BaseCharIterator charStart(read.seqBegin, bliss::ASCII2<Alphabet>());
          BaseCharIterator charEnd(read.seqEnd, bliss::ASCII2<Alphabet>());

          //== set up the kmer generating iterators.
          KmerIterType start(charStart, true);
          KmerIterType end(charEnd, false);

          // NOTE: need to call *start to actually evaluate.  question is whether ++ should be doing computation.
          for (; start != end; ++start)
          {
            index.insert(*start);  // right side either RVO assignment (if iterator/functor created the object) or copy (if returning a cached object)
                                  // then copy assign
          }
        }



    };



    /**
     * @class    bliss::index::KmerIndex
     * @brief    base type  for Kmer Index
     * @details
     *
     */
    template<unsigned int Kmer_Size, typename Alphabet, typename FileFormat = bliss::io::FASTQ, bool ThreadLocal = true>
    class KmerPositionIndex : public KmerIndex<bliss::index::distributed_multimap<bliss::Kmer<Kmer_Size, Alphabet>,
                                                                                  bliss::io::FASTQSequenceId,
                                                                                        bliss::io::CommunicationLayer<ThreadLocal>,
                                                                                        bliss::KmerSuffixHasher<bliss::Kmer<Kmer_Size, Alphabet> >,
                                                                                        bliss::KmerPrefixHasher<bliss::Kmer<Kmer_Size, Alphabet> > >,
                                                FileFormat, void, ThreadLocal>
    {
      public:
        /// DEFINE kmer index type, also used for reverse complement
        using KmerType = bliss::Kmer<Kmer_Size, Alphabet>;

        /// define the index storage type.  using an infix of kmer since prefix is used for distributing to processors.
        using MapType = bliss::index::distributed_multimap<KmerType,
            bliss::io::FASTQSequenceId,
            bliss::io::CommunicationLayer<ThreadLocal>,
            bliss::KmerSuffixHasher<KmerType>,
            bliss::KmerPrefixHasher<KmerType> >;

        using ValueType = typename MapType::value_type;
        using IdType = typename ValueType::second_type;

        using BaseIndexType = KmerIndex<MapType, FileFormat, void, ThreadLocal>;
        /// DEFINE kmer index type, also used for reverse complement

        using RangeType = typename BaseIndexType::RangeType;

        // once have ParserType, then get the sequence type and sequence id type
        using SeqType = typename BaseIndexType::SeqType;

        //==========below are to be redefined for each index type

        ///////////// INDEX TYPES
        /// converter from ascii to alphabet values
        using BaseCharIterator = bliss::iterator::transform_iterator<typename SeqType::IteratorType, bliss::ASCII2<Alphabet> >;

        /// kmer generation iterator
        typedef bliss::KmerGenerationIterator<BaseCharIterator, KmerType>             KmerIterType;

        /// kmer position iterator type
        using IdIterType = bliss::iterator::SequenceIdIterator<IdType>;

        /// combine kmer iterator and position iterator to create an index iterator type.
        using KmerIndexIterType = bliss::iterator::ZipIterator<KmerIterType, IdIterType>;


      public:
        /**
         * initializes the index
         * @param comm
         * @param comm_size
         */
        KmerPositionIndex(MPI_Comm _comm, int _comm_size, int num_threads = 1)  : BaseIndexType(_comm, _comm_size, num_threads)
        {
          this->index.setLookupAnswerCallback(std::function<void(std::pair<KmerType, IdType>&)>(&receivePositionAnswer<KmerType, IdType>));
          this->index.init();
        };
        virtual ~KmerPositionIndex() {};


      protected:


        /**
         * @brief processes a single FASTQ read or sequence into index, and insert (atomic)
         * @param read
         * @param index
         */
        virtual void buildForSequence(SeqType &read, MapType& index) {

          if (read.seqBegin == read.seqEnd) return;

          //== set up the kmer generating iterators.
          KmerIterType start(BaseCharIterator(read.seqBegin, bliss::ASCII2<Alphabet>()), true);
          KmerIterType end(BaseCharIterator(read.seqEnd, bliss::ASCII2<Alphabet>()), false);

          //== set up the position iterators
          IdIterType id_start(read.id);
          IdIterType id_end(read.id);

          // ==== set up the zip iterators
          KmerIndexIterType index_start(start, id_start);
          KmerIndexIterType index_end(end, id_end);

          for (; index_start != index_end; ++index_start)
          {
            index.insert(*index_start);  // right side either RVO assignment (if iterator/functor created the object) or copy (if returning a cached object)
                                  // then copy assign
          }
        }


    };


    /**
     * @class    bliss::index::KmerIndex
     * @brief    base type  for Kmer Index
     * @details
     *
     */
    template<unsigned int Kmer_Size, typename Alphabet, typename FileFormat = bliss::io::FASTQ, bool ThreadLocal = true>
    class KmerPositionAndQualityIndex : public KmerIndex<bliss::index::distributed_multimap<bliss::Kmer<Kmer_Size, Alphabet>,
                                                                                          std::pair<bliss::io::FASTQSequenceId, float>,
                                                                                          bliss::io::CommunicationLayer<ThreadLocal>,
                                                                                          bliss::KmerSuffixHasher<bliss::Kmer<Kmer_Size, Alphabet> >,
                                                                                          bliss::KmerPrefixHasher<bliss::Kmer<Kmer_Size, Alphabet> > >,
                                                       FileFormat, float, ThreadLocal>
    {
      public:
        /// DEFINE kmer index type, also used for reverse complement
        using KmerType = bliss::Kmer<Kmer_Size, Alphabet>;

        /// define the index storage type.  using an infix of kmer since prefix is used for distributing to processors.
        using MapType = bliss::index::distributed_multimap<KmerType,
            std::pair<bliss::io::FASTQSequenceId, float>,
            bliss::io::CommunicationLayer<ThreadLocal>,
            bliss::KmerSuffixHasher<KmerType>,
            bliss::KmerPrefixHasher<KmerType> >;

        using ValueType = typename MapType::value_type;

        using IdType = typename ValueType::second_type::first_type;
        using QualityType = typename ValueType::second_type::second_type;

        using BaseIndexType = KmerIndex<MapType, FileFormat, float, ThreadLocal>;
        /// DEFINE kmer index type, also used for reverse complement

        using RangeType = typename BaseIndexType::RangeType;

        // once have ParserType, then get the sequence type and sequence id type
        using SeqType = typename BaseIndexType::SeqType;

        //==========below are to be redefined for each index type

        ///////////// INDEX TYPES
        /// converter from ascii to alphabet values
        using BaseCharIterator = bliss::iterator::transform_iterator<typename SeqType::IteratorType, bliss::ASCII2<Alphabet> >;

        /// kmer generation iterator
        typedef bliss::KmerGenerationIterator<BaseCharIterator, KmerType>             KmerIterType;

        /// kmer position iterator type
        using IdIterType = bliss::iterator::SequenceIdIterator<IdType>;

        using QualIterType = bliss::index::QualityScoreGenerationIterator<typename SeqType::IteratorType, Kmer_Size, bliss::index::Illumina18QualityScoreCodec<QualityType> >;

        using KmerInfoIterType = bliss::iterator::ZipIterator<IdIterType, QualIterType>;
        using KmerInfoType = typename ValueType::second_type;

        /// combine kmer iterator and position iterator to create an index iterator type.
        using KmerIndexIterType = bliss::iterator::ZipIterator<KmerIterType, KmerInfoIterType>;


      public:
        /**
         * initializes the index
         * @param comm
         * @param comm_size
         */
        KmerPositionAndQualityIndex(MPI_Comm _comm, int _comm_size, int num_threads = 1) : BaseIndexType(_comm, _comm_size, num_threads)
        {

          this->index.setLookupAnswerCallback(std::function<void(std::pair<KmerType, KmerInfoType>&)>(&receivePositionAndQualityAnswer<KmerType, KmerInfoType>));
          this->index.init();
        };
        virtual ~KmerPositionAndQualityIndex() {};


      protected:



        /**
         * @brief processes a single FASTQ read or sequence into index, and insert (atomic)
         * @param read
         * @param index
         */
        virtual void buildForSequence(SeqType &read, MapType& index) {

          if (read.seqBegin == read.seqEnd || read.qualBegin == read.qualEnd) return;

          //== set up the kmer generating iterators.
          KmerIterType start(BaseCharIterator(read.seqBegin, bliss::ASCII2<Alphabet>()), true);
          KmerIterType end(BaseCharIterator(read.seqEnd, bliss::ASCII2<Alphabet>()), false);

          //== set up the position iterators
          IdIterType id_start(read.id);
          IdIterType id_end(read.id);

          // set up the quality iterator
          QualIterType qual_start(read.qualBegin);
          QualIterType qual_end(read.qualEnd);

          KmerInfoIterType info_start(id_start, qual_start);
          KmerInfoIterType info_end(id_end, qual_end);

          // ==== set up the zip iterators
          KmerIndexIterType index_start(start, info_start);
          KmerIndexIterType index_end(end, info_end);

          for (; index_start != index_end; ++index_start)
          {
            index.insert(*index_start);  // right side either RVO assignment (if iterator/functor created the object) or copy (if returning a cached object)
                                  // then copy assign
          }
        }


    };


    namespace retired {

    /**
     * @class    bliss::index::KmerIndex
     * @brief    base type  for Kmer Index
     * @details  uses simple kmer generator
     *
     */
    template<unsigned int Kmer_Size, typename Alphabet, typename FileFormat = bliss::io::FASTQ, bool ThreadLocal = true>
    class KmerCountIndexOld : public KmerIndex<bliss::index::distributed_counting_map<bliss::Kmer<Kmer_Size, Alphabet>,
                                                                                      bliss::io::CommunicationLayer<ThreadLocal>,
                                                                                      bliss::KmerSuffixHasher<bliss::Kmer<Kmer_Size, Alphabet> >,
                                                                                      bliss::KmerPrefixHasher<bliss::Kmer<Kmer_Size, Alphabet> > >,
                                              FileFormat, void, ThreadLocal>
    {

      public:
        using KmerType = bliss::Kmer<Kmer_Size, Alphabet>;

        using MapType = bliss::index::distributed_counting_map<KmerType,
            bliss::io::CommunicationLayer<ThreadLocal>,
            bliss::KmerSuffixHasher<KmerType>,
            bliss::KmerPrefixHasher<KmerType> >;

        using ValueType = typename MapType::value_type;
        using CountType = typename ValueType::second_type;

        using BaseIndexType = KmerIndex<MapType, FileFormat, void, ThreadLocal>;
        /// DEFINE kmer index type, also used for reverse complement

        using RangeType = typename BaseIndexType::RangeType;

        // once have ParserType, then get the sequence type and sequence id type
        using SeqType = typename BaseIndexType::SeqType;

        //==========below are to be redefined for each index type

        ///////////// INDEX TYPES
        /// define kmer iterator, and the functors for it.
        typedef bliss::index::generate_kmer_simple<SeqType, Alphabet, KmerType>                                  KmoleculeOpType;
        typedef bliss::iterator::buffered_transform_iterator<KmoleculeOpType, typename KmoleculeOpType::BaseIterType> KmoleculeIterType;
        typedef bliss::iterator::transform_iterator<KmoleculeIterType, bliss::utils::KmoleculeToKmerFunctor<KmerType> >       KmerIterType;

        typedef KmerIterType KmerIndexIterType;

      public:
        /**
         * initializes the index
         * @param comm
         * @param comm_size
         */
        KmerCountIndexOld(MPI_Comm _comm, int _comm_size, int num_threads = 1) : BaseIndexType(_comm, _comm_size, num_threads)
        {
          this->index.setLookupAnswerCallback(std::function<void(std::pair<KmerType, bliss::index::count_t>&)>(&receiveCountAnswer<KmerType>));
          this->index.init();
        };
        virtual ~KmerCountIndexOld() {};


      protected:

        /**
         * @brief processes a single FASTQ read or sequence into index, and insert (atomic)
         * @param read
         * @param index
         */
        virtual void buildForSequence(SeqType &read, MapType& index) {

          if (read.seqBegin == read.seqEnd) return;

          //== set up the kmer generating iterators.
          KmoleculeOpType kmer_op;
          bliss::utils::KmoleculeToKmerFunctor<KmerType> transform;
          KmerIndexIterType start(KmoleculeIterType(read.seqBegin, kmer_op), transform);
          KmerIndexIterType end(KmoleculeIterType(read.seqEnd, kmer_op), transform);


          // NOTE: need to call *start to actually evaluate.  question is whether ++ should be doing computation.
          unsigned int i = 0;
          for (; i < (Kmer_Size - 1) && start != end; ++start, ++i);  // compute but discard the first K - 1.

          for (; start != end; ++start, ++i)
          {
            index.insert(*start);  // right side either RVO assignment (if iterator/functor created the object) or copy (if returning a cached object)
                                  // then copy assign
          }
        }

    };




    /**
     * @class    bliss::index::KmerIndex
     * @brief    base type  for Kmer Index
     * @details
     *
     */
    template<unsigned int Kmer_Size, typename Alphabet, typename FileFormat = bliss::io::FASTQ, bool ThreadLocal = true>
    class KmerPositionIndexOld : public KmerIndex<bliss::index::distributed_multimap<bliss::Kmer<Kmer_Size, Alphabet>,
                                                                                     bliss::io::FASTQSequenceId,
                                                                                        bliss::io::CommunicationLayer<ThreadLocal>,
                                                                                        bliss::KmerSuffixHasher<bliss::Kmer<Kmer_Size, Alphabet> >,
                                                                                        bliss::KmerPrefixHasher<bliss::Kmer<Kmer_Size, Alphabet> > >,
                                                  FileFormat, void, ThreadLocal>
        {

        public:
        using KmerType = bliss::Kmer<Kmer_Size, Alphabet>;

        using MapType = bliss::index::distributed_multimap<KmerType,
                                                           bliss::io::FASTQSequenceId,
                                                            bliss::io::CommunicationLayer<ThreadLocal>,
                                                            bliss::KmerSuffixHasher<KmerType>,
                                                            bliss::KmerPrefixHasher<KmerType> >;

        using ValueType = typename MapType::value_type;
        using IdType = typename ValueType::second_type;


        using BaseIndexType = KmerIndex<MapType, FileFormat, void, ThreadLocal>;
        /// DEFINE kmer index type, also used for reverse complement

        using RangeType = typename BaseIndexType::RangeType;

        // once have ParserType, then get the sequence type and sequence id type
        using SeqType = typename BaseIndexType::SeqType;

        //==========below are to be redefined for each index type

        ///////////// INDEX TYPES
        /// define kmer iterator, and the functors for it.
        typedef bliss::index::generate_kmer_simple<SeqType, Alphabet, KmerType>                                  KmoleculeOpType;
        typedef bliss::iterator::buffered_transform_iterator<KmoleculeOpType, typename KmoleculeOpType::BaseIterType> KmoleculeIterType;
        typedef bliss::iterator::transform_iterator<KmoleculeIterType, bliss::utils::KmoleculeToKmerFunctor<KmerType> >       KmerIterType;

         /// kmer position iterator type
        using IdIterType = bliss::iterator::SequenceIdIterator<IdType>;

        /// combine kmer iterator and position iterator to create an index iterator type.
        using KmerIndexIterType = bliss::iterator::ZipIterator<KmerIterType, IdIterType>;

      public:
        /**
         * initializes the index
         * @param comm
         * @param comm_size
         */
        KmerPositionIndexOld(MPI_Comm _comm, int _comm_size, int num_threads = 1) : BaseIndexType(_comm, _comm_size, num_threads)
      {
          this->index.setLookupAnswerCallback(std::function<void(std::pair<KmerType, IdType>&)>(&receivePositionAnswer<KmerType, IdType>));
          this->index.init();
        };
        virtual ~KmerPositionIndexOld() {};


      protected:


        /**
         * @brief processes a single FASTQ read or sequence into index, and insert (atomic)
         * @param read
         * @param index
         */
        virtual void buildForSequence(SeqType &read, MapType& index) {

          if (read.seqBegin == read.seqEnd) return;

          //== set up the kmer generating iterators.
          KmoleculeOpType kmer_op;
          bliss::utils::KmoleculeToKmerFunctor<KmerType> transform;
          KmerIterType start(KmoleculeIterType(read.seqBegin, kmer_op), transform);
          KmerIterType end(KmoleculeIterType(read.seqEnd, kmer_op), transform);

          //== set up the position iterators
          IdIterType id_start(read.id);
          IdIterType id_end(read.id);

          // ==== set up the zip iterators
          KmerIndexIterType index_start(start, id_start);
          KmerIndexIterType index_end(end, id_end);


          unsigned int i = 0;
          for (; i < (Kmer_Size - 1) && index_start != index_end; ++index_start, ++i);  // compute but discard the first K - 1.

          for (; index_start != index_end; ++index_start, ++i)
          {
            index.insert(*index_start);  // right side either RVO assignment (if iterator/functor created the object) or copy (if returning a cached object)
                                  // then copy assign
          }
        }


    };






    /**
     * @class    bliss::index::KmerIndex
     * @brief    base type  for Kmer Index
     * @details
     *
     */
    template<unsigned int Kmer_Size, typename Alphabet, typename FileFormat = bliss::io::FASTQ, bool ThreadLocal = true>
    class KmerPositionAndQualityIndexOld : public KmerIndex<bliss::index::distributed_multimap<bliss::Kmer<Kmer_Size, Alphabet>,
                                                                                               std::pair<bliss::io::FASTQSequenceId, float>,
                                                                                               bliss::io::CommunicationLayer<ThreadLocal>,
                                                                                               bliss::KmerSuffixHasher<bliss::Kmer<Kmer_Size, Alphabet> >,
                                                                                               bliss::KmerPrefixHasher<bliss::Kmer<Kmer_Size, Alphabet> > >,
                                                            FileFormat, float, ThreadLocal>
      {

      public:
        using KmerType = bliss::Kmer<Kmer_Size, Alphabet>;

        using MapType = bliss::index::distributed_multimap<KmerType,
            std::pair<bliss::io::FASTQSequenceId, float>,
        bliss::io::CommunicationLayer<ThreadLocal>,
        bliss::KmerSuffixHasher<KmerType>,
        bliss::KmerPrefixHasher<KmerType> >;

        using ValueType = typename MapType::value_type;


        using IdType = typename ValueType::second_type::first_type;
        using QualityType = typename ValueType::second_type::second_type;

        using BaseIndexType = KmerIndex<MapType, FileFormat, float, ThreadLocal>;
        /// DEFINE kmer index type, also used for reverse complement

        using RangeType = typename BaseIndexType::RangeType;

        // once have ParserType, then get the sequence type and sequence id type
        using SeqType = typename BaseIndexType::SeqType;

        //==========below are to be redefined for each index type

        ///////////// INDEX TYPES
        /// define kmer iterator, and the functors for it.
        typedef bliss::index::generate_kmer_simple<SeqType, Alphabet, KmerType>                                  KmoleculeOpType;
        typedef bliss::iterator::buffered_transform_iterator<KmoleculeOpType, typename KmoleculeOpType::BaseIterType> KmoleculeIterType;
        typedef bliss::iterator::transform_iterator<KmoleculeIterType, bliss::utils::KmoleculeToKmerFunctor<KmerType> >       KmerIterType;


        /// kmer position iterator type
        using IdIterType = bliss::iterator::SequenceIdIterator<IdType>;


        using QualIterType = bliss::index::QualityScoreGenerationIterator<typename SeqType::IteratorType, Kmer_Size, bliss::index::Illumina18QualityScoreCodec<float> >;

        using KmerInfoIterType = bliss::iterator::ZipIterator<IdIterType, QualIterType>;
        using KmerInfoType = typename ValueType::second_type;


        /// combine kmer iterator and position iterator to create an index iterator type.
        using KmerIndexIterType = bliss::iterator::ZipIterator<KmerIterType, KmerInfoIterType>;


      public:
        /**
         * initializes the index
         * @param comm
         * @param comm_size
         */
        KmerPositionAndQualityIndexOld(MPI_Comm _comm, int _comm_size, int num_threads = 1) : BaseIndexType(_comm, _comm_size, num_threads)
        {

          this->index.setLookupAnswerCallback(std::function<void(std::pair<KmerType, KmerInfoType>&)>(&receivePositionAndQualityAnswer<KmerType, KmerInfoType>));
          this->index.init();
        };
        virtual ~KmerPositionAndQualityIndexOld() {};


      protected:


        /**
         * @brief processes a single FASTQ read or sequence into index, and insert (atomic)
         * @param read
         * @param index
         */
        virtual void buildForSequence(SeqType &read, MapType& index) {

          if (read.seqBegin == read.seqEnd || read.qualBegin == read.qualEnd) return;

          //== set up the kmer generating iterators.
          KmoleculeOpType kmer_op;
          bliss::utils::KmoleculeToKmerFunctor<KmerType> transform;
          KmerIterType start(KmoleculeIterType(read.seqBegin, kmer_op), transform);
          KmerIterType end(KmoleculeIterType(read.seqEnd, kmer_op), transform);

          //== set up the position iterators
          IdIterType id_start(read.id);
          IdIterType id_end(read.id);

          // set up the quality iterator
          QualIterType qual_start(read.qualBegin);
          QualIterType qual_end(read.qualEnd);

          KmerInfoIterType info_start(id_start, qual_start);
          KmerInfoIterType info_end(id_end, qual_end);


          // ==== set up the zip iterators
          KmerIndexIterType index_start(start, info_start);
          KmerIndexIterType index_end(end, info_end);

            // NOTE: need to call *start to actually evaluate.  question is whether ++ should be doing computation.
            unsigned int i = 0;
            for (; i < (Kmer_Size - 1) && index_start != index_end; ++index_start, ++i);  // compute but discard the first K - 1.

            for (; index_start != index_end; ++index_start, ++i)
            {
              index.insert(*index_start);  // right side either RVO assignment (if iterator/functor created the object) or copy (if returning a cached object)
                                    // then copy assign
            }
          }


    };






    } /* namespace retired */

  } /* namespace index */
} /* namespace bliss */

#endif /* KMERINDEX_HPP_ */
