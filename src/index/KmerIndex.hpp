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

namespace bliss
{
  namespace index
  {

    /**
     * @class    bliss::index::KmerIndex
     * @brief
     * @details
     *
     */
    template<unsigned int Kmer_Size, typename Alphabet, typename FileFormat, typename QualityType>
    class KmerIndex;

    template<unsigned int Kmer_Size, typename Alphabet>
    class KmerIndex<Kmer_Size, Alphabet, bliss::io::FASTQ, void>
    {
    };

    template<unsigned int Kmer_Size, typename Alphabet, typename QualityType>
    class KmerIndex<Kmer_Size, Alphabet, bliss::io::FASTQ, QualityType>
    {

      protected:
        static_assert(std::is_floating_point<QualityType>::value, "QualityType should be floating point.");

        ////////// COMMON TYPES

        /// define kmer index type, also used for reverse complement
        using KmerType = bliss::Kmer<Kmer_Size, bliss::AlphabetTraits<Alphabet>::BITS_PER_CHAR>;

        /// define quality type
        typedef bliss::io::FASTQSequenceId                                                IdType;       // hardcode for now

        /// define Range type
        typedef bliss::partition::range<size_t>                                 RangeType;

        /// define the input file type
        // raw data type :  use CharType
        typedef bliss::io::FASTQLoader<CharType, false, true>                   FileLoaderType;
        typedef typename FileLoaderType::L2BlockType::iterator                      FileBlockIterType;

        /// define read type
        typedef bliss::io::SequenceWithQuality<FileBlockIterType, Alphabet, QualityType>  SeqType;

        /// define the transform iterator type
        typedef bliss::io::FASTQParser<FileBlockIterType, Alphabet, QualityType>    ParserType;
        typedef bliss::io::SequencesIterator<ParserType, FileBlockIterType>             SeqIterType;

        /// define kmer quality GENERATOR types
        typedef bliss::index::SangerToLogProbCorrect<QualityType>               QualityEncoderType;
        typedef bliss::index::generate_qual<SeqType, KmerSize, QualityType, QualityEncoderType > QualOpType;


        ///////////// INDEX TYPES
        /// kmer index element type and corresponding generator type
        typedef bliss::index::KmerIndexElementWithIdAndQuality<KmerSize, KmerType, IdType, QualityType>   KmerIndexValueType;
        typedef bliss::index::generate_kmer<SeqType, KmerIndexValueType>                                  KmerOpType;

        /// define the index storage type
        typedef bliss::index::distributed_multimap<KmerType, KmerIndexValueType, bliss::io::CommunicationLayer> IndexType;
        /// define the KmerIndexElement Generator type.
        typedef bliss::index::KmerIndexGeneratorWithQuality<KmerOpType, IndexType, QualOpType>            KmerIndexComputeType;



      public:
        KmerIndex();
        virtual ~KmerIndex();


        template <typename Compute, typename Index>
        void buildIndex(FileLoaderType &loader, Index &index, const int& rank,
                        const int& nthreads, const int& chunkSize)
        {
          std::chrono::high_resolution_clock::time_point t1, t2;
          std::chrono::duration<double> time_span;
          t1 = std::chrono::high_resolution_clock::now();

          INFO("buildIndex rank.tid = " << rank << "." << omp_get_thread_num());

          // local variables for information only.
          int nReads = 0;    // read count
          int nChunks = 0;    // chunk count

          // uses the fastq iterator as the queue itself, instead of master/slave.
          //   at this point, no strong difference.
        #pragma omp parallel num_threads(nthreads) shared(loader, nthreads, index, rank) reduction(+:nReads, nChunks) OMP_SHARE_DEFAULT
          {
            /// instantiate a local parser
            ParserType parser;

            int tid = 0;
        #ifdef USE_OPENMP
            tid = omp_get_thread_num();
        #endif

            /// instantiate a local kmer generator
            Compute op;

            /// local variables for loop
            typename FileLoaderType::L2BlockType chunk;
            SeqType read;

            /// initialize the loop by getting the first chunk
            chunk = loader.getNextL2Block(tid);
            while (chunk.getRange().size() > 0) {
              /// get the chunk of data

              ///  and wrap the chunk inside an iterator that emits Reads.
              SeqIterType fastq_start(parser, chunk.begin(), chunk.end(), chunk.getRange().start);
              SeqIterType fastq_end(chunk.end());

              /// loop over the reads
              for (; fastq_start != fastq_end; ++fastq_start)
              {
                // first get read
                read = *fastq_start;

                // then compute and store into index (this will generate kmers and insert into index)
                op(read, index);
                ++nReads;

                // do a little status print.
                if (nReads % 20000 == 0)
                  INFO("buildIndex rank.tid=" << rank << "." <<  tid << " nReads=" << nReads);
              }

              /// get read for next loop iteration.
              chunk = loader.getNextL2Block(tid);
              ++nChunks;
            }

            INFO("buildIndex rank.tid=" << rank << "." << tid << " nChunks=" << nChunks);

          }  // compute threads parallel

          t2 = std::chrono::high_resolution_clock::now();
          time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
          INFO("buildIndex rank=" << rank << " nReads=" << nReads << " nChunks=" << nChunks << " elapsed time: " << time_span.count() << "s.");

          /// this MPI process is done.  now flush the index to all other nodes.
          t1 = std::chrono::high_resolution_clock::now();
          index.flush();
          t2 = std::chrono::high_resolution_clock::now();
          time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
          INFO("flushIndex rank=" << rank << " elapsed time: " << time_span.count() << "s.");

        }


    };

  } /* namespace index */
} /* namespace bliss */

#endif /* KMERINDEX_HPP_ */
