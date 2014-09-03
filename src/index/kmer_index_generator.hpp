/**
 * @file    kmer_index_generator.hpp
 * @ingroup
 * @author  tpan
 * @brief
 * @details
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef KMER_INDEX_GENERATOR_H_
#define KMER_INDEX_GENERATOR_H_

#include "iterators/buffered_transform_iterator.hpp"


namespace bliss
{
  namespace index
  {

    // TODO:  need to change the templates.
    // Sequence should be parameterized with base iterator and alphabet - DONE
    // given those + k, should be able to get kmer index element type  - DONE
    // given kmer index element type, and fastq sequence type, should be able to specialize kmer generator type.  DONE
    // given the kmer index element type, should be able to define IndexType type.  But won't.  IndexType may be of different types.
    // alignas(64)
//    struct countType
//    {
//        uint64_t c;
//    };

    // this can become a 1 to n transformer???
    template<typename KmerGenOp, typename IndexType>
    class KmerIndexGenerator {
      public:
        typedef typename KmerGenOp::SequenceType                    SequenceType;
        typedef typename KmerGenOp::KmerIndexType                   KmerIndexType;
        typedef typename KmerIndexType::KmerType                    KmerType;
        typedef typename KmerGenOp::OutputType                      KmerIndexPairType;
        typedef bliss::iterator::buffered_transform_iterator<KmerGenOp, typename KmerGenOp::BaseIterType> KmerIter;

        static constexpr unsigned int K = KmerType::size;

        KmerIndexGenerator() {}

        void operator()(SequenceType &read, IndexType &index) {
          KmerGenOp kmer_op(read.id);
          KmerIter start(read.seq, kmer_op);
          KmerIter end(read.seq_end, kmer_op);

          KmerIndexPairType kmer;
//          int kmerCount = 0;

//          int tid = omp_get_thread_num();

          // NOTE: need to call *start to actually evaluate.  question is whether ++ should be doing computation.
          int i = 0;
//          KmerType idx;
          for (; i < (K - 1) && start != end; ++start, ++i);  // compute but discard the first K - 1.

          for (; start != end; ++start, ++i)
          {
            //kmer = *start;
            index.insert(*start);  // right side either RVO assignment (if iterator/functor created the object) or copy (if returning a cached object)
                                  // then copy assign

            // some debugging output
            // printf("kmer send to %lx, key %lx, pos %d, qual %f\n", kmer.first, kmer.second.kmer, kmer.second.id.pos, kmer.second.qual);

            // sending the kmer.
            //printf("rank %d thread %d, staging to buffer %d\n", rank, tid, kmer.first % nprocs + (nprocs * tid) );

            // TODO: abstract this to hide details.
//            idx = this->_hash(kmer.first.kmer, kmer.second);
//            if (counts[this->_tid] % 100000 == 0)
//              printf("rank %d thread %d hashing to %lu of %lu buffers\n", this->_rank, this->_tid, index, buffers.size()); fflush(stdout);
            //printf("rank %d thread %d hashing to %lu, fill = %lu\n", this->_rank, tid, index, buffers[index].size());
//            index.insert(idx, kmer.first);
            //      printf("kmer send to %lx, key %lx, pos %d, qual %f\n", kmer.first, kmer.second.kmer, kmer.second.id.pos, kmer.second.qual);
//            ++kmerCount;


          }
//          counts[tid].c += kmerCount;
        }


    };

    template<typename KmerGenOp, typename IndexType, typename QualGenOp>
    class KmerIndexGeneratorWithQuality : public KmerIndexGenerator<KmerGenOp, IndexType>
    {
      public:
        typedef KmerIndexGenerator<KmerGenOp, IndexType>  BaseType;
        typedef typename KmerGenOp::SequenceType                    SequenceType;
        typedef typename KmerGenOp::KmerIndexType                   KmerIndexType;
        typedef typename KmerIndexType::KmerType                    KmerType;
        typedef typename KmerGenOp::OutputType                      KmerIndexPairType;
        typedef bliss::iterator::buffered_transform_iterator<KmerGenOp, typename KmerGenOp::BaseIterType> KmerIter;
        typedef bliss::iterator::buffered_transform_iterator<QualGenOp,
            typename QualGenOp::BaseIterType>                       QualIter;
        typedef typename QualGenOp::QualityType                     QualityType;

        static constexpr unsigned int K = KmerType::size;

        KmerIndexGeneratorWithQuality() {
          static_assert(std::is_same<SequenceType,
                        typename QualGenOp::SequenceType>::value,
                        "Kmer Generation and Quality Generation should use the same type");
        }


        void operator()(SequenceType &read, IndexType &index) {

          KmerGenOp kmer_op(read.id);
          KmerIter start(read.seqBegin, kmer_op);
          KmerIter end(read.seqEnd, kmer_op);

          QualGenOp qual_op;
          QualIter qstart(read.qualBegin, qual_op);
          QualIter qend(read.qualEnd, qual_op);

          KmerIndexPairType kmer;
//          int kmerCount = 0;

          // NOTE: need to get tid here. depending on how this is called, thread id may not be initialized correctly.
          //int tid = omp_get_thread_num();

          // NOTE: need to call *start to actually evaluate.  question is whether ++ should be doing computation.
//          KmerType idx;
          int i = 0;
          for (; i < (K - 1) && (start != end) && (qstart != qend); ++start, ++qstart, ++i);   // compute but discard first K-1 chars.

          for (; (start != end) && (qstart != qend); ++start, ++qstart, ++i)
          {

            kmer = *start;
            kmer.second.qual = *qstart;

            // some debugging output
            // printf("kmer send to %lx, key %lx, pos %d, qual %f\n", kmer.first, kmer.second.kmer, kmer.second.id.pos, kmer.second.qual);

            // sending the kmer.
            //printf("rank %d thread %d, staging to buffer %d\n", rank, tid, kmer.first % nprocs + (nprocs * tid) );
            if (fabs(kmer.second.qual) > std::numeric_limits<QualityType>::epsilon())
            {
              // sending the kmer.
              //printf("rank %d thread %d, staging to buffer %d\n", rank, tid, kmer.first % nprocs + (nprocs * tid) );

              index.insert(kmer);

              // TODO: abstract this to hide details.
          // idx = this->_hash(kmer.first.kmer, kmer.second);
//              if (counts[tid] % 100000 == 0)
//                printf("rank %d thread %d hashing to %lu of %lu kmer\n", this->_rank, tid, index, buffers.size()); fflush(stdout);
              //printf("rank %d thread %d hashing %lu to %lu, fill = %lu\n", this->_rank, tid, kmer.first.kmer, index, buffers[index].size());
          // index.insert(idx, kmer.first);
              //      printf("kmer send to %lx, key %lx, pos %d, qual %f\n", kmer.first, kmer.second.kmer, kmer.second.id.pos, kmer.second.qual);
//              ++kmerCount;
            }
            else
            {
              //      printf("BAD kmer quality.  key %lx, pos %d, qual %f\n", kmer.second.kmer, kmer.second.id.pos, kmer.second.qual);
            }

          }
//          counts[tid].c += kmerCount;
        }

    };


  } // namespace index

} // namespace bliss



#endif /* KMER_INDEX_GENERATOR_H_ */
