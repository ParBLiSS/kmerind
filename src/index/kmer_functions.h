/**
 * @file		kmer_functions.h
 * @ingroup
 * @author	tpan
 * @brief
 * @details
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef KMER_FUNCTIONS_H_
#define KMER_FUNCTIONS_H_

#include <vector>
#include <tuple>

namespace bliss
{
  namespace index
  {
    // TODO:  need to change the templates.
    // FASTQ_SEQUENCE should be parameterized with base iterator and alphabet - DONE
    // given those + k, should be able to get kmer index element type  - DONE
    // given kmer index element type, and fastq sequence type, should be able to specialize kmer generator type.  DONE
    // given the kmer index element type, should be able to define sendbuffer type.

    // this can become a 1 to n transformer???
    template<typename SendBuffer, typename KmerGenOp>
    class KmerIndexGenerator {
      public:
        typedef typename KmerGenOp::SequenceType		SequenceType;
        typedef std::vector<SendBuffer>  				    BufferType;
        typedef bliss::iterator::buffered_transform_iterator<KmerGenOp, typename KmerGenOp::BaseIterType> KmerIter;
        typedef typename KmerGenOp::KmerType 			  KmerType;

        KmerIndexGenerator(int nprocs, int rank, int tid) :
          _nprocs(nprocs), _rank(rank), _tid(tid)  {
          static_assert(std::is_same<SequenceType,
                        typename SendBuffer::ValueType>::value,
                        "Kmer Generation and Send Buffer should use the same type");
        }

      protected:
        int _nprocs;
        int _rank;
        int _tid;
        constexpr int K = KmerType::KmerSize::K;

        void impl(SeqType &read, int j, BufType &buffers, std::vector<size_t> &counts) {
          KmerGenOp kmer_op(read.id);
          KmerIter start(read.seq, kmer_op);
          KmerIter end(read.seq_end, kmer_op);

          std::pair<typename KmerType::kmer_type, KmerType> index_kmer;
          //uint64_t kmerCount;


          // NOTE: need to call *start to actually evaluate.  question is whether ++ should be doing computation.
          for (int i = 0; start != end; ++start, ++i)
          {

            if (i < (K - 1))
              continue;

            index_kmer = *start;

            // some debugging output
            // printf("kmer send to %lx, key %lx, pos %d, qual %f\n", index_kmer.first, index_kmer.second.kmer, index_kmer.second.id.components.pos, index_kmer.second.qual);

            // sending the kmer.
            //printf("rank %d thread %d, staging to buffer %d\n", rank, tid, index_kmer.first % nprocs + (nprocs * tid) );

            // TODO: abstract this to hide details.
            buffers[index_kmer.first % nprocs + (nprocs * tid) ].buffer(index_kmer.second);
            //      printf("kmer send to %lx, key %lx, pos %d, qual %f\n", index_kmer.first, index_kmer.second.kmer, index_kmer.second.id.components.pos, index_kmer.second.qual);
            //++kmerCount;
            counts[tid] += 1;
          }
        }
    };

    template<typename SendBuffer, typename KmerGenOp, typename QualGenOp>
    class KmerIndexGeneratorQuality : public KmerIndexGenerator<SendBuffer,
                                          KmerGenOp>
    {
      public:
        typedef KmerIndexGeneratorNoQuality<SendBuffer, KmerGenOp> BaseType;
        typedef typename BaseType::SequenceType SequenceType;
        typedef typename BaseType::BufferType BufferType;
        typedef typename BaseType::KmerIter KmerIter;
        typedef typename BaseType::KmerType KmerType;
        typedef bliss::iterator::buffered_transform_iterator<QualGenOp,
            typename QualGenOp::BaseIterType> QualIter;

        KmerIndexGeneratorQuality(int nprocs, int rank, int tid)
            : BaseType(nprocs, rank, tid)
        {
          static_assert(std::is_same<SequenceType, typename QualGenOp::SequenceType>::value, "Kmer Generation and Quality Generation should use the same type");
        }

      protected:

        void impl(SeqType &read, int j, BufType &buffers, std::vector<size_t> &counts) {

          KmerGenOp kmer_op(read.id);
          KmerIter start(read.seq, kmer_op);
          KmerIter end(read.seq_end, kmer_op);

          QualGenOp qual_op;
          QualIter qstart(read.qual, qual_op);
          QualIter qend(read.qual_end, qual_op);

          std::pair<typename KmerType::kmer_type, KmerType> index_kmer;
          //uint64_t kmerCount;


          // NOTE: need to call *start to actually evaluate.  question is whether ++ should be doing computation.
          for (int i = 0; (start != end) && (qstart != qend); ++start, ++qstart, ++i)
          {

            if (i < (K - 1))
              continue;

            index_kmer = *start;
            index_kmer.second.qual = *qstart;

            // some debugging output
            // printf("kmer send to %lx, key %lx, pos %d, qual %f\n", index_kmer.first, index_kmer.second.kmer, index_kmer.second.id.components.pos, index_kmer.second.qual);

            // sending the kmer.
            //printf("rank %d thread %d, staging to buffer %d\n", rank, tid, index_kmer.first % nprocs + (nprocs * tid) );
            if (fabs(index_kmer.second.qual) > std::numeric_limits<
                    typename KmerType::QualityType>::epsilon())
            {
              // sending the kmer.
              //printf("rank %d thread %d, staging to buffer %d\n", rank, tid, index_kmer.first % nprocs + (nprocs * tid) );

              // TODO: abstract this to hide details.
              buffers[index_kmer.first % nprocs + (nprocs * tid) ].buffer(index_kmer.second);
              //      printf("kmer send to %lx, key %lx, pos %d, qual %f\n", index_kmer.first, index_kmer.second.kmer, index_kmer.second.id.components.pos, index_kmer.second.qual);
              // ++kmerCount;
              counts[tid] += 1;
            }
            else
            {
              //      printf("BAD kmer quality.  key %lx, pos %d, qual %f\n", index_kmer.second.kmer, index_kmer.second.id.components.pos, index_kmer.second.qual);
            }

          }
        }

    };


  } // namespace index

} // namespace bliss



#endif /* KMER_FUNCTIONS_H_ */
