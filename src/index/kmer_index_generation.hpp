/**
 * @file		kmer_index_generation.hpp
 * @ingroup
 * @author	tpan
 * @brief
 * @details
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef KMER_INDEX_GENERATION_H_
#define KMER_INDEX_GENERATION_H_

#include <vector>
#include <tuple>

namespace bliss
{
  namespace index
  {

    /**
     * template Param:  T is type of data
     *                  P is number of bins to hash to.
     */
    template<typename T>
    struct XorModulus
    {
        T p;
        T offset;
        XorModulus(T _p, T _offset) : p(_p), offset(_offset) {
          assert(_p > 0);
          assert(_offset >= 0);
          printf("XorModulus nprocs %lu, offset %lu\n", p, offset); fflush(stdout);
        };

        T operator()(const T &v1, const T &v2) {
          return (p == 1 ? offset : (v1 ^ v2) % p + offset);
        }
    };

    // TODO:  need to change the templates.
    // FASTQ_SEQUENCE should be parameterized with base iterator and alphabet - DONE
    // given those + k, should be able to get kmer index element type  - DONE
    // given kmer index element type, and fastq sequence type, should be able to specialize kmer generator type.  DONE
    // given the kmer index element type, should be able to define sendbuffer type.  But won't.  SendBuffer may be of different types.

    // this can become a 1 to n transformer???
    template<typename KmerGenOp, typename SendBuffer, typename HashFunc>
    class KmerIndexGenerator {
      public:
        typedef typename KmerGenOp::SequenceType		SequenceType;
        typedef std::vector<SendBuffer>  				    BufferType;
        typedef bliss::iterator::buffered_transform_iterator<KmerGenOp, typename KmerGenOp::BaseIterType> KmerIter;
        typedef typename KmerGenOp::KmerValueType   KmerType;
        typedef typename KmerGenOp::OutputType      KmerIndexPairType;

        KmerIndexGenerator(int nprocs, int rank, int tid) :
          _nprocs(nprocs), _rank(rank), _tid(tid), h(nprocs, nprocs * tid)  {
          static_assert(std::is_same<typename KmerGenOp::KmerIndexType,
                        typename SendBuffer::ValueType>::value,
                        "Kmer Generation and Send Buffer should use the same type");

          printf("KmerIndexGenerator nprocs: %d rank %d threadid %d\n", nprocs, rank, tid);  fflush(stdout);
        }

        void operator()(SequenceType &read, int j, BufferType &buffers, std::vector<size_t> &counts) {
          KmerGenOp kmer_op(read.id);
          KmerIter start(read.seq, kmer_op);
          KmerIter end(read.seq_end, kmer_op);

          KmerIndexPairType index_kmer;
          //uint64_t kmerCount;


          // NOTE: need to call *start to actually evaluate.  question is whether ++ should be doing computation.
          for (int i = 0; start != end; ++start, ++i)
          {

            if (i < (KmerType::KmerSize::size - 1))
              continue;

            index_kmer = *start;

            // some debugging output
            // printf("kmer send to %lx, key %lx, pos %d, qual %f\n", index_kmer.first, index_kmer.second.kmer, index_kmer.second.id.components.pos, index_kmer.second.qual);

            // sending the kmer.
            //printf("rank %d thread %d, staging to buffer %d\n", rank, tid, index_kmer.first % nprocs + (nprocs * tid) );

            // TODO: abstract this to hide details.
            uint64_t index = this->h(index_kmer.first.kmer, index_kmer.second);
//            printf("rank %d thread %d hashing to %lu of %lu buffers\n", this->_rank, this->_tid, index, buffers.size()); fflush(stdout);
            printf("rank %d thread %d hashing to %lu, fill = %lu\n", this->_rank, this->_tid, index, buffers[index].size());
            buffers[index].buffer(index_kmer.first);
            //      printf("kmer send to %lx, key %lx, pos %d, qual %f\n", index_kmer.first, index_kmer.second.kmer, index_kmer.second.id.components.pos, index_kmer.second.qual);
            //++kmerCount;
            counts[this->_tid] += 1;

          }
        }

      protected:
        int _nprocs;
        int _rank;
        int _tid;
        HashFunc h;

    };

    template<typename KmerGenOp, typename SendBuffer, typename HashFunc, typename QualGenOp>
    class KmerIndexGeneratorWithQuality : public KmerIndexGenerator<KmerGenOp, SendBuffer, HashFunc>
    {
      public:
        typedef KmerIndexGenerator<KmerGenOp, SendBuffer, HashFunc> BaseType;
        typedef typename KmerGenOp::SequenceType                    SequenceType;
        typedef std::vector<SendBuffer>                             BufferType;
        typedef typename KmerGenOp::KmerIndexType                   KmerIndexType;
        typedef typename KmerIndexType::KmerType                    KmerType;
        typedef typename KmerIndexType::SizeType                    KmerSizeType;
        typedef typename KmerGenOp::OutputType                      KmerIndexPairType;
        typedef bliss::iterator::buffered_transform_iterator<KmerGenOp, typename KmerGenOp::BaseIterType> KmerIter;
        typedef bliss::iterator::buffered_transform_iterator<QualGenOp,
            typename QualGenOp::BaseIterType>                       QualIter;
        typedef typename QualGenOp::QualityType                     QualityType;

        static constexpr int K = KmerSizeType::size;

        KmerIndexGeneratorWithQuality(int nprocs, int rank, int tid)
          : BaseType(nprocs, rank, tid)
        {
          static_assert(std::is_same<SequenceType, typename QualGenOp::SequenceType>::value, "Kmer Generation and Quality Generation should use the same type");
          printf("KmerIndexGeneratorWithQuality nprocs: %d rank %d threadid %d\n", nprocs, rank, tid); fflush(stdout);
        }


        void operator()(SequenceType &read, int j, BufferType &buffers, std::vector<size_t> &counts) {

          KmerGenOp kmer_op(read.id);
          KmerIter start(read.seq, kmer_op);
          KmerIter end(read.seq_end, kmer_op);

          QualGenOp qual_op;
          QualIter qstart(read.qual, qual_op);
          QualIter qend(read.qual_end, qual_op);

          KmerIndexPairType index_kmer;
          //uint64_t kmerCount;


          // NOTE: need to call *start to actually evaluate.  question is whether ++ should be doing computation.
          for (int i = 0; (start != end) && (qstart != qend); ++start, ++qstart, ++i)
          {

            if (i < (K - 1))
              continue;

            index_kmer = *start;
            index_kmer.first.qual = *qstart;

            // some debugging output
            // printf("kmer send to %lx, key %lx, pos %d, qual %f\n", index_kmer.first, index_kmer.second.kmer, index_kmer.second.id.components.pos, index_kmer.second.qual);

            // sending the kmer.
            //printf("rank %d thread %d, staging to buffer %d\n", rank, tid, index_kmer.first % nprocs + (nprocs * tid) );
            if (fabs(index_kmer.first.qual) > std::numeric_limits<QualityType>::epsilon())
            {
              // sending the kmer.
              //printf("rank %d thread %d, staging to buffer %d\n", rank, tid, index_kmer.first % nprocs + (nprocs * tid) );

              // TODO: abstract this to hide details.
              uint64_t index = this->h(index_kmer.first.kmer, index_kmer.second);
//              printf("rank %d thread %d hashing to %lu of %lu buffers\n", this->_rank, this->_tid, index, buffers.size()); fflush(stdout);
              printf("rank %d thread %d hashing to %lu, fill = %lu\n", this->_rank, this->_tid, index, buffers[index].size());
              buffers[index].buffer(index_kmer.first);
              //      printf("kmer send to %lx, key %lx, pos %d, qual %f\n", index_kmer.first, index_kmer.second.kmer, index_kmer.second.id.components.pos, index_kmer.second.qual);
              // ++kmerCount;
              counts[this->_tid] += 1;
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



#endif /* KMER_INDEX_GENERATION_H_ */
