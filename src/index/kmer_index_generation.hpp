/**
 * @file    kmer_index_generation.hpp
 * @ingroup
 * @author  tpan
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
        const T nprocs;
        const T nthreads;
        XorModulus(T _procs, T _threads) : nprocs(_procs), nthreads(_threads) {
          assert(_procs > 0);
          assert(_threads > 0);
//          printf("XorModulus nprocs %lu, nthreads %lu\n", nprocs, nthreads); fflush(stdout);
        }

        /**
         *
         * @param v1
         * @param v2
         * @param tid     set this at run time so the same object can be used by multiple threads.
         * @return        processor to send this stuff to.
         */
        T operator()(const T &v1, const T &v2) {
          return (nprocs == 1 ? 0 : (v1 ^ v2) % nprocs);
        }
    };

    // TODO:  need to change the templates.
    // FASTQ_SEQUENCE should be parameterized with base iterator and alphabet - DONE
    // given those + k, should be able to get kmer index element type  - DONE
    // given kmer index element type, and fastq sequence type, should be able to specialize kmer generator type.  DONE
    // given the kmer index element type, should be able to define sendbuffer type.  But won't.  SendBuffer may be of different types.
    // alignas(64)
    struct countType
    {
        uint64_t c;
    };

    // this can become a 1 to n transformer???
    template<typename KmerGenOp, typename SendBuffer, typename HashFunc>
    class KmerIndexGenerator {
      public:
        typedef typename KmerGenOp::SequenceType    SequenceType;
        typedef std::vector<SendBuffer>              BufferType;
        typedef typename KmerGenOp::KmerIndexType                   KmerIndexType;
        typedef typename KmerIndexType::KmerType                    KmerType;
        typedef typename KmerIndexType::SizeType                    KmerSizeType;
        typedef typename KmerGenOp::OutputType                      KmerIndexPairType;
        typedef bliss::iterator::buffered_transform_iterator<KmerGenOp, typename KmerGenOp::BaseIterType> KmerIter;

        static constexpr int K = KmerSizeType::size;

        KmerIndexGenerator(int nprocs, int rank, int nthreads) :
          _nprocs(nprocs), _rank(rank), _hash(nprocs, nthreads)  {
          static_assert(std::is_same<typename KmerGenOp::KmerIndexType,
                        typename SendBuffer::ValueType>::value,
                        "Kmer Generation and Send Buffer should use the same type");
        }

        void operator()(SequenceType &read, BufferType &buffers) {
          KmerGenOp kmer_op(read.id);
          KmerIter start(read.seq, kmer_op);
          KmerIter end(read.seq_end, kmer_op);

          KmerIndexPairType index_kmer;
//          int kmerCount = 0;

          int tid = omp_get_thread_num();

          // NOTE: need to call *start to actually evaluate.  question is whether ++ should be doing computation.
          int i = 0;
          for (; i < (K - 1) && start != end; ++start, ++i);
          for (; start != end; ++start, ++i)
          {

            index_kmer = *start;  // right side either RVO assignment (if iterator/functor created the object) or copy (if returning a cached object)
                                  // then copy assign

            // some debugging output
            // printf("kmer send to %lx, key %lx, pos %d, qual %f\n", index_kmer.first, index_kmer.second.kmer, index_kmer.second.id.components.pos, index_kmer.second.qual);

            // sending the kmer.
            //printf("rank %d thread %d, staging to buffer %d\n", rank, tid, index_kmer.first % nprocs + (nprocs * tid) );

            // TODO: abstract this to hide details.
            uint64_t index = this->_hash(index_kmer.first.kmer, index_kmer.second);
//            if (counts[this->_tid] % 100000 == 0)
//              printf("rank %d thread %d hashing to %lu of %lu buffers\n", this->_rank, this->_tid, index, buffers.size()); fflush(stdout);
            //printf("rank %d thread %d hashing to %lu, fill = %lu\n", this->_rank, tid, index, buffers[index].size());
            buffers[index].buffer(index_kmer.first);
            //      printf("kmer send to %lx, key %lx, pos %d, qual %f\n", index_kmer.first, index_kmer.second.kmer, index_kmer.second.id.components.pos, index_kmer.second.qual);
//            ++kmerCount;


          }
//          counts[tid].c += kmerCount;
        }

      protected:
        int _nprocs;
        int _rank;
        HashFunc _hash;

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
          static_assert(std::is_same<SequenceType,
                        typename QualGenOp::SequenceType>::value,
                        "Kmer Generation and Quality Generation should use the same type");
        }


        void operator()(SequenceType &read, BufferType &buffers) {

          KmerGenOp kmer_op(read.id);
          KmerIter start(read.seq, kmer_op);
          KmerIter end(read.seq_end, kmer_op);

          QualGenOp qual_op;
          QualIter qstart(read.qual, qual_op);
          QualIter qend(read.qual_end, qual_op);

          KmerIndexPairType index_kmer;
//          int kmerCount = 0;

          // NOTE: need to get tid here. depending on how this is called, thread id may not be initialized correctly.
          //int tid = omp_get_thread_num();

          // NOTE: need to call *start to actually evaluate.  question is whether ++ should be doing computation.
          int i = 0;
          for (; i < (K - 1) && (start != end) && (qstart != qend); ++start, ++qstart, ++i);
          for (; (start != end) && (qstart != qend); ++start, ++qstart, ++i)
          {

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
              uint64_t index = this->_hash(index_kmer.first.kmer, index_kmer.second);
//              if (counts[tid] % 100000 == 0)
//                printf("rank %d thread %d hashing to %lu of %lu buffers\n", this->_rank, tid, index, buffers.size()); fflush(stdout);
              //printf("rank %d thread %d hashing %lu to %lu, fill = %lu\n", this->_rank, tid, index_kmer.first.kmer, index, buffers[index].size());
              buffers[index].buffer(index_kmer.first);
              //      printf("kmer send to %lx, key %lx, pos %d, qual %f\n", index_kmer.first, index_kmer.second.kmer, index_kmer.second.id.components.pos, index_kmer.second.qual);
//              ++kmerCount;
            }
            else
            {
              //      printf("BAD kmer quality.  key %lx, pos %d, qual %f\n", index_kmer.second.kmer, index_kmer.second.id.components.pos, index_kmer.second.qual);
            }

          }
//          counts[tid].c += kmerCount;
        }

    };


  } // namespace index

} // namespace bliss



#endif /* KMER_INDEX_GENERATION_H_ */
