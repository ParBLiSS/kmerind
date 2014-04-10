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
    template<typename Derived>
    struct KmerIndexGenerator {
        typedef typename Derived::SequenceType SeqType;
        typedef typename Derived::BufferType BufType;

        void operator()(SeqType &read, int nprocs, int rank,
                        int tid, int j, BufType &buffers, std::vector<size_t> &counts) {
          static_cast<Derived*>(this)->impl(read, nprocs, rank, tid, j, buffers, counts);
        }
    };


    template<typename SendBuffer, typename KmerGenOp>
    class KmerIndexGeneratorNoQuality : public KmerIndexGenerator<KmerIndexGeneratorNoQuality> {
      public:
        typedef typename KmerGenOp::SequenceType		SequenceType;
        typedef std::vector<SendBuffer>  				BufferType;
        typedef bliss::iterator::buffered_transform_iterator<KmerGenOp, typename KmerGenOp::BaseIterType> KmerIter;
        typedef typename KmerGenOp::KmerType 			KmerType;

        KmerIndexGeneratorNoQuality(int nprocs, int rank, int tid) :
            kmer_op(read.id), start(read.seq, kmer_op), end(read.seq_end, kmer_op) {
          static_assert(std::is_same<SequenceType, typename SendBuffer::ValueType>::value, "Kmer Generation and Send Buffer should use the same type");

        }

      protected:
        KmerGenOp kmer_op;
        KmerIter start;
        KmerIter end;
		
		int _nprocs;
		int _rank;
		int _tid;
		int _j;
		
        std::pair<typename kmer_struct_type::kmer_type, kmer_struct_type> index_kmer;
        uint64_t kmerCount;

    };

    template<typename SendBuffer, typename KmerGenOp, typename QualGenOp>
    class KmerIndexGeneratorQuality : public KmerIndexGenerator<KmerIndexGeneratorQuality>,
		KmerIndexGeneratorNoQuality<SendBuffer, KmerGenOp> {
	  public:
		typedef KmerIndexGeneratorNoQuality<SendBuffer, KmerGenOp> BaseType;
        typedef typename BaseType::SequenceType		SequenceType;
        typedef typename BaseType::BufferType		BufferType;
        typedef typename BaseType::KmerIter			KmerIter;
        typedef typename BaseType::KmerType  		KmerType;
        typedef bliss::iterator::buffered_transform_iterator<QualGenOp, typename QualGenOp::BaseIterType> QualIter;

        KmerIndexGeneratorQuality(SequenceType &read, int nprocs, int rank,
            int tid, int j, BufferType &buffers, std::vector<size_t> &counts) :
            KmerIndexGeneratorNoQuality(read, nprocs, rank, tid, j, buffers, counts),
			qual_op(), qstart(read.qual, qual_op), qend(read.qual_end, qual_op)
		{
			static_assert(std::is_same<SequenceType, typename QualGenOp::SequenceType>::value, "Kmer Generation and Quality Generation should use the same type");

        }

      protected:
        QualGenOp qual_op;
        QualIter qstart;
        QualIter qend;

    };


    void compute_kmer_index_from_read()
    {




      int i = -1;
      // NOTE: need to call *start to actually evaluate.  question is whether ++ should be doing computation.
      for (; (start != end) && (qstart != qend); ++start, ++qstart)
      {
        ++i;

        if (i < (K - 1))
          continue;
        //  kmers.push_back(*start);
        //  kmer ^= *start;
        //  qual = *qstart - qual;
        index_kmer = *start;
        index_kmer.second.qual = *qstart;

        // some debugging output
        // printf("kmer send to %lx, key %lx, pos %d, qual %f\n", index_kmer.first, index_kmer.second.kmer, index_kmer.second.id.components.pos, index_kmer.second.qual);

        if (fabs(index_kmer.second.qual) > std::numeric_limits<
                typename kmer_struct_type::qual_type>::epsilon())
        {
          // sending the kmer.
          //printf("rank %d thread %d, staging to buffer %d\n", rank, tid, index_kmer.first % nprocs + (nprocs * tid) );
          buffers[index_kmer.first % nprocs + (nprocs * tid) ].buffer(index_kmer.second);
          //      printf("kmer send to %lx, key %lx, pos %d, qual %f\n", index_kmer.first, index_kmer.second.kmer, index_kmer.second.id.components.pos, index_kmer.second.qual);
          ++kmerCount;
          counts[tid] += 1;
        }
        else
        {
          //      printf("BAD kmer quality.  key %lx, pos %d, qual %f\n", index_kmer.second.kmer, index_kmer.second.id.components.pos, index_kmer.second.qual);
        }
      }

    }


    void store_kmer_index() {

    }


  } // namespace index

} // namespace bliss



#endif /* KMER_FUNCTIONS_H_ */
