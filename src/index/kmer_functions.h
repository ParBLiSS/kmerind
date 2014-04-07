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
    // FASTQ_SEQUENCE should be parameterized with base iterator and alphabet
    // given those + k, should be able to get kmer index element type
    // given kmer index element type, and fastq sequence type, should be able to specialize kmer generator type.
    // given the kmer index element type, should be able to define sendbuffer type.

    // this can become a 1 to n transformer???
    template<typename SendBuffer, typename KmerGenOp, typename QualGenOp>
    void compute_kmer_index_from_read(bliss::iterator::fastq_sequence<char*> &read, int nprocs, int rank,
                 int tid, int j, std::vector<SendBuffer> &buffers, std::vector<size_t> &counts)
    {

      KmerGenOp kmer_op(read.id);
      typedef bliss::iterator::buffered_transform_iterator<KmerGenOp, typename KmerGenOp::BaseIterType> KmerIter;
      KmerIter start(read.seq, kmer_op);
      KmerIter end(read.seq_end, kmer_op);

      QualGenOp qual_op;
      typedef bliss::iterator::buffered_transform_iterator<QualGenOp, typename QualGenOp::BaseIterType> QualIter;
      QualIter qstart(read.qual, qual_op);
      QualIter qend(read.qual_end, qual_op);

      typedef typename KmerGenOp::KmerType kmer_struct_type;
      std::pair<typename kmer_struct_type::kmer_type, kmer_struct_type> index_kmer;

      uint64_t kmerCount = 0;

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
