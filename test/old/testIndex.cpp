/**
 *
 * test generating the indices
 *
 * 1,2,3,4 MPI procs, each reading a different part.
 *
 * step 1.  parse fastq into reads
 * setp 2.  for each read, iterate and generate kmers (and reverse complement)
 * step 3.  for each read, concurrently generate quality score
 * step 4.  bin the kmers for redistribution
 * step 5.  MPI send kmers.  (threading?)
 *
 */

#include "mpi.h"


#include <string>     // for std::string
#include <sys/stat.h>  // for stat
#include <cassert>
#include <iostream>
#include <chrono>

#include <string.h>
#include <vector>
#include <bitset>
#include <cmath>
#include <sstream>


#include "utils/logging.h"
#include "bliss-config.hpp"

#include "iterators/range.hpp"
#include "io/fastq_loader.hpp"

#include "common/alphabets.hpp"
#include "common/alphabet_traits.hpp"
#include "iterators/buffered_transform_iterator.hpp"

template<typename ALPHABET, typename Iterator, typename TO, int K>
struct generate_kmers
{
    TO forward;
    TO reverse;
    TO xored;
    TO xoredRecoverable;

    static constexpr BitSizeType nBits =
        bliss::AlphabetTraits<ALPHABET>::getBitsPerChar();
    static constexpr BitSizeType shift =
        bliss::AlphabetTraits<ALPHABET>::getBitsPerChar() * (K - 1);
    static constexpr AlphabetSizeType max =
        bliss::AlphabetTraits<ALPHABET>::getSize() - 1;

    static constexpr int word_size = sizeof(TO) * 8;
    static constexpr TO mask_reverse = ~(static_cast<TO>(0))
        >> (word_size - shift - nBits);
    static constexpr TO mask_lower_half = ~(static_cast<TO>(0))
        >> (word_size - nBits * (K + 1) / 2);

    generate_kmers()
        : forward()
    {
    }

    size_t operator()(Iterator &iter)
    {
      char val = ALPHABET::FROM_ASCII[static_cast<size_t>(*iter)];
      forward >>= nBits;
      forward |= (static_cast<TO>(val) << shift);

      char complement = max - val;
      reverse <<= nBits;
      reverse |= static_cast<TO>(complement);
      reverse &= mask_reverse;

      xored = forward ^ reverse;
//      INFO( "kmer:\t" << std::bitset<word_size>(forward) << std::endl << "\t" << std::bitset<word_size>(reverse) );
//      INFO( "\t" << std::bitset<word_size>(xored) );

      xoredRecoverable = (xored & ~mask_lower_half)
          | (forward & mask_lower_half);

      ++iter;
      return 1;
    }

    TO operator()()
    {
      return xoredRecoverable;
    }

};

/**
 * compute kmer quality based on phred quality score.
 *
 * Q = -10 log_10 P
 *
 * where P is probability of error.
 *
 * probability of kmer being correct is: (1-p_1)(1-p_2)...(1-p_k)
 * probability of kmer being incorrect then is 1 - (p(correct...))
 *
 * Phrad adds the Phred score, which amounts to p_1 * p_2 * ... * p_k because of the logarithm.
 *  (probability that all bases in the sequence are incorrect.)
 * Quake is computes a kmer score based on probability in the phred score.  not clear exactly how.
 *
 * value range is in ascii from !(33) to ~(126), using sanger encoding.
 *
 * // now computing k-correct with approximation gives
 * //    1 - sum_1..k(p_i) + 1/2 sum_1..k sum_1..k (p_i)(p_j) - 3rd order term + 4th order term ...
 * //
 * //   p(incorrect) is sum_1..k(p_i) - 1/2 sum_1..k sum_1..k (p_i)(p_j) + 3rd order term - 4th order term ...
 * //
 * // =====>>>  max_1..k(p_i) < p(incorrect) < sum_1..k(p_i).
 * //
 * // approximate with linear terms.  (i.e. probability of 2 or more bases being called incorrectly is low).  bound is not that tight.
 * //
 * // =====>>>  max_1..k(p_i) < p(
 * //
 * // max_1..k(p_i) ~= min_1..k(q_i)
 *
 *
 *  each term in (1-p_i) is (1- 10^(q/-10)) = (1 - e ^ ((log 10)(q/-10)))
 *  accumulate using sum(log(term))., expand using e^(sum(log(term))), calculate -10 log_10(1-cumulative).
 *
 */
template<typename ENCODING, typename Iterator, typename TO, int K>
struct generate_qual
{
    typedef typename std::iterator_traits<Iterator>::value_type TI;

    TO value;
    TO internal;
    TO terms[K];
    int pos;
    int zeroCount;

    static constexpr TO logE_10 = std::log(10);
    static_assert(std::is_floating_point<TO>::value, "generate_qual output needs to be floating point type");

    generate_qual()
        : value(1), pos(0)
    {
      for (int i = 0; i < K; ++i)
      {
        terms[i] = 0.0;
      }
      zeroCount = K;
    }

    TO compute(TI v)
    {
      if (v == 0)
        return 0.0;   // probability of being correct
      else
        return 1.0 - std::exp(static_cast<TO>(v) / -10.0f * logE_10); // prob of being correct.
    }

    size_t operator()(Iterator &iter)
    {

      // drop the old value
      TO oldval = terms[pos];

      // add the new value,       // update the position  - circular queue
      TO newval = compute(*iter - 33);  // this is for Sanger encoding.
      terms[pos] = newval;
      pos = (pos + 1) % K;

      // save the old zero count.
      int oldZeroCount = zeroCount;

      // update the zero count
      if (newval == 0.0)
      {
//        INFOF("ZERO!\n");
        ++zeroCount;
      }
      if (oldval == 0.0)
      {
        --zeroCount;
      }

      // if any is zero, then return 0.
      if (zeroCount > 0)
      {
        internal = 0.0;
        value = 0.0;
      }
      else
      {
        // else there is no zero,

        if (oldZeroCount == 1)
        {
          //INFOF("HAD ZEROS!\n");
          // removed a zero.  so recalculate.
          internal = 0.0;
          for (int i = 0; i < K; ++i)
          {
            internal += std::log(terms[i]);
          }
        }
        else
        {
          // there was not a zero.  so update.
          internal = internal + std::log(newval) - std::log(oldval);
        }

        TO v = std::exp(internal);  // prob of being correct.
        if (v == 1.0)
        {
          value = 1000;
          INFOF("confident kmer!\n");
        }
        else
        {
          value = -10.0 * std::log10(1.0 - std::exp(internal));
        }
      }

//      INFOF("qual %d oldval = %f, newval = %f, internal = %f, output val = %f\n", *iter, oldval, newval, internal, value);

      // move the iterator.
      ++iter;
      return value;
    }

    TO operator()()
    {
      return value;
    }

};

#define K 11
typedef bliss::io::fastq_loader<DNA, float> FileLoaderType;

typedef bliss::index::KmerSize<21> KmerSize;
typedef uint64_t KmerType;
typedef float QualityType;
typedef bliss::index::KmerIndexElementWithIdAndQuality<KmerSize, KmerType, bliss::common::ShortSequenceKmerId, QualityType> KmerIndexType;

typedef CharType* BaseIterType;

typedef DNA Alphabet;
typedef bliss::io::Sequence<BaseIterType>  SequenceType;

typedef bliss::index::generate_kmer<SequenceType, KmerIndexType> kmer_op_type;
typedef bliss::index::generate_qual<SequenceType, KmerSize, bliss::index::SangerToLogProbCorrect<double>, QualityType> qual_op_type;



int main(int argc, char* argv[])
{
  LOG_INIT();

  std::string filename("/home/tpan/src/bliss/test/data/test.fastq");
//  std::string filename("/mnt/data/1000genome/HG00096/sequence_read/SRR077487_1.filt.fastq");

  if (argc > 1)
  {
    filename.assign(argv[1]);
  }

  int rank = 0, nprocs = 1;
#ifdef USE_MPI
  // initialize MPI
  MPI_Init(&argc, &argv);

  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == 0)
    INFO( "USE_MPI is set" );
#endif

  // first thread gets the file size.
  uint64_t file_size = 0;
  if (rank == 0)
  {
    struct stat filestat;
    stat(filename.c_str(), &filestat);
    file_size = static_cast<uint64_t>(filestat.st_size);
    INFOF("block size is %ld\n", filestat.st_blksize);
    INFOF("sysconf block size is %ld\n", sysconf(_SC_PAGE_SIZE));
  }

#ifdef USE_MPI
  if (nprocs > 1) {
    // broadcast to all
    MPI_Bcast(&file_size, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
  }
#endif

  if (rank == nprocs - 1)
  {
    INFOF("file size is %ld\n", file_size);
  }
  /////////////// now try to open the file

  // real data:  mmap is better for large files and limited memory.
  //             preloading is better for smaller files and/or large amount of memory.
  // stream processing means data does not need to be buffered in memory - more efficient.

  // file access:  better to work with a few pages at a time, or to work with large block?

  // first generate an approximate partition.
  FileLoaderType::RangeType r =
      FileLoaderType::RangeType::block_partition(nprocs,
                                                          rank, 0, file_size);
  INFO( rank << " equipart: " << r );
  INFO( rank << " test block aligned: "
            << r.align_to_page(sysconf(_SC_PAGE_SIZE)) );

  // now open the file and begin reading
  std::chrono::high_resolution_clock::time_point t1, t2;
  std::chrono::duration<double> time_span1;
  std::chrono::duration<double> time_span2;
  std::chrono::duration<double> time_span3;
  std::chrono::duration<double> time_span;

  {
    t1 = std::chrono::high_resolution_clock::now();
    // now open the file
    FileLoaderType loader(filename, r, file_size);
//    FileLoaderType loader(filename, file_size);
    t2 = std::chrono::high_resolution_clock::now();
    time_span3 = std::chrono::duration_cast<std::chrono::duration<double>>(
        t2 - t1);
    INFO("MMap rank " << rank << " elapsed time: " << time_span3.count() << "s.");

    INFO( rank << " record adjusted " << loader.getRange() );

    t1 = std::chrono::high_resolution_clock::now();
    typename FileLoaderType::IteratorType fastq_start = loader.begin();
    typename FileLoaderType::IteratorType fastq_end = loader.end();

    uint64_t id = 0;

    uint64_t readCount = 0;
    for (; fastq_start != fastq_end; ++fastq_start)
    {

      id ^= (*fastq_start).id;
      ++readCount;
    }
    t2 = std::chrono::high_resolution_clock::now();
    time_span3 = std::chrono::duration_cast<std::chrono::duration<double>>(
        t2 - t1);
    INFO("reads " << readCount << " rank " << rank << " elapsed time: " << time_span3.count() << "s.");

    INFOF("avoid compiler optimizing out the ops %ld\n", id);

    t1 = std::chrono::high_resolution_clock::now();

    fastq_start = loader.begin();
    fastq_end = loader.end();

    typedef generate_kmers<DNA, char*, uint64_t, K> op_type;
    op_type kmer_op;
    typedef bliss::iterator::buffered_transform_iterator<op_type, char*> read_iter_type;

    uint64_t kmer = 0;
    uint64_t baseCount = 0;

    for (; fastq_start != fastq_end; ++fastq_start)
    {
      read_iter_type start((*fastq_start).seq, kmer_op);
      read_iter_type end((*fastq_start).seq_end, kmer_op);

      int i = -1;
      // NOTE: need to call *start to actually evaluate.  question is whether ++ should be doing computation.
      for (; (start != end); ++start)
      {
        ++i;
        ++baseCount;
        if (i < (K - 1))
          continue;
        kmer ^= *start;

      }
    }
    t2 = std::chrono::high_resolution_clock::now();
    time_span3 = std::chrono::duration_cast<std::chrono::duration<double>>(
        t2 - t1);
    INFO("bases " << baseCount << " kmer generation rank " << rank << " elapsed time: " << time_span3.count() << "s.");

    INFOF("avoid compiler optimizing out the ops %ld\n", kmer);

    t1 = std::chrono::high_resolution_clock::now();

    struct SANGER
    {
    };
    typedef generate_qual<SANGER, char*, double, K> qual_op_type;
    qual_op_type qual_op;
    typedef bliss::iterator::buffered_transform_iterator<qual_op_type, char*> qual_iter_type;

    kmer = 0;
    double qual = 0;
    uint64_t kmerCount = 0;

    fastq_start = loader.begin();
    fastq_end = loader.end();


    for (; fastq_start != fastq_end; ++fastq_start)
    {
      read_iter_type start((*fastq_start).seq, kmer_op);
      read_iter_type end((*fastq_start).seq_end, kmer_op);

      qual_iter_type qstart((*fastq_start).qual, qual_op);
      qual_iter_type qend((*fastq_start).qual_end, qual_op);

      int i = -1;
      // NOTE: need to call *start to actually evaluate.  question is whether ++ should be doing computation.
      for (; (start != end) && (qstart != qend); ++start, ++qstart)
      {
        ++i;

        if (i < (K - 1))
          continue;
//        kmers.push_back(*start);
        kmer ^= *start;
        qual = *qstart - qual;
        ++kmerCount;
      }

    }
    t2 = std::chrono::high_resolution_clock::now();
    time_span3 = std::chrono::duration_cast<std::chrono::duration<double>>(
        t2 - t1);
    INFO("kmer + qual " << kmerCount << " generation rank " << rank << " elapsed time: " << time_span3.count() << "s.");

    INFOF("avoid compiler optimizing out the ops %lx %lf\n", kmer, qual);

  }

// MEMORY INTENSIVE BELOW

//  {
//    t1 = std::chrono::high_resolution_clock::now();
//    // now open the file
//    bliss::io::fastq_loader loader(filename, r, file_size, true);
//    t2 = std::chrono::high_resolution_clock::now();
//    time_span3 = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
//    INFO("MMap preload " << "elapsed time: " << time_span3.count() << "s.");
//
//    INFO( rank << " record adjusted preload " << loader.getRange() );
//
//
//    t1 = std::chrono::high_resolution_clock::now();
//
//    bliss::io::fastq_loader::iterator fastq_start = loader.begin();
//    bliss::io::fastq_loader::iterator fastq_end = loader.end();
//
//    typedef generate_kmers<DNA, char*, uint64_t, 21>  op_type;
//    op_type kmer_op;
//    typedef bliss::iterator::buffered_transform_iterator<op_type, char*> read_iter_type;
//
//    // transform and generate kmers
////    std::vector<uint64_t> kmers;
//    uint64_t kmer;
//    for (; fastq_start != fastq_end; ++fastq_start) {
//      read_iter_type start((*fastq_start).seq, kmer_op);
//      read_iter_type end((*fastq_start).seq_end, kmer_op);
//
//      int i = -1;
//      std::stringstream ss;
//      for (; start != end; ++start) {
//        ++i;
//        ss.str(std::string());
//        if (i < 20) continue;
//  //      kmers.push_back(*start);
//
//        kmer = *start;
//
//        ss << kmer;
//      }
//
//    }
//    t2 = std::chrono::high_resolution_clock::now();
//    time_span3 = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
//    INFO("kmer generation preload " << "elapsed time: " << time_span3.count() << "s.");
//
//
//  }

#ifdef USE_MPI
    MPI_Finalize();
#endif

  return 0;

}
