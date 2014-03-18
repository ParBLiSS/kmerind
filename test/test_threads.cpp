/**
 * @file		test_threads.cpp
 * @ingroup
 * @author	tpan
 * @brief
 * @details
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */

#include <iostream>
#include <thread>

#include <cstdio>
#include <cmath>

#include "mpi.h"
#include "omp.h"

#include "utils/logging.h"
#include "config.hpp"

#include "iterators/range.hpp"
#include "io/fastq_loader.hpp"

#include "common/alphabets.hpp"
#include "common/AlphabetTraits.hpp"
#include "iterators/buffered_transform_iterator.hpp"

template<typename TO>
struct kmer_struct {
    TO kmer;
    TO revcomp;
    bliss::iterator::read_id id;
};



template<typename ALPHABET, typename Iterator, typename TO, int K>
struct generate_kmer
{
    kmer_struct<TO, TO> kmer;

    static constexpr BitSizeType nBits =
        bliss::AlphabetTraits<ALPHABET>::getBitsPerChar();
    static constexpr BitSizeType shift =
        bliss::AlphabetTraits<ALPHABET>::getBitsPerChar() * (K - 1);
    static constexpr AlphabetSizeType max =
        bliss::AlphabetTraits<ALPHABET>::getSize() - 1;

    static constexpr int word_size = sizeof(TO) * 8;
    static constexpr TO mask_reverse = ~(static_cast<TO>(0))
        >> (word_size - shift - nBits);

    generate_kmer() : kmer() {
      kmer.first = 0;
      kmer.second = 0;
    }

    size_t operator()(Iterator &iter)
    {
      char val = ALPHABET::FROM_ASCII[static_cast<size_t>(*iter)];
      kmer.first >>= nBits;
      kmer.first |= (static_cast<TO>(val) << shift);

      char complement = max - val;
      kmer.second <<= nBits;
      kmer.second |= static_cast<TO>(complement);
      kmer.second &= mask_reverse;

      ++iter;
      return 1;
    }

    std::pair<TO, TO> operator()()
    {
      return kmer;
    }

};



template<typename TO, int K>
struct xor_kmer
{
    TO operator()(const std::pair<TO, TO> &kmer_rev) {
      return kmer_rev.first ^ kmer_rev.second;
    }
};

template<typename ALPHABET, typename TO, int K>
struct xor_kmer_recoverable
{
    static constexpr BitSizeType nBits =
        bliss::AlphabetTraits<ALPHABET>::getBitsPerChar();
    static constexpr BitSizeType shift =
        bliss::AlphabetTraits<ALPHABET>::getBitsPerChar() * (K - 1);
    static constexpr int word_size = sizeof(TO) * 8;

    static constexpr TO mask_reverse = ~(static_cast<TO>(0))
        >> (word_size - shift - nBits);
    static constexpr TO mask_lower_half = ~(static_cast<TO>(0))
        >> (word_size - nBits * (K + 1) / 2);

    TO operator()(const std::pair<TO, TO> &kmer_rev) {
      TO xored = kmer_rev.first ^ kmer_rev.second;

      return (xored & ~mask_lower_half) | (kmer_rev.first & mask_lower_half);
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
//        printf("ZERO!\n");
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
          //printf("HAD ZEROS!\n");
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
          printf("confident kmer!\n");
        }
        else
        {
          value = -10.0 * std::log10(1.0 - std::exp(internal));
        }
      }

//      printf("qual %d oldval = %f, newval = %f, internal = %f, output val = %f\n", *iter, oldval, newval, internal, value);

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


typedef generate_kmer<DNA, char*, uint64_t, K> op_type;
op_type kmer_op;
typedef bliss::iterator::buffered_transform_iterator<op_type, char*> read_iter_type;

struct SANGER
{
};
typedef generate_qual<SANGER, char*, double, K> qual_op_type;
qual_op_type qual_op;
typedef bliss::iterator::buffered_transform_iterator<qual_op_type, char*> qual_iter_type;



void fileio() {
  // open the file and create a fastq iterator
    // each dereference has pointers.  choices:
      // a. directly use.  mmap may be jumping around
      // b. preload in fastqloader.  jumping around in memory
      // c. get the pointers and copy the data locally.

  // for each compute thread, accessing next element on file iterator returns a preloaded
  // data object.  with locking to update pointer.

  // return data object pointing to the local copy.
  printf("file io\n");

}

template<typename TO>
struct mpi_kmer_struct {
    TO kmer;
    bliss::iterator::read_id id;
};


// this can become a 1 to n transformer???
void compute(int rank, int pid, bliss::iterator::fastq_sequence<char*> &read, int j) {


  read_iter_type start(read.seq, kmer_op);
  read_iter_type end(read.seq_end, kmer_op);

  qual_iter_type qstart(read.qual, qual_op);
  qual_iter_type qend(read.qual_end, qual_op);

  uint64_t kmer = 0;
  double qual = 0;
  uint64_t kmerCount = 0;

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


// use a vector of MPIBuffers to manage
template<typename T>
class MPIBuffer {
    // need some default MPI buffer size, then pack in sizeof(T) blocks as many times as possible.

  public:
    MPIBuffer() {
      // initialize internal buffers (double buffering)
    }


    bool buffer(const uint32_t id, const T & val) {
      // synchronized

      // store value

      // if full, call send (block if other buffer is sending)
    }

    bool flush() {
      // no more coming in.  call by master thread only.
    }

  protected:
    // 2 buffers per target
    // active buffer content count
    // inactive buffer send status

    bool send(const uint32_t id) {
      // synchronized.

      // check inactive buffer status.
      // if being sent, wait for that to complete

      // swap active and inactive buffer


      // async send full inactive buffer.

    }
};

void networkwrite() {
  // instantiate bins (m per proc)

  // atomic add to bin. update bin count

  // if bin count is at limit, MPI send length, then send data
  //
  // final write
    // flush: sends length
      // if length > 0, send data
      // else length = 0, don't send anything.  return to waiting.

  // send termination signal (send empty message).
  printf("network write\n");

}

void networkread() {
  // initiate an array of senders, mark all as 1.

  // iprobe for a message.  if empty message, then don't need to listen on that node anymore

  // if has length,
    // if length > 0, then receive data and do something
  // else return to probe.

  printf("network read\n");
}


int main(int argc, char** argv) {

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
    std::cout << "USE_MPI is set" << std::endl;
#endif

  // first thread gets the file size.
  uint64_t file_size = 0;
  if (rank == 0)
  {
    struct stat filestat;
    stat(filename.c_str(), &filestat);
    file_size = static_cast<uint64_t>(filestat.st_size);
    fprintf(stderr, "block size is %ld\n", filestat.st_blksize);
    fprintf(stderr, "sysconf block size is %ld\n", sysconf(_SC_PAGE_SIZE));
  }

#ifdef USE_MPI
  if (nprocs > 1) {
    // broadcast to all
    MPI_Bcast(&file_size, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
  }
#endif

  if (rank == nprocs - 1)
  {
    fprintf(stderr, "file size is %ld\n", file_size);
  }
  /////////////// now try to open the file

  // real data:  mmap is better for large files and limited memory.
  //             preloading is better for smaller files and/or large amount of memory.
  // stream processing means data does not need to be buffered in memory - more efficient.

  // file access:  better to work with a few pages at a time, or to work with large block?

  // first generate an approximate partition.
  bliss::io::file_loader::range_type r =
      bliss::io::file_loader::range_type::block_partition(nprocs,
                                                          rank, 0, file_size);
  std::cout << rank << " equipart: " << r << std::endl;
  std::cout << rank << " test block aligned: "
            << r.align_to_page(sysconf(_SC_PAGE_SIZE)) << std::endl;

  // now open the file and begin reading
  std::chrono::high_resolution_clock::time_point t1, t2;
  std::chrono::duration<double> time_span;

  {
    t1 = std::chrono::high_resolution_clock::now();
    // now open the file
    bliss::io::fastq_loader loader(filename, r, file_size);
//    bliss::io::file_loader loader(filename, file_size);
    t2 = std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(
        t2 - t1);
    INFO("MMap rank " << rank << " elapsed time: " << time_span.count() << "s.");

    std::cout << rank << " record adjusted " << loader.getRange() << std::endl;





    std::cout << "starting threads" << std::endl;

  //  std::thread fileio_t(fileio);
    std::thread networkwrite_t(networkwrite);
    std::thread networkread_t(networkread);

    std::cout << "threads started" << std::endl;

    // do some work using openmp
    t1 = std::chrono::high_resolution_clock::now();
    bliss::io::fastq_loader::iterator fastq_start = loader.begin();
    bliss::io::fastq_loader::iterator fastq_end = loader.end();


    int threadcount = 4;

    bliss::iterator::fastq_sequence<char*> read;
    int i = 0;
#pragma omp parallel shared(fastq_start, fastq_end, i, read) num_threads(threadcount)
#pragma omp single nowait
    {
      for (; fastq_start != fastq_end; ++fastq_start, ++i) {
        // get data, and move to next for the other threads

//          // first get read
            read = *fastq_start;
//
//            ++fastq_start;
//            ++i;
#pragma omp task firstprivate(read, i)
        {
          // copy read.  not doing that right now.

          // then compute
          compute(rank, omp_get_thread_num(), read, i);

          // release resource

        }

        // next iteration will check to see if the iterator is at end,
        // and if not, get and compute.
      }
#pragma omp taskwait
    }

    t2 = std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(
        t2 - t1);
    INFO("computation rank " << rank << " elapsed time: " << time_span.count() << "s.");

    //fileio_t.join();
    networkwrite_t.join();
    networkread_t.join();
  }
  std::cout << "threads done" << std::endl;

#ifdef USE_MPI
    MPI_Finalize();
#endif

  return 0;
}
