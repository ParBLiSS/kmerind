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

#include "mpi.h"

#include <iostream>
//#include <thread>
#include <vector>
#include <unordered_map>
#include <chrono>

#include <string.h>
#include <cstdio>
#include <cmath>

#include "omp.h"

#include "utils/logging.h"
#include "config.hpp"

#include "iterators/range.hpp"
#include "io/fastq_loader.hpp"

#include "common/alphabets.hpp"
#include "common/AlphabetTraits.hpp"
#include "iterators/buffered_transform_iterator.hpp"


template<typename T1, typename T2=double>
struct kmer_struct {
    typedef T1 kmer_type;
    typedef T2 qual_type;

    bliss::iterator::read_id id;
    T1 kmer;
    T2 qual;
};



template<typename ALPHABET, typename Iterator, typename KMER, int K>
struct generate_kmer
{
    KMER kmer;
    typedef typename KMER::kmer_type TO;

    TO revcomp;
    uint16_t pos;

    static constexpr BitSizeType nBits =
        bliss::AlphabetTraits<ALPHABET>::getBitsPerChar();
    static constexpr BitSizeType shift =
        bliss::AlphabetTraits<ALPHABET>::getBitsPerChar() * (K - 1);
    static constexpr AlphabetSizeType max =
        bliss::AlphabetTraits<ALPHABET>::getSize() - 1;

    static constexpr int word_size = sizeof(TO) * 8;
    static constexpr TO mask_reverse = ~(static_cast<TO>(0))
        >> (word_size - shift - nBits);
//    static constexpr TO mask_lower_half = ~(static_cast<TO>(0))
//        >> (word_size - nBits * (K + 1) / 2);


    generate_kmer(const bliss::iterator::read_id &_rid) : kmer(), revcomp(0), pos(0) {
      kmer.kmer = 0;
      kmer.qual = 0;
      kmer.id = _rid;
    }

    size_t operator()(Iterator &iter)
    {
      // store the kmer information.
      char val = ALPHABET::FROM_ASCII[static_cast<size_t>(*iter)];
      kmer.kmer >>= nBits;
      kmer.kmer |= (static_cast<TO>(val) << shift);
      kmer.id.components.pos = pos;

      // generate the rev complement
      char complement = max - val;
      revcomp <<= nBits;
      revcomp |= static_cast<TO>(complement);
      revcomp &= mask_reverse;

      ++iter;
      ++pos;
      return 1;
    }

    std::pair<TO, KMER> operator()()
    {
      // compute the reverse complement
      TO xored = kmer.kmer ^ revcomp;

//      TO xored_recoverable = (xored & ~mask_lower_half) | (_kmer.forward & mask_lower_half);


      return std::pair<TO, KMER>(xored, kmer);
    }

};


/**
 * Phred Scores.  see http://en.wikipedia.org/wiki/FASTQ_format.
 *
 * do it smartly - use metaprogramming to initialize this stuff.
 */
template<typename T>
struct SANGER
{
    typedef T value_type;
    static constexpr T logE_10 = std::log(10.0);
    static_assert(std::is_floating_point<T>::value, "generate_qual output needs to be floating point type");

    static constexpr T offset = 33;

    constexpr T computeLogPCorrect(const int v)
    {
      return std::log(static_cast<T>(1.0) - std::exp(static_cast<T>(v) / static_cast<T>(-10.0) * logE_10));
    }

    static constexpr T logPCorr[94] =
    {
     std::numeric_limits<T>::max(),
     computeLogPCorrect(1),
     computeLogPCorrect(2),
     computeLogPCorrect(3),
     computeLogPCorrect(4),
     computeLogPCorrect(5),
     computeLogPCorrect(6),
     computeLogPCorrect(7),
     computeLogPCorrect(8),
     computeLogPCorrect(9),
     computeLogPCorrect(10),
     computeLogPCorrect(11),
     computeLogPCorrect(12),
     computeLogPCorrect(13),
     computeLogPCorrect(14),
     computeLogPCorrect(15),
     computeLogPCorrect(16),
     computeLogPCorrect(17),
    };
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
 *  TODO: this is currently the bottleneck?  switch to using a look up table instead of computing log all the time.
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
        return std::numeric_limits<TO>::max();   // probability of being correct
      else
        return std::log(static_cast<TO>(1.0) - std::exp(static_cast<TO>(v) / static_cast<TO>(-10.0) * logE_10)); // prob of being correct.
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
      if (newval == std::numeric_limits<TO>::max())
      {
//        printf("ZERO!\n");
        ++zeroCount;
      }
      if (oldval == std::numeric_limits<TO>::max())
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
        // else there is no zero, so valid values.

        if (oldZeroCount == 1)
        {
          //printf("HAD ZEROS!\n");
          // removed a zero.  so recalculate.
          internal = 0.0;
          for (int i = 0; i < K; ++i)
          {
            internal += terms[i];
          }
        }
        else
        {
          // there was not a zero.  so update.
          internal = internal + newval - oldval;
        }

        // compute prob of kmer being incorrect.
        if (fabs(internal) < std::numeric_limits<TO>::epsilon())
        {
          value = 1000.0;
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

#define K 21

typedef kmer_struct<uint64_t, double> kmer_struct_type;

typedef generate_kmer<DNA, char*, kmer_struct_type, K> kmer_op_type;
typedef bliss::iterator::buffered_transform_iterator<kmer_op_type, char*> read_iter_type;


typedef generate_qual<SANGER, char*, double, K> qual_op_type;
qual_op_type qual_op;
typedef bliss::iterator::buffered_transform_iterator<qual_op_type, char*> qual_iter_type;



//void fileio() {
//  // open the file and create a fastq iterator
//    // each dereference has pointers.  choices:
//      // a. directly use.  mmap may be jumping around
//      // b. preload in fastqloader.  jumping around in memory
//      // c. get the pointers and copy the data locally.
//
//  // for each compute thread, accessing next element on file iterator returns a preloaded
//  // data object.  with locking to update pointer.
//
//  // return data object pointing to the local copy.
//  printf("file io\n");
//
//}
//
//


// use a vector of MPIBuffers to manage process'
template<typename T>
class MPISendBuffer {
    // need some default MPI buffer size, then pack in sizeof(T) blocks as many times as possible.

  public:
    static const int END_TAG = 0;


    MPISendBuffer(MPI_Comm _comm, const int _targetRank, const size_t nbytes)
       : accepting(true), targetRank(_targetRank), comm(_comm), send_request(MPI_REQUEST_NULL), total(0)  {
      // initialize internal buffers (double buffering)
      capacity = nbytes / sizeof(T);
      //printf("created buffer for send, capacity = %ld\n", capacity);

      active_buf.reserve(capacity);
      inactive_buf.reserve(capacity);

#ifdef USE_OPENMP
      omp_init_lock(&lock);
#endif
    }

    virtual ~MPISendBuffer() {
      // buffer should be empty at this point.
      assert(active_buf.size() == 0);
      assert(inactive_buf.size() == 0);

      // all previous messages should be done.

#ifdef USE_OPENMP
      omp_destroy_lock(&lock);
#endif

    }

    void buffer(const T & val) {
      assert(accepting);   // if this assert fails, the calling thread's logic is faulty.

      // locking.
#ifdef USE_OPENMP
      omp_set_lock(&lock);
#endif
      // store value
      active_buf.push_back(val);
      // if full, call send (block if other buffer is sending)
      if (active_buf.size() == capacity) {
        send();
      }

#ifdef USE_OPENMP
      omp_unset_lock(&lock);
#endif
    }

    void flush() {
      int rank = 0;
      MPI_Comm_rank(comm, &rank);
      int pid = 0;
#ifdef USE_OPENMP
      pid = omp_get_thread_num();
//
//     assert(pid == 0);   // called by master thread only.
#endif
      accepting = false;

      // no more entries to buffer. send all remaining.
      send();

      // one more time to wait for prev send to finish.
      send();

      // send termination signal (send empty message).
      //printf("%d done sending to %d\n", rank, targetRank); fflush(stdout);
      MPI_Send(nullptr, 0, MPI_INT, targetRank, END_TAG, comm);
      printf("%d.%d  %ld flushed.\n", rank, pid, total); fflush(stdout);
    }

  protected:


    // 2 buffers per target
    std::vector<T> active_buf;
    std::vector<T> inactive_buf;

    // if buffer is still accepting stuff
    bool accepting;

    // rank of target entry.
    const int targetRank;
    size_t capacity;

    MPI_Comm comm;
    // inactive buffer send status
    MPI_Request send_request;

    size_t total;

#ifdef USE_OPENMP
      omp_lock_t lock;
#endif


    void send() {
      int rank = 0;
      MPI_Comm_rank(comm, &rank);
      int pid = 0;
#ifdef USE_OPENMP
      pid = omp_get_thread_num();
#endif
      // only called from within locked code block
      std::chrono::high_resolution_clock::time_point time1, time2;
      std::chrono::duration<double> time_span;

      // check inactive buffer is sending?
      if (send_request != MPI_REQUEST_NULL) {
        // if being sent, wait for that to complete
        //printf("Waiting for prev send %d -> %d to finish\n", rank, targetRank); fflush(stdout);
        //time1 = std::chrono::high_resolution_clock::now();Traffic

        MPI_Wait(&send_request, MPI_STATUS_IGNORE);
        //time2 = std::chrono::high_resolution_clock::now();
        //time_span = std::chrono::duration_cast<std::chrono::duration<double>>(
        //    time2 - time1);

        //printf("Waiting for prev send %d -> %d finished in %f\n", rank, targetRank, time_span.count()); fflush(stdout);
      }
      // clear the inactive buf
      inactive_buf.clear();
      assert(inactive_buf.size() == 0);
      int s = active_buf.size();
      // swap active and inactive buffer  (this swaps the contents)
      inactive_buf.swap(active_buf);
      assert(active_buf.size() == 0);
      assert(inactive_buf.size() == s);

      // async send inactive buffer if it's not empty.  Requires that remote side uses probe to get the size of data first.

      // TODO: send the data as bytes.  This assumes the receive end knows how to parse the elements
      if (inactive_buf.size() > 0) {
        total += inactive_buf.size();
        //printf("%d sending to %d\n", rank, targetRank);
        MPI_Isend(inactive_buf.data(), inactive_buf.size() * sizeof(T), MPI_UNSIGNED_CHAR, targetRank, pid + 1, comm, &send_request);
        //printf("SEND %d.%d send to %d, number of entries %ld, total %ld\n", rank, pid, targetRank, inactive_buf.size(), total); fflush(stdout);

      }

    }
};



// this can become a 1 to n transformer???
void compute(bliss::iterator::fastq_sequence<char*> &read, int nprocs, int rank, int pid, int j, std::vector<MPISendBuffer<kmer_struct_type> > &buffers ) {

  kmer_op_type kmer_op(read.id);
  read_iter_type start(read.seq, kmer_op);
  read_iter_type end(read.seq_end, kmer_op);

  qual_iter_type qstart(read.qual, qual_op);
  qual_iter_type qend(read.qual_end, qual_op);


  std::pair<kmer_struct_type::kmer_type, kmer_struct_type> index_kmer;

  uint64_t kmerCount = 0;

  int i = -1;
  // NOTE: need to call *start to actually evaluate.  question is whether ++ should be doing computation.
  for (; (start != end) && (qstart != qend); ++start, ++qstart)
  {
    ++i;

    if (i < (K - 1))
      continue;
//        kmers.push_back(*start);
//    kmer ^= *start;
//    qual = *qstart - qual;
    index_kmer = *start;
    index_kmer.second.qual = *qstart;

    // some debugging output
    //printf("kmer send to %lx, key %lx, pos %d, qual %f\n", index_kmer.first, index_kmer.second.kmer, index_kmer.second.id.components.pos, index_kmer.second.qual);

    // sending the kmer.
    buffers[index_kmer.first % nprocs].buffer(index_kmer.second);

    ++kmerCount;
  }

}


typedef std::unordered_multimap<kmer_struct_type::kmer_type, kmer_struct_type> kmer_map_type;
void networkread(MPI_Comm comm, const int nprocs, const int rank, const size_t buf_size, const int senders, kmer_map_type &kmers) {

  // track how many active senders remain.
  int n_senders = nprocs;  // hack.  each proc has multiple threads.
  MPI_Status status;

  size_t capacity = buf_size / sizeof(kmer_struct_type);
  kmer_struct_type *array = new kmer_struct_type[capacity];
  //printf("created temp storage for read, capacity = %ld\n", capacity); fflush(stdout);
  memset(array, 0, capacity * sizeof(kmer_struct_type));
  int received = 0;
  int count = 0;
  int src;
  int tag;

  while (n_senders > 0) {
    // probe for a message.  if empty message, then don't need to listen on that node anymore
    //printf("probing...\n");
    MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, comm, &status);

    src = status.MPI_SOURCE;
    tag = status.MPI_TAG;
    MPI_Get_count(&status, MPI_UNSIGNED_CHAR, &received);   // length is always > 0?

    if (tag == MPISendBuffer<kmer_struct_type>::END_TAG) {
      // end of messaging.
      printf("RECV %d receiving END signal %d from %d\n", rank, received, src); fflush(stdout);
      MPI_Recv(reinterpret_cast<unsigned char*>(array), received, MPI_UNSIGNED_CHAR, src, tag, comm, MPI_STATUS_IGNORE);
      --n_senders;
      continue;
    }

    //printf("RECV %d receiving %d bytes from %d.%d\n", rank, received, src, tag - 1); fflush(stdout);

    // TODO: FIX segfault here at iter 2.
    MPI_Recv(reinterpret_cast<unsigned char*>(array), received, MPI_UNSIGNED_CHAR, src, tag, comm, MPI_STATUS_IGNORE);

    // if message size is 0, then no work to do.
    if (received == 0)
      continue;

    assert(received % sizeof(kmer_struct_type) == 0);  // no partial messages.
    count = received / sizeof(kmer_struct_type);
    assert(count > 0);
    //printf("RECV+ %d receiving %d bytes or %d records from %d.%d\n", rank, received, count, src, tag - 1); fflush(stdout);


    // now process the array
    for (int i = 0; i < count; ++i) {
      kmers.insert(kmer_map_type::value_type(array[i].kmer, array[i]));
    }

//    memset(array, 0, capacity * sizeof(kmer_struct_type));
  }

  delete [] array;

  printf("network read done\n"); fflush(stdout);
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
  int provided;

#ifdef USE_OPENMP
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

  if (provided < MPI_THREAD_MULTIPLE) {
    printf("ERROR: The MPI Library Does not have full thread support.\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
#else
  MPI_Init(&argc, &argv);
#endif

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



    std::vector<MPISendBuffer<kmer_struct_type> > buffers;
    for (int i = 0; i < nprocs; ++i) {
      buffers.push_back(std::move(MPISendBuffer<kmer_struct_type>(MPI_COMM_WORLD, i, 8192*1024)));
    }


    std::unordered_multimap<kmer_struct_type::kmer_type, kmer_struct_type> index;


  //  std::thread fileio_t(fileio);
//    std::thread networkwrite_t(networkwrite);
    //std::thread networkread_t(networkread);



    // do some work using openmp
    t1 = std::chrono::high_resolution_clock::now();
    bliss::io::fastq_loader::iterator fastq_start = loader.begin();
    bliss::io::fastq_loader::iterator fastq_end = loader.end();

    int num_threads = 1;
#ifdef USE_OPENMP
    num_threads = omp_get_max_threads();
    assert(num_threads >= 3);
    fprintf(stderr, "rank %d MAX THRADS = %d\n", rank, num_threads);
    omp_set_nested(1);
    omp_set_dynamic(0);
#endif

    int senders = 0;
    MPI_Allreduce(&num_threads, &senders, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    bliss::iterator::fastq_sequence<char*> read;
    int i = 0;
#pragma omp parallel sections num_threads(2)
    {


#pragma omp section
    {
      networkread(MPI_COMM_WORLD, nprocs, rank, 8192*1024, senders, index);
      printf("kmer recording completed.  records: %ld\n", index.size());
    } // network read thread

#pragma omp section
    {  // compute threads
#pragma omp parallel num_threads(num_threads - 1)
#pragma omp single nowait
    {
      for (; fastq_start != fastq_end; ++fastq_start, ++i) {
        // get data, and move to next for the other threads

        // first get read
        read = *fastq_start;
//
//            ++fastq_start;
//            ++i;
#pragma omp task firstprivate(read, i)
        {
          // copy read.  not doing that right now.

          // then compute
          compute(read, nprocs, rank, omp_get_thread_num(), i, buffers);

          if (i % 100000 == 0)
            printf("%d\n", i);
          // release resource
        }

        // next iteration will check to see if the iterator is at end,
        // and if not, get and compute.
      }
#pragma omp taskwait


      t2 = std::chrono::high_resolution_clock::now();
      time_span = std::chrono::duration_cast<std::chrono::duration<double>>(
          t2 - t1);
      INFO("computation rank " << rank << " elapsed time: " << time_span.count() << "s.");


      // send the last part out.
      for (int i = 0; i < nprocs; ++i) {
        buffers[(i+ rank) % nprocs].flush();
      }

    }
    printf("compute completed\n");



    }  // compute threads

    } // outer parallel sections




    //fileio_t.join();
//    networkwrite_t.join();
    //networkread_t.join();
  }





#ifdef USE_MPI
    MPI_Finalize();
#endif

  return 0;
}
