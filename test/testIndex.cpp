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

#include <string>     // for std::string
#include <sys/stat.h>  // for stat
#include <cassert>
#include <iostream>
#include <chrono>

#include "mpi.h"

#include "utils/logging.h"
#include "config.hpp"

#include "iterators/range.hpp"
#include "io/file_loader.hpp"
#include "io/fastq_loader.hpp"

#include <string.h>
#include <vector>
#include <bitset>

#include "common/alphabets.hpp"
#include "common/AlphabetTraits.hpp"
//#include "iterators/transform_iterator.hpp"
#include "iterators/buffered_transform_iterator.hpp"

template<typename ALPHABET, typename Iterator, typename TO, int K >
struct generate_kmers {
    TO forward;
    TO reverse;
    TO xored;
    static constexpr BitSizeType nBits = bliss::AlphabetTraits<ALPHABET>::getBitsPerChar();
    static constexpr BitSizeType shift = bliss::AlphabetTraits<ALPHABET>::getBitsPerChar() * (K - 1);
    static constexpr AlphabetSizeType max = bliss::AlphabetTraits<ALPHABET>::getSize() - 1;

    static constexpr int word_size = sizeof(TO) * 8;
    static constexpr TO mask_reverse =  ~(static_cast<TO>(0))  >> (word_size - shift - nBits);

    generate_kmers() : forward()
    { }

    size_t operator()(Iterator &iter) {
      char val = ALPHABET::FROM_ASCII[static_cast<size_t>(*iter)];
      forward >>= nBits;
      forward |= (static_cast<TO>(val) << shift);

      char complement = max - val;
      reverse <<= nBits;
      reverse |= static_cast<TO>(complement);
      reverse &= mask_reverse;

      xored = forward ^ reverse;
//      std::cout << "kmer:\t" << std::bitset<word_size>(forward) << std::endl << "\t" << std::bitset<word_size>(reverse) << std::endl;
//      std::cout << "\t" << std::bitset<word_size>(xored) << std::endl;
      ++iter;
      return 1;
    }

    TO operator()() {
      return xored;
    }

//    TO operator()(CharType c) {
//      forward >>= nBits;
//      forward |= (static_cast<TO>(ALPHABET::FROM_ASCII[c]) << shift);
//      return forward;
//    }


};




int main(int argc, char* argv[]) {
  LOG_INIT();


  std::string filename("/home/tpan/src/bliss/test/test.fastq");

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
  }

#ifdef USE_MPI
  // broadcast to all
  MPI_Bcast(&file_size, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
#endif

  if (rank == nprocs - 1)
    fprintf(stderr, "file size is %ld\n", file_size);

  /////////////// now try to open the file


  // first generate an approximate partition.
  bliss::io::file_loader::range_type r = bliss::io::file_loader::range_type::block_partition(file_size, nprocs, rank);
  std::cout << rank << " equipart: " << r << std::endl;
  std::cout << rank << " test block aligned: " << r.align_to_page(sysconf(_SC_PAGE_SIZE)) << std::endl;

  // now open the file and begin reading
  std::chrono::high_resolution_clock::time_point t1, t2;
  std::chrono::duration<double> time_span1;
  std::chrono::duration<double> time_span2;
  std::chrono::duration<double> time_span3;
  std::chrono::duration<double> time_span;


  {
    t1 = std::chrono::high_resolution_clock::now();
    // now open the file
    bliss::io::fastq_loader loader(filename, r, file_size);
    t2 = std::chrono::high_resolution_clock::now();
    time_span3 = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    INFO("MMap " << "elapsed time: " << time_span3.count() << "s.");

    std::cout << rank << " record adjusted " << loader.getRange() << std::endl;


    t1 = std::chrono::high_resolution_clock::now();

    bliss::io::fastq_loader::iterator fastq_start = loader.begin();
    bliss::io::fastq_loader::iterator fastq_end = loader.end();

    typedef generate_kmers<DNA, char*, uint64_t, 21>  op_type;
    op_type kmer_op;
    typedef bliss::iterator::buffered_transform_iterator<op_type, char*> read_iter_type;

    // transform and generate kmers
    std::vector<uint64_t> kmers;
    for (; fastq_start != fastq_end; ++fastq_start) {
      read_iter_type start((*fastq_start).seq, kmer_op);
      read_iter_type end((*fastq_start).seq_end, kmer_op);

      int i = -1;
      // NOTE: need to call *start to actually evaluate.  question is whether ++ should be doing computation.
      for (; start != end; ++start) {
        ++i;
        if (i < 20) continue;
        kmers.push_back(*start);
      }

    }
    t2 = std::chrono::high_resolution_clock::now();
    time_span3 = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    INFO("kmer generation " << "elapsed time: " << time_span3.count() << "s.");

    int i = 0;
    for (auto kmer : kmers) {
      std::cout << i << " " << std::bitset<64>(kmer) << std::endl;
      ++i;
    }


  }


  {
    t1 = std::chrono::high_resolution_clock::now();
    // now open the file
    bliss::io::fastq_loader loader(filename, r, file_size, true);
    t2 = std::chrono::high_resolution_clock::now();
    time_span3 = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    INFO("MMap preload " << "elapsed time: " << time_span3.count() << "s.");

    std::cout << rank << " record adjusted preload " << loader.getRange() << std::endl;


    t1 = std::chrono::high_resolution_clock::now();

    bliss::io::fastq_loader::iterator fastq_start = loader.begin();
    bliss::io::fastq_loader::iterator fastq_end = loader.end();

    typedef generate_kmers<DNA, char*, uint64_t, 21>  op_type;
    op_type kmer_op;
    typedef bliss::iterator::buffered_transform_iterator<op_type, char*> read_iter_type;

    // transform and generate kmers
    std::vector<uint64_t> kmers;
    for (; fastq_start != fastq_end; ++fastq_start) {
      read_iter_type start((*fastq_start).seq, kmer_op);
      read_iter_type end((*fastq_start).seq_end, kmer_op);

      int i = -1;
      for (; start != end; ++start) {
        ++i;
        if (i < 20) continue;
        kmers.push_back(*start);
      }

    }
    t2 = std::chrono::high_resolution_clock::now();
    time_span3 = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    INFO("kmer generation preload " << "elapsed time: " << time_span3.count() << "s.");


  }
#ifdef USE_MPI
  MPI_Finalize();

#endif


  return 0;


}
