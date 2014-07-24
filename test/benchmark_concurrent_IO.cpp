//#define USE_OPENMP 1

#include "config.hpp"

#include <iostream> // cout
#include <iomanip>  // for setprecision

#ifdef USE_MPI
#include "mpi.h"
#endif

#include <cmath>      // log
#include <chrono>     // timing
#include <algorithm>  // for std::min
#include <fcntl.h>      // for open
#include <sys/stat.h>   // block size.
#include <unistd.h>     // sysconf
#include <sys/mman.h>   // mmap
#include <cstring>      // memcpy, strerror

#include "omp_patterns.hpp"
#include "io/file_loader.hpp"
#include "io/fastq_loader.hpp"
#include "common/alphabets.hpp"
#include "io/fastq_iterator.hpp"

template <typename OT, bool buffering = false, bool preloading = false>
struct readMMap {

    RangeType r;
    unsigned char* data;
    unsigned char* mapped_data;

    int file_handle;
    size_t page_size;

    size_t start;
    size_t chunkSize;

    readMMap(std::string filename, int nprocs, int rank, int nthreads, int _chunkSize) : data(nullptr), mapped_data(nullptr), chunkSize(_chunkSize) {

      /// get file size.
      struct stat filestat;
      int ret = stat(filename.c_str(), &filestat);
      if (ret < 0 ) {
        std::cerr << "ERROR in file open to get size." << std::endl;
        exit(-1);
      }
      size_t file_size = static_cast<size_t>(filestat.st_size);


      /// open the file and get a handle.
      file_handle = open(filename.c_str(), O_RDONLY);
      if (file_handle == -1)
      {
        int myerr = errno;
        std::cerr << "ERROR in file open: ["  << filename << "] error " << myerr << ": " << strerror(myerr);
        exit(-1);
      }

      /// get the block size
      page_size = sysconf(_SC_PAGE_SIZE);




      bliss::partition::BlockPartitioner<bliss::partition::range<size_t> > part;
      RangeType _r(0, file_size);
      part.configure(_r, nprocs);


      // mmap
      r = part.getNext(rank);
      r = r.align_to_page(page_size);
      mapped_data = (unsigned char*)mmap(nullptr, r.end - r.block_start,
                                           PROT_READ,
                                           MAP_PRIVATE, file_handle,
                                           r.block_start);
      data = mapped_data + (r.start - r.block_start);

      if (preloading)
      {
        data = new unsigned char[r.size()];
        memcpy(data, mapped_data + (r.start - r.block_start), r.size());

        munmap(mapped_data, r.end - r.block_start);
        mapped_data = data;
      }

      start = r.start;
    }

    ~readMMap() {
      // unmap
      if (preloading) {
        delete [] data;
      } else {
        munmap(mapped_data, r.end - r.block_start);
      }

      // close the file handle
      if (file_handle != -1) {
        close(file_handle);
        file_handle = -1;
      }
    }

    RangeType getRange() {
      return r;
    }
    size_t getChunkSize() {
      return chunkSize;
    }

    std::string getName() {
      return "readMMap";
    }

    void reset() {
      start = r.start;
    }


    bool operator()(int tid, size_t &count, OT &v) {
      size_t s;

#pragma omp atomic capture
      {
        s = start;
        start += chunkSize;
      }

      RangeType r1(s, s+chunkSize);
      RangeType r2(s, s+2*chunkSize);
      r1 &= r;
      r2 &= r;

      if (r2.size() == 0) {
        //printf("%d empty at %lu - %lu, in range %lu - %lu! chunkSize %lu\n", tid, r2.start, r2.end, r.start, r.end, chunkSize);
        return true;
      }

      unsigned char * ld = data  + (r1.start - r.start);
      if (buffering) {
        ld = new unsigned char[r2.size()];   // TODO: can preallocate.
        memcpy(ld, data + (r1.start - r.start), r2.size());
      }

      // simulate multiple traversals.
      unsigned char c = 0;
      for (size_t i = 0; i < r1.size(); ++i) {
        c = std::max(ld[i], c);
      }

      unsigned char d = 255;
      for (size_t i = 0; i < r1.size(); ++i) {
        d = std::min(ld[i], d);
      }

      // simulate kmer computation
      uint64_t km = c;
      km += d;
      size_t lcount = 0;
      for (size_t i = 0 ; i < r1.size(); ++i, ++lcount) {
        km <<= 8;
        km |= static_cast<uint64_t>(ld[i]);
      }

      // simulate quality score computation.
      OT tv = static_cast<OT>(km) / static_cast<OT>(std::numeric_limits<OT>::max() );
      for (size_t i = 0 ; i < r1.size(); ++i) {
        tv += log2(ld[i]);
      }

      if (buffering)
        delete [] ld;

      count += lcount;

      v += tv;
      return false;
    }
};


template <typename OT, bool buffering = false, bool preloading = false>
struct readFileLoader {

    typedef bliss::io::FileLoader<unsigned char, buffering, preloading> LoaderType;

    size_t page_size;
    typename LoaderType::DataType data;

    LoaderType loader;

    readFileLoader(std::string filename, int nprocs, int rank, int nthreads, int chunkSize) :
      loader(filename, nprocs, rank, nthreads, chunkSize)  {

      RangeType r = loader.getNextPartitionRange(rank);

      loader.load(r);

      /// get the block size
      page_size = sysconf(_SC_PAGE_SIZE);

//      printf("loader from thread %d range = %lu %lu\n", omp_get_thread_num(), _r.start, _r.end);

      data = loader.getData();
    }

    ~readFileLoader() {
      // unmap
      loader.unload();
    }

    RangeType getRange() {
      return data.getRange();
    }
    size_t getChunkSize() {
      return loader.getChunkSize();
    }


    std::string getName() {
      return "readFileLoader";
    }

    void reset() {
      loader.resetChunkRange();
    }

    bool operator()(int tid, size_t &count, OT &v) {
      RangeType r = data.getRange();
      RangeType r1 = loader.getNextChunkRange(tid);
      RangeType r2(r1.start, r1.end + loader.getChunkSize());

      r1 &= r;
      r2 &= r;
      // try copying the data.

      if (r2.size() == 0)
        return true;

      auto ld = data.begin() + (r1.start - r.start);

      // simulate multiple traversals.
      unsigned char c = 0;
      for (size_t i = 0; i < r1.size(); ++i) {
        c = std::max(ld[i], c);
      }

      unsigned char d = 255;
      for (size_t i = 0; i < r1.size(); ++i) {
        d = std::min(ld[i], d);
      }

      size_t lcount = 0;
      // simulate kmer computation
      uint64_t km = c;
      km += d;
      for (size_t i = 0 ; i < r1.size(); ++i, ++lcount) {
        km <<= 8;
        km |= static_cast<uint64_t>(ld[i]);
      }

      // simulate quality score computation.
      OT tv = static_cast<OT>(km) / static_cast<OT>(std::numeric_limits<OT>::max() );
      for (size_t i = 0 ; i < r1.size(); ++i) {
        tv += log2(ld[i]);
      }

      count += lcount;
      v += tv;
      return false;
    }
};


template <typename OT, bool buffering = false, bool preloading = false>
struct readFileLoaderAtomic {

    typedef bliss::io::FileLoader<unsigned char, buffering, preloading> LoaderType;


    size_t page_size;

    LoaderType loader;

    readFileLoaderAtomic(std::string filename, int nprocs, int rank, int nthreads, int chunkSize) :
      loader(filename, nprocs, rank, nthreads, chunkSize)  {

      RangeType r = loader.getNextPartitionRange(rank);

      loader.load(r);

      /// get the block size
      page_size = sysconf(_SC_PAGE_SIZE);

      // mmap
//      printf("loader from thread %d range = %lu %lu\n", omp_get_thread_num(), _r.start, _r.end);

    }

    ~readFileLoaderAtomic() {
      // unmap
      loader.unload();
    }

    RangeType getRange() {
      return loader.getData().getRange();
    }
    size_t getChunkSize() {
      return loader.getChunkSize();
    }

    std::string getName() {
      return "readFileLoaderAtomic";
    }

    void reset() {
      loader.resetChunkRange();
    }

    bool operator()(int tid, size_t &count, OT &v) {

      // try copying the data.
      RangeType r = loader.getNextChunkRange(tid);
      typename LoaderType::DataBlockType data = loader.getChunk(tid, r);

      if (data.begin() == data.end())
        return true;

      // simulate multiple traversals.
      unsigned char c = 0;
      for (auto iter = data.begin(); iter != data.end(); ++iter) {
        c = std::max(*iter, c);
      }

      unsigned char d = 255;

      for (auto iter = data.begin(); iter != data.end(); ++iter) {
        d = std::min(*iter, d);
      }
      size_t lcount = 0;

      // simulate kmer computation
      uint64_t km = c;
      km += d;
      for (auto iter = data.begin(); iter != data.end(); ++iter, ++lcount) {
        km <<= 8;
        km |= static_cast<uint64_t>(*iter);
      }


      // simulate quality score computation.
      OT tv = static_cast<OT>(km) / static_cast<OT>(std::numeric_limits<OT>::max() );
      for (auto iter = data.begin(); iter != data.end(); ++iter) {
        tv += log2(*iter);
      }

      count += lcount;
      v+= tv;
      return false;
    }
};



template <typename OT, bool buffering = false, bool preloading = false>
struct readFASTQ {
    typedef bliss::io::FASTQLoader<unsigned char, buffering, preloading> LoaderType;


    size_t page_size;

    LoaderType loader;

    readFASTQ(std::string filename, int nprocs, int rank, int nthreads, int chunkSize) :
      loader(filename, nprocs, rank, nthreads, chunkSize)  {

      RangeType r = loader.getNextPartitionRange(rank);

      loader.load(r);

//      printf("reading: %lu - %lu, values '%c' '%c'\n", r.start, r.end, *(loader.getData().begin()), *(loader.getData().end()));



      /// get the block size
      page_size = sysconf(_SC_PAGE_SIZE);

      // mmap
//      printf("loader from thread %d range = %lu %lu\n", omp_get_thread_num(), _r.start, _r.end);

    }

    ~readFASTQ() {
      // unmap
      loader.unload();
    }

    RangeType getRange() {
      return loader.getData().getRange();
    }
    size_t getChunkSize() {
      return loader.getChunkSize();
    }

    std::string getName() {
      return "readFASTQ";
    }

    void reset() {
      loader.resetChunkRange();
    }

    bool operator()(int tid, size_t &count, OT &v) {

      // try copying the data.
      RangeType r = loader.getNextChunkRange(tid);
      typename LoaderType::DataBlockType data = loader.getChunk(tid, r);

      if (data.begin() == data.end() )
        return true;

      // simulate multiple traversals.
      unsigned char c = 0;
      for (auto iter = data.begin(); iter != data.end(); ++iter) {
        c = std::max(*iter, c);
      }

      unsigned char d = 255;

      for (auto iter = data.begin(); iter != data.end(); ++iter) {
        d = std::min(*iter, d);
      }

      // simulate kmer computation
      size_t lcount = 0;
      uint64_t km = c;
      km += d;
      for (auto iter = data.begin(); iter != data.end(); ++iter, ++lcount) {
        km <<= 8;
        km |= static_cast<uint64_t>(*iter);
      }

      // simulate quality score computation.
      OT tv = static_cast<OT>(km) / static_cast<OT>(std::numeric_limits<uint64_t>::max() );
      for (auto iter = data.begin(); iter != data.end(); ++iter) {
        tv += log2(*iter);
      }

      count += lcount;
      v += tv;
      return false;
    }
};

template <typename OT, bool buffering = false, bool preloading = false>
struct FASTQIterator {
    typedef bliss::io::FASTQLoader<unsigned char, buffering, preloading> LoaderType;

    size_t page_size;


    LoaderType loader;

    typedef float QualityType;
    typedef DNA Alphabet;
    typedef typename LoaderType::BlockIteratorType BaseIterType;
    typedef bliss::io::fastq_sequence_quality<BaseIterType, Alphabet, QualityType>  SequenceType;

    typedef bliss::io::fastq_parser<BaseIterType, Alphabet, QualityType>  ParserType;
    typedef bliss::io::fastq_iterator<ParserType, BaseIterType>           IteratorType;

    FASTQIterator(std::string filename, int nprocs, int rank, int nthreads, int chunkSize) :
      loader(filename, nprocs, rank, nthreads, chunkSize)  {

      RangeType r = loader.getNextPartitionRange(rank);
      loader.load(r);

      /// get the block size
      page_size = sysconf(_SC_PAGE_SIZE);

      // mmap
//      printf("loader from thread %d range = %lu %lu\n", omp_get_thread_num(), _r.start, _r.end);

    }

    ~FASTQIterator() {
      // unmap
      loader.unload();
    }

    RangeType getRange() {
      return loader.getData().getRange();
    }
    size_t getChunkSize() {
      return loader.getChunkSize();
    }

    std::string getName() {
      return "FASTQIterator";
    }


    void reset() {
      loader.resetChunkRange();
    }

    bool operator()(int tid, size_t &count, OT &v) {

      // try copying the data.
      RangeType r = loader.getNextChunkRange(tid);
      typename LoaderType::DataBlockType data = loader.getChunk(tid, r);

      if (data.begin() ==  data.end())
        return true;

      // traverse using fastq iterator.
      ParserType parser;
      IteratorType fastq_start(parser, data.begin(), data.end(), data.getRange());
      IteratorType fastq_end(parser, data.end(), data.getRange());

      SequenceType read;
      unsigned char c = 0;
      unsigned char d = 255;
      uint64_t km = 0;
      OT tv = 0;
      size_t lcount = 0;

      for (; fastq_start != fastq_end; ++fastq_start, ++lcount)
      {
        read = *fastq_start;

        // now simulate the compute

        // simulate kmer computation
        // simulate multiple traversals.
        for (BaseIterType iter = read.seq; iter != read.seq_end; ++iter) {
          c = std::max(*iter, c);
        }
        for (BaseIterType iter = read.qual; iter != read.qual_end; ++iter) {
          c = std::max(*iter, c);
        }

        for (BaseIterType iter = read.seq; iter != read.seq_end; ++iter) {
          d = std::min(*iter, d);
        }
        for (BaseIterType iter = read.qual; iter != read.qual_end; ++iter) {
          d = std::min(*iter, d);
        }

        // simulate kmer computation
        for (BaseIterType iter = read.seq; iter != read.seq_end; ++iter) {
          km <<= 8;
          km |= static_cast<uint64_t>(*iter);
        }
        for (BaseIterType iter = read.qual; iter != read.qual_end; ++iter) {
          km <<= 8;
          km |= static_cast<uint64_t>(*iter);
        }

        // simulate quality score computation.
        tv += static_cast<OT>(km) / static_cast<OT>(std::numeric_limits<uint64_t>::max() );
        for (BaseIterType iter = read.seq; iter != read.seq_end; ++iter) {
          tv += log2(*iter);
        }
        for (BaseIterType iter = read.qual; iter != read.qual_end; ++iter) {
          tv += log2(*iter);
        }

      }

      count += lcount;
      v += tv;
      return false;
    }
};


template <typename OT, bool buffering = false, bool preloading = false>
struct FASTQIterator2 {
    typedef bliss::io::FASTQLoader<unsigned char, buffering, preloading> LoaderType;

    size_t page_size;


    LoaderType loader;

    typedef float QualityType;
    typedef DNA Alphabet;
    typedef typename LoaderType::BlockIteratorType BaseIterType;
    typedef bliss::io::fastq_sequence_quality<BaseIterType, Alphabet, QualityType>  SequenceType;

    typedef bliss::io::fastq_parser<BaseIterType, Alphabet, QualityType>  ParserType;
    typedef bliss::io::fastq_iterator<ParserType, BaseIterType>           IteratorType;

    FASTQIterator2(std::string filename, int nprocs, int rank, int nthreads, int chunkSize) :
      loader(filename, nprocs, rank, nthreads, chunkSize)  {

      RangeType r = loader.getNextPartitionRange(rank);
      loader.load(r);
//      printf("loader from thread %d loaded range = %lu %lu\n", omp_get_thread_num(), loader.getData().getRange().start, loader.getData().getRange().end);
//      printf("loader from thread %d file range = %lu %lu\n", omp_get_thread_num(), loader.getFileRange().start, loader.getFileRange().end);

      /// get the block size
      page_size = sysconf(_SC_PAGE_SIZE);

    }

    ~FASTQIterator2() {
      // unmap
      loader.unload();
    }

    RangeType getRange() {
      return loader.getData().getRange();
    }
    size_t getChunkSize() {
      return loader.getChunkSize();
    }

    std::string getName() {
      return "FASTQIterator2";
    }

    void reset() {
      loader.resetChunkRange();
    }

    bool operator()(int tid, size_t &count, OT &v) {

      // try copying the data.
      RangeType r = loader.getNextChunkRange(tid);
      typename LoaderType::DataBlockType data = loader.getChunk(tid, r);

      if (data.begin() ==  data.end())
        return true;

      // traverse using fastq iterator.
      ParserType parser;
      IteratorType fastq_start(parser, data.begin(), data.end(), data.getRange());
      IteratorType fastq_end(parser, data.end(), data.getRange());

      SequenceType read;
      uint64_t km = 0;
      OT tv = 0;
      size_t lcount = 0;

//      int i = 0, j = 0;

      BaseIterType iter;

      for (; fastq_start != fastq_end; ++fastq_start, ++lcount)
      {
        read = *fastq_start;

        // now simulate the compute

        // simulate kmer computation

        // simulate kmer computation
        iter = read.seq;
        for (;
            iter != read.seq_end;
            ++iter) {
          km <<= 8;
          km |= static_cast<uint64_t>(*iter);
//          ++i;
        }

        // simulate quality score computation.
        tv += static_cast<OT>(km) / static_cast<OT>(std::numeric_limits<uint64_t>::max() );
        iter = read.qual;
        for (;   // slow because of log2.
            iter != read.qual_end;
            ++iter) {
          tv += log2(*iter);
//          ++j;
        }

//        assert(i == j);

      }

      count += lcount;
      v += tv;
      return false;
    }
};


template <typename OT, bool buffering = false, bool preloading = false>
struct FASTQIteratorNoQual {
    typedef bliss::io::FASTQLoader<unsigned char, buffering, preloading> LoaderType;

    size_t page_size;


    LoaderType loader;

    typedef float QualityType;
    typedef DNA Alphabet;
    typedef typename LoaderType::BlockIteratorType BaseIterType;
    typedef bliss::io::fastq_sequence_quality<BaseIterType, Alphabet, QualityType>  SequenceType;

    typedef bliss::io::fastq_parser<BaseIterType, Alphabet, QualityType>  ParserType;
    typedef bliss::io::fastq_iterator<ParserType, BaseIterType>           IteratorType;

    FASTQIteratorNoQual(std::string filename, int nprocs, int rank, int nthreads, int chunkSize) :
      loader(filename, nprocs, rank, nthreads, chunkSize)  {

      RangeType r = loader.getNextPartitionRange(rank);
      loader.load(r);

      /// get the block size
      page_size = sysconf(_SC_PAGE_SIZE);

      // mmap
//      printf("loader from thread %d range = %lu %lu\n", omp_get_thread_num(), _r.start, _r.end);

    }

    ~FASTQIteratorNoQual() {
      // unmap
      loader.unload();
    }

    RangeType getRange() {
      return loader.getData().getRange();
    }
    size_t getChunkSize() {
      return loader.getChunkSize();
    }

    std::string getName() {
      return "FASTQIteratorNoQual";
    }
    void reset() {
      loader.resetChunkRange();
    }

    bool operator()(int tid, size_t &count, OT &v) {

      // try copying the data.
      RangeType r = loader.getNextChunkRange(tid);
      typename LoaderType::DataBlockType data = loader.getChunk(tid, r);
//      printf("thread %d getting block %lu-%lu, got block of length %ld\n", omp_get_thread_num(), start, end, (data.end() - data.begin()));


      if (data.begin() == data.end()) {
//        std::cerr << " range = " << start << "-" << end << std::endl;
        return true;
      }

//      std::cout << "Range: " << data.getRange() << std::endl;
//      std::cout << "Start: " << data.begin()[0] << std::endl;
//      std::cout << "End: "   << data.end()[0] << std::endl;
//      std::cout << "len: "   << (data.end() - data.begin())  << " len from range: " << (data.getRange().end - data.getRange().start) << std::endl;


      // traverse using fastq iterator.
      ParserType parser;
      IteratorType fastq_start(parser, data.begin(), data.end(), data.getRange());
      IteratorType fastq_end(parser, data.end(), data.getRange());

      SequenceType read;
      uint64_t km = 0;
      OT tv = 0;
      size_t lcount = 0;

      for (; fastq_start != fastq_end; ++fastq_start, ++lcount)
      {
        read = *fastq_start;

        // now simulate the compute

        // simulate kmer computation

        // simulate kmer computation
        for (BaseIterType iter = read.seq; iter != read.seq_end; ++iter) {
          km <<= 8;
          km |= static_cast<uint64_t>(*iter);
        }

        // simulate quality score computation.
        tv += static_cast<OT>(km) / static_cast<OT>(std::numeric_limits<uint64_t>::max() );

      }
      count += lcount;
      v += tv;
      return false;
    }
};



void printTiming(std::string tag, std::string name, int rank, int nprocs, int nthreads,
                 const std::chrono::duration<double>& time_span, int iter,
                 double v, size_t count)
{
  std::cout << name << "\t" << tag <<"\tMPI rank: " << rank << "/" << nprocs << "\tOMP "
            << nthreads << " threads\ttook " << std::fixed
            << std::setprecision(6) << time_span.count() / iter
            << "s,\tresult = " << v << " count = " << count << std::endl;
}

int main(int argc, char* argv[])
{

  int rank = 0, nprocs = 1;
#ifdef USE_MPI

  // initialize MPI
  MPI_Init(&argc, &argv);

  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == 0)
  std::cout << "USE_MPI is set" << std::endl;
#endif

  int nthreads = 1;
#ifdef USE_OPENMP
  if (rank == 0)
  std::cout << "USE_OPENMP is set" << std::endl;
  omp_set_nested(1);
  omp_set_dynamic(0);
  nthreads = omp_get_max_threads();
#endif


  if (argc > 1)
    nthreads = atoi(argv[1]);

  size_t step = 4096;
  if (argc > 2)
    step = atoi(argv[2]);

  /// set up the file.
  std::string filename(PROJ_SRC_DIR);
  filename.append("/test/data/test.fastq");
  if (argc > 3)
  {
    filename.assign(argv[3]);
  }

  int iter = 10;
  if (argc > 4)
    iter = atoi(argv[4]);




#if defined(TEST_OP_MMAP)
  typedef readMMap<             double, true, false> OpType;
#endif

#if defined(TEST_OP_FILELOADER)
  typedef readFileLoader<       double, true, false> OpType;
#endif

#if defined(TEST_OP_FILELOADER_ATOMIC)
  typedef readFileLoaderAtomic< double, true, false> OpType;
#endif

#if defined(TEST_OP_FASTQ)
  typedef readFASTQ<            double, true, false> OpType;
#endif

#if defined(TEST_OP_FASTQIter)
  typedef FASTQIterator<        double, true, false> OpType;
#endif

#if defined(TEST_OP_FASTQIter2)
  typedef FASTQIterator2<       double, true, false> OpType;
#endif

#if defined(TEST_OP_FASTQIterNoQual)
  typedef FASTQIteratorNoQual<  double, true, false> OpType;
#endif

  OpType op(filename, nprocs, rank, nthreads, step);

  double v = 0.0;
  size_t count = 0;

  std::chrono::high_resolution_clock::time_point t1, t2;
  std::chrono::duration<double> time_span;

  /// Workers only, critical
#if defined(USE_MPI)
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  t1 = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < iter; ++i) {
    op.reset();
    count = 0;
    v = P2P<OpType, double>(op, nthreads, count);
  }
#if defined(USE_MPI)
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  t2 = std::chrono::high_resolution_clock::now();
  time_span =
      std::chrono::duration_cast<std::chrono::duration<double>>(
          t2 - t1);
  printTiming("P2P critical:", op.getName(), rank, nprocs, nthreads, time_span, iter, v, count);


  /// master slave
#if defined(USE_MPI)
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  t1 = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < iter; ++i) {
    op.reset();
    count = 0;
    v= MasterSlave<OpType, double>(op, nthreads, count);
  }
#if defined(USE_MPI)
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  t2 = std::chrono::high_resolution_clock::now();
  time_span =
      std::chrono::duration_cast<std::chrono::duration<double>>(
          t2 - t1);
  printTiming("MS Wait:", op.getName(), rank, nprocs, nthreads, time_span, iter, v, count);

  /// master slave No Wait
#if defined(USE_MPI)
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  t1 = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < iter; ++i) {
    op.reset();
    count = 0;
    v= MasterSlaveNoWait<OpType, double>(op, nthreads, count);
  }
#if defined(USE_MPI)
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  t2 = std::chrono::high_resolution_clock::now();
  time_span =
      std::chrono::duration_cast<std::chrono::duration<double>>(
          t2 - t1);
  printTiming("MS NoWait:", op.getName(), rank, nprocs, nthreads, time_span, iter, v, count);

  /// parallel for
#if defined(USE_MPI)
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  t1 = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < iter; ++i) {
    op.reset();
    count = 0;
    v = ParFor<OpType, double>(op, nthreads, count);
  }
#if defined(USE_MPI)
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  t2 = std::chrono::high_resolution_clock::now();
  time_span =
      std::chrono::duration_cast<std::chrono::duration<double>>(
          t2 - t1);
  printTiming("PARFOR:\t", op.getName(), rank, nprocs, nthreads, time_span, iter, v, count);


  //// block parallel  for
#if defined(USE_MPI)
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  t1 = std::chrono::high_resolution_clock::now();
  count = 0;
  for (int i = 0; i < iter; ++i) {
    v = 0;
    count = 0;
#pragma omp parallel default(none) shared(nthreads, step, filename) num_threads(nthreads) reduction(+:v, count)
    {
      OpType op2(filename, nthreads, omp_get_thread_num(), 1, step);
      op2.reset();
      double v0 = Sequential<OpType, double>(op2, 1, count);
      //printf("%d processing range %lu %lu. result = %f\n", omp_get_thread_num(), op2.getRange().start, op2.getRange().end, v0);

      v += v0;
    }
  }
#if defined(USE_MPI)
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  t2 = std::chrono::high_resolution_clock::now();
  time_span =
      std::chrono::duration_cast<std::chrono::duration<double>>(
          t2 - t1);
  printTiming("BLOCK PARFOR:", op.getName(), rank, nprocs, nthreads, time_span, iter, v, count);



  //// serial for
#if defined(USE_MPI)
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  t1 = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < iter; ++i) {
    op.reset();
    count = 0;
    v = Sequential<OpType, double>(op, 1, count);
  }
#if defined(USE_MPI)
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  t2 = std::chrono::high_resolution_clock::now();
  time_span =
      std::chrono::duration_cast<std::chrono::duration<double>>(
          t2 - t1);
  printTiming("SEQFOR:\t", op.getName(), rank, nprocs, nthreads, time_span, iter, v, count);




#ifdef USE_MPI
  MPI_Finalize();

#endif

  return 0;

}

