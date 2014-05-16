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
#include "io/fastq_partition_helper.hpp"
#include "common/alphabets.hpp"
#include "io/fastq_iterator.hpp"

template <typename OT, bool buffering = false, bool preloading = false>
struct readMMap {

    RangeType r;
    unsigned char* data;
    unsigned char* mapped_data;

    int file_handle;
    size_t page_size;

    readMMap(std::string filename, RangeType _r) {
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

      // mmap
      r = _r.align_to_page(page_size);
      mapped_data = (unsigned char*)mmap(nullptr, r.end - r.block_start,
                                           PROT_READ,
                                           MAP_PRIVATE, file_handle,
                                           r.block_start);
      data = mapped_data + r.start - r.block_start;

      if (preloading)
      {
        data = new unsigned char[r.end - r.start];
        memcpy(data, mapped_data + r.start - r.block_start, r.end - r.start);

        munmap(mapped_data, r.end - r.block_start);
        mapped_data = data;
      }
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

    size_t getSeqSize() {
      return 1;
    }

    RangeType getRange() {
      return r;
    }

    std::string getName() {
      return "readMMap";
    }

    void reset() {
    }


    OT operator()(const size_t start, const size_t end, size_t &count) {
      size_t s = std::max(r.start, start);
      size_t e = std::min(r.end, end);
      size_t e2 = std::min(r.end, end + (end-start));
      // try copying the data.

      unsigned char * ld = data  + s - r.start;
      if (buffering) {
        ld = new unsigned char[e2-s];   // TODO: can preallocate.
        memcpy(ld, data + s - r.start, e2-s);
      }

      // simulate multiple traversals.
      unsigned char c = 0;
      for (size_t i = 0; i < (e2 - s); ++i) {
        c = std::max(ld[i], c);
      }

      unsigned char d = 255;
      for (size_t i = 0; i < (e2 - s); ++i) {
        d = std::min(ld[i], d);
      }

      // simulate kmer computation
      uint64_t km = c;
      km += d;
      size_t lcount = 0;
      for (size_t i = 0 ; i < e - s; ++i, ++lcount) {
        km <<= 8;
        km |= static_cast<uint64_t>(ld[i]);
      }

      // simulate quality score computation.
      OT tv = static_cast<OT>(km) / static_cast<OT>(std::numeric_limits<uint64_t>::max() );
      for (size_t i = e-s ; i < (e2 - s); ++i) {
        tv += log2(ld[i]);
      }

      if (buffering)
        delete [] ld;

      count += lcount;

      return tv;
    }
};


template <typename OT, bool buffering = false, bool preloading = false>
struct readFileLoader {

    typedef bliss::io::file_loader<unsigned char, buffering, preloading> LoaderType;

    size_t page_size;
    typename LoaderType::DataType data;

    LoaderType loader;

    readFileLoader(std::string filename, RangeType _r) :
      loader(filename)  {

      loader.setRange(_r);

      loader.load();

      /// get the block size
      page_size = sysconf(_SC_PAGE_SIZE);

//      printf("loader from thread %d range = %lu %lu\n", omp_get_thread_num(), _r.start, _r.end);

      data = loader.getData();
    }

    ~readFileLoader() {
      // unmap
      loader.unload();
    }

    size_t getSeqSize() {
      return 1;
    }

    RangeType getRange() {
      return data.getRange();
    }

    std::string getName() {
      return "readFileLoader";
    }

    void reset() {
      loader.resetChunks();
    }

    OT operator()(const size_t start, const size_t end, size_t &count) {
      size_t s = std::max(data.getRange().start, start);
      size_t e = std::min(data.getRange().end, end);
      size_t e2 = std::min(data.getRange().end, end + (end-start));
      // try copying the data.

      auto ld = data.begin() + (s - data.getRange().start);

      // simulate multiple traversals.
      unsigned char c = 0;
      for (size_t i = 0; i < (e2 - s); ++i) {
        c = std::max(ld[i], c);
      }

      unsigned char d = 255;
      for (size_t i = 0; i < (e2 - s); ++i) {
        d = std::min(ld[i], d);
      }

      // simulate kmer computation
      uint64_t km = c;
      km += d;
      for (size_t i = 0 ; i < e - s; ++i) {
        km <<= 8;
        km |= static_cast<uint64_t>(ld[i]);
      }

      size_t lcount = 0;
      // simulate quality score computation.
      OT tv = static_cast<OT>(km) / static_cast<OT>(std::numeric_limits<uint64_t>::max() );
      for (size_t i = e-s ; i < (e2 - s); ++i, ++lcount) {
        tv += log2(ld[i]);
      }

      count += lcount;
      return tv;
    }
};

template<typename Iterator, typename Range>
struct PartitionHelper {
    typedef typename Range::ValueType SizeType;
    typedef typename std::iterator_traits<Iterator>::value_type ValueType;
    typedef Iterator  IteratorType;

    const SizeType operator()(const Iterator &iter, const Range & parent, const Range &target) const {
      return target.start;
    }
};


template <typename OT, bool buffering = false, bool preloading = false>
struct readFileLoaderAtomic {

    typedef bliss::io::file_loader<unsigned char, buffering, preloading> LoaderType;


    size_t page_size;

    LoaderType loader;
    PartitionHelper<typename LoaderType::InputIteratorType, RangeType> helper;
    PartitionHelper<typename LoaderType::IteratorType, RangeType> chunkHelper;

    readFileLoaderAtomic(std::string filename, RangeType _r) :
      loader(filename)  {
      loader.setRange(_r);


      loader.adjustRange(helper);

      loader.load();

      /// get the block size
      page_size = sysconf(_SC_PAGE_SIZE);

      // mmap
//      printf("loader from thread %d range = %lu %lu\n", omp_get_thread_num(), _r.start, _r.end);

    }

    ~readFileLoaderAtomic() {
      // unmap
      loader.unload();
    }
    size_t getSeqSize() {
      return loader.getSeqSize(chunkHelper, 3);
    }

    RangeType getRange() {
      return loader.getData().getRange();
    }

    std::string getName() {
      return "readFileLoaderAtomic";
    }

    void reset() {
      loader.resetChunks();
    }

    OT operator()(const size_t start, const size_t end, size_t &count) {

      // try copying the data.

      typename LoaderType::DataBlockType data = loader.getNextChunkAtomic(chunkHelper, end - start);

      if (data.begin() == data.end())
        return 0;

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
      uint64_t km = c;
      km += d;
      for (auto iter = data.begin(); iter != data.end(); ++iter) {
        km <<= 8;
        km |= static_cast<uint64_t>(*iter);
      }

      size_t lcount = 0;

      // simulate quality score computation.
      OT tv = static_cast<OT>(km) / static_cast<OT>(std::numeric_limits<uint64_t>::max() );
      for (auto iter = data.begin(); iter != data.end(); ++iter, ++lcount) {
        tv += log2(*iter);
      }

      count += lcount;
      return tv;
    }
};



template <typename OT, bool buffering = false, bool preloading = false>
struct readFASTQ {
    typedef bliss::io::file_loader<unsigned char, buffering, preloading> LoaderType;


    size_t page_size;

    LoaderType loader;
    bliss::io::FASTQPartitionHelper<typename LoaderType::InputIteratorType, RangeType> helper;
    bliss::io::FASTQPartitionHelper<typename LoaderType::IteratorType, RangeType> chunkHelper;

    readFASTQ(std::string filename, RangeType _r) :
      loader(filename)  {

      loader.setRange(_r);

      loader.adjustRange(helper);

      loader.load();

      /// get the block size
      page_size = sysconf(_SC_PAGE_SIZE);

      // mmap
//      printf("loader from thread %d range = %lu %lu\n", omp_get_thread_num(), _r.start, _r.end);

    }

    ~readFASTQ() {
      // unmap
      loader.unload();
    }
    size_t getSeqSize() {
      return loader.getSeqSize(chunkHelper, 3);
    }

    RangeType getRange() {
      return loader.getData().getRange();
    }
    std::string getName() {
      return "readFASTQ";
    }

    void reset() {
      loader.resetChunks();
    }

    OT operator()(const size_t start, const size_t end, size_t &count) {

      // try copying the data.
      typename LoaderType::DataBlockType data = loader.getNextChunkAtomic(chunkHelper, end - start);

      if (data.begin() == data.end() )
        return 0;

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
      uint64_t km = c;
      km += d;
      for (auto iter = data.begin(); iter != data.end(); ++iter) {
        km <<= 8;
        km |= static_cast<uint64_t>(*iter);
      }

      size_t lcount = 0;
      // simulate quality score computation.
      OT tv = static_cast<OT>(km) / static_cast<OT>(std::numeric_limits<uint64_t>::max() );
      for (auto iter = data.begin(); iter != data.end(); ++iter, ++lcount) {
        tv += log2(*iter);
      }

      count += lcount;
      return tv;
    }
};

template <typename OT, bool buffering = false, bool preloading = false>
struct FASTQIterator {
    typedef bliss::io::file_loader<unsigned char, buffering, preloading> LoaderType;

    size_t page_size;


    LoaderType loader;
    bliss::io::FASTQPartitionHelper<typename LoaderType::InputIteratorType, RangeType> helper;
    bliss::io::FASTQPartitionHelper<typename LoaderType::IteratorType, RangeType> chunkHelper;

    typedef float QualityType;
    typedef DNA Alphabet;
    typedef typename LoaderType::BlockIteratorType BaseIterType;
    typedef bliss::io::fastq_sequence_quality<BaseIterType, Alphabet, QualityType>  SequenceType;

    typedef bliss::io::fastq_parser<BaseIterType, Alphabet, QualityType>  ParserType;
    typedef bliss::io::fastq_iterator<ParserType, BaseIterType>           IteratorType;

    FASTQIterator(std::string filename, RangeType _r) :
      loader(filename)  {

      loader.setRange(_r);
      loader.adjustRange(helper);
      loader.load();

      /// get the block size
      page_size = sysconf(_SC_PAGE_SIZE);

      // mmap
//      printf("loader from thread %d range = %lu %lu\n", omp_get_thread_num(), _r.start, _r.end);

    }

    ~FASTQIterator() {
      // unmap
      loader.unload();
    }
    size_t getSeqSize() {
      return loader.getSeqSize(chunkHelper, 3);
    }

    RangeType getRange() {
      return loader.getData().getRange();
    }
    std::string getName() {
      return "FASTQIterator";
    }


    void reset() {
      loader.resetChunks();
    }

    OT operator()(const size_t start, const size_t end, size_t &count) {

      // try copying the data.
      typename LoaderType::DataBlockType data = loader.getNextChunkAtomic(chunkHelper, end - start);

      if (data.begin() ==  data.end())
        return 0;

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
      return tv;
    }
};


template <typename OT, bool buffering = false, bool preloading = false>
struct FASTQIterator2 {
    typedef bliss::io::file_loader<unsigned char, buffering, preloading> LoaderType;

    size_t page_size;


    LoaderType loader;
    bliss::io::FASTQPartitionHelper<typename LoaderType::InputIteratorType, RangeType> helper;
    bliss::io::FASTQPartitionHelper<typename LoaderType::IteratorType, RangeType> chunkHelper;

    typedef float QualityType;
    typedef DNA Alphabet;
    typedef typename LoaderType::BlockIteratorType BaseIterType;
    typedef bliss::io::fastq_sequence_quality<BaseIterType, Alphabet, QualityType>  SequenceType;

    typedef bliss::io::fastq_parser<BaseIterType, Alphabet, QualityType>  ParserType;
    typedef bliss::io::fastq_iterator<ParserType, BaseIterType>           IteratorType;

    FASTQIterator2(std::string filename, RangeType _r) :
      loader(filename)  {

//      printf("loader from thread %d initial range = %lu %lu\n", omp_get_thread_num(), loader.getRange().start, loader.getRange().end);

      loader.setRange(_r);
//      printf("loader from thread %d requested range = %lu %lu\n", omp_get_thread_num(), loader.getRange().start, loader.getRange().end);
      loader.adjustRange(helper);

//      printf("loader from thread %d adjusted range = %lu %lu\n", omp_get_thread_num(), loader.getRange().start, loader.getRange().end);

      loader.load();
//      printf("loader from thread %d loaded range = %lu %lu\n", omp_get_thread_num(), loader.getData().getRange().start, loader.getData().getRange().end);
//      printf("loader from thread %d file range = %lu %lu\n", omp_get_thread_num(), loader.getFileRange().start, loader.getFileRange().end);

      /// get the block size
      page_size = sysconf(_SC_PAGE_SIZE);

    }

    ~FASTQIterator2() {
      // unmap
      loader.unload();
    }
    size_t getSeqSize() {
      return loader.getSeqSize(chunkHelper, 3);
    }
    RangeType getRange() {
      return loader.getData().getRange();
    }
    std::string getName() {
      return "FASTQIterator2";
    }

    void reset() {
      loader.resetChunks();
    }

    OT operator()(const size_t start, const size_t end, size_t &count) {

      // try copying the data.
      typename LoaderType::DataBlockType data = loader.getNextChunkAtomic(chunkHelper, end - start);

      if (data.begin() ==  data.end())
        return 0;

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
        for (BaseIterType iter = read.qual; iter != read.qual_end; ++iter) {
          tv += log2(*iter);
        }

      }

      count += lcount;
      return tv;
    }
};


template <typename OT, bool buffering = false, bool preloading = false>
struct FASTQIteratorNoQual {
    typedef bliss::io::file_loader<unsigned char, buffering, preloading> LoaderType;

    size_t page_size;


    LoaderType loader;
    bliss::io::FASTQPartitionHelper<typename LoaderType::InputIteratorType, RangeType> helper;
    bliss::io::FASTQPartitionHelper<typename LoaderType::IteratorType, RangeType> chunkHelper;

    typedef float QualityType;
    typedef DNA Alphabet;
    typedef typename LoaderType::BlockIteratorType BaseIterType;
    typedef bliss::io::fastq_sequence_quality<BaseIterType, Alphabet, QualityType>  SequenceType;

    typedef bliss::io::fastq_parser<BaseIterType, Alphabet, QualityType>  ParserType;
    typedef bliss::io::fastq_iterator<ParserType, BaseIterType>           IteratorType;

    FASTQIteratorNoQual(std::string filename, RangeType _r) :
      loader(filename)  {

      loader.setRange(_r);


      loader.adjustRange(helper);

      loader.load();

      /// get the block size
      page_size = sysconf(_SC_PAGE_SIZE);

      // mmap
//      printf("loader from thread %d range = %lu %lu\n", omp_get_thread_num(), _r.start, _r.end);

    }

    ~FASTQIteratorNoQual() {
      // unmap
      loader.unload();
    }

    size_t getSeqSize() {
      return loader.getSeqSize(chunkHelper, 3);
    }
    RangeType getRange() {
      return loader.getData().getRange();
    }
    std::string getName() {
      return "FASTQIteratorNoQual";
    }
    void reset() {
      loader.resetChunks();
    }

    OT operator()(const size_t start, const size_t end, size_t &count) {


      // try copying the data.
      //printf("thread %d getting block %lu-%lu\n", omp_get_thread_num(), start, end);

      typename LoaderType::DataBlockType data = loader.getNextChunkAtomic(chunkHelper, end - start);

//      printf("thread %d getting block %lu-%lu, got block of length %ld\n", omp_get_thread_num(), start, end, (data.end() - data.begin()));


      if (data.begin() == data.end()) {
//        std::cerr << " range = " << start << "-" << end << std::endl;
        return 0;
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
      return tv;
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

#ifdef USE_OPENMP
  if (rank == 0)
  std::cout << "USE_OPENMP is set" << std::endl;
  omp_set_nested(1);
  omp_set_dynamic(0);
#endif


  int nthreads = 1;
  if (argc > 1)
    nthreads = atoi(argv[1]);

  size_t step = 30;
  if (argc > 2)
    step = atoi(argv[2]);

  /// set up the file.
  std::string filename("/home/tpan/src/bliss/test/data/test.fastq");
  if (argc > 3)
  {
    filename.assign(argv[3]);
  }

  int iter = 10;
  if (argc > 4)
    iter = atoi(argv[4]);

  /// get file size
  size_t file_size = 0;
  if (rank == 0)
  {
    struct stat filestat;
    stat(filename.c_str(), &filestat);
    file_size = static_cast<size_t>(filestat.st_size);
  }

#if defined(USE_MPI)
  if (nprocs > 1) {
    /// broadcast filesize to all
    MPI_Bcast(&file_size, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
  }
#endif

  RangeType r = RangeType::block_partition(nprocs, rank, 0, file_size, 0);


//  typedef readMMap<             double, true , false> OpType;
//  typedef readFileLoader<       double, true , false> OpType;
//  typedef readFileLoaderAtomic< double, true , false> OpType;
//  typedef readFASTQ<            double, true , false> OpType;
//  typedef FASTQIterator<        double, true , false> OpType;
  typedef FASTQIterator2<       double, true , false> OpType;
//  typedef FASTQIteratorNoQual<  double, true , false> OpType;

  OpType op(filename, r);

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
    v = P2P(op, nthreads, op.getRange(), op.getSeqSize() * step, count);
  }
#if defined(USE_MPI)
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  t2 = std::chrono::high_resolution_clock::now();
  time_span =
      std::chrono::duration_cast<std::chrono::duration<double>>(
          t2 - t1);
  if (rank == 0)
  printTiming("P2P critical:", op.getName(), rank, nprocs, nthreads, time_span, iter, v, count);


  /// Workers only, atomic  - close to time for parfor.
#if defined(USE_MPI)
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  t1 = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < iter; ++i) {
    op.reset();
    count = 0;
    v = P2P_atomic(op, nthreads, op.getRange(), op.getSeqSize() * step, count);
  }
#if defined(USE_MPI)
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  t2 = std::chrono::high_resolution_clock::now();
  time_span =
      std::chrono::duration_cast<std::chrono::duration<double>>(
          t2 - t1);
  if (rank == 0)
  printTiming("P2P atomic:", op.getName(), rank, nprocs, nthreads, time_span, iter, v, count);


  /// master slave
#if defined(USE_MPI)
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  t1 = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < iter; ++i) {
    op.reset();
    count = 0;
    v= MasterSlave(op, nthreads, op.getRange(), op.getSeqSize() * step, count);
  }
#if defined(USE_MPI)
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  t2 = std::chrono::high_resolution_clock::now();
  time_span =
      std::chrono::duration_cast<std::chrono::duration<double>>(
          t2 - t1);
  if (rank == 0)
  printTiming("MS Wait:", op.getName(), rank, nprocs, nthreads, time_span, iter, v, count);

  /// master slave No Wait
#if defined(USE_MPI)
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  t1 = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < iter; ++i) {
    op.reset();
    count = 0;
    v= MasterSlaveNoWait(op, nthreads, op.getRange(), op.getSeqSize() * step, count);
  }
#if defined(USE_MPI)
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  t2 = std::chrono::high_resolution_clock::now();
  time_span =
      std::chrono::duration_cast<std::chrono::duration<double>>(
          t2 - t1);
  if (rank == 0)
  printTiming("MS NoWait:", op.getName(), rank, nprocs, nthreads, time_span, iter, v, count);

  /// parallel for
#if defined(USE_MPI)
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  t1 = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < iter; ++i) {
    op.reset();
    count = 0;
    v = ParFor(op, nthreads, op.getRange(), op.getSeqSize() * step, count);
  }
#if defined(USE_MPI)
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  t2 = std::chrono::high_resolution_clock::now();
  time_span =
      std::chrono::duration_cast<std::chrono::duration<double>>(
          t2 - t1);
  if (rank == 0)
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
#pragma omp parallel default(none) firstprivate(nthreads, step, filename, r) num_threads(nthreads) reduction(+:v, count)
    {
      RangeType r2 = r.block_partition(nthreads, omp_get_thread_num());
      OpType op2(filename, r2);
      v = Sequential(op2, nthreads, op2.getRange(), op2.getSeqSize() * step, count);
    }
  }
#if defined(USE_MPI)
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  t2 = std::chrono::high_resolution_clock::now();
  time_span =
      std::chrono::duration_cast<std::chrono::duration<double>>(
          t2 - t1);
  if (rank == 0)
  printTiming("BLOCK PARFOR:", op.getName(), rank, nprocs, nthreads, time_span, iter, v, count);



  //// serial for
#if defined(USE_MPI)
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  t1 = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < iter; ++i) {
    op.reset();
    count = 0;
    v = Sequential(op, 1, op.getRange(), op.getSeqSize() * step, count);
  }
#if defined(USE_MPI)
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  t2 = std::chrono::high_resolution_clock::now();
  time_span =
      std::chrono::duration_cast<std::chrono::duration<double>>(
          t2 - t1);
  if (rank == 0)
  printTiming("SEQFOR:\t", op.getName(), rank, nprocs, nthreads, time_span, iter, v, count);


#ifdef USE_MPI
  MPI_Finalize();

#endif

  return 0;

}

