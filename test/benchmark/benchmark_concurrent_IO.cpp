//#define USE_OPENMP 1

#include "bliss-config.hpp"

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
#include "io/sequence_iterator.hpp"

template <typename OT, bool buffering = false, bool preloading = false, bool populate = false, bool huge = false, bool advise = false>
struct readMMap {

    RangeType r;
    unsigned char* data;
    void *mapped_data;
    void *mapped;

    int file_handle;
    size_t file_size;
    size_t page_size;

    size_t start;
    size_t chunkSize;

    readMMap(std::string filename, int nprocs, int rank, int nthreads, int _chunkSize) : data(nullptr), mapped_data(nullptr), chunkSize(_chunkSize) {

      /// get file size.
      struct stat filestat;
      int ret = stat(filename.c_str(), &filestat);
      if (ret < 0 ) {
        ERROR( "ERROR in file open to get size." );
        exit(-1);
      }
      file_size = static_cast<size_t>(filestat.st_size);


      /// open the file and get a handle.
      file_handle = open(filename.c_str(), O_RDONLY);
      if (file_handle == -1)
      {
        int myerr = errno;
        FATAL("ERROR in file open: ["  << filename << "] error " << myerr << ": " << strerror(myerr));
        USED_BY_LOGGER_ONLY(myerr);

      }

      /// get the block size
      if (huge)
        page_size = 2 * 1024 * 1024;
      else
        page_size = sysconf(_SC_PAGE_SIZE);




      bliss::partition::BlockPartitioner<bliss::partition::range<size_t> > part;
      RangeType _r(0, file_size);
      part.configure(_r, nprocs);

      // configure mmap
      int advice =  MADV_SEQUENTIAL | MADV_WILLNEED;
      int flags = MAP_PRIVATE;
      if (populate) flags |= MAP_POPULATE;
      if (huge) {
        // ONLY OWKRS FOR PRIVATE, ANONYMOUS mappinig (not backed by file).
        flags |= MAP_HUGETLB;
        advice |= MADV_HUGEPAGE;
      }

      // mmap
      r = part.getNext(rank);
      typename RangeType::ValueType block_start = RangeType::align_to_page(r, page_size);
      typename RangeType::ValueType block_end = block_start + ((r.end - block_start + page_size - 1) / page_size) * page_size;
      block_end = std::min(file_size, block_end);

      mapped_data = mmap64(nullptr, block_end - block_start,
                                           PROT_READ,
                                           flags,
                                           file_handle,
                                           block_start);
      if (mapped_data == MAP_FAILED) {
        int myerr = errno;
        FATAL("ERROR in map: ["  << filename << "] error " << myerr << ": " << strerror(myerr) << "length " << block_end - block_start << " offset " << block_start);
        USED_BY_LOGGER_ONLY(myerr);

        exit(1);
      }

      data = reinterpret_cast<unsigned char*>(mapped_data) + (r.start - block_start);
      if (advise) madvise(mapped_data, block_end - block_start, advice);


      if (preloading)
      {
        auto rr = posix_memalign(&mapped, page_size, block_end - block_start);
        BLISS_UNUSED(rr);
        memcpy(mapped, mapped_data, block_end - block_start);  // copy some extra because of alignment issue.
        data = reinterpret_cast<unsigned char*>(mapped) + (r.start - block_start);

        munmap(mapped_data, block_end - block_start);
        mapped_data = mapped;
      }

      start = r.start;
    }

    ~readMMap() {

      // unmap
      if (preloading) {
        free(mapped);
      } else {

        typename RangeType::ValueType block_start = RangeType::align_to_page(r, page_size);
        typename RangeType::ValueType block_end = block_start + ((r.end - block_start + page_size - 1) / page_size) * page_size;
        block_end = std::min(file_size, block_end);

        munmap(mapped_data, block_end - block_start);
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
      { s = start; start += chunkSize; }

      typename RangeType::ValueType block_start = RangeType::align_to_page(r, page_size);
      typename RangeType::ValueType block_end = block_start + ((r.end - block_start + page_size - 1) / page_size) * page_size;
      block_end = std::min(file_size, block_end);

      RangeType r1(s, s+chunkSize);
      RangeType r2(s, s+2*chunkSize);
      r1.intersect(r);
      r2.intersect(r);


      if (r1.size() == 0) {
        //INFOF("%d empty at %lu - %lu, in range %lu - %lu! chunkSize %lu\n", tid, r2.start, r2.end, r.start, r.end, chunkSize);
        return true;
      }

      typename RangeType::ValueType aligned_s = RangeType::align_to_page(r2, page_size);
      typename RangeType::ValueType aligned_t = aligned_s + ((r2.end - aligned_s + page_size - 1) / page_size) * page_size;
      aligned_t = std::min(aligned_t, block_end);


      unsigned char * ld = data  + (r2.start - r.start);
      void *p;
      if (buffering) {
        auto rr = posix_memalign(&p, page_size, aligned_t - aligned_s);   // TODO: can preallocate.
        BLISS_UNUSED(rr);
        memcpy(p, reinterpret_cast<unsigned char*>(mapped_data) + (aligned_s - block_start), aligned_t - aligned_s);
        ld = reinterpret_cast<unsigned char*>(p) + (r2.start - aligned_s);
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
        free(p);

      count += lcount;

      v += tv;
      return false;
    }
};



template <typename OT, bool buffering = false, bool preloading = false>
struct readFileLoader {

    typedef bliss::io::FileLoader<unsigned char, 0, bliss::io::BaseFileParser, buffering, preloading> LoaderType;

    size_t page_size;
    typename LoaderType::L1BlockType data;

    LoaderType loader;

    readFileLoader(std::string filename, int nprocs, int rank, int nthreads, int chunkSize) :
      loader(filename, nprocs, rank, nthreads, chunkSize)  {

      loader.getNextL1Block();


      /// get the block size
      page_size = sysconf(_SC_PAGE_SIZE);

//      INFOF("loader from thread %d range = %lu %lu\n", omp_get_thread_num(), _r.start, _r.end);

      data = loader.getCurrentL1Block();
    }

    ~readFileLoader() {
      // unmap
    }

    RangeType getRange() {
      return data.getRange();
    }
    size_t getChunkSize() {
      return loader.getL2BlockSize();
    }


    std::string getName() {
      return "readFileLoader";
    }

    void reset() {
      loader.resetL2Partitioner();
    }

    bool operator()(int tid, size_t &count, OT &v) {
      RangeType r = data.getRange();
      RangeType r1 = loader.getNextL2Block(tid).getRange();

      r1.intersect(r);
      // try copying the data.

      if (r1.size() == 0)
        return true;   // finished

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

    typedef bliss::io::FileLoader<unsigned char, 0, bliss::io::BaseFileParser, buffering, preloading> LoaderType;


    size_t page_size;

    LoaderType loader;

    readFileLoaderAtomic(std::string filename, int nprocs, int rank, int nthreads, int chunkSize) :
      loader(filename, nprocs, rank, nthreads, chunkSize)  {

      loader.getNextL1Block();

      /// get the block size
      page_size = sysconf(_SC_PAGE_SIZE);

      // mmap
//      INFOF("loader from thread %d range = %lu %lu\n", omp_get_thread_num(), _r.start, _r.end);

    }

    ~readFileLoaderAtomic() {
      // unmap
    }

    RangeType getRange() {
      return loader.getCurrentL1Block().getRange();
    }
    size_t getChunkSize() {
      return loader.getL2BlockSize();
    }

    std::string getName() {
      return "readFileLoaderAtomic";
    }

    void reset() {
      loader.resetL2Partitioner();
    }

    bool operator()(int tid, size_t &count, OT &v) {

      // try copying the data.
      typename LoaderType::L2BlockType data = loader.getNextL2Block(tid);

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

      loader.getNextL1Block();

//      INFOF("reading: %lu - %lu, values '%c' '%c'\n", r.start, r.end, *(loader.getData().begin()), *(loader.getData().end()));



      /// get the block size
      page_size = sysconf(_SC_PAGE_SIZE);

      // mmap
//      INFOF("loader from thread %d range = %lu %lu\n", omp_get_thread_num(), _r.start, _r.end);

    }

    ~readFASTQ() {
    }

    RangeType getRange() {
      return loader.getCurrentL1Block().getRange();
    }
    size_t getChunkSize() {
      return loader.getL2BlockSize();
    }

    std::string getName() {
      return "readFASTQ";
    }

    void reset() {
      loader.resetL2Partitioner();
    }

    bool operator()(int tid, size_t &count, OT &v) {

      // try copying the data.
      typename LoaderType::L2BlockType data = loader.getNextL2Block(tid);

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
struct SequencesIterator {
    typedef bliss::io::FASTQLoader<unsigned char, buffering, preloading> LoaderType;

    size_t page_size;


    LoaderType loader;

    typedef float QualityType;
    typedef bliss::common::DNA Alphabet;
    typedef typename LoaderType::L2BlockType::iterator BaseIterType;

    typedef bliss::io::FASTQParser<BaseIterType>  ParserType;

    typedef bliss::io::SequencesIterator<BaseIterType, bliss::io::FASTQParser>           IteratorType;

    typedef typename ::std::iterator_traits<IteratorType>::value_type  SequenceType;

    SequencesIterator(std::string filename, int nprocs, int rank, int nthreads, int chunkSize) :
      loader(filename, nprocs, rank, nthreads, chunkSize)  {

      loader.getNextL1Block();

      /// get the block size
      page_size = sysconf(_SC_PAGE_SIZE);

      // mmap
//      INFOF("loader from thread %d range = %lu %lu\n", omp_get_thread_num(), _r.start, _r.end);

    }

    ~SequencesIterator() {
      // unmap
    }

    RangeType getRange() {
      return loader.getCurrentL1Block().getRange();
    }
    size_t getChunkSize() {
      return loader.getL2BlockSize();
    }

    std::string getName() {
      return "SequencesIterator";
    }


    void reset() {
      loader.resetL2Partitioner();
    }

    bool operator()(int tid, size_t &count, OT &v) {

      // try copying the data.
      typename LoaderType::L2BlockType data = loader.getNextL2Block(tid);

      if (data.begin() ==  data.end())
        return true;

      // traverse using fastq iterator.
      ParserType parser;
      IteratorType fastq_start(parser, data.begin(), data.end(), data.getRange().start);
      IteratorType fastq_end(data.end());

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
        for (BaseIterType iter = read.seqBegin; iter != read.seqEnd; ++iter) {
          c = std::max(*iter, c);
        }
        for (BaseIterType iter = read.qualBegin; iter != read.qualEnd; ++iter) {
          c = std::max(*iter, c);
        }

        for (BaseIterType iter = read.seqBegin; iter != read.seqEnd; ++iter) {
          d = std::min(*iter, d);
        }
        for (BaseIterType iter = read.qualBegin; iter != read.qualEnd; ++iter) {
          d = std::min(*iter, d);
        }

        // simulate kmer computation
        for (BaseIterType iter = read.seqBegin; iter != read.seqEnd; ++iter) {
          km <<= 8;
          km |= static_cast<uint64_t>(*iter);
        }
        for (BaseIterType iter = read.qualBegin; iter != read.qualEnd; ++iter) {
          km <<= 8;
          km |= static_cast<uint64_t>(*iter);
        }

        // simulate quality score computation.
        tv += static_cast<OT>(km) / static_cast<OT>(std::numeric_limits<uint64_t>::max() );
        for (BaseIterType iter = read.seqBegin; iter != read.seqEnd; ++iter) {
          tv += log2(*iter);
        }
        for (BaseIterType iter = read.qualBegin; iter != read.qualEnd; ++iter) {
          tv += log2(*iter);
        }

      }

      count += lcount;
      v += tv;
      return false;
    }
};


template <typename OT, bool buffering = false, bool preloading = false>
struct SequencesIterator2 {
    typedef bliss::io::FASTQLoader<unsigned char, buffering, preloading> LoaderType;

    size_t page_size;


    LoaderType loader;

    typedef float QualityType;
    typedef bliss::common::DNA Alphabet;
    typedef typename LoaderType::L2BlockType::iterator BaseIterType;

    typedef bliss::io::FASTQParser<BaseIterType>  ParserType;
    typedef bliss::io::SequencesIterator<BaseIterType, bliss::io::FASTQParser>           IteratorType;

    typedef typename ::std::iterator_traits<IteratorType>::value_type  SequenceType;


    SequencesIterator2(std::string filename, int nprocs, int rank, int nthreads, int chunkSize) :
      loader(filename, nprocs, rank, nthreads, chunkSize)  {

      loader.getNextL1Block();

//      INFOF("loader from thread %d loaded range = %lu %lu\n", omp_get_thread_num(), loader.getData().getRange().start, loader.getData().getRange().end);
//      INFOF("loader from thread %d file range = %lu %lu\n", omp_get_thread_num(), loader.getFileRange().start, loader.getFileRange().end);

      /// get the block size
      page_size = sysconf(_SC_PAGE_SIZE);

    }

    ~SequencesIterator2() {
      // unmap

    }

    RangeType getRange() {
      return loader.getCurrentL1Block().getRange();
    }
    size_t getChunkSize() {
      return loader.getL2BlockSize();
    }

    std::string getName() {
      return "SequencesIterator2";
    }

    void reset() {
      loader.resetL2Partitioner();
    }

    bool operator()(int tid, size_t &count, OT &v) {

      // try copying the data.
      typename LoaderType::L2BlockType data = loader.getNextL2Block(tid);

      if (data.begin() ==  data.end())
        return true;

      // traverse using fastq iterator.
      ParserType parser;
      IteratorType fastq_start(parser, data.begin(), data.end(), data.getRange().start);
      IteratorType fastq_end(data.end());

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
        iter = read.seqBegin;
        for (;
            iter != read.seqEnd;
            ++iter) {
          km <<= 8;
          km |= static_cast<uint64_t>(*iter);
//          ++i;
        }

        // simulate quality score computation.
        tv += static_cast<OT>(km) / static_cast<OT>(std::numeric_limits<uint64_t>::max() );
        iter = read.qualBegin;
        for (;   // slow because of log2.
            iter != read.qualEnd;
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
struct SequencesIteratorNoQual {
    typedef bliss::io::FASTQLoader<unsigned char, buffering, preloading> LoaderType;

    size_t page_size;


    LoaderType loader;

    typedef float QualityType;
    typedef bliss::common::DNA Alphabet;
    typedef typename LoaderType::L2BlockType::iterator BaseIterType;

    typedef bliss::io::FASTQParser<BaseIterType>  ParserType;
    typedef bliss::io::SequencesIterator<BaseIterType, bliss::io::FASTQParser>           IteratorType;
    typedef typename ::std::iterator_traits<IteratorType>::value_type  SequenceType;


    SequencesIteratorNoQual(std::string filename, int nprocs, int rank, int nthreads, int chunkSize) :
      loader(filename, nprocs, rank, nthreads, chunkSize)  {

      loader.getNextL1Block();

      /// get the block size
      page_size = sysconf(_SC_PAGE_SIZE);

      // mmap
//      INFOF("loader from thread %d range = %lu %lu\n", omp_get_thread_num(), _r.start, _r.end);

    }

    ~SequencesIteratorNoQual() {
      // unmap
    }

    RangeType getRange() {
      return loader.getCurrentL1Block().getRange();
    }
    size_t getChunkSize() {
      return loader.getL2BlockSize();
    }

    std::string getName() {
      return "SequencesIteratorNoQual";
    }
    void reset() {
      loader.resetL2Partitioner();
    }

    bool operator()(int tid, size_t &count, OT &v) {

      // try copying the data.
      typename LoaderType::L2BlockType data = loader.getNextL2Block(tid);
//      INFOF("thread %d getting block %lu-%lu, got block of length %ld\n", omp_get_thread_num(), start, end, (data.end() - data.begin()));


      if (data.begin() == data.end()) {
//        INFO( " range = " << start << "-" << end );
        return true;
      }

//      INFO( "Range: " << data.getRange() );
//      INFO( "Start: " << data.begin()[0] );
//      INFO( "End: "   << data.end()[0] );
//      INFO( "len: "   << (data.end() - data.begin())  << " len from range: " << (data.getRange().end - data.getRange().start) );


      // traverse using fastq iterator.
      ParserType parser;
      IteratorType fastq_start(parser, data.begin(), data.end(), data.getRange().start);
      IteratorType fastq_end(data.end());

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
        for (BaseIterType iter = read.seqBegin; iter != read.seqEnd; ++iter) {
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
  INFO( name << "\t" << tag <<"\tMPI rank: " << rank << "/" << nprocs << "\tOMP "
            << nthreads << " threads\ttook " << std::fixed
            << std::setprecision(6) << time_span.count() / iter
            << "s,\tresult = " << v << " count = " << count );
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
  INFO( "USE_MPI is set" );
#endif

  int nthreads = 1;
#ifdef USE_OPENMP
  if (rank == 0)
  INFO( "USE_OPENMP is set" );
  omp_set_nested(1);
  omp_set_dynamic(0);
  nthreads = omp_get_max_threads();
#endif


  if (argc > 1)
    nthreads = atoi(argv[1]);

  size_t step = 2 * 1024 * 1024;
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
  typedef readMMap<             double, true, false, false, false, false> OpType;
#endif

#if defined(TEST_OP_MMAP_ADVISE)
  typedef readMMap<             double, true, false, false, false, true> OpType;
#endif


#if defined(TEST_OP_MMAP_POPULATE)
  typedef readMMap<             double, true, false, true, false, false> OpType;
#endif

#if defined(TEST_OP_MMAP_POPULATE_ADVISE)
  typedef readMMap<             double, true, false, true, false, true> OpType;
#endif

// HUGETLB does not work with files.
//#if defined(TEST_OP_MMAP_HUGEPAGE)
//  typedef readMMap<             double, true, false, false, true, false> OpType;
//#endif
//
//#if defined(TEST_OP_MMAP_HUGEPAGE_ADVISE)
//  typedef readMMap<             double, true, false, false, true, true> OpType;
//#endif
//
//#if defined(TEST_OP_MMAP_POPULATE_HUGEPAGE)
//  typedef readMMap<             double, true, false, true, true, false> OpType;
//#endif
//
//#if defined(TEST_OP_MMAP_POPULATE_HUGEPAGE_ADVISE)
//  typedef readMMap<             double, true, false, true, true, true> OpType;
//#endif


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
  typedef ::SequencesIterator<        double, true, false> OpType;
#endif

#if defined(TEST_OP_FASTQIter2)
  typedef SequencesIterator2<       double, true, false> OpType;
#endif

#if defined(TEST_OP_FASTQIterNoQual)
  typedef SequencesIteratorNoQual<  double, true, false> OpType;
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
#pragma omp parallel OMP_SHARE_DEFAULT shared(nthreads, step, filename) num_threads(nthreads) reduction(+:v, count)
    {
      OpType op2(filename, nthreads, omp_get_thread_num(), 1, step);
      op2.reset();
      double v0 = Sequential<OpType, double>(op2, 1, count);
      //INFOF("%d processing range %lu %lu. result = %f\n", omp_get_thread_num(), op2.getRange().start, op2.getRange().end, v0);

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

