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

template <typename OT>
struct readMMap {

    RangeType r;
    unsigned char* data;
    unsigned char* mapped_data;

    int file_handle;
    size_t page_size;
    bool buffering;
    bool preloading;

    readMMap(std::string filename, RangeType _r, bool _preloading = false, bool _buffering = false) {
      /// open the file and get a handle.
      file_handle = open(filename.c_str(), O_RDONLY);
      if (file_handle == -1)
      {
        int myerr = errno;
        std::cerr << "ERROR in file open: ["  << filename << "] error " << myerr << ": " << strerror(myerr);
        exit(-1);
      }

      buffering = _buffering;
      /// get the block size
      page_size = sysconf(_SC_PAGE_SIZE);


      // mmap
      r = _r.align_to_page(page_size);
      mapped_data = (unsigned char*)mmap(nullptr, r.end - r.block_start,
                                           PROT_READ,
                                           MAP_PRIVATE, file_handle,
                                           r.block_start);
      data = mapped_data + r.start - r.block_start;

      preloading = _preloading;
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

    void reset() {
    }


    OT operator()(size_t start, size_t end) {
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
      for (size_t i = 0 ; i < e - s; ++i) {
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

      return tv;
    }
};


template <typename OT>
struct readFileLoader {

    RangeType r;
    unsigned char* data;
    size_t page_size;
    bool buffering;

    bliss::io::file_loader<unsigned char> loader;

    readFileLoader(std::string filename, RangeType _r, bool _preloading = false, bool _buffering = false) :
      buffering(_buffering), loader(filename)  {

      if (_preloading)
        loader.load(0.1f);
      else
        loader.load();

      /// get the block size
      page_size = sysconf(_SC_PAGE_SIZE);

      // mmap
      r = loader.getRange();
      data = loader.begin();
    }

    ~readFileLoader() {
      // unmap
      loader.unload();
      data = nullptr;
    }

    void reset() {
      loader.resetChunks();
    }

    OT operator()(size_t start, size_t end) {
      size_t s = std::max(r.start, start);
      size_t e = std::min(r.end, end);
      size_t e2 = std::min(r.end, end + (end-start));
      // try copying the data.

      unsigned char * ld = data + s - r.start;
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
      for (size_t i = 0 ; i < e - s; ++i) {
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

      return tv;
    }
};

template<typename CharT, typename SizeT = size_t>
struct PartitionHelper {
    typedef SizeT SizeType;

    SizeT operator()(const CharT* _data,
                            const SizeT &start, const SizeT &end) {
      return start;
    }
};


template <typename OT>
struct readFileLoaderAtomic {

    RangeType r;
    unsigned char* data;
    size_t page_size;
    bool buffering;

    bliss::io::file_loader<unsigned char> loader;
    PartitionHelper<unsigned char> helper;

    readFileLoaderAtomic(std::string filename, RangeType _r, bool _preloading = false, bool _buffering = false) :
      buffering(_buffering), loader(filename)  {

      if (_preloading)
        loader.load(0.1f);
      else
        loader.load();

      /// get the block size
      page_size = sysconf(_SC_PAGE_SIZE);

      // mmap
      r = loader.getRange();
      data = loader.begin();
    }

    ~readFileLoaderAtomic() {
      // unmap
      loader.unload();
      data = nullptr;
    }

    void reset() {
      loader.resetChunks();
    }

    OT operator()(size_t start, size_t end) {

      // try copying the data.

      unsigned char *startPtr, *endPtr;
      loader.getNextChunkAtomic(helper, startPtr, endPtr, end - start, buffering);

      if (startPtr == nullptr || endPtr == nullptr)
        std::cerr << " range = " << start << "-" << end << std::endl;


      // simulate multiple traversals.
      unsigned char c = 0;
      for (unsigned char *iter = startPtr; iter != endPtr; ++iter) {
        c = std::max(*iter, c);
      }

      unsigned char d = 255;

      for (unsigned char *iter = startPtr; iter != endPtr; ++iter) {
        d = std::min(*iter, d);
      }

      // simulate kmer computation
      uint64_t km = c;
      km += d;
      for (unsigned char *iter = startPtr; iter != endPtr; ++iter) {
        km <<= 8;
        km |= static_cast<uint64_t>(*iter);
      }

      // simulate quality score computation.
      OT tv = static_cast<OT>(km) / static_cast<OT>(std::numeric_limits<uint64_t>::max() );
      for (unsigned char *iter = startPtr; iter != endPtr; ++iter) {
        tv += log2(*iter);
      }

      if (buffering)
        delete [] startPtr;

      return tv;
    }
};


void printTiming(std::string tag, int rank, int nprocs, int nthreads,
                 const std::chrono::duration<double>& time_span, int iter,
                 double v)
{
  std::cout << tag << "\tMPI rank: " << rank << "/" << nprocs << "\tOMP "
            << nthreads << " threads\ttook " << std::fixed
            << std::setprecision(6) << time_span.count() / iter
            << "s,\tresult = " << v << std::endl;
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

  size_t step = 128;
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
  //readMMap<double> op(filename, r, false, true);
  //readFileLoader<double> op(filename, r, false, true);
  readFileLoaderAtomic<double> op(filename, r, false, true);
  double v = 0.0;

  std::chrono::high_resolution_clock::time_point t1, t2;
  std::chrono::duration<double> time_span;

  /// Workers only, critical
#if defined(USE_MPI)
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  t1 = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < iter; ++i) {
    op.reset();
    v = P2P(op, nthreads, r, step);
  }
#if defined(USE_MPI)
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  t2 = std::chrono::high_resolution_clock::now();
  time_span =
      std::chrono::duration_cast<std::chrono::duration<double>>(
          t2 - t1);
  if (rank == 0)
  printTiming("P2P critical:", rank, nprocs, nthreads, time_span, iter, v);


  /// Workers only, atomic  - close to time for parfor.
#if defined(USE_MPI)
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  t1 = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < iter; ++i) {
    op.reset();
    v = P2P_atomic(op, nthreads, r, step);
  }
#if defined(USE_MPI)
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  t2 = std::chrono::high_resolution_clock::now();
  time_span =
      std::chrono::duration_cast<std::chrono::duration<double>>(
          t2 - t1);
  if (rank == 0)
  printTiming("P2P atomic:", rank, nprocs, nthreads, time_span, iter, v);


  /// master slave
#if defined(USE_MPI)
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  t1 = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < iter; ++i) {
    op.reset();
    v= MasterSlave(op, nthreads, r, step);
  }
#if defined(USE_MPI)
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  t2 = std::chrono::high_resolution_clock::now();
  time_span =
      std::chrono::duration_cast<std::chrono::duration<double>>(
          t2 - t1);
  if (rank == 0)
  printTiming("MS Wait:", rank, nprocs, nthreads, time_span, iter, v);

  /// master slave No Wait
#if defined(USE_MPI)
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  t1 = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < iter; ++i) {
    op.reset();
    v= MasterSlaveNoWait(op, nthreads, r, step);
  }
#if defined(USE_MPI)
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  t2 = std::chrono::high_resolution_clock::now();
  time_span =
      std::chrono::duration_cast<std::chrono::duration<double>>(
          t2 - t1);
  if (rank == 0)
  printTiming("MS NoWait:", rank, nprocs, nthreads, time_span, iter, v);

  /// parallel for
#if defined(USE_MPI)
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  t1 = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < iter; ++i) {
    op.reset();
    v = ParFor(op, nthreads, r, step);
  }
#if defined(USE_MPI)
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  t2 = std::chrono::high_resolution_clock::now();
  time_span =
      std::chrono::duration_cast<std::chrono::duration<double>>(
          t2 - t1);
  if (rank == 0)
  printTiming("PARFOR:\t", rank, nprocs, nthreads, time_span, iter, v);

  //// serial for
#if defined(USE_MPI)
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  t1 = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < iter; ++i) {
    op.reset();
    v = Sequential(op, nthreads, r, step);
  }
#if defined(USE_MPI)
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  t2 = std::chrono::high_resolution_clock::now();
  time_span =
      std::chrono::duration_cast<std::chrono::duration<double>>(
          t2 - t1);
  if (rank == 0)
  printTiming("SEQFOR:\t", rank, nprocs, nthreads, time_span, iter, v);


#ifdef USE_MPI
  MPI_Finalize();

#endif

  return 0;

}

