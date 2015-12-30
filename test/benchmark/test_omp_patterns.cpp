#include "bliss-config.hpp"

#include <iostream> // cout
#include <iomanip>  // for setprecision

#ifdef USE_MPI
#include "mpi.h"
#endif

#include <cmath>      // log
#include <chrono>     // timing
#include <algorithm>  // for std::min

#include "omp_patterns.hpp"
#include "partition/partitioner.hpp"
#include "utils/logging.h"

template <typename OT>
struct compute {
    bool operator()(int tid, size_t count, OT &v) {
      v = 0;
      for (size_t i = 0 ; i < count; ++i) {
        v += log2(static_cast<double>(i) + 1);
        //ofs << j << std::endl;
      }
      return true;
    }

    RangeType getRange() {
      return RangeType(0, 100);
    }

    size_t getChunkSize() {
      return 100;
    }
};


void printTiming(std::string tag, int rank, int nprocs, int nthreads,
                 const std::chrono::duration<double>& time_span, int iter,
                 double v, size_t &count)
{
  BL_INFO( tag << "\tMPI rank: " << rank << "/" << nprocs << "\tOMP "
            << nthreads << " threads\ttook " << std::fixed
            << std::setprecision(6) << time_span.count() / iter
            << "s,\tresult = " << v << ",\tcount" << count );
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
  BL_INFO( "USE_MPI is set" );

#endif

#ifdef USE_OPENMP
  if (rank == 0)
  BL_INFO( "USE_OPENMP is set" );
  omp_set_nested(1);
  omp_set_dynamic(0);
#endif


  int nthreads = 1;
  if (argc > 1)
    nthreads = atoi(argv[1]);

  size_t max = 1000000;
  if (argc > 2)
    max = atoi(argv[2]);

  int iter = 10;
  if (argc > 3)
    iter = atoi(argv[3]);

  bliss::partition::BlockPartitioner<RangeType> part;
  part.configure(RangeType(0, max), nprocs);

  compute<double> op;
  double v = 0.0;
  size_t count = 0;

  std::chrono::high_resolution_clock::time_point t1, t2;
  std::chrono::duration<double> time_span;

  /// Workers only, critical
#if defined(USE_MPI)
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  t1 = std::chrono::high_resolution_clock::now();
  count = 0;
  for (int i = 0; i < iter; ++i)
    v = P2P<compute<double>, double>(op, nthreads, count);
#if defined(USE_MPI)
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  t2 = std::chrono::high_resolution_clock::now();
  time_span =
      std::chrono::duration_cast<std::chrono::duration<double>>(
          t2 - t1);
  if (rank == 0)
    printTiming("P2P critical:", rank, nprocs, nthreads, time_span, iter, v, count);



  /// master slave
#if defined(USE_MPI)
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  t1 = std::chrono::high_resolution_clock::now();
  count = 0;
  for (int i = 0; i < iter; ++i)
    v= MasterSlave<compute<double>, double>(op, nthreads, count);
#if defined(USE_MPI)
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  t2 = std::chrono::high_resolution_clock::now();
  time_span =
      std::chrono::duration_cast<std::chrono::duration<double>>(
          t2 - t1);
  if (rank == 0)
    printTiming("MS Wait:", rank, nprocs, nthreads, time_span, iter, v, count);

  /// master slave No Wait
#if defined(USE_MPI)
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  t1 = std::chrono::high_resolution_clock::now();
  count = 0;
  for (int i = 0; i < iter; ++i)
    v= MasterSlaveNoWait<compute<double>, double>(op, nthreads, count);
#if defined(USE_MPI)
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  t2 = std::chrono::high_resolution_clock::now();
  time_span =
      std::chrono::duration_cast<std::chrono::duration<double>>(
          t2 - t1);
  if (rank == 0)
    printTiming("MS NoWait:", rank, nprocs, nthreads, time_span, iter, v, count);

  /// parallel for
#if defined(USE_MPI)
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  t1 = std::chrono::high_resolution_clock::now();
  count = 0;
  for (int i = 0; i < iter; ++i)
    v = ParFor<compute<double>, double>(op, nthreads, count);
#if defined(USE_MPI)
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  t2 = std::chrono::high_resolution_clock::now();
  time_span =
      std::chrono::duration_cast<std::chrono::duration<double>>(
          t2 - t1);
  if (rank == 0)
    printTiming("PARFOR:\t", rank, nprocs, nthreads, time_span, iter, v, count);

  //// serial for
#if defined(USE_MPI)
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  t1 = std::chrono::high_resolution_clock::now();
  count = 0;
  for (int i = 0; i < iter; ++i)
    v = Sequential<compute<double>, double>(op, nthreads, count);
#if defined(USE_MPI)
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  t2 = std::chrono::high_resolution_clock::now();
  time_span =
      std::chrono::duration_cast<std::chrono::duration<double>>(
          t2 - t1);
  if (rank == 0)
    printTiming("SEQFOR:\t", rank, nprocs, nthreads, time_span, iter, v, count);

#ifdef USE_MPI
  MPI_Finalize();

#endif

  return 0;

}

