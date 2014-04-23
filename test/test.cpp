
#include <iostream>
#include <fstream>
#include <stdio.h>

#include "mpi.h"
#include <omp.h>

#include "utils/logging.h"
#include <cmath>
#include <chrono>

#define USE_OPENMP 1

int main(int argc, char* argv[])
{
  LOG_INIT();

  int rank = 0, nprocs = 1;

#ifdef USE_MPI
  std::cout << "USE_MPI is set" << std::endl;

  // initialize MPI
  MPI_Init(&argc, &argv);

  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

#endif

  int nthreads = 1;
  if (argc > 2) {
    nthreads = atoi(argv[2]);
  }

  int step = 128;
  if (argc > 3) {
    step = atoi(argv[3]);
  }

  std::chrono::high_resolution_clock::time_point t1, t2;
  std::chrono::duration<double> time_span;

  t1 = std::chrono::high_resolution_clock::now();

  size_t i = 0;
  size_t max = 1000000;
  if (argc > 1)
    max = atoi(argv[1]);
  bool done = false;
  double v = 0;


#ifdef USE_OPENMP
  std::cout << "USE_OPENMP is set" << std::endl;

#pragma omp parallel num_threads(nthreads) shared(max, done, i, step, v) default(none)
  {
#endif


//      printf("MPI rank: %d/%d OMP thread: %d/%d\n", (rank+1), nprocs, (tid+1), nthreads);

    double tv = 0;
    size_t j = 0;
    size_t bound = 0;

    do {

#ifdef USE_OPENMP
#pragma omp critical (assign)
#endif
      {
        j = i;
        i += step;
      }

      tv = 0;
      bound = (j + step > max ? max : j + step);
      for ( ; j < bound; ++j) {
        tv += log2(j + 1);
        //ofs << j << std::endl;
      }

      if (j >= max) {
        done = true;
#ifdef USE_OPENMP
#pragma omp flush(done)
#endif
      }


#pragma omp atomic
      v += tv;

    } while (!done);



#ifdef USE_OPENMP
  }
#endif


  t2 = std::chrono::high_resolution_clock::now();
  time_span =
      std::chrono::duration_cast<std::chrono::duration<double>>(
          t2 - t1);

  INFO(
      "P2P: MPI rank: " << rank << "/" << nprocs << " OMP took " << time_span.count() << "s, result = " << v);


  // finished peer to peer

  /// start master slave
  t1 = std::chrono::high_resolution_clock::now();
  v = 0;
#ifdef USE_OPENMP
#pragma omp parallel num_threads(nthreads) shared(max, v, step) default(none)
  {
#endif

#pragma omp single
    {
#pragma omp task untied
      for (size_t i = 0; i < max; i += step)
      {
#pragma omp task
        {
          double tv = 0;
          size_t bound = (i +step < max ? i+step : max);
          size_t j = i;
          for ( ; j < bound; ++j )
            tv += log2(j+1);
#pragma omp atomic
          v += tv;
//          printf("%d %lu %lf\n", omp_get_thread_num(), j, tv);
        }

      } //for
    } // omp single


#ifdef USE_OPENMP
  }
#endif


  t2 = std::chrono::high_resolution_clock::now();
  time_span =
      std::chrono::duration_cast<std::chrono::duration<double>>(
          t2 - t1);

  INFO(
      "MS: MPI rank: " << rank << "/" << nprocs << " OMP took " << time_span.count() << "s, result = " << v);


  // finished peer to peer

  /// start master slave
  t1 = std::chrono::high_resolution_clock::now();

  v = 0;
#ifdef USE_OPENMP
#pragma omp parallel for num_threads(nthreads) schedule(dynamic,1024) default(none) shared(max, rank, nprocs) reduction(+ : v)
#endif
for (size_t j = 0; j < max; ++j)
  {


//      printf("MPI rank: %d/%d OMP thread: %d/%d\n", (rank+1), nprocs, (tid+1), nthreads);


        v += log2(j + 1);
        //ofs << j << std::endl;
  }
  t2 = std::chrono::high_resolution_clock::now();
  time_span =
      std::chrono::duration_cast<std::chrono::duration<double>>(
          t2 - t1);

  INFO(
      "PARFOR MPI rank: " << rank << "/" << nprocs << " OMP took " << time_span.count() << "s, result = " << v);



  //// serial

  // finished peer to peer

  /// start master slave
  t1 = std::chrono::high_resolution_clock::now();

  v = 0;
  for (size_t j = 0; j < max; ++j)
  {


//      printf("MPI rank: %d/%d OMP thread: %d/%d\n", (rank+1), nprocs, (tid+1), nthreads);


        v += log2(j + 1);
        //ofs << j << std::endl;
  }
  t2 = std::chrono::high_resolution_clock::now();
  time_span =
      std::chrono::duration_cast<std::chrono::duration<double>>(
          t2 - t1);

  INFO(
      "SEQFOR MPI rank: " << rank << "/" << nprocs << " OMP took " << time_span.count() << "s, result = " << v);



#ifdef USE_MPI
  MPI_Finalize();

#endif

  return 0;

}

