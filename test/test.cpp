
#include "config.hpp"

#include <iostream>
#include <stdio.h>

#include "mpi.h"
#include <omp.h>

#include "utils/logging.h"

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

    int nthreads = 1, tid = 0;

#ifdef USE_OPENMP
    std::cout << "USE_OPENMP is set" << std::endl;

    // initialize openmp
#pragma omp parallel private(nthreads, tid)
    {
      nthreads = omp_get_num_threads();
      tid = omp_get_thread_num();
#endif

      INFO("MPI rank: " << (rank + 1) << "/" << nprocs << " OMP Thread: " << (tid+1) << "/" << nthreads);

//      printf("MPI rank: %d/%d OMP thread: %d/%d\n", (rank+1), nprocs, (tid+1), nthreads);


#ifdef USE_OPENMP
    }
#endif



#ifdef USE_MPI
    MPI_Finalize();

#endif


    return 0;

}

