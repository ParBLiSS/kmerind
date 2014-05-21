/**
 * @file		mpi_runner.hpp
 * @ingroup
 * @author	tpan
 * @brief   runs some MPI work.  has a source, a compute algo, and a sink.
 * @details can either initialize directly, or use an existing MPI communicator.
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef MPIRUNNER_HPP_
#define MPIRUNNER_HPP_

#include "mpi.h"
#include "cassert"

#include "config.hpp"

namespace bliss
{
  namespace concurrent
  {

    /**
     * @class			bliss::concurrent::MPIRunner
     * @brief
     * @details
     *
     */
    template <bool threaded>
    class MPIRunner
    {
      protected:
        MPI_Comm comm;
        int nprocs;
        int rank;
        bool commWasCreated;


      public:
        MPIRunner(int argc, char** argv) {

#if defined(USE_MPI)
          if (threaded) {

            int provided;
            MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

            if (provided < MPI_THREAD_MULTIPLE) {
              printf("ERROR: The MPI Library Does not have full thread support.\n");
              MPI_Abort(MPI_COMM_WORLD, 1);
            }
          } else {
            MPI_Init(&argc, &argv);
          }

          comm = MPI_COMM_WORLD;

          MPI_Comm_size(comm, &nprocs);
          MPI_Comm_rank(comm, &rank);

          commWasCreated = true;

          if (rank == 0)
            std::cout << "USE_MPI is set" << std::endl;


#else
          //TODO:  need to support no MPI.
          fprintf(stderr, "ERROR: need MPI support\n");
          exit(1);
#endif
        };

        // TODO: a splitting comm constructor


        // comm is now handled by the MPIRunner.
        MPIRunner(MPI_Comm &_comm) {
          comm = _comm;
          MPI_Comm_size(comm, &nprocs);
          MPI_Comm_rank(comm, &rank);
          commWasCreated = false;

        };

        virtual ~MPIRunner() {
          if (commWasCreated) {
            MPI_Finalize();
          }

        };
    };

  } /* namespace concurrent */
} /* namespace bliss */

#endif /* MPIRUNNER_HPP_ */
