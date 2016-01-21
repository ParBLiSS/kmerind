/**
 * @file    mpi_runner.hpp
 * @ingroup wip
 * @author  Tony Pan <tpan7@gatech.edu>
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

#include "bliss-config.hpp"
#include "taskrunner/runner.hpp"

namespace bliss
{
namespace concurrent
{

/**
 * @class     bliss::concurrent::MPIRunner
 * @brief     runs a set of tasks on MPI processes
 * @details   note that each MPI process is in its own memory address
 *            space. We need to assign a single task to it based on its
 *            rank.
 *
 *            NOT USED.
 *
 */
template <bool threaded>
class MPIRunner : public Runner
{
    static_assert(false, "no good usecase for this class.  do not use.");

  protected:
    Runnable r;

  protected:
    MPI_Comm comm;
    bool commWasCreated;

  public:
    MPIRunner(int argc, char** argv) : Runner()
    {
#ifdef USE_MPI
      if (threaded) {
        int provided;
        MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

        // check if there is MPI thread support
        if (provided < MPI_THREAD_MULTIPLE) {
          ERRORF("ERROR: The MPI Library Does not have full thread support.\n");
          MPI_Abort(MPI_COMM_WORLD, 1);
        }
      } else {
        MPI_Init(&argc, &argv);
      }

      // get size and rank
      comm = MPI_COMM_WORLD;
      MPI_Comm_size(comm, &groupSize);
      MPI_Comm_rank(comm, &id);

      commWasCreated = true;

      if (id == 0)
        // TODO: replace with logging
        INFO( "USE_MPI is set" );
#else
      static_assert(false, "MPIRunner Used when compilation is not set to use MPI");
#endif
    }

    // TODO: a splitting comm constructor


    // comm is now handled by the MPIRunner.
    MPIRunner(MPI_Comm &_comm) : Runner()
    {
      comm = _comm;
      MPI_Comm_size(comm, &groupSize);
      MPI_Comm_rank(comm, &id);
      commWasCreated = false;
    }

    virtual ~MPIRunner()
    {
      if (commWasCreated) {
        MPI_Finalize();
      }
    }

    MPI_Comm & getComm()
    {
      return comm;
    }

    virtual void addTask(Runnable &t)
    {
      r = t;
    }

    virtual void run()
    {
      r.run();
    }

    virtual void synchronize()
    {
      MPI_Barrier(comm);
    }



};

} /* namespace concurrent */
} /* namespace bliss */

#endif /* MPIRUNNER_HPP_ */
