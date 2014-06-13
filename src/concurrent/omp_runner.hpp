/**
 * @file    omp_runner.hpp
 * @ingroup
 * @author  tpan
 * @brief
 * @details
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef OMP_RUNNER_HPP_
#define OMP_RUNNER_HPP_

#include "omp.h"
#include <cassert>

#include "config.hpp"
#include "concurrent/runner.hpp"

namespace bliss
{
namespace concurrent
{

/**
 * @class      bliss::concurrent::OMPRunner
 * @brief
 * @details
 *
 */
class OMPRunner : public Runner
{
  public:

    /**
     * @brief Creates a new OMPRunner instance with the given number of
     *        threads
     *
     * @param nThreads The number of OpenMP threads to use.
     */
    OMPRunner(int nThreads) : groupSize(nThreads) {
#ifdef USE_OPENMP
      omp_set_num_threads(nThreads);
      id = omp_get_thread_num();
#else
      static_assert(false, "OMPRunner used although compilation is not set to use OpenMP");
#endif
    }

    /**
     * @brief The destructor of this class.
     */
    virtual ~OMPRunner() {};

    /**
     * @brief Adds a task to this Runner.
     *
     * @param t
     */
    virtual void addTask(Runnable &t) {
      r = t;
    }

    /**
     * @brief Runs all tasks.
     */
    virtual void run() {
#pragma omp parallel num_threads(groupSize) default(none) firstprivate(r)
      r.run();
    }

    /**
     * @brief Synchronizes all threads using a barrier.
     */
    virtual void synchronize() {
#pragma omp barrier
    }

  protected:
    /// The task
    Runnable r;
};

} /* namespace concurrent */
} /* namespace bliss */

#endif /* OMP_RUNNER_HPP_ */
