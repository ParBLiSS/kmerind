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
  protected:
    /// The task
    Runnable r;


  public:

    /**
     * @brief Creates a new OMPRunner instance with the given number of
     *        threads
     *
     * @param nThreads The number of OpenMP threads to use.
     */
#ifdef USE_OPENMP
    OMPRunner(int nThreads) : Runner(omp_get_thread_num(), nThreads) {
      omp_set_num_threads(nThreads);
    }
#else
    OMPRunner(int nThreads) : Runner() {
      static_assert(false, "OMPRunner used although compilation is not set to use OpenMP");
    }
#endif

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
      this->r = t;
    }

    /**
     * @brief Runs all tasks.
     */
    virtual void run() {
      Runnable lr = r;
#pragma omp parallel num_threads(groupSize) default(none) firstprivate(lr)
      this->r.run();
    }

    /**
     * @brief Synchronizes all threads using a barrier.
     */
    virtual void synchronize() {
#pragma omp barrier
    }

};

} /* namespace concurrent */
} /* namespace bliss */

#endif /* OMP_RUNNER_HPP_ */
