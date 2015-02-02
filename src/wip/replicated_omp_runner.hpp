/**
 * @file    replicated_omp_runner.hpp
 * @ingroup
 * @author  tpan
 * @brief
 * @details
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef REPLICATED_OMP_RUNNER_HPP_
#define REPLICATED_OMP_RUNNER_HPP_

#include "omp.h"
#include <cassert>

#include "config.hpp"
#include "wip/runner.hpp"
#include <list>

namespace bliss
{
namespace concurrent
{

/**
 * @class      bliss::concurrent::ReplicatedOMPRunner
 * @brief       runner that performs the same task for each thread.
 * @details     if the task is compound, a sequential runner can be used.
 *
 */
class ReplicatedOMPRunner : public Runner
{
  protected:
    /// The task
    std::list<std::unique_ptr<Runnable> > q;
    const int nThreads;

  public:

    /**
     * @brief Creates a new ReplicatedOMPRunner instance with the given number of
     *        threads
     *
     * @param nThreads The number of OpenMP threads to use.
     */
    ReplicatedOMPRunner(int num_threads) : Runner(), nThreads(num_threads) {
#ifdef USE_OPENMP
      omp_set_nested(1);
//      omp_set_num_threads(num_threads);
#else
      static_assert(false, "ReplicatedOMPRunner used although compilation is not set to use OpenMP");
#endif
    }

    /**
     * @brief The destructor of this class.
     */
    virtual ~ReplicatedOMPRunner() {};

    /**
     * @brief Runs all tasks.
     */
    void operator()() {
#pragma omp parallel num_threads(nThreads) default(none)
      {
        // iterator for list is valid during append (no insertion in middle, no deletion).
        // no deletion since this is shared.
        auto iter = q.begin();
        while (iter != q.end()) {
          (*iter)->operator()();
          ++iter;
        }
      }
    }

    /**
     * @brief Adds a task to this Runner.
     *
     * @param t
     */
    virtual bool addTask(std::unique_ptr<Runnable> &&t) {
      if (!blocked.load(std::memory_order_relaxed)) {
        q.push_back(std::forward<std::unique_ptr<Runnable> >(t));
        return true;
      } else
        return false;
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

#endif /* REPLICATED_OMP_RUNNER_HPP_ */
