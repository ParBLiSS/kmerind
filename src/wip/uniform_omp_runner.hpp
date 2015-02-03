/**
 * @file    uniform_omp_runner.hpp
 * @ingroup
 * @author  tpan
 * @brief
 * @details
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef UNIFORM_OMP_RUNNER_HPP_
#define UNIFORM_OMP_RUNNER_HPP_

#include "omp.h"
#include <cassert>

#include "config.hpp"
#include "wip/runner.hpp"
#include "concurrent/lockfree_queue.hpp"
#include <vector>

namespace bliss
{
namespace concurrent
{

/**
 * @class      bliss::concurrent::UniformOMPRunner
 * @brief       runner that performs the same task(s) for each thread.
 * @details     if the task is compound, a sequential runner can be used.
 *
 */
class UniformOMPRunner : public Runner
{
  protected:
    /// The task
    std::vector<bliss::concurrent::ThreadSafeQueue<Runnable* > > qs;

    const int nThreads;

  public:

    /**
     * @brief Creates a new UniformOMPRunner instance with the given number of
     *        threads
     *
     * @param nThreads The number of OpenMP threads to use.
     */
    UniformOMPRunner(int num_threads) : Runner(), qs(num_threads), nThreads(num_threads) {
#ifdef USE_OPENMP
      if (omp_get_nested() <= 0) omp_set_nested(1);
      if (omp_get_dynamic() <= 0) omp_set_dynamic(1);
#else
      static_assert(false, "UniformOMPRunner used although compilation is not set to use OpenMP");
#endif

      qs.clear();
      for (int i = 0; i < nThreads; ++i) {
        qs.push_back(std::move(bliss::concurrent::ThreadSafeQueue<Runnable* >()));
      }
    }

    /**
     * @brief The destructor of this class.
     */
    virtual ~UniformOMPRunner() {};

    /**
     * @brief Runs all tasks.
     */
    void operator()() {

#pragma omp parallel num_threads(nThreads) default(none)
      {
        // iterator for list is valid during append (no insertion in middle, no deletion).
        // no deletion since this is shared.
        bliss::concurrent::ThreadSafeQueue<Runnable* >& q = qs[omp_get_thread_num()];
        while (q.canPop())
        {
          auto v = std::move(q.waitAndPop());
          if (v.first) {
            (v.second)->operator()();
  //              } else {
  //                printf("blocked nothing in queue\n");
          }
        }
      }
    }

    /**
     * @brief Adds a task to this Runner.
     *
     * @param t
     */
    virtual bool addTask(Runnable* t) {
      bool out = true;
      int id;
      for (int i = 0; i < nThreads; ++i) {
        id = (omp_get_thread_num() + i) % nThreads;
        bliss::concurrent::ThreadSafeQueue<Runnable* >& q = qs[id];
        auto result = q.waitAndPush(std::forward<Runnable* >(t));

        out &= result.first;
      }
      return out;

    }

    virtual void disableAdd() {
      int id;
      for (int i = 0; i < nThreads; ++i) {
        id = (omp_get_thread_num() + i) % nThreads;
        bliss::concurrent::ThreadSafeQueue<Runnable* >& q = qs[id];
        q.disablePush();
      }
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

#endif /* UNIFORM_OMP_RUNNER_HPP_ */
