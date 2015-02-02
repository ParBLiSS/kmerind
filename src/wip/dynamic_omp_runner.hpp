/**
 * @file    dynamic_omp_runner.hpp
 * @ingroup
 * @author  tpan
 * @brief
 * @details
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef DYNAMIC_OMP_RUNNER_HPP_
#define DYNAMIC_OMP_RUNNER_HPP_

#include "omp.h"
#include <cassert>
#include "concurrent/lockfree_queue.hpp"

#include "config.hpp"
#include "wip/runner.hpp"


namespace bliss
{
namespace concurrent
{

/**
 * @class      bliss::concurrent::DynamicOMPRunner
 * @brief       OpenMP task runner, processing a queue of tasks until queue is done (blocked)
 * @details     used for continuously changing queue, threads process tasks in demand driven way.
 *              also use for one thread per task
 *
 */
class DynamicOMPRunner : public Runner
{
  protected:
    bliss::concurrent::ThreadSafeQueue<std::unique_ptr<Runnable> > q;

    const int nThreads;

  public:
    DynamicOMPRunner(int num_threads) : Runner(), nThreads(num_threads)
    {
#ifdef USE_OPENMP
//      omp_set_num_threads(num_threads);
      omp_set_nested(1);
#else
      static_assert(false, "OMPRunner Used When compilation is not set to use OpenMP");
#endif
    }


    virtual ~DynamicOMPRunner() {};

    void operator()()
    {
#pragma omp parallel num_threads(nThreads) default(none)
      {
#pragma omp single nowait
        {
          while (!q.isEmpty())
          {
#pragma omp task
            {
              auto v = std::move(q.waitAndPop());
              if (v.first)
                (v.second)->operator()();
            }
          }
        }
      }
    }

    virtual bool addTask(std::unique_ptr<Runnable> &&t)
    {
      auto result = q.waitAndPush(std::forward<std::unique_ptr<Runnable> >(t));
      return result.first;
    }

    /// flush currently queued tasks and disallow further new tasks
    void blockAdd() {
      q.disablePush();
    };

    /// flush currently queued tasks and disallow further new tasks
    void unblockAdd() {
      q.enablePush();
    };

    virtual void synchronize()
    {
#pragma omp barrier
    }

};

} /* namespace concurrent */
} /* namespace bliss */

#endif /* DYNAMIC_OMP_RUNNER_HPP_ */
