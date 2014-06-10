/**
 * @file    mixed_omp_runner.hpp
 * @ingroup
 * @author  tpan
 * @brief
 * @details
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef MIXED_OMP_RUNNER_HPP_
#define MIXED_OMP_RUNNER_HPP_

#include "omp.h"
#include <cassert>

#include "config.hpp"
#include "concurrent/runner.hpp"

namespace bliss
{
namespace concurrent
{

/**
 * @class      bliss::concurrent::MixedOMPRunner
 * @brief
 * @details
 *
 */
class MixedOMPRunner : public Runner
{
  public:
    MixedOMPRunner(int nThreads) : groupSize(nThreads)
    {
#ifdef USE_OPENMP
      omp_set_num_threads(nThreads);
      id = omp_get_thread_num();
#else
      static_assert(false, "OMPRunner Used When compilation is not set to use OpenMP");
#endif

    }
    virtual ~MixedOMPRunner() {};

    virtual void addTask(Runnable &t)
    {
      q.push_back(t);
    }

    virtual void run()
    {
#pragma omp parallel num_threads(groupSize) default(none) shared(q)
      {
#pragma omp single nowait
        {
          for (Runnable r : q)
          {
#pragma omp task
            {
              r.run();
            }
          }
        }
      }
    }

    virtual void synchronize()
    {
#pragma omp barrier
    }

  protected:
    std::vector<Runnable> q;
};

} /* namespace concurrent */
} /* namespace bliss */

#endif /* MIXED_OMP_RUNNER_HPP_ */
