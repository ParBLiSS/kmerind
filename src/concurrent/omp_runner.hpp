/**
 * @file		omp_runner.hpp
 * @ingroup
 * @author	tpan
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
     * @class			bliss::concurrent::OMPRunner
     * @brief
     * @details
     *
     */
    class OMPRunner : public Runner
    {
      public:
        OMPRunner(int nThreads) : groupSize(nThreads) {
#ifdef USE_OPENMP
          omp_set_num_threads(nThreads);
          id = omp_get_thread_num();
#else
          static_assert(false, "OMPRunner Used When compilation is not set to use OpenMP");
#endif

        }
        virtual ~OMPRunner() {};

        virtual void addTask(Runnable &t) {
          r = t;
        }

        virtual void run() {

#pragma omp parallel num_threads(groupSize) default(none) firstprivate(r)
          r.run();
        }

        virtual void synchronize() {
#pragma omp barrier
        }

      protected:
        Runnable r;

    };

  } /* namespace concurrent */
} /* namespace bliss */

#endif /* OMP_RUNNER_HPP_ */
