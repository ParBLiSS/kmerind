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

#include "config.hpp"
#include "concurrent/lockfree_queue.hpp"

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
    bliss::concurrent::ThreadSafeQueue<std::shared_ptr<Runnable> > q;

    const int nThreads;

  public:
    DynamicOMPRunner(const int num_threads) : Runner(), nThreads(num_threads)
    {
#ifdef USE_OPENMP
      if (omp_get_nested() <= 0) omp_set_nested(1);
      if (omp_get_dynamic() <= 0) omp_set_dynamic(1);
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
        {  // one thread to do this
//#pragma omp task untied default(none)                     // this slows it down by ORDERS OF MAGNITUDE.  DO NOT USE
//          {  // untied so can move to different threads
            while (q.canPop())  // flush the queue now.
            {
#pragma omp task default(none)
              {
                auto v = std::move(q.waitAndPop());
                if (v.first) {
                  (v.second)->operator()();
//                  counter2.fetch_add(1);

                }
              } // omp task
            } // while

            printf("tid %d done generating tasks.\n", omp_get_thread_num());
//          } // omp task untied.
        } // omp single
      } // omp parallel
      printf("Dynamic runner completed.\n");

    }  // operator()



    virtual bool addTask(std::shared_ptr<Runnable> &&t)
    {
//      counter.fetch_add(1);

        auto result = q.waitAndPush(std::forward<std::shared_ptr<Runnable> >(t));
//        printf("add to Dynamic runner.  size %lu, disabled %s\n", q.getSize(), (q.canPush() ? "n" : "y"));
        return result.first;
    }

    virtual size_t getTaskCount() {
      return q.getSize();
    }

    virtual bool isAddDisabled() {
      return !(q.canPush());
    }

    /// flush currently queued tasks and disallow further new tasks
    virtual void disableAdd() {
      q.disablePush();
    };

    virtual void synchronize()
    {
#pragma omp barrier
    }

};

} /* namespace concurrent */
} /* namespace bliss */

#endif /* DYNAMIC_OMP_RUNNER_HPP_ */
