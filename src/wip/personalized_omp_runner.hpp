/**
 * @file    personalized_omp_runner.hpp
 * @ingroup
 * @author  tpan
 * @brief
 * @details
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef PERSONALIZED_OMP_RUNNER_HPP_
#define PERSONALIZED_OMP_RUNNER_HPP_


#include <cassert>

#include "concurrent/lockfree_queue.hpp"

#include "wip/runner.hpp"

namespace bliss
{
namespace concurrent
{

/**
 * @class     bliss::concurrent::PersonalizedOMPRunner
 * @brief     runs a sequence of tasks
 * @details   can have multiple tasks, run in order of adding.
 *            MORE than just a compound task, since this does NOT just have a single input and output.
 *
 */
class PersonalizedOMPRunner : public Runner
{
  protected:
    // threadsafe queue because addTask may be called from a different thread.
    bliss::concurrent::ThreadSafeQueue<std::shared_ptr<Runnable> > q;

    const int nThreads;

  public:
    PersonalizedOMPRunner(const int num_threads) : Runner(), nThreads(num_threads) {
#ifdef USE_OPENMP
      if (omp_get_nested() <= 0) omp_set_nested(1);
      if (omp_get_dynamic() <= 0) omp_set_dynamic(1);
#else
      static_assert(false, "OMPRunner Used When compilation is not set to use OpenMP");
#endif
    };

    virtual ~PersonalizedOMPRunner() {};

    void operator()()
    {
      this->disableAdd();   // do not allow more tasks to be added.

      // get the size, use it for parallel for.
      size_t count = q.getSize();

#pragma omp parallel for num_threads(nThreads) default(none) schedule(dynamic) shared(count)
      for (size_t i = 0; i < count; ++i)
      {
        auto v = std::move(q.tryPop());
        if (v.first) {
          (v.second)->operator()();
//          counter2.fetch_add(1);

        }
      }
      printf("Personalized runner completed.\n");
    }

    virtual bool addTask(std::shared_ptr<Runnable> &&t)
    {
   //   counter.fetch_add(1);

      auto result = q.waitAndPush(std::forward<std::shared_ptr<Runnable> >(t));
//      printf("add to Personalized runner.  size %lu, disabled %s\n", q.getSize(), (q.canPush() ? "n" : "y"));

      return result.first;
    }

    virtual size_t getTaskCount() {
      return q.getSize();
    }

    virtual bool isAddDisabled() {
      return !(q.canPush());
    }

    virtual void disableAdd() {
      q.disablePush();
    }

    virtual void synchronize()
    {
#pragma omp barrier
    }

};

} /* namespace concurrent */
} /* namespace bliss */

#endif /* PERSONALIZED_OMP_RUNNER_HPP_ */
