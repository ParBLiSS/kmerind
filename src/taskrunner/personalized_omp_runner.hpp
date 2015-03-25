/**
 * @file    personalized_omp_runner.hpp
 * @ingroup taskrunner
 * @author  Tony Pan <tpan7@gatech.edu>
 * @brief	an OMP task execution engine for a fixed list of tasks
 * @details	the threads block decompose or round robin the task execution
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef PERSONALIZED_OMP_RUNNER_HPP_
#define PERSONALIZED_OMP_RUNNER_HPP_


#include "omp.h"
#include <cassert>
#include "config.hpp"

#include "concurrent/mutexlock_queue.hpp"
#include "taskrunner/runner.hpp"

namespace bliss
{
namespace concurrent
{

/**
 * @class     bliss::concurrent::PersonalizedOMPRunner
 * @brief     runs a sequence of tasks, one per thread.  does not allow adding tasks during execution
 * @details   can have multiple tasks, run in order of adding.
 *            the tasks are assigned to threads in a round robin manner or block decomposing manner,
 *            depending on the OMP scheduler.
 *
 *            This runner differs from the dynamic OMP runner in that when
 *            task engine starts, the task queue no longer accepts new tasks.
 *
 */
class PersonalizedOMPRunner : public Runner
{
  protected:
    // threadsafe queue because addTask may be called from a different thread.
    bliss::concurrent::ThreadSafeQueue<std::shared_ptr<Runnable> , bliss::concurrent::LockType::MUTEX> q;

    const int nThreads;

  public:
    /// constructor
    PersonalizedOMPRunner(const int num_threads) : Runner(), nThreads(num_threads) {
#ifdef USE_OPENMP
      if (omp_get_nested() <= 0) omp_set_nested(1);
      if (omp_get_dynamic() <= 0) omp_set_dynamic(1);
#else
      static_assert(false, "OMPRunner Used When compilation is not set to use OpenMP");
#endif
    };

    /// destructor
    virtual ~PersonalizedOMPRunner() {};

    /**
     * @brief executes the tasks
     * @details   This method first blocks the queue from accepting new tasks
     * 	then it iterates over all tasks and assign them to threads in blocks
     * 	in round robin fashion.
     */
    void operator()()
    {
      this->disableAdd();   // do not allow more tasks to be added.

      // get the size, use it for parallel for.
      size_t count = q.getSize();

      size_t proc = 0;

#pragma omp parallel for num_threads(nThreads) default(none) schedule(dynamic) shared(count) reduction(+: proc)
      for (size_t i = 0; i < count; ++i)
      {
        auto v = std::move(q.tryPop());
        if (v.first) {
        	// run the task
          (v.second)->operator()();
          ++proc;
        }
      }
      INFOF("Personalized runner generated %lu tasks and completed %lu tasks.\n", count, proc);
    }

    /// add a new task to queue
    virtual bool addTask(std::shared_ptr<Runnable> &&t)
    {
      DEBUGF("add to Personalized runner.  size %lu, disabled %s\n", q.getSize(), (q.canPush() ? "n" : "y"));
      return q.waitAndPush(std::forward<std::shared_ptr<Runnable> >(t)).first;
    }

    /// count the number of pending tasks
    virtual size_t getTaskCount() {
      return q.getSize();
    }

    /// check to see if new tasks can be added to the queue
    virtual bool isAddDisabled() {
      return !(q.canPush());
    }

    /// flush currently queued tasks and disallow further new tasks
    virtual void disableAdd() {
      q.disablePush();
    }

    /// allows synchronization (not used right now)
    virtual void synchronize()
    {
#pragma omp barrier
    }

};

} /* namespace concurrent */
} /* namespace bliss */

#endif /* PERSONALIZED_OMP_RUNNER_HPP_ */
