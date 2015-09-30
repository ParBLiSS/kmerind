/**
 * @file    uniform_omp_runner.hpp
 * @ingroup taskrunner
 * @author  Tony Pan <tpan7@gatech.edu>
 * @brief   an OMP task execution engine where all threads execute the same list of tasks
 * @details allows adding tasks during execution
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef UNIFORM_OMP_RUNNER_HPP_
#define UNIFORM_OMP_RUNNER_HPP_

#include "omp.h"
#include <vector>
#include <cassert>

#include "bliss-config.hpp"
#include "taskrunner/runner.hpp"
//#include "concurrent/mutexlock_queue.hpp"
#include "concurrent/lockfree_queue.hpp"


#define LockType bliss::concurrent::LockType::LOCKFREE

namespace bliss
{
namespace concurrent
{

/**
 * @class      bliss::concurrent::UniformOMPRunner
 * @brief       runner that performs the same task(s) for each thread.
 * @details     for multiple tasks, either tasks can be added directly
 * 				or a sequential runner can be used.
 */
class UniformOMPRunner : public Runner
{
  protected:



    // since each thread executes the same tasks we need an array of queues
	// using threadsafe queue because task addition may be multithreaded.
    std::vector<bliss::concurrent::ThreadSafeQueue<std::shared_ptr<Runnable> , LockType> > qs;

    const int nThreads;

  public:

    /**
     * @brief Creates a new UniformOMPRunner instance with the given number of
     *        threads
     *
     * @param nThreads The number of OpenMP threads to use.
     */
    UniformOMPRunner(const int num_threads) : Runner(), qs(num_threads), nThreads(num_threads) {
#ifdef USE_OPENMP
      if (omp_get_nested() <= 0) omp_set_nested(1);
      if (omp_get_dynamic() <= 0) omp_set_dynamic(1);
#else
      static_assert(false, "UniformOMPRunner used although compilation is not set to use OpenMP");
#endif

      qs.clear();
      for (int i = 0; i < nThreads; ++i) {
        qs.push_back(bliss::concurrent::ThreadSafeQueue<std::shared_ptr<Runnable>, LockType >());
      }
    }

    /**
     * @brief The destructor of this class.
     */
    virtual ~UniformOMPRunner() {};

    /**
     * @brief Runs all tasks until queue is empty.  each thread performs the same tasks.
     */
    void operator()() {
      size_t proc = 0;
#pragma omp parallel num_threads(nThreads) default(none) reduction(+: proc)
      {
        // iterator for list is valid during append (no insertion in middle, no deletion).
        // no deletion since this is shared.
        bliss::concurrent::ThreadSafeQueue<std::shared_ptr<Runnable>, LockType >& q = qs[omp_get_thread_num()];
        while (q.canPop())
        {
          auto v = q.waitAndPop();
          if (v.first) {
            (v.second)->operator()();
            ++proc;
          }
        }
      }
      INFOF("Uniform Runner completed %lu tasks.", proc);

    }

    /**
     * @brief Adds a task to this Runner.
     * @details  each thread has its own queue.  adding a task adds to all queues.
     * @param t
     */
    virtual bool addTask(std::shared_ptr<Runnable> &&t) {
      bool out = true;
      int id;
      int i = 0;
      // make copy and insert
      for (; i < nThreads-1; ++i) {
        id = (omp_get_thread_num() + i) % nThreads;
        bliss::concurrent::ThreadSafeQueue<std::shared_ptr<Runnable>, LockType >& q = qs[id];
        auto result = q.waitAndPush(std::shared_ptr<Runnable>(t));

        out &= result.first;
      }
      // insert the real thing.
      id = (omp_get_thread_num() + nThreads-1) % nThreads;
      bliss::concurrent::ThreadSafeQueue<std::shared_ptr<Runnable>, LockType >& q = qs[id];
      auto result = q.waitAndPush(std::forward<std::shared_ptr<Runnable> >(t));
      out &= result.first;

      DEBUGF("add to Uniform runner.  size %lu, disabled %s\n", q.getSize(), (q.canPush() ? "n" : "y"));

      return out;

    }

    /// count the number of pending tasks
    virtual size_t getTaskCount() {
      size_t max = 0;
      int id;
      for (int i = 0; i < nThreads; ++i) {
        id = (omp_get_thread_num() + i) % nThreads;
        max = std::max(max, qs[id].getSize());
      }
      return max;
    }

    /// check to see if new tasks can be added to the queue
    virtual bool isAddDisabled() {
      bool canPush = true;
      int id;
      for (int i = 0; i < nThreads; ++i) {
        id = (omp_get_thread_num() + i) % nThreads;
        canPush &= qs[id].canPush();
      }
      return !canPush;
    }

    /// flush currently queued tasks and disallow further new tasks
    virtual void disableAdd() {
      int id;
      for (int i = 0; i < nThreads; ++i) {
        id = (omp_get_thread_num() + i) % nThreads;
        qs[id].disablePush();
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
