/**
 * @file    sequential_runner.hpp
 * @ingroup taskrunner
 * @author  Tony Pan <tpan7@gatech.edu>
 * @brief	a sequential task engine that executes all tasks in a queue.
 * @details
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef SEQUENTIAL_RUNNER_HPP_
#define SEQUENTIAL_RUNNER_HPP_


#include <cassert>

//#include "concurrent/mutexlock_queue.hpp"
#include "concurrent/lockfree_queue.hpp"

#include "taskrunner/runner.hpp"

#define LockType bliss::concurrent::LockType::LOCKFREE


namespace bliss
{
namespace concurrent
{

/**
 * @class     bliss::concurrent::SequentialRunner
 * @brief     runs a sequence of tasks
 * @details   can have multiple tasks, run in order of adding.
 *            MORE than just a compound task, since this does NOT just have a single input and output.
 *			  tasks can be added during execution.
 */
class SequentialRunner : public Runner
{
  protected:
    // using thread safe queue because other threads could be calling addTask.
    bliss::concurrent::ThreadSafeQueue<std::shared_ptr<Runnable> , LockType> q;

  public:
    /// constructor
    SequentialRunner() : Runner() {};

    /// default destructor
    virtual ~SequentialRunner() {};

    /**
	 * @brief executes the tasks
	 * @details   executes all tasks in a queue until it is empty or disbled.
	 */
    void operator()()
    {
      size_t proc = 0;
      while (q.canPop())
      {
        auto v = q.waitAndPop();
        if (v.first) {
          (v.second)->operator()();
          ++proc;
        }
      }
      INFOF("Sequential runner completed %lu tasks.", proc);

    }
    /// add a new task to queue
    virtual bool addTask(std::shared_ptr<Runnable> &&t)
    {
      DEBUGF("add to Sequential runner.  size %lu, disabled %s", q.getSize(), (q.canPush() ? "n" : "y"));
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

    /// allows synchronization.  does not apply to sequential runner. (not used right now)
    virtual void synchronize()
    {
    }

};

} /* namespace concurrent */
} /* namespace bliss */

#endif /* SEQUENTIAL_RUNNER_HPP_ */
