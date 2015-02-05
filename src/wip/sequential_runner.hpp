/**
 * @file    sequential_runner.hpp
 * @ingroup
 * @author  tpan
 * @brief
 * @details
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef SEQUENTIAL_RUNNER_HPP_
#define SEQUENTIAL_RUNNER_HPP_


#include <cassert>

#include "concurrent/lockfree_queue.hpp"

#include "wip/runner.hpp"

namespace bliss
{
namespace concurrent
{

/**
 * @class     bliss::concurrent::SequentialRunner
 * @brief     runs a sequence of tasks
 * @details   can have multiple tasks, run in order of adding.
 *            MORE than just a compound task, since this does NOT just have a single input and output.
 *
 */
class SequentialRunner : public Runner
{
  protected:
    // using thread safe queue because other threads could be calling addTask.
    bliss::concurrent::ThreadSafeQueue<std::shared_ptr<Runnable> > q;

  public:
    SequentialRunner() : Runner() {};

    virtual ~SequentialRunner() {};

    void operator()()
    {
      while (q.canPop())
      {
        auto v = std::move(q.waitAndPop());
        if (v.first) {
          (v.second)->operator()();
//          counter2.fetch_add(1);
        }
      }
      printf("Sequential runner completed.\n");

    }

    virtual bool addTask(std::shared_ptr<Runnable> &&t)
    {
//      counter.fetch_add(1);

      auto result = q.waitAndPush(std::forward<std::shared_ptr<Runnable> >(t));
//      printf("add to Sequential runner.  size %lu, disabled %s\n", q.getSize(), (q.canPush() ? "n" : "y"));

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
    }

};

} /* namespace concurrent */
} /* namespace bliss */

#endif /* SEQUENTIAL_RUNNER_HPP_ */
