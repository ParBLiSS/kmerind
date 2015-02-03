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
    bliss::concurrent::ThreadSafeQueue<Runnable* > q;

  public:
    SequentialRunner() : Runner() {};

    virtual ~SequentialRunner() {};

    void operator()()
    {
      // should q be allowed to be changed?

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

    virtual bool addTask(Runnable* t)
    {
      auto result = q.waitAndPush(std::forward<Runnable* >(t));

      return result.first;
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
