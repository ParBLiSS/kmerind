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
#include <deque>

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
    std::deque<std::unique_ptr<Runnable> > q;

  public:
    SequentialRunner() : Runner() {};

    virtual ~SequentialRunner() {

    };

    void operator()()
    {
      while (!q.empty())
      {
        q.front()->operator()();
        q.pop_front();
      }
    }

    virtual bool addTask(std::unique_ptr<Runnable> &&t)
    {
      bool b = !blocked;
      if (b) {
        q.push_back(std::forward<std::unique_ptr<Runnable> >(t));
      }
      return b;
    }

    virtual void synchronize()
    {
    }

};

} /* namespace concurrent */
} /* namespace bliss */

#endif /* SEQUENTIAL_RUNNER_HPP_ */
