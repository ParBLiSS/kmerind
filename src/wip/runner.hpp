/**
 * @file    runner.hpp
 * @ingroup
 * @author  tpan
 * @brief
 * @details
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef RUNNER_HPP_
#define RUNNER_HPP_

#include "wip/runnable.hpp"
#include <atomic>
#include <memory>

namespace bliss
{
namespace concurrent
{

/**
 * @class     bliss::concurrent::Runner
 * @brief     abstract base class for executing some tasks
 * @details   uses data Src, data Sink, and Runnable for computation.
 *
 *            multiple sources can be represented as a single one with some
 *            non-conflicting mapping
 *            - One-to-One: (block or cyclic partition), demand driven
 *                          (locking), or random (locking)
 *            - many to many: a mapping function (hash) followed by demand
 *                            driven (locking)
 *            task may be independent of source (constant source or zero source)
 *            multiple sink assignment (hash).  single sink requires locking.
 *
 *            there is always a sink.
 *
 *            coarse grain partition.
 *
 *            inherit from Runnable to allow nesting.
 *
 * TODO:           not subclassed yet: streaming runner.
 */

class Runner : public Runnable
{

  protected:
    std::atomic<bool> blocked;

  public:
    Runner() : blocked(false) {};
    virtual ~Runner() {
      blocked = true;
    };

    // function to run
    virtual void operator()() = 0;

    virtual bool addTask(std::unique_ptr<Runnable> &&_t) = 0;

    /// flush currently queued tasks and disallow further new tasks
    void blockAdd() {
      blocked.store(true, std::memory_order_relaxed);
    };

    /// flush currently queued tasks and disallow further new tasks
    void unblockAdd() {
      blocked.store(false, std::memory_order_relaxed);
    };

    virtual void synchronize() = 0;


};



} /* namespace concurrent */
} /* namespace bliss */

#endif /* RUNNER_HPP_ */
