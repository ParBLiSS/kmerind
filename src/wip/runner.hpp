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

#include <atomic>
#include <memory>

#include "wip/runnable.hpp"

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

  public:
    Runner() {};
    virtual ~Runner() {
    };

    // function to run
    virtual void operator()() = 0;

    virtual bool addTask(Runnable* _t) = 0;

    /// flush currently queued tasks and disallow further new tasks
    virtual void disableAdd() = 0;

    virtual void synchronize() = 0;


};



} /* namespace concurrent */
} /* namespace bliss */

#endif /* RUNNER_HPP_ */
