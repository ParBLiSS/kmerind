/**
 * @file    runner.hpp
 * @ingroup taskrunner
 * @author  Tony Pan
 * @brief	an abstract class to represent a task engine, "runner"
 * @details derived from Runnable so that Runner can be treated as tasks as well.
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef RUNNER_HPP_
#define RUNNER_HPP_

#include <atomic>
#include <memory>

#include "taskrunner/runnable.hpp"

namespace bliss
{
namespace concurrent
{

/**
 * @class     bliss::concurrent::Runner
 * @brief     abstract base class for executing some tasks
 * @details   represents coarse grain partition.
 *
 *            inherit from Runnable to allow nesting of Runners as tasks
 *
 *            specific task implementation allows the following flow of
 *            control patterns to be constructed:
 *            1. execute once:  task perform computation and terminate
 *            2. looping:  task reinsert into the same task queue at end of computation
 *            3. dependent computation:  task instantiate additional task of different type
 *                  and reinsert into same task queue
 *            4. conditional computation:  task generate additional task depending
 *                  on some state, and reinsert into same task queue
 *            5. signaling: task generate additional task and insert into a different queue,
 *                  e.g. master thread queue.
 *            6. schedule on specific resource: task generate additional tasks and
 *                  insert into queue tied to particular resource (e.g. io, accelerator, etc)
 *
 *
 *            supported concurrency patterns:
 *            1. sequential:  single thread executes all tasks in queue in order
 *            2. concurrent, uniform execution:  multiple threads, each execute all
 *                tasks in queue.  implementation provided as uniform_omp_runner, or
 *                alternatively can be done via sequential runner nested in personalized_omp_runner
 *            3. concurrent, personalized execution:  multiple threads, each execute
 *                one task in queue (no dynamic scheduling).  tasks may be heterogeneous.  implementation provided
 *                via personalized_omp_runner.  Note that once the execution starts, no
 *                additional tasks may be added, which differs from other implementations.
 *            4. concurrent, demand driven execution:  multiple threads.  tasks are
 *                scheduled into the omp task queue, and executed by available free threads.
 *
 *
 *
 * @note      OpenMP has its own omp task queue.  The omp task queue is populated
 *            by a single thread in the parallel section.  this precludes
 *            1. defining tasks or a sequence of tasks before entering the
 *                parallel block (minor issue)
 *              1a. user has to define all possible compute tasks, as well as a master
 *                  task that populates the task queue.
 *            2. switching to pthread or other threading backend that does not
 *                provide a native task queue and thread pool (design issue).
 *            3. supporting all concurrency patterns:
 *                a. sequential can use omp task construct with 2 threads, 1 for creating task and 1 to execute.
 *                b. uniform parallel execution requires p copies of the task, omp for with static schedule clause
 *                c. personalized parallel execution requires omp for with static schedule clause
 *                d. dynamic execution can be represented as omp task queue.
 *
 *            for 1a, if the master task also has computation, such as sending/receiving
 *            over communication or io, master task becomes complex and potentially
 *            cannot dynamically schedule its own subtasks on demand.
 *
 *            for 3a, we need 2 threads instead of 1.
 *            for 3b and 3c, since we need for loop, we require a list/vector and
 *              cannot dynamically insert anyways.
 *
 *            if omp task queue were to be used directly without a middleware level
 *            task queue, then we need a master task, and would still need an api
 *            like the one defined here to encapsulate the parallel blocks.
 *            at the same time, dynamically adding tasks becomes hard or impossible
 *            which would need to be done by master thread anyway, not by application
 *            code directly
 *
 *            summary:  A middleware task queue simplifies logic and enable certain
 *              concurrency patterns, when compared to using omp task queue directly.
 *
 */
class Runner : public Runnable
{
  public:
	/// default constructor
    Runner() {};
    /// default destructor
    virtual ~Runner() {};

    /// interface function, for executes the tasks
    virtual void operator()() = 0;

    /// interface for adding new tasks to queue
    virtual bool addTask(std::shared_ptr<Runnable> &&_t) = 0;

    /// interface function for getting the number of pending tasks
    virtual size_t getTaskCount() = 0;

    /// interface function for checking if queue is blocked from future task addition
    virtual bool isAddDisabled() = 0;

    /// flush currently queued tasks and disallow further new tasks
    virtual void disableAdd() = 0;

    /// interface function for synchronizing between runners.
    virtual void synchronize() = 0;


};



} /* namespace concurrent */
} /* namespace bliss */

#endif /* RUNNER_HPP_ */
