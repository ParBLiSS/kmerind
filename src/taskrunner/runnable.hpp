/**
 * @file    runnable.hpp
 * @ingroup taskrunner
 * @author  Tony Pan
 * @brief	abstract base class for Runnable objects, including tasks and task engine
 * @details	common base class allows task engine to be treated as tasks.
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef RUNNABLE_HPP_
#define RUNNABLE_HPP_

#include <memory>
#include "utils/logging.h"

namespace bliss
{
namespace concurrent
{

/**
 * @class   bliss::concurrent::Runnable
 * @brief   not instantiable. An abstract concept of a runnable task.
 *          could be multithreaded, mpi based, or a single threaded task.
 * @details subclasses include runner and is subclasses, and user defined
 *          tasks.
 *
 *         	This interface allows Runners to be treated as tasks as well, so
 *         	nesting Runners is possible.
 *
 *          since we allow mix of different types of Runnables (tasks),
 *          we cannot use CRTP and static polymorphism, which requires
 *          all types to be known at compile time.
 *
 *          Use dynamic polymorphism and restrict task to be coarse grain
 *          (NOT inside inner loop).
 */
class Runnable : public std::enable_shared_from_this<Runnable>
{
  public:
	/// constructor
    Runnable() {};

    /// destructor
    virtual ~Runnable() {};

    /// abstract "run" operator.
    virtual void operator()() = 0;
};

} /* namespace concurrent */
} /* namespace bliss */

#endif /* RUNNABLE_HPP_ */
