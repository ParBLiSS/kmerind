/**
 * @file    runnable.hpp
 * @ingroup
 * @author  tpan
 * @brief
 * @details
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
    Runnable() {};

    virtual ~Runnable() {};

    virtual void operator()() = 0;
};

} /* namespace concurrent */
} /* namespace bliss */

#endif /* RUNNABLE_HPP_ */
