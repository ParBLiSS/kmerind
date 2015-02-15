/**
 * @file    task.hpp
 * @ingroup taskrunner
 * @author  Tony Pan
 * @brief	API interface for a task.  essentially a functor
 * @details
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef TASK_HPP_
#define TASK_HPP_

#include "taskrunner/runnable.hpp"

namespace bliss
{
namespace concurrent
{

/**
 * @class     bliss::concurrent::Task
 * @brief	  a conceptual unit of work, or a task.
 * @details   Defined in the form of a functor.
 * 			  It is the user's responsibility to define the input and output
 * 			  and the computation.
 *
 * 			  Derived from Runnable so that the task engine can provide a
 * 			  standard task management interface (add, run, etc).
 *
 */
  class Task : public Runnable
  {
    public:
	  /// default construtor
	  Task() {};
	  /// default destructor
      virtual ~Task() {};
      /// interface function for running a task.  essentially a functor.
      virtual void operator()() = 0;
  };

} /* namespace concurrent */
} /* namespace bliss */

#endif /* TASK_HPP_ */
