/**
 * @file		runner.hpp
 * @ingroup
 * @author	tpan
 * @brief
 * @details
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef RUNNER_HPP_
#define RUNNER_HPP_

namespace bliss
{
  namespace concurrent
  {

    /**
     * @class			bliss::concurrent::Runner
     * @brief     abstract base class for executing some tasks
     * @details   uses data Src, data Sink, and Task for computation.
     *
     *            multiple sources can be represented as a single one with some non-conflicting mapping
     *              One-to-One (block or cyclic partition), demand driven (locking), or random (locking)
     *              many to many is a mapping function (hash) followed by demand driven (locking)
     *            task may be independent of source (constant source or zero source)
     *            multiple sink assignment (hash).  single sink requires locking.
     *
     *            there is always a sink.
     */
    template <typename Task, typename Sink, typename TaskToSinkMap, typename Source = void, typename TaskToSourceMap = void >
    class Runner
    {
      public:
        Runner();
        virtual ~Runner();

        void run();
    };



  } /* namespace concurrent */
} /* namespace bliss */

#endif /* RUNNER_HPP_ */
