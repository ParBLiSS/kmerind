/**
 * @file		sequential_runner.hpp
 * @ingroup
 * @author	tpan
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
#include <vector>

#include "config.hpp"
#include "concurrent/runner.hpp"

namespace bliss
{
  namespace concurrent
  {

    /**
     * @class			bliss::concurrent::SequentialRunner
     * @brief
     * @details   can have multiple tasks, run in order of adding.
     *
     */
    class SequentialRunner : public Runner
    {
      public:
        SequentialRunner() : Runner() {};
        virtual ~SequentialRunner() {};

        virtual void addTask(Runnable &t) {
          q.push_back(t);
        }

        virtual void run() {
          for (Runnable t : q) {
            t.run();
          }
        }

        virtual void synchronize() {
        }
      protected:
        std::vector<Runnable> q;

    };

  } /* namespace concurrent */
} /* namespace bliss */

#endif /* SEQUENTIAL_RUNNER_HPP_ */
