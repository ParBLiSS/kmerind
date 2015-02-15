/**
 * bliss-Release@bliss copyright 2014 Georgia Institute of Technology
 *
 * @file    thread_pool.hpp
 * @ingroup wip
 *
 * @date:  Mar 20, 2014
 * @author: Tony Pan <tpan7@gatech.edu>
 * @description
 *
 * inspired by https://github.com/greyfade/workqueue/blob/master/threadpool.hpp
 * modifications:  none yet.
 */
#ifndef THREADPOOL_H_
#define THREADPOOL_H_

#include <future>
#include <thread>
#include <deque>
#include <vector>
#include <utility>
#include <chrono>
#include <functional>
#include <type_traits>


namespace bliss
{
  namespace concurrent
  {

    /*
    * This file is licensed under the zlib/libpng license, included in this
    * distribution in the file COPYING.
    */

    class thread_pool {
        std::vector<std::thread> pool;
        bool stop;

        std::mutex access;
        std::condition_variable cond;
        std::deque<std::function<void()>> tasks;

      public:
        explicit thread_pool(int nr = 1) : stop(false) {
          while(--nr > 0) {
            add_worker();
          }
        }
        ~thread_pool() {
          stop = true;
          for(std::thread &t : pool) {
            t.join();
          }
          pool.clear();
        }

        template<class Rt>
        auto add_task(std::packaged_task<Rt()>& pt) -> std::future<Rt> {
            std::unique_lock<std::mutex> lock(access);

            auto ret = pt.get_future();
            tasks.push_back([&pt]{pt();});

            cond.notify_one();

            return ret;
        }

      private:
        void add_worker() {
          std::thread t([this]() {
            while(!stop) {
              std::function<void()> task;
              {
                std::unique_lock<std::mutex> lock(access);
                if(tasks.empty()) {
                  cond.wait_for(lock, std::chrono::duration<int, std::milli>(5));
                  continue;
                }
                task = std::move(tasks.front());
                tasks.pop_front();
              }
              task();
            }
          });
          pool.push_back(std::move(t));
        }
    };


  } /* namespace concurrent */
} /* namespace bliss */

#endif /* THREADPOOL_H_ */
