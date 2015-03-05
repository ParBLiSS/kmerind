/**
 * @file		test_runners.cpp
 * @ingroup
 * @author	Tony Pan <tpan7@gatech.edu>
 * @brief
 * @details
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */

#include <iostream>
#include <unistd.h>  // sleep
#include <chrono>

//#include "concurrent/mpi_runner.hpp"
#include "taskrunner/uniform_omp_runner.hpp"
#include "taskrunner/personalized_omp_runner.hpp"
#include "taskrunner/dynamic_omp_runner.hpp"
#include "taskrunner/sequential_runner.hpp"
#include "taskrunner/task.hpp"

static const int iter = 100000;
std::atomic<int> cc;
std::atomic<int> cc2;


class Test : public bliss::concurrent::Task
{
  public:
    Test(std::string msg) : bliss::concurrent::Task(), content(msg) {};
    virtual ~Test() {};

    virtual void operator()() {
      //INFOF("tid %d: %s\n", omp_get_thread_num(), content.c_str());
    }

  protected:
    std::string content;

};


class Test2 : public bliss::concurrent::Task
{
  protected:
    std::string content;
    std::shared_ptr<bliss::concurrent::Runner> parent;
  public:
    Test2(std::string msg,
          const std::shared_ptr<bliss::concurrent::Runner>& _parent) :
            bliss::concurrent::Task(), content(msg), parent(_parent) {};
    virtual ~Test2() {};

    virtual void operator()() {
      //INFOF("tid %d: %s\n", omp_get_thread_num(), content.c_str());
      parent->addTask(this->shared_from_this());
//      if (parent.addTask(this))
//        INFOF("tid %d: %s reinserted\n", omp_get_thread_num(), content.c_str());
    }


};

class Test3 : public bliss::concurrent::Task
{
  protected:
    std::string content;
    std::shared_ptr<bliss::concurrent::Runner> other;
  public:
    Test3(std::string msg, const std::shared_ptr<bliss::concurrent::Runner>& _other) : bliss::concurrent::Task(), content(msg), other(_other) {};
    virtual ~Test3() {};

    virtual void operator()() {
      INFOF("tid %d: %s start\n", omp_get_thread_num(), content.c_str());

      usleep(1000);

      INFOF("tid %d: %s stop\n", omp_get_thread_num(), content.c_str());
      other->disableAdd();
      INFOF("tid %d: %s stopped\n", omp_get_thread_num(), content.c_str());
    }

};



class Test4 : public bliss::concurrent::Task
{
  protected:
    std::shared_ptr<bliss::concurrent::Runner> work_r;

    std::string content;
  public:
    Test4(std::string msg,
           const std::shared_ptr<bliss::concurrent::Runner>& work) :
             work_r(work), content(msg) {};

    virtual ~Test4() {};

    virtual void operator()() {
      auto id = cc.fetch_add(1, std::memory_order_relaxed);

      if (id < iter) {
        //INFOF("tid %d: %s, count %d\n", omp_get_thread_num(), content.c_str(), i);
        //INFOF("%dR%d ", omp_get_thread_num(), i);

        work_r->addTask(this->shared_from_this());
      } else {
        work_r->disableAdd();
      }
    }

};
template<class Receiver>
class Sender : public bliss::concurrent::Task
{
  protected:
    std::shared_ptr<bliss::concurrent::Runner> work_r;
    std::shared_ptr<bliss::concurrent::Runner> comm_r;
    std::string content;

    int id;


  public:
    Sender(std::string msg,
           const std::shared_ptr<bliss::concurrent::Runner>& comm,
           const std::shared_ptr<bliss::concurrent::Runner>& work, int _id) :
             work_r(work), comm_r(comm), content(msg), id(_id) {};

    virtual ~Sender() {};

    virtual void operator()() {
      //INFOF("tid %d: %s\n", omp_get_thread_num(), content.c_str());
      //INFOF("%dS ", omp_get_thread_num());
      _mm_pause();
//      INFOF("tid %d: %s add recv %d to comm runner size %lu, disabled %s \n", omp_get_thread_num(), content.c_str(), id, comm_r->getTaskCount(), (comm_r->isAddDisabled() ? "y" : "n"));
      comm_r->addTask(std::move(std::shared_ptr<Runnable>(new Receiver("recv", comm_r, work_r, id))));

    }

};

template<class SENDER>
class Node : public bliss::concurrent::Task
{
  protected:
    std::shared_ptr<bliss::concurrent::Runner> work_r;
    std::shared_ptr<bliss::concurrent::Runner> comm_r;
    std::string content;
    int id;

  public:
    Node(std::string msg,
         const std::shared_ptr<bliss::concurrent::Runner>& comm,
         const std::shared_ptr<bliss::concurrent::Runner>& work, int _id) :
           work_r(work), comm_r(comm), content(msg), id(_id) {};

    virtual ~Node() {};

    virtual void operator()() {

      _mm_pause();
      _mm_pause();
      _mm_pause();
      _mm_pause();
      _mm_pause();
      _mm_pause();
      _mm_pause();
      _mm_pause();
      _mm_pause();
      _mm_pause();

//      INFOF("tid %d: %s add send %d to comm runner size %lu, disabled %s \n", omp_get_thread_num(), content.c_str(), id, comm_r->getTaskCount(), (comm_r->isAddDisabled() ? "y" : "n"));
      comm_r->addTask(std::move(std::shared_ptr<Runnable>(new SENDER("sender", comm_r, work_r, id))));

    }

};

template<class SENDER>
class Source : public bliss::concurrent::Task
{
  protected:
    std::shared_ptr<bliss::concurrent::Runner> work_r;
    std::shared_ptr<bliss::concurrent::Runner> comm_r;
    std::string content;
  public:
    Source(std::string msg,
           const std::shared_ptr<bliss::concurrent::Runner>& comm,
           const std::shared_ptr<bliss::concurrent::Runner>& work) :
             work_r(work), comm_r(comm), content(msg) {};

    virtual ~Source() {};

    virtual void operator()() {
      auto id = cc.fetch_add(1, std::memory_order_relaxed);

        //INFOF("tid %d: %s, count %d\n", omp_get_thread_num(), content.c_str(), i);
        //INFOF("%dR%d ", omp_get_thread_num(), i);
        _mm_pause();
        _mm_pause();
        _mm_pause();

//        INFOF("tid %d: %s add send %d to comm runner size %lu, disabled %s \n", omp_get_thread_num(), content.c_str(), id, comm_r->getTaskCount(), (comm_r->isAddDisabled() ? "y" : "n"));
//        INFOF("tid %d: %s add self %d to work runner size %lu, disabled %s \n", omp_get_thread_num(), content.c_str(), id, work_r->getTaskCount(), (work_r->isAddDisabled() ? "y" : "n"));
        comm_r->addTask(std::move(std::shared_ptr<Runnable>(new SENDER("sender", comm_r, work_r, id))));

      if (id+1 < iter) {
        work_r->addTask(this->shared_from_this());
      } else {

        if (id %1000 == 0) INFOF("tid %d: %s, count %d\n", omp_get_thread_num(), content.c_str(), id);
      }
    }

};

class Sink : public bliss::concurrent::Task
{
  protected:
    std::shared_ptr<bliss::concurrent::Runner> work_r;
    std::shared_ptr<bliss::concurrent::Runner> comm_r;
    std::string content;

    int id;

  public:
    Sink(std::string msg,
         const std::shared_ptr<bliss::concurrent::Runner>& comm,
         const std::shared_ptr<bliss::concurrent::Runner>& work, int _id) :
           work_r(work), comm_r(comm), content(msg), id(_id) {};

    virtual ~Sink() {};

    virtual void operator()() {
      auto i = cc2.fetch_add(1, std::memory_order_relaxed);

      if (i+1 >= iter) {
        INFOF("tid %d: disabled runners at %d, item %d, iter %d\n", omp_get_thread_num(), id, i, iter);
        work_r->disableAdd();
        comm_r->disableAdd();
      } else {
        //INFOF("%dW ", omp_get_thread_num());
        _mm_pause();
        _mm_pause();
        _mm_pause();
        _mm_pause();
        _mm_pause();
        _mm_pause();
        _mm_pause();
        _mm_pause();
        _mm_pause();
        _mm_pause();

//        INFOF("tid %d: sink id %d, item %d, iter %d, work runner size %lu, disabled %s\n", omp_get_thread_num(), id ,i, iter, work_r->getTaskCount(), (work_r->isAddDisabled() ? "y" : "n"));
      }

    }

};

template<class WRITER>
class Receiver : public bliss::concurrent::Task
{
  protected:
    std::shared_ptr<bliss::concurrent::Runner> work_r;
    std::shared_ptr<bliss::concurrent::Runner> comm_r;
    std::string content;
    int id;

  public:
    Receiver(std::string msg,
             const std::shared_ptr<bliss::concurrent::Runner>& comm,
             const std::shared_ptr<bliss::concurrent::Runner>& work, int _id) :
               work_r(work), comm_r(comm),  content(msg), id(_id) {};

    virtual ~Receiver() {};

    virtual void operator()() {
        //INFOF("%dC%d ", omp_get_thread_num(), i);
        _mm_pause();
        _mm_pause();

//        INFOF("tid %d: %s add write %d to work runner size %lu, disabled %s \n", omp_get_thread_num(), content.c_str(), id, work_r->getTaskCount(), (work_r->isAddDisabled() ? "y" : "n"));

        work_r->addTask(std::move(std::shared_ptr<Runnable>(new WRITER("write", comm_r, work_r, id))));

    }

};



int main(int argc, char** argv) {

  std::chrono::high_resolution_clock::time_point t1, t2;
  std::chrono::duration<double> time_span;


  {
    t1 = std::chrono::high_resolution_clock::now();

    std::shared_ptr<bliss::concurrent::Runner> sr(new bliss::concurrent::SequentialRunner());
    //// test sequential.
    INFOF("serial runner\n");
    sr->addTask(std::move(std::shared_ptr<bliss::concurrent::Runnable>(new Test("task 1"))));
    sr->addTask(std::move(std::shared_ptr<bliss::concurrent::Runnable>(new Test("task 2"))));
    sr->addTask(std::move(std::shared_ptr<bliss::concurrent::Runnable>(new Test("task 3"))));
    sr->addTask(std::move(std::shared_ptr<bliss::concurrent::Runnable>(new Test("task 4"))));
    sr->disableAdd();

    t2 = std::chrono::high_resolution_clock::now();
    time_span =
        std::chrono::duration_cast<std::chrono::duration<double>>(
            t2 - t1);
    INFO( "setup: time " << time_span.count() );

    t1 = std::chrono::high_resolution_clock::now();

    sr->operator()();

    t2 = std::chrono::high_resolution_clock::now();
    time_span =
        std::chrono::duration_cast<std::chrono::duration<double>>(
            t2 - t1);
    INFO( "run: time " << time_span.count() );


  }


  {
    t1 = std::chrono::high_resolution_clock::now();

    cc.store(0, std::memory_order_relaxed);

    std::shared_ptr<bliss::concurrent::Runner> sr(new bliss::concurrent::SequentialRunner());
    //// test sequential.
    INFOF("serial runner\n");
    sr->addTask(std::move(std::shared_ptr<bliss::concurrent::Runnable>(new Test4("task 1", sr))));


    t2 = std::chrono::high_resolution_clock::now();
    time_span =
        std::chrono::duration_cast<std::chrono::duration<double>>(
            t2 - t1);
    INFO( "setup: time " << time_span.count() );

    t1 = std::chrono::high_resolution_clock::now();

    sr->operator()();

    t2 = std::chrono::high_resolution_clock::now();
    time_span =
        std::chrono::duration_cast<std::chrono::duration<double>>(
            t2 - t1);
    INFO( "run: time " << time_span.count() );


  }



  {
    t1 = std::chrono::high_resolution_clock::now();

    std::shared_ptr<bliss::concurrent::Runner> sr(new bliss::concurrent::PersonalizedOMPRunner(3));
    //// test sequential.
    INFOF("personalized OMP runner\n");
    sr->addTask(std::move(std::shared_ptr<bliss::concurrent::Runnable>(new Test("task 1"))));
    sr->addTask(std::move(std::shared_ptr<bliss::concurrent::Runnable>(new Test("task 2"))));
    sr->addTask(std::move(std::shared_ptr<bliss::concurrent::Runnable>(new Test("task 3"))));
    sr->addTask(std::move(std::shared_ptr<bliss::concurrent::Runnable>(new Test("task 4"))));
    sr->disableAdd();

    t2 = std::chrono::high_resolution_clock::now();
    time_span =
        std::chrono::duration_cast<std::chrono::duration<double>>(
            t2 - t1);
    INFO( "setup: time " << time_span.count() );

    t1 = std::chrono::high_resolution_clock::now();

    sr->operator()();

    t2 = std::chrono::high_resolution_clock::now();
    time_span =
        std::chrono::duration_cast<std::chrono::duration<double>>(
            t2 - t1);
    INFO( "run: time " << time_span.count() );


  }


  {
    t1 = std::chrono::high_resolution_clock::now();

    cc.store(0, std::memory_order_relaxed);

    std::shared_ptr<bliss::concurrent::Runner> sr(new bliss::concurrent::PersonalizedOMPRunner(4));
    //// test sequential.
    INFOF("personalized OMP runner\n");
    sr->addTask(std::move(std::shared_ptr<bliss::concurrent::Runnable>(new Test4("task 1", sr))));
    sr->addTask(std::move(std::shared_ptr<bliss::concurrent::Runnable>(new Test4("task 2", sr))));
    sr->addTask(std::move(std::shared_ptr<bliss::concurrent::Runnable>(new Test4("task 3", sr))));
    sr->addTask(std::move(std::shared_ptr<bliss::concurrent::Runnable>(new Test4("task 4", sr))));

    t2 = std::chrono::high_resolution_clock::now();
    time_span =
        std::chrono::duration_cast<std::chrono::duration<double>>(
            t2 - t1);
    INFO( "setup: time " << time_span.count() );

    t1 = std::chrono::high_resolution_clock::now();

    sr->operator()();

    t2 = std::chrono::high_resolution_clock::now();
    time_span =
        std::chrono::duration_cast<std::chrono::duration<double>>(
            t2 - t1);
    INFO( "run: time " << time_span.count() );


  }



  {
    t1 = std::chrono::high_resolution_clock::now();

    std::shared_ptr<bliss::concurrent::Runner> sr(new bliss::concurrent::UniformOMPRunner(3));
    //// test sequential.
    INFOF("uniform OMP runner\n");
    sr->addTask(std::move(std::shared_ptr<bliss::concurrent::Runnable>(new Test("task 1"))));
    sr->addTask(std::move(std::shared_ptr<bliss::concurrent::Runnable>(new Test("task 2"))));
    sr->addTask(std::move(std::shared_ptr<bliss::concurrent::Runnable>(new Test("task 3"))));
    sr->addTask(std::move(std::shared_ptr<bliss::concurrent::Runnable>(new Test("task 4"))));
    sr->disableAdd();

    t2 = std::chrono::high_resolution_clock::now();
    time_span =
        std::chrono::duration_cast<std::chrono::duration<double>>(
            t2 - t1);
    INFO( "setup: time " << time_span.count() );

    t1 = std::chrono::high_resolution_clock::now();

    sr->operator()();

    t2 = std::chrono::high_resolution_clock::now();
    time_span =
        std::chrono::duration_cast<std::chrono::duration<double>>(
            t2 - t1);
    INFO( "run: time " << time_span.count() );


  }

  {
    t1 = std::chrono::high_resolution_clock::now();

    cc.store(0, std::memory_order_relaxed);

    std::shared_ptr<bliss::concurrent::Runner> sr(new bliss::concurrent::UniformOMPRunner(4));
    //// test sequential.
    INFOF("uniform OMP runner\n");
    sr->addTask(std::move(std::shared_ptr<bliss::concurrent::Runnable>(new Test4("task 1", sr))));
    sr->addTask(std::move(std::shared_ptr<bliss::concurrent::Runnable>(new Test4("task 2", sr))));
    sr->addTask(std::move(std::shared_ptr<bliss::concurrent::Runnable>(new Test4("task 3", sr))));
    sr->addTask(std::move(std::shared_ptr<bliss::concurrent::Runnable>(new Test4("task 4", sr))));

    t2 = std::chrono::high_resolution_clock::now();
    time_span =
        std::chrono::duration_cast<std::chrono::duration<double>>(
            t2 - t1);
    INFO( "setup: time " << time_span.count() );

    t1 = std::chrono::high_resolution_clock::now();

    sr->operator()();

    t2 = std::chrono::high_resolution_clock::now();
    time_span =
        std::chrono::duration_cast<std::chrono::duration<double>>(
            t2 - t1);
    INFO( "run: time " << time_span.count() );


  }

  {
    t1 = std::chrono::high_resolution_clock::now();

    std::shared_ptr<bliss::concurrent::Runner> sr(new bliss::concurrent::DynamicOMPRunner(3));
    //// test sequential.
    INFOF("dynamic OMP runner\n");
    sr->addTask(std::move(std::shared_ptr<bliss::concurrent::Runnable>(new Test("task 1"))));
    sr->addTask(std::move(std::shared_ptr<bliss::concurrent::Runnable>(new Test("task 2"))));
    sr->addTask(std::move(std::shared_ptr<bliss::concurrent::Runnable>(new Test("task 3"))));
    sr->addTask(std::move(std::shared_ptr<bliss::concurrent::Runnable>(new Test("task 4"))));
    sr->disableAdd();

    t2 = std::chrono::high_resolution_clock::now();
    time_span =
        std::chrono::duration_cast<std::chrono::duration<double>>(
            t2 - t1);
    INFO( "setup: time " << time_span.count() );

    t1 = std::chrono::high_resolution_clock::now();

    sr->operator()();

    t2 = std::chrono::high_resolution_clock::now();
    time_span =
        std::chrono::duration_cast<std::chrono::duration<double>>(
            t2 - t1);
    INFO( "run: time " << time_span.count() );

  }


  {
    t1 = std::chrono::high_resolution_clock::now();

    cc.store(0, std::memory_order_relaxed);
    std::shared_ptr<bliss::concurrent::Runner> sr(new bliss::concurrent::DynamicOMPRunner(4));

    //// test sequential.
    INFOF("dynamic OMP runner \n");
    sr->addTask(std::move(std::shared_ptr<bliss::concurrent::Runnable>(new Test4("task 1", sr))));
    sr->addTask(std::move(std::shared_ptr<bliss::concurrent::Runnable>(new Test4("task 2", sr))));
    sr->addTask(std::move(std::shared_ptr<bliss::concurrent::Runnable>(new Test4("task 3", sr))));
    sr->addTask(std::move(std::shared_ptr<bliss::concurrent::Runnable>(new Test4("task 4", sr))));

    t2 = std::chrono::high_resolution_clock::now();
    time_span =
        std::chrono::duration_cast<std::chrono::duration<double>>(
            t2 - t1);
    INFO( "setup: time " << time_span.count() );

    t1 = std::chrono::high_resolution_clock::now();

    sr->operator()();

    t2 = std::chrono::high_resolution_clock::now();
    time_span =
        std::chrono::duration_cast<std::chrono::duration<double>>(
            t2 - t1);
    INFO( "run: time " << time_span.count() );

  }

  {
    t1 = std::chrono::high_resolution_clock::now();

    std::shared_ptr<bliss::concurrent::Runner> sr(new bliss::concurrent::DynamicOMPRunner(3));

    //// test sequential.
    INFOF("dynamic OMP runner with self-adding tasks and blocked queue\n");
    sr->addTask(std::move(std::shared_ptr<bliss::concurrent::Runnable>(new Test2("task 1", sr))));
    sr->addTask(std::move(std::shared_ptr<bliss::concurrent::Runnable>(new Test2("task 2", sr))));
    sr->addTask(std::move(std::shared_ptr<bliss::concurrent::Runnable>(new Test2("task 3", sr))));
    sr->addTask(std::move(std::shared_ptr<bliss::concurrent::Runnable>(new Test2("task 4", sr))));
    sr->disableAdd();

    t2 = std::chrono::high_resolution_clock::now();
    time_span =
        std::chrono::duration_cast<std::chrono::duration<double>>(
            t2 - t1);
    INFO( "setup: time " << time_span.count() );

    t1 = std::chrono::high_resolution_clock::now();

    sr->operator()();

    t2 = std::chrono::high_resolution_clock::now();
    time_span =
        std::chrono::duration_cast<std::chrono::duration<double>>(
            t2 - t1);
    INFO( "run: time " << time_span.count() );

  }

    {
      t1 = std::chrono::high_resolution_clock::now();

      INFOF("sequential OMP with separate thread to stop\n");

      std::shared_ptr<bliss::concurrent::Runner> sr(new bliss::concurrent::SequentialRunner());

      sr->addTask(std::move(std::shared_ptr<bliss::concurrent::Runnable>(new Test2("task 1", sr))));
      sr->addTask(std::move(std::shared_ptr<bliss::concurrent::Runnable>(new Test2("task 2", sr))));
      sr->addTask(std::move(std::shared_ptr<bliss::concurrent::Runnable>(new Test2("task 3", sr))));
      sr->addTask(std::move(std::shared_ptr<bliss::concurrent::Runnable>(new Test2("task 4", sr))));
      sr->addTask(std::move(std::shared_ptr<bliss::concurrent::Runnable>(new Test2("task 5", sr))));


      std::shared_ptr<bliss::concurrent::Runner> sr2(new bliss::concurrent::DynamicOMPRunner(2));
      sr2->addTask(std::move(std::shared_ptr<bliss::concurrent::Runnable>(new Test3("control", sr))));
      sr2->addTask(std::move(sr->shared_from_this()));
      sr2->disableAdd();

      t2 = std::chrono::high_resolution_clock::now();
      time_span =
          std::chrono::duration_cast<std::chrono::duration<double>>(
              t2 - t1);
      INFO( "setup: time " << time_span.count() );

      t1 = std::chrono::high_resolution_clock::now();

      sr2->operator()();

      t2 = std::chrono::high_resolution_clock::now();
      time_span =
          std::chrono::duration_cast<std::chrono::duration<double>>(
              t2 - t1);
      INFO( "run: time " << time_span.count() );
    }

  {
    t1 = std::chrono::high_resolution_clock::now();

    INFOF("dynamic OMP with separate thread to stop, no reinsert\n");

    std::shared_ptr<bliss::concurrent::Runner> sr(new bliss::concurrent::DynamicOMPRunner(3));

    sr->addTask(std::move(std::shared_ptr<bliss::concurrent::Runnable>(new Test("task 1"))));
    sr->addTask(std::move(std::shared_ptr<bliss::concurrent::Runnable>(new Test("task 2"))));
    sr->addTask(std::move(std::shared_ptr<bliss::concurrent::Runnable>(new Test("task 3"))));
    sr->addTask(std::move(std::shared_ptr<bliss::concurrent::Runnable>(new Test("task 4"))));
    sr->addTask(std::move(std::shared_ptr<bliss::concurrent::Runnable>(new Test("task 5"))));

    std::shared_ptr<bliss::concurrent::Runner> sr2(new bliss::concurrent::DynamicOMPRunner(2));
    sr2->addTask(std::move(std::shared_ptr<bliss::concurrent::Runnable>(new Test3("control", sr))));
    sr2->addTask(std::move(sr->shared_from_this()));
    sr2->disableAdd();


    t2 = std::chrono::high_resolution_clock::now();
    time_span =
        std::chrono::duration_cast<std::chrono::duration<double>>(
            t2 - t1);
    INFO( "setup: time " << time_span.count() );

    t1 = std::chrono::high_resolution_clock::now();

    sr2->operator()();

    t2 = std::chrono::high_resolution_clock::now();
    time_span =
        std::chrono::duration_cast<std::chrono::duration<double>>(
            t2 - t1);
    INFO( "run: time " << time_span.count() );
  }


    {
      t1 = std::chrono::high_resolution_clock::now();

      INFOF("dynamic OMP with separate thread to stop\n");

      std::shared_ptr<bliss::concurrent::Runner> sr(new bliss::concurrent::DynamicOMPRunner(3));

      sr->addTask(std::move(std::shared_ptr<bliss::concurrent::Runnable>(new Test2("task 1", sr))));
      sr->addTask(std::move(std::shared_ptr<bliss::concurrent::Runnable>(new Test2("task 2", sr))));
      sr->addTask(std::move(std::shared_ptr<bliss::concurrent::Runnable>(new Test2("task 3", sr))));
      sr->addTask(std::move(std::shared_ptr<bliss::concurrent::Runnable>(new Test2("task 4", sr))));
      sr->addTask(std::move(std::shared_ptr<bliss::concurrent::Runnable>(new Test2("task 5", sr))));


      std::shared_ptr<bliss::concurrent::Runner> sr2(new bliss::concurrent::DynamicOMPRunner(2));

      sr2->addTask(std::move(std::shared_ptr<bliss::concurrent::Runnable>(new Test3("control", sr))));
      sr2->addTask(std::move(sr->shared_from_this()));
      sr2->disableAdd();

      t2 = std::chrono::high_resolution_clock::now();
      time_span =
          std::chrono::duration_cast<std::chrono::duration<double>>(
              t2 - t1);
      INFO( "setup: time " << time_span.count() );

      t1 = std::chrono::high_resolution_clock::now();

      sr2->operator()();

      t2 = std::chrono::high_resolution_clock::now();
      time_span =
          std::chrono::duration_cast<std::chrono::duration<double>>(
              t2 - t1);
      INFO( "run: time " << time_span.count() );
    }

    {
      t1 = std::chrono::high_resolution_clock::now();

      cc.store(0, std::memory_order_relaxed);
      cc2.store(0, std::memory_order_relaxed);

      INFOF("index build pattern\n");

      std::shared_ptr<bliss::concurrent::Runner> commRunner(new bliss::concurrent::SequentialRunner());
      std::shared_ptr<bliss::concurrent::Runner> workRunner(new bliss::concurrent::DynamicOMPRunner(3));

      // initially put as many senders as there are threads.  hook them to receivers
      workRunner->addTask(std::move(std::shared_ptr<bliss::concurrent::Task>(new Source<Sender<Receiver<Sink> > >("reader 1", commRunner, workRunner))));
      workRunner->addTask(std::move(std::shared_ptr<bliss::concurrent::Task>(new Source<Sender<Receiver<Sink> > >("reader 2", commRunner, workRunner))));
      workRunner->addTask(std::move(std::shared_ptr<bliss::concurrent::Task>(new Source<Sender<Receiver<Sink> > >("reader 3", commRunner, workRunner))));



      bliss::concurrent::PersonalizedOMPRunner appRunner(2);
      appRunner.addTask(std::move(commRunner->shared_from_this()));
      appRunner.addTask(std::move(workRunner->shared_from_this()));
      appRunner.disableAdd();


      t2 = std::chrono::high_resolution_clock::now();
      time_span =
          std::chrono::duration_cast<std::chrono::duration<double>>(
              t2 - t1);
      INFO( "setup: time " << time_span.count() );

      t1 = std::chrono::high_resolution_clock::now();

      appRunner();

      t2 = std::chrono::high_resolution_clock::now();
      time_span =
          std::chrono::duration_cast<std::chrono::duration<double>>(
              t2 - t1);
      INFO( "run: time " << time_span.count() );
    }

    {
      t1 = std::chrono::high_resolution_clock::now();

      cc.store(0, std::memory_order_relaxed);
      cc2.store(0, std::memory_order_relaxed);

      INFOF("index query pattern\n");

      std::shared_ptr<bliss::concurrent::Runner> commRunner(new bliss::concurrent::SequentialRunner());
      std::shared_ptr<bliss::concurrent::Runner> workRunner(new bliss::concurrent::DynamicOMPRunner(3));


      // initially put as many senders as there are threads.  hook them to receivers
      workRunner->addTask(std::move(std::shared_ptr<bliss::concurrent::Task>(new Source<Sender<Receiver<Node<Sender<Receiver<Sink> > > > > >("query 1", commRunner, workRunner))));
      workRunner->addTask(std::move(std::shared_ptr<bliss::concurrent::Task>(new Source<Sender<Receiver<Node<Sender<Receiver<Sink> > > > > >("query 2", commRunner, workRunner))));
      workRunner->addTask(std::move(std::shared_ptr<bliss::concurrent::Task>(new Source<Sender<Receiver<Node<Sender<Receiver<Sink> > > > > >("query 3", commRunner, workRunner))));
      
      

      bliss::concurrent::PersonalizedOMPRunner appRunner(2);

      appRunner.addTask(std::move(commRunner->shared_from_this()));
      appRunner.addTask(std::move(workRunner->shared_from_this()));
      appRunner.disableAdd();


      t2 = std::chrono::high_resolution_clock::now();
      time_span =
          std::chrono::duration_cast<std::chrono::duration<double>>(
              t2 - t1);
      INFO( "setup: time " << time_span.count() );

      t1 = std::chrono::high_resolution_clock::now();

      appRunner();

      t2 = std::chrono::high_resolution_clock::now();
      time_span =
          std::chrono::duration_cast<std::chrono::duration<double>>(
              t2 - t1);
      INFO( "run: time " << time_span.count() );
    }


}
