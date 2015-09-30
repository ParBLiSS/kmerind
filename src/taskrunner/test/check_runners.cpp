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
//#include <unistd.h>  // sleep
#include <chrono>

//#include "concurrent/mpi_runner.hpp"
#include "taskrunner/uniform_omp_runner.hpp"
#include "taskrunner/personalized_omp_runner.hpp"
#include "taskrunner/dynamic_omp_runner.hpp"
#include "taskrunner/sequential_runner.hpp"
#include "taskrunner/task.hpp"

#include <xmmintrin.h>

static const int iter = 100000;
std::atomic<int> cc;
std::atomic<int> cc2;


class BasicTask : public bliss::concurrent::Task
{
  public:
    BasicTask(std::string msg) : bliss::concurrent::Task(), content(msg) {};
    virtual ~BasicTask() {};

    virtual void operator()() {
      //INFOF("tid %d: %s", omp_get_thread_num(), content.c_str());
    }

  protected:
    std::string content;

};


class SelfAddingTask : public bliss::concurrent::Task
{
  protected:
    std::string content;
    std::shared_ptr<bliss::concurrent::Runner> parent;
  public:
    SelfAddingTask(std::string msg,
          const std::shared_ptr<bliss::concurrent::Runner>& _parent) :
            bliss::concurrent::Task(), content(msg), parent(_parent) {};
    virtual ~SelfAddingTask() {};

    virtual void operator()() {
      //INFOF("tid %d: %s", omp_get_thread_num(), content.c_str());
      parent->addTask(this->shared_from_this());
//      if (parent.addTask(this))
//        INFOF("tid %d: %s reinserted", omp_get_thread_num(), content.c_str());
    }


};

class DelayedDisableTask : public bliss::concurrent::Task
{
  protected:
    std::string content;
    std::shared_ptr<bliss::concurrent::Runner> other;
  public:
    DelayedDisableTask(std::string msg, const std::shared_ptr<bliss::concurrent::Runner>& _other) : bliss::concurrent::Task(), content(msg), other(_other) {};
    virtual ~DelayedDisableTask() {};

    virtual void operator()() {
      INFOF("tid %d: %s start", omp_get_thread_num(), content.c_str());

      for (int i = 0; i < 1000; ++i)
        _mm_pause();

      INFOF("tid %d: %s stop", omp_get_thread_num(), content.c_str());
      other->disableAdd();
      INFOF("tid %d: %s stopped", omp_get_thread_num(), content.c_str());
    }

};



class FixedIterationDisablingTask : public bliss::concurrent::Task
{
  protected:
    std::shared_ptr<bliss::concurrent::Runner> work_r;

    std::string content;
  public:
    FixedIterationDisablingTask(std::string msg,
           const std::shared_ptr<bliss::concurrent::Runner>& work) :
             work_r(work), content(msg) {};

    virtual ~FixedIterationDisablingTask() {};

    virtual void operator()() {
      auto id = cc.fetch_add(1, std::memory_order_relaxed);

      if (id < iter) {
        //INFOF("tid %d: %s, count %d", omp_get_thread_num(), content.c_str(), i);
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
      //INFOF("tid %d: %s", omp_get_thread_num(), content.c_str());
      //INFOF("%dS ", omp_get_thread_num());
      _mm_pause();
//      INFOF("tid %d: %s add recv %d to comm runner size %lu, disabled %s ", omp_get_thread_num(), content.c_str(), id, comm_r->getTaskCount(), (comm_r->isAddDisabled() ? "y" : "n"));
      comm_r->addTask(std::shared_ptr<Runnable>(new Receiver("recv", comm_r, work_r, id)));

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
//      _mm_pause();
//      _mm_pause();
//      _mm_pause();
//      _mm_pause();
//      _mm_pause();
//      _mm_pause();
//      _mm_pause();
//      _mm_pause();
//      _mm_pause();

//      INFOF("tid %d: %s add send %d to comm runner size %lu, disabled %s ", omp_get_thread_num(), content.c_str(), id, comm_r->getTaskCount(), (comm_r->isAddDisabled() ? "y" : "n"));
      comm_r->addTask(std::shared_ptr<Runnable>(new SENDER("node", comm_r, work_r, id)));

    }

};

template<class SENDER>
class Source : public bliss::concurrent::Task
{
  protected:
    std::shared_ptr<bliss::concurrent::Runner> work_r;
    std::shared_ptr<bliss::concurrent::Runner> comm_r;
    std::string content;

    bool justFinished;

  public:
    Source(std::string msg,
           const std::shared_ptr<bliss::concurrent::Runner>& comm,
           const std::shared_ptr<bliss::concurrent::Runner>& work) :
             work_r(work), comm_r(comm), content(msg), justFinished(false) {};

    virtual ~Source() {};

    virtual void operator()() {
      auto id = cc.fetch_add(1, std::memory_order_relaxed);

      if (id < iter) {
          //INFOF("tid %d: %s, count %d", omp_get_thread_num(), content.c_str(), i);
          //INFOF("%dR%d ", omp_get_thread_num(), i);
          _mm_pause();
//          _mm_pause();
//          _mm_pause();

//        INFOF("tid %d: %s add send %d to comm runner size %lu, disabled %s ", omp_get_thread_num(), content.c_str(), id, comm_r->getTaskCount(), (comm_r->isAddDisabled() ? "y" : "n"));
//        INFOF("tid %d: %s add self %d to work runner size %lu, disabled %s ", omp_get_thread_num(), content.c_str(), id, work_r->getTaskCount(), (work_r->isAddDisabled() ? "y" : "n"));
        comm_r->addTask(std::shared_ptr<Runnable>(new SENDER("sender", comm_r, work_r, id)));

        work_r->addTask(this->shared_from_this());
      } else {
        if (!justFinished) {
          justFinished = true;
          INFOF("tid %d: %s finished, count %d", omp_get_thread_num(), content.c_str(), id);
        }
      }
    }

};

class Sink : public bliss::concurrent::Task
{
  protected:
    std::shared_ptr<bliss::concurrent::Runner> work_r;
    std::shared_ptr<bliss::concurrent::Runner> comm_r;
    std::string content;

    bool justFinished;
    int id;

  public:
    Sink(std::string msg,
         const std::shared_ptr<bliss::concurrent::Runner>& comm,
         const std::shared_ptr<bliss::concurrent::Runner>& work, int _id) :
           work_r(work), comm_r(comm), content(msg), justFinished(false), id(_id) {};

    virtual ~Sink() {};

    virtual void operator()() {
      auto i = cc2.fetch_add(1, std::memory_order_relaxed);


      if (i+1 >= iter) {
        INFOF("tid %d: disabled runners at %d, item %d, iter %d", omp_get_thread_num(), id, i, iter);

        if (!justFinished) {
          justFinished = true;
          work_r->disableAdd();
          comm_r->disableAdd();
          INFOF("tid %d: sink id %d, item %d, iter %d, work runner size %lu, disabled %s", omp_get_thread_num(), id ,i, iter, work_r->getTaskCount(), (work_r->isAddDisabled() ? "y" : "n"));
      }

      } else {
        //INFOF("%dW ", omp_get_thread_num());
        _mm_pause();
//        _mm_pause();
//        _mm_pause();
//        _mm_pause();
//        _mm_pause();
//        _mm_pause();
//        _mm_pause();
//        _mm_pause();
//        _mm_pause();
//        _mm_pause();
//        INFOF("tid %d: sink id %d, item %d, iter %d, work runner size %lu, disabled %s", omp_get_thread_num(), id ,i, iter, work_r->getTaskCount(), (work_r->isAddDisabled() ? "y" : "n"));

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
//        _mm_pause();

//        INFOF("tid %d: %s add write %d to work runner size %lu, disabled %s ", omp_get_thread_num(), content.c_str(), id, work_r->getTaskCount(), (work_r->isAddDisabled() ? "y" : "n"));

        work_r->addTask(std::shared_ptr<Runnable>(new WRITER("write", comm_r, work_r, id)));

    }

};



int main(int argc, char** argv) {

  std::chrono::high_resolution_clock::time_point t1, t2;
  std::chrono::duration<double> time_span;


  {
    t1 = std::chrono::high_resolution_clock::now();

    std::shared_ptr<bliss::concurrent::Runner> sr(new bliss::concurrent::SequentialRunner());
    //// test sequential.
    INFOF("\nserial runner");
    sr->addTask(std::shared_ptr<bliss::concurrent::Runnable>(new BasicTask("task 1")));
    sr->addTask(std::shared_ptr<bliss::concurrent::Runnable>(new BasicTask("task 2")));
    sr->addTask(std::shared_ptr<bliss::concurrent::Runnable>(new BasicTask("task 3")));
    sr->addTask(std::shared_ptr<bliss::concurrent::Runnable>(new BasicTask("task 4")));
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
    INFOF("\nserial runner");
    sr->addTask(std::shared_ptr<bliss::concurrent::Runnable>(new FixedIterationDisablingTask("task 1", sr)));


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
    INFOF("\npersonalized OMP runner");
    sr->addTask(std::shared_ptr<bliss::concurrent::Runnable>(new BasicTask("task 1")));
    sr->addTask(std::shared_ptr<bliss::concurrent::Runnable>(new BasicTask("task 2")));
    sr->addTask(std::shared_ptr<bliss::concurrent::Runnable>(new BasicTask("task 3")));
    sr->addTask(std::shared_ptr<bliss::concurrent::Runnable>(new BasicTask("task 4")));
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
    INFOF("\npersonalized OMP runner");
    sr->addTask(std::shared_ptr<bliss::concurrent::Runnable>(new FixedIterationDisablingTask("task 1", sr)));
    sr->addTask(std::shared_ptr<bliss::concurrent::Runnable>(new FixedIterationDisablingTask("task 2", sr)));
    sr->addTask(std::shared_ptr<bliss::concurrent::Runnable>(new FixedIterationDisablingTask("task 3", sr)));
    sr->addTask(std::shared_ptr<bliss::concurrent::Runnable>(new FixedIterationDisablingTask("task 4", sr)));

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
    INFOF("\nuniform OMP runner");
    sr->addTask(std::shared_ptr<bliss::concurrent::Runnable>(new BasicTask("task 1")));
    sr->addTask(std::shared_ptr<bliss::concurrent::Runnable>(new BasicTask("task 2")));
    sr->addTask(std::shared_ptr<bliss::concurrent::Runnable>(new BasicTask("task 3")));
    sr->addTask(std::shared_ptr<bliss::concurrent::Runnable>(new BasicTask("task 4")));
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
    INFOF("\nuniform OMP runner");
    sr->addTask(std::shared_ptr<bliss::concurrent::Runnable>(new FixedIterationDisablingTask("task 1", sr)));
    sr->addTask(std::shared_ptr<bliss::concurrent::Runnable>(new FixedIterationDisablingTask("task 2", sr)));
    sr->addTask(std::shared_ptr<bliss::concurrent::Runnable>(new FixedIterationDisablingTask("task 3", sr)));
    sr->addTask(std::shared_ptr<bliss::concurrent::Runnable>(new FixedIterationDisablingTask("task 4", sr)));

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
    INFOF("\ndynamic OMP runner");
    sr->addTask(std::shared_ptr<bliss::concurrent::Runnable>(new BasicTask("task 1")));
    sr->addTask(std::shared_ptr<bliss::concurrent::Runnable>(new BasicTask("task 2")));
    sr->addTask(std::shared_ptr<bliss::concurrent::Runnable>(new BasicTask("task 3")));
    sr->addTask(std::shared_ptr<bliss::concurrent::Runnable>(new BasicTask("task 4")));
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
    INFOF("\ndynamic OMP runner ");
    sr->addTask(std::shared_ptr<bliss::concurrent::Runnable>(new FixedIterationDisablingTask("task 1", sr)));
    sr->addTask(std::shared_ptr<bliss::concurrent::Runnable>(new FixedIterationDisablingTask("task 2", sr)));
    sr->addTask(std::shared_ptr<bliss::concurrent::Runnable>(new FixedIterationDisablingTask("task 3", sr)));
    sr->addTask(std::shared_ptr<bliss::concurrent::Runnable>(new FixedIterationDisablingTask("task 4", sr)));

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
    INFOF("\ndynamic OMP runner with self-adding tasks and blocked queue");
    sr->addTask(std::shared_ptr<bliss::concurrent::Runnable>(new SelfAddingTask("task 1", sr)));
    sr->addTask(std::shared_ptr<bliss::concurrent::Runnable>(new SelfAddingTask("task 2", sr)));
    sr->addTask(std::shared_ptr<bliss::concurrent::Runnable>(new SelfAddingTask("task 3", sr)));
    sr->addTask(std::shared_ptr<bliss::concurrent::Runnable>(new SelfAddingTask("task 4", sr)));
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

      INFOF("\nsequential OMP with separate thread to stop");

      std::shared_ptr<bliss::concurrent::Runner> sr(new bliss::concurrent::SequentialRunner());

      sr->addTask(std::shared_ptr<bliss::concurrent::Runnable>(new SelfAddingTask("task 1", sr)));
      sr->addTask(std::shared_ptr<bliss::concurrent::Runnable>(new SelfAddingTask("task 2", sr)));
      sr->addTask(std::shared_ptr<bliss::concurrent::Runnable>(new SelfAddingTask("task 3", sr)));
      sr->addTask(std::shared_ptr<bliss::concurrent::Runnable>(new SelfAddingTask("task 4", sr)));
      sr->addTask(std::shared_ptr<bliss::concurrent::Runnable>(new SelfAddingTask("task 5", sr)));


      std::shared_ptr<bliss::concurrent::Runner> sr2(new bliss::concurrent::DynamicOMPRunner(2));
      sr2->addTask(std::shared_ptr<bliss::concurrent::Runnable>(new DelayedDisableTask("control", sr)));
      sr2->addTask(sr->shared_from_this());
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

    INFOF("\ndynamic OMP with separate thread to stop, no reinsert");

    std::shared_ptr<bliss::concurrent::Runner> sr(new bliss::concurrent::DynamicOMPRunner(3));

    sr->addTask(std::shared_ptr<bliss::concurrent::Runnable>(new BasicTask("task 1")));
    sr->addTask(std::shared_ptr<bliss::concurrent::Runnable>(new BasicTask("task 2")));
    sr->addTask(std::shared_ptr<bliss::concurrent::Runnable>(new BasicTask("task 3")));
    sr->addTask(std::shared_ptr<bliss::concurrent::Runnable>(new BasicTask("task 4")));
    sr->addTask(std::shared_ptr<bliss::concurrent::Runnable>(new BasicTask("task 5")));

    std::shared_ptr<bliss::concurrent::Runner> sr2(new bliss::concurrent::DynamicOMPRunner(2));
    sr2->addTask(std::shared_ptr<bliss::concurrent::Runnable>(new DelayedDisableTask("control", sr)));
    sr2->addTask(sr->shared_from_this());
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

      INFOF("\ndynamic OMP with separate thread to stop");

      std::shared_ptr<bliss::concurrent::Runner> sr(new bliss::concurrent::DynamicOMPRunner(3));

      sr->addTask(std::shared_ptr<bliss::concurrent::Runnable>(new SelfAddingTask("task 1", sr)));
      sr->addTask(std::shared_ptr<bliss::concurrent::Runnable>(new SelfAddingTask("task 2", sr)));
      sr->addTask(std::shared_ptr<bliss::concurrent::Runnable>(new SelfAddingTask("task 3", sr)));
      sr->addTask(std::shared_ptr<bliss::concurrent::Runnable>(new SelfAddingTask("task 4", sr)));
      sr->addTask(std::shared_ptr<bliss::concurrent::Runnable>(new SelfAddingTask("task 5", sr)));


      std::shared_ptr<bliss::concurrent::Runner> sr2(new bliss::concurrent::DynamicOMPRunner(2));

      sr2->addTask(std::shared_ptr<bliss::concurrent::Runnable>(new DelayedDisableTask("control", sr)));
      sr2->addTask(sr->shared_from_this());
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

      INFOF("\nindex build pattern.");

      std::shared_ptr<bliss::concurrent::Runner> commRunner(new bliss::concurrent::SequentialRunner());
      std::shared_ptr<bliss::concurrent::Runner> workRunner(new bliss::concurrent::DynamicOMPRunner(3));

      // initially put as many senders as there are threads.  hook them to receivers
      workRunner->addTask(std::shared_ptr<bliss::concurrent::Task>(new Source<Sender<Receiver<Sink> > >("reader 1", commRunner, workRunner)));
      workRunner->addTask(std::shared_ptr<bliss::concurrent::Task>(new Source<Sender<Receiver<Sink> > >("reader 2", commRunner, workRunner)));
      workRunner->addTask(std::shared_ptr<bliss::concurrent::Task>(new Source<Sender<Receiver<Sink> > >("reader 3", commRunner, workRunner)));



      bliss::concurrent::PersonalizedOMPRunner appRunner(2);
      appRunner.addTask(commRunner->shared_from_this());
      appRunner.addTask(workRunner->shared_from_this());
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

      INFOF("\nindex query pattern");

      std::shared_ptr<bliss::concurrent::Runner> commRunner(new bliss::concurrent::SequentialRunner());
      std::shared_ptr<bliss::concurrent::Runner> workRunner(new bliss::concurrent::DynamicOMPRunner(3));


      // initially put as many senders as there are threads.  hook them to receivers
      workRunner->addTask(std::shared_ptr<bliss::concurrent::Task>(new Source<Sender<Receiver<Node<Sender<Receiver<Sink> > > > > >("query 1", commRunner, workRunner)));
      workRunner->addTask(std::shared_ptr<bliss::concurrent::Task>(new Source<Sender<Receiver<Node<Sender<Receiver<Sink> > > > > >("query 2", commRunner, workRunner)));
      workRunner->addTask(std::shared_ptr<bliss::concurrent::Task>(new Source<Sender<Receiver<Node<Sender<Receiver<Sink> > > > > >("query 3", commRunner, workRunner)));
      
      

      bliss::concurrent::PersonalizedOMPRunner appRunner(2);

      appRunner.addTask(commRunner->shared_from_this());
      appRunner.addTask(workRunner->shared_from_this());
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
