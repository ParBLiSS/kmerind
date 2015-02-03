/**
 * @file		test_runners.cpp
 * @ingroup
 * @author	tpan
 * @brief
 * @details
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */

#include <iostream>
#include <unistd.h>  // sleep

//#include "concurrent/mpi_runner.hpp"
#include "wip/uniform_omp_runner.hpp"
#include "wip/dynamic_omp_runner.hpp"
#include "wip/sequential_runner.hpp"
#include "wip/task.hpp"

static const int iter = 1000000;

class Test : public bliss::concurrent::Task
{
  public:
    Test(std::string msg) : bliss::concurrent::Task(), content(msg) {};
    virtual ~Test() {};

    virtual void operator()() {
      //printf("tid %d: %s\n", omp_get_thread_num(), content.c_str());
    }

  protected:
    std::string content;

};


class Test2 : public bliss::concurrent::Task
{
  protected:
    std::string content;
    bliss::concurrent::Runner& parent;
  public:
    Test2(std::string msg, bliss::concurrent::Runner& _parent) : bliss::concurrent::Task(), content(msg), parent(_parent) {};
    virtual ~Test2() {};

    virtual void operator()() {
      //printf("tid %d: %s\n", omp_get_thread_num(), content.c_str());
      parent.addTask(this);
//      if (parent.addTask(this))
//        printf("tid %d: %s reinserted\n", omp_get_thread_num(), content.c_str());
    }


};

class Test3 : public bliss::concurrent::Task
{
  protected:
    std::string content;
    bliss::concurrent::Runner& other;
  public:
    Test3(std::string msg, bliss::concurrent::Runner& _other) : bliss::concurrent::Task(), content(msg), other(_other) {};
    virtual ~Test3() {};

    virtual void operator()() {
      printf("tid %d: %s start\n", omp_get_thread_num(), content.c_str());

      usleep(1000);

      printf("tid %d: %s stop\n", omp_get_thread_num(), content.c_str());
      other.disableAdd();
      printf("tid %d: %s stopped\n", omp_get_thread_num(), content.c_str());
    }

};

class Sender : public bliss::concurrent::Task
{
  protected:
    bliss::concurrent::Runner& r;
    std::string content;



  public:
    Sender(std::string msg, bliss::concurrent::Runner& _other) : r(_other), content(msg) {};

    virtual ~Sender() {};

    virtual void operator()() {
      //printf("tid %d: %s\n", omp_get_thread_num(), content.c_str());
      //printf("%dS ", omp_get_thread_num());
      _mm_pause();
      //r.addTask(this);

    }

};

template<class SENDER>
class Reader : public bliss::concurrent::Task
{
  protected:
    volatile std::atomic<int>& cc;
    bliss::concurrent::Runner& r;
    bliss::concurrent::Runner& comm_r;
    std::string content;
  public:
    Reader(std::string msg, bliss::concurrent::Runner& mine, bliss::concurrent::Runner& _other, std::atomic<int>& count) : cc(count), r(mine), comm_r(_other), content(msg) {};

    virtual ~Reader() {};

    virtual void operator()() {
      auto i = cc.fetch_add(1, std::memory_order_relaxed);

      if (i < iter) {
        //printf("tid %d: %s, count %d\n", omp_get_thread_num(), content.c_str(), i);
        //printf("%dR%d ", omp_get_thread_num(), i);
        _mm_pause();
        _mm_pause();
        _mm_pause();

        comm_r.addTask(new SENDER("sender", comm_r));
        r.addTask(this);
      }
    }

};

class Writer : public bliss::concurrent::Task
{
  protected:
    std::string content;
  public:
    Writer(std::string msg) : content(msg) {};

    virtual ~Writer() {};

    virtual void operator()() {
//      printf("tid %d: %s\n", omp_get_thread_num(), content.c_str());
      //printf("%dW ", omp_get_thread_num());
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
    }

};

template<class WRITER>
class Receiver : public bliss::concurrent::Task
{
  protected:
    std::atomic<int>& cc;
    bliss::concurrent::Runner& r;
    bliss::concurrent::Runner& work_r;
    std::string content;
  public:
    Receiver(std::string msg, bliss::concurrent::Runner& mine, bliss::concurrent::Runner& _other, std::atomic<int>& count) : cc(count), r(mine), work_r(_other), content(msg) {};

    virtual ~Receiver() {};

    virtual void operator()() {
      auto i = cc.load(std::memory_order_relaxed);
//      printf("tid %d: %s.  count %d\n", omp_get_thread_num(), content.c_str(), i);

      if (i < iter) {
        //printf("%dC%d ", omp_get_thread_num(), i);
        _mm_pause();

        work_r.addTask(new WRITER("write"));

         r.addTask(this);
      } else {
        work_r.disableAdd();
        r.disableAdd();
      }
    }

};



int main(int argc, char** argv) {

  {
    Test t1("task 1");
    Test t2("task 2");
    Test t3("task 3");
    Test t4("task 4");

    bliss::concurrent::SequentialRunner sr;
    //// test sequential.
    printf("serial runner\n");
    sr.addTask(&t1);
    sr.addTask(&t2);
    sr.addTask(&t3);
    sr.addTask(&t4);
    sr.disableAdd();
    sr();

  }

  {
    Test t1("task 1");
    Test t2("task 2");
    Test t3("task 3");
    Test t4("task 4");

    bliss::concurrent::UniformOMPRunner sr(3);
    //// test sequential.
    printf("uniform OMP runner\n");
    sr.addTask(&t1);
    sr.addTask(&t2);
    sr.addTask(&t3);
    sr.addTask(&t4);
    sr.disableAdd();
    sr();


  }



  {
    Test t1("task 1");
    Test t2("task 2");
    Test t3("task 3");
    Test t4("task 4");

    bliss::concurrent::DynamicOMPRunner sr(3);
    //// test sequential.
    printf("dynamic OMP runner\n");
    sr.addTask(&t1);
    sr.addTask(&t2);
    sr.addTask(&t3);
    sr.addTask(&t4);
    sr.disableAdd();
    sr();

  }


  {
    bliss::concurrent::DynamicOMPRunner sr(3);
    Test2 t1("task 1", sr);
    Test2 t2("task 2", sr);
    Test2 t3("task 3", sr);
    Test2 t4("task 4", sr);

    //// test sequential.
    printf("dynamic OMP runner with self-adding tasks and blocked queue\n");
    sr.addTask(&t1);
    sr.addTask(&t2);
    sr.addTask(&t3);
    sr.addTask(&t4);
    sr.disableAdd();
    sr();

  }

    {
      printf("sequential OMP with separate thread to stop\n");

      bliss::concurrent::SequentialRunner sr;

      Test2 t1("task 1", sr);
      Test2 t2("task 2", sr);
      Test2 t3("task 3", sr);
      Test2 t4("task 4", sr);
      Test2 t5("task 5", sr);

      sr.addTask(&t1);
      sr.addTask(&t2);
      sr.addTask(&t3);
      sr.addTask(&t4);
      sr.addTask(&t5);


      Test3 t("control", sr);

      bliss::concurrent::DynamicOMPRunner sr2(2);
      sr2.addTask(&t);
      sr2.addTask(&sr);
      sr2.disableAdd();
      sr2();
    }

  {
    printf("dynamic OMP with separate thread to stop, no reinsert\n");

    bliss::concurrent::DynamicOMPRunner sr(3);

    Test t1("task 1");
    Test t2("task 2");
    Test t3("task 3");
    Test t4("task 4");
    Test t5("task 5");

    sr.addTask(&t1);
    sr.addTask(&t2);
    sr.addTask(&t3);
    sr.addTask(&t4);
    sr.addTask(&t5);


    Test3 t("control", sr);

    bliss::concurrent::DynamicOMPRunner sr2(2);
    sr2.addTask(&t);
    sr2.addTask(&sr);
    sr2.disableAdd();
    sr2();
  }


    {
      printf("dynamic OMP with separate thread to stop\n");

      bliss::concurrent::DynamicOMPRunner sr(3);

      Test2 t1("task 1", sr);
      Test2 t2("task 2", sr);
      Test2 t3("task 3", sr);
      Test2 t4("task 4", sr);
      Test2 t5("task 5", sr);

      sr.addTask(&t1);
      sr.addTask(&t2);
      sr.addTask(&t3);
      sr.addTask(&t4);
      sr.addTask(&t5);


      Test3 t("control", sr);

      bliss::concurrent::DynamicOMPRunner sr2(2);
      sr2.addTask(&t);
      sr2.addTask(&sr);
      sr2.disableAdd();
      sr2();
    }


    {
      printf("index build pattern\n");

      bliss::concurrent::DynamicOMPRunner workRunner(4);
      bliss::concurrent::SequentialRunner commRunner;

      std::atomic<int> count(0);
      // initially put as many senders as there are threads.  hook them to receivers
      Reader<Sender> t1("reader 1", commRunner, workRunner, count);
      Reader<Sender> t2("reader 2", commRunner, workRunner, count);
      Reader<Sender> t3("reader 3", commRunner, workRunner, count);

      workRunner.addTask(&t1);
      workRunner.addTask(&t2);
      workRunner.addTask(&t3);


      Receiver<Writer> recvr("receiver", commRunner, workRunner, count);

      commRunner.addTask(&recvr);


      bliss::concurrent::DynamicOMPRunner appRunner(2);
      appRunner.addTask(&commRunner);
      appRunner.addTask(&workRunner);
      appRunner.disableAdd();
      appRunner();
    }



}
