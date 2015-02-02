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

//#include "concurrent/mpi_runner.hpp"
#include "wip/replicated_omp_runner.hpp"
#include "wip/dynamic_omp_runner.hpp"
#include "wip/sequential_runner.hpp"
#include "wip/task.hpp"


class Test : public bliss::concurrent::Task
{
  public:
    Test(std::string msg) : bliss::concurrent::Task(), content(msg) {};
    virtual ~Test() {};

    virtual void operator()() {
      printf("tid %d: %s\n", omp_get_thread_num(), content.c_str());
    }

  protected:
    std::string content;

};



int main(int argc, char** argv) {

  {
    std::unique_ptr<Test> t1(new Test("task 1"));
    std::unique_ptr<Test> t2(new Test("task 2"));
    std::unique_ptr<Test> t3(new Test("task 3"));
    std::unique_ptr<Test> t4(new Test("task 4"));

    bliss::concurrent::SequentialRunner sr;
    //// test sequential.
    printf("serial runner\n");
    sr.addTask(std::move(t1));
    sr.addTask(std::move(t2));
    sr.addTask(std::move(t3));
    sr.addTask(std::move(t4));
    sr();
  }

  {
    std::unique_ptr<Test> t1(new Test("task 1"));
    std::unique_ptr<Test> t2(new Test("task 2"));
    std::unique_ptr<Test> t3(new Test("task 3"));
    std::unique_ptr<Test> t4(new Test("task 4"));

    bliss::concurrent::ReplicatedOMPRunner sr(3);
    //// test sequential.
    printf("replicated OMP runner\n");
    sr.addTask(std::move(t1));
    sr.addTask(std::move(t2));
    sr.addTask(std::move(t3));
    sr.addTask(std::move(t4));
    sr();
  }

  {
    std::unique_ptr<Test> t1(new Test("task 1"));
    std::unique_ptr<Test> t2(new Test("task 2"));
    std::unique_ptr<Test> t3(new Test("task 3"));
    std::unique_ptr<Test> t4(new Test("task 4"));

    bliss::concurrent::DynamicOMPRunner sr(3);
    //// test sequential.
    printf("dynamic OMP runner\n");
    sr.addTask(std::move(t1));
    sr.addTask(std::move(t2));
    sr.addTask(std::move(t3));
    sr.addTask(std::move(t4));
    sr.blockAdd();
    sr();
  }



//
//
//
//  //// test mixing different runners
//  printf("mpi runner with omp child runners\n");
//  bliss::concurrent::MPIRunner<false> mpir2(mpir.getComm());
//  if (mpir2.getId() == 0) mpir2.addTask(t1);
//  if (mpir2.getId() == 1) mpir2.addTask(sr);
//  mpir2.run();
//
//  if (mpir2.getId() == 0) mpir2.addTask(ompr);
//  mpir2.run();
//
//  if (mpir2.getId() == 0) mpir2.addTask(mompr);
//  mpir2.run();
//
//
//  if (mpir.getId() == 0) {
//    printf("omp runner with serial child runners\n");
//    ompr.addTask(sr);
//    ompr.run();
//  }




  //// test mpi
//  printf("test MPI runner\n");
//  bliss::concurrent::MPIRunner<false> mpir(argc, argv);
//  if (mpir.getId() == 0) mpir.addTask(t1);
//  if (mpir.getId() == 1) mpir.addTask(t2);
//  if (mpir.getId() == 2) mpir.addTask(t3);
//  if (mpir.getId() == 3) mpir.addTask(t4);
//  mpir.run();

}
