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

#include "concurrent/mpi_runner.hpp"
#include "concurrent/omp_runner.hpp"
#include "concurrent/mixed_omp_runner.hpp"
#include "concurrent/sequential_runner.hpp"
#include "concurrent/runnable.hpp"


class Test : public bliss::concurrent::Runnable
{
  public:
    test(std::string msg) : content(msg) {};
    virtual ~test() {};

    virtual void run() {
      printf("%s\n", content.c_str());
    }

  protected:
    std::string content;

};



int main(int argc, char** argv) {

  Test t1("task 1");
  Test t2("task 2");
  Test t3("task 3");
  Test t4("task 4");


  //// test mpi
  printf("test MPI runner\n");
  bliss::concurrent::MPIRunner mpir(argc, argv);
  if (mpir.getId() == 0) mpir.addTask(t1);
  if (mpir.getId() == 1) mpir.addTask(t2);
  if (mpir.getId() == 2) mpir.addTask(t3);
  if (mpir.getId() == 3) mpir.addTask(t4);
  mpir.run();

  bliss::concurrent::SequentialRunner sr;
  bliss::concurrent::OMPRunner ompr(3);
  bliss::concurrent::MixedOMPRunner mompr(3);

  if (mpir.getId() == 0) {

    //// test sequential.
    printf("serial runner\n");
    sr.addTask(t1);
    sr.addTask(t2);
    sr.addTask(t3);
    sr.addTask(t4);
    sr.run();


    //// test omp next
    printf("mpi runner\n");
    ompr.addTask(t3);
    ompr.run();


    //// test mixed omp next
    printf("mixed mpi runner\n");
    mompr.addTask(t1);
    mompr.addTask(t2);
    mompr.addTask(t3);
    mompr.addTask(t4);
    mompr.addTask(t4);
    mompr.run;
  }

  //// test mixing different runners
  printf("mpi runner with omp child runners\n");
  bliss::concurrent::MPIRunner mpir2(mpir.getComm());
  if (mpir2.getId() == 0) mpir2.addTask(t1);
  if (mpir2.getId() == 1) mpir2.addTask(sr);
  mpir2.run();

  if (mpir2.getId() == 0) mpir2.addTask(ompr);
  mpir2.run();

  if (mpir2.getId() == 0) mpir2.addTask(mompr);
  mpir2.run();


  if (mpir.getId() == 0) {
    printf("omp runner with serial child runners\n");
    ompr.addTask(sr);
    ompr.run();
  }






}
