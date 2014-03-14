/**
 * @file		test_threads.cpp
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
#include <thread>


void fileio() {
  // open the file and create a fastq iterator
    // each dereference has pointers.  choices:
      // a. directly use.  mmap may be jumping around
      // b. preload in fastqloader.  jumping around in memory
      // c. get the pointers and copy the data locally.

  // for each compute thread, accessing next element on file iterator returns a preloaded
  // data object.  with locking to update pointer.

  // return data object pointing to the local copy.
  printf("file io\n");

}

void compute(int i) {
   // do something.

  printf("compute %d\n", i);

}

void networkwrite() {
  // instantiate bins (m per proc)

  // atomic add to bin. update bin count

  // if bin count is at limit, MPI send length, then send data
  //
  // final write
    // flush: sends length
      // if length > 0, send data
      // else length = 0, don't send anything.  return to waiting.

  // send termination signal (send empty message).
  printf("network write\n");

}

void networkread() {
  // initiate an array of senders, mark all as 1.

  // iprobe for a message.  if empty message, then don't need to listen on that node anymore

  // if has length,
    // if length > 0, then receive data and do something
  // else return to probe.

  printf("network read\n")
}


int main(int argc, char** argv) {

  std::cout << "starting threads" << std::endl;

  std::thread fileio_t(fileio);
  std::thread networkwrite_t(networkwrite);
  std::thread networkread_t(networkread);

  std::cout << "threads started" << std::endl;
  // do some work using openmp

  int threadcount = 8;

#pragma omp parallel for
  for (int i = 0; i < 8; ++i) {
    compute(i);
  }

  std::cout << "computation done" << std::endl;

  fileio_t.join();
  networkwrite_t.join();
  networkread_t.join();

  std::cout << "threads done" << std::endl;

  return 0;
}
