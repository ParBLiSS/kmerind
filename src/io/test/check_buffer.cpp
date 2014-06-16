/**
 * @file		check_buffer.cpp
 * @ingroup
 * @author	tpan
 * @brief
 * @details
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */


#include "io/Buffer.hpp"


int main(int argc, char** argv) {


  int x = 123456;


  bliss::io::Buffer<false> threadlocalBuffer(8192);

  // check insertion.
  for (int i = 0; i < 9000; i += 4) {
    if (! threadlocalBuffer.append(&x, sizeof(int))) {
      printf("canb't insert!\n");
    }
  }

  // check clear
  threadlocalBuffer.clear();

  // check isFull
  for (int i = 0; i < 9000; i += 4) {
    threadlocalBuffer.append(&x, sizeof(int));
    if (! threadlocalBuffer.isFull()) {
      printf("isFULL!\n");
    }
  }



  bliss::io::Buffer<true>  threadsafeBuffer(8192);







}
