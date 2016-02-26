/**
 * @file    icc-autodeduce.cpp
 * @ingroup
 * @author  tpan
 * @brief
 * @details
 *
 * Copyright (c) 2016 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */

#include <cstdio>

namespace tester {
  template <typename T, unsigned long l>
  void test(T (&x)[l]) {
    printf ("good");

  }
}

int main(int argc, char** argv) {

  unsigned char x[20];


  tester::test(x);

}
