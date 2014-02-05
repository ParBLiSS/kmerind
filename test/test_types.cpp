/**
 * test_types.cpp
 *
 *  Created on: Feb 5, 2014
 *      Author: tpan
 */

#include <cstdio>
#include <iostream>
#include <typeinfo>

template<typename T>
void print() {
  printf("unknown\n");
}

template<>
void print<int>() {
  printf("int\n");
}
template<>
void print<float>() {
  printf("float\n");
}

template<>
void print<double>() {
  printf("double\n");
}

template<typename INPUT, typename OUTPUT>
struct test {
    OUTPUT operator()(INPUT v) {
      return static_cast<OUTPUT>(v);
    }
};

int main(int argc, char* argv[]){

  test<int, float> t1;
  std::cout << t1(1) << std::endl;

  typedef decltype(t1(1)) tt;

  print<tt>();

  typedef struct test<int, double> TEST;
  TEST t2;
  typedef decltype(t2(1)) tt2;
  printf("type name: %s\n", typeid(tt2).name());
}
