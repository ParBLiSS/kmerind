/**
 * test_types.cpp
 *
 *  Created on: Feb 5, 2014
 *      Author: tpan
 */

#include <cstdio>
#include <iostream>
#include <typeinfo>
#include <type_traits>
#include <utility>
#include <chrono>

#include <vector>
#include <functional>

#include <iterators/transform_iterator.hpp>


template<typename INPUT, typename OUTPUT>
struct testStruct {
    INPUT u;

    OUTPUT operator()(INPUT v) {
      return static_cast<OUTPUT>(v);
    }
    OUTPUT foo(INPUT v) {
      return static_cast<OUTPUT>(v);
    }
};

template<typename INPUT, typename OUTPUT>
class testClass {
  public:
    INPUT u;

    OUTPUT operator()(INPUT v) {
      return static_cast<OUTPUT>(v);
    }
    OUTPUT foo(INPUT v) {
      return static_cast<OUTPUT>(v);
    }

    static OUTPUT bar(INPUT v) {
      return static_cast<OUTPUT>(v);
    }
};

template<typename INPUT, typename OUTPUT>
struct testClassWConstructor {
    testClassWConstructor(INPUT v) {
    }

    OUTPUT operator()(INPUT v) {
      return static_cast<OUTPUT>(v);
    }
    OUTPUT foo(INPUT v) {
      return static_cast<OUTPUT>(v);
    }
};

template<typename INPUT, typename OUTPUT>
class testClassWDConstructor {
  public:
    testClassWDConstructor() {
    }

    OUTPUT operator()(INPUT v) {
      return static_cast<OUTPUT>(v);
    }
    OUTPUT foo(INPUT v) {
      return static_cast<OUTPUT>(v);
    }
};

template<typename INPUT, typename OUTPUT>
class testClassWPConstructor {
    testClassWPConstructor(INPUT v) {
    }

    OUTPUT operator()(INPUT v) {
      return static_cast<OUTPUT>(v);
    }
    OUTPUT foo(INPUT v) {
      return static_cast<OUTPUT>(v);
    }
};

template<typename INPUT, typename OUTPUT>
class testClassMultiOps {
  public:
    INPUT u;

    INPUT operator()() {
      return 0;
    }
    OUTPUT operator()(INPUT v) {
      return static_cast<OUTPUT>(v);
    }
    INPUT operator()(INPUT v, INPUT u) {
      return v + u;
    }
    OUTPUT foo(INPUT v) {
      return static_cast<OUTPUT>(v);
    }
};

template<typename INPUT, typename OUTPUT>
OUTPUT testFunc(INPUT v) {
  return static_cast<OUTPUT>(v);
};

template<typename INPUT, typename OUTPUT>
static OUTPUT testFuncStatic(INPUT v) {
  return static_cast<OUTPUT>(v);
};

template<typename T1, typename T2>
struct addConstFunctor {
    T2 operator()(T1 v) {
      return static_cast<T2>(v) + static_cast<T2>(1);
    }
};

template<typename T1, typename T2>
struct addConstMemberFunction {
    T2 foo(T1 v) {
      return static_cast<T2>(v) + static_cast<T2>(1);
    }
};

template<typename T1, typename T2>
T2 addConstFunction(T1 v) {
  return static_cast<T2>(v) + static_cast<T2>(1);
}

template<typename T1, typename T2>
struct squareFunctor {
    T2 operator()(T1 v) {
      return static_cast<T2>(v) * static_cast<T2>(v);
    }
};

template<typename T1, typename T2>
struct squareMemberFunction {
    T2 foo(T1 v) {
      return static_cast<T2>(v) * static_cast<T2>(v);
    }
};

template<typename T1, typename T2>
T2 squareFunction(T1 v) {
  return static_cast<T2>(v) * static_cast<T2>(v);
}


// this approaches is very similar to http://en.cppreference.com/w/cpp/types/result_of impl.
// all these are the same problems as with the traits classes, and is essentially the result_of implementation sans declval.
// this is unlikely to work for member functions, unless class type is specified as the first parameter of the function.
// declval is compile time only.
//template<typename T, typename A1>
//struct resultType {
//    static constexpr T* obj = nullptr;
//    static constexpr A1* par = nullptr;
//    typedef decltype((*obj)(*par)) return_type;
//};


int main(int argc, char* argv[]){

  // types to test:
  //  function pointer
  //  member function pointer
  //  static function pointer
  //  static member function pointer
  //  functor
  //  lambda function
  //  std::function
  //  std::bind output

  // NOTE:
  // & operator gets us the "address of".  e.g. address of function, == function pointer.
  //   in void (*f)(int) = somefunction,  f is a pointer that is being dereferenced when somefunction is assigned to it.
  //   so f == &somefunction
  // std::is_function operates on a func object, not a func pointer.
  // result_of uses declval internally (declval is compile time only), so need to be supplied with type info.
  //    further, result_of is implemented with decltype((F instance)(args instance)), requires F to be type of the function pointer.
  //    the object to pointer type conversion can be done in the constructor of the iterator and/or in the traits class.
  // std::function acts like a functor instance, so don't deference it.
  // lambda function also acts like a functor instance.  (it's own type)
  // implicit conversion from lambda function to function pointer (function pointer is a pointer type)


  // define some lambda functions
  std::function<double(double)> lf0    = [](double x){return 1;};
  auto                          lf1    = [](double x){return x;};
  auto                          lf2    = [](double x){return x*x;};
  double (* lf2_ptr)(double)           = lf2;



  // test result_of with: struct, class, struct with default constructor and defined constructor, class with constructors, member operator, member function, and standalone functions.
  static_assert(std::is_same<std::result_of< testStruct<int, float>(int) >::type, float>::value, "1.1");                                           // operator
  static_assert(std::is_same<std::result_of< decltype(&testStruct<int, float>::foo)(testStruct<int, float>, int) >::type, float>::value, "1.2");  // member function - use function pointer syntax
                                                                                                                                                // need decltype for the member function pointer
  static_assert(std::is_same<std::result_of< decltype(&testStruct<int, float>::operator())(testStruct<int, float>, int)>::type, float>::value, "1.3"); // member operator

  typedef decltype(&testStruct<int, float>::operator()) opType;
  static_assert(std::is_same<std::result_of< opType(testStruct<int, float>, int)>::type, float>::value, "1.3"); // member operator

  static_assert(std::is_same<std::result_of< testClass<int, float>(int)>::type, float>::value, "4.1");                                            // private constructor fails, public constructor compiles
  static_assert(std::is_same<std::result_of< decltype(&testClass<int, float>::operator())(testClass<int, float>, int)>::type, float>::value, "4.2");
  static_assert(std::is_same<std::result_of< decltype(&testClass<int, float>::foo)(testClass<int, float>, int) >::type, float>::value, "4.3");  // member function - use function pointer syntax
                                                                                                                                                // need decltype for the member function pointer
  static_assert(std::is_same<std::result_of< testClassWConstructor<int, float>(int)>::type, float>::value, "5.1");                               // private constructor fails, public constructor compiles
  static_assert(std::is_same<std::result_of< decltype(&testClassWConstructor<int, float>::operator())(testClassWConstructor<int, float>, int)>::type, float>::value, "5.2");
  static_assert(std::is_same<std::result_of< decltype(&testClassWConstructor<int, float>::foo)(testClassWConstructor<int, float>, int) >::type, float>::value, "5.3");  // member function - use function pointer syntax
                                                                                                                                                // need decltype for the member function pointer
  static_assert(std::is_same<std::result_of< testClassWDConstructor<int, float>(int)>::type, float>::value, "6.1");                               // private constructor fails, public constructor compiles
  static_assert(std::is_same<std::result_of< decltype(&testClassWDConstructor<int, float>::operator())(testClassWDConstructor<int, float>, int)>::type, float>::value, "6.2");
  static_assert(std::is_same<std::result_of< decltype(&testClassWDConstructor<int, float>::foo)(testClassWDConstructor<int, float>, int) >::type, float>::value, "6.3");  // member function - use function pointer syntax

  // private constructor. fails.
  //static_assert(std::is_same<std::result_of< testClassWPConstructor<int, float>(int)>::type, float>::value, "7.1");                               // private constructor fails, public constructor compiles
  static_assert(std::is_same<std::result_of< decltype(&testClassWPConstructor<int, float>::operator())(testClassWPConstructor<int, float>, int)>::type, float>::value, "6.2");
  static_assert(std::is_same<std::result_of< decltype(&testClassWPConstructor<int, float>::foo)(testClassWPConstructor<int, float>, int) >::type, float>::value, "6.3");  // member function - use function pointer syntax

  // need decltype for the member function pointer
  static_assert(std::is_same<std::result_of< decltype(&testFunc<int, float>)(int)>::type, float>::value, "8.1");

  // test static function
  static_assert(std::is_same<std::result_of< decltype(&testFuncStatic<int, float>)(int)>::type, float>::value, "9.1");
  static_assert(std::is_same<std::result_of< decltype(&testClass<int, float>::bar)(int) >::type, float>::value, "9.2");  // member function - use function pointer syntax

  // test lambda function
  static_assert(std::is_same<std::result_of< decltype(lf0)(double)>::type, double>::value, "8.2");
  static_assert(std::is_same<std::result_of< decltype(lf1)(double)>::type, double>::value, "8.3");  // member function - use function pointer syntax
  static_assert(std::is_same<std::result_of< decltype(lf2)(double)>::type, double>::value, "8.4");
  static_assert(std::is_same<std::result_of< decltype(lf2_ptr)(double)>::type, double>::value, "8.5");  // member function - use function pointer syntax



  // summary:  when use result_of, can only use it on functions or function pointers, or callable type (e.g. operator() defined).  with operator() defined,
  //           result_of appears to be using the function pointer appropriately, but requires public constructor if a class.
  //           struct is all public.  but best to be explicit and specify ::operator()
  // note that member function and regular functions require different parameters - member function also need type of the containing type.  so need to check.


  //////////// test using decltype without instantiating (actually, with "instantiation" using declval and default constructor)
  // decltype evaluates the entity or expression for its type.. . this is useful when dealing with types with different constructors.
  static_assert(std::is_same< decltype( std::declval< testStruct<int, float> >().operator()(0)), float >::value, "11.1");
  static_assert(std::is_same< decltype( std::declval< testStruct<int, float> >()(0)), float>::value, "11.2");
  static_assert(std::is_same< decltype( std::declval< testStruct<int, float> >().foo(0)), float>::value, "11.3");

  static_assert(std::is_same< decltype( std::declval< testClass<int, float> >().operator()(0)), float >::value, "14.1");
  static_assert(std::is_same< decltype( std::declval< testClass<int, float> >()(0)), float>::value, "14.2");
  static_assert(std::is_same< decltype( std::declval< testClass<int, float> >().foo(0)), float>::value, "14.3");

  static_assert(std::is_same< decltype( std::declval< testClassWConstructor<int, float> >().operator()(0)), float >::value, "15.1");
  static_assert(std::is_same< decltype( std::declval< testClassWConstructor<int, float> >()(0)), float>::value, "15.2");
  static_assert(std::is_same< decltype( std::declval< testClassWConstructor<int, float> >().foo(0)), float>::value, "15.3");

  static_assert(std::is_same< decltype( std::declval< testClassWDConstructor<int, float> >().operator()(0)), float >::value, "16.1");
  static_assert(std::is_same< decltype( std::declval< testClassWDConstructor<int, float> >()(0)), float>::value, "16.2");
  static_assert(std::is_same< decltype( std::declval< testClassWDConstructor<int, float> >().foo(0)), float>::value, "16.3");

  static_assert(std::is_same< decltype( std::declval< testClassWPConstructor<int, float> >().operator()(0)), float >::value, "18.1");
  static_assert(std::is_same< decltype( std::declval< testClassWPConstructor<int, float> >()(0)), float>::value, "18.2");
  static_assert(std::is_same< decltype( std::declval< testClassWPConstructor<int, float> >().foo(0)), float>::value, "18.3");

  static_assert(std::is_same< decltype( testFunc<int, float>(0)), float>::value, "17");


  // summary:
  // this removes the runtime instance requirement.  also, there is no requirement for visibility of constructor.
  // it may not select the correct overloaded function, however.  also, we still would need to instantiate the parameter.


  // test with constructed struct and classes
  testStruct<int, float> t1;
  static_assert(std::is_same<decltype(t1(1)), float>::value, "21.1");
  static_assert(std::is_same<decltype(t1.operator()(1)), float>::value, "21.2");
  static_assert(std::is_same<decltype(t1.foo(1)), float>::value, "21.3");
  testClass<int, float> t2;
  static_assert(std::is_same<decltype(t2(1)), float>::value, "22.1");
  static_assert(std::is_same<decltype(t2.operator()(1)), float>::value, "22.2");
  static_assert(std::is_same<decltype(t2.foo(1)), float>::value, "22.3");
  testClassWConstructor<int, float> tc2(1);
  static_assert(std::is_same<decltype(tc2(1)), float>::value, "24.1");
  static_assert(std::is_same<decltype(tc2.operator()(1)), float>::value, "24.2");
  static_assert(std::is_same<decltype(tc2.foo(1)), float>::value, "24.3");
  testClassWDConstructor<int, float> tdc2;
  static_assert(std::is_same<decltype(tdc2(1)), float>::value, "26.1");
  static_assert(std::is_same<decltype(tdc2.operator()(1)), float>::value, "26.2");
  static_assert(std::is_same<decltype(tdc2.foo(1)), float>::value, "26.3");

  // function pointer IS a constructed object.
  static_assert(std::is_same<decltype(testFunc<int, float>(1)), float>::value, "27.1");
  // summary
  // if we have an instance of a struct/class, or a function pointer, works.  but this may not be possible during type inference for templating.  hence the declval.
  // because of instance, cannot have private constructor.


  // overall summary:
  // use result_of avoids explicit instatiation of data types.

  // overloading tests.


  // assess the properties of a struct, member function, and function pointer.
  printf("testStruct is class: %s\n", std::is_class< testStruct<int, float> >::value ? "yes" : "no");
  printf("testStruct is constructible: %s\n", std::is_constructible< testStruct<int, float> >::value ? "yes" : "no");
  printf("testStruct is function: %s\n", std::is_function< testStruct<int, float> >::value ? "yes" : "no");
  printf("testStruct is member function: %s\n", std::is_member_function_pointer< testStruct<int, float> >::value ? "yes" : "no");

  // does not work when using decltype( std::declval<testStruct<int, float> >().operator() )  because whats' in decltype is not a pointer.
  // can use &testStruct<int, float>::operator().  can't use a bound member function address (i.e. an object's) to form a pointer to member function.
  printf("testStruct.operator() is class: %s\n", std::is_class< decltype(&testStruct<int, float>::operator() ) >::value ? "yes" : "no");
  printf("testStruct.operator() is constructible: %s\n", std::is_constructible< decltype(&testStruct<int, float>::operator()) >::value ? "yes" : "no");
  printf("testStruct.operator() is function: %s\n", std::is_function< decltype(&testStruct<int, float>::operator() ) >::value ? "yes" : "no");
  printf("testStruct.operator() is member function: %s\n", std::is_member_function_pointer< decltype(&testStruct<int, float>::operator()) >::value ? "yes" : "no");

  // does not work when using decltype( std::declval<testStruct<int, float> >().operator() )  because whats' in decltype is not a pointer.
  // can use &testStruct<int, float>::operator().  can't use a bound member function address (i.e. an object's) to form a pointer to member function.
  printf("testStruct.foo() is class: %s\n", std::is_class< decltype(&testStruct<int, float>::foo ) >::value ? "yes" : "no");
  printf("testStruct.foo() is constructible: %s\n", std::is_constructible< decltype(&testStruct<int, float>::foo) >::value ? "yes" : "no");
  printf("testStruct.foo() is function: %s\n", std::is_function< decltype(&testStruct<int, float>::foo ) >::value ? "yes" : "no");
  printf("testStruct.foo() is member function: %s\n", std::is_member_function_pointer< decltype(&testStruct<int, float>::foo) >::value ? "yes" : "no");

  printf("testFunc is class: %s\n", std::is_class< decltype(testFunc<int, float> ) >::value ? "yes" : "no");
  printf("testFunc is constructible: %s\n", std::is_constructible< decltype(testFunc<int, float>) >::value ? "yes" : "no");
  printf("testFunc is function: %s\n", std::is_function< decltype(testFunc<int, float> ) >::value ? "yes" : "no");
  printf("testFunc is member function: %s\n", std::is_member_function_pointer< decltype(testFunc<int, float>) >::value ? "yes" : "no");

  printf("&testFunc is class: %s\n", std::is_class< decltype(&testFunc<int, float> ) >::value ? "yes" : "no");
  printf("&testFunc is constructible: %s\n", std::is_constructible< decltype(&testFunc<int, float>) >::value ? "yes" : "no");
  printf("&testFunc is function: %s\n", std::is_function< decltype(&testFunc<int, float> ) >::value ? "yes" : "no");
  printf("&testFunc is member function: %s\n", std::is_member_function_pointer< decltype(&testFunc<int, float>) >::value ? "yes" : "no");

  printf("testFuncStatic is class: %s\n", std::is_class< decltype(testFuncStatic<int, float> ) >::value ? "yes" : "no");
  printf("testFuncStatic is constructible: %s\n", std::is_constructible< decltype(testFuncStatic<int, float>) >::value ? "yes" : "no");
  printf("testFuncStatic is function: %s\n", std::is_function< decltype(testFuncStatic<int, float> ) >::value ? "yes" : "no");
  printf("testFuncStatic is member function: %s\n", std::is_member_function_pointer< decltype(testFuncStatic<int, float>) >::value ? "yes" : "no");

  printf("&testFuncStatic is class: %s\n", std::is_class< decltype(&testFuncStatic<int, float> ) >::value ? "yes" : "no");
  printf("&testFuncStatic is constructible: %s\n", std::is_constructible< decltype(&testFuncStatic<int, float>) >::value ? "yes" : "no");
  printf("&testFuncStatic is function: %s\n", std::is_function< decltype(&testFuncStatic<int, float> ) >::value ? "yes" : "no");
  printf("&testFuncStatic is member function: %s\n", std::is_member_function_pointer< decltype(&testFuncStatic<int, float>) >::value ? "yes" : "no");


  printf("testStruct functor_trait class is (%s) %s\n",
         typeid(testStruct<int, double>).name(),
         typeid(bliss::iterator::func_traits<testStruct<int, double>, true, int >::class_type).name());
  printf("testStruct.operator() functor_trait class is (%s) %s\n",
         typeid(testStruct<int, long>).name(),
         typeid(bliss::iterator::func_traits<decltype(&testStruct<int, long>::operator()), false, int>::class_type).name());
  printf("testStruct.foo functor_trait class is (%s) %s\n",
         typeid(testStruct<int, int>).name(),
         typeid(bliss::iterator::func_traits<decltype(&testStruct<int, int>::foo), false, int>::class_type).name());
  printf("testFun functor_trait class is (nullptr) %s\n",
         typeid(bliss::iterator::func_traits<decltype(&testFunc<int, char>), false, int>::class_type).name());


  printf("testStruct functor_trait result is (d) %s\n", typeid(bliss::iterator::func_traits<testStruct<int, double>, true, int >::return_type).name());
  printf("testStruct.operator() functor_trait result is (l) %s\n", typeid(bliss::iterator::func_traits<decltype(&testStruct<int, long>::operator()), false, int>::return_type).name());
  printf("testStruct.foo functor_trait result is (i) %s\n", typeid(bliss::iterator::func_traits<decltype(&testStruct<int, int>::foo), false, int>::return_type).name());
  printf("testFun functor_trait result is (c) %s\n", typeid(bliss::iterator::func_traits<decltype(&testFunc<int, char>), false, int>::return_type).name());

  // should cause compile error
  // printf("testFunc.u functor_trait result is %s\n", typeid(bliss::iterator::func_traits<decltype(&testStruct<int, double>::u)>::return_type).name());

  // with overloaded operators.  standard decltype will fail.  should use exact operator.
  printf("testClassMultiOps.operator()() functor_trait result is (d) %s\n", typeid(bliss::iterator::func_traits<testClassMultiOps<int, double>, true>::return_type).name());
  printf("testClassMultiOps.operator()(int) functor_trait result is (d) %s\n", typeid(bliss::iterator::func_traits<testClassMultiOps<int, double>, true, int>::return_type).name());
  printf("testClassMultiOps.operator()(int,int) functor_trait result is (d) %s\n", typeid(bliss::iterator::func_traits<testClassMultiOps<int, double>, true, int, int>::return_type).name());


  printf("type of testStruct is %s\n", typeid(testStruct<int, float>).name());
  printf("type of testStruct output is %s\n", typeid(std::result_of< testClassWConstructor<int, float>(int)>::type).name());
  printf("type of testStruct.operator() output is %s\n", typeid(std::result_of< decltype(&testClassWConstructor<int, float>::operator())(testClassWConstructor<int, float>, int)>::type).name());
  printf("type of testStruct.foo output is %s\n", typeid(std::result_of< decltype(&testClassWConstructor<int, float>::foo)(testClassWConstructor<int, float>,int)>::type).name());
  printf("type of testClass is %s\n", typeid(testClass<int, float>).name());
  printf("type of testFunc is %s\n", typeid(testFunc<int, float>).name());
  printf("type of testFunc* is %s\n", typeid(testFunc<int, float>).name());

  std::cout << "type of &testFunc decltype is " << typeid( decltype(&testFunc<int, float>) ).name() << std::endl;
  printf("type of &testFunc output is %s\n", typeid( std::result_of< decltype(&testFunc<int, float>)(int)>::type ).name());

  printf("\n");

  std::chrono::high_resolution_clock::time_point time1, time2;
  std::chrono::duration<double> time_span;

  double gold_sum_squares = 0;
  double gold_sum = 0;
  // init the data.
  std::vector<int> data;
  for (int i = 0; i < 1000000; ++i) {
    data.push_back(i);
    gold_sum += double(i + 1);
    gold_sum_squares += double(i) * double(i);
  }
  printf("gold:  sum %f with add const\n", gold_sum);
  printf("gold:  sum squares %f \n", gold_sum_squares);

  double sum = 0;
  double sum_squares = 0;

  // CORRECTNESS and performance TESTING.  naive works, right?
  std::vector<int>::iterator baseIter = data.begin();
  std::vector<int>::iterator baseEnd = data.end();

  time1 = std::chrono::high_resolution_clock::now();
  for (; baseIter != baseEnd ; ++baseIter) {
    sum += static_cast<double>(*baseIter) + 1.0;
  }
  time2 = std::chrono::high_resolution_clock::now();
  time_span = std::chrono::duration_cast<std::chrono::duration<double>>(time2 - time1);
  printf("no functor:  sum %f with addConst elapsed time: %f s\n", sum, time_span.count());

  baseIter = data.begin();
  baseEnd = data.end();

  time1 = std::chrono::high_resolution_clock::now();
  for (; baseIter != baseEnd ; ++baseIter) {
    sum_squares += static_cast<double>(*baseIter) * static_cast<double>(*baseIter);
  }
  time2 = std::chrono::high_resolution_clock::now();
  time_span = std::chrono::duration_cast<std::chrono::duration<double>>(time2 - time1);
  printf("no functor:  sum_squares %f with square elapsed time: %f s\n", sum_squares, time_span.count());


  // test with functor
  addConstFunctor<int, double> a;
  bliss::iterator::transform_iterator<addConstFunctor<int, double>, std::vector<int>::iterator > iter1(data.begin(), a, a);
  bliss::iterator::transform_iterator<addConstFunctor<int, double>, std::vector<int>::iterator > end1(data.end(), a, a);

  sum = 0;
  time1 = std::chrono::high_resolution_clock::now();
  for (size_t i = 0; i < data.size(); ++i) {
    sum += iter1[i];
  }
  time2 = std::chrono::high_resolution_clock::now();
  time_span = std::chrono::duration_cast<std::chrono::duration<double>>(time2 - time1);
  printf("functor:  sum %f with addConst and [] elapsed time: %f s\n", sum, time_span.count());


  sum = 0;
  iter1 = bliss::iterator::transform_iterator<addConstFunctor<int, double>, std::vector<int>::iterator >(data.begin(), a, a);
  time1 = std::chrono::high_resolution_clock::now();
  for (; iter1 != end1 ; ++iter1) {
    sum += *iter1;
  }
  time2 = std::chrono::high_resolution_clock::now();
  time_span = std::chrono::duration_cast<std::chrono::duration<double>>(time2 - time1);
  printf("functor: sum %f with addConst elapsed time: %f s\n", sum, time_span.count());


  squareFunctor<int, double> b;
  bliss::iterator::transform_iterator<squareFunctor<int, double>, std::vector<int>::iterator > iter2(data.begin(), b, b);
  bliss::iterator::transform_iterator<squareFunctor<int, double>, std::vector<int>::iterator > end2(data.end(), b, b);

  sum_squares = 0;
  time1 = std::chrono::high_resolution_clock::now();
  for (; iter2 != end2 ; ++iter2) {
    sum_squares += *iter2;
  }
  time2 = std::chrono::high_resolution_clock::now();
  time_span = std::chrono::duration_cast<std::chrono::duration<double>>(time2 - time1);
  printf("functor:  sum_squares %f with square elapsed time: %f s\n", sum_squares, time_span.count());



  // function pointer
  addConstMemberFunction<int, double> c;
  bliss::iterator::transform_iterator<decltype(&addConstMemberFunction<int, double>::foo), std::vector<int>::iterator > iter3(data.begin(), &addConstMemberFunction<int, double>::foo, c);
  bliss::iterator::transform_iterator<decltype(&addConstMemberFunction<int, double>::foo), std::vector<int>::iterator > end3(data.end(), &addConstMemberFunction<int, double>::foo, c);

  sum = 0;
  time1 = std::chrono::high_resolution_clock::now();
  for (; iter3 != end3 ; ++iter3) {
    sum += *iter3;
  }
  time2 = std::chrono::high_resolution_clock::now();
  time_span = std::chrono::duration_cast<std::chrono::duration<double>>(time2 - time1);
  printf("member func:  SLOW sum %f with addConst elapsed time: %f s\n", sum, time_span.count());


  squareMemberFunction<int, double> d;
  bliss::iterator::transform_iterator<decltype(&squareMemberFunction<int, double>::foo), std::vector<int>::iterator > iter4(data.begin(),&squareMemberFunction<int, double>::foo, d);
  bliss::iterator::transform_iterator<decltype(&squareMemberFunction<int, double>::foo), std::vector<int>::iterator > end4(data.end(),&squareMemberFunction<int, double>::foo, d);

  sum_squares = 0;
  time1 = std::chrono::high_resolution_clock::now();
  for (; iter4 != end4 ; ++iter4) {
    sum_squares += *iter4;
  }
  time2 = std::chrono::high_resolution_clock::now();
  time_span = std::chrono::duration_cast<std::chrono::duration<double>>(time2 - time1);
  printf("mem func:  SLOW sum_squares %f with square elapsed time: %f s\n", sum_squares, time_span.count());


  // test function
  bliss::iterator::transform_iterator<decltype(&addConstFunction<int, double>), std::vector<int>::iterator > iter5(data.begin(), &addConstFunction<int, double>, 0);
  bliss::iterator::transform_iterator<decltype(&addConstFunction<int, double>), std::vector<int>::iterator > end5(data.end(), &addConstFunction<int, double>, 0);

  sum = 0;
  time1 = std::chrono::high_resolution_clock::now();
  for (; iter5 != end5 ; ++iter5) {
    sum += *iter5;
  }
  time2 = std::chrono::high_resolution_clock::now();
  time_span = std::chrono::duration_cast<std::chrono::duration<double>>(time2 - time1);
  printf("func:  sum %f with addConst elapsed time: %f s\n", sum, time_span.count());


  bliss::iterator::transform_iterator<decltype(&squareFunction<int, double>), std::vector<int>::iterator > iter6(data.begin(), &squareFunction<int, double>, 0);
  bliss::iterator::transform_iterator<decltype(&squareFunction<int, double>), std::vector<int>::iterator > end6(data.end(), &squareFunction<int, double>, 0);

  sum_squares = 0;
  time1 = std::chrono::high_resolution_clock::now();
  for (; iter6 != end6 ; ++iter6) {
    sum_squares += *iter6;
  }
  time2 = std::chrono::high_resolution_clock::now();
  time_span = std::chrono::duration_cast<std::chrono::duration<double>>(time2 - time1);
  printf("func:  sum_squares %f with square elapsed time: %f s\n", sum_squares, time_span.count());



  // member function pointer

//  typedef decltype(&addConst2::foo) cf_type;
//  addConst2 ac2;
//  cf_type fptr2 = &addConst2::foo;
//  float f2 = (ac2.*fptr2)(4);

//  f_type fptr = &addC;
//  float f = (*fptr)(4);



}
