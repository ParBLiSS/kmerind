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

#include <string.h>

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

    OUTPUT operator()() {
      return 0;
    }
    OUTPUT operator()(INPUT v) {
      return static_cast<OUTPUT>(v);
    }
    OUTPUT operator()(INPUT v, INPUT u) {
      return v + u;
    }
    OUTPUT foo(INPUT v) {
      return static_cast<OUTPUT>(v);
    }
};

template<typename INPUT, typename OUTPUT>
class testClassMultiFuncs {
  public:
    INPUT u;

    OUTPUT operator()(INPUT v) {
      return static_cast<OUTPUT>(v);
    }
    OUTPUT foo(INPUT v, INPUT u) {
      return v + u;
    }
    OUTPUT foo(INPUT v) {
      return static_cast<OUTPUT>(v);
    }
    OUTPUT foo() {
      return 0;
    }
};


template<typename INPUT, typename OUTPUT>
OUTPUT testFunc(INPUT v) {
  return static_cast<OUTPUT>(v);
};

template<typename INPUT, typename OUTPUT>
class testClassConstRef {
  public:
    OUTPUT testFunc(INPUT v) {
      return static_cast<OUTPUT>(v);
    }
    OUTPUT testConstFunc(const INPUT v) {
      return static_cast<OUTPUT>(v);
    }
    OUTPUT testRefFunc(INPUT & v) {
      return static_cast<OUTPUT>(v);
    }
    OUTPUT testConstRefFunc(const INPUT & v) {
      return static_cast<OUTPUT>(v);
    }
    OUTPUT testFuncConst(INPUT v) const{
      return static_cast<OUTPUT>(v);
    }
    OUTPUT testConstFuncConst(const INPUT v) const{
      return static_cast<OUTPUT>(v);
    }
    OUTPUT testRefFuncConst(INPUT & v) const {
      return static_cast<OUTPUT>(v);
    }
    OUTPUT testConstRefFuncConst(const INPUT & v) const {
      return static_cast<OUTPUT>(v);
    }
    OUTPUT& testRefFuncRef(INPUT & v) {
      return reinterpret_cast<OUTPUT>(v);
    }
    OUTPUT& testRefFuncConstRef(INPUT & v) const {
      return reinterpret_cast<OUTPUT>(v);
    }
    OUTPUT& testConstRefFuncConstRef(const INPUT & v) const {
      return reinterpret_cast<OUTPUT>(v);
    }
};

template<typename INPUT, typename OUTPUT>
struct testFunctorXXXX {
    OUTPUT operator()(INPUT v) {
      return static_cast<OUTPUT>(v);
    }
};
template<typename INPUT, typename OUTPUT>
struct testFunctorXXCX {
    OUTPUT operator()(INPUT v) const{
      return static_cast<OUTPUT>(v);
    }
};
template<typename INPUT, typename OUTPUT>
struct testFunctorCXXX {
    OUTPUT operator()(const INPUT v) {
      return static_cast<OUTPUT>(v);
    }
};
template<typename INPUT, typename OUTPUT>
struct testFunctorCXCX {
    OUTPUT operator()(const INPUT v) const{
      return static_cast<OUTPUT>(v);
    }
};
template<typename INPUT, typename OUTPUT>
struct testFunctorXRXX {
    OUTPUT operator()(INPUT & v) {
      return static_cast<OUTPUT>(v);
    }
};
template<typename INPUT, typename OUTPUT>
struct testFunctorXRXR {
    OUTPUT& operator()(INPUT & v) {
      return reinterpret_cast<OUTPUT>(v);
    }
};
template<typename INPUT, typename OUTPUT>
struct testFunctorXRCX {
    OUTPUT operator()(INPUT & v) const {
      return static_cast<OUTPUT>(v);
    }
};
template<typename INPUT, typename OUTPUT>
struct testFunctorXRCR {
    OUTPUT& operator()(INPUT & v) const {
      return reinterpret_cast<OUTPUT>(v);
    }
};
template<typename INPUT, typename OUTPUT>
struct testFunctorCRXX {
    OUTPUT operator()(const INPUT & v) {
      return static_cast<OUTPUT>(v);
    }
};
template<typename INPUT, typename OUTPUT>
struct testFunctorCRCX {
    OUTPUT operator()(const INPUT & v) const {
      return static_cast<OUTPUT>(v);
    }
};
template<typename INPUT, typename OUTPUT>
struct testFunctorCRCR {
    OUTPUT& operator()(const INPUT & v) const {
      return reinterpret_cast<OUTPUT>(v);
    }
};



template<typename INPUT, typename OUTPUT>
OUTPUT testOverloadedFunc(INPUT v) {
  return static_cast<OUTPUT>(v);
};
template<typename INPUT, typename OUTPUT>
OUTPUT testOverloadedFunc(INPUT v, INPUT v2) {
  return static_cast<OUTPUT>(v);
};
template<typename INPUT, typename OUTPUT>
OUTPUT testOverloadedFunc() {
  return static_cast<OUTPUT>(1);
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
// also, need to specify T to account for overloading.
// declval is compile time only.
// declval does not require public default constructors.
// declval on member function pointer looks like???
//template<typename T, typename A1>
//struct resultType {
//    static constexpr T* obj = nullptr;
//    static constexpr A1* par = nullptr;
//    typedef decltype((*obj)(*par)) return_type;
//};



template<typename F, typename C, typename R, typename... Args>
void testFunctionTraits(std::string name) {

  printf("NAME: %s\n", name.c_str());
  printf("\ttypeid.name %s\n", typeid(F).name());

  printf("\tis class: %s\n", std::is_class< F >::value ? "yes" : "no");
  printf("\tis constructible: %s\n", std::is_constructible< F >::value ? "yes" : "no");
  printf("\tis function: %s\n", std::is_function< F >::value ? "yes" : "no");
  printf("\tis member function: %s\n", std::is_member_function_pointer< F >::value ? "yes" : "no");

  printf("\tfunctor_trait class type is %s (== %s ? %s)\n",
         typeid(typename bliss::iterator::func_traits<F, std::is_class<F>::value, Args... >::class_type).name(),
         typeid(C).name(),
         strcmp(typeid(typename bliss::iterator::func_traits<F, std::is_class<F>::value, Args... >::class_type).name(),
                typeid(C).name()) == 0 ? "yes" : "no");

  printf("\tfunctor_trait result type is %s (== %s ? %s)\n",
         typeid(typename bliss::iterator::func_traits<F, std::is_class<F>::value, Args... >::return_type).name(),
         typeid(R).name(),
         strcmp(typeid(typename bliss::iterator::func_traits<F, std::is_class<F>::value, Args... >::return_type).name(),
                typeid(R).name()) == 0 ? "yes" : "no");

  // does not work for functor as this. type not defined.  don't know why. ....
//  printf("\tresult_of result type is %s (== %s ?)\n",
//         typeid(typename std::result_of<F(Args...)>::type).name(),
//         typeid(R).name());

  printf("\n");

}

template<typename F, typename V>
void testTransformIterator(F f,
                           typename bliss::iterator::func_traits<F,
                                                                 std::is_class<F>::value,
                                                                 typename std::iterator_traits<typename V::iterator>::value_type
                                                                 >::class_type c,
                           V& data, std::string name, std::string func_name)
{
  // test with functor
  typedef bliss::iterator::transform_iterator<F, typename V::iterator> t_iter_type;
  t_iter_type iter1(data.begin(), f, c);
  t_iter_type end1(data.end(), f, c);

  std::chrono::high_resolution_clock::time_point time1, time2;
  std::chrono::duration<double> time_span;


  double sum = 0;
  time1 = std::chrono::high_resolution_clock::now();
  for (size_t i = 0; i < data.size(); ++i) {
    sum += iter1[i];
  }
  time2 = std::chrono::high_resolution_clock::now();
  time_span = std::chrono::duration_cast<std::chrono::duration<double>>(time2 - time1);
  std::cout << name << " " << func_name << ":  sum " << sum << " and [] elapsed time: " << time_span.count() << " s" << std::endl;


  sum = 0;
  iter1 = t_iter_type(data.begin(), f, c);
  time1 = std::chrono::high_resolution_clock::now();
  for (; iter1 != end1 ; ++iter1) {
    sum += *iter1;
  }
  time2 = std::chrono::high_resolution_clock::now();
  time_span = std::chrono::duration_cast<std::chrono::duration<double>>(time2 - time1);
  std::cout << name << " " << func_name << ":  sum " << sum << "  elapsed time: " << time_span.count() << " s" << std::endl;

}


int main(int argc, char* argv[]){

  // types to test:
  //  done function pointer
  //  done member function pointer
  //  done static function pointer
  //  done static member function pointer
  //  done functor
  //  done lambda function
  //  done std::function (as lambda function)
  //  done constness and references for input, reference for output
  //  done const member function
  //  std::bind output

  // NOTE:
  // & operator gets us the "address of".  e.g. address of function, == function pointer.
  //   in void (*f)(int) = somefunction,  f is a pointer that is being dereferenced when somefunction is assigned to it.
  //   so f == &somefunction
  // std::is_function operates on a func object, not a func pointer.
  // result_of uses declval internally (declval is compile time only), so need to be supplied with type info (decltype on function pointer).
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

  // private constructor. fails. but declval version works (used in result_of)
  //static_assert(std::is_same<std::result_of< testClassWPConstructor<int, float>(int)>::type, float>::value, "7.1");                               // private constructor fails, public constructor compiles

  // need decltype for the member function pointer
  static_assert(std::is_same<std::result_of< decltype(&testFunc<int, float>)(int)>::type, float>::value, "8.1");

  // test lambda function
  static_assert(std::is_same<std::result_of< decltype(lf0)(double)>::type, double>::value, "8.2");
  static_assert(std::is_same<std::result_of< decltype(lf1)(double)>::type, double>::value, "8.3");  // member function - use function pointer syntax
  static_assert(std::is_same<std::result_of< decltype(lf2)(double)>::type, double>::value, "8.4");
  static_assert(std::is_same<std::result_of< decltype(lf2_ptr)(double)>::type, double>::value, "8.5");  // member function - use function pointer syntax



  // summary:  when use result_of, can only use it on functions or function pointers, or callable type (e.g. operator() defined).  with operator() defined,
  //           result_of appears to be using the function pointer, but requires public constructor if a class.
  //           struct is all public.  but best to be explicit and specify ::operator()
  // note that member function and regular functions require different parameters - member function also need type of the containing type.  so need to check.


  //////////// test using decltype without instantiating (actually, with "instantiation" using declval and default constructor)
  // decltype evaluates the entity or expression for its type.. . this is useful when dealing with types with different constructors.
  // this is similar to the internals of std::result_of, except the declval is called on the class, then an explicit function is used wtih decltype.
  static_assert(std::is_same< decltype( std::declval< testStruct<int, float> >().operator()(0)), float >::value, "11.1");
  static_assert(std::is_same< decltype( std::declval< testStruct<int, float> >()(0)), float>::value, "11.2");
  static_assert(std::is_same< decltype( std::declval< testStruct<int, float> >().foo(0)), float>::value, "11.3");

  // doesn't work for member function this way.
  //static_assert(std::is_same< decltype( std::declval< &testStruct<int, float>::foo >()(0)), float>::value, "11.4");

  // testFunc is already a pointer.
  static_assert(std::is_same< decltype( testFunc<int, float>(0)), float>::value, "17");
  // this method does not seem to work for overloaded operators.

  // summary:
  // this removes the runtime instance requirement.  also, there is no requirement for visibility of constructor.
  // it may not select the correct overloaded function, however.  also, we still would need to instantiate the parameter.


  // test with constructed struct and classes.  this does actual instantiation.
  testStruct<int, float> t1;
  static_assert(std::is_same<decltype(t1(1)), float>::value, "21.1");
  static_assert(std::is_same<decltype(t1.operator()(1)), float>::value, "21.2");
  static_assert(std::is_same<decltype(t1.foo(1)), float>::value, "21.3");

  // function pointer IS a constructed object.
  static_assert(std::is_same<decltype(testFunc<int, float>(1)), float>::value, "27.1");
  // summary
  // if we have an instance of a struct/class, or a function pointer, works.  but this may not be possible during type inference for templating.  hence the declval.
  // because of instance, cannot have private constructor.


  // overall summary:
  // use result_of avoids explicit instatiation of data types.


  // testing:
  // assess the properties of a struct, member function, and function pointer.
  testFunctionTraits<testStruct<int, float>, testStruct<int, float>, float, int>("testStruct<int, float>");
  testFunctionTraits<decltype(&testStruct<int, float>::operator()), testStruct<int, float>, float, int>("decltype(&testStruct<int, float>::operator())");
  testFunctionTraits<decltype(&testStruct<int, float>::foo), testStruct<int, float>, float, int>("decltype(&testStruct<int, float>::foo)");
  testFunctionTraits<decltype(&testFunc<int, float>), std::nullptr_t, float, int>("decltype(&testFunc<int, float>)");
  testFunctionTraits<decltype(&testFuncStatic<int, float>), std::nullptr_t, float, int>("decltype(&testFuncStatic<int, float>) - static");
  // static member function is same as static function.
  testFunctionTraits<decltype(&testClass<int, float>::bar), std::nullptr_t, float, int>("decltype(&testClass<int, float>::bar) - static");

  // class with overloaded operators.  standard decltype will fail with operator name and functions.  functors are okay.
  testFunctionTraits<testClassMultiOps<int, double>, testClassMultiOps<int, double>, int>("testClassMultiOps<int, double>())");
  testFunctionTraits<testClassMultiOps<int, double>, testClassMultiOps<int, double>, double, int>("testClassMultiOps<int, double>(int)");
  testFunctionTraits<testClassMultiOps<int, double>, testClassMultiOps<int, double>, int, int, int>("testClassMultiOps<int, double>(int, int))");
  // alternatively use typedef Ret ([class::]*blah)(args...); then use blah as the type, or use static cast.
  testFunctionTraits<decltype(static_cast<float (*)()>(&testOverloadedFunc<int, float>)), std::nullptr_t, float>("&testOverloadedFunc<int, double>()) - static cast");
  testFunctionTraits<decltype(static_cast<float (*)(int)>(&testOverloadedFunc<int, float>)), std::nullptr_t, float, int>("&testOverloadedFunc<int, double>(int)) - static cast");
  testFunctionTraits<decltype(static_cast<float (*)(int, int)>(&testOverloadedFunc<int, float>)), std::nullptr_t, float, int, int>("&testOverloadedFunc<int, double>(int, int)) - static cast");

  // overloaded member function pointer.
  testFunctionTraits<decltype(static_cast<float (testClassMultiFuncs<int, float>::*)()>(&testClassMultiFuncs<int, float>::foo)), testClassMultiFuncs<int, float>, float>("&testClassMultiFuncs<int, float>::foo()) - static cast");
  testFunctionTraits<decltype(static_cast<float (testClassMultiFuncs<int, float>::*)(int)>(&testClassMultiFuncs<int, float>::foo)), testClassMultiFuncs<int, float>, float, int>("&testClassMultiFuncs<int, float>::foo(int)) - static cast");
  testFunctionTraits<decltype(static_cast<float (testClassMultiFuncs<int, float>::*)(int, int)>(&testClassMultiFuncs<int, float>::foo)), testClassMultiFuncs<int, float>, float, int, int>("&testClassMultiFuncs<int, float>::foo(int, int)) - static cast");

  // does not work when using decltype( std::declval<testStruct<int, float> >().operator() )  because whats' in decltype is not a pointer.
  //printf("decltype( std::declval<testStruct<int, float> >().operator() ) is class: %s\n", std::is_class< decltype( std::declval<testStruct<int, float> >().operator() ) >::value ? "yes" : "no");

  // can't use a bound member function address (i.e. an object's) to form a pointer to member function.
  // printf("testFunc.u functor_trait result is %s\n", typeid(bliss::iterator::func_traits<decltype(&testStruct<int, double>::u)>::return_type).name());

  // lambda functions.
  testFunctionTraits<decltype(lf0), std::function<double(double)>, double, double>("decltype(lf0)");
  testFunctionTraits<decltype(lf1), std::nullptr_t, double, double>("decltype(lf1)");
  testFunctionTraits<decltype(lf2), std::nullptr_t, double, double>("decltype(lf2)");
  testFunctionTraits<decltype(lf2_ptr), std::nullptr_t, double, double>("decltype(lf2_ptr)");

  // const and ref functions
  testFunctionTraits<decltype(&testClassConstRef<int, float>::testFunc),                 testClassConstRef<int, float>, float, int>("decltype(&testClassConstRef<int, float>::testFunc)");
  testFunctionTraits<decltype(&testClassConstRef<int, float>::testConstFunc),            testClassConstRef<int, float>, float, int>("decltype(&testClassConstRef<int, float>::testConstFunc)");
  testFunctionTraits<decltype(&testClassConstRef<int, float>::testRefFunc),              testClassConstRef<int, float>, float, int>("decltype(&testClassConstRef<int, float>::testRefFunc)");
  testFunctionTraits<decltype(&testClassConstRef<int, float>::testConstRefFunc),         testClassConstRef<int, float>, float, int>("decltype(&testClassConstRef<int, float>::testConstRefFunc)");
  testFunctionTraits<decltype(&testClassConstRef<int, float>::testFuncConst),            testClassConstRef<int, float>, float, int>("decltype(&testClassConstRef<int, float>::testFuncConst)");
  testFunctionTraits<decltype(&testClassConstRef<int, float>::testConstFuncConst),       testClassConstRef<int, float>, float, int>("decltype(&testClassConstRef<int, float>::testConstFuncConst)");
  testFunctionTraits<decltype(&testClassConstRef<int, float>::testRefFuncConst),         testClassConstRef<int, float>, float, int>("decltype(&testClassConstRef<int, float>::testRefFuncConst)");
  testFunctionTraits<decltype(&testClassConstRef<int, float>::testConstRefFuncConst),    testClassConstRef<int, float>, float, int>("decltype(&testClassConstRef<int, float>::testConstRefFuncConst)");
  testFunctionTraits<decltype(&testClassConstRef<int, float>::testRefFuncRef),           testClassConstRef<int, float>, int,   int>("decltype(&testClassConstRef<int, float>::testRefFuncRef)");
  testFunctionTraits<decltype(&testClassConstRef<int, float>::testRefFuncConstRef),      testClassConstRef<int, float>, int,   int>("decltype(&testClassConstRef<int, float>::testRefFuncConstRef)");
  testFunctionTraits<decltype(&testClassConstRef<int, float>::testConstRefFuncConstRef), testClassConstRef<int, float>, int,   int>("decltype(&testClassConstRef<int, float>::testConstRefFuncConstRef)");

  // const and ref functors
  testFunctionTraits<testFunctorXXXX<int, float>, testFunctorXXXX<int, float>, float, int>("testFunctorXXXX<int, float>");
  testFunctionTraits<testFunctorXXCX<int, float>, testFunctorXXCX<int, float>, float, int>("testFunctorXXCX<int, float>");
  testFunctionTraits<testFunctorCXXX<int, float>, testFunctorCXXX<int, float>, float, int>("testFunctorCXXX<int, float>");
  testFunctionTraits<testFunctorCXCX<int, float>, testFunctorCXCX<int, float>, float, int>("testFunctorCXCX<int, float>");
  testFunctionTraits<testFunctorXRXX<int, float>, testFunctorXRXX<int, float>, float, int>("testFunctorXRXX<int, float>");
  testFunctionTraits<testFunctorXRXR<int, float>, testFunctorXRXR<int, float>, float, int>("testFunctorXRXR<int, float>");
  testFunctionTraits<testFunctorXRCX<int, float>, testFunctorXRCX<int, float>, float, int>("testFunctorXRCX<int, float>");
  testFunctionTraits<testFunctorXRCR<int, float>, testFunctorXRCR<int, float>, float, int>("testFunctorXRCR<int, float>");
  testFunctionTraits<testFunctorCRXX<int, float>, testFunctorCRXX<int, float>, float, int>("testFunctorCRXX<int, float>");
  testFunctionTraits<testFunctorCRCX<int, float>, testFunctorCRCX<int, float>, float, int>("testFunctorCRCX<int, float>");
  testFunctionTraits<testFunctorCRCR<int, float>, testFunctorCRCR<int, float>, float, int>("testFunctorCRCR<int, float>");


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
  testTransformIterator(a, a, data, "functor", "addConstFunctor<int, double>");


  squareFunctor<int, double> b;
  testTransformIterator(b, b, data, "functor", "squareFunctor<int, double>");



  // function pointer
  addConstMemberFunction<int, double> c;
  testTransformIterator(&addConstMemberFunction<int, double>::foo, c, data, "mem func", "&addConstMemberFunction<int, double>::foo");

  squareMemberFunction<int, double> d;
  testTransformIterator(&squareMemberFunction<int, double>::foo, d, data, "mem func", "&squareMemberFunction<int, double>::foo");

  // test function
  testTransformIterator(&addConstFunction<int, double>, nullptr, data, "func", "&addConstFunction<int, double>");

  testTransformIterator(&squareFunction<int, double>, nullptr, data, "func", "&squareFunction<int, double>");


  // member function pointer format.

//  typedef decltype(&addConst2::foo) cf_type;
//  addConst2 ac2;
//  cf_type fptr2 = &addConst2::foo;
//  float f2 = (ac2.*fptr2)(4);

//  f_type fptr = &addC;
//  float f = (*fptr)(4);



}