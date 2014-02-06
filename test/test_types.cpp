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



//template<typename T>
//void print() {
//  printf("unknown\n");
//}
//
//template<>
//void print<int>() {
//  printf("int\n");
//}
//template<>
//void print<float>() {
//  printf("float\n");
//}
//
//template<>
//void print<double>() {
//  printf("double\n");
//}

template<typename INPUT, typename OUTPUT>
struct testStruct {
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
    OUTPUT operator()(INPUT v) {
      return static_cast<OUTPUT>(v);
    }
    OUTPUT foo(INPUT v) {
      return static_cast<OUTPUT>(v);
    }
};


template<typename INPUT, typename OUTPUT>
struct testStructWDConstructor {
    testStructWDConstructor() {
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
struct testStructWConstructor {
    testStructWConstructor(INPUT v) {
    }

    OUTPUT operator()(INPUT v) {
      return static_cast<OUTPUT>(v);
    }
    OUTPUT foo(INPUT v) {
      return static_cast<OUTPUT>(v);
    }
};

template<typename INPUT, typename OUTPUT>
class testClassWConstructor {
  public:
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
OUTPUT testFunc(INPUT v) {
  return static_cast<OUTPUT>(v);
};


// 3 types:  global function, member function, functor
//template<typename FUNC, typename OT>
//struct function_traits {
//    typedef OT type;
//};
//
//
//template<typename FUNC>
//struct function_traits<FUNC, std::enable_if< std::is_class<FUNC>::value && std::is_constructible<FUNC>::value,
//                                             decltype( std::declval< FUNC >().operator()(0) ) > >
//{
//    typedef decltype( std::declval< FUNC >().operator()(0)) type;
//};
//
//template<typename FUNC>
//struct function_traits<FUNC, std::enable_if< std::is_member_function_pointer<FUNC>::value || std::is_function<FUNC>::value,
//                                             decltype( FUNC(0)) > >
//{
//    typedef decltype( FUNC(0)) type;
//};

template<typename FUNC>
struct function_traits
{
    typedef typename std::conditional< std::is_class<FUNC>::value && std::is_constructible<FUNC>::value,
                                        decltype( std::declval< FUNC >().operator()(0) ),
                                        std::enable_if< std::is_member_function_pointer<FUNC>::value || std::is_function<FUNC>::value,
                                                        decltype( FUNC(0) )
                                                      >
                                      >::type type;
//    typedef typename std::enable_if< std::is_class<FUNC>::value && std::is_constructible<FUNC>::value,
//                                        decltype( std::declval< FUNC >().operator()(0) ) > type2;
//
//    typedef typename std::enable_if< std::is_member_function_pointer<FUNC>::value || std::is_function<FUNC>::value,
//                                        decltype( FUNC(0) ) >::type type2;
};




//template<class FUNC, class T>
//class container {
//    // define bar using types of FUNC
//    typename std::enable_if<std::is_constructible<FUNC>::value,
//                            decltype(std::declval<FUNC>())>::type
//    callme(T v) {
//      return FUNC(v);
//    }
//
//    typename std::enable_if<!std::is_constructible<FUNC>::value,
//                             std::result_of<FUNC(T)>::type>::type
//    callme(T v)
//    {
//      return FUNC(v);
//    }
//};

//template<class FUNC, class BASE_VECTOR>
//class myvec : public std::vector<
//  typename std::enable_if<std::is_constructible<FUNC>::value,
//                            decltype(std::declval<FUNC>())>::type >
//{
//};
//
//template<class FUNC, class BASE_VECTOR>
//class myvec : public std::vector<
//  typename std::enable_if<!std::is_constructible<FUNC>::value,
//  std::result_of<FUNC(BASE_VECTOR::value_type)>::type>::type >
//{
//};



int main(int argc, char* argv[]){


  // test result_of with: struct, class, struct with default constructor and defined constructor, class with constructors, member operator, member function, and standalone functions.
  static_assert(std::is_same<std::result_of< testStruct<int, float>(int)>::type, float>::value, "1.1");                                           // operator
  //static_assert(std::is_same<std::result_of< testStruct<int, float>::foo(int)>::type, float>::value, "2");     // bad syntax.  testStruct<int, float>::foo fail compilation
  //static_assert(std::is_same<std::result_of< &(testStruct<int, float>::foo)(testStruct<int, float>, int)>::type, float>::value, "2");           // bad syntax, function pointer...
  static_assert(std::is_same<std::result_of< decltype(&testStruct<int, float>::foo)(testStruct<int, float>, int) >::type, float>::value, "1.2");  // member function - use function pointer syntax
                                                                                                                                                // need decltype for the member function pointer
  static_assert(std::is_same<std::result_of< decltype(&testStruct<int, float>::operator())(testStruct<int, float>, int)>::type, float>::value, "1.3"); // member operator

  typedef decltype(&testStruct<int, float>::operator()) opType;
  static_assert(std::is_same<std::result_of< opType(testStruct<int, float>, int)>::type, float>::value, "1.3"); // member operator


  static_assert(std::is_same<std::result_of< testStructWConstructor<int, float>(int)>::type, float>::value, "2.1");                               // operator?
  static_assert(std::is_same<std::result_of< decltype(&testStructWConstructor<int, float>::foo)(testStructWConstructor<int, float>, int) >::type, float>::value, "2.2");  // member function - use function pointer syntax
                                                                                                                                                // need decltype for the member function pointer
  static_assert(std::is_same<std::result_of< decltype(&testStructWConstructor<int, float>::operator())(testStructWConstructor<int, float>, int)>::type, float>::value, "2.3");  // member operator

  static_assert(std::is_same<std::result_of< testStructWDConstructor<int, float>(int)>::type, float>::value, "3.1");                               // operator, right?
  static_assert(std::is_same<std::result_of< decltype(&testStructWDConstructor<int, float>::foo)(testStructWDConstructor<int, float>, int) >::type, float>::value, "3.2");  // member function - use function pointer syntax
                                                                                                                                                // need decltype for the member function pointer
  static_assert(std::is_same<std::result_of< decltype(&testStructWDConstructor<int, float>::operator())(testStructWDConstructor<int, float>, int)>::type, float>::value, "3.3");  // member operator


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
                                                                                                                                                // need decltype for the member function pointer

  static_assert(std::is_same<std::result_of< decltype(&testFunc<int, float>)(int)>::type, float>::value, "7.1");
  // static_assert(std::is_same<std::result_of< testFunc<int, float>(int) >::type, float>::value, "7.2");       // fails compile

 // static_assert(std::is_same<std::result_of< testClassWPConstructor<int, float>(int)>::type, float>::value, "8.1");                               // private constructor fails, public constructor compiles
  static_assert(std::is_same<std::result_of< decltype(&testClassWPConstructor<int, float>::operator())(testClassWPConstructor<int, float>, int)>::type, float>::value, "6.2");
  static_assert(std::is_same<std::result_of< decltype(&testClassWPConstructor<int, float>::foo)(testClassWPConstructor<int, float>, int) >::type, float>::value, "6.3");  // member function - use function pointer syntax
                                                                                                                                                // need decltype for the member function pointer


  // summary:  when use result_of, can only use it on functions or function pointers, or callable type (e.g. operator() defined).  with operator() defined,
  //           result_of appears to be using the function pointer appropriately, but requires public constructor if a class.
  //           struct is more forgiving.  but best to be explicit and specify ::operator()
  // note that member function and regular functions require different parameters - member function also need type of the containing type.  so need to check.

   // test using decltype without instantiating (actually, with "instantiation" using declval and default constructor)
   // decltype evaluates the entity or expression for its type.. . this is useful when dealing with types with different constructors.

  static_assert(std::is_same< decltype( std::declval< testStruct<int, float> >().operator()(0)), float >::value, "11.1");
  static_assert(std::is_same< decltype( std::declval< testStruct<int, float> >()(0)), float>::value, "11.2");
  static_assert(std::is_same< decltype( std::declval< testStruct<int, float> >().foo(0)), float>::value, "11.3");

  static_assert(std::is_same< decltype( std::declval< testStructWConstructor<int, float> >().operator()(0)), float >::value, "12.1");
  static_assert(std::is_same< decltype( std::declval< testStructWConstructor<int, float> >()(0)), float>::value, "12.2");
  static_assert(std::is_same< decltype( std::declval< testStructWConstructor<int, float> >().foo(0)), float>::value, "12.3");

  static_assert(std::is_same< decltype( std::declval< testStructWDConstructor<int, float> >().operator()(0)), float >::value, "13.1");
  static_assert(std::is_same< decltype( std::declval< testStructWDConstructor<int, float> >()(0)), float>::value, "13.2");
  static_assert(std::is_same< decltype( std::declval< testStructWDConstructor<int, float> >().foo(0)), float>::value, "13.3");

  static_assert(std::is_same< decltype( std::declval< testClass<int, float> >().operator()(0)), float >::value, "14.1");
  static_assert(std::is_same< decltype( std::declval< testClass<int, float> >()(0)), float>::value, "14.2");
  static_assert(std::is_same< decltype( std::declval< testClass<int, float> >().foo(0)), float>::value, "14.3");

  static_assert(std::is_same< decltype( std::declval< testClassWConstructor<int, float> >().operator()(0)), float >::value, "15.1");
  static_assert(std::is_same< decltype( std::declval< testClassWConstructor<int, float> >()(0)), float>::value, "15.2");
  static_assert(std::is_same< decltype( std::declval< testClassWConstructor<int, float> >().foo(0)), float>::value, "15.3");

  static_assert(std::is_same< decltype( std::declval< testClassWDConstructor<int, float> >().operator()(0)), float >::value, "16.1");
  static_assert(std::is_same< decltype( std::declval< testClassWDConstructor<int, float> >()(0)), float>::value, "16.2");
  static_assert(std::is_same< decltype( std::declval< testClassWDConstructor<int, float> >().foo(0)), float>::value, "16.3");

  static_assert(std::is_same< decltype( testFunc<int, float>(0)), float>::value, "17");

  static_assert(std::is_same< decltype( std::declval< testClassWPConstructor<int, float> >().operator()(0)), float >::value, "18.1");
  static_assert(std::is_same< decltype( std::declval< testClassWPConstructor<int, float> >()(0)), float>::value, "18.2");
  static_assert(std::is_same< decltype( std::declval< testClassWPConstructor<int, float> >().foo(0)), float>::value, "18.3");

  // summary:
  // this removes the runtime instance requirement.  also, there is no requirement for visibility of constructor.
  // however, there is still a difference between


  // test with constructed struct and classes
  testStruct<int, float> t1;
  static_assert(std::is_same<decltype(t1(1)), float>::value, "21.1");
  static_assert(std::is_same<decltype(t1.operator()(1)), float>::value, "21.2");
  static_assert(std::is_same<decltype(t1.foo(1)), float>::value, "21.3");
  testClass<int, float> t2;
  static_assert(std::is_same<decltype(t2(1)), float>::value, "22.1");
  static_assert(std::is_same<decltype(t2.operator()(1)), float>::value, "22.2");
  static_assert(std::is_same<decltype(t2.foo(1)), float>::value, "22.3");
  testStructWConstructor<int, float> tc1(1);
  static_assert(std::is_same<decltype(tc1(1)), float>::value, "23.1");
  static_assert(std::is_same<decltype(tc1.operator()(1)), float>::value, "23.2");
  static_assert(std::is_same<decltype(tc1.foo(1)), float>::value, "23.3");
  testClassWConstructor<int, float> tc2(1);
  static_assert(std::is_same<decltype(tc2(1)), float>::value, "24.1");
  static_assert(std::is_same<decltype(tc2.operator()(1)), float>::value, "24.2");
  static_assert(std::is_same<decltype(tc2.foo(1)), float>::value, "24.3");
  testStructWDConstructor<int, float> tdc1;
  static_assert(std::is_same<decltype(tdc1(1)), float>::value, "25.1");
  static_assert(std::is_same<decltype(tdc1.operator()(1)), float>::value, "25.2");
  static_assert(std::is_same<decltype(tdc1.foo(1)), float>::value, "25.3");
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
  // use declval to create instance if constructible, then use decltype.  if function pointer, just use decltype.  "call" with some dummy value.

  printf("testStruct is class: %s\n", std::is_class< testStruct<int, float> >::value ? "yes" : "no");
  printf("testStruct is constructible: %s\n", std::is_constructible< testStruct<int, float> >::value ? "yes" : "no");

  printf("testClass is class: %s\n", std::is_class< testClass<int, float> >::value ? "yes" : "no");
  printf("testClass is constructible: %s\n", std::is_constructible< testClass<int, float> >::value ? "yes" : "no");

  printf("testFunc is class: %s\n", std::is_class< decltype(testFunc<int, float> ) >::value ? "yes" : "no");
  printf("testFunc is constructible: %s\n", std::is_constructible< decltype(testFunc<int, float>) >::value ? "yes" : "no");

//  printf("testFunc is function: %s\n", std::is_function< decltype(testFunc<int, float>) >::value ? "yes" : "no");
//  printf("testFunc is memberfunction pointer: %s\n", std::is_constructible< typename testClass<int, float> >::value ? "yes" : "no");
//

//  typedef function_traits<  testStruct<int, float> > tt;
  //static_assert(std::is_same< function_traits< testFunc<int, float> >::type, float>::value, "28");
//  typename testStruct<int, float> x;
//  static_assert(std::is_same< function_traits< x.operator() >::type, float>::value, "29");
//  static_assert(std::is_same< function_traits< testStruct<int, float> >::type, float>::value, "30");

}
