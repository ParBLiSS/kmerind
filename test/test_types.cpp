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

#include <vector>

#include <iterators/transform_iterator.hpp>


//#include <ext/functor/Functor.h>

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


template<typename T1, typename T2>
struct Transform {
    T2 operator()(T1 const & in) {

    };
};


//template<typename FUNC>
//struct function_traits
//{
//    typedef FUNC type2;
//
//    typedef typename std::conditional< std::is_class<FUNC>::value && std::is_constructible<FUNC>::value,
//                                        decltype( std::declval< FUNC >().operator()(0) ),
//                                        std::conditional< std::is_function< typename FUNC >::value,
//                                                          decltype( FUNC(0) ),
//                                                          std::enable_if< std::is_member_function_pointer< decltype(typename FUNC) >::value,
//                                                                          decltype( FUNC(0) )
//                                                                          >
//                                                          >
//                                        >::type type;
//};


//template<typename Functor,
//         typename Base_Iterator >
//struct transform_iterator_functor
//  : public std::iterator<typename std::iterator_traits<Base_Iterator>::iterator_category,
//                         typename Functor::ResultType,
//                         typename std::iterator_traits<Base_Iterator>::difference_type,
//                         typename Functor::ResultType*,
//                         typename Functor::ResultType&
//                         >
//  {
//     typedef T value_type;
//  };
struct addConst {
    float operator()(int v) {
      return static_cast<float>(v) + 4.0;
    }
};

struct addConst2 {
    float foo(int v) {
      return static_cast<float>(v) + 20.0;
    }
};

float addC(int v) {
  return static_cast<float>(v) + 10.0;
}


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
  //static_assert(std::is_same<std::result_of< testFunc<int, float>(int) >::type, float>::value, "7.2");       // fails compile

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
  // it may not select the correct overloaded function, however


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
  printf("testStruct is function: %s\n", std::is_function< testStruct<int, float> >::value ? "yes" : "no");
  printf("testStruct is member function: %s\n", std::is_member_function_pointer< testStruct<int, float> >::value ? "yes" : "no");

  printf("testClass is class: %s\n", std::is_class< testClass<int, float> >::value ? "yes" : "no");
  printf("testClass is constructible: %s\n", std::is_constructible< testClass<int, float> >::value ? "yes" : "no");

  printf("testFunc is class: %s\n", std::is_class< decltype(testFunc<int, float> ) >::value ? "yes" : "no");
  printf("testFunc is constructible: %s\n", std::is_constructible< decltype(testFunc<int, float>) >::value ? "yes" : "no");
  printf("testFunc is function: %s\n", std::is_function< decltype(testFunc<int, float> ) >::value ? "yes" : "no");
  printf("testFunc is member function: %s\n", std::is_member_function_pointer< decltype(testFunc<int, float>) >::value ? "yes" : "no");

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


  // test instantiate some containers  This works.
  //transform_iterator< testStruct<int, float>, std::vector<int>::iterator, float> base;
  //transform_iterator< testStruct<int, float>, std::vector<int>::iterator> iter3;

  // for function pointers - decltype returns "T (* func)(T2)".  we know T2.  can we get T's type?
  // also, how to have compiler distinguish between function pointer and type?



  printf("type of testStruct is %s\n", typeid(testStruct<int, float>).name());
  printf("type of testSTruct output is %s\n", typeid(std::result_of< testClassWConstructor<int, float>(int)>::type).name());
  printf("type of testClass is %s\n", typeid(testClass<int, float>).name());
  printf("type of testFunc is %s\n", typeid(testFunc<int, float>).name());
  printf("type of testFunc decltype is %s\n", typeid( decltype(&testFunc<int, float>)).name());
  printf("type of testFunc output is %s\n", typeid( std::result_of< decltype(&testFunc<int, float>)(int)>::type ).name());


//  trans_iter< testFunc<int, float>, std::vector<int>::iterator> iter4;

  std::vector<int> data = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};

  typedef bliss::iterator::transform_iterator<addConst, std::vector<int>::iterator > iterType;
  addConst a;

  iterType iter(data.begin(), a);
  iterType end(data.end(), a);

  printf("created transform iterator\n");
  // iterType iter2 = iter;
  // iterType end = iter + n;
  for (; iter != end ; ++iter) {
    printf("%d -> %f, ", *(iter.getBaseIterator()), *(iter));
  }
  printf("\n");


  printf("type of addC is %s\n", typeid(decltype(&addConst2::foo)).name());
  typedef decltype(&addConst2::foo) cf_type;
  printf("type of testFunc decltype is %s\n", typeid( cf_type ).name());
  printf("type of testFunc output is %s\n", typeid( std::result_of< cf_type(addConst2, int)>::type ).name());
  printf("testFunc is class: %s\n", std::is_class< cf_type >::value ? "yes" : "no");
  printf("testFunc is constructible: %s\n", std::is_constructible< cf_type >::value ? "yes" : "no");
  printf("testFunc is function: %s\n", std::is_function< cf_type >::value ? "yes" : "no");
  printf("testFunc is member function: %s\n", std::is_member_function_pointer< cf_type >::value ? "yes" : "no");

  addConst2 ac2;
  cf_type fptr2 = &addConst2::foo;
  float f2 = (ac2.*fptr2)(4);

  printf("result of calling addC2 with 4 is %f\n", f2);


  typedef bliss::iterator::transform_iterator<decltype(&addConst2::foo), std::vector<int>::iterator > iterType2;

  iterType2 iter3(data.begin(), ac2, &addConst2::foo);
  iterType2 end3(data.end(), ac2, &addConst2::foo);

  printf("created transform iterator\n");
  // iterType iter3 = iter;
  // iterType end = iter + n;
  for (; iter3 != end3 ; ++iter3) {
    printf("%d -> %f, ", *(iter3.getBaseIterator()), *(iter3));

  }
  printf("\n");


  printf("type of addC is %s\n", typeid(decltype(&addC)).name());
  typedef decltype(&addC) f_type;
  typedef decltype(addC) f_type2;
  printf("type of testFunc decltype is %s\n", typeid( f_type ).name());
  printf("type of testFunc output is %s\n", typeid( std::result_of< f_type(int)>::type ).name());
  printf("testFunc is class: %s\n", std::is_class< f_type >::value ? "yes" : "no");
  printf("testFunc is constructible: %s\n", std::is_constructible< f_type >::value ? "yes" : "no");
  printf("testFunc is function: %s\n", std::is_function< f_type >::value ? "yes" : "no");
  printf("testFunc is member function: %s\n", std::is_member_function_pointer< f_type >::value ? "yes" : "no");
  printf("testFunc is class: %s\n", std::is_class< f_type2 >::value ? "yes" : "no");
  printf("testFunc is constructible: %s\n", std::is_constructible< f_type2 >::value ? "yes" : "no");
  printf("testFunc is function: %s\n", std::is_function< f_type2 >::value ? "yes" : "no");
  printf("testFunc is member function: %s\n", std::is_member_function_pointer< f_type2 >::value ? "yes" : "no");

  f_type fptr = &addC;
  float f = (*fptr)(4);

  printf("result of calling addC with 4 is %f\n", f);
//
//  printf("result %f\n", f);

//
  typedef bliss::iterator::transform_iterator< decltype(&addC), std::vector<int>::iterator> iter2Type;

  iter2Type iter2(data.begin(), &addC);
  iter2Type end2(data.end(), &addC);
  printf("created transform iterator 2\n");

  for (; iter2 != end2 ; ++iter2) {
    printf("%d -> %f, ", *(iter2.getBaseIterator()), *(iter2));

  }
  printf("\n");


}
