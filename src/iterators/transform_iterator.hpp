/**
 * transform_iterator.hpp
 *
 *  Created on: Feb 5, 2014
 *      Author: tpan
 */

#ifndef TRANSFORM_ITERATOR_HPP_
#define TRANSFORM_ITERATOR_HPP_

#include <iterator>

namespace bliss
{



namespace iterator
{


  template <typename func>
  void foo(func f)
  {
    std::vector<int> vec;
    for (auto i : vec)
    {
      f(i);
    }
  }


  template<typename T1, typename T2>
  struct Transform {
      T2 operator()(T1 const & in) {

      };
  };

  template <typename FuncType, typename InputType>
  struct GetFuncOutputType
  {
      static constexpr InputType t;
      typename decltype(FuncType(t)) value_type;
  };



template<class Functor,
  class Base_Iterator >
class transform_iterator
    : public std::iterator<std::random_access_iterator_tag,
          typename decltype(Functor()),
          Base_Iterator::difference_type,
          typename iterator_traits<_Iterator>::pointer,
          typename iterator_traits<_Iterator>::reference>

{

}

    void main()
    {
        int bar;
        int baz;
        int max = -1;
        foo([&max](int x) {return if (x > max) max = x;});

    }

} // iterator
} // bliss
#endif /* TRANSFORM_ITERATOR_HPP_ */
