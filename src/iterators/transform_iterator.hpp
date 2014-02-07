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

  template<typename T1, typename T2>
  struct Transform {
      T2 operator()(T1 const & in) {

      };
  };

template<class Functor,
         class Base_Iterator >
class transform_iterator
    : public std::iterator<
           std::iterator_traits<Base_Iterator>::iterator_category,
          std::result_of<Functor(typename Base_Iterator::value_type)>::type,
          Base_Iterator::difference_type,
          std::add_pointer<std::result_of<Functor(typename Base_Iterator::value_type)>::type>::type,
          std::add_reference<std::result_of<Functor(typename Base_Iterator::value_type)>::type>::type

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
