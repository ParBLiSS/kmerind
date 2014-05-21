/**
 * @file		iterator_functor.hpp
 * @ingroup
 * @author	tpan
 * @brief
 * @details
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef ITERATOR_FUNCTOR_HPP_
#define ITERATOR_FUNCTOR_HPP_

namespace bliss
{
  namespace iterator
  {

    /**
     * base class for functors.  not implemented.  check with Patrick.
     *
     * @class			bliss::iterator::IteratorFunctor
     * @brief
     * @details
     *
     */
    template <typename Iterator, typename Output, typename Range, typename Derived>
    class IteratorFunctor
    {
      public:
        IteratorFunctor();
        virtual ~IteratorFunctor();


        Iterator increment(const Iterator &iter, const Iterator &end, const Range &r);

        Output dereference(const Iterator &iter);

    };

  } /* namespace iterator */
} /* namespace bliss */

#endif /* ITERATOR_FUNCTOR_HPP_ */
