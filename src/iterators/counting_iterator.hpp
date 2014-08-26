/**
 * @file    counting_iterator.hpp
 * @ingroup bliss::iterators
 * @author  tpan
 * @brief   contains a counting iterator, which iterates over a range from start to end.
 * @details
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef COUNTING_ITERATOR_HPP_
#define COUNTING_ITERATOR_HPP_

#include <iterator>

#include "partition/range.hpp"


namespace bliss
{
  namespace iterator
  {

    /**
     * @class    bliss::iterator::CountingIterator
     * @brief    iterator that iterates from start to end in specified steps size
     * @details  counting iterator counts from the start to the end in regular steps.
     *           this is useful to establish index values for some other array.
     * @tparam RangeType  the type of range to iterate over.
     */
    template<typename RangeType>
    class CountingIterator :  public std::iterator<std::random_access_iterator_tag,
                                                   typename RangeType::ValueType,
                                                   decltype(std::declval<RangeType>().size())>
    {
      protected:
        /// output value type
        using T = RangeType::ValueType;

        /// the range within which to iterate over
        const RangeType range;

        /// the step for each iteration.
        const T step;

      public:

        /**
         * constructor
         * @param _range
         * @param step
         */
        CountingIterator(const RangeType& _range, const T &step) : range(_range), step(_step) {};

        virtual ~CountingIterator();
    };

  } /* namespace iterator */
} /* namespace bliss */

#endif /* COUNTING_ITERATOR_HPP_ */
