/**
 * range.hpp
 *
 *  Created on: Feb 18, 2014
 *      Author: tpan
 */

#ifndef RANGE_HPP_
#define RANGE_HPP_

#include <cassert>

namespace bliss
{
  namespace iterator
  {

    /**
     * Range specified with offset, length, and overlap.  specific for 1D.
     */
    template<typename T>
    struct range {

        T offset;   // starting position of range.
        T length;   // width of range
        T overlap;  // overlap with the next range
        T stride;   // stride is the distance between each successive elements

        range(T const& _offset,
              T const& _length,
              T const& _overlap = 0,
              T const& _stride = 0)
          : offset(_offset),
            length(_length),
            stride(_stride),
            overlap(_overlap)
        {}

        range(range<T> const &other)
          : offset(other.offset),
            length(other.length),
            stride(other.stride),
            overlap(other.overlap)
        {}

        range()
        : offset(0),
          length(0),
          stride(0),
          overlap(0)
        {}




        // TODO: comparators
        // TODO: +/-, extend, etc.

        /**
         *  block partitioning
         */
        static range<T> block_partition(T   const& total,
                                        int const& np,
                                        int const& pid,
                                        T   const& overlap = 0,
                                        T   const& stride = 0)
        {
          assert(total > 0);
          assert(overlap >= 0);
          assert(stride >= 0);
          assert(np > 0);
          assert(pid >= 0 && pid < np);

          range<T> output;
          output.overlap = overlap;
          output.stride = stride;

          if (np == 1)
          {
            output.offset = 0;
            output.length = total;
            return output;
          }

          T div = total / static_cast<T>(np);
          T rem = total % static_cast<T>(np);
          if (static_cast<T>(pid) < rem)
          {
            output.offset = static_cast<T>(pid) * (div + 1);
            output.length = div + 1 + overlap;
          }
          else
          {
            output.offset = static_cast<T>(pid) * div + rem;
            output.length = div + overlap;
          }

          assert(output.offset < total);
          if ((output.offset + output.length) >= total)
            output.length = total - output.offset;

          return output;
        }

        /**
         * align the range to page boundaries.
         */
        static range<T> align_to_page(range<T> const & input,
                                      T const &page_size,
                                      T const & total,
                                      int const & np)
        {
          range<T> output(input);

          // change offset to align by page size.  extend range length
          output.overlap = input.overlap;

          if (np == 1) {
            output.offset = input.offset;
            output.length = input.length;
            return output;
          }

          output.offset = (input.offset / page_size) * page_size;
          output.length = ((input.offset + input.length + page_size - 1) / page_size) * page_size - output.offset;

          assert(output.offset < total);
          if ((output.offset + output.length) >= total)
            output.length = total - output.offset;

          return output;
        }
    };

  } /* namespace functional */
} /* namespace bliss */
#endif /* RANGE_HPP_ */
