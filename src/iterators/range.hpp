/**
 * range.hpp
 *
 *  Created on: Feb 18, 2014
 *      Author: tpan
 */

#ifndef RANGE_HPP_
#define RANGE_HPP_

#include <cassert>
#include <iostream>
#include <limits>

namespace bliss
{
  namespace iterator
  {
    /**
     * Range specified with offset, length, and overlap.  specific for 1D.
     */
    template<typename T>
    struct range
    {

        T block_start;   // starting position of range aligned to block
        T start; // offset from the beginning of block, if range is block aligned.
        T end;   // length of range, starting from offset_in_block.
//        T step;   // stride is the distance between each successive elements
        T overlap;  // amount of overlap at each end

        range(const T &_start, const T &_end, const T &_overlap = 0)
            : block_start(_start), start(_start), end(_end),
              overlap(_overlap)
        {
          assert(_start <= _end);
        }

        range(range<T> const &other)
            : block_start(other.block_start), start(other.start),
              end(other.end), overlap(other.overlap)
        {
        }

        range()
            : block_start(0), start(0), end(0), overlap(0)
        {
        }

        range<T>& operator=(range<T> const & other)
        {
          block_start = other.block_start;
          start = other.start;
          end = other.end;
          overlap = other.overlap;
          return *this;
        }

        bool operator==(const range<T> &other)
        {
          // same if the data range is identical and step is same.
          // not comparing overlap or block start.

          return (start == other.start) && (end == other.end);
        }

        // TODO: when needed: comparators
        // TODO: when needed: +/-
        // TODO: when needed: set operations

        /**
         *  block partitioning
         *  allows start and end to be negative. but require that end > start
         */
        static range<T> block_partition(const size_t &np, const size_t &pid,
                                        const T &start, const T &end,
                                        const T &_overlap = 0)
        {
          assert( start <= end);
          assert(_overlap >= 0);
          assert(pid < np);  // both non-negative.

          range<T> output(start, end, _overlap);

          if (np == 1)
            return output;

          T div = (end - start) / static_cast<T>(np);
          T rem = (end - start) % static_cast<T>(np);
          if (static_cast<T>(pid) < rem)
          {
            output.start += static_cast<T>(pid) * (div + 1);
            output.end = output.start + (div + 1) + _overlap;
          }
          else
          {
            output.start += static_cast<T>(pid) * div + rem;
            output.end = output.start + div + _overlap;
          }

          if (output.end > end)
            output.end = end;
          assert(output.start < output.end);

          output.block_start = output.start;
          return output;
        }

        static range<T> block_partition(const size_t &np, const size_t &pid,
                                 const range<T> &other)
        {
          return block_partition(np, pid, other.start, other.end, other.overlap);
        }

        range<T> block_partition(const size_t &np, const size_t &pid)
        {
          return block_partition(np, pid, this->start, this->end, this->overlap);
        }


        /**
         * align the range to page boundaries.
         */
        range<T> align_to_page(const size_t &page_size) const
        {
          assert(page_size > 0);

          range<T> output(*this);

          //printf("page size: %ld\n", static_cast<T>(page_size));

          // change start to align by page size.  extend range start.
          // note that if output.start is negative, it will put block_start at a bigger address than the start.
          output.block_start = (output.start / page_size) * page_size;

          if (output.block_start > output.start)  // only enters if start is negative.
          {

//            printf("block start: %ld\n", static_cast<size_t>(output.block_start));

            // if near lowest possible value, then we can't align further..  assert this situation.
            assert(
                (output.block_start - std::numeric_limits<T>::lowest()) > page_size);

            // deal with negative start position.
            output.block_start = output.block_start - page_size;

          }
          // leave end as is.

          return output;
        }

        bool is_page_aligned(const size_t &page_size) const
        {
          assert(page_size > 0);
          return (this->block_start % page_size) == 0;
        }

        /**
         * non-overlapping block partitions.
         *
         * takes into account the block size (e.g. page size) to force alignment
         *
         * uses overlap.  kernel handles bringing whole page in for the last (not full) page.
         *
         * result is the returned ranges are block aligned, ranges are mutually exclusive.
         */
        // deprecated
        //static range<T> block_partition_page_aligned(const T1& total, const T1& overlap, const T1& blocksize, const T2& np, const T2& pid)
        //{
        //  assert(total > 0);
        //  assert(blocksize > 0);
        //  assert(overlap >= 0);
        //  assert(np > 0);
        //  assert(pid >= 0 && pid < np);
        //
        //  RangeType<T1> output;
        //  output.overlap = overlap;
        //
        //  if (np == 1)
        //  {
        //    output.offset = 0;
        //    output.length = total;
        //    return output;
        //  }
        //
        //  // spread the number of blocks first.
        //  T1 nblock = total / blocksize;
        //
        //  T1 div = nblock / static_cast<T1>(np);
        //  T1 rem = nblock % static_cast<T1>(np);
        //  if (static_cast<T1>(pid) < rem)
        //  {
        //    output.offset = static_cast<T1>(pid) * (div + 1) * blocksize;
        //    output.length = (div + 1) * blocksize + overlap;
        //  }
        //  else
        //  {
        //    output.offset = (static_cast<T1>(pid) * div + rem) * blocksize;
        //    output.length = div * blocksize + overlap;
        //  }
        //
        //  assert(output.offset < total);
        //  if ((output.offset + output.length) >= total)
        //    output.length = total - output.offset;
        //
        //  return output;
        //}

    };

    template<typename T>
    std::ostream& operator<<(std::ostream& ost, const range<T>& r)
    {
      ost << "range: block@" << r.block_start << " [" << r.start << ":" << r.end << ") overlap " << r.overlap;
      return ost;
    }

  } /* namespace functional */
} /* namespace bliss */
#endif /* RANGE_HPP_ */
