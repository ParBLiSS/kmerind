/**
 * @file    range.hpp
 * @ingroup iterators
 * @author  Tony Pan
 * @brief   Generic representation of an interval.
 * @details Represents an interval with start, end, and overlap length.
 *   Also contains a block_start to mark beginning of a underlying data block, such as paging size.
 *
 * @copyright BLISS: Copyright (c) Georgia Institute of Technology.  All Rights Preserved.
 *
 * TODO add Licence
 */

#ifndef RANGE_HPP_
#define RANGE_HPP_

#include <cassert>
#include <iostream>
#include <limits>

namespace bliss
{
  /**
   * @namespace iterator
   */
  namespace iterator
  {
    /**
     * @class range
     * @brief Range specified with offsets and overlap.  specific for 1D.
     * @details
     *
     * @tparam T  data type used for the start and end offsets and overlap.
     */
    template<typename T>
    struct range
    {
        /**
         * @var   block_start
         * @brief starting position of range aligned to underlying block boundary
         */
        T block_start;
        /**
         * @var   start
         * @brief starting position of a range in absolute coordinates
         */
        T start;
        /**
         * @var     end
         * @brief   end position of a range in absolute coordinates.
         * @details End points to 1 position past the last element in the range
         */
        T end;
        /**
         * @var   overlap
         * @brief amount of overlap between adjacent ranges.
         */
        T overlap;

        /**
         * @brief   construct directly from start and end offsets and overlap
         * @details _start should be less than or equal to _end
         *
         * @param[in] _start    starting position of range in absolute coordinates
         * @param[in] _end      ending position of range in absoluate coordinates.
         * @param[in] _overlap  amount of overlap between adjacent ranges.  optional
         */
        range(const T &_start, const T &_end, const T &_overlap = 0)
            : block_start(_start), start(_start), end(_end),
              overlap(_overlap)
        {
          assert(_start <= _end);
        }

        /**
         * @brief   copy construct from the field values of another range
         * @param[in] other   the range object to copy from
         */
        range(const range<T> &other)
            : block_start(other.block_start), start(other.start),
              end(other.end), overlap(other.overlap)
        {
        }

        /**
         * @brief   default constructor.  construct an empty range, with start and end initialized to 0.
         *
         */
        range()
            : block_start(0), start(0), end(0), overlap(0)
        {
        }


        /**
         * @brief assignment operator.  copy the field values from the operand range.
         *
         * @param[in] other   the range object to copy from
         * @return  the calling range object, with field values copied from the operand range object
         */
        range<T>& operator=(const range<T> & other)
        {
          block_start = other.block_start;
          start = other.start;
          end = other.end;
          overlap = other.overlap;
          return *this;
        }

        /**
         * @brief equals operator.  compares 2 ranges' start and end positions only.
         *
         * @param[in] other   The range to compare to
         * @return  true if 2 ranges are same.  false otherwise.
         */
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
         * @brief   static function.  block partitioning of a range
         * @details Given the number of partition, the partition element desired, and a start and end range, deterministically compute the subrange.
         *    start end end are allowed to be negative, but end >= start is required.
         *
         * @param[in] np        number of partitions
         * @param[in] pid       id of specific partition desired
         * @param[in] start     the starting offset of the range to be partitioned.
         * @param[in] end       the ending offset of the range to be partitioned
         * @param[in] _overlap  the overlap between partitions
         * @return              computed subrange
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

        /**
         * @brief static function.  block partitioning of a range
         * @details Given the number of partition, the partition element desired, and a start and end range, deterministically compute the subrange.
         *    start end end are allowed to be negative, but end >= start is required.
         *
         * @param[in] np        number of partitions
         * @param[in] pid       id of specific partition desired
         * @param[in] other     range object to be partitioned
         * @return              computed subrange
         */
        static range<T> block_partition(const size_t &np, const size_t &pid,
                                 const range<T> &other)
        {
          return block_partition(np, pid, other.start, other.end, other.overlap);
        }

        /**
         * @brief static function.  block partitioning of a range
         * @details Given the number of partition, the partition element desired, deterministically compute a subrange for the current range object.
         *    start end end are allowed to be negative, but end >= start is required.
         *
         * @param[in] np        number of partitions
         * @param[in] pid       id of specific partition desired
         * @return              computed subrange
         */
        range<T> block_partition(const size_t &np, const size_t &pid)
        {
          return block_partition(np, pid, this->start, this->end, this->overlap);
        }


        /**
         * @brief   align the range to underlying block boundaries, e.g. disk page size
         * @details range is aligned to underlying block boundaries by moving the block_start variable back towards minimum
         *    if range start is too close to the data type's minimum, then assertion is thrown.
         *
         * @param[in] page_size   the size of the underlying block.
         * @return                updated range
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

        /**
         * @brief     check to see if the range has been aligned to underlying block boundary.
         *
         * @param[in] page_size   the size of the underlying block.
         * @return    true if range is block aligned, false otherwise.
         */
        bool is_page_aligned(const size_t &page_size) const
        {
          assert(page_size > 0);
          return (this->block_start % page_size) == 0;
        }

        /*
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

    /**
     * @brief << operator to write out range object's fields.
     * @param[in/out] ost   output stream to which the content is directed.
     * @param[in]     r     range object to write out
     * @return              output stream object
     */
    template<typename T>
    std::ostream& operator<<(std::ostream& ost, const range<T>& r)
    {
      ost << "range: block@" << r.block_start << " [" << r.start << ":" << r.end << ") overlap " << r.overlap;
      return ost;
    }

  } /* namespace functional */
} /* namespace bliss */
#endif /* RANGE_HPP_ */
