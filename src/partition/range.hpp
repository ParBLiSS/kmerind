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
 *
 * TODO: tests for intersection, union, complement, and shift operators.
 */

#ifndef RANGE_HPP_
#define RANGE_HPP_

#include <cassert>
#include <iostream>
#include <limits>
#include <algorithm>

#include <partition/partitioner.hpp>

namespace bliss
{
  /**
   * @namespace iterator
   */
  namespace partition
  {
    /**
     * @class range
     * @brief Range specified with offsets and overlap.  specific for 1D.  overlap is on END side only, and is included in the END.
     * @details
     *
     * @tparam T  data type used for the start and end offsets and overlap.
     */
    template<typename T>
    struct range
    {
        typedef T ValueType;

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
         * @brief   end position of a range in absolute coordinates.  DOES include overlap.
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
          assert(_overlap >= 0);
        }

        /**
         * @brief   copy construct from the field values of another range
         * @param[in] other   the range object to copy from
         */
        range(const range<T> &other)
            : block_start(other.block_start), start(other.start),
              end(other.end), overlap(other.overlap)
        {}

        /**
         * @brief   default constructor.  construct an empty range, with start and end initialized to 0.
         *
         */
        range()
            : block_start(0), start(0), end(0), overlap(0)
        {}


        /**
         * @brief assignment operator.  copy the field values from the operand range.
         *
         * @param[in] other   the range object to copy from
         * @return  the calling range object, with field values copied from the operand range object
         */
        range<T>& operator =(const range<T> & other)
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
        bool operator ==(const range<T> &other) const
        {
          // same if the data range is identical and step is same.
          // not comparing overlap or block start.

          return (start == other.start) && (end == other.end);
        }

        // TODO: when needed: comparators  - what does that mean?
        // TODO: when needed: +/-

        /**
         * @brief union of range (as |=)  NOTE: result could include previously unincluded ranges
         * @param other
         * @return
         */
        range<T>& operator |=(const range<T>& other) {
          block_start =
                start = std::min(start,       other.start);
          end   =       std::max(end,         other.end);
          overlap =     std::max(overlap,     other.overlap);

          return *this;
        }
        /**
         * @brief union of range (as |)  NOTE: result could include previously unincluded ranges
         * @param other
         * @return
         */
        range<T> operator |(const range<T>& other) const
        {
          range<T> output(*this);
          output |= other;
          return output;
        }

        /**
         * @brief intersection of range (as &=)
         * @param other
         * @return
         */
        range<T>& operator &=(const range<T>& other)
        {
          start =       std::max(start,       other.start);
          end   =       std::min(end,         other.end);
          overlap =     std::max(overlap,     other.overlap);

          // in case the ranges do not intersect
          block_start =
                start = std::min(start, end);

          return *this;
        }
        /**
         * @brief intersection of range (as |)
         * @param other
         * @return
         */
        range<T> operator &(const range<T>& other) const
        {
          range<T> output(*this);
          output &= other;
          return output;
        }
        /**
         * @brief complement operation (as -=).  order matters
         * @param other
         * @return
         */
        range<T>& operator -=(const range<T>& other)
        {
          // cases: other.start < start < end  : output should be other.start <-> other.start
          //        start < other.start < end  : output should be       start <-> other.start
          //        start < end < other.start  : output should be       start <-> end
          block_start =
                start = std::min(start,       other.start);
                end   = std::min(end,         other.start);
          return *this;
        }
        /**
         * @brief complement operation (as -). order matters
         * @param other
         * @return
         */
        range<T> operator -(const range<T>& other) const
        {
          range<T> output(*this);
          output -= other;
          return output;
        }

        /**
         * @brief right shift operation (as >>=).
         * @param amount
         * @return
         */
        range<T>& operator >>=(const T& amount)
        {
          start += amount;
          end   += amount;
          block_start = start;
          return *this;
        }
        /**
         * @brief right shift operation (as >>).
         * @param amount
         * @return
         */
        range<T> operator >>(const T& amount) const
        {
          range<T> output(*this);
          output >>= amount;
          return output;
        }
        /**
         * @brief left shift operation (as <<=).
         * @param amount
         * @return
         */
        range<T>& operator <<=(const T& amount)
        {
          start -= amount;
          end   -= amount;
          block_start = start;
          return *this;
        }
        /**
         * @brief left shift operation (as <<).
         * @param amount
         * @return
         */
        range<T> operator <<(const T& amount)
        {
          range<T> output(*this);
          output <<= amount;
          return output;
        }



        /**
         * @brief   static function.  block partitioning of a range
         * @details Given the number of partition, the partition element desired, and a start and end range, deterministically compute the subrange.
         *    start and end are allowed to be negative, but end >= start is required.
         *    OVERLAP IS PART OF END.
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
          std::cerr << "range::block_partition  TESTING ONLY.  to be deprecated." << std::endl;

          BlockPartitioner<T> partitioner(range<T>(start, end, _overlap), np);
          return partitioner.getNext(pid);

//          assert( start <= end);
//          assert(_overlap >= 0);
//          assert(pid >= 0);
//          assert(pid < np);  // both non-negative.
//
//          range<T> output(start, end, _overlap);
//
//          if (np == 1)
//            return output;
//
//          size_t length = end - start;
//
//          size_t div = length / np;
//          size_t rem = length % np;
//          if (pid < rem)
//          {
//            output.start += static_cast<T>(pid * (div + 1));
//            output.end = output.start + static_cast<T>(div + 1) + _overlap;
//          }
//          else
//          {
//            output.start += static_cast<T>(pid * div + rem);
//            output.end = output.start + static_cast<T>(div) + _overlap;
//          }
//
//          // last entry.  no overlap.
//          if (output.end > end) {
//            output.end = end;
//            output.overlap = 0;
//          }
//          assert(output.start <= output.end);
//
//          output.block_start = output.start;
//          return output;
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
        range<T> block_partition(const size_t &np, const size_t &pid) const
        {
          return block_partition(np, pid, this->start, this->end, this->overlap);
        }

//        /**
//
//         * @brief   static function.  page aligned, non-overlapping block partitioning of a range
//         * @details Given the number of partition, the partition element desired, and a start and end range, deterministically compute the subrange.
//         *    start and end are allowed to be negative, but end >= start is required.
//         *    OVERLAP IS PART OF END.
//         *
//         * @param[in] np        number of partitions
//         * @param[in] pid       id of specific partition desired
//         * @param[in] page_size the size of the underlying block.
//         * @param[in] start     the starting offset of the range to be partitioned.
//         * @param[in] end       the ending offset of the range to be partitioned
//         * @param[in] _overlap  the overlap between partitions
//         * @return              computed subrange
//         */
//        static range<T> paged_block_partition(const size_t &np, const size_t &pid, const size_t& page_size,
//                                              const T &start, const T &end,
//                                              const T &_overlap = 0)
//        {
//          assert( start <= end);
//          assert(_overlap >= 0);
//          assert(pid >= 0)
//          assert(pid < np);  // both non-negative.
//          assert(page_size > 0);
//
//          range<T> output(start, end, _overlap);
//
//          if (np == 1)
//          {
//            return output;
//          }
//
//
//          // construct in units of page.
//          range<T> pages;
//          pages.start = start / page_size;
//          pages.end = (end + page_size - 1) / page_size;
//          pages.overlap = (_overlap + page_size - 1) / page_size;
//
//          // now partition the pages
//          output = block_partition(np, pid, pages);
//
//          // expand the new output
//          output.start *= page_size;
//          output.end *= page_size;
//          output.overlap = 0;
//
//          // bound by original start and end.
//          output.start = std::max(output.start, start);
//          output.end = std::min(output.end, end);
//
//          // last check
//          assert(output.start <= output.end);
//          output.block_start = output.start;
//
//          return output;
//        }
//
//
//
//        /**
//         * @brief static function.  page aligned, non-overlapping block partitioning of a range
//         * @details Given the number of partition, the partition element desired, and a start and end range, deterministically compute the subrange.
//         *    start end end are allowed to be negative, but end >= start is required.
//         *
//         * @param[in] np        number of partitions
//         * @param[in] pid       id of specific partition desired
//         * @param[in] other     range object to be partitioned
//         * @return              computed subrange
//         */
//        static range<T> paged_block_partition(const size_t &np, const size_t &pid,
//                                              const size_t &page_size,
//                                              const range<T> &other)
//        {
//          return paged_block_partition(np, pid, page_size, other.start, other.end, other.overlap);
//        }
//
//        /**
//         * @brief static function.  page aligned, non-overlapping block partitioning of a range
//         * @details Given the number of partition, the partition element desired, deterministically compute a subrange for the current range object.
//         *    start end end are allowed to be negative, but end >= start is required.
//         *
//         * @param[in] np        number of partitions
//         * @param[in] pid       id of specific partition desired
//         * @return              computed subrange
//         */
//        range<T> paged_block_partition(const size_t &np, const size_t &pid, const size_t &page_size) const
//        {
//          return paged_block_partition(np, pid, page_size, this->start, this->end, this->overlap);
//        }


        /**
         * @brief   align the range to underlying block boundaries, e.g. disk page size
         * @details range is aligned to underlying block boundaries by moving the block_start variable back towards minimum
         *    if range start is too close to the data type's minimum, then assertion is thrown.
         *
         * @param[in] page_size   the size of the underlying block.
         * @return                updated range
         */
        range<T>& align_to_page(const size_t &page_size)
        {
          assert(page_size > 0);

          //printf("page size: %ld\n", static_cast<T>(page_size));

          // change start to align by page size.  extend range start.
          // note that if output.start is negative, it will put block_start at a bigger address than the start.
          block_start = (start / page_size) * page_size;

          if (block_start > start)  // only enters if start is negative.
          {

//            printf("block start: %ld\n", static_cast<size_t>(block_start));

            // if near lowest possible value, then we can't align further..  assert this situation.
            assert(
                (block_start - std::numeric_limits<T>::lowest()) > page_size);

            // deal with negative start position.
            block_start = block_start - page_size;

          }
          // leave end as is.

          return *this;
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


        size_t size() const
        {
          return end - start;
        }

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
