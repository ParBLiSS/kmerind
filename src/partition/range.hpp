/**
 * @file    range.hpp
 * @ingroup bliss::partition
 * @author  Tony Pan
 * @brief   Generic representation of an interval on a 1D data structure
 * @details Represents an interval with start, end, and overlap length.
 *   Also contains a block_start to mark beginning of a underlying data block, such as paging size.
 *
 * @copyright BLISS: Copyright (c) Georgia Institute of Technology.  All Rights Preserved.
 *
 * TODO add Licence
 *
 */

#ifndef RANGE_HPP_
#define RANGE_HPP_

#include <stdexcept>
#include <iostream>       // for printing to ostream
#include <limits>         // for numeric limits
#include <algorithm>      // for min/max

//#include <partition/partitioner.hpp>

namespace bliss
{
  /**
   * @namespace partition
   */
  namespace partition
  {

    /**
     * @class range
     * @brief   Generic representation of an interval on a 1D data structure
     * @details Range specified with offsets and overlap.  specific for 1D.  overlap is on END side only, and is included in the END.
     *          Works for continuous value ranges (floating points).
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
          if (_end < _start)
            throw std::invalid_argument("ERROR: range constructor: end is less than start");
          if (_overlap < 0) throw std::invalid_argument("ERROR: range constructor: overlap is less than 0");
        }

        /**
         * @brief   copy construct from the field values of another range
         * @param[in] other   the range object to copy from
         */
        range(const range<T> &other)
            : block_start(other.block_start), start(other.start),
              end(other.end), overlap(other.overlap)
        {}

        ////// move constructor and move assignment operator
        ///  NOTE: these are NOT defined because they would take more ops than the copy constructor/assignment operator.


        /**
         * @brief   default constructor.  construct an empty range, with start and end initialized to 0.
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
        static bool equal(const range<T> &self, const range<T> &other) const
        {
          // same if the data range is identical and step is same.
          // not comparing overlap or block start.

          return (self.start == other.start) && (self.end == other.end);
        }

        /**
         * @brief equals operator.  compares 2 ranges' start and end positions only.
         *
         * @param[in] other   The range to compare to
         * @return  true if 2 ranges are same.  false otherwise.
         */
        bool equal(const range<T> &other) const
        {
          // same if the data range is identical and step is same.
          // not comparing overlap or block start.

          return (start == other.start) && (end == other.end);
        }


        /**
         * @brief union of range (as |=)  NOTE: result could include previously unincluded ranges
         * @param other
         * @return
         */
        range<T>& range_union(const range<T>& other) {
          block_start = std::min(block_start, other.block_start);
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
        static range<T> range_union(const range<T>& self, const range<T>& other) const
        {
          range<T> output(*this);
          output.range_union(other);
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
         * @brief intersection of range (as &)
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
         * @brief determines if this range contains the other range.
         * @param other   The range object that may be inside this one.
         * @return    bool, true if other is inside this range.
         */
        bool contains(const range<T> &other) const {
          return (other.start >= this->start) && (other.end <= this->end);
        }

        /**
         * @brief determines if this range overlaps the other range.
         * @param other   The range object that may be overlapping this one.
         * @return    bool, true if other overlaps this range.
         */
        bool overlaps(const range<T> &other) const {
          return (other & *this).size() > 0;
        }


        /**
         * @brief   align the range to underlying block boundaries, e.g. disk page size.  only for integral types
         * @details range is aligned to underlying block boundaries by moving the block_start variable back towards minimum
         *    if range start is too close to the data type's minimum, then assertion is thrown.
         * @tparam TT							type of values for start/end.  used to check if type is integral.
         * @param[in] page_size   the size of the underlying block.
         * @return                updated range
         */
        template<typename TT = T>
        typename std::enable_if<std::is_integral<TT>::value, range<TT> >::type& align_to_page(const size_t &page_size)
        {
          if (page_size == 0) throw std::invalid_argument("ERROR: range align_to_page: page size specified as 0.");

          // change start to align by page size.  extend range start.
          // note that if output.start is negative, it will put block_start at a bigger address than the start.
          block_start = (start / page_size) * page_size;

          if (block_start > start)  // only enters if start is negative.
          {

            // if near lowest possible value, then we can't align further.  assert this situation.
            if ((block_start - std::numeric_limits<TT>::lowest()) < page_size)
              throw std::range_error("ERROR: range align_to_page: start is within a single page size of a signed data type minimum. cannot align page.");

            // deal with negative start position.
            block_start = block_start - page_size;
          }
          // leave end as is.

          return *this;
        }

        /**
         * @brief     check to see if the range has been aligned to underlying block boundary.   only for integral types
         *
				 * @tparam		TT					type of start/end.  used to check if type is integral
         * @param[in] page_size   the size of the underlying block.
         * @return    true if range is block aligned, false otherwise.
         */
        template<typename TT = T>
        typename std::enable_if<std::is_integral<TT>::value, bool >::type is_page_aligned(const size_t &page_size) const
        {
          return (this->block_start % page_size) == 0;
        }


        /**
         * @brief   get the integral size of the range between [start, end)
         * @return  unsigned length of the range.
         */
        template<typename TT=T>
        typename std::enable_if<std::is_integral<TT>::value, size_t>::type size() const
        {
          return end - start;
        }

        /**
         * @brief   get the floating point size of the range between [start, end)
         * @return  length of the range.
         */
        template<typename TT=T>
        typename std::enable_if<std::is_floating_point<TT>::value, TT>::type size() const
        {
          return end - start;
        }

    };

    /**
     * @brief << operator to write out range object's fields.  Signed integral version
     * @tparam  T           type of values used by Range internally.  This is type deduced by compiler.
     * @param[in/out] ost   output stream to which the content is directed.
     * @param[in]     r     range object to write out
     * @return              output stream object
     */
    template<typename T>
    typename std::enable_if<std::is_signed<T>::value and std::is_integral<T>::value, std::ostream>::type& operator<<(std::ostream& ost, const range<T>& r)
    {
      ost << "range: block@" << static_cast<int64_t>(r.block_start) << " [" << static_cast<int64_t>(r.start) << ":" << static_cast<int64_t>(r.end) << ") overlap " << static_cast<int64_t>(r.overlap);
      return ost;
    }

    /**
     * @brief << operator to write out range object's fields.  Unsigned integral version
     * @tparam  T           type of values used by Range internally.  This is type deduced by compiler.
     * @param[in/out] ost   output stream to which the content is directed.
     * @param[in]     r     range object to write out
     * @return              output stream object
     */
    template<typename T>
    typename std::enable_if<!std::is_signed<T>::value and std::is_integral<T>::value, std::ostream>::type& operator<<(std::ostream& ost, const range<T>& r)
    {
      ost << "range: block@" << static_cast<uint64_t>(r.block_start) << " [" << static_cast<uint64_t>(r.start) << ":" << static_cast<uint64_t>(r.end) << ") overlap " << static_cast<uint64_t>(r.overlap);
      return ost;
    }

    /**
     * @brief << operator to write out range object's fields.  floating point version
     * @tparam  T           type of values used by Range internally.  This is type deduced by compiler.
     * @param[in/out] ost   output stream to which the content is directed.
     * @param[in]     r     range object to write out
     * @return              output stream object
     */
    template<typename T>
    typename std::enable_if<std::is_floating_point<T>::value, std::ostream>::type& operator<<(std::ostream& ost, const range<T>& r)
    {
      ost << "range: block@" << static_cast<double>(r.block_start) << " [" << static_cast<double>(r.start) << ":" << static_cast<double>(r.end) << ") overlap " << static_cast<double>(r.overlap);
      return ost;
    }


  } /* namespace partition */
} /* namespace bliss */
#endif /* RANGE_HPP_ */
