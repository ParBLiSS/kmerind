/*
 * Copyright 2015 Georgia Institute of Technology
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/**
 * @file    range.hpp
 * @ingroup partition
 * @author  Tony Pan <tpan7@gatech.edu>
 * @brief   Generic representation of an interval on a 1D data structure
 * @details Represents an interval with start, end, and overlap length.
 *   Convenience functions are provided to compute an page aligned starting position that is
 *      suitable for the underlying data storage device, such as disk. however, page aligned start value
 *      is not stored here.
 *
 *
 */

#ifndef RANGE_HPP_
#define RANGE_HPP_

#include <stdexcept>
#include <iostream>       // for printing to ostream
#include <limits>         // for numeric limits
#include <algorithm>      // for min/max

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
     * @details Range specified with offsets.  specific for 1D.
     *          Also works for continuous value ranges (floating points).
     *
     * @note    Overlap concept is dropped since it is not used at the moment
     *
     *            vs std::partition:
     *            	std::partition is defined as std::pair<iterator, iterator>.  this class also supports other types.
     *            	std::partition is heavier weight than Range's start and end values.
     *				Finally, std::partition cannot support floating point ranges.
     *
     *				If T is an iterator, it must be a random access iterators due to need for comparators.
     *
     *
     * @tparam T  data type used for the start and end offsets and overlap.
     */
    template<typename T>
    struct range
    {
    	/// DEFINE value type of the range.
        typedef T ValueType;

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
         * @brief   construct directly from start and end offsets and overlap
         * @details _start should be less than or equal to _finish
         *
         * @param[in] _start    starting position of range in absolute coordinates
         * @param[in] _finish      ending position of range in absoluate coordinates.
         * @param[in] _overlap    amount of overlap between adjacent ranges.  optional
         */
        range(const T &_start, const T &_finish)
            : start(_start), end(_finish)
        {
          if (_finish < _start)
            throw std::invalid_argument("ERROR: range constructor: end is less than start");
        }

        /**
         * @brief   copy construct from the field values of another range
         * @param[in] other   the range object to copy from
         */
        range(const range<T> &other)
            : start(other.start), end(other.end)
        {}

        //============= move constructor and move assignment operator
        //  NOTE: these are NOT defined because they would take more ops than the copy constructor/assignment operator.

        /**
         * @brief   default constructor.  construct an empty range, with start and end initialized to 0.
         */
        range() : start(), end(){}


        /**
         * @brief assignment operator.  copy the field values from the operand range.
         *
         * @param[in] other   the range object to copy from
         * @return            the calling range object, with field values copied from the operand range object
         */
        range<T>& operator =(const range<T> & other)
        {
          start = other.start;
          end = other.end;
          return *this;
        }

        /**
         * @brief static equals function.  compares 2 ranges' start and end positions.
         *
         * @param[in] other   The range to compare to
         * @return            true if 2 ranges are same.  false otherwise.
         */
        static bool equal(const range<T> &lhs, const range<T> &rhs)
        {
          // same if the data range is identical and step is same.

          return (lhs.start == rhs.start) && (lhs.end == rhs.end);
        }

        /**
         * @brief equals function.  compares this range's start and end positions to the "other" range. (excluding overlap region)
         *
         * @param[in] other   The range to compare to
         * @return            true if 2 ranges are same.  false otherwise.
         */
        bool equal(const range<T> &other) const
        {
          // same if the data range is identical and step is same.

          return range<T>::equal(*this, other);
        }

        /**
         * @brief equals operator.  compares 2 ranges' start and end positions only. (excluding overlap region)
         *
         * @param[in] other   The range to compare to
         * @return            true if 2 ranges are same.  false otherwise.
         */
        bool operator ==(const range<T> &other) const
        {
          // same if the data range is identical and step is same.

          return range<T>::equal(*this, other);
        }



        /**
         * @brief Union of "this" range with "other" range, updates "this" range in place.
         * @details   Given 2 ranges R1 and R2, each with start s, and end e,
         *            then union is defined as R = [min(R1.s, R2.s), max(R1.e, R2.e))
         *            choice of end also chooses the overlap.
         * @note      this call requires that the ranges not to be disjoint.
         *            this is a subset of the set-union definition.
         *
         * @param other   the range to form the union with
         */
        void merge(const range<T>& other) {
          if (this->is_disjoint(other) && !this->is_adjacent(other))
             throw std::invalid_argument("ERROR: Range merge() with ranges that are not overlapping or adjacent");


          // merge defined by start and overlap_finish
          start =       (start <= other.start            ) ? start : other.start;
          end   =       (end >= other.end                ) ? end : other.end;
        }

        /**
         * @brief Static function for merging of 2 range, creating a new range object along the way..
         * @details   Given 2 ranges R1 and R2, each with start s, and end e,
         *            then union is defined as R = [min(R1.s, R2.s), max(R1.e, R2.e))
         *            choice of end also chooses the overlap.
         * @note      this call requires that the ranges not to be disjoint.
         *            this is a subset of the set-union definition.
         *
         * @param first   the first range to form the union with
         * @param second  the second range to form the union with
         * @return        a new range object containing the union of this and "other".
         */
        static range<T> merge(const range<T>& first, const range<T>& second)
        {
          range<T> output(first);
          output.merge(second);
          return output;
        }

        /**
         * @brief Intersection of "this" range with "other" range, updates "this" range in place.
         * @details   Given 2 ranges R1 and R2, each with start s, and end e,
         *            then intersection is defined as R = [max(R1.s, R2.s), min(R1.e, R2.e))
         *            choice of end also chooses the overlap.
         * @note      result may be an empty range, which will contain start = end = min(R1.e, R2.e)
         *			  this operation is not symmetric anymore.  if not intersecting,
         *			  	will set range to one of the ends of *this.
         * @param other   the range to form the intersection with
         */
        void intersect(const range<T>& other)
        {
          // intersect when there are no intersections between even overlap region
          auto lstart =       (start >= other.start            ) ? start : other.start;
          auto lend =         (end <= other.end                ) ? end : other.end;

          // check new start and end to see if valid (lstart <= lend).  if so, use them.
          // else set both start and end to end closer to other range.
          auto closer = (lend == end) ? end : start;
          start = (lstart <= lend) ? lstart : closer;  // second half = invalid intersection, then choose start to be closer to the other range.
          end = (lstart <= lend) ? lend : closer;
        }

        /**
         * @brief Static function for Intersection of 2 ranges, return a new Range object..
         * @details   Given 2 ranges R1 and R2, each with start s, and end e,
         *            then intersection is defined as R = [max(R1.s, R2.s), min(R1.e, R2.e))
         *            choice of end also chooses the overlap.
         * @note      result may be an empty range, which will contain start = end = min(R1.e, R2.e)
         *			  operation is not symmetric.  if not intersecting, both ends set to
         *			  	one of the ends of first range.
         * @param first   the first range to form the intersection with
         * @param second   the second range to form the intersection with
         * @return        updated current range containing the intersection of this and "other".
         */
        static range<T> intersect(const range<T>& first, const range<T>& second)
        {
          range<T> output(first);
          output.intersect(second);
          return output;
        }


//        /**
//         * @brief set-theoretic difference between "this" and "other" ranges.  updates "this" range object.
//         * @details following standard definition of set-theoretic difference, also known as "relative complement".
//         *             Given 2 ranges R1 and R2, each with start s, and end e,
//         *             set-theoretic difference of R1 and R2 is denoted as R1 \ R2, defined as part of R1 not in R2.
//         *                then range_difference is then defined as R = [min(R1.s, R2.s), min(R1.e, R2.s))
//         *    cases: other.start < other.end < start < end  : output should be start     <-> end
//         *           other.start < start < other.end < end  : output should be other.end <-> end
//         *           other.start < start < end < other.end  : output should be start     <-> start        or    end <-> end
//         *           start < other.start < other.end < end  : output should be start     <-> other.start  and  other.end <-> end
//         *           start < other.start < end < other.end  : output should be start     <-> other.start
//         *           start < end < other.start < other.end  : output should be start     <-> end
//         *
//         * NOT DEFINED.  DOCUMENTATION LEFT HERE FOR FUTURE REFERENCE.
//         *
//         * @param other
//         * @return
//         */

        /**
         * @brief Shift range to the right (positive direction on number line) by a specified amount.
         * @details   both the start and end position of the range are incremented by the specified amount
         *            the size does not change.
         * @param amount    the amount to right shift the range by
         */
        void shiftRight(const T& amount)
        {
          start += amount;
          end   += amount;
        }
        /**
         * @brief Shift range to the right (positive direction on number line) by a specified amount.
         * @details   both the start and end position of the range are incremented by the specified amount
         *            the size does not change.
         * @param amount    the amount to right shift the range by
         * @return          new range with updated start and end
         */
        static range<T> shiftRight(const range<T>& r, const T& amount)
        {
          range<T> output(r);
          output.shiftRight(amount);
          return output;
        }

        /**
         * @brief Shift range to the right (positive direction on number line) by a specified amount.
         * @details   both the start and end position of the range are incremented by the specified amount
         *            the size does not change.
         * @param amount    the amount to right shift the range by
         * @return          updated range object.
         */
        range<T>& operator +=(const T& amount)
        {
          this->shiftRight(amount);
          return *this;
        }
        /**
         * @brief Shift range to the right (positive direction on number line) by a specified amount.
         * @details   both the start and end position of the range are incremented by the specified amount
         *            the size does not change.
         * @param amount    the amount to right shift the range by
         * @return          new range with updated start and end
         */
        range<T> operator +(const T& amount) const
        {
          range<T> output(*this);
          output.shiftRight(amount);
          return output;
        }


        /**
         * @brief Shift range to the left (negative direction on number line) by a specified amount.
         * @details   both the start and end position of the range are decremented by the specified amount
         *            the size does not change.
         * @param amount    the amount to left shift the range by
         */
        void shiftLeft(const T& amount)
        {
          start -= amount;
          end   -= amount;
        }
        /**
         * @brief Shift range to the left (negative direction on number line) by a specified amount.
         * @details   both the start and end position of the range are incremented by the specified amount
         *            the size does not change.
         * @param amount    the amount to left shift the range by
         * @return          new range with updated start and end
         */
        static range<T> shiftLeft(const range<T>& r, const T& amount)
        {
          range<T> output(r);
          output.shiftLeft(amount);
          return output;
        }

        /**
         * @brief Shift range to the left (negative direction on number line) by a specified amount.
         * @details   both the start and end position of the range are decremented by the specified amount
         *            the size does not change.
         * @param amount    the amount to left shift the range by
         * @return          updated range object.
         */
        range<T>& operator -=(const T& amount)
        {
          this->shiftLeft(amount);
          return *this;
        }
        /**
         * @brief Shift range to the left (negative direction on number line) by a specified amount.
         * @details   both the start and end position of the range are incremented by the specified amount
         *            the size does not change.
         * @param amount    the amount to left shift the range by
         * @return          new range with updated start and end
         */
        range<T> operator -(const T& amount)
        {
          range<T> output(*this);
          output.shiftLeft(amount);
          return output;
        }

        /**
         * @brief Determines if the other range is completely inside this range, and the other range is not a zero-length one.
         * @details       Given 2 ranges R1 and R2, each with start s, and end e,
         *                R1 contains R2 if R1.s <= R2.s, and R1.e >= R2.e
         *
         * @note          Note that this is not communicative.
         *                overlap regions are included in the comparison.
         * @param other   The range object that may be inside this one.
         * @return        bool, true if other is inside this range.
         */
        bool contains(const range<T> &other) const {
          return (other.start >= this->start) && (other.end <= this->end);
        }
        bool contains(T const & v) const {
          return (v >= this->start) && (v < this->end);
        }

        /**
         * @brief Determines if this range is disjoint from the other range.  adjacent == disjoint
         * @details       Given 2 ranges R1 and R2, each with start s, and end e,
         *                R1 is joint to R2 if R1.e < R2.s || R2.e < R1.s.
         *                This includes the "overlap" within each range.
         *                Note that this is communicative.
         *
         *
         *
         * @param other   The range object that may be disjoint from this one.
         * @return        bool, true if other is disjoint from this range.
         */
        bool is_disjoint(const range<T> &other) const {
          return (this->start >= other.end) || (this->end <= other.start);
        }


        /**
         * @brief Determines if this range overlaps the other range.  adjacent == !overlaps.
         * @details       Given 2 ranges R1 and R2, each with start s, and end e,
         *                R1 overlaps R2 if the intersection of R1 and R2 has non-zero size.
         * @note          this is communicative.
         *                overlap region is included.
         * @param other   The range object that may be overlapping this one.
         * @return        bool, true if other overlaps this range.
         */
        bool overlaps(const range<T> &other) const {
          return !is_disjoint(other);
        }


        /**
         * @brief Determines if this range is adjacent to the other range.
         * @details       Given 2 ranges R1 and R2, each with start s, and end e,
         *                R1 is adjacent to R2 if R1.s == R2.e || R2.s == R1.e.
         *                This includes the "overlap" within each range.
         * @note          this is communicative.
         *
         * @param other   The range object that may be adjacent to this one.
         * @return        bool, true if other is adjacent to this range.
         */
        bool is_adjacent(const range<T> &other) const {
          return (this->start == other.end) || (this->end == other.start);
        }


        /**
         * @brief   Static function.  align the range to underlying block boundaries, e.g. disk page size.  only for integral types
         * @details range is aligned to underlying block boundaries by moving the block_start variable back towards minimum
         *    if range start is too close to the data type's minimum, then assertion is thrown.
         * @tparam TT							type of values for start/end.  used to check if type is integral.
         * @param[in] start       the start position of the range to be aligned.
         * @param[in] page_size   the size of the underlying block.
         * @return                the page aligned start position of the range.
         */
        template<typename TT = T>
        static typename std::enable_if<std::is_integral<TT>::value, TT >::type align_to_page(const TT& start, const size_t &page_size)
        {
          if (page_size == 0) throw std::invalid_argument("ERROR: range align_to_page: page size specified as 0.");

          // change start to align by page size.  extfinish range start.
          // note that if output.start is negative, it will put block_start at a bigger address than the start.
          TT block_start = (start / page_size) * page_size;

          if (block_start > start)  // only enters if start is negative.
          {

            // if near lowest possible value, then we can't align further.  assert this situation.
            if (static_cast<size_t>(block_start - std::numeric_limits<TT>::lowest()) < page_size)
              throw std::range_error("ERROR: range align_to_page: start is within a single page size of a signed data type minimum. cannot align page.");

            // deal with negative start position.
            block_start = block_start - page_size;
          }
          // leave end as is.

          return block_start;
        }

        /**
         * @brief   Static function.  align the range to underlying block boundaries, e.g. disk page size.  only for integral types
         * @details range is aligned to underlying block boundaries by moving the block_start variable back towards minimum
         *    if range start is too close to the data type's minimum, then assertion is thrown.
         *
         *   static so Range does not need to be modified during alignment.
         *
         * @tparam TT             type of values for start/end.  used to check if type is integral.
         * @param[in] r           the range to be aligned.
         * @param[in] page_size   the size of the underlying block.
         * @return                the page aligned start position of the range.
         */
        template<typename TT = T>
        static typename std::enable_if<std::is_integral<TT>::value, TT >::type align_to_page(const range<TT>& r, const size_t &page_size)
        {
          return range<T>::align_to_page(r.start, page_size);
        }

        /**
         * @brief     Static function. check to see if the start position has been aligned to underlying block boundary.   only for integral types
         *
				 * @tparam		TT					type of start/end.  used to check if type is integral
         * @param[in] start       the start position of the range to be cehcked for alignment.
         * @param[in] page_size   the size of the underlying block.
         * @return                true if range is block aligned, false otherwise.
         */
        template<typename TT = T>
        static typename std::enable_if<std::is_integral<TT>::value, bool >::type is_page_aligned(const TT& start, const size_t &page_size)
        {
          return (start % page_size) == 0;
        }

        /**
         * @brief     Static function. check to see if the start position has been aligned to underlying block boundary.   only for integral types
         *
         * @tparam    TT          type of start/end.  used to check if type is integral
         * @param[in] r           the range to be checked for alignment.
         * @param[in] page_size   the size of the underlying block.
         * @return                true if range is block aligned, false otherwise.
         */
        template<typename TT = T>
        static typename std::enable_if<std::is_integral<TT>::value, bool >::type is_page_aligned(const range<TT>& r, const size_t &page_size)
        {
          return range<T>::is_page_aligned(r.start, page_size);
        }

        /**
         * @brief   get the integral size of the range between [start, end), including the overlap region
         * @return  unsigned length of the range.
         */
        template<typename TT=T>
        typename std::enable_if<!std::is_floating_point<TT>::value, size_t>::type size() const
        {
          return end - start;   // should allow both integral and iterator TT.
        }

        /**
         * @brief   get the floating point size of the range between [start, end), including the overlap region
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
      ost << "range: block [" << static_cast<int64_t>(r.start) << ":" << static_cast<int64_t>(r.end) << ")";
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
      ost << "range: block [" << static_cast<uint64_t>(r.start) << ":" << static_cast<uint64_t>(r.end) << ")";
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
      ost << "range: block [" << static_cast<double>(r.start) << ":" << static_cast<double>(r.end) << ")";
      return ost;
    }


  } /* namespace partition */
} /* namespace bliss */
#endif /* RANGE_HPP_ */
