/**
 * @file		data_block.hpp
 * @ingroup
 * @author	tpan
 * @brief
 * @details
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef DATA_BLOCK_HPP_
#define DATA_BLOCK_HPP_

#include <cassert>

#include <iterator>
#include <type_traits>
#include <vector>
#include <algorithm>

#include "iterators/container_traits.hpp"
#include "iterators/range.hpp"

namespace bliss
{
  namespace io
  {

    /**
     *  issue with reference:  need to initialize this object with references, so harder to do for member variable instances.

     *  design constraints:
     *   1. support iterators  (pointer, normal iterators) for input, iterators for output
     *   2. managers the buffer memory internally.
     *        a. not thread safe.  okay.
     *        b. always available, not always used.
     *   3. reuse the buffer memory.
     *        a. have "assign" and "clear" functions.
     *   4. choice of use of buffer is at runtime
     *        a. assign is given the directive to buffer or not.  clear always clears.
     *        b. need to keep state so know whether buffering has occurred.
     *   4. hide the fact that data has (not) been buffered in memory
     *        a. interface has "begin", "end" functions, pointing to either original data, or buffered data.
     *        b. this means that there are 2 different type of iterators that could be returned.
     *          how to present a unified "begin" and "end" interface given the input and output iterators are different?
     *          can add member function with template parameter to select "True" or "False" for using the buffer.  Providing 2 types BUFFERING and NO_BUFFERING for this.
     *
     *  use "assign" to return one of 2 (buffering/non-buffering) internal objects,  have caller use "auto" for the object type, and "auto" for begin/end output type.
     *
     *  The DataBlock class is used for both having a backing store and not, decidable at run time.
     *
     *  FORGET ABOUT HAVING 1 CLASS THAT ALLOWS BOTH RUNTIME BUFFERING AND NOT.  complicating the usage and API.
     *    DO TEMPLATING.
     *    Type is templated.  Interface is identical.
     *
     */


    template<typename Derived, typename Iterator, typename Range, typename OutputIterator>
    class DataBlock {
      public:
        typedef typename std::iterator_traits<Iterator>::value_type           ValueType;
        static_assert(!(std::is_same<ValueType, void>::value), "Iterator is NOT valid.");
        typedef OutputIterator                                    iterator;
        typedef typename std::add_const<OutputIterator>::type     const_iterator;

        const Range& getRange() const {
          return range;
        }
        bool isEmpty() {
          return empty;
        }

        bool hasBuffer() {
          return static_cast<Derived*>(this)->hasBufferImpl();
        }
        void assign(const Iterator &_start, const Iterator &_end, const Range &_range) {
          empty = _range.end <= _range.start;
          range = _range;
          static_cast<Derived*>(this)->clearImpl();
          if (!empty) static_cast<Derived*>(this)->assignImpl(_start, _end);
        }
        void clear() {
          empty = true;
          range = Range();
          static_cast<Derived*>(this)->clearImpl();
        }

      protected:
        Range range;
        bool empty;

    };




    template<typename Iterator, typename Range, typename Container = std::vector<typename std::iterator_traits<Iterator>::value_type> >
    class BufferedDataBlock : public DataBlock<BufferedDataBlock<Iterator, Range, Container>, Iterator, Range, typename Container::iterator>
    {
      public:
        typedef DataBlock<BufferedDataBlock<Iterator, Range, Container>, Iterator, Range, typename Container::iterator> SuperType;
        static_assert(is_container<Container>::value, "Container is not valid.  Should support begin() and end() at the least.");
        static_assert(std::is_same<typename SuperType::ValueType, typename Container::value_type>::value, "Iterator and Container should have the same element types");

        bool hasBufferImpl() {
          return true;
        }
        void assignImpl(const Iterator &_start, const Iterator &_end) {
          std::copy(_start, _end, std::inserter(buffer, buffer.begin()));
        }
        void clearImpl() {
          buffer.clear();
        }

        /**
         *
         * @return
         */
        typename Container::iterator begin()
        {
          return buffer.begin();
        }

        /**
         *
         * @return
         */
        typename Container::const_iterator cbegin() const
        {
          return buffer.cbegin();
        }


        /**
         *
         * @return
         */
        typename Container::iterator end()
        {
          return buffer.end();
        }

        /**
         *
         * @return
         */
        typename Container::const_iterator cend() const
        {
          return buffer.cend();
        }


      protected:
        Container buffer;
    };

    /**
     * @brief << operator to write out range object's fields.
     * @param[in/out] ost   output stream to which the content is directed.
     * @param[in]     r     range object to write out
     * @return              output stream object
     */
    template<typename Iterator, typename Range, typename Container>
    std::ostream& operator<<(std::ostream& ost, const BufferedDataBlock<Iterator, Range, Container>& db)
    {
      std::ostream_iterator<typename Container::value_type> oit(ost);
      std::copy(db.cbegin(), db.cend(), oit);
      return ost;
    }



    template<typename Iterator, typename Range>
    class UnbufferedDataBlock :
        public DataBlock<UnbufferedDataBlock<Iterator, Range>, Iterator, Range, typename std::remove_const<typename std::remove_reference<Iterator>::type>::type>
    {
      public:
        typedef DataBlock<UnbufferedDataBlock<Iterator, Range>, Iterator, Range, typename std::remove_const<typename std::remove_reference<Iterator>::type>::type>  SuperType;

        bool hasBufferImpl() {
          return false;
        }
        void assignImpl(const Iterator &_start, const Iterator &_end) {
          startIter = _start;
          endIter = _end;
        }
        void clearImpl() {
          startIter = endIter = Iterator();
        }

        typename std::remove_const<typename std::remove_reference<Iterator>::type>::type begin()
        {
          return startIter;
        }
        typename std::add_lvalue_reference<typename std::add_const<Iterator>::type>::type cbegin() const
        {
          return startIter;
        }
        typename std::remove_const<typename std::remove_reference<Iterator>::type>::type end()
        {
          return endIter;
        }
        typename std::add_lvalue_reference<typename std::add_const<Iterator>::type>::type cend() const
        {
          return endIter;
        }

      protected:
        Iterator startIter;   // if buffering, kept for future use.
        Iterator endIter;     // kept for future use.
    };

    /**
     * @brief << operator to write out range object's fields.
     * @param[in/out] ost   output stream to which the content is directed.
     * @param[in]     r     range object to write out
     * @return              output stream object
     */
    template<typename Iterator, typename Range>
    std::ostream& operator<<(std::ostream& ost, const UnbufferedDataBlock<Iterator, Range>& db)
    {
      std::ostream_iterator<typename std::iterator_traits<Iterator>::value_type > oit(ost);
      std::copy(db.cbegin(), db.cend(), oit);
      return ost;
    }



  } /* namespace io */
} /* namespace bliss */

#endif /* DATA_BLOCK_HPP_ */
