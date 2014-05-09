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
#include <vector>
#include <algorithm>

#include "iterators/container_traits.hpp"

namespace bliss
{
  namespace io
  {

    /**
     *  design constraints:
     *   1. hide the fact that data has (not) been replicated in memory (buffered)
     *        a. interface has "begin", "end" functions, pointing to either original data, or buffered data.
     *   2. reuse caller defined buffer memory when supplied
     *        a. need to handle data not owned by self - thread life time
     *        b. (perhaps we should just make it clear that dataBlock IS the owner of the buffer memory, and this is NOT threadsafe?)
     *        c. iterator does NOT allow clearing the previous data?
     *            see 2b
     *   3. avoid excessive instantiations.
     *        a. don't want internal "container", even if empty, if possible.
     *            See 2b.  can have internal container, as long as lifetime of data_block is not too short.
     *        b. pointer to container, or reference to container.
     *            see 2b
     *   4. support iterators for input, containers/iterators for output
     *        output iterator should be limited to pointers - container's iterators are for reading or modifying existing stuff, not for inserting.
     *          can use back_insert_iterator adaptor, but would need to clear the container first.
     *            see 2b.
     *   5. choose the buffering or none buffering variant at runtime.
     *        member variable of calling class should be of a base type that is NOT templated based on buffering requirement.
     *          this is hard to achieve because runtime information is needed.
     *          can template the caller, but this gives less flexibiltiy.
     *          can make this class choose at runtime whether to buffer or not.  may be expensive.
     *          can add member function with template parameter to select "True" or "False" for using the buffer.  Providing 2 types BUFFERING and NO_BUFFERING for this.
     *
     *  issue with reference:  need to initialize this object with references, so harder to do for member variable instances.
     *
     *  have ctor in 2 forms:  with buffer size, and without buffer size.
     *     with buffer size -> construct a local buffer.
     *     without buffer size -> no local buffer.
     *  have function for "assign" - without buffer, store the start and end of iterator.  with buffer, copy of content into local data buffer.
     *    this means that the member variables can NOT be const refs.
     *  rest is same: begin, end, getRange (returns global coord for the data block.
     *
     *  container can be specified as T*.  internally it will use a std::vector in that case.
     *
     *  The DataBlock class is used for both having a backing store and not, decidable at run time.
     *
     *
     *  member function specialization:
     *
     */
    struct BUFFER_CHOICE {};
    struct BUFFER_ON : public BUFFER_CHOICE {};
    struct BUFFER_OFF : public BUFFER_CHOICE {};

    /**
     * @class     bliss::io::DataBlock
     * @brief     abstraction to represent a block of data.
     * @details   container type defaults to Iterator type.
     *            Iterator:  can be pointer or an iterator class.
     *            Container:  defaults to Iterator.  container is pointer or a regular iterator,
     *              then default to a std::vector.
     */
    template<typename Iterator, typename Range,
             typename Container = std::vector<typename std::iterator_traits<Iterator>::value_type> >
    class DataBlock
    {
      public:

        typedef typename std::iterator_traits<Iterator>::value_type           ValueType;

      protected:
        static constexpr bool iteratorIsValid = !(std::is_same<ValueType, void>::value);

        static_assert(iteratorIsValid, "Iterator is NOT valid.");
        static_assert(is_container<Container>::value, "Container is not valid.  Should support begin() and end() at the least.");
        static_assert(std::is_same<ValueType, typename Container::value_type>::value, "Iterator and Container should have the same element types");



      public:

        /**
         *
         */

        // begin and end functions need to choose the return iterator type depending on information available at time of "calling".
        DataBlock() : range(), startIter(), endIter(), buffer(), buffered(false)
        {
        }

        /**
         * constructor for non-buffering
         * @param _start
         * @param _end
         * @param _range
         */
        template<typename buffering>
        void assign(const Iterator &_start, const Iterator &_end, const Range &_range, buffering) {
          range = _range;
          startIter = _start;
          endIter = _end;

          buffer.clear();
          buffered = std::is_same<buffering, BUFFER_ON>::value;
          if (buffered)
            std::copy(startIter, endIter, std::inserter(buffer, buffer.begin()));
        }

        /**
         *
         */
        virtual ~DataBlock() {
          buffer.clear();
        }

//        /**
//         *
//         * @return
//         */
//        template<typename buffering>
//        typename std::enable_if<std::is_same<buffering, BUFFER_ON>::value, typename Container::iterator>::type begin(buffering = BUFFER_ON())
//        {
//          assert(buffered);
//          return buffer.begin();
//        }
//        template<typename buffering>
//        typename std::enable_if<std::is_same<buffering, BUFFER_OFF>::value, Iterator>::type begin(buffering = BUFFER_OFF())
//        {
//          return startIter;
//        }

        // TODO:  use auto return type.
        auto begin(BUFFER_ON) {  // was typename Container::iterator
          assert(buffered);
          return buffer.begin();
        }
        auto begin(BUFFER_OFF) {  // was Iterator
          return startIter;
        }
        /**
         *
         * @return
         */
        template<typename buffering>
        typename std::enable_if<std::is_same<buffering, BUFFER_ON>::value, typename std::add_lvalue_reference<typename Container::const_iterator>::type >::type cbegin(buffering = BUFFER_ON()) const
        {
          assert(buffered);
          return buffer.cbegin();
        }

        template<typename buffering>
        typename std::enable_if<std::is_same<buffering, BUFFER_OFF>::value, typename std::add_lvalue_reference<typename std::add_const<Iterator>::type>::type >::type cbegin(buffering = BUFFER_OFF()) const
        {
          return startIter;
        }

        /**
         *
         * @return
         */
        template<typename buffering>
        typename std::enable_if<std::is_same<buffering, BUFFER_ON>::value, typename Container::iterator>::type end(buffering = BUFFER_ON())
        {
          assert(buffered);
          return buffer.end();
        }
        template<typename buffering>
        typename std::enable_if<std::is_same<buffering, BUFFER_OFF>::value, Iterator>::type end(buffering = BUFFER_OFF())
        {
          return endIter;
        }
        /**
         *
         * @return
         */
        template<typename buffering>
        typename std::enable_if<std::is_same<buffering, BUFFER_ON>::value, typename std::add_lvalue_reference<typename Container::const_iterator>::type >::type cend(buffering = BUFFER_ON()) const
        {
          assert(buffered);
          return buffer.cend();
        }

        template<typename buffering>
        typename std::enable_if<std::is_same<buffering, BUFFER_OFF>::value, typename std::add_lvalue_reference<typename std::add_const<Iterator>::type>::type >::type cend(buffering = BUFFER_OFF()) const
        {
          return endIter;
        }
        /**
         *
         * @return
         */
        const Range& getRange() const {
          return range;
        }

      protected:
        Range range;
        Iterator startIter;   // if buffering, kept for future use.
        Iterator endIter;     // kept for future use.

        Container buffer;
        bool buffered;
    };


  } /* namespace io */
} /* namespace bliss */

#endif /* DATA_BLOCK_HPP_ */
