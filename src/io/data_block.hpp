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
     *          can add member function with template parameter to select "True" or "False" for using the buffer.
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
     *
     *
     */




    /**
     * @class     bliss::io::DataBlock
     * @brief     abstraction to represent a block of data.
     * @details
     */
    template<typename Iterator, typename Range, typename Container = void>
    class DataBlock
    {
      public:
        constexpr bool hasContainer = !std::is_same<Container, void>::value;

        typedef typename std""

        constexpr bool containerIsPointer = std::is_pointer<Container>::value;

        // if Container is a pointer type, use it directly as iterator.  else use it directly.
        typedef typename std::conditional<containerIsPointer,
                                          Container,
                                          typename Container::iterator>::type   InternalIterT;
        typedef typename std::iterator_traits<InternalIterT>::value_type        OutputValueType;
        typedef typename std::conditional<containerIsPointer,
                                          std::vector<OutputValueType>,
                                          Container>                            StoreType;
//        typedef Container                                                       StoreType;
        typedef typename StoreType::iterator                                    OutputIteratorType;

        DataBlock() : range(), startIter(), endIter(), buffer()
        {
        }

        /**
         * constructor for non-buffering
         * @param _start
         * @param _end
         * @param _range
         */
        template <bool copy>

        void assign(const Iterator &_start, const Iterator &_end, const Range &_range)
        {
          range = _range;
          startIter = _start;
          endIter = _end;

          buffer.clear();
          std::copy(startIter, endIter, std::inserter(buffer, buffer.begin()));
        }

        virtual ~DataBlock() {
        }

        template<bool buffering>
        typename std::enable_if<buffering, OutputIteratorType>::type
        begin() {
          return buffer.begin();
        }
        OutputIteratorType end() {
          return buffer.end();
        }

        const OutputIteratorType& cbegin() const {
          return buffer.cbegin();
        }
        const OutputIteratorType& cend() const {
          return buffer.cend();
        }

        const Range& getRange() const {
          return range;
        }

      protected:
        Range range;
        Iterator startIter;   // if buffering, kept for future use.
        Iterator endIter;     // kept for future use.
                                        // reference so will update with caller's values
        StoreType buffer;



    };

    /**
     * @class     bliss::io::DataBlock
     * @brief     abstraction to represent a block of data.
     * @details   buffer has to be allocated externally to the right size.  This is to allow reusing a buffer by the calling code.
     *            constructor determines whether buffering or not.  (could have done it as template, but that is not as flexible during runtime.
     *            constructors have _end in case Iterator is not a random access iterator so can't do _start + n
     */
    template<typename Iterator, typename Range>
    class DataBlock<Iterator, Range, void>
    {
      public:
        typedef Iterator OutputIterator;

        DataBlock() : range(), startIter(), endIter()
        {
        }


        /**
         * constructor for non-buffering
         * @param _start
         * @param _end
         * @param _range
         */
        void assign(const Iterator &_start, const Iterator &_end, const Range &_range)
        {
          range = _range;
          startIter = _start;
          endIter = _end;
        }

        virtual ~DataBlock() {
        }

        OutputIterator begin() {
          return startIter;
        }
        OutputIterator end() {
          return endIter;
        }

        const OutputIterator& cbegin() const {
          return startIter;
        }
        const OutputIterator& cend() const {
          return endIter;
        }
        const Range& getRange() const {
          return range;
        }

      protected:
        Range range;
        Iterator startIter;   // if buffering, kept for future use.
        Iterator endIter;     // kept for future use.
                                        // reference so will update with caller's values
    };



  } /* namespace io */
} /* namespace bliss */

#endif /* DATA_BLOCK_HPP_ */
