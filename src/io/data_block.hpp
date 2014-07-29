/**
 * @file		data_block.hpp
 * @ingroup bliss::io
 * @author	tpan
 * @brief   contains classes for managing an in memory 1D address range, with or without buffering
 * @details presents a standardized interface for buffered and unbuffered data blocks.
 *
 *          DataBlock is the base class.  from which there are 2 subclasses:  BufferedDataBlock and UnbufferedDataBlock
 *          both subclasses share the same API set with different return types:
 *            begin/end
 *            cbegin/cend
 *            assign/clear
 *            hasBuffer
 *            isEmpty (common)
 *            getRange (common)
 *
 *          5 design approaches were considered:
 *          1. a single DataBlock class with runtime choice of buffering
 *              problem:  would require iterator accessors to take a parameter to enable overloading (for the different return types)
 *                  logic would be complicated
 *              benefit:  runtime choice of overloading.
 *          2. a single DataBlock class with template specialization to choose the iterator accessor return types.
 *              problem:  use enable_if - needs a boolean value that is based on a user choice
 *                  either a boolean template paramter will be needed, or the container type of void means no buffering
 *                  a large percentage of functions (first 4 sets above) will need to be specialized with enable_if, thus replicating declarations
 *          3. template specialization with an empty default, and 2 specialized classes
 *              problem:  same requirement that there is some way to distinguish between buffering and not buffering
 *                  replicated functions in different specializations
 *              benefit:  fewer enable_if.
 *                  easily specialize for new container types (e.g. raw arrays)
 *          4. template specialization with default being the UnbufferedDataBlock type
 *              problem and benefits :  same as 3
 *              problem:  new parameter types (e.g. container type) may be incorrectly handled by default.
 *          5. CRTP with DataBlock being the base and Buffered and Unbuffered as derived classes
 *              problem: replicated APIs and some code.  more complex hierarchy.  harder to read and understand.
 *              benefit: common member and methods are in the base class.
 *
 *
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
#include <utility>

#include "utils/container_traits.hpp"
#include "partition/range.hpp"

namespace bliss
{
  namespace io
  {

    /**
     * @class   DataBlock
     * @brief   templated, high level abstraction of a subrange of a sequence of data elements.
     * @details Given a range and the corresponding raw data, the goal of this class is to
     *          present a consistent iterator-based interface for traversing the data,
     *          with or without buffering.  The challenge is when buffering, data is copied into a
     *          container, which will have a different data type than the source data's iterator type.
     *
     *          An example usage is to represent portions of a memory mapped file.  When buffered, the iterator
     *          type is no longer pointer, but rather std::vector::iterator.
     *
     *          DataBlock's design requirements are listed below:
     *           1. stores iterators (pointer, normal iterators) for input and output traversal of the data
     *           2. manages the buffer memory internally and transparently.
     *                a. Okay not to be thread safe.
     *                a. iterators for either original data, or buffered data accessed through begin()/end(),etc.
     *           3. reuse the buffer memory when a BufferedBlock is remapped to another range.
     *                a. have "assign" and "clear" functions.
     *           4. store the conceptual "range" of the data, e.g. offsets of the DataBlock within the entire available data sequence.
     *
     *
     * @tparam  Derived         Derived class's type, used by CRTP to call the correct impl method for the public interfaces in thsi class
     * @tparam  Iterator        Iterator type for the source data
     * @tparam  Range           Range within the overall data element sequence that this block is wrapping.
     * @tparam  OutputIterator  Iterator type for the output
     */
    template<typename Derived, typename Iterator, typename Range, typename OutputIterator>
    class DataBlock {
      public:
        /**
         * @typedef   ValueType
         * @brief     type of the data elements in the DataBlock
         */
        typedef typename std::iterator_traits<Iterator>::value_type           ValueType;
        // assert that ValueType is defined.
        static_assert(!(std::is_same<ValueType, void>::value), "Iterator is NOT valid.");

        /**
         * @typedef iterator
         * @brief   type of output iterator
         */
        typedef OutputIterator                                    iterator;

        /**
         * @typedef const_iterator
         * @brief   const version of output iterator
         */
        typedef typename std::add_const<OutputIterator>::type     const_iterator;

        /**
         *  @brief  default constructor
         */
        DataBlock() : range(), empty(true) {};

        /**
         * @brief         move constructor.  moves range and empty flag.
         * @param other   source DataBlock to move
         */
        DataBlock(DataBlock<Derived, Iterator, Range, OutputIterator>&& other) : range(other.range), empty(other.empty) {
          other.range = Range();
          other.empty = true;
        }

        /**
         * @brief         move assignment operator.  moves range and empty flag.
         * @param other   source DataBlock to move
         * @return        self, updated
         */
        DataBlock<Derived, Iterator, Range, OutputIterator>& operator=(DataBlock<Derived, Iterator, Range, OutputIterator>&& other) {
          if (this != &other) {
            empty = other.empty;  other.empty = true;
            range = other.range;  other.range = Range();
          }
          return *this;
        }

        /**
         * @brief         copy constructor.  copies range and empty flag.
         * @param other   source DataBlock to copy
         */
        DataBlock(const DataBlock<Derived, Iterator, Range, OutputIterator>& other) : range(other.range), empty(other.empty){
        }

        /**
         * @brief         copy assignment operator.  copies range and empty flag.
         * @param other   source DataBlock to copy
         * @return        self, updated
         */
        DataBlock<Derived, Iterator, Range, OutputIterator>& operator=(const DataBlock<Derived, Iterator, Range, OutputIterator>& other) {
          if (this != &other) {
            empty = other.empty;
            range = other.range;
          }
          return *this;
        }

        /**
         * @brief default destructor
         */
        virtual ~DataBlock() {};

        /**
         * @brief   accessor to get the range (offset start/end) represented by this DataBlock
         * @return  Range of the DataBlock
         */
        const Range& getRange() const {
          return range;
        }

        /**
         * @brief   indicate if the DataBlock is empty (e.g. newly allocated)
         * @return  bool, true if DataBlock is empty.
         */
        bool isEmpty() {
          return empty;
        }

        /**
         * @brief   indicate if the DataBlock is buffering. calls the subclass' impl method.
         * @return  bool, true if DataBlock is buffering.
         */
        bool hasBuffer() {
          return static_cast<Derived*>(this)->hasBufferImpl();
        }

        /**
         * @brief    assigns some data to the DataBlock.  check if empty.  if not, then call derived class to assign.
         * @param _start    begin iterator  (input)
         * @param _end      end iterator    (input)
         * @param _range    the corresponding range (offset values)
         */
        void assign(const Iterator &_start, const Iterator &_end, const Range &_range) {
          empty = _range.end <= _range.start;
          range = _range;
          static_cast<Derived*>(this)->clearImpl();
          if (!empty) static_cast<Derived*>(this)->assignImpl(_start, _end);
        }

        /**
         * @brief   clears the DataBlock by setting it to empty, clear the Range, and call the derived class's clear implementation.
         */
        void clear() {
          empty = true;
          range = Range();
          static_cast<Derived*>(this)->clearImpl();
        }

      protected:
        /**
         * @var range
         * @brief   range object describing the start and end of the data block, in 1D coordinates (i.e. offsets)
         */
        Range range;

        /**
         * @var empty
         * @brief   boolean flag indicating that the datablock does not contain valid data.
         */
        bool empty;

    };



    /**
     * @class BufferedDataBlock
     * @brief   BufferedDataBlock uses an internally allocated container to mirror the source data.
     * @details Source data type may be pointer or a proper iterator type.  The container type then determine the output iterator type
     *          The class derives from the base DataBlock data type using CRTP.
     * @tparam Iterator   Source data iterator type.
     * @tparam Range      Range data type
     * @tparam Container  container type for buffer.  defaults to std::vector.
     */
    template<typename Iterator, typename Range, typename Container = std::vector<typename std::iterator_traits<Iterator>::value_type> >
    class BufferedDataBlock : public DataBlock<BufferedDataBlock<Iterator, Range, Container>, Iterator, Range, typename Container::iterator>
    {
      protected:
        /**
         * @typedef SuperType
         * @brief superclass's type
         */
        typedef DataBlock<BufferedDataBlock<Iterator, Range, Container>, Iterator, Range, typename Container::iterator> SuperType;
        static_assert(std::is_same<typename SuperType::ValueType, typename Container::value_type>::value, "Iterator and Container should have the same element types");

        // check to make sure data types are correct.
        static_assert(bliss::utils::container_traits<Container>::isAssignable, "ERROR: Container needs to be a sequence container with 'assign' method.");
        static_assert(bliss::utils::container_traits<Container>::isConstIterable, "ERROR: Container needs to be a sequence container with 'cbegin' and 'cend' methods.");
        static_assert(bliss::utils::container_traits<Container>::isIterable, "ERROR: Container needs to be a sequence container with 'begin' and 'end' methods.");

      public:

        /**
         * @brief default constructor
         */
        BufferedDataBlock() : SuperType() {};

        /**
         * @brief         move constructor.  calls baseclass' move constructor, and also moves the buffered data
         * @param other   source BufferedDataBlock to move
         */
        BufferedDataBlock(BufferedDataBlock<Iterator, Range, Container>&& other) :
          SuperType(std::forward(other)), buffer(std::move(other.buffer)) {
        }

        /**
         * @brief         move assignment operator.  calls baseclass' move constructor, and also moves the buffered data
         * @param other   source BufferedDataBlock to move
         * @return        self, updated
         */
        BufferedDataBlock<Iterator, Range, Container>& operator=(BufferedDataBlock<Iterator, Range, Container>&& other) {
          if (this != &other) {
            this->SuperType::operator=(std::forward(other));
            buffer = std::move(other.buffer);
          }
          return *this;
        }

        /**
         * @brief         copy constructor.  calls baseclass' copy constructor, and also copys the buffered data
         * @param other   source BufferedDataBlock to copy
         */
        BufferedDataBlock(const BufferedDataBlock<Iterator, Range, Container>& other) : SuperType(other), buffer(other.buffer) {
        }

        /**
         * @brief         copy assignment operator.  calls baseclass' copy constructor, and also copys the buffered data
         * @param other   source BufferedDataBlock to copy
         * @return        self, updated
         */
        BufferedDataBlock<Iterator, Range, Container>& operator=(const BufferedDataBlock<Iterator, Range, Container>& other) {
          if (this != &other) {
            this->SuperType::operator=(other);
            buffer = other.buffer;
          }
          return *this;
        }

        /**
         * @brief default destructor
         */
        virtual ~BufferedDataBlock() {};


        /**
         * @brief   indicates if this class is buffering.  BufferDataBlock specific implementation.
         * @return  bool.  always true for BufferedDataBlock
         */
        bool hasBufferImpl() {
          return true;
        }

        /**
         * @brief  copies data between start and end into container.
         * @param _start    iterator pointing to beginning of data to copy
         * @param _end      iterator pointing to end of data to copy.
         */
        void assignImpl(const Iterator &_start, const Iterator &_end) {
          // using assign (and let buffer grow, do not explicitly clear in between) sped things up TREMENDOUSLY.
          buffer.assign(_start, _end);
        }

        /**
         * @brief   clears the buffer.  BufferDataBlock specific implementation - calls "clear" on the container
         */
        void clearImpl() {
          buffer.clear();
        }

        /**
         * @brief   return the start of the buffered datablock.  points to the start of the internal container.
         * @return  iterator pointing to beginning of buffer.
         */
        typename Container::iterator begin()
        {
          return buffer.begin();
        }

        /**
         * @brief   return the start of the buffered datablock.  points to the start of the internal container.
         * @return  const iterator pointing to beginning of buffer.
         */
        typename Container::const_iterator cbegin() const
        {
          return buffer.cbegin();
        }


        /**
         * @brief   return the end of the buffered datablock.  points to the end of the internal container.
         * @return  iterator pointing to end of buffer.
         */
        typename Container::iterator end()
        {
          return buffer.end();
        }

        /**
         * @brief   return the end of the buffered datablock.  points to the end of the internal container.
         * @return  const iterator pointing to end of buffer.
         */
        typename Container::const_iterator cend() const
        {
          return buffer.cend();
        }


      protected:
        /**
         * @var buffer
         * @brief container used as buffer to hold a copy of the input data
         */
        Container buffer;
    };

    /**
     * @brief << operator to write out BufferedDataBlock object's fields.
     * @tparam Iterator   Source data iterator type.
     * @tparam Range      Range data type
     * @tparam Container  container type for buffer.  defaults to std::vector.
     * @param[in/out] ost   output stream to which the content is directed.
     * @param[in]     db    BufferedDataBlock object to write out
     * @return              output stream object
     */
    template<typename Iterator, typename Range, typename Container>
    std::ostream& operator<<(std::ostream& ost, const BufferedDataBlock<Iterator, Range, Container>& db)
    {
      std::ostream_iterator<typename Container::value_type> oit(ost);
      std::copy(db.cbegin(), db.cend(), oit);
      return ost;
    }


    /**
     * @class UnbufferedDataBlock
     * @brief   UnbufferedDataBlock uses .
     * @details Source data type may be pointer or a proper iterator type.  The container type then determine the output iterator type
     *          The class derives from the base DataBlock data type using CRTP.
     * @tparam Iterator   Source data iterator type.
     * @tparam Range      Range data type
     */
    template<typename Iterator, typename Range>
    class UnbufferedDataBlock :
        public DataBlock<UnbufferedDataBlock<Iterator, Range>, Iterator, Range, typename std::remove_const<typename std::remove_reference<Iterator>::type>::type>
    {
      public:
        typedef DataBlock<UnbufferedDataBlock<Iterator, Range>, Iterator, Range, typename std::remove_const<typename std::remove_reference<Iterator>::type>::type>  SuperType;

        /// constructor is default
        UnbufferedDataBlock() : SuperType() {};

        /// move constructor and assignment operator
        UnbufferedDataBlock(UnbufferedDataBlock<Iterator, Range>&& other) :
          SuperType(std::forward(other)), startIter(other.startIter), endIter(other.endIter) {
          other.startIter = Iterator();
          other.endIter = Iterator();
        }

        UnbufferedDataBlock<Iterator, Range>& operator=(UnbufferedDataBlock<Iterator, Range>&& other) {
          if (this != &other) {
            this->SuperType::operator=(std::forward(other));
            startIter = other.startIter;      other.startIter = Iterator();
            endIter = other.endIter;          other.endIter = Iterator();
          }
          return *this;
        }

        /// copy constructor and assignemtn operator
        UnbufferedDataBlock(const UnbufferedDataBlock<Iterator, Range>& other) : SuperType(other), startIter(other.startIter), endIter(other.endIter) {
        }

        UnbufferedDataBlock<Iterator, Range>& operator=(const UnbufferedDataBlock<Iterator, Range>& other) {
          if (this != &other) {
            this->SuperType::operator=(other);
            startIter = other.startIter;
            endIter = other.endIter;
          }
          return *this;
        }

        virtual ~UnbufferedDataBlock() {};



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

        typename SuperType::iterator begin()
        {
          return startIter;
        }
        typename SuperType::const_iterator cbegin() const
        {
          return startIter;
        }
        typename SuperType::iterator end()
        {
          return endIter;
        }
        typename SuperType::const_iterator cend() const
        {
          return endIter;
        }

      protected:
        Iterator startIter;   // if buffering, kept for future use.
        Iterator endIter;     // kept for future use.
    };

    /**
     * @brief << operator to write out UnbufferedDataBlock object's fields.
     * @tparam Iterator   Source data iterator type.
     * @tparam Range      Range data type
     * @tparam Container  container type for buffer.  defaults to std::vector.
     * @param[in/out] ost   output stream to which the content is directed.
     * @param[in]     db    BufferedDataBlock object to write out
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
