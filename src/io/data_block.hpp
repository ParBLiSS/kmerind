/**
 * @file		data_block.hpp
 * @ingroup io
 * @author	Tony Pan <tpan7@gatech.edu>
 * @brief   contains classes for managing an in-memory 1D address range, along with iterators pointing to the actual data.  with or without buffering
 * @details Presents a standardized interface for buffered and unbuffered data blocks.
 *
 *          Given a range and the corresponding raw data, the goal of DataBlock is to
 *          present a consistent iterator-based interface for traversing the data,
 *          with or without buffering.  The challenge is when buffering, data is copied into a
 *          container, which will have a different data type than the source data's iterator type.
 *
 *          DataBlock's design requirements are listed below:
 *           1. stores iterators (pointer, normal iterators) for input and output traversal of the data
 *           2. consistent API for iterator accessor
 *                a. iterators for either original data, or buffered data accessed through begin()/end(),etc.
 *           3. manages the buffer memory internally and transparently, if any
 *                a. Okay not to be thread safe.
 *                b. reuse the buffer memory when a BufferedBlock is remapped to another range.
 *           4. store the conceptual "range" of the data, e.g. offsets within the entire available data sequence.
 *
 *          An example usage is to represent portions of a memory mapped file.  When buffered, the iterator
 *          type is no longer pointer, but rather std::vector::iterator.
 *
 *          All derived/specialized classes for buffering and unbuffering DataBlock need to support the following functions
 *            begin/end
 *            cbegin/cend
 *            assign/clear
 *            hasBuffer
 *            isEmpty (common)
 *            getRange (common)
 *
 *          5 design approaches were considered to support flexible container choice, consistent api, and clear semantic meaning of classes.
 *          Issue to consider:  when a container is used for buffering, the iterator type probably is not the same as the input iterator.
 *          1. a single DataBlock class with runtime choice of buffering
 *              problem:  would require iterator accessors to take a parameter to enable overloading (for the different return types)
 *                  logic would be complicated
 *              benefit:  runtime choice of overloading.
 *          2. a single DataBlock class with template specialization for buffer on/off.
 *              problem:  use enable_if to choose the iterator accessor return types
 *                  either a boolean template paramter will be needed, or the container type of void means no buffering
 *                  a large percentage of functions (first 4 sets above) will need to be specialized with enable_if, thus replicating declarations
 *          3. template specialization with an empty default, and 2 specialized classes
 *              problem:  same requirement that there is some way to distinguish between buffering and not buffering
 *                  replicated functions in different specializations
 *              benefit:  fewer enable_if.
 *                  easily specialize for new container types (e.g. raw arrays)
 *          4. template specialization with base being the BufferedDataBlock type and specialization being UnBufferedDataBlock
 *              problem and benefits :  same as 3
 *              problem:  new parameter types (e.g. container type) may be incorrectly handled by default.
 *          5. CRTP with DataBlock being the base and BufferedDataBlock and UnbufferedDataBlock as derived classes, container type as template
 *              problem: more complex hierarchy.  harder to read and understand.
 *              benefit: common member and methods are in the base class, enforcing consistent method naming, etc.
 *
 *  CHOSE: 4, with buffering indicated by the type of buffer container, and default buffer container is std::vector (using alias to make this clear)
 *    Use type aliasing to make semantically easy-to-understand DataBlock class names: BufferedDataBlock and UnBufferedDataBlock
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef DATA_BLOCK_HPP_
#define DATA_BLOCK_HPP_

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
     * @brief   BufferedDataBlock uses an internally allocated container to mirror the source data.
     * @details Source data type may be pointer or a proper iterator type.  The container type then determine the output iterator type
     *          Provides iterator accessor functions that point to the container's iterator accessor functions.
     * @tparam Iterator   Source data iterator type.
     * @tparam Range      Range Type within the overall data element sequence that this block is wrapping.
     * @tparam Container  container type for buffer.  defaults to std::vector.
     */
    template<typename Iterator, typename Range, typename Container>
    class DataBlock {
      protected:
        // assert that ValueType is defined.
        static_assert(!(std::is_same<typename std::iterator_traits<Iterator>::value_type, void>::value), "Iterator is NOT valid.");

        // assert that container is valid.
        static_assert(bliss::utils::container_traits<Container>::isAssignable, "ERROR: Container needs to be a sequence container with 'assign' method.");
        static_assert(bliss::utils::container_traits<Container>::isConstIterable, "ERROR: Container needs to be a sequence container with 'cbegin' and 'cend' methods.");
        static_assert(bliss::utils::container_traits<Container>::isIterable, "ERROR: Container needs to be a sequence container with 'begin' and 'end' methods.");

        /**
         * @typedef   ValueType
         * @brief     type of the data elements in the DataBlock
         */
        typedef typename std::iterator_traits<Iterator>::value_type           ValueType;

      public:
        /**
         * @typedef iterator
         * @brief   type of output iterator
         */
        typedef decltype(std::declval<Container>().begin())             iterator;

        /**
         * @typedef const_iterator
         * @brief   const version of output iterator
         */
        typedef decltype(std::declval<Container>().cbegin())            const_iterator;

        // assert that the template parameter's iterator (input) and the local iterator (output) have same value type.
        static_assert(std::is_same<typename std::iterator_traits<Iterator>::value_type,
                      typename std::iterator_traits<iterator>::value_type>::value,
                      "Iterator and Container use different value types.");


        /**
         *  @brief  default constructor
         */
        DataBlock() : range(), empty(true), buffer() {};

        /**
         * @brief         move constructor.  moves range and empty flag.
         * @param other   source DataBlock to move
         */
        DataBlock(DataBlock<Iterator, Range, Container>&& other) : range(other.range), empty(other.empty), buffer(std::move(other.buffer)) {
          other.range = Range();
          other.empty = true;
          other.buffer.clear();
        };

        /**
         * @brief         move assignment operator.  moves range and empty flag.
         * @param other   source DataBlock to move
         * @return        self, updated
         */
        DataBlock<Iterator, Range, Container>& operator=(DataBlock<Iterator, Range, Container>&& other) {
          if (this != &other) {
            range = other.range;  other.range = Range();
            empty = other.empty;  other.empty = true;
            buffer.swap(other.buffer);
            other.buffer.clear();
          }

          return *this;
        };

        /**
         * @brief         copy constructor.  copies range and empty flag.
         * @param other   source DataBlock to copy
         */
        DataBlock(const DataBlock<Iterator, Range, Container>& other) : range(other.range), empty(other.empty), buffer(other.buffer) {};

        /**
         * @brief         copy assignment operator.  copies range and empty flag.
         * @param other   source DataBlock to copy
         * @return        self, updated
         */
        DataBlock<Iterator, Range, Container>& operator=(const DataBlock<Iterator, Range, Container>& other) {
          if (this != &other) {
            range = other.range;
            empty = other.empty;
            buffer = other.buffer;
          }
          return *this;
        };


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
          return true;
        }

        /**
         * @brief    assigns some data to the DataBlock.  check if empty.  if not, then call derived class to assign.
         * @param _start    begin iterator  (input)
         * @param _end      end iterator    (input)
         * @param _range    the corresponding range (offset values)
         */
        void assign(const Iterator &_start, const Iterator &_end, const Range &_range) {
          if (_start == _end || _range.size() == 0) {
            clear();
          } else {
            empty = false;
            range = _range;
            buffer.assign(_start, _end);
          }
        }

        /**
         * @brief   clears the DataBlock by setting it to empty, clear the Range, and call the derived class's clear implementation.
         */
        void clear() {
          empty = true;
          range = Range();
          buffer.clear();
        }

        /**
         * @brief   begin iterator accessor
         * @details points to the start of the internal container.
         *          same signature as for std containers.
         * @return  start iterator for the data in memory
         */
        iterator begin()
        {
          return buffer.begin();
        }

        /**
         * @brief   constant begin iterator accessor
         * @details points to the start of the internal container.
         *          same signature as for std containers.
         * @return  const start iterator for the data in memory
         */
        const_iterator cbegin() const
        {
          return buffer.cbegin();
        }

        /**
         * @brief   constant end iterator accessor
         * @details points to the end of the internal container.
         *          same signature as for std containers.
         * @return  const end iterator for the data in memory
         */
        iterator end()
        {
          return buffer.end();
        }

        /**
         * @brief   constant end iterator accessor
         * @details points to the end of the internal container.
         *          same signature as for std containers.
         * @return  const end iterator for the data in memory
         */
        const_iterator cend() const
        {
          return buffer.cend();
        }


      protected:
        /// range object describing the start and end of the data block, in 1D coordinates (i.e. offsets)
        Range range;

        /// boolean flag indicating that the datablock does not contain valid data.
        bool empty;

        /// container used as buffer to hold a copy of the input data
        Container buffer;

    };


    /**
     * @typedef VectorBuffer
     * @brief   aliased to std vector.  when this is given to DataBlock, vector is used for buffering.
     */
    template<class T, class A = std::allocator<T> > using VectorBuffer = std::vector<T, A>;


    /**
     * @typedef BufferedDataBlock
     * @brief   aliased type name for Buffered DataBlock, for convenience.
     */
    template<typename Iterator, typename Range> using BufferedDataBlock = DataBlock<Iterator, Range, VectorBuffer<typename std::iterator_traits<Iterator>::value_type> >;





    /**
     * @typedef NoBuffer
     * @brief   empty type.  when this is given to DataBlock, no buffering is done.
     */
    struct NoBuffer {};

    /**
     * @class   DataBlock
     * @brief   UnBuffered DataBlock provides iterator accessors to the input iterators.
     * @details Source data type may be pointer or a proper iterator type.  The input iterators are passed through as output iterators.
     *
     * @tparam  Iterator        Iterator type for the source data
     * @tparam  Range           Range within the overall data element sequence that this block is wrapping.
     */
    template<typename Iterator, typename Range>
    class DataBlock<Iterator, Range, NoBuffer> {
      protected:
        // assert that ValueType is defined.
        static_assert(!(std::is_same<typename std::iterator_traits<Iterator>::value_type, void>::value), "Iterator is NOT valid.");

        /**
         * @typedef   ValueType
         * @brief     type of the data elements in the DataBlock
         */
        typedef typename std::iterator_traits<Iterator>::value_type           ValueType;

      public:
        /**
         * @typedef iterator
         * @brief   type of output iterator
         */
        typedef typename std::remove_reference<Iterator>::type                iterator;

        /**
         * @typedef const_iterator
         * @brief   const version of output iterator
         */
        typedef typename std::add_const<iterator>::type                       const_iterator;

        /**
         *  @brief  default constructor
         */
        DataBlock() : range(), empty(true), startIter(), endIter() {};

        /**
         * @brief         move constructor.  moves range and empty flag.
         * @param other   source DataBlock to move
         */
        DataBlock(DataBlock<Iterator, Range, NoBuffer>&& other) = default;

        /**
         * @brief         move assignment operator.  moves range and empty flag.
         * @param other   source DataBlock to move
         * @return        self, updated
         */
        DataBlock<Iterator, Range, NoBuffer>& operator=(DataBlock<Iterator, Range, NoBuffer>&& other) = default;

        /**
         * @brief         copy constructor.  copies range and empty flag.
         * @param other   source DataBlock to copy
         */
        DataBlock(const DataBlock<Iterator, Range, NoBuffer>& other) = default;

        /**
         * @brief         copy assignment operator.  copies range and empty flag.
         * @param other   source DataBlock to copy
         * @return        self, updated
         */
        DataBlock<Iterator, Range, NoBuffer>& operator=(const DataBlock<Iterator, Range, NoBuffer>& other) = default;

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
          return false;
        }

        /**
         * @brief    assigns some data to the DataBlock.  check if empty.  if not, then call derived class to assign.
         * @param _start    begin iterator  (input)
         * @param _end      end iterator    (input)
         * @param _range    the corresponding range (offset values)
         */
        void assign(const Iterator &_start, const Iterator &_end, const Range &_range) {
          if (_start == _end || _range.size() == 0) {
            clear();
          } else {
            empty = false;
            range = _range;

            startIter = _start;
            endIter = _end;
          }
        }

        /**
         * @brief   clears the DataBlock by setting it to empty, clear the Range, and call the derived class's clear implementation.
         */
        void clear() {
          empty = true;
          range = Range();
          startIter = endIter = Iterator();
        }

        /**
         * @brief   begin iterator accessor
         * @details pass through the input iterator because we are not buffering
         *          same signature as for std containers.
         * @return  start iterator for the data in memory
         */
        iterator begin()
        {
          return startIter;
        }

        /**
         * @brief   constant begin iterator accessor
         * @details pass through the input iterator because we are not buffering
         *          same signature as for std containers.
         * @return  const start iterator for the data in memory
         */
        const_iterator cbegin() const
        {
          return startIter;
        }

        /**
         * @brief   constant end iterator accessor
         * @details pass through the input iterator because we are not buffering
         *          same signature as for std containers.
         * @return  const end iterator for the data in memory
         */
        iterator end()
        {
          return endIter;
        }

        /**
         * @brief   constant end iterator accessor
         * @details pass through the input iterator because we are not buffering
         *          same signature as for std containers.
         * @return  const end iterator for the data in memory
         */
        const_iterator cend() const
        {
          return endIter;
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

        /**
         * @var startIter
         * @brief the input begin iterator
         */
        Iterator startIter;

        /**
         * @var endIter
         * @brief the input end iterator
         */
        Iterator endIter;
    };

    /**
     * @typedef UnbufferedDataBlock
     * @brief   aliased type name for Unbuffered DataBlock, for convenience.
     */
    template<typename Iterator, typename Range> using UnbufferedDataBlock = DataBlock<Iterator, Range, NoBuffer>;





    /**
     * @brief << operator to write out DataBlock object's actual data.
     * @tparam Iterator   Source data iterator type.
     * @tparam Range      Range data type
     * @tparam Container  container type for buffer.  defaults to std::vector.
     * @param[in/out] ost   output stream to which the content is directed.
     * @param[in]     db    BufferedDataBlock object to write out
     * @return              output stream object
     */
    template<typename Iterator, typename Range, typename BufferContainer>
    std::ostream& operator<<(std::ostream& ost, const DataBlock<Iterator, Range, BufferContainer>& db)
    {
      std::ostream_iterator<typename std::iterator_traits<decltype(db.cbegin())>::value_type > oit(ost);
      std::copy(db.cbegin(), db.cend(), oit);
      return ost;
    }



  } /* namespace io */
} /* namespace bliss */

#endif /* DATA_BLOCK_HPP_ */
