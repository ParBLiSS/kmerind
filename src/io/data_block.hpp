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
     * @class     bliss::io::DataBlock
     * @brief     abstraction to represent a block of data.
     * @details   buffer has to be allocated externally to the right size.  This is to allow reusing a buffer by the calling code.
     *            constructor determines whether buffering or not.  (could have done it as template, but that is not as flexible during runtime.
     *            constructors have _end in case Iterator is not a random access iterator so can't do _start + n
     */
    template<typename Iterator, typename Range>
    class DataBlock
    {
      public:
        /**
         * constructor for non-buffering
         * @param _start
         * @param _end
         * @param _range
         */
        DataBlock(const Iterator &_start, const Iterator &_end, const Range &_range) :
          range(_range), srcStartIter(_start), srcEndIter(_end), destStartIter(_start), destEndIter(_end), buffering(false)
        {
        }

        /**
         * constructor for buffering.  buffer is from caller.  needs to be allocated properly.
         * @param _start
         * @param _end
         * @param _range
         * @param _destStart
         * @param destLen
         */
        DataBlock(const Iterator &_start, const Iterator &_end, const Range &_range, Iterator _destStart, const size_t destLen) :
          range(_range), srcStartIter(_start), srcEndIter(_end), destStartIter(_destStart), buffering(true) {

          assert((range.end - range.start) <= destLen);

          destEndIter = std::copy(srcStartIter, srcEndIter, destStartIter);
        }

        virtual ~DataBlock() {
        }

        const Iterator& begin() const {
          return destStartIter;
        }
        const Iterator& end() const {
          return destEndIter;
        }
        const Range& getRange() const {
          return range;
        }

      protected:
        const Range range;
        const Iterator &srcStartIter;   // if buffering, kept for future use.
        const Iterator &srcEndIter;     // kept for future use.
                                        // reference so will update with caller's values

        Iterator destStartIter;
        Iterator destEndIter;

        const bool buffering;
    };




  } /* namespace io */
} /* namespace bliss */

#endif /* DATA_BLOCK_HPP_ */
