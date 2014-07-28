/**
 * @file		Partitioner.hpp
 * @ingroup bliss::partition
 * @author	tpan
 * @brief   contains several class that provide different logic for partitioning a range
 * @details contains block, cyclic, and demand driven (THREAD SAFE) partitioners
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef PARTITIONER_HPP_
#define PARTITIONER_HPP_

#include <stdexcept>
#include <atomic>
#include <config.hpp>
#include <partition/range.hpp>

namespace bliss
{
  /**
   * @namespace partition
   */
  namespace partition
  {

    /**
     * @class			Partitioner
     * @brief     base class to partition a range.
     * @details   operates without knowledge of the data that the range refers to,  so no "search" type of partitioning
     *
     *            each partition consists of 1 or more chunks.  chunks in a partition share the same conceptual partition id.
     *
     *            Simple types (essentially functinoids, use default copy constructor and assignment operator.  no need to move versions.
     *            uses the Curiously Recursive Template Pattern to enforce consistent API from base to derived classes.
     *
     *            uses default constructor, copy constructor, and move constructor.
     *
     *          Should work for continuous value ranges, but this is not tested.
     *
     *
     * @tparam  Range     Range object to be partitioned.
     * @tparam  Derived   Subclass, used to specialize the Base class for deciding the private impl of public api.
     */
    template<typename Range, typename Derived>
    class Partitioner {
      protected:

        /**
         * @typedef ValueType
         * @brief   value type used by the range's start/end/overlap
         */
        using ValueType = typename Range::ValueType;

        /**
         * @typedef chunkSizeType
         * @brief   value type used for chunkSize.  if integral ype for valueType, then size_t.  else (floating point) same as ValueType
         */
        using ChunkSizeType = typename std::conditional<std::is_integral<ValueType>::value,
                                                        size_t, ValueType>::type;

        /**
         * @var src
         * @brief   range to be partitioned
         */
        Range src;

        /**
         * @var end
         * @brief   end of the range to be partitioned.
         */
        Range end;

        /**
         * @var nPartitions
         * @brief   number of partitions to divide the range into.  chunks in a partition have the same partition id.  1 or more partition is assigned to a caller.
         */
        size_t nPartitions;

        /**
         * @var chunkSize
         * @brief   size of each partition.  computed for block partitioner, and user specified for cyclic and demand driven.
         */
        ChunkSizeType chunkSize;

        static const uint8_t BEFORE = 0;
        static const uint8_t DURING = 1;
        static const uint8_t AFTER = 2;


      public:
        /**
         * @brief default destructor
         */
        virtual ~Partitioner() {};

        /**
         * @brief configures the partitioner with the source range, number of partitions.
         * @param _src          range object to be partitioned.
         * @param _nPartitions the number of partitions to divide this range into
         * @param _chunkSize   the number size of each partition.  default is 0, computed later.
         */
        void configure(const Range &_src, const size_t &_nPartitions, const ChunkSizeType &_chunkSize) {
          src = _src;
          end = _src;
          end.start = end.end;
          nPartitions = _nPartitions;
          chunkSize = _chunkSize;

        }

        /**
         * @brief get the next sub range (chunk) within a partition.
         * @details the function calls the specific subclass implementation.
         * @param partId    the id of the partition to retrieve
         * @return          range object containing the start and end of the next chunk in the partition.
         */
        inline Range& getNext(const size_t &partId) {
          // both non-negative.
          if (partId >= nPartitions) throw std::invalid_argument("ERROR: partitioner.getNext called with partition id larger than number of partitions.");

          return static_cast<Derived*>(this)->getNextImpl(partId);
        }

        /**
         * @brief   resets the partitioner to the initial condition.
         * @details useful for cyclic and demand-driven partitioner, where there is an internal counter pointing to the next chunk range.
         */
        void reset() {
          static_cast<Derived*>(this)->resetImpl();
        }

    };


    /**
     * @class BlockPartitioner
     * @brief A partition that creates equal sized partitions during division of the range.
     * @details   Each partition has a size that's guaranteed to be within 1 of each other in size.
     *            inherits from base Partitioner using CRTP, and implements the detail impl for getNext and reset
     * @tparam Range  the range type used.
     */
    template<typename Range>
    class BlockPartitioner : public Partitioner<Range, BlockPartitioner<Range> >
    {
      protected:
        /**
         * @var curr
         * @brief the result range.  cached here for speed
         */
        Range curr;

        /**
         * @var done
         * @brief boolean indicating that partition has already been done.
         */
        bool done;

        /**
         * @typedef ValueType
         */
        using ValueType = typename Range::ValueType;

        /**
         * @typedef ChunkSizeType
         */
        using ChunkSizeType = decltype(std::declval<Range>().size());

        /**
         * @var rem
         * @brief the left-over that is to be spread out to the first rem partitions.
         */
        ValueType rem;

        /**
         * @typedef BaseClassType
         * @brief   the superclass type.
         */
        using BaseClassType = Partitioner<Range, BlockPartitioner<Range> >;

      public:

        /**
         * @brief default destructor
         */
        virtual ~BlockPartitioner() {};

        /**
         * @brief configures the partitioner with the source range, number of partitions
         * @param _src          range object to be partitioned.
         * @param _nPartitions the number of partitions to divide this range into
         * @param _chunkSize  size of each chunk.  here it should be 0 so they will be auto computed.
         */
        void configure(const Range &_src, const size_t &_nPartitions, const ChunkSizeType &_chunkSize = 0) {
          this->BaseClassType::configure(_src, _nPartitions, _chunkSize);

          // compute the size of partitions and number of partitions that have size larger than others by 1.
          // play fpr chhunkSize to be 0 - rem would be more than 0, so the first rem partitions will have size of 1.
          this->chunkSize = this->src.size() / static_cast<ChunkSizeType>(this->nPartitions);
          // not using modulus since we ValueType may be float.
          rem = this->src.size() - this->chunkSize * static_cast<ChunkSizeType>(this->nPartitions);

          resetImpl();
        }

        /**
         * @brief       get the next chunk in the partition.  for BlockPartition, there is only 1 chunk.
         * @details     first call to this function returns the partitioned range.  subsequent calls get the end object (start == end)
         *              to call again and get meaningful result, reset() must be called.
         *              no need to really store much - everything can be recomputed relatively fast.
         * @param partId   partition id for the sub range.
         * @return      range of the partition
         */
        inline Range& getNextImpl(const size_t& partId) {
          // param validation already done.

          // if previously computed, then return end object, signalling no more chunks.
          if (done) return this->end;  // can only call this once.

          // if just 1 partition, return.
          if (this->nPartitions == 1) return this->src;

          // compute the subrange's start and end, spreading out the remainder to the first rem chunks/partitions.
          if (partId < rem)
          {
            // each chunk with partition id < rem gets 1 extra.
            curr.start = this->src.start + static_cast<ValueType>(partId) * (this->chunkSize + 1);
            curr.end = (this->src.end - curr.start > static_cast<ValueType>(this->chunkSize + 1) + this->src.overlap) ?
                curr.start + static_cast<ValueType>(this->chunkSize + 1) + this->src.overlap :
                this->src.end;
          }
          else
          {
            // first rem chunks have 1 extra element than chunkSize.  the rest have chunk size number of elements.
            curr.start = this->src.start + static_cast<ValueType>(partId) *  this->chunkSize + rem;
            curr.end = (this->src.end - curr.start > static_cast<ValueType>(this->chunkSize) + this->src.overlap) ?
                curr.start + static_cast<ValueType>(this->chunkSize) + this->src.overlap :
                this->src.end;
          }

          // last entry.  no overlap.  done via intersection
          //curr &= this->src;

          curr.block_start = curr.start;

          done = true;
          return curr;
        }

        /**
         * @brief reset size to full range, mark the done to false.
         */
        void resetImpl() {
          curr = this->src;
          done = false;
        }

    };

    /**
     * @class CyclicPartitioner
     * @brief cyclically partition the range into chunkSize blocks.
     *        uses CRTP pattern to enforce consistent API.
     * @tparam Range  type of range to be partitioned.
     */
    template<typename Range>
    class CyclicPartitioner : public Partitioner<Range, CyclicPartitioner<Range> >
    {

      protected:
        /**
         * @typedef BaseClassType
         * @brief   the superclass type.
         */
        using BaseClassType = Partitioner<Range, CyclicPartitioner<Range> >;

        /**
         * @typedef ChunkSizeType
         */
        using ChunkSizeType = typename BaseClassType::ChunkSizeType;


        /**
         * @var done
         * @brief An array of "state", one for each partition.
         */
        uint8_t *state;

        /**
         * @var curr
         * @brief An array of "curr" subranges, one for each partition.  updated as getNext is called.
         */
        Range *curr;

        /**
         * @var nChunks
         * @brief   number of chunks in the src range.
         * @details by comparing nChunks to nPartitions (user specified) we can tell if some partitions will not be receiving a chunk.
         *          this is a mechanism to deal with nPartition that are too large.
         */
        size_t nChunks;

        /**
         * @var stride
         * @brief size of a stride, which is chunkSize *  number of partitions (stride for successive call to getNext)
         */
        ChunkSizeType stride;



      public:

        /**
         * @brief default destructor.  cleans up arrays for callers.
         */
        virtual ~CyclicPartitioner() {
          if (state) delete [] state;
          if (curr) delete [] curr;
        }

        /**
         * @brief configures the partitioner with the source range, number of partitions
         * @param _src          range object to be partitioned.
         * @param _nPartitions the number of partitions to divide this range into
         * @param _chunkSize    Size of the chunks
         *
         */
        void configure(const Range &_src, const size_t &_nPartitions, const ChunkSizeType &_chunkSize) {
          this->BaseClassType::configure(_src, _nPartitions, _chunkSize);

          // get the maximum number of chunks in the source range.
          nChunks = std::floor((this->chunkSize - 1 + this->src.size()) / this->chunkSize);

          // stride:  if there are more chunks than partition, then set real stride, else make it src size.
          stride = (nChunks > this->nPartitions) ? this->chunkSize * static_cast<ChunkSizeType>(this->nPartitions) : this->src.size();

          state = new uint8_t[std::min(nChunks, this->nPartitions)];
          curr = new Range[std::min(nChunks, this->nPartitions)];

          resetImpl();
        }

        /**
         * @brief       get the next chunk in the partition.  for CyclickPartition, keep getting until done.
         * @details     each call to this function gets the next chunk belonging to this partition.
         *              for cyclicPartitioner, the chunks are separated by chunkSize * nPartitions
         * @param partId   partition id for the sub range.
         * @return      range of the partition
         */
        inline Range& getNextImpl(const size_t& partId) {

          // if nChunks < nPartitions, then each partId will only get a chunk once.
          // in that case, if partId is >= nChunks, then end is returned.
          if (partId >= nChunks) return this->end;

          // if this partition is done, return the last entry (end)
          if (state[partId] == BaseClassType::AFTER)  return this->end;

          /// first iteration, use initialized value
          if (state[partId] == BaseClassType::BEFORE) {
            state[partId] = BaseClassType::DURING;
            return curr[partId];
          }
          // else not the first and not last, so increment.

          /// comparing to amount of range available - trying to avoid data type overflow.

          // shift starting position by stride length
          if (this->src.end - curr[partId].start >= stride) {
            // has room.  so shift
            curr[partId].start += stride;
          } else {
            state[partId] = BaseClassType::AFTER;         // if outside range, done
            return this->end;
          }

          // shift end position by stride length - start is NOT outside range.
          // overlap is already part of curr[partId].end set from resetImpl.
          if (this->src.end - curr[partId].end > stride) {
            // end is not outside parent range.
            curr[partId].end += stride;
          } else {
            // end is outside parent range
            curr[partId].end = this->src.end;
            state[partId] = BaseClassType::AFTER;
          }

          return curr[partId];
        }

        /**
         * @brief resets the partitioner by resetting the internal done and subrange arrays.
         * @details this function also serves to initialize the arrays.
         */
        void resetImpl()
        {
          size_t s = std::min(nChunks, this->nPartitions);
          for (size_t i = 0; i < s; ++i) {
            state[i] = BaseClassType::BEFORE;
            curr[i].start = this->src.start + static_cast<ChunkSizeType>(i) * this->chunkSize;
            // if it's the last chunk (s-1), and number of chunks is less or equal to number of partitions
            // so if i = nChunks - 1
            // do this check to avoid overflow.
            curr[i].end = (i == nChunks - 1) ?
                this->src.end :
                curr[i].start + this->chunkSize + this->src.overlap;
            curr[i].block_start = curr[i].start;
            curr[i].overlap = this->src.overlap;

//            printf ("values start-end: %lu %lu\n", curr[i].start, curr[i].end);

          }
        }


    };


    /**
     * @class DemandDrivenPartitioner
     * @brief a partitioner that assigns chunks to partition in the order that the getNext function is called.
     * @tparam Range  type of the range object to be partitioned.
     */
    template<typename Range>
    class DemandDrivenPartitioner : public Partitioner<Range, DemandDrivenPartitioner<Range> >
    {

      protected:
        /**
         * @typedef BaseClassType
         * @brief   the superclass type.
         */
        using BaseClassType = Partitioner<Range, DemandDrivenPartitioner<Range> >;

        /**
         * @typedef ChunkSizeType
         */
        using ChunkSizeType = typename BaseClassType::ChunkSizeType;

        /**
         * @var chunkOffset
         * @brief the offset for the next chunk to be returned.  TREHAD SAFE
         */
        std::atomic<ChunkSizeType> chunkOffset;

        /**
         * @var id of current chunk
         * @brief  used for tracking the chunk begin returned.  useful if nChunk < nPartition
         */
        std::atomic<size_t> chunkId;

        /**
         * @var nChunks
         * @brief   number of chunks in the src range.
         * @details by comparing nChunks to nPartitions (user specified) we can tell if some partitions will not be receiving a chunk.
         *          this is a mechanism to deal with nPartition that are too large.
         */
        size_t nChunks;

        /**
         * @var curr
         * @brief internal cache of subrange objects that each partition last retrieved.
         */
        Range *curr;

        /**
         * @done
         * @brief boolean indicating that the range has been exhausted by the traversal.  THREAD SAFE
         */
        std::atomic<bool> done;

        /**
         * @typedef ValueType
         * @brief   type for the start/end/overlap
         */
        using ValueType = typename Range::ValueType;


      public:
        /**
         * @brief default destructor
         */
        virtual ~DemandDrivenPartitioner() {
          if (curr != nullptr) delete [] curr;
        }


        /**
         * @brief configures the partitioner with the source range, number of partitions
         * @param _src          range object to be partitioned.
         * @param _nPartitions  the number of partitions to divide this range into
         * @param _chunkSize  size of each chunk for the partitioning.
         *
         */
        void configure(const Range &_src, const size_t &_nPartitions, const ChunkSizeType &_chunkSize) {
          this->BaseClassType::configure(_src, _nPartitions, _chunkSize);

          // get the maximum number of chunks in the source range.
          nChunks = std::floor((this->chunkSize - 1 + this->src.size()) / this->chunkSize);

          // if there are more partitions than chunks, then the array represents mapping from chunkId to subrange
          // else if there are more chunks than partitions, then the array represents the most recent chunk assigned to a partition.
          curr = new Range[std::min(nChunks, this->nPartitions)];

          resetImpl();
        };

///        TODO:

        template<typename T = ChunkSizeType>
        typename std::enable_if<std::is_integral<T>, ValueType>::type getNextOffset() {

        }
        template<typename T = ChunkSizeType>
        typename std::enable_if<std::is_floating_point<T>, ValueType>::type getNextOffset() {

        }

        /**
         * @brief       get the next chunk in the partition.  for DemandDrivenPartition, keep getting until done.
         * @details     each call to this function gets the next chunk in the range and assign it to a partition.
         *              the sequence of partition ids depends on call order.
         *
         *              ASSUMPTION:  no 2 concurrent callers will be requesting the same partition Id.
         *
         *              THREAD SAFE
         * @param partId   partition id for the sub range.
         * @return      range of the partition
         */
        inline Range& getNextImpl(const size_t& partId) {
          // all done, so return end
          if (done.load(std::memory_order_consume)) return this->end;

          ValueType s = chunkOffset.fetch_add(this->chunkSize, std::memory_order_acq_rel);

          /// identify the location in array to store the result
          // first get the id of the chunk we are returning.
          size_t id = chunkId.fetch_add(1, std::memory_order_acq_rel);
          // if there are more partitions than chunks, then the array represents mapping from chunkId to subrange
          // else if there are more chunks than partitions, then the array represents the most recent chunk assigned to a partition.
          id = (nChunks < this->nPartitions ? id : partId);

          if (s >= this->src.end) {
            done.store(true, std::memory_order_release);
            return this->end;
          } else {

            curr[id].block_start =
                curr[id].start = s;
            curr[id].end = (this->src.end - s > this->chunkSize + this->src.overlap) ?
                s + this->chunkSize + this->src.overlap : this->src.end;

            return curr[id];
          }
        }

        /**
         * @brief resets the partitioner by resetting the internal subrange arrays.  also reset the offset to the start of the src range.
         * @details this function also serves to initialize the subrange array.
         */
        void resetImpl() {
          printf("RESET: number of chunks, partitions: %lu, %lu\n", nChunks, this->nPartitions);


          // these 2 calls probably should be synchronized together.
          chunkOffset.store(this->src.start, std::memory_order_release);
          done.store(false, std::memory_order_release);

          size_t s = std::min(this->nPartitions, nChunks);
          for (size_t i = 0; i < s; ++i) {
            curr[i] = Range(this->src.start, this->src.start, this->src.overlap);
          }
        }


    };


  } /* namespace partition */
} /* namespace bliss */

#endif /* PARTITIONER_HPP_ */
