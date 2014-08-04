/**
 * @file		Partitioner.hpp
 * @ingroup bliss::partition
 * @author	tpan
 * @brief   contains several class that provide different logic for partitioning a range
 * @details contains block, cyclic, and demand driven (THREAD SAFE) partitioners
						logic implementation uses comparison to avoid overflows and implicit casting where needed.
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
     *            PARTIONER essentially divide the RANGE into CHUNKS, then assign PARTITION id to the chunks.
     *            each partition consists of 0 or more chunks.  chunks in a partition share the same (implicit) partition id.
     *
     *            Simple types (essentially functinoids), use default copy constructor and assignment operator.  no need to move versions.
     *            uses the Curiously Recursive Template Pattern to enforce consistent API from base to derived classes and to provide a path
     *            for deriving new classes

     *            uses default constructor, copy constructor, and move constructor.
     *
     *            Works for continuous value ranges.

     *            Partitioners are designed for multithreaded usage.  For Block and Cyclic partitioners, the chunks are computed
     *            independent of other threads.  For DemandDriven partitioner, the chunks are computed using atomc operation to provide
     *            memory fence/synchronization.
     *
     * @note      The purpose of partitioning is to divide the data for computation.  Partitioner will therefore divide the source range,
     *            including the overlap region (see range class documentation), into equal parts.
     *            A new overlap region length can be specified during partitioning.  When overlap region length is specified, all subranges
     *            from the partitioner, except the last subrange, will have its "overlap" variable set to the overlap region size.  The last
     *            subrange has overlap region of size 0.
     *
     *            Rationale for NOT implementing Iterator interface:
     *              Partitioners support multi-threaded access.  For Block and Cyclic, no interthread communication is required.
     *              ThreadId is needed for caching but we can instantiate one partitioner per thread with the appropriate threadId
     *
     *              For DemandDriven, interthread communication is needed via memory fence/synchronization.  We therefore need a
     *              single instance that is shared between all threads, else coordination is not possible.  While computing and returning nextChunk
     *              is not a problem, threadId will need to be provided per iterator call to allow caching, which significantly
     *              deviates from std iterator interface.  Alternatively, the partitioner will need to use the thread library's call
     *              to get thread id, per call to the increment function in the iterator interface.
     *
     *              it MAY be possible to use std::bind to avoid this, but it probably will be complicated.
     *
     *              For consistency, simplicity of interface, and performance, Partitioner interface therefore does NOT follow std::iterator convention.
     *
     *
     * @tparam  Range     Range object to be partitioned.
     * @tparam  Derived   Subclass, used to specialize the Base class for deciding the private impl of public api.
     */
    template<typename Range, typename Derived>
    class Partitioner {
      protected:

        /**
         * @typedef RangeValueType
         * @brief   value type used by the range's start/end/overlap
         */
        using RangeValueType = typename Range::ValueType;

        /**
         * @typedef chunkSizeType
         * @brief   value type used for chunkSize.  if integral ype for valueType, then size_t.  else (floating point) same as RangeValueType
         */
        using ChunkSizeType = typename std::conditional<std::is_integral<RangeValueType>::value,
                                                        size_t, RangeValueType>::type;

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
         * @note    this value excludes overlap region size.
         */
        ChunkSizeType chunkSize;

        /**
         * @var     overlapSize
         * @brief   the overlapSize that each subrange should contain.  Parameter of this partitioning.
         */
        RangeValueType overlapSize;

        /**
         * @brief   computes the number of chunks in a src range.  for Integral type only.
         * @details computed using the nonoverlapping part of each subrange, length in chunkSize
         * @tparam  R  range type.  used for enable_if.
         * @return  number of chunks in a range.
         */
        template <typename R = Range>
        typename std::enable_if<std::is_integral<typename R::ValueType>::value, size_t>::type computeNumberOfChunks() {
          return std::floor((this->chunkSize - 1 + this->src.size()) / this->chunkSize);
        }
        /**
         * @brief   computes the number of chunks in a src range, excluding the overlap region.  for floating type only.
         * @tparam  R  range type.  used for enable_if.
         * @return  number of chunks in a range, based on chunkSize.
         */
        template <typename R = Range>
        typename std::enable_if<std::is_floating_point<typename R::ValueType>::value, size_t>::type computeNumberOfChunks() {
          return this->src.size() / this->chunkSize;
        }

        /**
         * @brief   internal function to compute the range for a chunk of the parent range.
         * @details The chunk range is computed by first setting the start to parent range start + start offset (e.g. rem from BlockPartitioner)
         * @note    chunkId must not exceed maximum number of chunks - enforced outside of this function.
         *
         * @param[out] r            the range to compute.  initialized in the function to pr.start + start_offset
         * @param[in] pr            the parent range
         * @param[in] start_offset  offset to add to the start position
         * @param[in] chunkId       the chunkId of the current range.
         * @param[in] chunk_size    the size of the chunk        (precomputed or user specified)
         * @param[in] overlap_size    the size of the overlap region (user specified)
         */
        void computeRangeForChunkId(Range & r, const Range & pr, const RangeValueType& start_offset,
                                    const size_t& chunkId, const ChunkSizeType& chunk_size, const RangeValueType& overlap_size) {
          // compute start
          r.start = pr.start + start_offset;
          if (pr.end - r.start > static_cast<RangeValueType>(chunkId) *  chunk_size) {
            r.start += static_cast<RangeValueType>(chunkId) *  chunk_size;
          } else {
            r.start = pr.end;
          }

          if (pr.end - r.start > static_cast<RangeValueType>(chunk_size) + overlap_size) {
            // far away from parent range's end
            r.end = r.start + static_cast<RangeValueType>(chunk_size) + overlap_size;
            r.overlap = overlap_size;
          } else
            r.end = pr.end;

            if (pr.end - r.start <= static_cast<RangeValueType>(chunk_size)) {
              // parent range's end is within the target subrange's chunk region
              r.overlap = 0;
            } else {
              // parent range's end is within the target subrange's overlap region
              r.overlap = pr.end - r.start - static_cast<RangeValueType>(chunk_size);
            }
        }

      public:
        Partitioner() : src(), end(), nPartitions(0), chunkSize(0), overlapSize(0) {};

        /**
         * @brief default destructor
         */
        virtual ~Partitioner() {};

        /**
         * @brief configures the partitioner with the source range, number of partitions.
         * @param _src          range object to be partitioned.
         * @param _nPartitions the number of partitions to divide this range into
         * @param _chunkSize   the size of each partition.  default is 0 (in subclasses), computed later.
         * @param _overlapSize   the size of the overlap overlap region.
         */
        void configure(const Range &_src, const size_t &_nPartitions, const ChunkSizeType &_chunkSize, const RangeValueType &_overlapSize) {

          if (_nPartitions == 0)
            throw std::invalid_argument("ERROR: partitioner c'tor: nPartitions is 0");
          nPartitions = _nPartitions;

          if (_overlapSize < 0)
            throw std::invalid_argument("ERROR: partitioner c'tor: overlapSize is < 0");
          overlapSize = _overlapSize;

          if (_chunkSize < 0)
            throw std::invalid_argument("ERROR: partitioner c'tor: chunkSize is < 0");
          chunkSize = _chunkSize;

          src = _src;  // leaving the overlap region size values same as before for src.
          end = Range(src.end, src.end, 0);
        }

        /**
         * @brief get the next sub range (chunk) within a partition.
         * @details the function calls the specific subclass implementation.
         * @param partId    the id of the partition to retrieve
         * @return          range object containing the start and end of the next chunk in the partition.
         */
        inline Range& getNext(const size_t &partId) {
          // both non-negative.
          if (partId >= nPartitions)
            throw std::invalid_argument("ERROR: partitioner.getNext called with partition id larger than number of partitions.");

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
         * @typedef BaseClassType
         * @brief   the superclass type.
         */
        using BaseClassType = Partitioner<Range, BlockPartitioner<Range> >;

        /**
         * @typedef ChunkSizeType
         */
        using ChunkSizeType = decltype(std::declval<Range>().size());

        /**
         * @typedef ValueType
         */
        using RangeValueType = typename BaseClassType::RangeValueType;

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
         * @var rem
         * @brief the left-over that is to be spread out to the first rem partitions.
         */
        RangeValueType rem;


      public:

        BlockPartitioner() : BaseClassType(), curr(), done(false), rem(0) {};


        /**
         * @brief default destructor
         */
        virtual ~BlockPartitioner() {};

        /**
         * @brief configures the partitioner with the source range, number of partitions
         * @param _src          range object to be partitioned.
         * @param _nPartitions the number of partitions to divide this range into
         * @param _chunkSize  size of each chunk.  here it should be 0 so they will be auto computed.
         * @param _overlapSize   the size of the overlap overlap region.
         *
         */
        void configure(const Range &_src, const size_t &_nPartitions,
                       const ChunkSizeType &_chunkSize = 0, const RangeValueType &_overlapSize = 0) {
          this->BaseClassType::configure(_src, _nPartitions, _chunkSize, _overlapSize);

          // compute the size of partitions and number of partitions that have size larger than others by 1.
          // if chunkSize to be 0 - rem would be more than 0, so the first rem partitions will have size of 1.

          // first compute the chunkSize excluding the overlap region
          this->chunkSize = this->src.size() / static_cast<ChunkSizeType>(this->nPartitions);

          // next figure out the remainder using chunkSize that excludes overlap region size.
          // not using modulus since we RangeValueType may be float.
          rem = this->src.size() - this->chunkSize * static_cast<ChunkSizeType>(this->nPartitions);

          resetImpl();
        }

        /**
         * @brief       get the next chunk in the partition.  for BlockPartition, there is only 1 chunk.  INTEGRAL Type Only
         * @details     first call to this function returns the partitioned range.  subsequent calls get the end object (start == end)
         *              to call again and get meaningful result, reset() must be called.
         *              no need to really store much - everything can be recomputed relatively fast.
         * @param partId   partition id for the sub range.
         * @return      range of the partition
         */
        template<typename R = Range>
        typename std::enable_if<std::is_integral<typename R::ValueType>::value, Range>::type& getNextImpl(const size_t& partId) {
          // param validation already done.

          // if previously computed, then return end object, signalling no more chunks.
          if (done) return this->end;  // can only call this once.

          // if just 1 partition, return.
          if (this->nPartitions == 1) return this->src;

          // compute the subrange's start and end, spreading out the remainder to the first rem chunks/partitions.
          ChunkSizeType cs = this->chunkSize;
         RangeValueType startOffset = 0;
          if (partId < rem)
          {
            // each chunk with partition id < rem gets 1 extra.
            cs += 1;
          }
          else
          {
            // first rem chunks have 1 extra element than chunkSize.  the rest have chunk size number of elements.
           startOffset = rem;
          }

          // compute the new range
          BaseClassType::computeRangeForChunkId(curr, this->src, startOffset, partId, cs, this->overlapSize);

          // block partitioning only allows calling this once.
          done = true;
          return curr;
        }


        /**
         * @brief       get the next chunk in the partition.  for BlockPartition, there is only 1 chunk.  Floating Point only
         * @details     first call to this function returns the partitioned range.  subsequent calls get the end object (start == end)
         *              to call again and get meaningful result, reset() must be called.
         *              no need to really store much - everything can be recomputed relatively fast.
         * @param partId   partition id for the sub range.
         * @return      range of the partition
         */
        template<typename R = Range>
        typename std::enable_if<std::is_floating_point<typename R::ValueType>::value, Range>::type& getNextImpl(const size_t& partId) {
          // param validation already done.

          // if previously computed, then return end object, signalling no more chunks.
          if (done) return this->end;  // can only call this once.

          // if just 1 partition, return.
          if (this->nPartitions == 1) return this->src;

          // compute the subrange's start and end, spreading out the remainder to the first rem chunks/partitions.
          BaseClassType::computeRangeForChunkId(curr, this->src, 0, partId, this->chunkSize, this->overlapSize);

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
         * @typedef RangeValueType
         */
        using RangeValueType = typename BaseClassType::RangeValueType;

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

				/**
				 * @var BEFORE
         * @brief	static variable indicating we have not yet called "getNext"
			   */
        static const uint8_t BEFORE = 0;
				
				/**
				 * @var DURING
         * @brief	static variable indicating we have called getNext at least once, but are not at the end yet
			   */
        static const uint8_t DURING = 1;
        
        /**
				 * @var BEFORE
         * @brief	static variable indicating we have reached the end of the range during calls to getNext
			   */
        static const uint8_t AFTER = 2;

      public:

        CyclicPartitioner() : BaseClassType(), state(nullptr), curr(nullptr), nChunks(0), stride(0) {};

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
        void configure(const Range &_src, const size_t &_nPartitions, const ChunkSizeType &_chunkSize, const RangeValueType &_overlapSize = 0) {
          this->BaseClassType::configure(_src, _nPartitions, _chunkSize, _overlapSize);

          // get the maximum number of chunks in the source range.
          // TODO: make this compatible with floating point.
          nChunks = BaseClassType::computeNumberOfChunks();

          // stride:  if there are less chunks than partition, we can only walk through once, so stride is src size.
          stride = (nChunks > this->nPartitions) ? 
          	this->chunkSize * static_cast<ChunkSizeType>(this->nPartitions) : 
          	this->src.size();

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

          /// if nChunks < nPartitions, then each partId will only get a chunk once.
          // in that case, if partId is >= nChunks, then end is returned.
          if (partId >= nChunks) return this->end;

          // if this partition is done, return the last entry (end)
          if (state[partId] == AFTER)  return this->end;

          /// first iteration, use initialized value
          if (state[partId] == BEFORE) {
            state[partId] = DURING;
            return curr[partId];
          }
          // else not the first and not last, so increment.

         /// comparing to amount of available range and set the start, end, overlap - trying to avoid data type overflow.

          if (this->src.end - curr[partId].end > stride) {
            // has room.  so shift  starting and end position by stride length
            curr[partId].start += stride;
            curr[partId].end += stride;
            // overlap doesn't need to change
          } else if (this->src.end - curr[partId].start <= stride) {
            // adding stride would make subrange start after parent range's end.
            state[partId] = AFTER;
            return this->end;
          } else {
            // parent range's end falls in the overlap (overlap) region
            // start can be moved.
            curr[partId].start += stride;
            // end is definitely at parent ranges' end
            curr[partId].end = this->src.end;
            // overlap needs to be computed
       	    curr[partId].overlap = (this->src.end - curr[partId].start <= this->chunkSize) ?
					    0 : this->src.end - curr[partId].start - this->chunkSize;
            // may have more for this partition id - if nPartition is 1, so not at end yet.
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
            state[i] = BEFORE;
            BaseClassType::computeRangeForChunkId(curr[i], this->src, 0, i, this->chunkSize, this->overlapSize);
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
         * @typedef RangeValueType
         * @brief   type for the start/end/overlap
         */
        using RangeValueType = typename BaseClassType::RangeValueType;

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
         * @brief     internal method to atomically increment the chunkOffset.  this is the version for integral type
         * @return    old value.  chunkOffset incremented.
         */
        template<typename T = ChunkSizeType>
        typename std::enable_if<std::is_integral<T>::value, RangeValueType>::type getNextOffset() {
          return chunkOffset.fetch_add(this->chunkSize, std::memory_order_acq_rel);
        }
        /**
         * @brief     internal method to atomically increment the chunkOffset for floating point types
         * @details   std::atomics does not support fetch_add on non-integral types.
         * @return    old value.  chunkOffset incremented.
         */
        template<typename T = ChunkSizeType>
        typename std::enable_if<std::is_floating_point<T>::value, RangeValueType>::type getNextOffset() {
          RangeValueType origval = chunkOffset.load(std::memory_order_consume);
          RangeValueType newval;
          do {
            newval = origval + this->chunkSize;
          } while (!chunkOffset.compare_exchange_weak(origval, newval, std::memory_order_acq_rel, std::memory_order_acquire));
          return origval;
        }

      public:
        DemandDrivenPartitioner() : BaseClassType(), chunkOffset(0), chunkId(0), nChunks(0), curr(nullptr), done(false) {};

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
        void configure(const Range &_src, const size_t &_nPartitions, const ChunkSizeType &_chunkSize, const RangeValueType &_overlapSize = 0) {
          this->BaseClassType::configure(_src, _nPartitions, _chunkSize, _overlapSize);

          // get the maximum number of chunks in the source range.
          // TODO: make this compatible with floating point.
          nChunks = BaseClassType::computeNumberOfChunks();

          // if there are more partitions than chunks, then the array represents mapping from chunkId to subrange
          // else if there are more chunks than partitions, then the array represents the most recent chunk assigned to a partition.
          curr = new Range[std::min(nChunks, this->nPartitions)];

          resetImpl();
        };




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

          // call internal function (so integral and floating point types are handled properly)
          RangeValueType s = getNextOffset<ChunkSizeType>();

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

            curr[id].start = s;
                // use comparison to avoid overflow.
            if (this->src.end - s > this->chunkSize + this->overlapSize) {
              // enough room for overlapSize
              curr[id].end = s + this->chunkSize + this->overlapSize;
              // range's overlap is left as default
            } else {
              // not enough room for even the chunk
              curr[id].end = this->src.end;
              if (this->src.end - s <= this->chunkSize) {
                curr[id].overlap = 0;
              } else {
                curr[id].overlap = this->src.end - s - this->chunkSize;
              }
            }

            return curr[id];
          }
        }

        /**
         * @brief resets the partitioner by resetting the internal subrange arrays.  also reset the offset to the start of the src range, chunk Id, and "done".
         * @details this function also serves to initialize the subrange array.          chunkOffset.store(this->src.start, std::memory_order_release);
          chunkId.store(0, )
         *
         */
        void resetImpl() {
          // these 3 calls probably should be synchronized together.
          chunkOffset.store(this->src.start, std::memory_order_release);
          chunkId.store(0, std::memory_order_release);
          done.store(false, std::memory_order_release);

          // range will be completely recomputed, so don't have to do much here.
        }


    };


  } /* namespace partition */
} /* namespace bliss */

#endif /* PARTITIONER_HPP_ */
