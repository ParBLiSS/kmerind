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
 * @file		partitioner.hpp
 * @ingroup partition
 * @author	Tony Pan <tpan7@gatech.edu>
 * @brief   contains several class that provide different logic for partitioning a range
 * @details contains block, cyclic, and demand driven (THREAD SAFE) partitioners
						logic implementation uses comparison to avoid overflows and implicit casting where needed.

 */
#ifndef PARTITIONER_HPP_
#define PARTITIONER_HPP_

#include <stdexcept>
#include <atomic>
#include "bliss-config.hpp"
#include <type_traits>
#include "partition/range.hpp"

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
     *
     *            uses default constructor, copy constructor, and move constructor.
     *
     *            Works for continuous value ranges.
     *
     *            Partitioners are designed for multithreaded usage.  For Block and Cyclic partitioners, the chunks are computed
     *            independent of other threads.  For DemandDriven partitioner, the chunks are computed using atomc operation to provide
     *            memory fence/synchronization.
     *
     * @note      The purpose of partitioning is to divide the range of data for computation.  Partitioner will therefore divide the source range into parts.
     *            An overlap region length can be specified during partitioning.  When overlap region length is specified, all subranges
     *            from the partitioner will have its start and previous subrange's end overlap.
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
         * @typedef SizeType
         * @brief   value type used for size.  if integral type for valueType, then size_t.  else (floating point) same as RangeValueType
         */
        using SizeType = typename std::conditional<std::is_floating_point<RangeValueType>::value,
                                                    RangeValueType, size_t>::type;

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
         * @var nChunks
         * @brief    number of chunks present in the src range.
         */
        size_t nChunks;

        /**
         * @var non_overlap_size
         * @brief   size of each partition excluding the overlap region.  computed for block partitioner, and user specified for cyclic and demand driven.
         */
        SizeType non_overlap_size;

        /**
         * @var     overlap_size
         * @brief   the overlap_size that each subrange should contain.  Parameter for this partitioning.
         */
        SizeType overlap_size;

        /**
         * @brief   computes the number of chunks in a src range.  for Integral type only.
         * @details computed using the nonoverlapping part of each subrange, length in non_overlap_size
         * @tparam  R  range type.  used for enable_if.
         * @return  number of chunks in a range.
         */
        template <typename R = Range>
        inline typename std::enable_if<!std::is_floating_point<typename R::ValueType>::value, size_t>::type computeNumberOfChunks() {
          return static_cast<std::size_t>( (this->src.size() + this->non_overlap_size - 1) / this->non_overlap_size);
        }
        /**
         * @brief   computes the number of chunks in a src range, excluding the overlap region.  for floating type only.
         * @tparam  R  range type.  used for enable_if.
         * @return  number of chunks in a range, based on non_overlap_size.
         */
        template <typename R = Range>
        inline typename std::enable_if<std::is_floating_point<typename R::ValueType>::value, size_t>::type computeNumberOfChunks() {
          return static_cast<std::size_t>(std::ceil(this->src.size() / this->non_overlap_size));
        }

        /**
         * @brief   internal function to compute the range for a chunk of the parent range.
         * @details The chunk range is computed by first setting the start to parent range start + start offset (e.g. rem from BlockPartitioner)
         * @note    chunkId must not exceed maximum number of chunks - enforced outside of this function.
         *
         * @param[out] r            the range to compute.  initialized in the function to pr.start + start_offset
         * @param[in] pr            the parent range
         * @param[in] start_offset  offset to add to the start position.  only used for block partitioner to spread the remainders to first blocks.
         * @param[in] chunkId       the chunkId of the current range.
         * @param[in] non_overlap_size    the size of the chunk excluding the overlapped portion        (precomputed or user specified)
         */
        Range computeRangeForChunkId(const Range & pr, const SizeType& start_offset, const size_t& chunkId) {

          Range r(pr.start + start_offset, pr.end);
          // compute start
          if (static_cast<SizeType>(pr.end - r.start) > static_cast<SizeType>(chunkId) * this->non_overlap_size) {
            r.start += static_cast<SizeType>(chunkId) * this->non_overlap_size;
          } else {
            r.start = pr.end;
          }

          // compute end
          if (static_cast<SizeType>(pr.end - r.start) > (this->non_overlap_size + this->overlap_size)) {
            // enough room
            r.end = r.start + this->non_overlap_size + this->overlap_size;
          } else {
            r.end = pr.end;
          }
          return r;
        }

      public:
        Partitioner() : src(), end(), nPartitions(1), nChunks(1), non_overlap_size(0), overlap_size(0) {};

        /**
         * @brief default destructor
         */
        virtual ~Partitioner() {};

        /**
         * @brief configures the partitioner with the source range, number of partitions.
         * @param _src          range object to be partitioned.
         * @param _nPartitions the number of partitions to divide this range into
         * @param _non_overlap_size   the size of each partition.  default is 0 (in subclasses), computed later.
         * @param _overlap_size   the size of the overlap region.
         * @return updated non_overlap_size.
         */
        SizeType configure(const Range &_src, const size_t &_nPartitions, const SizeType &_non_overlap_size, const SizeType &_overlap_size) {

//          if (_src.size() == 0)
//            throw std::invalid_argument("ERROR: partitioner c'tor: src range size is 0");

          if (_nPartitions == 0)
            throw std::invalid_argument("ERROR: partitioner c'tor: nPartitions is 0");
          this->nPartitions = _nPartitions;

          if (_overlap_size < 0)
            throw std::invalid_argument("ERROR: partitioner c'tor: overlap_size is < 0");
          this->overlap_size = _overlap_size;

//          if (_non_overlap_size <= 0)
//            throw std::invalid_argument("ERROR: partitioner c'tor: non_overlap_size is <= 0");
          // check elsewhere.
          this->non_overlap_size = _non_overlap_size;

          this->src = _src;  // leaving the overlap region size values same as before for src.
          this->end = Range(this->src.end, this->src.end);


//          std::cout << "partitioner configured for range " << this->src << " with " << nPartitions <<
//              " parts, size " << non_overlap_size << " + " << overlap_size << std::endl;

          return this->non_overlap_size;
        }

        /**
         * @brief get the next sub range (chunk) within a partition.
         * @details the function calls the specific subclass implementation.
         * @param partId    the id of the partition to retrieve
         * @return          range object containing the start and end of the next chunk in the partition.
         */
        inline Range getNext(const size_t &partId) {
          // both non-negative.
          if (partId >= this->nPartitions)
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
        friend Partitioner<Range, BlockPartitioner<Range> >;


      protected:
        /**
         * @typedef BaseClassType
         * @brief   the superclass type.
         */
        using BaseClassType = Partitioner<Range, BlockPartitioner<Range> >;

        /**
         * @typedef SizeType
         */
        using SizeType = typename BaseClassType::SizeType;

        /**
         * @typedef ValueType
         */
        using RangeValueType = typename BaseClassType::RangeValueType;

        /**
         * @var done
         * @brief boolean indicating that partition has already been done. THREAD SAFE
         */
        bool* done;

        /**
         * @var rem
         * @brief the left-over that is to be spread out to the first rem partitions.
         */
        SizeType rem;


      public:

        BlockPartitioner() : BaseClassType(), done(nullptr), rem(0) {};

        /**
         * @brief default destructor
         */
        virtual ~BlockPartitioner() {
          if (done != nullptr) {
            delete [] done;
            done = nullptr;
          }
        };

        /**
         * @brief configures the partitioner with the source range, number of partitions
         * @param _src          range object to be partitioned.
         * @param _nPartitions the number of partitions to divide this range into
         * @param _non_overlap_size  size of each chunk.  here it should be 0 so they will be auto computed.
         * @param _overlap_size   the size of the overlap region.
         * @return updated non_overlap_size.
         *
         */
        SizeType configure(const Range &_src, const size_t &_nPartitions,
                       const SizeType &_non_overlap_size = 1, const SizeType &_overlap_size = 0) {
          if (_nPartitions == 0)
            throw std::invalid_argument("ERROR: partitioner: number of partitions is 0");


          SizeType chunk_size = (_src.size() / static_cast<SizeType>(_nPartitions));

          // okay for non_overlap_size to be 0 here.
          this->BaseClassType::configure(_src, _nPartitions, chunk_size, _overlap_size);

          //Array of booleans for each partition
          if (done) delete [] done;
          done = new bool[this->nPartitions];

          // next figure out the remainder using non_overlap_size
          // not using modulus since we SizeType may be float.
          this->rem = _src.size() - this->non_overlap_size * static_cast<SizeType>(this->nPartitions);

          this->nChunks = (chunk_size == 0 && !std::is_floating_point<SizeType>::value) ? this->rem : _nPartitions;

          resetImpl();

          return this->non_overlap_size;
        }

      protected:

        /**
         * @brief       get the next chunk in the partition.  for BlockPartition, there is only 1 chunk.  INTEGRAL Type Only
         * @details     first call to this function returns the partitioned range.  subsequent calls get the end object (start == end)
         *              to call again and get meaningful result, reset() must be called.
         *              no need to really store much - everything can be recomputed relatively fast.
         * @param partId   partition id for the sub range.
         * @return      range of the partition
         */
        template<typename R = Range>
        typename std::enable_if<!std::is_floating_point<typename R::ValueType>::value, Range>::type
        getNextImpl(const size_t& partId) {
          // param validation already done.

          // if previously computed, then return end object, signaling no more chunks.
          if (done[partId]) return this->end;  // can only call this once.

          done[partId] = true;

          // if just 1 partition, return.
          if (this->nPartitions == 1) {

            return this->src;
          } else {

            // compute the subrange's start and end, spreading out the remainder to the first rem chunks/partitions.
            // each chunk with partition id < rem gets 1 extra.  after that offset is rem.
            SizeType offset = (partId < this->rem) ? partId : this->rem;

            // compute the new range
            Range r = BaseClassType::computeRangeForChunkId(this->src, offset, partId);
            if (partId < this->rem) r.end += 1;  // each chunk with partition id < rem gets 1 extra.

//            std::cout << "computed range: " << r << std::endl;

            return r;
          }
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
        typename std::enable_if<std::is_floating_point<typename R::ValueType>::value, Range>::type
        getNextImpl(const size_t& partId) {
          // param validation already done.

          // if previously computed, then return end object, signalling no more chunks.
          if (done[partId]) return this->end;  // can only call this once.

          done[partId] = true;

          // if just 1 partition, return.
          if (this->nPartitions == 1) {
            return this->src;
          } else {
            // compute the subrange's start and end, spreading out the remainder to the first rem chunks/partitions.
            return BaseClassType::computeRangeForChunkId(this->src, 0.0, partId);
          }
        }


        /**
         * @brief reset size to full range, mark the done to false.
         */
        void resetImpl() {
          if (done)
            for(size_t i=0 ; i < this->nPartitions; i++)
            {
              done[i] = false;
            }
        }

    };

    /**
     * @class CyclicPartitioner
     * @brief cyclically partition the range into non_overlap_size blocks.
     *        uses CRTP pattern to enforce consistent API.
     * @tparam Range  type of range to be partitioned.
     */
    template<typename Range>
    class CyclicPartitioner : public Partitioner<Range, CyclicPartitioner<Range> >
    {
        friend Partitioner<Range, CyclicPartitioner<Range> >;


      protected:
        /**
         * @typedef BaseClassType
         * @brief   the superclass type.
         */
        using BaseClassType = Partitioner<Range, CyclicPartitioner<Range> >;

        /**
         * @typedef SizeType
         */
        using SizeType = typename BaseClassType::SizeType;

        /**
         * @typedef RangeValueType
         */
        using RangeValueType = typename BaseClassType::RangeValueType;

        /**
         * @var done
         * @brief An array of "state", one for each partition.
         */
        size_t *n_strides;


      public:

        CyclicPartitioner() : BaseClassType(), n_strides(nullptr) {};

        /**
         * @brief default destructor.  cleans up arrays for callers.
         */
        virtual ~CyclicPartitioner() {
          if (n_strides) {
            delete [] n_strides;
            n_strides = nullptr;
          }
        }

        /**
         * @brief configures the partitioner with the source range, number of partitions
         * @param _src          range object to be partitioned.
         * @param _nPartitions the number of partitions to divide this range into
         * @param _non_overlap_size    Size of the chunks
         * @param _overlap_size   the size of the overlap region.
         * @return updated non_overlap_size.
         *
         */
        SizeType configure(const Range &_src, const size_t &_nPartitions, const SizeType &_non_overlap_size, const SizeType &_overlap_size = 0) {
          if (_non_overlap_size <= 0)
            throw std::invalid_argument("ERROR: partitioner c'tor: non_overlap_size is <= 0");

          this->BaseClassType::configure(_src, _nPartitions, _non_overlap_size, _overlap_size);

          if (n_strides) delete [] n_strides;
          n_strides = new size_t[this->nPartitions];

          this->nChunks = this->computeNumberOfChunks();

          resetImpl();

          return this->non_overlap_size;
        }

      protected:
        /**
         * @brief       get the next chunk in the partition.  for CyclickPartition, keep getting until done.
         * @details     each call to this function gets the next chunk belonging to this partition.
         *              for cyclicPartitioner, the chunks are separated by non_overlap_size * nPartitions
         * @note        if nChunks < nPartitions, then each partId will only get a chunk once.
         *                in that case, if partId is >= nChunks, then end is returned.
         * @param partId   partition id for the sub range.
         * @return      range of the partition.  if no more, return "end", which has start and end position set to end.
         */
        inline Range getNextImpl(const size_t& partId) {

          // convert partId to the current chunk_id based on n_strides (number of iterations);
          size_t chunk_id = partId + this->nPartitions * n_strides[partId];

          // if chunk_id is greater than total number of nChunks, return empty range.
          if (chunk_id >= this->nChunks) return this->end;

          // else can return a real chunk.
          else {
            n_strides[partId] += 1;
            return BaseClassType::computeRangeForChunkId(this->src, 0, chunk_id);

          }

        }

        /**
         * @brief resets the partitioner by resetting the internal done and subrange arrays.
         * @details this function also serves to initialize the arrays.
         */
        void resetImpl()
        {
          if (n_strides)
            for (size_t i = 0; i < this->nPartitions; ++i) {
              n_strides[i] = 0;
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

        friend Partitioner<Range, DemandDrivenPartitioner<Range> >;

      protected:
        /**
         * @typedef BaseClassType
         * @brief   the superclass type.
         */
        using BaseClassType = Partitioner<Range, DemandDrivenPartitioner<Range> >;

        /**
         * @typedef SizeType
         */
        using SizeType = typename BaseClassType::SizeType;

        /**
         * @typedef RangeValueType
         * @brief   type for the start/end/overlap
         */
        using RangeValueType = typename BaseClassType::RangeValueType;

        /**
         * @var chunk_id
         * @brief  the id for the next chunk to be returned.
         */
        std::atomic<size_t> chunk_id;


        /**
         * @done
         * @brief boolean indicating that the range has been exhausted by the traversal.  THREAD SAFE
         * @details this variable is here because if we only check chunk_id, it might get deleted by one thread while another is incrementing it, resulting in a value that is still in a valid range.
         *    using an atomic bool is less likely for this to happen.
         */
        std::atomic<bool> done;


      public:
        DemandDrivenPartitioner() : BaseClassType(), chunk_id(0), done(false) {};

        /**
         * @brief default destructor
         */
        virtual ~DemandDrivenPartitioner() {}


        /**
         * @brief configures the partitioner with the source range, number of partitions
         * @note  should be called by single thread.
         * @param _src          range object to be partitioned.
         * @param _nPartitions  the number of partitions to divide this range into
         * @param _non_overlap_size  size of each chunk for the partitioning.
         * @return updated non_overlap_size.
         *
         */
        SizeType configure(const Range &_src, const size_t &_nPartitions, const SizeType &_non_overlap_size, const SizeType &_overlap_size = 0) {
          if (_non_overlap_size <= 0)
            throw std::invalid_argument("ERROR: partitioner c'tor: non_overlap_size is <= 0");

          this->BaseClassType::configure(_src, _nPartitions, _non_overlap_size, _overlap_size);

          this->nChunks = this->computeNumberOfChunks();

          resetImpl();

          return this->non_overlap_size;
        };


      protected:

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
         inline Range getNextImpl(const size_t& partId) {

          // all done, so return end
          if (done.load(std::memory_order_seq_cst)) return this->end;

          // convert partId to the current chunk_id based on n_strides (number of iterations);
          size_t chunk_id = this->chunk_id.fetch_add(1, std::memory_order_seq_cst);

          // if chunk_id is greater than total number of nChunks, return empty range.
          if (chunk_id >= this->nChunks) {
            done.store(true, std::memory_order_seq_cst);
            return this->end;
          }

          // else can return a real chunk.
          else
            return BaseClassType::computeRangeForChunkId(this->src, 0, chunk_id);

        }

        /**
         * @brief resets the partitioner by resetting the internal subrange arrays.  also reset the offset to the start of the src range, chunk Id, and "done".
         * @details this function also serves to initialize the subrange array.          chunkOffset.store(this->src.start, std::memory_order_seq_cst);
          chunkId.store(0, )
         *
         */
        void resetImpl() {
          chunk_id.store(0, std::memory_order_seq_cst);
          done.store(false, std::memory_order_seq_cst);

        }


    };


  } /* namespace partition */
} /* namespace bliss */

#endif /* PARTITIONER_HPP_ */
