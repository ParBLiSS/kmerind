/**
 * @file		Partitioner.hpp
 * @ingroup
 * @author	tpan
 * @brief
 * @details
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef PARTITIONER_HPP_
#define PARTITIONER_HPP_

#include <cassert>
#include <config.hpp>
#include <partition/range.hpp>

namespace bliss
{
  namespace partition
  {

    /**
     * @class			bliss::partition::Partitioner
     * @brief
     * @details   partition a range.  (without knowledge of the data specified by the range,  so no "search" type of partitioning)
     *
     *            Simple types (essentially functinoids, use default copy constructor and assignment operator.  no need to move versions.
     *
     */
    template<typename Range, typename Derived>
    class Partitioner {
      protected:
        typedef typename Range::ValueType T;

        Range r;
        Range end;
        T length;
        int nParts;
        size_t chunkSize;

      public:
        virtual ~Partitioner() {};

        void configure(const Range &src, const int &numPartitions, const size_t &_chunkSize) {

          assert(numPartitions > 0);
          assert(_chunkSize > 0);

          r = src;
          end.start = end.end = src.end;
          length = src.end - src.start;
          nParts = numPartitions;
          chunkSize = _chunkSize;
        }

        inline Range& getNext(const int pid) {
//          if (pid >= nParts) printf("pid: %d, parts %d", pid, nParts);

          assert(pid >= 0);
          assert(pid < nParts);  // both non-negative.

          return static_cast<Derived*>(this)->getNextImpl(pid);
        }

        void reset() {
          static_cast<Derived*>(this)->resetImpl();
        }

    };


    template<typename Range>
    class BlockPartitioner : public Partitioner<Range, BlockPartitioner<Range> >
    {
      protected:
        Range curr;
        bool done;

        typedef typename Range::ValueType T;

        T div;
        T rem;

      public:
        virtual ~BlockPartitioner() {};

        void configure(const Range &src, const int &numPartitions, const size_t &_chunkSize) {
          this->Partitioner<Range, BlockPartitioner<Range> >::configure(src, numPartitions, _chunkSize);

          div = this->length / numPartitions;
          rem = this->length % numPartitions;

          resetImpl();
        }


        inline Range& getNextImpl(const int pid) {

          if (done) return this->end;  // can only call this once.

          if (this->nParts == 1) return this->r;

          if (pid < rem)
          {
            curr.start += static_cast<T>(pid * (div + 1));
            curr.end = curr.start + static_cast<T>(div + 1) + this->r.overlap;
          }
          else
          {
            curr.start += static_cast<T>(pid * div + rem);
            curr.end = curr.start + static_cast<T>(div) + this->r.overlap;
          }

          // last entry.  no overlap.
          curr &= this->r;

          curr.block_start = curr.start;

          done = true;
          return curr;
        }

        void resetImpl() {
          curr = this->r;
          done = false;
        }

    };


    template<typename Range>
    class CyclicPartitioner : public Partitioner<Range, CyclicPartitioner<Range> >
    {

      protected:
        bool *done;
        Range *curr;
        size_t iterSize;

      public:
        CyclicPartitioner() : done(nullptr), curr(nullptr), iterSize(0) {};

        virtual ~CyclicPartitioner() {
          if (done) delete [] done;
          if (curr) delete [] curr;
        }

        void configure(const Range &src, const int &numPartitions, const size_t &_chunkSize) {
          this->Partitioner<Range, CyclicPartitioner<Range> >::configure(src, numPartitions, _chunkSize);

          done = new bool[numPartitions];
          curr = new Range[numPartitions];

          resetImpl();

          iterSize = _chunkSize * numPartitions;
        }

        inline Range& getNextImpl(const int pid) {

          if (done[pid])  return this->end;                                // if done, return.

          curr[pid].start += iterSize;  // compute starting positio
          if (curr[pid].start > this->r.end) {
            done[pid] = true;         // if outside range, done
            return this->end;
          }

          curr[pid].end = curr[pid].start + this->chunkSize;
          curr[pid] &= this->r;  // keep within range

          return curr[pid];
        }

        void resetImpl()
        {
          for (size_t i = 0; i < this->nParts; ++i) {
            done[i] = false;
            curr[i] = Range(i * this->chunkSize, (i+1) * this->chunkSize, this->r.overlap) & this->r;
          }
        }


    };



    template<typename Range>
    class DemandDrivenPartitioner : public Partitioner<Range, DemandDrivenPartitioner<Range> >
    {

      protected:
        size_t chunkOffset;
        Range *curr;
        bool done;

        typedef typename Range::ValueType T;

      public:
        DemandDrivenPartitioner() : chunkOffset(0), curr(nullptr), done(false) {};

        virtual ~DemandDrivenPartitioner() {
          if (curr != nullptr) delete [] curr;
        }
        void configure(const Range &src, const int &numPartitions, const size_t &_chunkSize) {
          this->Partitioner<Range, DemandDrivenPartitioner<Range> >::configure(src, numPartitions, _chunkSize);

          curr = new Range[numPartitions];

          resetImpl();
        };

        inline Range& getNextImpl(const int pid) {
          // TODO:  what is MPI equivalent?
          if (done) return this->end;

          T s = 0;

#pragma omp atomic capture
          {
            s = chunkOffset;
            chunkOffset += this->chunkSize;
          }

          if (s >= this->r.end) {
            done = true;
            return this->end;
          }

          curr[pid].start = s;
          curr[pid].end = s + this->chunkSize;
          curr[pid] &= this->r;

          return curr[pid];
        }

        void resetImpl() {
          chunkOffset = this->r.start;
          done = false;
          for (size_t i = 0; i < this->nParts; ++i) {
            curr[i] = Range();
          }
        }


    };


  } /* namespace partition */
} /* namespace bliss */

#endif /* PARTITIONER_HPP_ */
