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
     * @details
     *
     */
    template<typename T, typename Derived>
    class Partitioner {
      public:
        Partitioner(const range<T> &src, const size_t &numPartitions) : r(src), end(src.end, src.end), nParts(numPartitions) {
          assert( r.start <= r.end);
          assert( r.overlap >= 0);
        };
        virtual ~Partitioner() {};

        inline range<T>& getNext(const size_t &pid) {
          assert(pid >= 0);
          assert(pid < nParts);  // both non-negative.

          return static_cast<Derived*>(this)->getNextImpl(pid);
        }

      protected:
        range<T> r;
        range<T> end;
        size_t nParts;

    };


    template<typename T>
    class BlockPartitioner : public Partitioner<T, BlockPartitioner<T> >
    {
      public:
        BlockPartitioner(const range<T> &src, const size_t &numPartitions) :
          Partitioner<T, BlockPartitioner<T> >(src, numPartitions), curr(src), done(false) {

        }
        virtual ~BlockPartitioner() {};

      protected:
        range<T> curr;
        bool done;

        inline range<T>& getNextImpl(const size_t &pid) {

          if (done) return end;  // can only call this once.

          if (nParts == 1) return r;

          T length = curr.end - curr.start;

          T div = length / nParts;
          T rem = length % nParts;
          if (pid < rem)
          {
            curr.start += static_cast<T>(pid * (div + 1));
            curr.end = curr.start + static_cast<T>(div + 1) + curr.overlap;
          }
          else
          {
            curr.start += static_cast<T>(pid * div + rem);
            curr.end = curr.start + static_cast<T>(div) + curr.overlap;
          }

          // last entry.  no overlap.
          curr &= r;

          curr.block_start = curr.start;

          done = true;
          return curr;
        }
    };


    template<typename T>
    class CyclicPartitioner : public Partitioner<T, CyclicPartitioner<T> >
    {
      public:
        CyclicPartitioner(const range<T> &src, const size_t &numPartitions, const size_t &_chunkSize) :
          Partitioner<T, CyclicPartitioner<T> >(src, numPartitions), chunkSize(_chunkSize) {

          done = new bool[numPartitions];
          curr = new range<T>[numPartitions];
          for (size_t i = 0; i < numPartitions; ++i) {
            done[i] = false;
            curr[i] = range<T>(i * _chunkSize, (i+1) * _chunkSize, src.overlap) & r;
          }
          iterSize = _chunkSize * numPartitions;
        };
        virtual ~CyclicPartitioner() {
          delete [] done;
          delete [] curr;
        }

      protected:
        bool *done;
        range<T> *curr;
        size_t chunkSize;
        size_t iterSize;

        inline range<T>& getNextImpl(const size_t &pid) {

          if (done[pid])  return end;                                // if done, return.

          curr[pid].start += iterSize;  // compute starting positio
          if (curr[pid].start > r.end) {
            done[pid] = true;         // if outside range, done
            return end;
          }

          curr[pid].end = curr[pid].start + chunkSize;
          curr[pid] &= r;  // keep within range

          return curr[pid];
        }
    };



    template<typename T>
    class DemandDrivenPartitioner : public Partitioner<T, DemandDrivenPartitioner<T> >
    {
      public:
        DemandDrivenPartitioner(const range<T> &src, const size_t &numPartitions, const size_t &_chunkSize) :
          Partitioner<T, DemandDrivenPartitioner<T> >(src, numPartitions), chunkOffset(src.start), chunkSize(_chunkSize), done(false) {

          curr = new range<T>[numPartitions];
          for (size_t i = 0; i < numPartitions; ++i) {
            curr[i] = range<T>();
          }
        };
        virtual ~DemandDrivenPartitioner() {
          delete [] curr;
        };

      protected:
        size_t chunkOffset;
        size_t chunkSize;
        range<T> *curr;
        bool done;

        inline range<T>& getNextImpl(const size_t &pid) {
          // TODO:  what is MPI equivalent?
          if (done) return end;

          T s = 0;

#pragma omp atomic capture
          {
            s = chunkOffset;
            chunkPos += chunkSize;
          }

          if (s >= r.end) {
            done = true;
            return end;
          }

          curr[pid].start = s;
          curr[pid].end = s + chunkSize;
          curr[pid] &= r;

          return curr[pid];
        }
    };


  } /* namespace partition */
} /* namespace bliss */

#endif /* PARTITIONER_HPP_ */
