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
     */
    template<typename T, typename Derived>
    class Partitioner {
      public:
        Partitioner(const range<T> &src, const int &numPartitions, const size_t &_chunkSize) : r(src), end(src.end, src.end), nParts(numPartitions), chunkSize(_chunkSize) {
          length = src.end - src.start;
        };
        virtual ~Partitioner() {};

        inline range<T>& getNext(const int pid) {
          assert(pid >= 0);
          assert(pid < nParts);  // both non-negative.

          return static_cast<Derived*>(this)->getNextImpl(pid);
        }

        void reset() {
          static_cast<Derived*>(this)->resetImpl();
        }

      protected:
        range<T> r;
        range<T> end;
        T length;
        int nParts;
        size_t chunkSize;

    };


    template<typename T>
    class BlockPartitioner : public Partitioner<T, BlockPartitioner<T> >
    {
      public:
        BlockPartitioner(const range<T> &src, const int &numPartitions, const size_t &chunkSize) :
          Partitioner<T, BlockPartitioner<T> >(src, numPartitions), curr(src), done(false) {

          div = length / nParts;
          rem = length % nParts;
        }
        virtual ~BlockPartitioner() {};

      protected:
        range<T> curr;
        bool done;
        T div;
        T rem;

        inline range<T>& getNextImpl(const int pid) {

          if (done) return end;  // can only call this once.

          if (nParts == 1) return r;

          if (pid < rem)
          {
            curr.start += static_cast<T>(pid * (div + 1));
            curr.end = curr.start + static_cast<T>(div + 1) + r.overlap;
          }
          else
          {
            curr.start += static_cast<T>(pid * div + rem);
            curr.end = curr.start + static_cast<T>(div) + r.overlap;
          }

          // last entry.  no overlap.
          curr &= r;

          curr.block_start = curr.start;

          done = true;
          return curr;
        }

        void resetImpl() {
          curr = r;
          done = false;
        }
    };


    template<typename T>
    class CyclicPartitioner : public Partitioner<T, CyclicPartitioner<T> >
    {
      public:
        CyclicPartitioner(const range<T> &src, const int &numPartitions, const size_t &_chunkSize) :
          Partitioner<T, CyclicPartitioner<T> >(src, numPartitions) {

          done = new bool[numPartitions];
          curr = new range<T>[numPartitions];

          resetImpl();

          iterSize = _chunkSize * numPartitions;
        };
        virtual ~CyclicPartitioner() {
          delete [] done;
          delete [] curr;
        }

      protected:
        bool *done;
        range<T> *curr;
        size_t iterSize;

        inline range<T>& getNextImpl(const int pid) {

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

        void resetImpl()
        {
          for (size_t i = 0; i < numPartitions; ++i) {
            done[i] = false;
            curr[i] = range<T>(i * chunkSize, (i+1) * chunkSize, r.overlap) & r;
          }
        }
    };



    template<typename T>
    class DemandDrivenPartitioner : public Partitioner<T, DemandDrivenPartitioner<T> >
    {
      public:
        DemandDrivenPartitioner(const range<T> &src, const int &numPartitions, const size_t &_chunkSize) :
          Partitioner<T, DemandDrivenPartitioner<T> >(src, numPartitions) {

          curr = new range<T>[numPartitions];

          resetImpl();
        };
        virtual ~DemandDrivenPartitioner() {
          delete [] curr;
        };

      protected:
        size_t chunkOffset;
        range<T> *curr;
        bool done;

        inline range<T>& getNextImpl(const int pid) {
          // TODO:  what is MPI equivalent?
          if (done) return end;

          T s = 0;

#pragma omp atomic capture
          {
            s = chunkOffset;
            chunkOffset += chunkSize;
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

        void resetImpl() {
          chunkOffset = r.start;
          done = false;
          for (size_t i = 0; i < numPartitions; ++i) {
            curr[i] = range<T>();
          }
        }
    };


  } /* namespace partition */
} /* namespace bliss */

#endif /* PARTITIONER_HPP_ */
