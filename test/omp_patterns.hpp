/**
 * @file		omp_patterns.hpp
 * @ingroup
 * @author	Tony Pan <tpan7@gatech.edu>
 * @brief
 * @details
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef OMP_PATTERNS_HPP_
#define OMP_PATTERNS_HPP_

#ifdef USE_OPENMP
#include <omp.h>
#endif

#include "partition/range.hpp"
typedef bliss::partition::range<size_t> RangeType;


////// TODO:  need to make the pipelines better.  DO THESE MAKE SENSE?
///**
// * transforms data retrieved from InputOp, using ComputeOp, and put results in OutputOp.
// * InputOp and OutputOp are operators, not just data structures.  examples are MPI receiver/senders
// *
// * Map or Reduce pattern.
// *
// * @param in
// * @param comp
// * @param out
// * @param nthreads
// * @param r
// * @param step
// */
//template<typename InputOp, typename ComputeOp, typename OutputOp>
//void TransformingPipeline(InputOp &in, ComputeOp &comp, OutputOp &out, const int &nthreads, const RangeType &r, const size_t step) {
//#pragma omp parallel sections num_threads(3) OMP_SHARE_DEFAULT
//  {
//#pragma omp section
//    {
//      out(nthreads, r, step);
//    }
//
//#pragma omp section
//    {
//      in(nthreads, r, step);
//    }
//
//#pragma omp section
//    {
//      comp(nthreads, r, step);
//    }
//  }
//}
//
//
///**
// * generate data using ComputeOp and put results in OutputOp
// * OutputOp is an operator, not jsut data structure.  example is MPI Sender
// *
// * Map or Reduce pattern.
// *
// * @param comp
// * @param out
// * @param nthreads
// * @param r
// * @param step
// */
//template<typename ComputeOp, typename OutputOp>
//void GeneratingPipeline(ComputeOp &comp, OutputOp &out, const int &nthreads, const RangeType &r, const size_t step) {
//#pragma omp parallel sections num_threads(2) OMP_SHARE_DEFAULT
//  {
//#pragma omp section
//    {
//      out(nthreads, r, step);
//    }
//
//#pragma omp section
//    {
//      comp(nthreads, r, step);
//    }
//  }
//}
//
//
///**
// * retrieve data from Input and compute some final output using ComputeOp
// * Input is an operator, not just data structure.  example is a MPI Receiver
// *
// * Map or Reduce Pattern.
// *
// * @param in
// * @param comp
// * @param nthreads
// * @param r
// * @param step
// */
//template<typename InputOp, typename ComputeOp>
//void ConsumingPipeline(InputOp &in, ComputeOp &comp, const int &nthreads, const RangeType &r, const size_t step) {
//#pragma omp parallel sections num_threads(numOps) OMP_SHARE_DEFAULT
//  {
//#pragma omp section
//    {
//      in(nthreads, r, step);
//    }
//
//#pragma omp section
//    {
//      comp(nthreads, r, step);
//    }
//  }
//}
//
//////// TODO:  need to make the pipelines better.




template<typename OP, typename OT>
OT P2P(OP &op, const int &nthreads, size_t &_count) {

  OT v = 0;
  size_t count = 0;


#pragma omp parallel num_threads(nthreads) shared(op) OMP_SHARE_DEFAULT reduction(+:v, count)
  {
//      INFOF("MPI rank: %d/%d OMP thread: %d/%d\n", (rank+1), nprocs, (tid+1), nthreads);
//    INFOF("thread %d range: %lu, %lu\n", omp_get_thread_num(), r.start, r.end);

//    OT tv = 0;
    bool done = false;

    int tid = 0;
#ifdef USE_OPENMP
    tid = omp_get_thread_num();
#endif
    do {

      // not sure if this is efficient or should a thread local variable be used?
      done = op(tid, count, v);
    } while (!done);
  }
  _count = count;
  return v;
}



template<typename OP, typename OT>
OT MasterSlave(OP &op, const int &nthreads, size_t &_count) {

  OT *v = new OT[nthreads];
  size_t *counts = new size_t[nthreads];
#pragma omp parallel for num_threads(nthreads) OMP_SHARE_DEFAULT shared(nthreads, v, counts) schedule(static)
  for (int i = 0; i < nthreads; ++i) {
    v[i] = 0;
    counts[i] = 0;
  }
  RangeType r = op.getRange();
  size_t step = op.getChunkSize();


#pragma omp parallel num_threads(nthreads) shared(r, v, counts, op, step) OMP_SHARE_DEFAULT
  {
//    INFOF("thread %d range: %lu, %lu\n", omp_get_thread_num(), r.start, r.end);
// wait or nowait?
#pragma omp single
    {
//#pragma omp task untied
      for (size_t i = r.start; i < r.end; i += step)
      {
#pragma omp task
        {
          int tid = 0;
#ifdef USE_OPENMP
          tid = omp_get_thread_num();
#endif
          op(tid, counts[tid], v[tid]);
        }

      } //for
    } // omp single
  } // omp parallel

  OT vo = 0;
  size_t count = 0;
#pragma omp parallel for num_threads(nthreads) OMP_SHARE_DEFAULT reduction(+:vo, count) shared(nthreads, v, counts) schedule(static)
  for (int i = 0; i < nthreads; ++i) {
    vo += v[i];
    count += counts[i];
  }

  delete [] v;
  delete [] counts;
  _count = count;
  return vo;
}

// TODO:  explicit queueing of input data read from disk.


template<typename OP, typename OT>
OT MasterSlaveNoWait(OP &op, const int &nthreads, size_t &_count) {

  OT *v = new OT[nthreads];
  size_t *counts = new size_t[nthreads];
#pragma omp parallel for num_threads(nthreads) OMP_SHARE_DEFAULT shared(nthreads, v, counts) schedule(static)
  for (int i = 0; i < nthreads; ++i) {
    v[i] = 0;
    counts[i] = 0;
  }
  RangeType r = op.getRange();
  size_t step = op.getChunkSize();


#pragma omp parallel num_threads(nthreads) shared(r, v, counts, op, step) OMP_SHARE_DEFAULT
  {
//    INFOF("thread %d range: %lu, %lu\n", omp_get_thread_num(), r.start, r.end);

// wait or nowait?
#pragma omp single nowait
    {
//#pragma omp task untied
      for (size_t i = r.start; i < r.end; i += step)
      {
#pragma omp task
        {
          int tid = 0;
#ifdef USE_OPENMP
          tid = omp_get_thread_num();
#endif
          op(tid, counts[tid], v[tid]);

//          INFOF("%d %lu %lf\n", omp_get_thread_num(), j, tv);
        }

      } //for
    } // omp single
  } // omp parallel

  OT vo = 0;
  size_t count = 0;
#pragma omp parallel for num_threads(nthreads) OMP_SHARE_DEFAULT reduction(+:vo, count) shared(nthreads, v, counts) schedule(static)
  for (int i = 0; i < nthreads; ++i) {
    vo += v[i];
    count += counts[i];
  }

  delete [] v;
  delete [] counts;
_count = count;
  return vo;
}


template<typename OP, typename OT>
OT ParFor(OP &op, const int &nthreads, size_t &_count) {

  OT v = 0;
  size_t count = 0;
  RangeType r = op.getRange();
  size_t step = op.getChunkSize();

#pragma omp parallel num_threads(nthreads) OMP_SHARE_DEFAULT shared(r, v, count, op, step)
  {
//    INFOF("thread %d range: %lu, %lu\n", omp_get_thread_num(), r.start, r.end);

    #pragma omp for schedule(guided) reduction(+ : v, count)
      for (size_t i = r.start; i < r.end; i+=step)
      {
        op(omp_get_thread_num(), count, v);
      }
  }
  _count = count;
  return v;
}


template<typename OP, typename OT>
OT Sequential(OP &op, const int &nthreads, size_t &_count) {

  OT v = 0;
  size_t count = 0;

//  INFOF("thread %d range: %lu, %lu\n", omp_get_thread_num(), r.start, r.end);
  RangeType r = op.getRange();
  size_t step = op.getChunkSize();


  for (size_t i = r.start; i < r.end; i+=step)
  {
    op(0, count, v);
  }
  _count = count;
  return v;
}



#endif /* OMP_PATTERNS_HPP_ */
