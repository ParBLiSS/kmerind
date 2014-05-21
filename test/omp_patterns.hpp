/**
 * @file		omp_patterns.hpp
 * @ingroup
 * @author	tpan
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

#include <iterators/range.hpp>
typedef bliss::iterator::range<size_t> RangeType;


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
//#pragma omp parallel sections num_threads(3) default(none)
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
//#pragma omp parallel sections num_threads(2) default(none)
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
//#pragma omp parallel sections num_threads(numOps) default(none)
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




template<typename OP, typename OT = typename std::result_of<OP(size_t, size_t, size_t&)>::type >
OT P2P(OP &op, const int &nthreads, const RangeType &r, const size_t step, size_t &_count) {

  size_t i = r.start;
  OT v = 0;
  size_t count = 0;


#pragma omp parallel num_threads(nthreads) shared(r, i, op) default(none) reduction(+:v, count)
  {
//      printf("MPI rank: %d/%d OMP thread: %d/%d\n", (rank+1), nprocs, (tid+1), nthreads);
//    printf("thread %d range: %lu, %lu\n", omp_get_thread_num(), r.start, r.end);

//    OT tv = 0;
    size_t j = 0;
    bool done = false;

    do {

#pragma omp critical (assign)
      {
        j = i;
        i += step;
      }
      if (j >= r.end) {
        done = true;
      } else {

        // not sure if this is efficient or should a thread local variable be used?
        v += op(j, std::min(j + step, r.end), count);
      }
    } while (!done);
  }
  _count = count;
  return v;
}


template<typename OP, typename OT = typename std::result_of<OP(size_t, size_t, size_t&)>::type>
OT P2P_atomic(OP &op, const int &nthreads, const RangeType &r, const size_t step, size_t &_count) {


  size_t i = r.start;
  OT v = 0;
  size_t count = 0;


#pragma omp parallel num_threads(nthreads) shared(r, i, op) default(none) reduction(+:v, count)
  {
//      printf("MPI rank: %d/%d OMP thread: %d/%d\n", (rank+1), nprocs, (tid+1), nthreads);
//    printf("thread %d range: %lu, %lu\n", omp_get_thread_num(), r.start, r.end);

    size_t j = 0;
    bool done = false;

    do {

#pragma omp atomic capture
      {
        j = i;
        i += step;
      }

      if (j >= r.end) {
        done = true;
      } else {
        v += op(j, std::min(j + step, r.end), count);
      }
    } while (!done);

  }
  _count = count;

  return v;
}


template<typename OP, typename OT = typename std::result_of<OP(size_t, size_t, size_t&)>::type>
OT MasterSlave(OP &op, const int &nthreads, const RangeType &r, const size_t step, size_t &_count) {

  OT *v = new OT[nthreads];
  size_t *counts = new size_t[nthreads];
#pragma omp parallel for num_threads(nthreads) default(none) shared(nthreads, v, counts) schedule(static)
  for (int i = 0; i < nthreads; ++i) {
    v[i] = 0;
    counts[i] = 0;
  }

#pragma omp parallel num_threads(nthreads) shared(r, v, counts, op) default(none)
  {
//    printf("thread %d range: %lu, %lu\n", omp_get_thread_num(), r.start, r.end);
// wait or nowait?
#pragma omp single
    {
#pragma omp task untied
      for (size_t i = r.start; i < r.end; i += step)
      {
#pragma omp task
        {
#ifdef USE_OPENMP
          v[omp_get_thread_num()] += op(i, std::min(i +step, r.end), counts[omp_get_thread_num()]);
#else
          v[0] += op(i, std::min(i + step, r.end), counts[0]);
#endif
//          printf("%d %lu %lf\n", omp_get_thread_num(), j, tv);
        }

      } //for
    } // omp single
  } // omp parallel

  OT vo = 0;
  size_t count = 0;
#pragma omp parallel for num_threads(nthreads) default(none) reduction(+:vo, count) shared(nthreads, v, counts) schedule(static)
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


template<typename OP, typename OT = typename std::result_of<OP(size_t, size_t, size_t&)>::type>
OT MasterSlaveNoWait(OP &op, const int &nthreads, const RangeType &r, const size_t step, size_t &_count) {

  OT *v = new OT[nthreads];
  size_t *counts = new size_t[nthreads];
#pragma omp parallel for num_threads(nthreads) default(none) shared(nthreads, v, counts) schedule(static)
  for (int i = 0; i < nthreads; ++i) {
    v[i] = 0;
    counts[i] = 0;
  }

#pragma omp parallel num_threads(nthreads) shared(r, v, counts, op) default(none)
  {
//    printf("thread %d range: %lu, %lu\n", omp_get_thread_num(), r.start, r.end);

// wait or nowait?
#pragma omp single nowait
    {
#pragma omp task untied
      for (size_t i = r.start; i < r.end; i += step)
      {
#pragma omp task
        {
#ifdef USE_OPENMP
          int tid = omp_get_thread_num();
          v[tid] += op(i, std::min(i + step, r.end), counts[tid]);

#else
          v[0] += op(i, std::min(i + step, r.end), counts[0]);
#endif
//          printf("%d %lu %lf\n", omp_get_thread_num(), j, tv);
        }

      } //for
    } // omp single
  } // omp parallel

  OT vo = 0;
  size_t count = 0;
#pragma omp parallel for num_threads(nthreads) default(none) reduction(+:vo, count) shared(nthreads, v, counts) schedule(static)
  for (int i = 0; i < nthreads; ++i) {
    vo += v[i];
    count += counts[i];
  }

  delete [] v;
  delete [] counts;
_count = count;
  return vo;
}


template<typename OP, typename OT = typename std::result_of<OP(size_t, size_t, size_t&)>::type>
OT ParFor(OP &op, const int &nthreads, const RangeType &r, const size_t step, size_t &_count) {

  OT v = 0;
  size_t count = 0;

#pragma omp parallel num_threads(nthreads) default(none) shared(r, v, count, op)
  {
//    printf("thread %d range: %lu, %lu\n", omp_get_thread_num(), r.start, r.end);

    #pragma omp for schedule(guided) reduction(+ : v, count)
      for (size_t i = r.start; i < r.end; i+=step)
      {
        v += op(i, std::min(i + step, r.end), count);
      }
  }
  _count = count;
  return v;
}


template<typename OP, typename OT = typename std::result_of<OP(size_t, size_t, size_t&)>::type>
OT Sequential(OP &op, const int &nthreads, const RangeType &r, const size_t step, size_t &_count) {

  OT v = 0;
  size_t count = 0;

//  printf("thread %d range: %lu, %lu\n", omp_get_thread_num(), r.start, r.end);

  for (size_t i = r.start; i < r.end; i+=step)
  {
    v += op(i, std::min(i + step, r.end), count);
  }
  _count = count;
  return v;
}



#endif /* OMP_PATTERNS_HPP_ */
