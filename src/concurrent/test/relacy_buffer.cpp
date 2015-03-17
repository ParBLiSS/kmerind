/**
 * @file		check_buffer.cpp
 * @ingroup
 * @author	Tony Pan <tpan7@gatech.edu>
 * @brief
 * @details
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */


#include "concurrent/buffer.hpp"

#include <vector>

#include "utils/iterator_test_utils.hpp"
#include "utils/logging.h"




// note that relacy simulate threads.
template<typename T, bliss::concurrent::LockType TS, int nProducers, int64_t MDSize, int64_t CAP, int64_t FILL>
struct testAppend : rl::test_suite<
    testAppend<T, TS, nProducers, MDSize, CAP, FILL>,
    nProducers>
{
  bliss::io::Buffer<TS, CAP, MDSize> b1;
  static constexpr int64_t nelems = CAP / sizeof(T);
  static constexpr int64_t target = FILL / sizeof(T);
  static constexpr int64_t attempts = ((target / nProducers) * nProducers);

  std::array<int64_t, nProducers> lsuccess;
  std::array<int64_t, nProducers> lfail;
  std::array<int64_t, nProducers> lswap;
  std::array<std::vector<T>, nProducers > lgold;



  static constexpr int64_t expectedSuccess = (nelems > attempts ? attempts : nelems);
  static constexpr int64_t expectedFailure = (nelems > attempts ? 0 : attempts - nelems);
  static constexpr int64_t expectedSwap = (CAP > (attempts * sizeof(T)) ? 0 : 1);


    void before()
    {
      for (int i = 0; i < nProducers; ++i)
      {
        lsuccess[i] = 0;
        lfail[i] = 0;
        lswap[i] = 0;
      }

      b1.clear_and_unblock_writes();
    }

    void thread(unsigned thread_index)
    {
      T v;
      int r;
      for (int64_t i = 0; i < target / nProducers; ++i)
      {
        v = static_cast<T>(i * nProducers + thread_index);
        r = b1.append(&v, sizeof(T));
        if ((r & 0x1) > 0) {
          ++lsuccess[thread_index];

          lgold[thread_index].push_back(v);

        } else
          ++lfail[thread_index];

        if ((r & 0x2) > 0) {
          ++lswap[thread_index];
        }
      }
    }

    void after()
    {
      int64_t success = 0;
      int64_t fail = 0;
      int64_t swap = 0;
      std::vector<T> gold;

      for (int i = 0; i < nProducers; ++i)
      {
        success += lsuccess[i];
        fail += lfail[i];
        swap += lswap[i];

        gold.insert(gold.end(), lgold[i].begin(), lgold[i].end());
      }

      bool same = compareUnorderedSequences(b1.operator T*(), gold.begin(), success);


      if ((success != expectedSuccess) || !same || (fail != expectedFailure) || (swap != expectedSwap) || ((swap > 0) != b1.is_read_only()))
        INFOF(
            "Test testAppend: LT %d nProd %d cap %ld fill %ld (%ld elems), attempts %ld, actual %ld/%lu, failed %ld/%lu, swap %ld/%lu",
            TS, nProducers, CAP, FILL, target, attempts,
            success, expectedSuccess, fail, expectedFailure, swap, expectedSwap);
      assert(success == expectedSuccess);
      assert(same);
      assert(fail == expectedFailure);
      assert(swap == expectedSwap);
      assert((swap > 0) == b1.is_read_only());

    }

    void invariant()
    {

    }

};



// note that relacy simulate threads.
template<typename T, bliss::concurrent::LockType TS, int nProducers, int64_t MDSize, int64_t CAP, int64_t BlockAt, int64_t UnblockAt>
struct testBlockedAppend : rl::test_suite<
    testBlockedAppend<T, TS, nProducers, MDSize, CAP, BlockAt, UnblockAt>,
    nProducers>
{
  bliss::io::Buffer<TS, CAP, MDSize> b1;
  static constexpr int64_t nelems = CAP / sizeof(T);
  static constexpr int64_t target1 = BlockAt / sizeof(T);
  static constexpr int64_t target2 = UnblockAt / sizeof(T);
  static constexpr int64_t attempts = ((nelems / nProducers) * nProducers);

  std::array<int64_t, nProducers> lsuccess;
  std::array<int64_t, nProducers> lfail;
  std::array<int64_t, nProducers> lswap;
  std::array<std::vector<T>, nProducers > lgold;


  // because the threads may interleave in various ways. it's possible that we have 0 failure,
  // (all other threads finish, 1 thread blocks near end, and a second thread immediately unblocks, both finish remainders)
  //  but it's not possible to have 0 success, because the blocking thread has to reach the block point from
  //  unblocked state.
  static constexpr int64_t expectedSuccessMin = target1 / nProducers;   // all threads fail except for part of 1
  static constexpr int64_t expectedFailureMin = 0;
  static constexpr int64_t expectedSuccessMax = attempts - expectedFailureMin;  // all threads succeed except for part of 1.
  static constexpr int64_t expectedFailureMax = attempts - expectedSuccessMin;


    void before()
    {
      for (int i = 0; i < nProducers; ++i)
      {
        lsuccess[i] = 0;
        lfail[i] = 0;
        lswap[i] = 0;
      }

      b1.clear_and_unblock_writes();
    }

    void thread(unsigned thread_index)
    {
      T v;
      int r;
      int64_t eid;
      for (int64_t i = 0; i < nelems / nProducers; ++i)
      {
        eid = i * nProducers + thread_index;

        if (eid == target1) {
          b1.block_writes();
        }
        if (eid == target2) {
          b1.unblock_writes();
        }

        v = static_cast<T>(eid);
        r = b1.append(&v, sizeof(T));
        if ((r & 0x1) > 0) {
          ++lsuccess[thread_index];

          lgold[thread_index].push_back(v);

        } else
          ++lfail[thread_index];

        if ((r & 0x2) > 0) {
          ++lswap[thread_index];
        }
      }
    }

    void after()
    {
      int64_t success = 0;
      int64_t fail = 0;
      int64_t swap = 0;
      std::vector<T> gold;

      for (int i = 0; i < nProducers; ++i)
      {
        success += lsuccess[i];
        fail += lfail[i];
        swap += lswap[i];

        gold.insert(gold.end(), lgold[i].begin(), lgold[i].end());
      }

      bool same = compareUnorderedSequences(b1.operator T*(), gold.begin(), success);


      if ((success < expectedSuccessMin || success > expectedSuccessMax) || !same || (fail < expectedFailureMin || fail > expectedFailureMax) )
        INFOF(
            "Test testAppend: LT %d nProd %d cap %ld (%ld elems), blockat %ld, unblockat %lu, attempts %ld, success %ld/[%lu-%lu], failed %ld/[%lu-%lu], swap %ld",
            TS, nProducers, CAP, nelems, BlockAt, UnblockAt, attempts,
            success, expectedSuccessMin, expectedSuccessMax, fail, expectedFailureMin, expectedFailureMax, swap);
      assert(success >= expectedSuccessMin && success <= expectedSuccessMax);
      assert(same);
      assert(fail >= expectedFailureMin && fail <= expectedFailureMax);

    }

    void invariant()
    {

    }

};



template<bliss::concurrent::LockType LT, int nProducers>
void simulate(rl::test_params p)
{
  rl::simulate<testAppend<int, LT, nProducers, 0, 128, 64> >(p);
  rl::simulate<testAppend<int, LT, nProducers, 0, 128, 128> >(p);
  rl::simulate<testAppend<int, LT, nProducers, 0, 128, 128 + sizeof(int)> >(p);

  rl::simulate<testAppend<int, LT, nProducers, 0, 127, 64> >(p);
  rl::simulate<testAppend<int, LT, nProducers, 0, 127, 127> >(p);
  rl::simulate<testAppend<int, LT, nProducers, 0, 127, 128> >(p);

  rl::simulate<testAppend<int, LT, nProducers, 1,   128, 128> >(p);
  rl::simulate<testAppend<int, LT, nProducers, 2,   128, 128> >(p);
  rl::simulate<testAppend<int, LT, nProducers, 4,   128, 128> >(p);
  rl::simulate<testAppend<int, LT, nProducers, 7,   128, 128> >(p);
  rl::simulate<testAppend<int, LT, nProducers, 8,   128, 128> >(p);
  rl::simulate<testAppend<int, LT, nProducers, 16,  128, 128> >(p);
  rl::simulate<testAppend<int, LT, nProducers, 200, 128, 128> >(p);

  rl::simulate<testBlockedAppend<int, LT, nProducers, 0, 128,  0, 64> >(p);
  rl::simulate<testBlockedAppend<int, LT, nProducers, 0, 128,  0, 128> >(p);
  rl::simulate<testBlockedAppend<int, LT, nProducers, 0, 128, 32, 32> >(p);
  rl::simulate<testBlockedAppend<int, LT, nProducers, 0, 128, 32, 64> >(p);
  rl::simulate<testBlockedAppend<int, LT, nProducers, 0, 128, 32, 128> >(p);

}



int main(int argc, char** argv) {



  int iterations = 1000;
  if (argc > 1)
  {
    iterations = atoi(argv[1]);
  }

  rl::test_params p;
  p.search_type = rl::sched_random;
  p.iteration_count = iterations;
  p.execution_depth_limit = 4000;



#if defined( BLISS_MUTEX)
  constexpr bliss::concurrent::LockType lt = bliss::concurrent::LockType::MUTEX;
#elif defined(BLISS_SPINLOCK)
  constexpr bliss::concurrent::LockType lt = bliss::concurrent::LockType::SPINLOCK;
#elif defined(BLISS_LOCKFREE)
  constexpr bliss::concurrent::LockType lt = bliss::concurrent::LockType::LOCKFREE;
#else //if defined(BLISS_NONE)
  constexpr bliss::concurrent::LockType lt = bliss::concurrent::LockType::NONE;
#endif


  simulate<lt, 1>(p);

#if defined(BLISS_MUTEX) || defined(BLISS_SPINLOCK) || defined(BLISS_LOCKFREE)

  simulate<lt, 2>(p);
  simulate<lt, 3>(p);
  simulate<lt, 4>(p);
  simulate<lt, 8>(p);

#endif




}
