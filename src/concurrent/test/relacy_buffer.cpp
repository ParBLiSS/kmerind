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
template<typename T, bliss::concurrent::LockType TS, int n_threads, int64_t MDSize, int64_t CAP, int64_t FILL>
struct testAppend : rl::test_suite<
    testAppend<T, TS, n_threads, MDSize, CAP, FILL>,
    n_threads>
{
  bliss::io::Buffer<TS, CAP, MDSize> b1;
  static constexpr int64_t cap_in_elems = CAP / sizeof(T);
  static constexpr int64_t nelems = FILL / sizeof(T);

  std::array<int64_t, n_threads> lsuccess;
  std::array<int64_t, n_threads> lfail;
  std::array<int64_t, n_threads> lswap;
  std::array<std::vector<T>, n_threads > lgold;



  static constexpr int64_t expectedSuccess = (cap_in_elems > nelems ? nelems : cap_in_elems);
  static constexpr int64_t expectedFailure = (cap_in_elems > nelems ? 0 : nelems - cap_in_elems);
  static constexpr int64_t expectedSwap = (CAP > (nelems * sizeof(T)) ? 0 : 1);


    void before()
    {
      for (int i = 0; i < n_threads; ++i)
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
      for (int64_t i = thread_index; i < nelems; i += n_threads)
      {
        v = static_cast<T>(i);
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

      for (int i = 0; i < n_threads; ++i)
      {
        success += lsuccess[i];
        fail += lfail[i];
        swap += lswap[i];

        gold.insert(gold.end(), lgold[i].begin(), lgold[i].end());
      }

      bool same = compareUnorderedSequences(b1.operator T*(), gold.begin(), success);


      if ((success != expectedSuccess) || !same || (fail != expectedFailure) || (swap != expectedSwap) || ((swap > 0) != b1.is_read_only()))
        INFOF(
            "Test testAppend: LT %d nProd %d cap %ld fill %ld (%ld elems), actual %ld/%lu, failed %ld/%lu, swap %ld/%lu",
            TS, n_threads, CAP, FILL, nelems,
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
template<typename T, bliss::concurrent::LockType TS, int n_threads, int64_t MDSize, int64_t CAP, int64_t FILL, int64_t BlockAt, int64_t UnblockAt>
struct testBlockedAppend : rl::test_suite<
    testBlockedAppend<T, TS, n_threads, MDSize, CAP, FILL, BlockAt, UnblockAt>,
    n_threads>
{
  bliss::io::Buffer<TS, CAP, MDSize> b1;
  static constexpr int64_t cap_in_elems = CAP / sizeof(T);
  static constexpr int64_t nelems = FILL / sizeof(T);
  static constexpr int64_t nelems1 = BlockAt / sizeof(T);
  static constexpr int64_t nelems2 = UnblockAt / sizeof(T);

  std::array<int64_t, n_threads> lsuccess;
  std::array<int64_t, n_threads> lfail;
  std::array<int64_t, n_threads> lswap;
  std::array<std::vector<T>, n_threads > lgold;


  // because the threads may interleave in various ways. it's possible that we have fully inserted buffer,
  // (all other threads finish, 1 thread blocks near end, and a second thread immediately unblocks, both finish remainders)
  //  but it's not possible to have 0 success, because the blocking thread has to reach the block point from
  //  unblocked state.
  static constexpr int64_t expectedSuccessMin = (nelems < nelems1 / n_threads) ? nelems : nelems1 / n_threads;   // all threads fail except for part of 1
  static constexpr int64_t expectedFailureMin = (nelems > cap_in_elems ? nelems - cap_in_elems : 0);
  static constexpr int64_t expectedSuccessMax = cap_in_elems;  // all threads succeed except for part of 1.
  static constexpr int64_t expectedFailureMax = nelems - expectedSuccessMin;


    void before()
    {
      for (int i = 0; i < n_threads; ++i)
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

      for (int64_t i = thread_index; i < nelems; i+= n_threads)
      {

        if (i == nelems1) {
          b1.block_writes();
        }
        if (i == nelems2) {
          b1.unblock_writes();
        }

        v = static_cast<T>(i);
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

      for (int i = 0; i < n_threads; ++i)
      {
        success += lsuccess[i];
        fail += lfail[i];
        swap += lswap[i];

        gold.insert(gold.end(), lgold[i].begin(), lgold[i].end());
      }

      bool same = compareUnorderedSequences(b1.operator T*(), gold.begin(), success);


      if ((success < expectedSuccessMin || success > expectedSuccessMax) || !same || (fail < expectedFailureMin || fail > expectedFailureMax) )
        INFOF(
            "Test testAppend: LT %d nProd %d cap %ld (%ld elems), blockat %ld, unblockat %lu, total elems %ld, success %ld/[%lu-%lu], failed %ld/[%lu-%lu], swap %ld",
            TS, n_threads, CAP, cap_in_elems, BlockAt, UnblockAt, nelems,
            success, expectedSuccessMin, expectedSuccessMax, fail, expectedFailureMin, expectedFailureMax, swap);
      assert(success >= expectedSuccessMin && success <= expectedSuccessMax);
      assert(same);
      assert(fail >= expectedFailureMin && fail <= expectedFailureMax);

    }

    void invariant()
    {

    }

};



template<bliss::concurrent::LockType LT, int n_threads>
void simulate(rl::test_params p)
{
  rl::simulate<testAppend<int, LT, n_threads, 0, 128, 64> >(p);
  rl::simulate<testAppend<int, LT, n_threads, 0, 128, 128> >(p);
  rl::simulate<testAppend<int, LT, n_threads, 0, 128, 128 + sizeof(int)> >(p);

  rl::simulate<testAppend<int, LT, n_threads, 0, 127, 64> >(p);
  rl::simulate<testAppend<int, LT, n_threads, 0, 127, 127> >(p);
  rl::simulate<testAppend<int, LT, n_threads, 0, 127, 128> >(p);

  rl::simulate<testAppend<int, LT, n_threads, 1,   128, 128> >(p);
  rl::simulate<testAppend<int, LT, n_threads, 2,   128, 128> >(p);
  rl::simulate<testAppend<int, LT, n_threads, 4,   128, 128> >(p);
  rl::simulate<testAppend<int, LT, n_threads, 7,   128, 128> >(p);
  rl::simulate<testAppend<int, LT, n_threads, 8,   128, 128> >(p);
  rl::simulate<testAppend<int, LT, n_threads, 16,  128, 128> >(p);
  rl::simulate<testAppend<int, LT, n_threads, 200, 128, 128> >(p);

  rl::simulate<testBlockedAppend<int, LT, n_threads, 0, 128, 128,  0, 64> >(p);
  rl::simulate<testBlockedAppend<int, LT, n_threads, 0, 128, 128,  0, 128> >(p);
  rl::simulate<testBlockedAppend<int, LT, n_threads, 0, 128, 128, 32, 32> >(p);
  rl::simulate<testBlockedAppend<int, LT, n_threads, 0, 128, 128, 32, 64> >(p);
  rl::simulate<testBlockedAppend<int, LT, n_threads, 0, 128, 128, 32, 128> >(p);

  rl::simulate<testBlockedAppend<int, LT, n_threads, 0, 128, 16, 32, 64> >(p);
  rl::simulate<testBlockedAppend<int, LT, n_threads, 0, 128, 48, 32, 64> >(p);
  rl::simulate<testBlockedAppend<int, LT, n_threads, 0, 128, 96, 32, 64> >(p);
  rl::simulate<testBlockedAppend<int, LT, n_threads, 0, 128, 144, 32, 64> >(p);
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
