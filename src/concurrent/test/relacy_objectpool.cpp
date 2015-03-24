/**
 * @file		check_bufferpool.cpp
 * @ingroup
 * @author	Tony Pan <tpan7@gatech.edu>
 * @brief
 * @details
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */

#include <vector>
#include <cstdlib>   // for rand

#include "utils/logging.h"

#include "utils/iterator_test_utils.hpp"
#include "concurrent/buffer.hpp"
#include "concurrent/referenced_object_pool.hpp"
#include "concurrent/unreferenced_object_pool.hpp"



// note that relacy simulate threads.
template<typename T, typename PoolType, int nThreads, int64_t Iter, int64_t PCAP = std::numeric_limits<int64_t>::max() >
struct testBufferSwap : rl::test_suite<
  testBufferSwap<T, PoolType, nThreads, Iter, PCAP >,
    nThreads>
{
  using BufferType = typename PoolType::ObjectType;


  std::array<int64_t, nThreads> lsuccess;
  std::array<int64_t, nThreads> lfail;
  std::array<int64_t, nThreads> lswap;
  std::array<std::vector<T>, nThreads > lgold;
  std::array<std::vector<T>, nThreads > lstored;

  PoolType pool;
  VAR_T(BufferType*) buf_ptr;


  testBufferSwap() : pool(PCAP), buf_ptr(nullptr) {}

    void before()
    {
      for (int i = 0; i < nThreads; ++i)
      {
        lsuccess[i] = 0;
        lfail[i] = 0;
        lswap[i] = 0;
      }

      VAR(buf_ptr) = pool.tryAcquireObject();
      if (VAR(buf_ptr)) VAR(buf_ptr)->clear_and_unblock_writes();
    }

    void thread(unsigned thread_index)
    {
      T v;
      unsigned int r;
      for (int64_t i = thread_index; i < Iter; i+= nThreads)
      {

        if (VAR(buf_ptr)) {
          v = static_cast<T>(i);
          r = VAR(buf_ptr)->append(&v, sizeof(T));
        } else {
          r = 0x0;
        }

        if ((r & 0x1) > 0) {
          ++lsuccess[thread_index];

          lgold[thread_index].push_back(v);

        } else
          ++lfail[thread_index];

        if ((r & 0x2) > 0) {
          ++lswap[thread_index];

          auto new_buf_ptr = pool.tryAcquireObject();
          if (new_buf_ptr) new_buf_ptr->clear_and_unblock_writes();
          auto lptr = VAR(buf_ptr);
          VAR(buf_ptr) = new_buf_ptr;

          if (lptr) lstored[thread_index].insert(lstored[thread_index].end(),
                                                lptr->operator T*(), lptr->operator T*() + lptr->getSize() / sizeof(T));


          pool.releaseObject(lptr);

        }
      }
    }

    void after()
    {
      int64_t success = 0;
      int64_t fail = 0;
      int64_t swap = 0;
      std::vector<T> gold;
      std::vector<T> stored;

      for (int i = 0; i < nThreads; ++i)
      {
        success += lsuccess[i];
        fail += lfail[i];
        swap += lswap[i];

        gold.insert(gold.end(), lgold[i].begin(), lgold[i].end());
        stored.insert(stored.end(), lstored[i].begin(), lstored[i].end());
      }

      if (VAR(buf_ptr)) {
        VAR(buf_ptr)->block_and_flush();
        stored.insert(stored.end(),
                      VAR(buf_ptr)->operator T*(), VAR(buf_ptr)->operator T*() + VAR(buf_ptr)->getSize() / sizeof(T));

        pool.releaseObject(VAR(buf_ptr));
      }


      bool same = compareUnorderedSequences(stored.begin(), gold.begin(), success);
      assert(same);
      assert(success == stored.size());
      assert(success != 0);

    }

    void invariant()
    {

    }

};




template<typename PoolType, int nThreads, int64_t Iter, int64_t PoolCap = std::numeric_limits<int64_t>::max()>
struct testPoolAcquire : rl::test_suite<
  testPoolAcquire<PoolType, nThreads, Iter, PoolCap >,
    nThreads>
{

  int64_t expectedSuccess, expectedFailure;

  std::array<int64_t, nThreads> lsuccess;
  std::array<int64_t, nThreads> lfail;

  PoolType pool;
  std::array<std::vector<typename PoolType::ObjectType *>, nThreads > objs;

  testPoolAcquire() : pool(PoolCap) {}

    void before()
    {
      for (int i = 0; i < nThreads; ++i)
      {
        lsuccess[i] = 0;
        lfail[i] = 0;
      }

      expectedSuccess = pool.isUnlimited() ? Iter :
          (Iter > pool.getCapacity() ? pool.getCapacity() : Iter);
      expectedFailure = pool.isUnlimited() ? 0 :
          (Iter > pool.getCapacity() ? Iter - pool.getCapacity() : 0);
    }

    void thread(unsigned thread_index)
    {
      for (int64_t i = thread_index; i < Iter ; i += nThreads)
      {
        auto ptr = pool.tryAcquireObject();
        if (ptr) {
          ++lsuccess[thread_index];
          objs[thread_index].push_back(ptr);
        } else {
          ++lfail[thread_index];
        }
      }
    }

    void after()
    {
      int64_t success = 0;
      int64_t fail = 0;

      for (int i = 0; i < nThreads; ++i)
      {
        success += lsuccess[i];
        fail += lfail[i];

        for (int j = 0; j < objs[i].size(); ++j) {
          pool.releaseObject(objs[i][j]);
        }
      }

      if (success != expectedSuccess || fail != expectedFailure)
    	  ERRORF("success %ld != expected %ld?  failed %ld != expectedFailure %ld ?", success, expectedSuccess, fail, expectedFailure);
      assert(success == expectedSuccess);
      assert(fail == expectedFailure);

      pool.reset();
    }

    void invariant()
    {

    }
};


template<typename PoolType, int nThreads, int64_t Iter, int64_t PoolCap = std::numeric_limits<int64_t>::max()>
struct testPoolRelease : rl::test_suite<
  testPoolRelease<PoolType, nThreads, Iter, PoolCap>,
    nThreads>
{
  int64_t expectedSuccess, expectedFailure;

  std::array<int64_t, nThreads> lsuccess;
  std::array<int64_t, nThreads> lfail;

  PoolType pool;
  std::array<std::vector<typename PoolType::ObjectType *>, nThreads > objs;
  int count;

  testPoolRelease() : pool(PoolCap), count(0) {}

    void before()
    {
      for (int i = 0; i < nThreads; ++i)
      {
        lsuccess[i] = 0;
        lfail[i] = 0;
      }

      count = 0;
      expectedSuccess = pool.isUnlimited() ? Iter :
          (Iter > pool.getCapacity() ? pool.getCapacity() : Iter);
      expectedFailure = expectedSuccess;

      for (int i = 0; i < Iter; ++i ) {
        auto ptr = pool.tryAcquireObject();
        if (ptr) {
          objs[i % nThreads].push_back(ptr);
          objs[(i + 1) % nThreads].push_back(ptr);
          ++count;
        }
      }
      assert(count == expectedSuccess);
    }

    void thread(unsigned thread_index)
    {
      for ( auto v : objs[thread_index] )
      {
        if (pool.releaseObject(v)) {
          ++lsuccess[thread_index];
        } else {
          ++lfail[thread_index];
        }
      }
    }

    void after()
    {
      int64_t success = 0;
      int64_t fail = 0;

      for (int i = 0; i < nThreads; ++i)
      {
        success += lsuccess[i];
        fail += lfail[i];

        objs[i].clear();
      }

      if (success != expectedSuccess || fail != expectedFailure)
    	  ERRORF("success %ld != expected %ld?  failed %ld != expectedFailure %ld ?", success, expectedSuccess, fail, expectedFailure);
      assert(success == expectedSuccess);
      assert(fail == expectedFailure);

      pool.reset();
    }

    void invariant()
    {

    }
};


template<bliss::concurrent::LockType LT, int nThreads,
	typename std::enable_if<nThreads == 1 &&
	(LT == bliss::concurrent::LockType::MUTEX ||
	 LT == bliss::concurrent::LockType::SPINLOCK), int>::type = 0>
void simulate(rl::test_params p)
{
		using PoolType = bliss::concurrent::ObjectPool< bliss::io::Buffer<bliss::concurrent::LockType::NONE, 8192>, LT, true >;

	  rl::simulate<testPoolAcquire<PoolType, nThreads, 200> >(p);
	  rl::simulate<testPoolAcquire<PoolType, nThreads, 100, 150> >(p);
	  rl::simulate<testPoolAcquire<PoolType, nThreads, 200, 150> >(p);
	  rl::simulate<testPoolRelease<PoolType, nThreads, 200> >(p);
	  rl::simulate<testPoolRelease<PoolType, nThreads, 100, 150> >(p);
	  rl::simulate<testPoolRelease<PoolType, nThreads, 200, 150> >(p);

	  rl::simulate<testBufferSwap<int, PoolType, nThreads, 200 > >(p);
	  rl::simulate<testBufferSwap<int, PoolType, nThreads, 100, 150 > >(p);
	  rl::simulate<testBufferSwap<int, PoolType, nThreads, 200, 150 > >(p);
}

template<bliss::concurrent::LockType LT, int nThreads,
	typename std::enable_if<nThreads == 1 &&
	(LT == bliss::concurrent::LockType::LOCKFREE ||
	 LT == bliss::concurrent::LockType::NONE), int>::type = 0>
void simulate(rl::test_params p)
{
using PoolType = bliss::concurrent::ObjectPool< bliss::io::Buffer<bliss::concurrent::LockType::NONE, 8192>, LT, false >;

  	  rl::simulate<testPoolAcquire<PoolType, nThreads, 200> >(p);
	  rl::simulate<testPoolAcquire<PoolType, nThreads, 100, 150> >(p);
	  rl::simulate<testPoolAcquire<PoolType, nThreads, 200, 150> >(p);

	  rl::simulate<testPoolRelease<PoolType, nThreads, 200> >(p);
	  rl::simulate<testPoolRelease<PoolType, nThreads, 100, 150> >(p);
	  rl::simulate<testPoolRelease<PoolType, nThreads, 200, 150> >(p);

	  rl::simulate<testBufferSwap<int, PoolType, nThreads, 200 > >(p);
	  rl::simulate<testBufferSwap<int, PoolType, nThreads, 100, 150 > >(p);
	  rl::simulate<testBufferSwap<int, PoolType, nThreads, 200, 150 > >(p);
}


template<bliss::concurrent::LockType LT, int nThreads,
	typename std::enable_if<nThreads != 1 &&
							(LT == bliss::concurrent::LockType::MUTEX ||
							 LT == bliss::concurrent::LockType::SPINLOCK), int>::type = 0>
void simulate(rl::test_params p)
{
	using PoolType = bliss::concurrent::ObjectPool< bliss::io::Buffer<bliss::concurrent::LockType::LOCKFREE, 8192>, LT, true >;
  rl::simulate<testPoolAcquire<PoolType, nThreads, 200> >(p);
  rl::simulate<testPoolAcquire<PoolType, nThreads, 100, 150> >(p);
  rl::simulate<testPoolAcquire<PoolType, nThreads, 200, 150> >(p);

  rl::simulate<testPoolRelease<PoolType, nThreads, 200> >(p);
  rl::simulate<testPoolRelease<PoolType, nThreads, 100, 150> >(p);
  rl::simulate<testPoolRelease<PoolType, nThreads, 200, 150> >(p);

  rl::simulate<testBufferSwap<int, PoolType, nThreads, 200 > >(p);
  rl::simulate<testBufferSwap<int, PoolType, nThreads, 100, 150 > >(p);
  rl::simulate<testBufferSwap<int, PoolType, nThreads, 200, 150 > >(p);
}


template<bliss::concurrent::LockType LT, int nThreads,
	typename std::enable_if<nThreads != 1 &&
							(LT == bliss::concurrent::LockType::LOCKFREE ||
							 LT == bliss::concurrent::LockType::NONE), int>::type = 0>
void simulate(rl::test_params p)
{
	using PoolType = bliss::concurrent::ObjectPool< bliss::io::Buffer<bliss::concurrent::LockType::LOCKFREE, 8192>, LT, false >;

  rl::simulate<testPoolAcquire<PoolType, nThreads, 200> >(p);
  rl::simulate<testPoolAcquire<PoolType, nThreads, 100, 150> >(p);
  rl::simulate<testPoolAcquire<PoolType, nThreads, 200, 150> >(p);

  rl::simulate<testPoolRelease<PoolType, nThreads, 200> >(p);
  rl::simulate<testPoolRelease<PoolType, nThreads, 100, 150> >(p);
  rl::simulate<testPoolRelease<PoolType, nThreads, 200, 150> >(p);

  rl::simulate<testBufferSwap<int, PoolType, nThreads, 200 > >(p);
  rl::simulate<testBufferSwap<int, PoolType, nThreads, 100, 150 > >(p);
  rl::simulate<testBufferSwap<int, PoolType, nThreads, 200, 150 > >(p);
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
#elif defined( BLISS_LOCKFREE )
  constexpr bliss::concurrent::LockType lt = bliss::concurrent::LockType::LOCKFREE;
#elif defined(BLISS_NONE)
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
