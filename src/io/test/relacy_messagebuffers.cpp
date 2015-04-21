

#include "config/relacy_config.hpp"
#include "utils/logging.h"
#include "utils/iterator_test_utils.hpp"
#include "io/message_buffers.hpp"
#include "concurrent/concurrent.hpp"

#include "concurrent/buffer.hpp"
#include "concurrent/referenced_object_pool.hpp"

#include <sstream>



// note that relacy simulate threads.
template<typename T, typename MBType, int commSize, int nThreads, int64_t nelems>
struct mbAppend : rl::test_suite<
  mbAppend<T, MBType, commSize, nThreads, nelems>,
    nThreads>
{
  using BufferType = typename MBType::BufferType;
  using PoolType = typename MBType::BufferPoolType;


  std::array<int64_t, nThreads> lsuccess;
  std::array<int64_t, nThreads> lfail;
  std::array<int64_t, nThreads> lsswap;
  std::array<int64_t, nThreads> lfswap;
  std::array<int64_t, nThreads> lupdating;

  std::array< std::vector<T>, nThreads > lgold;
  std::array< std::vector<VAR_T(BufferType*)>, nThreads > lstored;

  PoolType pool;

  MBType buffers;


  mbAppend() : pool(), buffers(pool, commSize, nThreads) {}

    void before()
    {
      for (int i = 0; i < nThreads; ++i)
      {
        lsuccess[i] = 0;
        lfail[i] = 0;
        lfswap[i] = 0;
        lsswap[i] = 0;
        lupdating[i] = 0;
      }

    }

    void thread(unsigned thread_index)
    {
      T data;
      bool op_suc;
      BufferType* ptr;

      for (int64_t i = thread_index; i < nelems; i+= nThreads)
      {
        data = static_cast<T>(i);
        T * data_remain = nullptr;
        size_t count_remain = 0;

        std::tie(op_suc, ptr) = buffers.append(&data, 1, data_remain, count_remain, 0, thread_index);

        if (op_suc) {
          ++lsuccess[thread_index];

          lgold[thread_index].push_back(data);

          if (ptr) {
            ++lsswap[thread_index];
            if ( ptr->is_writing() ) ++lupdating[thread_index];

            lstored[thread_index].push_back(ptr);


//            std::stringstream ss;
//            T* vals = ptr->operator T*();
//            int count = ptr->getSize()/sizeof(T);
//            std::copy(vals, vals + count, std::ostream_iterator<T>(ss, ", "));
//            DEBUGF("swapped buffer: [%s]", ss.str().c_str());
          }

        } else {
          ++lfail[thread_index];

          if (ptr) {
            ++lfswap[thread_index];
            if ( ptr->is_writing() ) ++lupdating[thread_index];


            lstored[thread_index].push_back(ptr);

//            std::stringstream ss;
//            T* vals = ptr->operator T*();
//            int count = ptr->getSize()/sizeof(T);
//            std::copy(vals, vals + count, std::ostream_iterator<T>(ss, ", "));
//            DEBUGF("swapped buffer: [%s]", ss.str().c_str());
          }

        }
      }
    }

    void after()
    {
      int64_t success = 0;
      int64_t fail = 0;
      int64_t sswap = 0;
      int64_t fswap = 0;
      int64_t updating = 0;
      std::vector<T> gold;
      gold.clear();
      std::vector<T> stored;
      stored.clear();

      int fullCount = 0;

      for (int i = 0; i < nThreads; ++i)
      {
        success += lsuccess[i];
        fail += lfail[i];
        sswap += lsswap[i];
        fswap += lfswap[i];
        updating += lupdating[i];

        gold.insert(gold.end(), lgold[i].begin(), lgold[i].end());
        fullCount += lstored[i].size();


        for ( auto ptr : lstored[i]) {
          if (VAR(ptr)) {
            T* vals = VAR(ptr)->operator T*();
            int count = VAR(ptr)->getSize()/sizeof(T);
            stored.insert(stored.end(), vals, vals+count);
            pool.releaseObject(VAR(ptr));
          }
        }
      }
//      std::stringstream ss;
//      std::copy(gold.begin(), gold.end(), std::ostream_iterator<T>(ss, ", "));
//      DEBUGF("gold: [%s]", ss.str().c_str());


      std::vector<BufferType*> finals = buffers.flushBufferForRank(0);
      for (auto ptr : finals) {
        if (ptr) {
          T* vals = ptr->operator T*();
          int count = ptr->getSize()/sizeof(T);
          stored.insert(stored.end(), vals, vals+count);
          pool.releaseObject(ptr);
        }
      }


      int expectedElemsInFullBuffer = (BufferType::getCapacity()/sizeof(T));
      bool exactFill = (BufferType::getCapacity() % sizeof(T)) == 0;

      int expectedFullMin = success / expectedElemsInFullBuffer - nThreads;
      int expectedFullMax = success / expectedElemsInFullBuffer;

      //INFOF("success: %ld, fail %ld, sswap %ld, fswap %ld, updating %ld", success, fail, sswap, fswap, updating);

      bool same = compareUnorderedSequences(stored.begin(), gold.begin(), success);
      assert(same);
      assert(success == stored.size());

      if (exactFill)
        assert((fail == sswap) && (sswap == fullCount) && (sswap >= expectedFullMin)  && (sswap <= expectedFullMax) && (fswap == 0));
      else
        assert((fail == fswap) && (fswap == fullCount) && (fswap >= expectedFullMin)  && (fswap <= expectedFullMax) && (sswap == 0));
      assert(updating == 0);
    }

    void invariant()
    {

    }

};



// note that relacy simulate threads.
template<typename T, typename MBType, int commSize, int nThreads, int64_t nelems>
struct mbTest : rl::test_suite<
  mbTest<T, MBType, commSize, nThreads, nelems>,
    nThreads>
{
  using BufferType = typename MBType::BufferType;
  using PoolType = typename MBType::BufferPoolType;


  std::array<int64_t, nThreads> lsuccess;
  std::array<int64_t, nThreads> ldatacount;
  std::array<int64_t, nThreads> lfail;
  std::array<int64_t, nThreads> lsswap;
  std::array<int64_t, nThreads> lfswap;
  std::array<int64_t, nThreads> lupdating;
  std::array< std::vector<T>, nThreads > lgold;
  std::array< std::vector<T>, nThreads > lstored;

  PoolType pool;

  MBType buffers;


  mbTest() : pool(), buffers(pool, commSize, nThreads) {}

    void before()
    {
      for (int i = 0; i < nThreads; ++i)
      {
        lsuccess[i] = 0;
        ldatacount[i] = 0;
        lfail[i] = 0;
        lfswap[i] = 0;
        lsswap[i] = 0;
        lupdating[i] = 0;
      }

    }

    void thread(unsigned thread_index)
    {
      T data;
      bool op_suc;
      BufferType* ptr;

      for (int64_t i = thread_index; i < nelems; i+= nThreads)
      {
        data = static_cast<T>(i);
        T * data_remain = nullptr;
        size_t count_remain = 0;

        std::tie(op_suc, ptr) = buffers.append(&data, 1, data_remain, count_remain, 0, thread_index);

        if (op_suc) {
          ++lsuccess[thread_index];

          lgold[thread_index].push_back(data);

          if (ptr) {
            ++lsswap[thread_index];
            if ( ptr->is_writing() ) ++lupdating[thread_index];

            int count = ptr->getSize() / sizeof(T);
            ldatacount[thread_index] += count;
            T* d = ptr->operator T*();
            lstored[thread_index].insert(lstored[thread_index].end(), d, d + count);

            pool.releaseObject(ptr);
          }

        } else {
          ++lfail[thread_index];

          if (ptr) {
            ++lfswap[thread_index];
            if ( ptr->is_writing() ) ++lupdating[thread_index];

            int count = ptr->getSize() / sizeof(T);
            ldatacount[thread_index] += count;
            T* d = ptr->operator T*();
            lstored[thread_index].insert(lstored[thread_index].end(), d, d + count);

            pool.releaseObject(ptr);
          }

        }
      }
    }

    void after()
    {
      int64_t success = 0;
      int64_t datacount = 0;
      int64_t fail = 0;
      int64_t sswap = 0;
      int64_t fswap = 0;
      int64_t updating = 0;

      std::vector<T> gold;
      std::vector<T> stored;


      for (int i = 0; i < nThreads; ++i)
      {
        success += lsuccess[i];
        datacount += ldatacount[i];
        fail += lfail[i];
        sswap += lsswap[i];
        fswap += lfswap[i];
        updating += lupdating[i];

        gold.insert(gold.end(), lgold[i].begin(), lgold[i].end());

        stored.insert(stored.end(), lstored[i].begin(), lstored[i].end());
      }


      std::vector<BufferType*> finals = buffers.flushBufferForRank(0);
      for (auto ptr : finals) {
        if (ptr) {
          T* vals = ptr->operator T*();
          int count = ptr->getSize()/sizeof(T);
          datacount += count;
          stored.insert(stored.end(), vals, vals+count);
          pool.releaseObject(ptr);
        }

      }

      int expectedElemsInFullBuffer = (BufferType::getCapacity()/sizeof(T));
      bool exactFill = (BufferType::getCapacity() % sizeof(T)) == 0;

      int expectedFullMin = success / expectedElemsInFullBuffer - nThreads;
      int expectedFullMax = success / expectedElemsInFullBuffer;
      //INFOF("success: %ld, fail %ld, sswap %ld, fswap %ld, updating %ld", success, fail, sswap, fswap, updating);

      bool same = compareUnorderedSequences(stored.begin(), gold.begin(), success);
      assert(same);
      assert(datacount == success);
      assert(success == stored.size());
      //assert(success == nelems);
      if (exactFill)
        assert((fail == sswap)  && (sswap >= expectedFullMin)  && (sswap <= expectedFullMax) && (fswap == 0));
      else
        assert((fail == fswap) && (fswap >= expectedFullMin)  && (fswap <= expectedFullMax) && (sswap == 0));
      assert(updating == 0);

    }

    void invariant()
    {

    }

};



// note that relacy simulate threads.
template<typename T, typename MBType, int commSize, int nThreads, int64_t nelems>
struct mbEnsureInsertTest : rl::test_suite<
  mbEnsureInsertTest<T, MBType, commSize, nThreads, nelems>,
    nThreads>
{
  using BufferType = typename MBType::BufferType;
  using PoolType = typename MBType::BufferPoolType;


  std::array<int64_t, nThreads> lsuccess;
  std::array<int64_t, nThreads> ldatacount;
  std::array<int64_t, nThreads> lsswap;
  std::array<int64_t, nThreads> lupdating;
  std::array< std::vector<T>, nThreads > lgold;
  std::array< std::vector<T>, nThreads > lstored;

  PoolType pool;

  MBType buffers;


  mbEnsureInsertTest() : pool(), buffers(pool, commSize, nThreads) {}

    void before()
    {
      for (int i = 0; i < nThreads; ++i)
      {
        lsuccess[i] = 0;
        ldatacount[i] = 0;
        lsswap[i] = 0;
        lupdating[i] = 0;
      }

    }

    void thread(unsigned thread_index)
    {
      T data;
      bool op_suc;
      BufferType* ptr;

      for (int64_t i = thread_index; i < nelems; i+= nThreads)
      {
        data = static_cast<T>(i);

        do {
          T * data_remain = nullptr;
          size_t count_remain = 0;

          std::tie(op_suc, ptr) = buffers.append(&data, 1, data_remain, count_remain, 0, thread_index);

          if (ptr) {
            ++lsswap[thread_index];
            if ( ptr->is_writing() ) ++lupdating[thread_index];

            int count = ptr->getSize() / sizeof(T);
            ldatacount[thread_index] += count;
            T* d = ptr->operator T*();
            lstored[thread_index].insert(lstored[thread_index].end(), d, d + count);

            pool.releaseObject(ptr);
          }

        } while ( !op_suc );

        if (op_suc) {
          ++lsuccess[thread_index];

          lgold[thread_index].push_back(data);
        }


      }
    }

    void after()
    {
      int64_t success = 0;
      int64_t datacount = 0;
      int64_t sswap = 0;
      int64_t updating = 0;

      std::vector<T> gold;
      std::vector<T> stored;


      for (int i = 0; i < nThreads; ++i)
      {
        success += lsuccess[i];
        datacount += ldatacount[i];
        sswap += lsswap[i];
        updating += lupdating[i];

        gold.insert(gold.end(), lgold[i].begin(), lgold[i].end());

        stored.insert(stored.end(), lstored[i].begin(), lstored[i].end());
      }


      std::vector<BufferType*> finals = buffers.flushBufferForRank(0);
      for (auto ptr : finals) {
        if (ptr) {
          T* vals = ptr->operator T*();
          int count = ptr->getSize()/sizeof(T);
          stored.insert(stored.end(), vals, vals+count);
          pool.releaseObject(ptr);
        }

      }

      int expectedElemsInFullBuffer = (BufferType::getCapacity()/sizeof(T));

      int expectedFullMin = success / expectedElemsInFullBuffer - nThreads;
      int expectedFullMax = success / expectedElemsInFullBuffer;

      //INFOF("success: %ld, sswap %ld, updating %ld", success, sswap, updating);

      bool same = compareUnorderedSequences(stored.begin(), gold.begin(), success);
      assert(same);
      assert(datacount == success);
      assert(success == stored.size());
      assert(success == nelems);
      assert((sswap >= expectedFullMin)  && (sswap <= expectedFullMax) );
      assert(updating == 0);

    }

    void invariant()
    {

    }

};


template<bliss::concurrent::LockType LT, int commSize, int nThreads>
void simulate(rl::test_params parm)
{
  using BufferType = bliss::io::Buffer<bliss::concurrent::LockType::NONE, 2047, 0>;
  using PoolType = bliss::concurrent::ObjectPool<BufferType, LT>;
  using MessageBuffersType = bliss::io::SendMessageBuffers<bliss::concurrent::LockType::THREADLOCAL, PoolType>;

  rl::simulate<mbAppend<int, MessageBuffersType, commSize, nThreads,  nThreads * 1536> >(parm);
//  rl::simulate<mbTest<int, MessageBuffersType, commSize,  nThreads,  nThreads * 1536> >(parm);
//  rl::simulate<mbEnsureInsertTest<int, MessageBuffersType, commSize,  nThreads,  nThreads * 1536> >(parm);

}



int main(int argc, char *argv[])
{

  int iterations = 1000;
  if (argc > 1)
  {
    iterations = atoi(argv[1]);
  }

  rl::test_params parm;
  parm.search_type = rl::sched_random;
  parm.iteration_count = iterations;
  parm.execution_depth_limit = 4000;


#if defined( BLISS_MUTEX)
  constexpr bliss::concurrent::LockType LT = bliss::concurrent::LockType::MUTEX;
#elif defined(BLISS_SPINLOCK)
  constexpr bliss::concurrent::LockType LT= bliss::concurrent::LockType::SPINLOCK;
#endif



  // set up MPI


  simulate<LT, 1, 1>(parm);
  simulate<LT, 1, 2>(parm);
  simulate<LT, 1, 3>(parm);
  simulate<LT, 1, 4>(parm);
  simulate<LT, 1, 8>(parm);
  simulate<LT, 1, 16>(parm);

  simulate<LT, 2, 4>(parm);
  simulate<LT, 3, 4>(parm);
  simulate<LT, 4, 4>(parm);

  // finalize MPI
  return 0;
}

