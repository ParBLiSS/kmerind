/**
 * @file		buffer_pool.hpp
 * @ingroup bliss::io
 * @author	tpan
 * @brief   this file defines in-memory Buffer Pool.
 * @details
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef BUFFERPOOL_HPP_
#define BUFFERPOOL_HPP_

#include <cassert>

#include <mutex>
#include <unordered_set>
#include <vector>
#include <stdexcept>

#include "utils/logging.h"
#include "io/buffer.hpp"
#include "concurrent/concurrent.hpp"



namespace bliss
{
  namespace io
  {
    // TODO: move constructor and assignment operators between BufferPools of different thread safeties.

    /**
     * @class     BufferPool
     * @brief     BufferPool manages a set of reusable, in-memory buffers.  In particular, it can limit the amount of memory usage.
     * @details   The class is templated to provide thread-safe or unsafe BufferPools, managing thread-safe or unsafe Buffers
     *
     *            Each buffer is a block of preallocated memory that can be appended into.
     *            Each Buffer instance is referenced by id using the [] or at() accessor functions
     *
     *            The caller can acquire a new buffer from the pool and release an old buffer back into the pool, by Id.
     *            Released buffers are marked as empty and available.
     *
     *            When the pool is exhausted, Acquire function returns false.
     *
     *      life cycle:  pool acquire():  buffer unblock().  -> allow rw.
     *              application:  buffer block() -> allow ro
     *              pool release():  buffer clear().
     *
     *
     * @note      Each buffer should be acquired by a single thread and released by a single thread.
     *            Note that Buffer does not actualy clear the memory, just resets the size, which is also the pointer for insertion.
     *            While the release method handles duplicate releases, there may be race conditions when multithreading, including
     *              likely loss of data (thread 1 appends data and thread 2 releases  it.
     *
     *            the pools RETAINS OWNERSHIP of the buffers.
     *
     *
     *
     *
     *  @tparam PoolTS    The thread safety property for the pool
     *  @tparam BufferTS  The thread safety property of each Buffer.
     */
    template<bliss::concurrent::ThreadSafety PoolTS, bliss::concurrent::ThreadSafety BufferTS = PoolTS>
    class BufferPool
    {
      public:
        /**
         * @brief     Type of Buffer in the BufferPool (thread safe or not)
         */
        typedef bliss::io::Buffer<BufferTS>               BufferType;
        /**
         * @brief     Index Type used to reference the Buffers in the BufferPool.
         */
        typedef int                                       IdType;

        /// constant, represent when a pool element or Buffer is not present.
        static constexpr IdType ABSENT = -1;

      protected:
        /**
         * @brief     capacity of the BufferPool (in number of Buffers)
         */
        IdType capacity;
        /**
         * @brief     capacity of the individual Buffers (in bytes)
         */
        size_t buffer_capacity;

        /**
         * @brief     Internal Set of available Buffers for immediate use.
         * @details   Using set instead of ThreadSafeQueue to ensure uniqueness of Buffer Ids.
         */
        std::unordered_set<IdType>                          available;

        /**
         * @brief     Internal Vector of all Buffers.  May be growable if the Pool was specified as such.
         */
        std::vector<BufferType> buffers;

        /**
         * @brief     boolean flag indicating whether the Pool is fixed size or growable.
         */
        bool fixedSize;

        /**
         * @brief     mutex to control access.
         */
        mutable std::mutex mutex;


        /**
         * @brief     Create a new Buffer object and put into the BufferPool internal vector.
         * @note      THREAD SAFE
         * @return    Id of the newly created Buffer.
         */
        template<bliss::concurrent::ThreadSafety TS = PoolTS>
        typename std::enable_if<TS == bliss::concurrent::THREAD_SAFE, IdType >::type createNewBuffer()
        {
          std::lock_guard<std::mutex> lock(mutex);
          printf("creating new Thread Safe Buffer\n");
          IdType bufferId = static_cast<IdType>(buffers.size());
          buffers.push_back(std::move(BufferType(buffer_capacity)));

          return bufferId;
        }
        /**
         * @brief     Get the Id of the next available Buffer.
         * @note      THREAD SAFE
         * @return    Id of the next available Buffer.
         */
        template<bliss::concurrent::ThreadSafety TS = PoolTS>
        typename std::enable_if<TS == bliss::concurrent::THREAD_SAFE, std::pair<bool, IdType> >::type getNextAvailable() {
          std::pair<bool, IdType> result;
          result.first = false;
          std::unique_lock<std::mutex> lock(mutex);
          if (!available.empty()) {
            result.first = true;
            result.second = *(available.begin());
            available.erase(result.second);
          }
          lock.unlock();
          return std::move(result);
        }

        /**
         * @brief     Create a new Buffer object and put into the BufferPool internal vector.
         * @note      THREAD UNSAFE
         * @return    Id of the newly created Buffer.
         */
        template<bliss::concurrent::ThreadSafety TS = PoolTS>
        typename std::enable_if<TS == bliss::concurrent::THREAD_UNSAFE, IdType >::type createNewBuffer()
        {
          printf("creating new Thread Unsafe Buffer\n");
          IdType bufferId = buffers.size();
          buffers.push_back(std::move(BufferType(buffer_capacity)));

          return bufferId;
        }
        /**
         * @brief     Get the Id of the next available Buffer.
         * @note      THREAD UNSAFE
         * @return    Id of the next available Buffer.
         */
        template<bliss::concurrent::ThreadSafety TS = PoolTS>
        typename std::enable_if<TS == bliss::concurrent::THREAD_UNSAFE, std::pair<bool, IdType> >::type getNextAvailable() {
          std::pair<bool, IdType> result;
          result.first = false;
          if (!available.empty()) {
            result.first = true;
            result.second = *(available.begin());
            available.erase(result.second);
          }
          return std::move(result);
        }


        /**
         * @brief     Private move constructor with mutex lock.
         * @param other   Source BufferPool object to move
         * @param l       Mutex lock on the Source BufferPool object
         */
        BufferPool(BufferPool<PoolTS, BufferTS>&& other, const std::lock_guard<std::mutex>& l) :
          capacity(other.capacity), buffer_capacity(other.buffer_capacity), available(std::move(other.available)),
          buffers(std::move(other.buffers)), fixedSize(other.fixedSize) {
          other.capacity = 0;
          other.buffer_capacity = 0;
          other.available.clear();
        };

//        /**
//         * @brief verifies that the id is valid - throws exception.
//         * @param id    Id of the buffer to check.
//         */
//        void checkId(const IdType& id) const {
//          if (id < 0 || id >= static_cast<IdType>(buffers.size()))
//            throw std::out_of_range("Error:  buffer id was specified with out of range value");
//        }


      public:
        /**
         * @brief     construct a Buffer Pool with buffers of capacity _buffer_capacity.  the number of buffers in the pool is set to _pool_capacity.
         *
         * @param _pool_capacity       number of buffers in the pool
         * @param _buffer_capacity     size of the individual buffers
         */
        BufferPool(const IdType _pool_capacity, const size_t _buffer_capacity) :
          capacity(_pool_capacity), buffer_capacity(_buffer_capacity), available(),   // ThreadSafeQueue is not size bound.
          buffers(), fixedSize(_pool_capacity != std::numeric_limits<IdType>::max()) {

          // get an estimated size first, so don't have to keep growing the vector
          IdType size_hint = fixedSize ? capacity : 4;
          //IdType size_hint = capacity;

          // reserve the buffer, and configure.  this part is not thread safe
          buffers.reserve(size_hint);
          //DEBUGF("RESERVING: %d for buffer pool", size_hint);
          for (IdType i = 0; i < size_hint; ++i) {
            buffers.push_back(std::move(BufferType(buffer_capacity)));
            buffers[i].block();
            releaseBuffer<PoolTS>(i);
          }
        };

        /**
         * @brief     default constructor is deleted.
         */
        BufferPool() = delete;

        /**
         * @brief     construct a Buffer pool with buffers of capacity _buffer_capacity.  the number of buffers in the pool is unbounded.
         *
         * @param _buffer_capacity    size of the individual buffers
         */
        explicit BufferPool(const size_t _buffer_capacity) :
            BufferPool<PoolTS, BufferTS>(std::numeric_limits<IdType>::max(), _buffer_capacity) {};

        /**
         * @brief     default copy constructor is deleted.
         * @param other   source BufferPool object to copy from.
         */
        explicit BufferPool(const BufferPool<PoolTS, BufferTS>& other) = delete;

        /**
         * @brief      Move constructor.  Delegates to the Move constructor that locks access on the source.
         * @param other    source BufferPool object to move.
         */
        explicit BufferPool(BufferPool<PoolTS, BufferTS>&& other) :
          BufferPool<PoolTS, BufferTS>(std::move(other), std::lock_guard<std::mutex>(other.mutex)) {};

        /**
         * @brief     default copy assignment operator is deleted.
         * @param other   source BufferPool object to copy from.
         * @return        self, with member variables copied from other.
         */
        BufferPool<PoolTS, BufferTS>& operator=(const BufferPool<PoolTS, BufferTS>& other) = delete;

        /**
         * @brief     move assignment operator.
         * @param other   source BufferPool object to move from.
         * @return        self, with member variables moved from other.
         */
        BufferPool<PoolTS, BufferTS>& operator=(BufferPool<PoolTS, BufferTS>&& other) {
          std::unique_lock<std::mutex> mylock(mutex, std::defer_lock),
                                        otherlock(other.mutex, std::defer_lock);
          std::lock(mylock, otherlock);

          capacity = other.capacity; other.capacity = 0;
          buffer_capacity = other.buffer_capacity; other.buffer_capacity = 0;
          fixedSize = other.fixedSize;
          available = std::move(other.available);
          buffers = std::move(other.buffers);
          return *this;
        };


        /**
         * @brief     default destructor
         */
        virtual ~BufferPool() {};

        /**
         * @brief     Current size of the BufferPool
         * @return  size, type IdType (aka int).
         */
        const IdType getSize() const  {
          return buffers.size();
        }

        /**
         * @brief     Current capacity of the BufferPool
         * @return    capacity, type IdTyype (aka int)
         */
        const IdType getCapacity() const {
          return capacity;
        }

        /**
         * @brief     Each buffer's maximum capacity.
         * @return  each buffer's maximum capacity.
         */
        const size_t getBufferCapacity() const {
          return buffer_capacity;
        }

        /**
         * @brief     whether the BufferPool is growable
         * @return    bool indicating if BufferPool is growable.
         */
        const bool isFixedSize() const {
          return fixedSize;
        }

        /**
         * @brief     Resets the entire BufferPool: all Buffers in pool are  marked as released and available.
         *
         * @note      This is not entirely thread safe.  the available set is cleared by a single thread, but other
         *            threads may be acquiring buffers while they are being released.
         *
         *            It is envisioned that this function should be called from a single thread.
         *            however, care must be taken to ensure that all other threads have completed their work, else data loss is likely.
         */
        void reset() {
          std::unique_lock<std::mutex> lock(mutex);
          available.clear();
          lock.unlock();
          for (IdType i = 0; i < static_cast<IdType>(buffers.size()); ++i) {
            this->buffers[i].block();
            releaseBuffer<PoolTS>(i);
          }
        }

        /**
         * check if the buffer specified by id is available for reuse.
         * @param id
         * @return
         */
        bool isAvailable(const IdType& id) {
          return available.find(id) != available.end();
        }


        /**
         * @brief     Get the next available Buffer by id.  if none available, return false in the first argument of the pair.
         * @return    std::pair with bool and Id,  first indicate whether acquire was successful, second indicate the BufferId if successful.
         */
        std::pair<bool, IdType> tryAcquireBuffer() {
          std::pair<bool, IdType> output = std::move(getNextAvailable<PoolTS>());
          //DEBUGF("acquired? %d, available count %ld, buffers count %ld", (output.first ? output.second : -1), available.size(), buffers.size());
          if (!output.first && !fixedSize) {
            // can grow and none available.  so allocate a new one and return.
            output.first = true;
            output.second = createNewBuffer<PoolTS>();

            //DEBUGF("acquired final %d, available count %ld, buffers count %ld", output.second, available.size(), buffers.size());
          }
          // clear just before use.  buffer also gets unblocked.  This allows a 3 stage life cycle: actively used (rw), blocked (ro), released (no access)
          // clear transitions from released to active use.  rw->ro should be done by application logic.  blocked to released is also done by application logic.
          this->buffers[output.second].unblock();
          assert(this->buffers[output.second].isEmpty());

          return output;
        }

        /*
         *  waitAndAcquireBuffer does not really make sense for non-threadsafe.  not going to have it for threadsafe version either
         *   - if there is nothing available in a single thread, waiting
         *  for the same thread to update does not make sense.
         * @param bufferId
         */
//        /**
//         *
//         * @param bufferId
//         */
//        IdType waitAndAcquireBuffer() {
//          IdType bufferId;
//          if (fixedSize) {
//            // have to wait for available
//            bufferId = available.waitAndPop();
//
//          } else {  // can grow, but first see if we can reuse one.
//
//            // check if there are some available.
//            auto newBuf = available.tryPop();
//            if (!newBuf.first) {
//              std::lock_guard<std::mutex> lock(mutex);
//              bufferId = createNewBuffer();
//            } else {
//              bufferId = newBuf.second;
//            }
//          }
//          return bufferId;
//        }

        /**
         * @brief     Release a buffer back to pool, by id.  if id is incorrect, throw std exception.  clears buffer before release back into available.
         * @note      THREAD SAFE
         * @param bufferId    The id of the Buffer to be released.
         */
        template<bliss::concurrent::ThreadSafety TS = PoolTS>
        typename std::enable_if<TS == bliss::concurrent::THREAD_SAFE, void>::type releaseBuffer(const IdType& id) {
          //checkId(id);

          std::lock_guard<std::mutex> lock(mutex);
          assert(this->buffers[id].isBlocked());
          this->buffers.at(id).clear();
          available.insert(id);

          //DEBUGF("buffer id %d is available for reuse.", id);
        }


        /**
         * @brief     Release a buffer back to pool, by id.  if id is incorrect, throw std exception.  clears buffer before release back into available.
         * @note      THREAD UNSAFE
         * @param bufferId    The id of the Buffer to be released.
         */
        template<bliss::concurrent::ThreadSafety TS = PoolTS>
        typename std::enable_if<TS == bliss::concurrent::THREAD_UNSAFE, void>::type releaseBuffer(const IdType& id) {
          //checkId(id);

          assert(this->buffers[id].isBlocked());
          this->buffers.at(id).clear();
          available.insert(id);

          //DEBUGF("buffer id %d is available for reuse.", id);
        }

        /**
         * @brief     Access the Buffer by id.  if Id is outside the range of the current BufferPool size, throw exception.
         *
         * @note  multiple threads can access the same id, and one thread can access multiple Ids.
         *
         * @param bufferId    Id of the buffer to be access (get from Acquire)
         * @return            const reference to the Buffer
         */
        const BufferType& operator[](const IdType& id) const {
          // can only access what has been allocated.  does not increase size of the pool.
          //checkId(id);

          return buffers.at(id);               // multiple threads can access the same buffer...
                                                // bounds checking is done here.
        }

        /**
         * @brief     Access the Buffer by id.  if Id is incorrect, throw exception.
         *
         * @note      multiple threads can access the same id, and one thread can access multiple ids.
         *
         * @param bufferId    Id of the buffer to be access (get from Acquire)
         * @return            const reference to the Buffer
         */
        const BufferType& at(const IdType& id) const {
          return buffers.at(id);
        }
    };

    template<bliss::concurrent::ThreadSafety PoolTS, bliss::concurrent::ThreadSafety BufferTS>
        constexpr typename BufferPool<PoolTS, BufferTS>::IdType BufferPool<PoolTS, BufferTS>::ABSENT;

  } /* namespace io */
} /* namespace bliss */

#endif /* BUFFERPOOL_HPP_ */
