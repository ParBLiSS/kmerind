/**
 * @file		MessageBuffers.hpp
 * @ingroup bliss::io
 * @author	tpan
 * @brief   MessageBuffers base class and SendMessageBuffers subclass for buffering data for MPI send/receive
 * @details SendMessageBuffers is a collection of in-memory Buffers that containing data that
 *            1. will be sent to remote destinations, such as a MPI process
 *            2. was received from remote destinations, such as a MPI process
 *
 *          In contrast to BufferPool, whose purpose is to prevent local process from blocking waiting for a buffer to clear,
 *          MessageBuffer's purpose is to batch the data to be transmitted, using BufferPool to prevent local process from blocking.
 *
 *          As there may be several different patterns of interactions between MessageBuffers and other MessageBuffers, MessageBuffers and
 *          local thread/process, there is a base class MessageBuffers that acts as proxy to the underlying BufferPool (memory management).
 *          Two subclasses are planned:  SendMessageBuffers and RecvMessageBuffers, but only SendMessageBuffers has been implemented.
 *
 *          Envisioned patterns of interaction:
 *            Send:     MessageBuffers -> mpi send
 *            Recv:     mpi recv -> MessageBuffers
 *            Forward:  mpi recv -> mpi send
 *            Local:    MessageBuffers -> MessageBuffers
 *
 *          The SendMessageBuffers accumulates data in its internal buffer, and notifies the caller when a buffer is full.
 *
 *          The (unimplemented) RecvMessageBuffers requires an event notification mechanism that works with callbacks to handle the received messages.
 *          This is currently implemented in CommunicationLayer.
 *
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef MESSAGEBUFFERS_HPP_
#define MESSAGEBUFFERS_HPP_

#include <vector>
#include <type_traits>
#include <typeinfo>
#include <iostream>
#include <sstream>
#include <mutex>
#include <condition_variable>
#include <xmmintrin.h>

#include "config.hpp"
#include "concurrent/buffer.hpp"
#include "concurrent/object_pool.hpp"
#include "concurrent/copyable_atomic.hpp"

#include <atomic>

#include <iterator> // for ostream_iterator

#include "omp.h"

namespace bliss
{
  namespace io
  {

    /**
     * @class			MessageBuffers
     * @brief     a data structure to buffering/batching messages for communication
     * @details   Templated base class for a collection of buffers that are used for communication.
     *
     *
     *            The class contains a BufferPool, which contains an array of reusable Buffer objects for in-memory storage.
     *
     * @tparam LockType   Indicates whether the class should be thread safe or not.
     */
    template<bliss::concurrent::LockType PoolLT, bliss::concurrent::LockType BufferLT = PoolLT, int64_t BufferCapacity = 8192>
    class MessageBuffers
    {
        // TODO: move consturctor and assignment operator to copy between thread safety levels.
      public:

        typedef Buffer<BufferLT, BufferCapacity> BufferType;

      protected:

        /// Internal BufferPool Type.  typedefed only to shorten the usage.
        typedef typename bliss::io::ObjectPool<PoolLT, BufferType>  BufferPoolType;

      public:
        /// IdType of a Buffer, aliased from BufferPool
        typedef typename BufferPoolType::ObjectPtrType            BufferPtrType;

      protected:
        /// a pool of in-memory Buffers for storage.
        BufferPoolType pool;

        /**
         * Gets the capacity of the Buffer instances.
         * @return    capacity of each buffer.
         */
        const size_t getBufferCapacity() const {
          return BufferCapacity;
        }

        /**
         * @brief Constructor.  creates a buffer given a specified per-buffer capacity, and the BufferPool capacity.
         * @note  Protected so base class function is not called directly.
         *
         *        limit BufferPool to have unlimited capacity, so as to ensure that we don't get nullptrs and thus don't need to
         *        check in code.
         *
         *
         * @param buffer_capacity   the Buffer's maximum capacity.  default to 8192.
         * @param pool_capacity     the BufferPool's capacity.  default to unlimited.
         */
        explicit MessageBuffers() : pool() {};  // unlimited size pool

        /**
         * @brief default copy constructor.  deleted.  since internal BufferPool does not allow copy construction/assignment.
         * @note Protected so base class function is not called directly.
         *
         * @param other     source MessageBuffers to copy from
         */
        explicit MessageBuffers(const MessageBuffers<PoolLT, BufferLT, BufferCapacity>& other) = delete;

        /**
         * @brief default copy assignment operator.  deleted.   since internal BufferPool does not allow copy construction/assignment.
         * @note Protected so base class function is not called directly.
         *
         * @param other     source MessageBuffers to copy from
         * @return          self.
         */
        MessageBuffers<PoolLT, BufferLT, BufferCapacity>& operator=(const MessageBuffers<PoolLT, BufferLT, BufferCapacity>& other) = delete;


        /**
         * @brief default move constructor.
         * @note Protected so base class function is not called directly.
         *
         * @param other     source MessageBuffers to move from
         */
        explicit MessageBuffers(MessageBuffers<PoolLT, BufferLT, BufferCapacity>&& other) : pool(std::move(other.pool)) {};


        /**
         * @brief default move assignment operator.
         * @note Protected so base class function is not called directly.
         *
         * @param other     source MessageBuffers to move from
         * @return          self.
         */
        MessageBuffers<PoolLT, BufferLT, BufferCapacity>& operator=(MessageBuffers<PoolLT, BufferLT, BufferCapacity>&& other) {
          pool = std::move(other.pool);
          return *this;
        }


      public:
        /**
         * @brief default destructor
         */
        virtual ~MessageBuffers() {};

        /**
         * @brief Releases a Buffer by it's BufferPool id, after the buffer is no longer needed.
         * @note  this should be called on a buffer that is not being used by MessageBuffers, i.e. flush does not satisfy this.
         *
         * @note  This to be called after the communication logic is done with the Send or Recv buffer.
         *
         * @param id    Buffer's id (assigned from BufferPool)
         */
        virtual void releaseBuffer(BufferPtrType ptr) {
          pool.releaseObject(ptr);
        }

      protected:
        /**
         * @brief Convenience method to release and clear all current Buffers in the pool.
         */
        virtual void reset() = 0;
    };



    /**
     * @class SendMessageBuffers
     * @brief SendMessageBuffers is a subclass of MessageBuffers designed to manage the actual buffering of data for a set of messaging targets.
     * @details
     *
     *  The class is designed with the following requirements in mind:
     *  1. data is destined for some remote process identifiable by a process id (e.g. rank)
     *  2. data is appended to buffer incrementally and with minimal blocking
     *  3. calling thread has a mechanism to process a full buffer. e.g. send it
     *  4. all data in a SendMessageBuffers are homogeneously typed.
     *
     *  The class is implemented to support the requirement:
     *  1. SendMessageBuffers is not aware of its data's type or metadata specifying its type.
     *  2. SendMessageBuffers uses the base MessageBuffers class' internal BufferPool to reuse memory
     *  3. SendMessageBuffers stores a vector of Buffer Ids (as assigned from BufferPool), thus mapping from process Rank to BufferId.
     *  4. SendMessageBuffers provides an Append function to incrementally add data to the Buffer object for the targt process rank.
     *  5. When a particular buffer is full, return the Buffer Id to the calling thread for it to process (send), and swap in an empty
     *      Buffer from BufferPool.
     *
     *  An internal vector is used to map from target process rank to Buffer pointer.  For multithreading flavors of message buffers,
     *    atomic<BufferType*> is used.
     *
     *  @note
     *  pool size is unlimited so we would not get nullptr, and therefore do not need to check it.
     *
     * @note
     *  For Buffers that are full, the calling thread will need to track the full buffers' Ids as this class evicts a full Buffer's id from it's
     *    vector of Buffer Ids.  The Calling Thread can do t his via a queue, as in the case of CommunicationLayer.
     *  Note that a Buffer is "in use" from the first call to the Buffer's append function to the completion of the send (e.g. via MPI Isend)
     *
     *  @tparam LockType    determines if this class should be thread safe
     */
    template<bliss::concurrent::LockType PoolLT, bliss::concurrent::LockType BufferLT = PoolLT, int64_t BufferCapacity = 8192>
    class SendMessageBuffers;

    /// DISABLED BECAUSE SWAPPING BUFFER PTR IS NOT SAFE in multiple threaded code.
//    template<bliss::concurrent::LockType PoolLT, bliss::concurrent::LockType BufferLT = PoolLT, int64_t BufferCapacity = 8192>
//    class SendMessageBuffers : public MessageBuffers<PoolLT, BufferLT, BufferCapacity>
//    {
//      protected:
//        using BaseType = MessageBuffers<PoolLT, BufferLT, BufferCapacity>;
//        using MyType = SendMessageBuffers<PoolLT, BufferLT, BufferCapacity>;
//
//      public:
//        /// Id type of the Buffers
//        using BufferType = typename BaseType::BufferType;
//        using BufferPtrType = typename BaseType::BufferPtrType;
//
//        /// the BufferPoolType from parent class
//        using BufferPoolType = typename BaseType::BufferPoolType;
//
//      protected:
//
//        /// Vector of pointers (atomic).  Provides a mapping from process id (0 to vector size), to bufferPtr (from BufferPool)
//        /// using vector of atomic pointers allows atomic swap of each pointer, without array of mutex or atomic_flags.
//        typedef typename std::conditional<PoolLT == bliss::concurrent::LockType::NONE, BufferType*, std::atomic<BufferType*> >::type BufferPtrTypeInternal;
//        std::vector< BufferPtrTypeInternal > buffers;
//
//        /// for synchronizing access to buffers (and pool).
//        mutable std::mutex mutex;
//        mutable std::atomic_flag spinlock;
//
//      private:
//
//        // private copy contructor
//        SendMessageBuffers(MyType&& other, const std::lock_guard<std::mutex>&) :
//          BaseType(std::move(other)), buffers(std::move(other.buffers))
//          {};
//
//        /**
//         * @brief default copy constructor.  deleted.
//         * @param other   source SendMessageBuffers to copy from.
//         */
//        explicit SendMessageBuffers(const MyType &other) = delete;
//
//        /**
//         * @brief default copy assignment operator, deleted.
//         * @param other   source SendMessageBuffers to copy from.
//         * @return        self
//         */
//        MyType& operator=(const MyType &other) = delete;
//
//
//        /**
//         * @brief Default constructor, deleted
//         */
//        SendMessageBuffers() = delete;
//
//      public:
//        /**
//         * @brief Constructor.
//         * @param numDests         The number of messaging targets/destinations
//         * @param bufferCapacity   The capacity of the individual buffers.  default 8192.
//         * @param poolCapacity     The capacity of the pool.  default unbounded.
//         */
//        explicit SendMessageBuffers(const int & numDests, int numThreads = 0) :
//          BaseType(), buffers(numDests)
//        {
//          this->reset();
//        };
//
//
//
//        /**
//         * default move constructor.  calls superclass move constructor first.
//         * @param other   source SendMessageBuffers to move from.
//         */
//        explicit SendMessageBuffers(MyType && other) :
//            SendMessageBuffers(std::forward<MyType >(other), std::lock_guard<std::mutex>(other.mutex)) {};
//
//
//        /**
//         * @brief default move assignment operator.
//         * @param other   source SendMessageBuffers to move from.
//         * @return        self
//         */
//        MyType& operator=(MyType && other) {
//          std::unique_lock<std::mutex> mine(mutex, std::defer_lock),
//                                          hers(other.mutex, std::defer_lock);
//          std::lock(mine, hers);
//
//          this->pool = std::move(other.pool);
//          buffers.clear(); std::swap(buffers, other.buffers);
//
//          return *this;
//        }
//
//        /**
//         * @brief default destructor
//         */
//        virtual ~SendMessageBuffers() {
//          // when pool destructs, all buffers will be destroyed.
//          buffers.clear();
//        };
//
//        /**
//         * @brief get the number of buffers.  should be same as number of targets for messages
//         * @return    number of buffers
//         */
//        const size_t getSize() const {
//          return buffers.size();
//        }
//
//        /**
//         * @brief get a raw pointer to the the BufferId that a messaging target maps to.
//         * @details for simplicity, not distinguishing between thread safe and unsafe versions
//         *
//         * @param targetRank   the target id for the messages.
//         * @return      reference to the unique_ptr.  when swapping, the content of the unique ptrs are swapped.
//         */
////        template<bliss::concurrent::LockType LT = PoolLT>
////        const typename std::enable_if<LT != bliss::concurrent::LockType::NONE, BufferType*>::type at(const int targetRank) const {
////          return buffers.at(targetRank).load();   // can use reference since we made pool unlimited, so never get nullptr
////        }
//        BufferType* at(const int targetRank) {
//          return (BufferType*)(buffers.at(targetRank));   // if atomic pointer, then will do load()
//        }
//
//        BufferType* operator[](const int targetRank) {
//          return at(targetRank);
//        }
//
//
//        /**
//         * @brief Reset the current MessageBuffers instance by first clearing its list of Buffer Ids, then repopulate it from the pool.
//         * @note  One thread only should call this.
//         */
//        virtual void reset() {
//
//          std::lock_guard<std::mutex> lock(mutex);
//          int count = buffers.size();
//          // release all buffers back to pool
//          for (int i = 0; i < count; ++i) {
//            if (this->at(i)) {
//              this->at(i)->block_and_flush();
//              this->pool.releaseObject(this->at(i));
//              buffers.at(i) = nullptr;
//            }
//          }
//          // reset the pool. local vector should contain a bunch of nullptrs.
//          this->pool.reset();
////          buffers.clear();
////          buffers.resize(count);
//
//          // populate the buffers from the pool
//          for (int i = 0; i < count; ++i) {
//            //printf("initializing buffer %p \n", this->at(i));
//            swapInEmptyBuffer<PoolLT>(i);
//            //printf("initialized buffer %p blocked? %s, empty? %s\n", this->at(i), this->at(i)->is_read_only()? "y" : "n", this->at(i)->isEmpty() ? "y" : "n");
//
//          }
//        }
//
//        /** thread safety:
//         * called by commlayer's send/recv thread, so single threaded there
//         * called by sendMessageBuffer's swapInEmptyBuffer
//         *    if threadsafe version of sendMessageBuffer, then the pointer to buffer to be released will be unique.
//         *    if not threadsafe version, then the sendMessageBuffer should be used by a single threaded program, so buffer to release should be unique.
//         *
//         *    overall, do not need to worry about multiple threads releasing the same buffer.
//         */
//        virtual void releaseBuffer(BufferPtrType ptr) {
//          if (ptr) {
//            ptr->block_and_flush();  // if already blocked and flushed, no op.
//
//            BaseType::releaseBuffer(ptr);
//          }
//        }
//
//        /**
//         * @brief Appends data to the target message buffer.
//         * @details   internally will try to swap in a new buffer when full, and notify the caller of the full buffer's id.
//         * need to return 2 things:
//         *  1. success/failure of current insert
//         *  2. indicator that there is a full buffer (the full buffer's id).
//         *
//         * Notes:
//         *  on left side -
//         *    targetBufferId is obtained outside of CAS, so when CAS is called to update bufferIdForProcId[dest], it may be different already  (on left side)
//         *  on right side -
//         *  fullBufferId and bufferId[dest] are set atomically, but newTargetBufferId is obtained outside of CAS, so the id of the TargetBufferId
//         *  may be more up to date than the CAS function results, which is okay.  Note that a fullBuffer will be marked as blocked to prevent further append.
//         *
//         * Table of values and compare exchange behavior (called when a buffer if full and an empty buffer is swapped in.)
//         *  targetBufferId,     bufferIdForProcId[dest],    newBufferId (from pool) ->  fullBufferId,     bufferId[dest],   newTargetBufferId
//         *  x                   x                   y                           x                 y                 y
//         *  x                   x                   -1                          x                 -1                -1
//         *  x                   z                   y                           -1                z                 z
//         *  x                   z                   -1                          -1                z                 z
//         *  x                   -1                  y                           -1                -1                -1
//         *  x                   -1                  -1                          -1                -1                -1
//         *  -1                  -1                  y                           -1                y                 y
//         *  -1                  -1                  -1                          -1                -1                -1
//         *  -1                  z                   y                           -1                -1                -1
//         *  -1                  z                   -1                          -1                -1                -1
//         *
//         * @note              if failure, data to be inserted would need to be handled by the caller.
//         *
//         * @param[in] data    data to be inserted, as byteArray
//         * @param[in] count   number of bytes to be inserted in the the Buffer
//         * @param[in] dest    messaging target for the data, decides which buffer to append into.
//         * @return            std::pair containing the status of the append (boolean success/fail), and the id of a full buffer if there is one.
//         */
//        std::pair<bool, BufferType*> append(const void* data, const size_t count, const int targetProc) {
//
//          //== if data being set is null, throw error
//          if (data == nullptr || count <= 0) {
//            throw (std::invalid_argument("ERROR: calling MessageBuffer append with nullptr"));
//          }
//
//          //== if there is not enough room for the new data in even a new buffer, LOGIC ERROR IN CODE: throw exception
//          if (count > this->getBufferCapacity()) {
//            throw (std::invalid_argument("ERROR: messageBuffer append with count exceeding Buffer capacity"));
//          }
//
//          //== if targetProc is outside the valid range, throw an error
//          if (targetProc < 0 || targetProc > getSize()) {
//            throw (std::invalid_argument("ERROR: messageBuffer append with invalid targetProc"));
//          }
//
//
//
//          // NOTE: BufferPool is unlimited in size, so don't need to check for nullptr, can just append directly.   is append really atomic?
//          // question is what happens in unique_ptr when dereferencing the internal pointer - if the dereference happens before a unique_ptr swap
//          // from swapInEmptyBuffer, then append would use the old object.  however, it was probably swapped because it's full,
//          // or flushing, so lock_read() would have been set for the full case.  now it depends on what happens in append.
//          // if flush, then block is set just before swap, so it's more likely this thread enters append without a block.
//
//          unsigned int appendResult = 0x0;
//          auto ptr = this->at(targetProc);
//
//          if (ptr) {
////            // DEBUGGING ONLY - for testCommLayer only.  test if call from CommLayer appended the wrong message to MessageBuffers.  Buffer has test for before and after as well.
////            int m = *((int*)data);
////            if ((m / 1000) % 10 == 1) {
////              if (m / 100000 != targetProc + 1) ERRORF("ERROR: DEBUG: MessageBuffers Append wrong input message: %d to proc %d", m, targetProc);
////            }
////            else {
////              if (m % 1000 != targetProc + 1) ERRORF("ERROR: DEBUG: MessageBuffers Append wrong input message: %d to proc %d", m, targetProc);
////            }
//            void * result = nullptr;
//
//            appendResult = ptr->append(data, count, result);  //DEBUGGING FORM
//
////            // DEBUGGING ONLY - for testCommLayer only.  test if call from CommLayer appended the wrong message to MessageBuffers.  Buffer has test for before and after as well.
////            if (result != nullptr) {
////              m = *((int*)result);
////              if ((m / 1000) % 10 == 1) {
////                if (m / 100000 != targetProc + 1) ERRORF("ERROR: DEBUG: MessageBuffers Append wrong outputmessage: %d to proc %d", m, targetProc);
////              }
////              else {
////                if (m % 1000 != targetProc + 1) ERRORF("ERROR: DEBUG: MessageBuffers Append wrong output message: %d to proc %d", m, targetProc);
////              }
////            } else {
////              if (appendResult & 0x1) {
////                ERRORF("ERROR: successful append but result ptr is null!");
////              }
////            }
//          }
//          else {
//            ERRORF("ERROR: Append: shared Buffer ptr is null and no way to swap in a different one.");
//            throw std::logic_error("ERROR: Append: Shared Buffer ptr is null and no way to swap in a different one.");
//          }
//
//
////          unsigned int appendResult = this->at(targetProc)->append(data, count);
//
////          printf("buffer blocked? %s, empty? %s\n", this->at(targetProc)->is_read_only()? "y" : "n", this->at(targetProc)->isEmpty() ? "y" : "n");
//
//          // now if appendResult is false, then we return false, but also swap in a new buffer.
//          // conditions are either full buffer, or blocked buffer.
//          // this call will mark the fullBuffer as blocked to prevent multiple write attempts while flushing the buffer
//          // ideally, we want it to be a reference that gets updated for all threads.
//
//          if (appendResult & 0x2) { // only 1 thread gets 0x2 result for a buffer.  all other threads either writes successfully or fails.
//
//            return std::move(std::make_pair(appendResult & 0x1, swapInEmptyBuffer<PoolLT>(targetProc) ) );
//          } else {
//            //if (preswap != postswap) printf("ERROR: NOSWAP: preswap %p and postswap %p are not the same.\n", preswap, postswap);
//            return std::move(std::make_pair(appendResult & 0x1, nullptr));
//          }
//
//        }
//
//        /**
//         * flushes the buffer at rank targetProc.
//         *
//         * @note  only 1 thread should call this per proc.
//         *        multiple calls to this may be needed to get all data that are written to the new empty buffer.
//         *
//         * @param targetProc
//         * @return
//         */
//        std::vector<BufferType*> flushBufferForRank(const int targetProc) {
//          //== if targetProc is outside the valid range, throw an error
//          if (targetProc < 0 || targetProc > getSize()) {
//            throw (std::invalid_argument("ERROR: messageBuffer append with invalid targetProc"));
//          }
//
//          std::vector<BufferType*> result;
////          // block the old buffer (when append full, will block as well).  each proc has a unique buffer.
////          this->at(targetProc)->block_and_flush();  // this part can be concurrent.
////
//          // passing in getBufferIdForRank result to ensure atomicity.
//          // may return ABSENT if there is no available buffer to use.
//          auto old = swapInEmptyBuffer<PoolLT>(targetProc);
//          if (old) {
//            old->block_and_flush();
//            result.push_back(old);
//          }
//          return result;
//        }
//
//        /**
//         * flushes the buffer at rank targetProc.
//         *
//         * @note  only 1 thread should call this per proc.
//         *        multiple calls to this may be needed to get all data that are written to the new empty buffer.
//         *
//         * @param targetProc
//         * @return
//         */
//        BufferType* threadFlushBufferForRank(const int targetProc, int thread_id = 0) {
//          //== if targetProc is outside the valid range, throw an error
//          if (targetProc < 0 || targetProc > getSize()) {
//            throw (std::invalid_argument("ERROR: messageBuffer append with invalid targetProc"));
//          }
//
////          // block the old buffer (when append full, will block as well).  each proc has a unique buffer.
////          this->at(targetProc)->block_and_flush();  // this part can be concurrent.
////
//          // passing in getBufferIdForRank result to ensure atomicity.
//          // may return ABSENT if there is no available buffer to use.
//          auto old = swapInEmptyBuffer<PoolLT>(targetProc);
//          if (old) old->block_and_flush();
//          return old;
//        }
//
//
//
//      protected:
//
//        /**
//         * @brief  Swap in an empty Buffer from BufferPool at the dest location in MessageBuffers.  The old buffer is returned.  allows swapping in a nullptr if pool is full.
//         * @details The new buffer may be BufferPoolType::ABSENT (invalid buffer, when there is no available Buffer from pool)
//         *
//         * Effect:  bufferIdForProcId[dest] gets a new Buffer Id, full Buffer is set to "blocked" to prevent further append.
//         *
//         * This is the THREAD SAFE version.   Uses std::compare_exchange_strong to ensure that another thread hasn't swapped in a new Empty one already.
//         *
//         * @note we make the caller get the new bufferIdForProcId[dest] directly instead of returning by reference in the method parameter variable
//         *  because this way there is less of a chance of race condition if a shared "old" variable was accidentally used.
//         *    and the caller will likely prefer the most up to date bufferIdForProcId[dest] anyway.
//         *
//         * @note  below is now relevant any more now that message buffer hold buffer ptrs and pool is unlimited.. we only have line 5 here.
//         * 	old is 		old assignment	action
//         *    ABSENT      available     not a valid case
//         *    ABSENT      at dest       swap, return ABS
//         *    ABSENT      not at dest   dest is already swapped.  no op, return ABS
//         *    not ABS     available     not used.  block old, no swap. return ABS
//         *    not ABS     at dest       swap, block old.  return old
//         *    not ABS     not at dest   being used.  no op. return ABS
//         *
//         *
//         * @note:  exchanges buffer at dest atomically.
//         * @note:  because only thread that just fills a buffer will do exchange, only 1 thread does exchange at a time.  can rely on pool's Lock to get one new buffer
//         *          can be assured that each acquired new buffer will be  used, and can ascertain that no 2 threads will replace the same full buffer.
//         *
//         * @tparam        used to choose thread safe vs not verfsion of the method.
//         * @param dest    position in bufferIdForProcId array to swap out
//         * @return        the BufferId that was swapped out.
//         */
//        template<bliss::concurrent::LockType LT = PoolLT>
//        typename std::enable_if<LT != bliss::concurrent::LockType::NONE, BufferType*>::type swapInEmptyBuffer(const int dest) {
//
//          auto ptr = this->pool.tryAcquireObject();
//          int i = 1;
//          while (!ptr) {
//            _mm_pause();
//            ptr = this->pool.tryAcquireObject();
//            ++i;
//          }
//          if (i > 200) WARNINGF("NOTICE: Concurrent Pool shared Buffer ptr took %d iterations to acquire, %d threads.", i, omp_get_num_threads());
//
//          if (ptr) {
//            ptr->clear_and_unblock_writes();
//            memset(ptr->operator int*(), 0, ptr->getCapacity());
//          }
//
//          auto oldbuf = buffers.at(dest).exchange(ptr);
//          if (oldbuf && oldbuf->isEmpty()) {
////            printf("oldbuf %p, blocked? %s\n", oldbuf, oldbuf->is_read_only() ? "y" : "n");
//            releaseBuffer(oldbuf);
//            return nullptr;
//          }
//          else return oldbuf;
//        }
//
//
//        /**
//         * @brief Swap in an empty Buffer from BufferPool at the dest location in MessageBuffers.  The old buffer is returned.
//         * @details The new buffer may be BufferPoolType::ABSENT (invalid buffer, when there is no available Buffer from pool)
//         *
//         * effect:  bufferIdForProcId[dest] gets a new Buffer Id, full Buffer is set to "blocked" to prevent further append.
//         *
//         * this is the THREAD UNSAFE version
//
//         * @note we make the caller get the new bufferIdForProcId[dest] directly instead of returning by reference in the method parameter variable
//         *  because this way there is less of a chance of race condition if a shared "old" variable was accidentally used.
//         *    and the caller will likely prefer the most up to date bufferIdForProcId[dest] anyway.
//         *
//         * @note
//         *    ABSENT      available     not a valid case
//         *    ABSENT      at dest       swap, return ABS
//         *    ABSENT      not at dest   dest is already swapped.  no op, return ABS
//         *    not ABS     available     not used.  block old, no swap. return ABS
//         *    not ABS     at dest       swap, block old.  return old
//         *    not ABS     not at dest   being used.  no op. return ABS
//         *
//         * @note:  requires that the queue at dest to be blocked.
//         *
//         *
//         * @tparam LT     used to choose thread safe vs not verfsion of the method.
//         * @param dest    position in bufferIdForProcId array to swap out
//         * @return        the BufferId that was swapped out.
//         */
//        template<bliss::concurrent::LockType LT = PoolLT>
//        typename std::enable_if<LT == bliss::concurrent::LockType::NONE, BufferType*>::type swapInEmptyBuffer(const int dest) {
//
//          auto oldbuf = this->at(dest);
//          auto ptr = this->pool.tryAcquireObject();
//          int i = 1;
//          while (!ptr) {
//            _mm_pause();
//            ptr = this->pool.tryAcquireObject();
//            ++i;
//          }
//          if (i > 200) WARNINGF("NOTICE: non Concurrent Pool shared Buffer ptr took %d iterations to acquire, %d threads.", i, omp_get_num_threads());
//
//          if (ptr) {
//            ptr->clear_and_unblock_writes();
//            memset(ptr->operator int*(), 0, ptr->getCapacity());
//          }
//
//          buffers.at(dest) = ptr;
//          if (oldbuf && oldbuf->isEmpty()) {
////            printf("oldbuf %p, empty? %s\n", oldbuf, oldbuf->isEmpty() ? "y" : "n");
//
//            releaseBuffer(oldbuf);
//            return nullptr;
//          }
//          else return oldbuf;  // swap the pointer to Buffer object, not Buffer's internal "data" pointer
//        }
//
//
//    };



    /**
     * @class SendMessageBuffers
     * @brief SendMessageBuffers is a subclass of MessageBuffers designed to manage the actual buffering of data for a set of messaging targets.
     * @details
     *
     *  The class is designed with the following requirements in mind:
     *  1. data is destined for some remote process identifiable by a process id (e.g. rank)
     *  2. data is appended to buffer incrementally and with minimal blocking
     *  3. calling thread has a mechanism to process a full buffer. e.g. send it
     *  4. all data in a SendMessageBuffers are homogeneously typed.
     *
     *  The class is implemented to support the requirement:
     *  1. SendMessageBuffers is not aware of its data's type or metadata specifying its type.
     *  2. SendMessageBuffers uses the base MessageBuffers class' internal BufferPool to reuse memory
     *  3. SendMessageBuffers stores a vector of Buffer Ids (as assigned from BufferPool), thus mapping from process Rank to BufferId.
     *  4. SendMessageBuffers provides an Append function to incrementally add data to the Buffer object for the targt process rank.
     *  5. When a particular buffer is full, return the Buffer Id to the calling thread for it to process (send), and swap in an empty
     *      Buffer from BufferPool.
     *
     *  An internal vector is used to map from target process rank to Buffer pointer.  For multithreading flavors of message buffers,
     *    atomic<BufferType*> is used.
     *
     *  @note
     *  pool size is unlimited so we would not get nullptr, and therefore do not need to check it.
     *
     * @note
     *  For Buffers that are full, the calling thread will need to track the full buffers' Ids as this class evicts a full Buffer's id from it's
     *    vector of Buffer Ids.  The Calling Thread can do t his via a queue, as in the case of CommunicationLayer.
     *  Note that a Buffer is "in use" from the first call to the Buffer's append function to the completion of the send (e.g. via MPI Isend)
     *
     *  @tparam LockType    determines if this class should be thread safe
     */
    template<bliss::concurrent::LockType PoolLT, int64_t BufferCapacity>
    class SendMessageBuffers<PoolLT, bliss::concurrent::LockType::NONE, BufferCapacity> :
    public MessageBuffers<PoolLT, bliss::concurrent::LockType::NONE, BufferCapacity>
    {
      protected:
        using BaseType = MessageBuffers<PoolLT, bliss::concurrent::LockType::NONE, BufferCapacity>;
        using MyType = SendMessageBuffers<PoolLT, bliss::concurrent::LockType::NONE, BufferCapacity>;

      public:
        /// Id type of the Buffers
        using BufferType = typename BaseType::BufferType;
        using BufferPtrType = typename BaseType::BufferPtrType;

        /// the BufferPoolType from parent class
        using BufferPoolType = typename BaseType::BufferPoolType;

      protected:

        /// Vector of pointers (atomic).  Provides a mapping from process id (0 to vector size), to bufferPtr (from BufferPool)
        /// using vector of atomic pointers allows atomic swap of each pointer, without array of mutex or atomic_flags.
        typedef typename std::conditional<PoolLT == bliss::concurrent::LockType::NONE, BufferType*,
            std::atomic<BufferType*> >::type BufferPtrTypeInternal;

        /// buffers, one set per target, and in each set, as many as there are local threads.
        std::vector< std::vector< BufferPtrTypeInternal > > buffers;

        /// for synchronizing access to buffers (and pool).
        mutable std::mutex mutex;
        mutable std::atomic_flag spinlock;

      private:

        // private copy contructor
        SendMessageBuffers(SendMessageBuffers&& other, const std::lock_guard<std::mutex>&) :
          BaseType(std::move(other)), buffers(std::move(other.buffers))
          {};

        /**
         * @brief default copy constructor.  deleted.
         * @param other   source SendMessageBuffers to copy from.
         */
        explicit SendMessageBuffers(const MyType &other) = delete;

        /**
         * @brief default copy assignment operator, deleted.
         * @param other   source SendMessageBuffers to copy from.
         * @return        self
         */
        MyType& operator=(const MyType &other) = delete;


        /**
         * @brief Default constructor, deleted
         */
        SendMessageBuffers() = delete;

      public:
        /**
         * @brief Constructor.
         * @param numDests         The number of messaging targets/destinations
         * @param bufferCapacity   The capacity of the individual buffers.  default 8192.
         * @param poolCapacity     The capacity of the pool.  default unbounded.
         */
        explicit SendMessageBuffers(const int & numDests, int numThreads = 0) :
          BaseType(), buffers()
        {
          int nThreads = (numThreads == 0 ? omp_get_max_threads() : numThreads);

          for (int i = 0; i < numDests; ++i) {
            buffers.push_back(std::vector<BufferPtrTypeInternal>(nThreads));
          }

          this->reset();
        };



        /**
         * default move constructor.  calls superclass move constructor first.
         * @param other   source SendMessageBuffers to move from.
         */
        explicit SendMessageBuffers(MyType && other) :
            SendMessageBuffers(std::forward<MyType >(other), std::lock_guard<std::mutex>(other.mutex)) {};


        /**
         * @brief default move assignment operator.
         * @param other   source SendMessageBuffers to move from.
         * @return        self
         */
        MyType& operator=(MyType && other) {
          std::unique_lock<std::mutex> mine(mutex, std::defer_lock),
                                          hers(other.mutex, std::defer_lock);
          std::lock(mine, hers);

          this->pool = std::move(other.pool);
          buffers.clear(); buffers = std::move(other.buffers);

          return *this;
        }

        /**
         * @brief default destructor
         */
        virtual ~SendMessageBuffers() {
          // when pool destructs, all buffers will be destroyed.

          buffers.clear();
        };

        /**
         * @brief get the number of buffers.  should be same as number of targets for messages
         * @return    number of buffers
         */
        const size_t getSize() const {
          return buffers.size();
        }

        /**
         * @brief get a raw pointer to the the BufferId that a messaging target maps to (FOR THE CALLING THREAD)
         * @details for simplicity, not distinguishing between thread safe and unsafe versions
         *
         * @param targetRank   the target id for the messages.
         * @return      reference to the unique_ptr.  when swapping, the content of the unique ptrs are swapped.
         */
//        template<bliss::concurrent::LockType LT = PoolLT>
//        const typename std::enable_if<LT != bliss::concurrent::LockType::NONE, BufferType*>::type at(const int targetRank) const {
//          return buffers.at(targetRank).load();   // can use reference since we made pool unlimited, so never get nullptr
//        }
        BufferType* at(const int targetRank, int srcThread = -1) {
          // openmp - thread ids are sequential from 0 to nThreads - 1, so can use directly as vector index.
          int tid = (srcThread == -1 ? omp_get_thread_num() : srcThread);
          // insert a new one if none there.   note that the entry is initialize to nullptr

          return (BufferType*)(buffers.at(targetRank).at(tid));   // if atomic pointer, then will do load()
        }


        /**
         * @brief get a raw pointer to the the BufferId that a messaging target maps to (FOR THE CALLING THREAD)
         */
        BufferType* operator[](const int targetRank) {
          return at(targetRank);
        }


        /**
         * @brief Reset the current MessageBuffers instance by first clearing its list of Buffer Ids, then repopulate it from the pool.
         * @note  One thread only should call this.
         */
        virtual void reset() {

          std::lock_guard<std::mutex> lock(mutex);
          int count = buffers.size();
          BufferType* ptr;
          int tid, nthreads;
          // release all buffers back to pool
          for (int i = 0; i < count; ++i) {

            nthreads = buffers.at(i).size();

            for (tid = 0; tid < nthreads; ++tid) {
              // get pointer to it (not using iter's copy of it, avoid changing the iterator)
              ptr = this->at(i, tid);
              if (ptr) {
                ptr->block_and_flush();
                this->pool.releaseObject(ptr);
                buffers.at(i).at(tid) = nullptr;
              }
            }
          }
          // reset the pool. local vector should contain a bunch of nullptrs.
          this->pool.reset();

          // populate the buffers from the pool
          for (int i = 0; i < count; ++i) {
            nthreads = buffers.at(i).size();

            for (tid = 0; tid < nthreads; ++tid) {
              swapInEmptyBuffer<PoolLT>(i, tid);
            }
          }
        }

        /** thread safety:
         * called by commlayer's send/recv thread, so single threaded there
         * called by sendMessageBuffer's swapInEmptyBuffer
         *    if threadsafe version of sendMessageBuffer, then the pointer to buffer to be released will be unique.
         *    if not threadsafe version, then the sendMessageBuffer should be used by a single threaded program, so buffer to release should be unique.
         *
         *    overall, do not need to worry about multiple threads releasing the same buffer.
         */
        virtual void releaseBuffer(BufferPtrType ptr) {
          if (ptr) {
            ptr->block_and_flush();  // if already blocked and flushed, no op.

            BaseType::releaseBuffer(ptr);
          }
        }

        /**
         * @brief Appends data to the target message buffer.
         * @details   internally will try to swap in a new buffer when full, and notify the caller of the full buffer's id.
         * need to return 2 things:
         *  1. success/failure of current insert
         *  2. indicator that there is a full buffer (the full buffer's id).
         *
         * Notes:
         *  on left side -
         *    targetBufferId is obtained outside of CAS, so when CAS is called to update bufferIdForProcId[dest], it may be different already  (on left side)
         *  on right side -
         *  fullBufferId and bufferId[dest] are set atomically, but newTargetBufferId is obtained outside of CAS, so the id of the TargetBufferId
         *  may be more up to date than the CAS function results, which is okay.  Note that a fullBuffer will be marked as blocked to prevent further append.
         *
         * Table of values and compare exchange behavior (called when a buffer if full and an empty buffer is swapped in.)
         *  targetBufferId,     bufferIdForProcId[dest],    newBufferId (from pool) ->  fullBufferId,     bufferId[dest],   newTargetBufferId
         *  x                   x                   y                           x                 y                 y
         *  x                   x                   -1                          x                 -1                -1
         *  x                   z                   y                           -1                z                 z
         *  x                   z                   -1                          -1                z                 z
         *  x                   -1                  y                           -1                -1                -1
         *  x                   -1                  -1                          -1                -1                -1
         *  -1                  -1                  y                           -1                y                 y
         *  -1                  -1                  -1                          -1                -1                -1
         *  -1                  z                   y                           -1                -1                -1
         *  -1                  z                   -1                          -1                -1                -1
         *
         * @note              if failure, data to be inserted would need to be handled by the caller.
         *
         * @param[in] data    data to be inserted, as byteArray
         * @param[in] count   number of bytes to be inserted in the the Buffer
         * @param[in] dest    messaging target for the data, decides which buffer to append into.
         * @return            std::pair containing the status of the append (boolean success/fail), and the id of a full buffer if there is one.
         */
        std::pair<bool, BufferType*> append(const void* data, const size_t count, const int targetProc, int thread_id = -1) {

          //== if data being set is null, throw error
          if (data == nullptr || count <= 0) {
            throw (std::invalid_argument("ERROR: calling MessageBuffer append with nullptr"));
          }

          //== if there is not enough room for the new data in even a new buffer, LOGIC ERROR IN CODE: throw exception
          if (count > this->getBufferCapacity()) {
            throw (std::invalid_argument("ERROR: messageBuffer append with count exceeding Buffer capacity"));
          }

          //== if targetProc is outside the valid range, throw an error
          if (targetProc < 0 || targetProc > getSize()) {
            throw (std::invalid_argument("ERROR: messageBuffer append with invalid targetProc"));
          }

          int tid = (thread_id < 0 ? omp_get_thread_num() : thread_id);


          // NOTE: BufferPool is unlimited in size, so don't need to check for nullptr, can just append directly.   is append really atomic?
          // question is what happens in unique_ptr when dereferencing the internal pointer - if the dereference happens before a unique_ptr swap
          // from swapInEmptyBuffer, then append would use the old object.  however, it was probably swapped because it's full,
          // or flushing, so lock_read() would have been set for the full case.  now it depends on what happens in append.
          // if flush, then block is set just before swap, so it's more likely this thread enters append without a block.

          unsigned int appendResult = 0x0;
          BufferType* ptr = this->at(targetProc, thread_id);

          if (ptr) {
//            // DEBUGGING ONLY - for testCommLayer only.  test if call from CommLayer appended the wrong message to MessageBuffers.  Buffer has test for before and after as well.
//            int m = *((int*)data);
//            if ((m / 1000) % 10 == 1) {
//              if (m / 100000 != targetProc + 1) ERRORF("ERROR: DEBUG: MessageBuffers Append wrong input message: %d to proc %d", m, targetProc);
//            }
//            else {
//              if (m % 1000 != targetProc + 1) ERRORF("ERROR: DEBUG: MessageBuffers Append wrong input message: %d to proc %d", m, targetProc);
//            }
            void * result = nullptr;

            appendResult = ptr->append(data, count, result);  //DEBUGGING FORM

//            // DEBUGGING ONLY - for testCommLayer only.  test if call from CommLayer appended the wrong message to MessageBuffers.  Buffer has test for before and after as well.
//            if (result != nullptr) {
//              m = *((int*)result);
//              if ((m / 1000) % 10 == 1) {
//                if (m / 100000 != targetProc + 1) ERRORF("ERROR: DEBUG: MessageBuffers Append wrong output message: %d to proc %d", m, targetProc);
//              }
//              else {
//                if (m % 1000 != targetProc + 1) ERRORF("ERROR: DEBUG: MessageBuffers Append wrong output message: %d to proc %d", m, targetProc);
//              }
//            } else {
//              if (appendResult & 0x1) {
//                ERRORF("ERROR: successful append but result ptr is null!");
//              }
//            }
          }
          else {
            // do NOT swap in one here, since access here is concurrent, multiple threads may try to swap.
            ERRORF("ERROR: Append: threadlocal Buffer ptr is null and no way to swap in a different one.");
            throw std::logic_error("ERROR: Append: threadlocal Buffer ptr is null and no way to swap in a different one.");
          }

//          printf("buffer blocked? %s, empty? %s\n", this->at(targetProc)->is_read_only()? "y" : "n", this->at(targetProc)->isEmpty() ? "y" : "n");

          // now if appendResult is false, then we return false, but also swap in a new buffer.
          // conditions are either full buffer, or blocked buffer.
          // this call will mark the fullBuffer as blocked to prevent multiple write attempts while flushing the buffer
          // ideally, we want it to be a reference that gets updated for all threads.

          if (appendResult & 0x2) { // only 1 thread gets 0x2 result for a buffer.  all other threads either writes successfully or fails.

            return std::move(std::make_pair(appendResult & 0x1, swapInEmptyBuffer<PoolLT>(targetProc, tid) ) );

          } else {
            //if (preswap != postswap) printf("ERROR: NOSWAP: preswap %p and postswap %p are not the same.\n", preswap, postswap);
            return std::move(std::make_pair(appendResult & 0x1, nullptr));
          }

        }

        /**
         * flushes the buffer at rank targetProc.
         *
         * @note  only 1 thread should call this per proc.
         *        multiple calls to this may be needed to get all data that are written to the new empty buffer.
         *
         * @param targetProc
         * @return
         */
        std::vector<BufferType*> flushBufferForRank(const int targetProc) {
          //== if targetProc is outside the valid range, throw an error
          if (targetProc < 0 || targetProc > getSize()) {
            throw (std::invalid_argument("ERROR: messageBuffer append with invalid targetProc"));
          }

          BufferType* old;
          std::vector<BufferType*> result;
          int tid = omp_get_thread_num();
          int nthreads = buffers.at(targetProc).size();
          for (int t = 0; t < nthreads; ++t) {
            old = swapInEmptyBuffer<PoolLT>(targetProc, (t + tid) % nthreads);
            if (old) {
              old->block_and_flush();
              result.push_back(old);
            }
          }

          return result;
        }

        BufferType* threadFlushBufferForRank(const int targetProc, int thread_id = -1) {
          //== if targetProc is outside the valid range, throw an error
          if (targetProc < 0 || targetProc > getSize()) {
            throw (std::invalid_argument("ERROR: messageBuffer append with invalid targetProc"));
          }

          int tid = (thread_id < 0 ? omp_get_thread_num() : thread_id);

          BufferType* old = swapInEmptyBuffer<PoolLT>(targetProc, tid);
          if (old) {
            old->block_and_flush();
            return old;
          }
        }

      protected:

        /**
         * @brief  Swap in an empty Buffer from BufferPool at the dest location in MessageBuffers.  The old buffer is returned.
         * @details The new buffer may be BufferPoolType::ABSENT (invalid buffer, when there is no available Buffer from pool)
         *
         * Effect:  bufferIdForProcId[dest] gets a new Buffer Id, full Buffer is set to "blocked" to prevent further append.
         *
         * This is the THREAD SAFE version.   Uses std::compare_exchange_strong to ensure that another thread hasn't swapped in a new Empty one already.
         *
         * @note we make the caller get the new bufferIdForProcId[dest] directly instead of returning by reference in the method parameter variable
         *  because this way there is less of a chance of race condition if a shared "old" variable was accidentally used.
         *    and the caller will likely prefer the most up to date bufferIdForProcId[dest] anyway.
         *
         * @note  below is now relevant any more now that message buffer hold buffer ptrs and pool is unlimited.. we only have line 5 here.
         *  old is    old assignment  action
         *    ABSENT      available     not a valid case
         *    ABSENT      at dest       swap, return ABS
         *    ABSENT      not at dest   dest is already swapped.  no op, return ABS
         *    not ABS     available     not used.  block old, no swap. return ABS
         *    not ABS     at dest       swap, block old.  return old
         *    not ABS     not at dest   being used.  no op. return ABS
         *
         *
         * @note:  exchanges buffer at dest atomically.
         * @note:  because only thread that just fills a buffer will do exchange, only 1 thread does exchange at a time.  can rely on pool's Lock to get one new buffer
         *          can be assured that each acquired new buffer will be  used, and can ascertain that no 2 threads will replace the same full buffer.
         *
         * @tparam        used to choose thread safe vs not verfsion of the method.
         * @param dest    position in bufferIdForProcId array to swap out
         * @return        the BufferId that was swapped out.
         */
        template<bliss::concurrent::LockType LT = PoolLT>
        typename std::enable_if<LT != bliss::concurrent::LockType::NONE, BufferType*>::type swapInEmptyBuffer(const int dest, int thread_id = -1) {

          int tid = thread_id < 0 ? omp_get_thread_num() : thread_id;

          auto ptr = this->pool.tryAcquireObject();
          int i = 1;
          while (!ptr) {
            _mm_pause();
            ptr = this->pool.tryAcquireObject();
            ++i;
          }
          if (i > 200) WARNINGF("NOTICE: Concurrent Pool threadlocal Buffer ptr took %d iterations to acquire, %d threads.", i, omp_get_num_threads());

          if (ptr) ptr->clear_and_unblock_writes();


          BufferType* oldbuf = buffers.at(dest).at(tid).exchange(ptr);
          //printf("INFO: swap old %p, new %p, for thread %d: target rank %d of %lu, curr size %lu\n", oldbuf, ptr, tid, dest, buffers.size(), buffers.at(dest).size());

          if (oldbuf && oldbuf->isEmpty()) {
//            printf("oldbuf %p, blocked? %s\n", oldbuf, oldbuf->is_read_only() ? "y" : "n");
            releaseBuffer(oldbuf);
            oldbuf = nullptr;
          }
          return oldbuf;
        }


        /**
         * @brief Swap in an empty Buffer from BufferPool at the dest location in MessageBuffers.  The old buffer is returned.
         * @details The new buffer may be BufferPoolType::ABSENT (invalid buffer, when there is no available Buffer from pool)
         *
         * effect:  bufferIdForProcId[dest] gets a new Buffer Id, full Buffer is set to "blocked" to prevent further append.
         *
         * this is the THREAD UNSAFE version

         * @note we make the caller get the new bufferIdForProcId[dest] directly instead of returning by reference in the method parameter variable
         *  because this way there is less of a chance of race condition if a shared "old" variable was accidentally used.
         *    and the caller will likely prefer the most up to date bufferIdForProcId[dest] anyway.
         *
         * @note
         *    ABSENT      available     not a valid case
         *    ABSENT      at dest       swap, return ABS
         *    ABSENT      not at dest   dest is already swapped.  no op, return ABS
         *    not ABS     available     not used.  block old, no swap. return ABS
         *    not ABS     at dest       swap, block old.  return old
         *    not ABS     not at dest   being used.  no op. return ABS
         *
         * @note:  requires that the queue at dest to be blocked.
         *
         *
         * @tparam LT     used to choose thread safe vs not verfsion of the method.
         * @param dest    position in bufferIdForProcId array to swap out
         * @return        the BufferId that was swapped out.
         */
        template<bliss::concurrent::LockType LT = PoolLT>
        typename std::enable_if<LT == bliss::concurrent::LockType::NONE, BufferType*>::type swapInEmptyBuffer(const int dest, int thread_id = -1) {

          int tid = thread_id < 0 ? omp_get_thread_num() : thread_id;

          BufferType* oldbuf = buffers.at(dest).at(tid);
          BufferType* ptr = this->pool.tryAcquireObject();
          int i = 1;
          while (!ptr) {
            _mm_pause();
            ptr = this->pool.tryAcquireObject();
            ++i;
          }
          if (i > 200) WARNINGF("NOTICE: non-Concurrent Pool threadlocal Buffer ptr took %d iterations to acquire, %d threads.", i, omp_get_num_threads());


          if (ptr) {
            ptr->clear_and_unblock_writes();
            //memset(ptr->operator int*(), 0, ptr->getCapacity());
          }

          buffers.at(dest).at(tid) = ptr;
          if (oldbuf && oldbuf->isEmpty()) {
            releaseBuffer(oldbuf);
            oldbuf = nullptr;
          }
          return oldbuf;
        }


    };

  } /* namespace io */
} /* namespace bliss */

#endif /* MESSAGEBUFFERS_HPP_ */
