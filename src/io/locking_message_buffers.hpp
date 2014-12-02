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
#include "io/locking_buffer.hpp"
#include "io/locking_object_pool.hpp"
#include <atomic>

#include <iterator> // for ostream_iterator

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
        virtual void reset() {
          pool.reset();
        }
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
    class SendMessageBuffers : public MessageBuffers<PoolLT, BufferLT, BufferCapacity>
    {
      public:
        /// Id type of the Buffers
        using BufferType = typename MessageBuffers<PoolLT, BufferLT, BufferCapacity>::BufferType;
        //using BufferPtrType = typename MessageBuffers<PoolLT, BufferLT, BufferCapacity>::BufferPtrType;

        /// the BufferPoolType from parent class
        using BufferPoolType = typename MessageBuffers<PoolLT, BufferLT, BufferCapacity>::BufferPoolType;

      protected:

        /// Vector of pointers (atomic).  Provides a mapping from process id (0 to vector size), to bufferPtr (from BufferPool)
        /// using vector of atomic pointers allows atomic swap of each pointer, without array of mutex or atomic_flags.
        typedef typename std::conditional<PoolLT == bliss::concurrent::LockType::NONE, BufferType*, std::atomic<BufferType*> >::type BufferPtrTypeInternal;
        std::vector< BufferPtrTypeInternal > buffers;

        /// for synchronizing access to buffers (and pool).
        mutable std::mutex mutex;
        mutable std::atomic_flag spinlock;

      private:

        // private copy contructor
        SendMessageBuffers(SendMessageBuffers&& other, const std::lock_guard<std::mutex>&) :
          MessageBuffers<PoolLT, BufferLT, BufferCapacity>(std::move(other)), buffers(std::move(other.buffers))
          {};

        /**
         * @brief default copy constructor.  deleted.
         * @param other   source SendMessageBuffers to copy from.
         */
        explicit SendMessageBuffers(const SendMessageBuffers<PoolLT, BufferLT, BufferCapacity> &other) = delete;

        /**
         * @brief default copy assignment operator, deleted.
         * @param other   source SendMessageBuffers to copy from.
         * @return        self
         */
        SendMessageBuffers<PoolLT, BufferLT, BufferCapacity>& operator=(const SendMessageBuffers<PoolLT, BufferLT, BufferCapacity> &other) = delete;


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
        explicit SendMessageBuffers(const int & numDests) :
          MessageBuffers<PoolLT, BufferLT, BufferCapacity>(), buffers(numDests)
        {
          this->reset();
        };



        /**
         * default move constructor.  calls superclass move constructor first.
         * @param other   source SendMessageBuffers to move from.
         */
        explicit SendMessageBuffers(SendMessageBuffers<PoolLT, BufferLT, BufferCapacity> && other) :
            SendMessageBuffers(std::forward<SendMessageBuffers<PoolLT, BufferLT, BufferCapacity> >(other), std::lock_guard<std::mutex>(other.mutex)) {};


        /**
         * @brief default move assignment operator.
         * @param other   source SendMessageBuffers to move from.
         * @return        self
         */
        SendMessageBuffers<PoolLT, BufferLT, BufferCapacity>& operator=(SendMessageBuffers<PoolLT, BufferLT, BufferCapacity> && other) {
          std::unique_lock<std::mutex> mine(mutex, std::defer_lock),
                                          hers(other.mutex, std::defer_lock);
          std::lock(mine, hers);

          this->pool = std::move(other.pool);
          buffers.clear(); std::swap(buffers, other.buffers);

          return *this;
        }

        /**
         * @brief default destructor
         */
        virtual ~SendMessageBuffers() {
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
         * @brief get a raw pointer to the the BufferId that a messaging target maps to.
         * @details for simplicity, not distinguishing between thread safe and unsafe versions
         *
         * @param targetRank   the target id for the messages.
         * @return      reference to the unique_ptr.  when swapping, the content of the unique ptrs are swapped.
         */
//        template<bliss::concurrent::LockType LT = PoolLT>
//        const typename std::enable_if<LT != bliss::concurrent::LockType::NONE, BufferType*>::type at(const int targetRank) const {
//          return buffers.at(targetRank).load();   // can use reference since we made pool unlimited, so never get nullptr
//        }
        BufferType* at(const int targetRank) const {
          return (BufferType*)(buffers.at(targetRank));   // if atomic pointer, then will do load()
        }

        BufferType* operator[](const int targetRank) const {
          return at(targetRank);
        }


        /**
         * @brief Reset the current MessageBuffers instance by first clearing its list of Buffer Ids, then repopulate it from the pool.
         * @note  One thread only should call this.
         */
        virtual void reset() {

          std::lock_guard<std::mutex> lock(mutex);
          // release all buffers back to pool
          for (int i = 0; i < buffers.size(); ++i) {
            if (this->at(i)) {
              this->at(i)->block_and_flush();
              this->pool.releaseObject(this->at(i));
            }
          }
          // reset the pool. local vector should contain a bunch of nullptrs.
          this->pool.reset();

          // populate the buffers from the pool
          for (int i = 0; i < buffers.size(); ++i) {
            swapInEmptyBuffer<PoolLT>(i);
          }
        }

        virtual void releaseBuffer(BufferType * ptr) {

          ptr->block_and_flush();  // if already blocked and flushed, no op.

          MessageBuffers<PoolLT, BufferLT, BufferCapacity>::releaseBuffer(ptr);
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
        std::pair<bool, BufferType*> append(const void* data, const size_t count, const int targetProc) {

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



          // NOTE: BufferPool is unlimited in size, so don't need to check for nullptr, can just append directly.   is append really atomic?
          // question is what happens in unique_ptr when dereferencing the internal pointer - if the dereference happens before a unique_ptr swap
          // from swapInEmptyBuffer, then append would use the old object.  however, it was probably swapped because it's full,
          // or flushing, so lock_read() would have been set for the full case.  now it depends on what happens in append.
          // if flush, then block is set just before swap, so it's more likely this thread enters append without a block.

          unsigned int appendResult = this->at(targetProc)->append(data, count);


          // now if appendResult is false, then we return false, but also swap in a new buffer.
          // conditions are either full buffer, or blocked buffer.
          // this call will mark the fullBuffer as blocked to prevent multiple write attempts while flushing the buffer
          // ideally, we want it to be a reference that gets updated for all threads.

          if (appendResult & 0x2) { // only 1 thread gets 0x2 result for a buffer.  all other threads either writes successfully or fails.

            return std::move(std::make_pair(appendResult & 0x1, swapInEmptyBuffer<PoolLT>(targetProc) ) );


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
        BufferType* flushBufferForRank(const int targetProc) {
          //== if targetProc is outside the valid range, throw an error
          if (targetProc < 0 || targetProc > getSize()) {
            throw (std::invalid_argument("ERROR: messageBuffer append with invalid targetProc"));
          }

//          // block the old buffer (when append full, will block as well).  each proc has a unique buffer.
//          this->at(targetProc)->block_and_flush();  // this part can be concurrent.
//
          // passing in getBufferIdForRank result to ensure atomicity.
          // may return ABSENT if there is no available buffer to use.
          auto old = swapInEmptyBuffer<PoolLT>(targetProc);
          if (old) old->block_and_flush();
          return old;
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
         * 	old is 		old assignment	action
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
        typename std::enable_if<LT != bliss::concurrent::LockType::NONE, BufferType*>::type swapInEmptyBuffer(const int dest) {

          auto ptr = this->pool.tryAcquireObject();
          ptr->clear_and_unblock_writes();
          auto oldbuf = buffers.at(dest).exchange(ptr);
          if (oldbuf.isEmpty()) {
            releaseBuffer(oldbuf);
            return nullptr;
          }
          else return oldbuf;
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
        typename std::enable_if<LT == bliss::concurrent::LockType::NONE, BufferType*>::type swapInEmptyBuffer(const int dest) {

          auto oldbuf = this->at(dest);
          buffers.at(dest) = this->pool.tryAcquireObject();
          buffers.at(dest)->clear_and_unblock_writes();
          if (oldbuf.isEmpty()) {
            releaseBuffer(oldbuf);
            return nullptr;
          }
          else return oldbuf;  // swap the pointer to Buffer object, not Buffer's internal "data" pointer
        }


    };

  } /* namespace io */
} /* namespace bliss */

#endif /* MESSAGEBUFFERS_HPP_ */
