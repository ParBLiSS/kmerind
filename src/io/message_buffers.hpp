/**
 * @file		MessageBuffers.hpp
 * @ingroup
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

#include <mutex>
#include <vector>
#include <type_traits>
#include <typeinfo>
#include <iostream>
#include <sstream>

#include "config.hpp"
#include "io/buffer.hpp"
#include "io/buffer_pool.hpp"

#include <iterator> // for ostream_iterator

namespace bliss
{
  namespace io
  {

    /**
     * @class			bliss::io::MessageBuffers
     * @brief     a data structure to buffering/batching messages for communication
     * @details   Templated base class for a collection of buffers that are used for communication.
     *
     *            Templated to support thread safe and unsafe operations
     *
     *            The class contains a BufferPool, which contains an array of Buffer objects for in-memory storage.
     *
     */
    template<bliss::concurrent::ThreadSafety ThreadSafety>
    class MessageBuffers
    {
        // TODO: move consturctor and assignment operator to copy between thread safety levels.
      protected:
        /**
         * Internal BufferPool Type.  typedefed only to shorten the usage.
         */
        typedef typename bliss::io::BufferPool<ThreadSafety>  BufferPoolType;

        /**
         * IdType of a Buffer, aliased from BufferPool
         */
        typedef typename BufferPoolType::IdType    BufferIdType;

        /**
         * a pool of in-memory Buffers for storage.
         */
        BufferPoolType pool;

        /**
         * Gets the capacity of the Buffer instances.
         * @return    capacity of each buffer.
         */
        const size_t getBufferCapacity() const {
          return pool.getBufferCapacity();
        }

        /**
         * Protected so base class function is not called directly.
         *
         * Constructor.  creates a buffer given a specified per-buffer capacity, and the BufferPool capacity.
         *
         * @param buffer_capacity   the Buffer's maximum capacity.  default to 8192.
         * @param pool_capacity     the BufferPool's capacity.  default to unlimited.
         */
        explicit MessageBuffers(const size_t & _buffer_capacity = 8192, const BufferIdType & pool_capacity = std::numeric_limits<BufferIdType>::max()) :
          pool(pool_capacity, _buffer_capacity) {

          if (pool_capacity < 0) {
            throw std::invalid_argument("ERROR: pool_capacity set to less than 0");
          }
        };

        /**
         * Protected so base class function is not called directly.
         *
         * default copy constructor.  deleted.  since internal BufferPool does not allow copy construction/assignment.
         * @param other     source MessageBuffers to copy from
         */
        explicit MessageBuffers(const MessageBuffers<ThreadSafety>& other) = delete;

        /**
         * Protected so base class function is not called directly.
         *
         * default copy assignment operator.  deleted.   since internal BufferPool does not allow copy construction/assignment.
         * @param other     source MessageBuffers to copy from
         * @return          self.
         */
        MessageBuffers<ThreadSafety>& operator=(const MessageBuffers<ThreadSafety>& other) = delete;


        /**
         * Protected so base class function is not called directly.
         *
         * default move constructor.
         * @param other     source MessageBuffers to move from
         */
        explicit MessageBuffers(MessageBuffers<ThreadSafety>&& other) : pool(std::move(other.pool)) {};


        /**
         * Protected so base class function is not called directly.
         *
         * default move assignment operator.
         * @param other     source MessageBuffers to move from
         * @return          self.
         */
        MessageBuffers<ThreadSafety>& operator=(MessageBuffers<ThreadSafety>&& other) {
          pool = std::move(other.pool);
          return *this;
        }


      public:
        /**
         * default destructor
         */
        virtual ~MessageBuffers() {};

        /**
         * Accesses the Buffer object inside the BufferPool by BufferPool assigned Id
         * @param id      Id of the Buffer to be accessed.
         * @return        const reference to the Buffer object.
         */
        const Buffer<ThreadSafety> & getBackBuffer(const BufferIdType& id) const {
          return this->pool[id];
        }


        /**
         * Releases a Buffer by it's BufferPool id, after the buffer is no longer needed.
         *
         * This to be called after the communication logic is done with the Send or Recv buffer.
         *
         * @param id    Buffer's id (assigned from BufferPool)
         */
        void releaseBuffer(const BufferIdType &id) {
          pool.releaseBuffer(id);
        }

        /**
         * Convenience method to release and clear all current Buffers in the pool.
         */
        virtual void reset() {
          pool.reset();
        }
    };



    /**
     * SendMessageBuffers is a subclass of MessageBuffers designed to manage the actual buffering of data for a set of messaging targets.
     *
     * The class is designed with the following requirements in mind:
     *  1. data is destined for some remote process identifiable by a process id (e.g. rank)
     *  2. data is appended to buffer incrementally and with minimal blocking
     *  3. calling thread has a mechanism to process a full buffer.
     *  4. data is typed, but all data in a SendMessageBuffers are homogeneously typed.
     *
     *  The class is implemented to support the requirement:
     *  1. SendMessageBuffers is not aware of its data type or metadata specifying its type.
     *  2. SendMessageBuffers uses the base MessageBuffers class' internal BufferPool to reuse memory
     *  3. SendMessageBuffers stores a vector of Buffer Ids (as assigned from BufferPool), thus mapping from process Rank to BufferId.
     *  4. SendMessageBuffers provides an Append function to incrementally add data to the Buffer object for the targt process rank.
     *  5. When a particular buffer is full, return the Buffer Id to the calling thread for it to process (send), and swap in an empty
     *      Buffer from BufferPool.
     *
     *
     *  For Buffers that are full, the calling thread will need to track the Ids as this class evicts a full Buffer's id from it's
     *    vector of Buffer Ids.  This can be done via a queue, as in the case of CommunicationLayer.
     *  Note that a Buffer is "in use" from the first call to the Buffer's append function to the completion of the send (e.g. via MPI Isend)
     *
     */
    template<bliss::concurrent::ThreadSafety ThreadSafety>
    class SendMessageBuffers : public MessageBuffers<ThreadSafety>
    {
      public:
        /**
         * Id type of the Buffers
         */
        using BufferIdType = typename MessageBuffers<ThreadSafety>::BufferIdType;

        /**
         * the BufferPoolType from parent class
         */
        using BufferPoolType = typename MessageBuffers<ThreadSafety>::BufferPoolType;

      protected:
        /**
         * Id type of the Buffers, but here depending on the specified ThreadSafety template parameter,
         * either BufferIdType or atomic version of BufferIdType
         */
        typedef typename std::conditional<ThreadSafety, std::atomic<BufferIdType>, BufferIdType>::type   IdType;

        /**
         * Vector of IdType (atomic or not).  Provides a mapping from process id (0 to vector size), to buffer Ids (from BufferPool)
         */
        std::vector< IdType > bufferIds;

      public:
        /**
         * Constructor.
         * @param numDests         The number of messaging targets/destinations
         * @param bufferCapacity   The capacity of the individual buffers.  default 8192.
         * @param poolCapacity     The capacity of the pool.  default unbounded.
         */
        explicit SendMessageBuffers(const int & numDests, const size_t & bufferCapacity = 8192, const BufferIdType & poolCapacity = std::numeric_limits<BufferIdType>::max()) :
          MessageBuffers<ThreadSafety>(bufferCapacity, poolCapacity), bufferIds(numDests)
        {
          /// initialize the bufferIds, acquiring them from the pool by calling reset.
          this->reset();
        };

        /**
         * Default constructor, deleted
         */
        SendMessageBuffers() :  MessageBuffers<ThreadSafety>() {};

        /**
         * default copy constructor.  deleted.
         * @param other   source SendMessageBuffers to copy from.
         */
        explicit SendMessageBuffers(const SendMessageBuffers<ThreadSafety> &other) = delete;

        /**
         * default move constructor.  calls superclass move constructor first.
         * @param other   source SendMessageBuffers to move from.
         */
        explicit SendMessageBuffers(SendMessageBuffers<ThreadSafety> && other) :
            MessageBuffers<ThreadSafety>(std::move(other)), bufferIds(std::move(other.bufferIds)) {};

        /**
         * default copy assignment operator, deleted.
         * @param other   source SendMessageBuffers to copy from.
         * @return        self
         */
        SendMessageBuffers<ThreadSafety>& operator=(const SendMessageBuffers<ThreadSafety> &other) = delete;


        /**
         * default move assignment operator.
         * @param other   source SendMessageBuffers to move from.
         * @return        self
         */
        SendMessageBuffers<ThreadSafety>& operator=(SendMessageBuffers<ThreadSafety> && other) {
          bufferIds = std::move(other.bufferIds);
          this->pool = std::move(other.pool);
          return *this;
        }

        /**
         * default destructor
         */
        virtual ~SendMessageBuffers() {};

        /**
         * get the number of buffers.  should be same as number of targets for messages
         * @return    number of buffers
         */
        const size_t getSize() const {
          return bufferIds.size();
        }

        /**
         * get the BufferId that a messaging target max to.
         * for simplicity, not distinguishing between thread safe and unsafe versions
         *
         * @param idx   the target id for the messages.
         * @return      the buffer id that the target id maps to.
         */
        const BufferIdType getBufferId(const int idx) const {
          return BufferIdType(bufferIds[idx]);
        }

        /**
         * get the list of active Buffer ids.  (for debugging)
         * @return  vector containing the active Buffer ids.
         */
        const std::vector<IdType>& getBufferIds() const {
          return bufferIds;
        }

        /**
         * convenience method to convert the list of Active Buffer Ids to string.  (for debugging)
         * @return   string with list of active Buffer Ids
         */
        const std::string bufferIdsToString() const {
          std::stringstream ss;
          std::ostream_iterator<IdType> ost(ss, ",");
          std::copy(bufferIds.begin(), bufferIds.end(), ost);
          ss << std::endl;
          return ss.str();
        }

        /**
         * Reset the current MessageBuffers instance by first clearing its list of Buffer Ids, then repopulate it from
         * the pool.
         */
        virtual void reset() {
          MessageBuffers<ThreadSafety>::reset();
          for (int i = 0; i < this->bufferIds.size(); ++i) {
            bufferIds[i] = this->pool.tryAcquireBuffer().second;
          }
        }

        /**
         * Appends data to the target message buffer.   internally will try to swap in a new buffer when full, and notify the caller of the full buffer's id.
         * need to return 2 things:
         *  success/failure of current insert
         *  indicator that there is a full buffer (the full buffer's id).
         *
         *  table of values and compare exchange behavior (called when a buffer if full and an empty buffer is swapped in.)
         *  targetBufferId is obtained outside of CAS, so when CAS is called, bufferIds[dest] may be different already  (on left side)
         *  fullBufferId and bufferId[dest] are set atomically, but newTargetBufferId is obtained outside of CAS, so the id of the TargetBufferId
         *  may be more up to date than the CAS function results, which is okay.  Note that a fullBuffer will be marked as blocked to prevent further append.
         *
         *  targetBufferId,     bufferIds[dest],    newBufferId (from pool) ->  fullBufferId,     bufferId[dest],   newTargetBufferId
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
         *
         * @param[in] data    data to be inserted, as byteArray
         * @param[in] count   number of bytes to be inserted in the the Buffer
         * @param[in] dest    messaging target for the data, decides which buffer to append into.
         * @return            std::pair containing the status of the append (boolean success/fail), and the id of a full buffer if there is one.
         */
        std::pair<bool, BufferIdType> append(const void* data, const size_t &count, const int &targetProc) {
          /// if count is 0, no write and this succeeds right away.
          if (count == 0)
            return std::pair<bool, BufferIdType>(true, BufferPoolType::INVALID);

          /// if there is not enough room for the new data in even a new buffer, LOGIC ERROR IN CODE: throw exception
          if (count > this->getBufferCapacity()) {
            throw (std::invalid_argument("ERROR: messageBuffer append with count exceeding Buffer capacity"));
          }

          /// if targetProc is outside the valid range, throw an error
          if (targetProc < 0 || targetProc > getSize()) {
            throw (std::invalid_argument("ERROR: messageBuffer append with invalid targetProc"));
          }

          /// if data being set is null, throw error
          if (data == nullptr) {
            throw (std::invalid_argument("ERROR: calling MessageBuffer append with nullptr"));
          }

          /// default fullBufferId return value.
          BufferIdType fullBufferId = BufferPoolType::INVALID;

          /// get the current Buffer's Id
          BufferIdType targetBufferId = getBufferId(targetProc);

          /// if targetBufferId is an invalid buffer, or if new data can't fit, need to replace bufferIds[dest]
          if ((targetBufferId == BufferPoolType::INVALID) ||
              this->pool[targetBufferId].getSize() > (this->getBufferCapacity() - count)) {
            // at this point, the local variables may be out of date already, and targetBufferId may already marked as full.

            // this call will mark the fullBuffer as blocked to prevent multiple write attempts while flushing the buffer
            fullBufferId = swapInEmptyBuffer<ThreadSafety>(targetProc, targetBufferId);    // swap in an empty buffer
            targetBufferId = getBufferId(targetProc);               // and get the updated targetBuffer id
          }


          if (targetBufferId == BufferPoolType::INVALID) {
            // we don't have a buffer to insert into
            // then don't insert
            return std::pair<bool, BufferIdType>(false, fullBufferId);
          } else {
            // else we have a buffer to insert into, so try insert.
            // if targetBufferId buffer is now full, or marked as full, will get a false.  the next thread to call this function will swap the buffer out.
            return std::pair<bool, BufferIdType>(this->pool[targetBufferId].append(data, count), fullBufferId);
          }
        }

      protected:
        /**
         * Swap in an empty Buffer from BufferPool at the dest location in MessageBuffers.  The old buffer is returned.
         * The new buffer may be BufferPoolType::INVALID (invalid buffer, when there is no available Buffer from pool)
         *
         * effect:  bufferIds[dest] gets a new Buffer Id, full Buffer is set to "blocked" to prevent further append.
         *
         * THREAD SAFE.   Uses std::compare_exchange_strong to ensure that another thread hasn't swapped in a new Empty one already.
         *
         * NOTE: we make the caller get the new BufferIds[dest] directly instead of returning by reference in the method parameter variable
         *  because this way there is less of a chance of race condition if a shared "old" variable was accidentally used.
         *    and the caller will likely prefer the most up to date bufferIds[dest] anyway.
         *
         * @param dest    position in BufferIds array to swap out
         * @return        the BufferId that was swapped out.
         */
        template<bliss::concurrent::ThreadSafety TS = ThreadSafety>
        typename std::enable_if<TS, BufferIdType>::type swapInEmptyBuffer(const int& dest, const BufferIdType& old) {

          /// get a new buffer and try to replace the existing.
          auto newBuffer = this->pool.tryAcquireBuffer();

          bool hasNewBuffer = newBuffer.first;
          BufferIdType newBufferId = newBuffer.second;               // if acquire fails, we have BufferPoolType::INVALID here.
          BufferIdType targetBufferId = old;

          /// now try to set the bufferIds[dest] to the new buffer id (valid, or BufferPoolType::INVALID if can't get a new one)
          // again, targetBufferId is BufferPoolType::INVALID or pointing to a full buffer.
          // use compare_exchange to ensure we only replace if bufferIds[dest] is still targetBufferId (no other thread has changed it)
          if (bufferIds[dest].compare_exchange_strong(targetBufferId, newBufferId, std::memory_order_acq_rel)) {
            // successful exchange.  bufferIds[dest] now has newBufferId value (will be accessed by caller). targetBufferId still has old value (full buffer, to be returned)

            /// block the old buffer.
            if (old != BufferPoolType::INVALID) this->pool[old].block();

            return old;   // value could be full buffer, or could be BufferPoolType::INVALID.
          } else {
            // failed exchange, another thread had already updated bufferIds[dest].
            // bufferIds[dest] will be accessed by caller later, should return BufferPoolType::INVALID as the replaced buffer id.

            // return the unused buffer id to pool.
            if (hasNewBuffer) {
              this->pool.releaseBuffer(newBufferId);
            }

            return BufferPoolType::INVALID;
          }

        }


        /**
         * Swap in an empty Buffer from BufferPool at the dest location in MessageBuffers.  The old buffer is returned.
         * The new buffer may be BufferPoolType::INVALID (invalid buffer, when there is no available Buffer from pool)
         *
         * effect:  bufferIds[dest] gets a new Buffer Id, full Buffer is set to "blocked" to prevent further append.
         *
         * THREAD UNSAFE.

         * NOTE: we make the caller get the new BufferIds[dest] directly instead of returning by reference in the method parameter variable
         *  because this way there is less of a chance of race condition if a shared "old" variable was accidentally used.
         *    and the caller will likely prefer the most up to date bufferIds[dest] anyway.
         *
         * @param dest    position in BufferIds array to swap out
         * @return        the BufferId that was swapped out.
         */
        template<bliss::concurrent::ThreadSafety TS = ThreadSafety>
        typename std::enable_if<!TS, BufferIdType>::type swapInEmptyBuffer(const int& dest, const BufferIdType& old) {


          /// get a new buffer to replace the existing.  set it whether tryAcquire succeeds or fails (id = BufferPoolType::INVALID)
          bufferIds[dest] = this->pool.tryAcquireBuffer().second;

          /// blocks the old buffer
          if (old != BufferPoolType::INVALID) this->pool[old].block();

          return old;
        }



    };

  } /* namespace io */
} /* namespace bliss */

#endif /* MESSAGEBUFFERS_HPP_ */
