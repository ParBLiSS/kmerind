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
     * @class			MessageBuffers
     * @brief     a data structure to buffering/batching messages for communication
     * @details   Templated base class for a collection of buffers that are used for communication.
     *
     *
     *            The class contains a BufferPool, which contains an array of reusable Buffer objects for in-memory storage.
     *
     * @tparam ThreadSafety   Indicates whether the class should be thread safe or not.
     */
    template<bliss::concurrent::ThreadSafety ThreadSafety>
    class MessageBuffers
    {
        // TODO: move consturctor and assignment operator to copy between thread safety levels.
      protected:

        /// Internal BufferPool Type.  typedefed only to shorten the usage.
        typedef typename bliss::io::BufferPool<ThreadSafety>  BufferPoolType;

        /// IdType of a Buffer, aliased from BufferPool
        typedef typename BufferPoolType::IdType    BufferIdType;

        /// a pool of in-memory Buffers for storage.
        BufferPoolType pool;

        /**
         * Gets the capacity of the Buffer instances.
         * @return    capacity of each buffer.
         */
        const size_t getBufferCapacity() const {
          return pool.getBufferCapacity();
        }

        /**
         * @brief Constructor.  creates a buffer given a specified per-buffer capacity, and the BufferPool capacity.
         * @note  Protected so base class function is not called directly.
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
         * @brief default copy constructor.  deleted.  since internal BufferPool does not allow copy construction/assignment.
         * @note Protected so base class function is not called directly.
         *
         * @param other     source MessageBuffers to copy from
         */
        explicit MessageBuffers(const MessageBuffers<ThreadSafety>& other) = delete;

        /**
         * @brief default copy assignment operator.  deleted.   since internal BufferPool does not allow copy construction/assignment.
         * @note Protected so base class function is not called directly.
         *
         * @param other     source MessageBuffers to copy from
         * @return          self.
         */
        MessageBuffers<ThreadSafety>& operator=(const MessageBuffers<ThreadSafety>& other) = delete;


        /**
         * @brief default move constructor.
         * @note Protected so base class function is not called directly.
         *
         * @param other     source MessageBuffers to move from
         */
        explicit MessageBuffers(MessageBuffers<ThreadSafety>&& other) : pool(std::move(other.pool)) {};


        /**
         * @brief default move assignment operator.
         * @note Protected so base class function is not called directly.
         *
         * @param other     source MessageBuffers to move from
         * @return          self.
         */
        MessageBuffers<ThreadSafety>& operator=(MessageBuffers<ThreadSafety>&& other) {
          pool = std::move(other.pool);
          return *this;
        }


      public:
        /**
         * @brief default destructor
         */
        virtual ~MessageBuffers() {};

        /**
         * @brief Accesses the Buffer object inside the BufferPool by BufferPool assigned Id
         * @param id      Id of the Buffer to be accessed.
         * @return        const reference to the Buffer object.
         */
        const Buffer<ThreadSafety> & getBackBuffer(const BufferIdType& id) const {
          return this->pool.at(id);
        }

        /**
         * @brief Releases a Buffer by it's BufferPool id, after the buffer is no longer needed.
         * @note  this should be called on a buffer that is not being used by MessageBuffers, i.e. flush does not satisfy this.
         *
         * @note  This to be called after the communication logic is done with the Send or Recv buffer.
         *
         * @param id    Buffer's id (assigned from BufferPool)
         */
        void releaseBuffer(const BufferIdType &id) {
          pool.releaseBuffer(id);
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
     * @note
     *  For Buffers that are full, the calling thread will need to track the full buffers' Ids as this class evicts a full Buffer's id from it's
     *    vector of Buffer Ids.  The Calling Thread can do t his via a queue, as in the case of CommunicationLayer.
     *  Note that a Buffer is "in use" from the first call to the Buffer's append function to the completion of the send (e.g. via MPI Isend)
     *
     *  @tparam ThreadSafety    determines if this class should be thread safe
     */
    template<bliss::concurrent::ThreadSafety ThreadSafety>
    class SendMessageBuffers : public MessageBuffers<ThreadSafety>
    {
      public:
        /// Id type of the Buffers
        using BufferIdType = typename MessageBuffers<ThreadSafety>::BufferIdType;

        /// the BufferPoolType from parent class
        using BufferPoolType = typename MessageBuffers<ThreadSafety>::BufferPoolType;

      protected:
        /// Id type of the Buffers, but here depending on the specified ThreadSafety template parameter, either BufferIdType or atomic version of BufferIdType
        using IdType = typename std::conditional<ThreadSafety, std::atomic<BufferIdType>, BufferIdType>::type;

        /// Vector of IdType (atomic or not).  Provides a mapping from process id (0 to vector size), to buffer Ids (from BufferPool)
        std::vector< IdType > bufferIdForProcId;

//        /// map from IdType to process id.  for quickly checking if a bufferId is in use.
//        std::unordered_map<IdType, int> procIdForBufferId;


        std::mutex mutex;


      public:
        /**
         * @brief Constructor.
         * @param numDests         The number of messaging targets/destinations
         * @param bufferCapacity   The capacity of the individual buffers.  default 8192.
         * @param poolCapacity     The capacity of the pool.  default unbounded.
         */
        explicit SendMessageBuffers(const int & numDests, const size_t & bufferCapacity = 8192, const BufferIdType & poolCapacity = std::numeric_limits<BufferIdType>::max()) :
          MessageBuffers<ThreadSafety>(bufferCapacity, poolCapacity), bufferIdForProcId(numDests)
        {
          //== initialize the bufferIdForProcId, acquiring them from the pool by calling reset.
          this->reset();
        };

        /**
         * @brief Default constructor, deleted
         */
        SendMessageBuffers() :  MessageBuffers<ThreadSafety>() {};

        /**
         * @brief default copy constructor.  deleted.
         * @param other   source SendMessageBuffers to copy from.
         */
        explicit SendMessageBuffers(const SendMessageBuffers<ThreadSafety> &other) = delete;

        /**
         * default move constructor.  calls superclass move constructor first.
         * @param other   source SendMessageBuffers to move from.
         */
        explicit SendMessageBuffers(SendMessageBuffers<ThreadSafety> && other) :
            MessageBuffers<ThreadSafety>(std::move(other)), bufferIdForProcId(std::move(other.bufferIdForProcId)) {};

        /**
         * @brief default copy assignment operator, deleted.
         * @param other   source SendMessageBuffers to copy from.
         * @return        self
         */
        SendMessageBuffers<ThreadSafety>& operator=(const SendMessageBuffers<ThreadSafety> &other) = delete;


        /**
         * @brief default move assignment operator.
         * @param other   source SendMessageBuffers to move from.
         * @return        self
         */
        SendMessageBuffers<ThreadSafety>& operator=(SendMessageBuffers<ThreadSafety> && other) {
          bufferIdForProcId = std::move(other.bufferIdForProcId);
          this->pool = std::move(other.pool);
          return *this;
        }

        /**
         * @brief default destructor
         */
        virtual ~SendMessageBuffers() {};

        /**
         * @brief get the number of buffers.  should be same as number of targets for messages
         * @return    number of buffers
         */
        const size_t getSize() const {
          return bufferIdForProcId.size();
        }

        /**
         * @brief get the BufferId that a messaging target maps to.
         * @details for simplicity, not distinguishing between thread safe and unsafe versions
         *
         * @param targetRank   the target id for the messages.
         * @return      the buffer id that the target id maps to.
         */
        inline const BufferIdType getBufferIdForRank(const int targetRank) const {
          return BufferIdType(bufferIdForProcId.at(targetRank));
        }

//        /**
//         * @brief get the  messaging target that a bufferId maps to.
//         * @details for simplicity, not distinguishing between thread safe and unsafe versions
//         *
//         * @param targetRank   the target id for the messages.
//         * @return      the buffer id that the target id maps to.
//         */
//        inline const int getRankForBufferId(const BufferIdType targetBufferId) const {
//          if (procIdForBufferId.find(targetBufferId) == procIdForBufferId.end()) return -1;
//          else return procIdForBufferId.at(targetBufferId);
//        }


        /**
         * @brief get the list of active Buffer ids.  (for debugging)
         * @return  vector containing the active Buffer ids.
         */
        const std::vector<IdType>& getBufferIdsForAllRanks() const {
          return bufferIdForProcId;
        }

        /**
         * @brief convenience method to convert the list of Active Buffer Ids to string.  (for debugging)
         * @return   string with list of active Buffer Ids
         */
        const std::string bufferIdsToString() const {
          std::stringstream ss;
          std::ostream_iterator<IdType> ost(ss, ",");
          std::copy(bufferIdForProcId.begin(), bufferIdForProcId.end(), ost);
          ss << std::endl;
          return ss.str();
        }

        /**
         * @brief Reset the current MessageBuffers instance by first clearing its list of Buffer Ids, then repopulate it from the pool.
         */
        virtual void reset() {
          MessageBuffers<ThreadSafety>::reset();
          for (int i = 0; i < this->bufferIdForProcId.size(); ++i) {
            auto newBuffer = std::move(this->pool.tryAcquireBuffer());
            if (!newBuffer.first)
              throw std::logic_error("not enough buffers in BufferPool to support this messagebuffers object.");

            bufferIdForProcId.at(i) = newBuffer.second;
//            procIdForBufferId[newBuffer.second] = i;
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
        std::pair<bool, BufferIdType> append(const void* data, const size_t count, const int targetProc) {

          std::unique_lock<std::mutex> lock(mutex);

          //== if count is 0, no write and this succeeds right away.
          if (count == 0)
            return std::move(std::pair<bool, BufferIdType>(true, BufferPoolType::ABSENT));

          //== if there is not enough room for the new data in even a new buffer, LOGIC ERROR IN CODE: throw exception
          if (count > this->getBufferCapacity()) {
            throw (std::invalid_argument("ERROR: messageBuffer append with count exceeding Buffer capacity"));
          }

          //== if targetProc is outside the valid range, throw an error
          if (targetProc < 0 || targetProc > getSize()) {
            throw (std::invalid_argument("ERROR: messageBuffer append with invalid targetProc"));
          }

          //== if data being set is null, throw error
          if (data == nullptr) {
            throw (std::invalid_argument("ERROR: calling MessageBuffer append with nullptr"));
          }

          bool appendResult = false;
          //== get the current Buffer's Id. local var because we need to check value AND get back buffer.
          // this means that the targetBufferId may not be associated with targetProc by the time append is called due to threading.
          // TODO:  fix:  use a dummy buffer that does not allow append.  getBackBuffer with ABSENT buffer id returns the dummy buffer.
          // this way, the actual buffer to insert is resolved in 1 single chained call with minimal branching.
          BufferIdType targetBufferId = getBufferIdForRank(targetProc);
          // now try to insert into the current bufferId.
          if (targetBufferId != BufferPoolType::ABSENT) {

            if (MessageBuffers<ThreadSafety>::getBackBuffer(targetBufferId).getCapacity() == 1) {
            	printf("ERROR 2!!!! targetProc %d maps to buffer id saved %d, current %d, pool size %d\n", targetProc, targetBufferId, getBufferIdForRank(targetProc), this->pool.getSize());
            	fflush(stdout);
              assert(MessageBuffers<ThreadSafety>::getBackBuffer(targetBufferId).getCapacity() != 1);
            }
            // else we have a buffer to insert into, so try insert.
            // if fails, then the buffer is full or blocked (sending or released)
            appendResult = MessageBuffers<ThreadSafety>::getBackBuffer(targetBufferId).append(data, count);

          } // else we don't have a buffer to insert into, so result is false.

          //== default fullBufferId return value.
          BufferIdType fullBufferId = BufferPoolType::ABSENT;
          // if insert result is false, then get a new buffer for the next time to insert
          if (!appendResult) {
              //== if targetBufferId is an invalid buffer, or if new data can't fit, need to replace bufferIdForProcId[dest]
              if (targetBufferId != BufferPoolType::ABSENT && MessageBuffers<ThreadSafety>::getBackBuffer(targetBufferId).getCapacity() == 1) {
                printf("ERROR 1!!!! targetProc %d maps to buffer id saved %d, current %d, pool size %d\n", targetProc, targetBufferId, getBufferIdForRank(targetProc), this->pool.getSize());
                fflush(stdout);
                assert(MessageBuffers<ThreadSafety>::getBackBuffer(targetBufferId).getCapacity() != 1);
              }

              // conditions are either no buffer to insert into, or full buffer, or blocked buffer.
    			// this call will mark the fullBuffer as blocked to prevent multiple write attempts while flushing the buffer
                fullBufferId = swapInEmptyBuffer<ThreadSafety>(targetProc, targetBufferId);    // swap in an empty buffer
                //targetBufferId = getBufferIdForRank(targetProc);               // and get the updated targetBuffer id


          }
          lock.unlock();
          // don't try to reinsert right here.




//          //== get the current Buffer's Id. local var to ensure all three lines use the same value.
//          BufferIdType targetBufferId = getBufferIdForRank(targetProc);
//          //== if targetBufferId is an invalid buffer, or if new data can't fit, need to replace bufferIdForProcId[dest]
//          if (targetBufferId != BufferPoolType::ABSENT && MessageBuffers<ThreadSafety>::getBackBuffer(targetBufferId).getCapacity() == 1) {
//            printf("ERROR 1!!!! targetProc %d maps to buffer id saved %d, current %d, pool size %d\n", targetProc, targetBufferId, getBufferIdForRank(targetProc), this->pool.getSize());
//            fflush(stdout);
//            assert(MessageBuffers<ThreadSafety>::getBackBuffer(targetBufferId).getCapacity() != 1);
//          }
//
//          if ((targetBufferId == BufferPoolType::ABSENT) ||
//        		  MessageBuffers<ThreadSafety>::getBackBuffer(targetBufferId).getSize() > (this->getBufferCapacity() - count)) {
//            // at this point, the local variables may be out of date already, and targetBufferId may already marked as full.
//
//            // this call will mark the fullBuffer as blocked to prevent multiple write attempts while flushing the buffer
//            fullBufferId = swapInEmptyBuffer<ThreadSafety>(targetProc, targetBufferId);    // swap in an empty buffer
//            //targetBufferId = getBufferIdForRank(targetProc);               // and get the updated targetBuffer id
//          }
//
//          bool appendResult = false;
//          targetBufferId = getBufferIdForRank(targetProc);
//          // now try to insert into the current bufferId.
//          if (targetBufferId != BufferPoolType::ABSENT) {
//
//            if (MessageBuffers<ThreadSafety>::getBackBuffer(targetBufferId).getCapacity() == 1) {
//            	printf("ERROR 2!!!! targetProc %d maps to buffer id saved %d, current %d, pool size %d\n", targetProc, targetBufferId, getBufferIdForRank(targetProc), this->pool.getSize());
//            	fflush(stdout);
//              assert(MessageBuffers<ThreadSafety>::getBackBuffer(targetBufferId).getCapacity() != 1);
//            }
//            // else we have a buffer to insert into, so try insert.
//            // if targetBufferId buffer is now full, or marked as full, will get a false.  the next thread to call this function will swap the buffer out.
//            appendResult = MessageBuffers<ThreadSafety>::getBackBuffer(targetBufferId).append(data, count);
//          } // else we don't have a buffer to insert into, so result is false.
//

          return std::move(std::make_pair(appendResult, fullBufferId));

        }


        BufferIdType flushBufferForRank(const int targetProc) {

          DEBUG("flush buffer for rank!");

          //== if targetProc is outside the valid range, throw an error
          if (targetProc < 0 || targetProc > getSize()) {
            throw (std::invalid_argument("ERROR: messageBuffer append with invalid targetProc"));
          }

          // passing in getBufferIdForRank result to ensure atomicity.
          // may return ABSENT if there is no available buffer to use.
          return swapInEmptyBuffer<ThreadSafety>(targetProc, getBufferIdForRank(targetProc));

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
         * @note
         * 	old is 		old assignment	action
         *    ABSENT      available     not a valid case
         *    ABSENT      at dest       swap, return ABS
         *    ABSENT      not at dest   dest is already swapped.  no op, return ABS
         *    not ABS     available     not used.  block old, no swap. return ABS
         *    not ABS     at dest       swap, block old.  return old
         *    not ABS     not at dest   being used.  no op. return ABS
         *
         * @tparam        used to choose thread safe vs not verfsion of the method.
         * @param dest    position in bufferIdForProcId array to swap out
         * @return        the BufferId that was swapped out.
         */
        template<bliss::concurrent::ThreadSafety TS = ThreadSafety>
        typename std::enable_if<TS, BufferIdType>::type swapInEmptyBuffer(const int dest, const BufferIdType old) {

          //== block the old buffer if it's not ABSENT, is available.  this means that it's not being used by message buffer, so swap is not going to happen
          if ((old != BufferPoolType::ABSENT) && this->pool.isAvailable(old)) {
        	  MessageBuffers<ThreadSafety>::getBackBuffer(old).block();
            return BufferPoolType::ABSENT;
          }

          //== get a new buffer and try to replace the existing.
          auto newBufferInfo = this->pool.tryAcquireBuffer();

          BufferIdType newBufferId = newBufferInfo.second;               // if acquire fails, we have BufferPoolType::ABSENT here.
          BufferIdType targetBufferId = old;

          //== now try to set the bufferIdForProcId[dest] to the new buffer id (valid, or BufferPoolType::ABSENT if can't get a new one)
          // again, targetBufferId is BufferPoolType::ABSENT or pointing to a full buffer.
          // use compare_exchange to ensure we only replace if bufferIdForProcId[dest] is still targetBufferId (no other thread has changed it)
          if (bufferIdForProcId.at(dest).compare_exchange_strong(targetBufferId, newBufferId, std::memory_order_acq_rel)) {
            // successful exchange.  bufferIdForProcId[dest] now has newBufferId value (will be accessed by caller). targetBufferId still has old value (full buffer, to be returned)

            //== block the old buffer if it's not ABSENT, is being used by MessageBuffers, and is associated with "dest".  else another call had already blocked it.
            if (targetBufferId != BufferPoolType::ABSENT) MessageBuffers<ThreadSafety>::getBackBuffer(targetBufferId).block();


            // return old buffer id.  value could be full buffer, or could be BufferPoolType::ABSENT, or a partially full buffer if this is called manually.
            return targetBufferId;
          } else {
            // failed exchange, another thread had already updated bufferIdForProcId[dest], INCLUDING PROCESSING THE OLD FULL BUFFER.
            // bufferIdForProcId[dest] will be accessed by caller later, should return BufferPoolType::ABSENT as the replaced buffer id.

            // return the unused buffer id to pool.
            if (newBufferInfo.first) {
            	MessageBuffers<ThreadSafety>::getBackBuffer(newBufferId).block();
              this->pool.releaseBuffer(newBufferId);
            }

            // if failed exchange, the old buffer still needs to be processed, so need to return old data type.
            //(targetBufferId now holds bufferIdForProcId[dest], an active buffer id.)
            return BufferPoolType::ABSENT;

          }
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
         * @tparam TS     used to choose thread safe vs not verfsion of the method.
         * @param dest    position in bufferIdForProcId array to swap out
         * @return        the BufferId that was swapped out.
         */
        template<bliss::concurrent::ThreadSafety TS = ThreadSafety>
        typename std::enable_if<!TS, BufferIdType>::type swapInEmptyBuffer(const int dest, const BufferIdType old) {

          //== block the old buffer if it's not ABSENT, is available.  this means that it's not being used by message buffer, so swap is not going to happen
          if ((old != BufferPoolType::ABSENT) && (this->pool.isAvailable(old))) {
        	  MessageBuffers<ThreadSafety>::getBackBuffer(old).block();
            return BufferPoolType::ABSENT;
          }


          if (bufferIdForProcId.at(dest) == old) {
            //== get a new buffer to replace the existing.  set it whether tryAcquire succeeds or fails (id = BufferPoolType::ABSENT)

            auto newBufferInfo = this->pool.tryAcquireBuffer();
            bufferIdForProcId.at(dest) = newBufferInfo.second;
            //== blocks the old buffer
            if (old != BufferPoolType::ABSENT) MessageBuffers<ThreadSafety>::getBackBuffer(old).block();
            return old;
          } else {
            // buffer already replaced for dest.  don't need to do anything
            return BufferPoolType::ABSENT;
          }

        }



    };

  } /* namespace io */
} /* namespace bliss */

#endif /* MESSAGEBUFFERS_HPP_ */
