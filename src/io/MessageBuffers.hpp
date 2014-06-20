/**
 * @file		MessageBuffers.hpp
 * @ingroup
 * @author	tpan
 * @brief   buffering data for MPI send/receive
 * @details MessageBuffers provides buffering for MPI Send/Receive messages.
 *          Terminology:  Message is an MPI compatible message, serialized from buffer, including the appropriate metadata.
 *                        Element is a single piece of data from the application code.
 *                        Buffer is a collection of elements.
 *
 *
 *          provides functions for storing single and collection of elements into the buffer, and function to serialize
 *          the buffer into a message.
 *
 *          2 modes:  send:  application code populates the buffer either 1 element at a time or multiple elements at a time
 *                            then application code will request the message
 *                    recv:  application code gets a buffer region and receive from comm channel
 *                            application code will request the entire buffer.
 *
 *          each element is destined for a particular mpi target rank, along with a tag that denotes its type.  The add method uses a byte array along with size.
 *          buffer has a limited size.
 *
 *          This class contains a vector of buffers, each with a different id.  The tag is set as template parameter, since we need to know type corresponding to the
 *          tag at compile time anyways.
 *
 *          The class allows sharing between thread but it can certainly be used within a single thread.
 *
 *          send/recv have different semantics?  single class would be nice - easier to support different patterns
 *            buffer -> send = send
 *            recv -> buffer = recv
 *            recv -> send   = forward
 *            buffer -> buffer = local, same rank send.
 *
 *          event driven with observer pattern on full/receive with CommLayer as handler? - avoids looping to check if full, before enqueue.
 *            comm layer can do this work without being an observer since it has to have a send method - can check for full and send.
 *
 *
 *
 *
 *          buffer's lifetime is from the first call to the buffer function to the completion of the MPI send, or
 *              from the start of the MPI recv to the consumption of the buffer.
 *
 *          QUESTIONS: reuse buffer?  std::move or copy content?
 *          MAY NEED A BUFFER OF THESE.  Do that later.
 *
 *
 *
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef MESSAGEBUFFERS_HPP_
#define MESSAGEBUFFERS_HPP_

#include <cassert>

#include <mutex>
#include <vector>

#include "config.hpp"
#include "io/buffer.hpp"
#include "io/buffer_pool.hpp"


namespace bliss
{
  namespace io
  {

    /**
     * @class			bliss::io::MessageBuffers
     * @brief     a data structure to buffering/batching messages for communication
     * @details   THREAD_SAFE enables std::mutex for thread safety.  this may introduce contention
     *            this class uses Buffer class for the actual in memory storage.
     *            this class also uses BufferPool to reuse memory.
     *
     *            this class's purpose is to manage the actual buffering of the data for a set of message targets.
     *              there is no need to have the MessageBuffers be aware of the tags - user can manage that.
     *            when a Buffer is full, we swap in an empty Buffer from BufferPool.  The  full Buffer
     *              is added to a processing queue by Buffer Id (used to get access the buffer from BufferPool).
     *
     *            CommLayer manages the "inProgress" queues, although the actual memory used in the inProgress queues are
     *              retrieved from the Buffers used here.
     *            Inter-thread communication is handled by the ThreadSafeQueues
     *
     *
     *            2 potential subclasses: send and receive, because they have different flow of control.
     *              RecvMessageBuffers is not needed at the moment so not implemented.
     *
     */
    template<int TAG, bliss::concurrent::ThreadSafety ThreadSafety>
    class MessageBuffers
    {
      protected:
        bliss::io::BufferPool<ThreadSafety> pool;

        /**
         * superclass MessageBuffers is not to be instantiated directly.
         * default, internal pool has unlimited capacity.
         *
         * @param numDests          number of destinations.  e.g. mpi comm size.
         * @param buffer_capacity   individual buffer's size in bytes
         */
        MessageBuffers(const int & numDests, const int & buffer_capacity) : pool(buffer_capacity) {};

      public:
        virtual ~MessageBuffers() {};

//        virtual bool acquireBuffer(size_t &id);
        /**
         * after the buffer has been consumed, release it back to the pool
         * @param id
         */
        void releaseBuffer(const size_t &id) throw (bliss::io::IOException) {
          pool.releaseBuffer(id);
        }

      private:

    };

    /**
     * each SendMessageBuffers class contains an vector of Buffer Ids (reference buffers in BufferPool), one for each target id
     *
     * for ones that are full and are not in the vector, the user of this class will need to track the Id and use this class to
     * access the internal BufferPool by id.
     *
     */
    template<int TAG, bliss::concurrent::ThreadSafety ThreadSafety>
    class SendMessageBuffers;

    /// specializations
    template<int TAG>
    class SendMessageBuffers<TAG, bliss::concurrent::THREAD_SAFE> : public MessageBuffers<TAG, bliss::concurrent::THREAD_SAFE>
    {
      protected:
        typedef typename bliss::io::BufferPool<bliss::concurrent::THREAD_SAFE>::IdType    BufferIdType;

        std::vector<std::atomic<BufferIdType> > bufferIds;  // mapping from process id (0 to vector size), to buffer Ids (from BufferPool)

      public:
        SendMessageBuffers(const int & numDests, const int & buffer_capacity) :
          MessageBuffers<TAG, bliss::concurrent::THREAD_SAFE>(numDests, buffer_capacity),
          bufferIds(numDests) {
          /// initialize the bufferIds by acquiring them from the pool.

          BufferIdType t;
          for (int i = 0; i < numDests; ++i) {
            pool.tryAcquireBuffer(t);
            bufferIds[i] = t;
          }
        };

        virtual ~SendMessageBuffers() {};

        /**
         * function appends data to the target message buffer.   internally will try to swap in a new buffer when full, and notify the caller of the full buffer's id.
         * need to return 3 things:
         *  status of current insert (in case the bufferpool is fixed size, and there is no available place to buffer)
         *  indicator that there is a full buffer, and the full buffer's id.
         *
         * @param[in] data
         * @param[in] count
         * @param[in] dest
         * @param[out] fullBufferId   id of the full buffer, if full.  if not full, -1.
         * @return
         */
        bool append(const void* data, const size_t &count, const int &dest, BufferIdType& fullBufferId) throw (bliss::io::IOException) {

          /// initialize fullBufferId
          fullBufferId = -1;

          size_t cap = 0;

          /// get the current Buffer's Id
          BufferIdType targetBufferId = bufferIds[dest].load(std::memory_order_consume);
          bool canFit = (targetBufferId != -1);   // if targetBufferId == -1, definitely can't fit.

          if (canFit) {
            cap = pool[targetBufferId].getCapacity();
            canFit = ((cap - pool[targetBufferId].getSize()) > count);
          }

          /// if there is not enough room for the new data in even a new buffer, LOGIC ERROR IN CODE: throw exception
          if (count > cap) {
            std::stringstream ss;
            ss << "ERROR: MessageBuffer append with count " << count << " larger than buffer capacity " << this->pool[targetBufferId].getCapacity();
            throw (bliss::io::IOException(ss.str()));
          }


          /// need to replace bufferIds[dest] if it's -1, or if new data can't fit
          if (!canFit) {
            BufferIdType newBufferId = -1;
            // at this point, targetBufferId may be -1 or some value, and the local variables may be out of date already.

            /// get a new buffer and try to replace the existing.
            if (this->pool.tryAcquireBuffer(newBufferId)) {
              // successfully gotten a new buffer. so now set it. again, targetBufferId is -1 or pointing to a full buffer.
              if (bufferIds[dest].compare_exchange_strong(targetBufferId, newBufferId, std::memory_order_acq_rel)) {
                // successful exchange.  bufferIds[dest] now has newBufferId value. targetBufferId has old value.
                // only a success exchange should return targetBufferId as fullBufferId.
                fullBufferId = targetBufferId;  // could be full buffer, or could be -1.

                targetBufferId = newBufferId;

              } else {
                // failed exchange, then another thread had already exchanged it.  so return the unused buffer to pool and return false.
                pool.releaseBuffer(newBufferId);
                // targetBufferId on failed exchange will contain the up-to-date bufferIds[dest], which could be -1.

                // leave fullBufferId as -1.
              }

              // now try insert using teh new targetBufferId.
              if (targetBufferId == -1) {   // some thread may have failed to get a newBufferId from pool, and exchanged a full buffer with -1.
                return false;
              } else {
                return this->pool[targetBufferId].append(data, count);
              }

            } else {
              // no buffers are available. append fails with targetBufferId set to -1.
              // if no thread has changed bufferIds[dest] yet then set it to -1.
              if (targetBufferId != -1) {
                if (bufferIds[dest].compare_exchange_strong(targetBufferId, -1, std::memory_order_acq_rel)) {
                  // targetBufferId is the old full buffer
                  fullBufferId = targetBufferId;
                } // else some thread changed targetBufferId, so don't want to mark the buffer as full and return.
              }

              return false;
            }


          }

          /// if no target buffer currently, then try to set one.  multiple threads can be doing this, only 1 should succeed.
//          if (targetBufferId == -1) {  // no target buffer (previous full and being sent, and pool has no available ones.
//            // if no target buffer, try get a new one.
//            if (this->pool.tryAcquireBuffer(newBufferId)) {
//              // successfully gotten a new buffer. so now set it
//              if (bufferIds[dest].compare_exchange_strong(targetBufferId, newBufferId, std::memory_order_acq_rel)) {
//                // successful exchange.  bufferIds[dest] now has newBufferId value. targetBufferId has old value (-1).
//                targetBufferId = newBufferId;
//              } else {
//                // if failed, then another thread had already exchanged it.  so return the unused buffer.
//                pool.releaseBuffer(newBufferId);
//                // targetBufferId on failed exchange will contain the up-to-date bufferIds[dest].
//              }
//            } else {
//              // no buffers are available. append failed with targetBufferId left as is. (either -1, or another buffer had updated it.  caller should retry.
//              return false;
//            }
//          }
//
//          /// try inserting into the target buffer
//          if (this->pool[targetBufferId].append(data, count)) {
//            // success.  done, return true;
//            return true;
//          }



        }

    };

    template<int TAG>
    class SendMessageBuffers<TAG, bliss::concurrent::THREAD_UNSAFE> : public MessageBuffers<TAG, bliss::concurrent::THREAD_UNSAFE>
    {
      protected:
        typedef typename bliss::io::BufferPool<bliss::concurrent::THREAD_UNSAFE>::IdType    BufferIdType;

        std::vector<BufferIdType> bufferIds;  // mapping from process id (0 to vector size), to buffer Ids (from BufferPool)

      public:
        SendMessageBuffers(const int & numDests, const int & buffer_capacity) :
          MessageBuffers<TAG, bliss::concurrent::THREAD_UNSAFE>(numDests, buffer_capacity),
          bufferIds(numDests) {
          /// initialize the bufferIds by acquiring them from the pool.

          BufferIdType t;
          for (int i = 0; i < numDests; ++i) {
            pool.tryAcquireBuffer(t);
            bufferIds[i] = t;
          }
        };
        virtual ~SendMessageBuffers() {};

        /**
         * function appends data to the target message buffer.   internally will try to swap in a new buffer when full, and notify the caller of the full buffer's id.
         * need to return 3 things:
         *  status of current insert (in case the bufferpool is fixed size, and there is no available place to buffer)
         *  indicator that there is a full buffer, and the full buffer's id.
         *
         * @param[in] data
         * @param[in] count
         * @param[in] dest
         * @param[out] fullBufferId   id of the full buffer, if full.  if not full, -1.
         * @return
         */
        bool append(const void* data, const size_t &count, const int &dest, BufferIdType& fullBufferId) throw (bliss::io::IOException) {

          /// initialize fullBufferId to no buffer returned.
          fullBufferId = -1;

          /// save the current Buffer's Id
          BufferIdType targetBufferId = bufferIds[dest];

          /// if the the target Id is -1, then no active buffer is present, so get one.
          if (targetBufferId == -1) {   // no target buffer (previous full and being sent, and pool has no available ones.
            // if no target buffer, try get a new one.
            if (this->pool.tryAcquireBuffer(targetBufferId)) {
              bufferIds[dest] = targetBufferId;
            } else {
              // can't get a new buffer, so return append failed.
              return false;
            }
          }

          /// try inserting into the target buffer
          if (this->pool[targetBufferId].append(data, count)) {
            // success, so return true.

            return true;
          }


          /// if can't insert (buffer full)

          /// if there is not enough room for the new data in even a new buffer, throw exception
          if (count > this->pool[targetBufferId].getCapacity()) {
            std::stringstream ss;
            ss << "ERROR: MessageBuffer append with count " << count << " larger than buffer capacity " << this->pool[targetBufferId].getCapacity();
            throw (bliss::io::IOException(ss.str()));
          }

          /// save the targetId for caller processing
          fullBufferId = targetBufferId;

          /// try to get a new bufferId from the pool
          if (this->pool.tryAcquireBuffer(targetBufferId)) {
            /// if can get a new bufferId, set the new bufferId into the vector, and retry insert, and return result of insertion
            bufferIds[dest] = targetBufferId;

            return this->pool[targetBufferId].append(data, count);
          } else {
            /// if can't get a new bufferId, pool is empty, so set the bufferIds to -1, and return the full buffer anyways so can be cleared.
            bufferIds[dest] = -1;

            /// could not get a new buffer right now, so fail
            return false;
          }

        }

    };



    /**
     * each RecvMessageBuffers contains a collection of Buffers, one for each received buffer to be processed.
     *
     * not implemented.
     */
    template<int TAG, bliss::concurrent::ThreadSafety ThreadSafety>
    class RecvMessageBuffers;


  } /* namespace io */
} /* namespace bliss */

#endif /* MESSAGEBUFFERS_HPP_ */
