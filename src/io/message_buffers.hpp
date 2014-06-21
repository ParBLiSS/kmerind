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
        const int bufferCapacity;

        /**
         * superclass MessageBuffers is not to be instantiated directly.
         * default, internal pool has unlimited capacity.
         *
         * @param numDests          number of destinations.  e.g. mpi comm size.
         * @param buffer_capacity   individual buffer's size in bytes
         */
        MessageBuffers(const int & numDests, const int & _buffer_capacity) : pool(_buffer_capacity), bufferCapacity(_buffer_capacity) {};


        const int getBufferCapacity() {
          return bufferCapacity;
        }

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

        size_t getSize() {
          return bufferIds.size();
        }




        /**
         * function appends data to the target message buffer.   internally will try to swap in a new buffer when full, and notify the caller of the full buffer's id.
         * need to return 3 things:
         *  status of current insert (in case the bufferpool is fixed size, and there is no available place to buffer)
         *  indicator that there is a full buffer, and the full buffer's id.
         *
         *  TODO:  show a table for compare exchange values and behavior
         *
         * @param[in] data
         * @param[in] count
         * @param[in] dest
         * @param[out] fullBufferId   id of the full buffer, if full.  if not full, -1.
         * @return
         */
        bool append(const void* data, const size_t &count, const int &dest, BufferIdType& oldBufferId) throw (bliss::io::IOException) {

          /// if there is not enough room for the new data in even a new buffer, LOGIC ERROR IN CODE: throw exception
          if (count > bufferCapacity) {
            std::stringstream ss;
            ss << "ERROR: MessageBuffer append with count " << count << " larger than buffer capacity " << bufferCapacity;
            throw (bliss::io::IOException(ss.str()));
          }

          /// default fullBufferId return value.
          BufferIdType fullBufferId = -1;

          /// get the current Buffer's Id
          BufferIdType targetBufferId = bufferIds[dest].load(std::memory_order_consume);

          /// need to replace bufferIds[dest], if it's -1, or if new data can't fit
          if ((targetBufferId == -1) ||
              ((this->bufferCapacity - pool[targetBufferId].getSize()) < count)) {
            // at this point, targetBufferId may be -1 or some value, and the local variables may be out of date already.

            /// get a new buffer and try to replace the existing.
            BufferIdType newBufferId = -1;
            bool hasNewBuffer = this->pool.tryAcquireBuffer(newBufferId);
            // if acquire fails, we have -1 for newBufferId.

            /// now try to set the bufferIds[dest] to the new buffer id (valid, or -1 if can't get a new one)
            // again, targetBufferId is -1 or pointing to a full buffer.

            // use compare_exchange to ensure we only replace if targetBufferId hasn't changed.
            if (bufferIds[dest].compare_exchange_strong(targetBufferId, newBufferId, std::memory_order_acq_rel)) {
              // successful exchange.  bufferIds[dest] now has newBufferId value. targetBufferId has old value.

              // only a success exchange should return targetBufferId as fullBufferId.
              fullBufferId = targetBufferId;  // value could be full buffer, or could be -1.

              targetBufferId = newBufferId;   // set targetBufferId to the new buffer, prep for append.

            } else {
              // failed exchange, then another thread had already updated bufferIds[dest].

              // leave fullBufferId as -1.

              // targetBufferId on failed exchange will contain the up-to-date bufferIds[dest], which could be -1.

              // return the unused buffer id to pool.
              if (hasNewBuffer)
                pool.releaseBuffer(newBufferId);

            }
          }

          oldBufferId = fullBufferId;
          // now try insert using the new targetBufferId.
          if (targetBufferId == -1)   // some other thread may have set this to -1 because it could not get a buffer from pool.  try again.
            return false;
          else
            return this->pool[targetBufferId].append(data, count);

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

        size_t getSize() {
          return bufferIds.size();
        }

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
        bool append(const void* data, const size_t &count, const int &dest, BufferIdType& oldBufferId) throw (bliss::io::IOException) {

          /// if there is not enough room for the new data in even a new buffer, throw exception
          if (count > bufferCapacity) {
            std::stringstream ss;
            ss << "ERROR: MessageBuffer append with count " << count << " larger than buffer capacity " << this->pool[targetBufferId].getCapacity();
            throw (bliss::io::IOException(ss.str()));
          }


          /// initialize fullBufferId to no buffer returned.
          BufferIdType fullBufferId = -1;

          /// save the current Buffer's Id
          BufferIdType targetBufferId = bufferIds[dest];

          /// need to replace bufferIds[dest], if it's -1, or if new data can't fit
          if ((targetBufferId == -1) ||
              ((this->bufferCapacity - pool[targetBufferId].getSize()) < count)) {

            // save the old buffer/full buffer to return
            fullBufferId = targetBufferId;

            /// get a new buffer and try to replace the existing.
            targetBufferId = -1;
            this->pool.tryAcquireBuffer(targetBufferId);

            /// now try to set the bufferIds[dest] to the new buffer id (valid, or -1 if can't get a new one)
            bufferIds[dest] = targetBufferId;

          }

          // now try insert using the new targetBufferId.
          oldBufferId = fullBufferId;
          if (targetBufferId == -1)
            return false;
          else
            return this->pool[targetBufferId].append(data, count);

        }

    };



    /**
     * each RecvMessageBuffers contains a collection of Buffers, one for each received buffer to be processed.
     *
     * not implemented.
     */
//    template<int TAG, bliss::concurrent::ThreadSafety ThreadSafety>
//    class RecvMessageBuffers;


  } /* namespace io */
} /* namespace bliss */

#endif /* MESSAGEBUFFERS_HPP_ */
