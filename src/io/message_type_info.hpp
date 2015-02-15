/**
 * @file    message_type_info.hpp
 * @ingroup bliss::io
 * @author  Tony Pan <tpan7@gatech.edu>
 * @brief   a data structure to manage data and control message communications associated with a message type.
 * @details defines TaggedEpoch for presenting tag and id of distinct communication episodes
 *          defines MessageTypeInfo to track active epochs for a Tag.
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef MESSAGE_TYPE_INFO_HPP_
#define MESSAGE_TYPE_INFO_HPP_


#include <atomic>
#include <memory>

#include "concurrent/copyable_atomic.hpp"


namespace bliss
{
  namespace io
  {

    /// convenience function to extract tag from TaggedEpoch
    typedef int64_t TaggedEpoch;
    inline int32_t getTagFromTaggedEpoch(const TaggedEpoch & te) {
      return te >> 32;
    }

    /// convenience function to extract epoch from TaggedEpoch
    inline uint32_t getEpochFromTaggedEpoch(const TaggedEpoch & te) {
      return te & 0x00000000FFFFFFFF;
    }

    /// CONTOL TAG's epoch
    constexpr TaggedEpoch CONTROL_TAGGED_EPOCH = 0;

    /// CONTOL TAG's epoch
    constexpr int32_t CONTROL_TAG = 0;

    /**
     * @class MessageTypeInfo
     * @brief Used to manage active communication episodes for a message type.
     * @details
     *          This class handles the generation of epoch number, store the epoch number
     *          count of control messages with matching tag/epoch for each epoch.
     *
     *          It provides thread safe access to the TaggedEpoch to count mapping,
     *          and provides a mechanism for calling thread to wait for the count to
     *          reach 0 for a TaggedEpoch.  This is done via mutex and condition variables.
     *
     *          finally it manages a pointer to the current MessageBuffers for a message type
     *          via a shared pointer, which provides reference counting.
     *            since MessageBuffers is long running, we don't expect it to change often,
     *            and since MessageBuffers is thread safe, we don't lock MessageBuffer's access.
     *
     */
    template<typename MessageBuffersType>
    class MessageTypeInfo {
      public:

      protected:
        /// shared pointer to MessageBuffers for storing pending communication.
        std::shared_ptr<MessageBuffersType> bufferPtr;

        /// indicates that a Tag (message type) is finished and no future messages will be of this type
        std::atomic<bool>         finished;  // do we need this?

        /// id of the next Epoch
        std::atomic<TaggedEpoch>  nextEpoch;

        /// number of processes involved in communication epoch
        mutable int commSize;

        /// map from epoch to count to track active epochs.
        std::unordered_map<TaggedEpoch, int > epoch_capacities;

        /// for locking access to epoch_capacities
        mutable std::mutex mutex;
        /// for use to wait for an epoch to finish
        mutable std::condition_variable condVar;

        /// deleted copy constructor
        explicit MessageTypeInfo(const MessageTypeInfo& other) = delete;

        /// deleted copy assignment operator
        MessageTypeInfo& operator=(const MessageTypeInfo& other) = delete;

        /// mutex locked move constructor.
        MessageTypeInfo(MessageTypeInfo&& other, const std::lock_guard<std::mutex> &) :
          bufferPtr(std::move(other.bufferPtr)),
          commSize(other.commSize),
          epoch_capacities(std::move(other.epoch_capacities)) {

          finished.store(other.finished.exchange(true, std::memory_order_relaxed));
          nextEpoch.store(other.nextEpoch.exchange(-1, std::memory_order_relaxed));

          other.commSize = 0;
        }


      public:
        /// default constructor.  needed by some containers.
        MessageTypeInfo() : finished(false), nextEpoch(-1), commSize(0) {}

        /// constructor
        MessageTypeInfo(const int32_t tag, const std::shared_ptr<MessageBuffersType>& ptr, const int epoch_capacity) :
          bufferPtr(ptr), finished(false), nextEpoch(static_cast<TaggedEpoch>(tag) << 32), commSize(epoch_capacity) {}

        /// move constructor
        explicit MessageTypeInfo(MessageTypeInfo&& other) :
            MessageTypeInfo(std::forward<MessageTypeInfo>(other), std::lock_guard<std::mutex>(other.mutex)) {}

        /// move assignment operator
        MessageTypeInfo& operator=(MessageTypeInfo&& other) {
          std::unique_lock<std::mutex> lock(mutex, std::defer_lock);
          std::unique_lock<std::mutex> otherlock(other.mutex, std::defer_lock);
          std::lock(lock, otherlock);

          bufferPtr = std::move(other.bufferPtr);

          finished.store(other.finished.exchange(true, std::memory_order_relaxed));
          // -1 in64 is 0xFFFFFFFFFFFFFFFF, which is -1, -1 for 2 ints.
          nextEpoch.store(other.nextEpoch.exchange(-1, std::memory_order_relaxed));
          commSize = other.commSize; other.commSize = 0;
          epoch_capacities = std::move(other.epoch_capacities);

          return *this;
        }


        /// default constructor
        virtual ~MessageTypeInfo() {}

        //===== accessors

        /// mark a message type as finished and return previous value
        bool finish() {
          return finished.exchange(true, std::memory_order_acq_rel);
        }

        /// check if a message type is finished
        bool isFinished() {
          return finished.load(std::memory_order_relaxed);
        }

        /**
         * @brief  get the next epoch id, and insert the epoch into the epoch-capacity map.
         * @detail   note that the mapping may be inserted before the epoch was created.
         *        this can happen if a remote process sent out a control message with that epoch before local process acquired.
         */
        TaggedEpoch acquireEpoch(std::atomic<int>& totalCount) {
          TaggedEpoch te = nextEpoch.fetch_add(1, std::memory_order_relaxed);

          std::lock_guard<std::mutex> lock(mutex);
          if (epoch_capacities.count(te) == 0) {
            epoch_capacities[te] = commSize;
            totalCount.fetch_add(1, std::memory_order_release);
          }
          return te;
        }
        /**
         * @brief  reduce the capacity associated with an epoch by 1.
         * @detail   note that the mapping may not exist yet when a control MPI message is received,
         *        as the local process not yet acquired the epoch.
         *        In this case, a new mapping is created then decremented.
         *        this can happen if a remote process sent out a control message with that epoch before local process acquired.
         */
        int countdownEpoch(const TaggedEpoch & te, std::atomic<int>& totalCount) throw (std::out_of_range) {
          std::lock_guard<std::mutex> lock(mutex);
          if (epoch_capacities.count(te) == 0) {
            epoch_capacities[te] = commSize;
            totalCount.fetch_add(1, std::memory_order_release);
            WARNINGF("epoch does not exist to count down. created: %d %d", getTagFromTaggedEpoch(te), getEpochFromTaggedEpoch(te));
          }
          return --epoch_capacities.at(te);
        }

        /// check if an epoch is finished.
        bool isEpochFinished(const TaggedEpoch & te) throw (std::out_of_range) {
          std::lock_guard<std::mutex> lock(mutex);
          if (epoch_capacities.count(te) == 0) {
            ERRORF("epoch does not exist to check for Finished: %d %d", getTagFromTaggedEpoch(te), getEpochFromTaggedEpoch(te));
            throw std::out_of_range("ERROR: epoch does not exist.");
          }
          return epoch_capacities.at(te) <= 0;
        }

        /// check if an epoch is active.
        bool isEpochPresent(const TaggedEpoch & te) {
          std::lock_guard<std::mutex> lock(mutex);
          return epoch_capacities.count(te) > 0;
        }

        /// wait for the capacity associated with an epoch to go down to 0. (i.e. epoch's communication is done).  uses condition variable
        void waitForEpochRelease(const TaggedEpoch & te) {
          std::unique_lock<std::mutex> lock(mutex);
          while (epoch_capacities.count(te) > 0) {
            condVar.wait(lock);
          }
          lock.unlock();
        }

        /// release an epoch.  notifies any waiting thread that epoch was released.
        void releaseEpoch(const TaggedEpoch & te, std::atomic<int>& totalCount) {
          std::unique_lock<std::mutex> lock(mutex);
          epoch_capacities.erase(te);
          totalCount.fetch_sub(1, std::memory_order_release);
          lock.unlock();

          condVar.notify_all();
        }

        /// get the Message Buffers class associated with this message type.
        std::shared_ptr<MessageBuffersType> getBuffer() {
          return std::shared_ptr<MessageBuffersType>(bufferPtr);
        }

    };


  } /* namespace io */
} /* namespace bliss */

#endif /* MESSAGE_TYPE_INFO_HPP_ */
