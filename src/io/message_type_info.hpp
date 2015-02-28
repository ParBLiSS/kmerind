/**
 * @file    message_type_info.hpp
 * @ingroup bliss::io
 * @author  Tony Pan <tpan7@gatech.edu>
 * @brief   a data structure to manage data and control message communications associated with a message type.
 * @details defines MessageTypeInfo for presenting tag and id of distinct communication episodes
 *          defines EpochInfo to track active epochs for all Tag.
 *
 *
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
#include "io/mpi_utils.hpp"

namespace bliss
{
  namespace io
  {


    /**
     * @class MessageTypeInfo
     * @brief Used to manage message buffers for a message type.
     * @details
     *          This class manages a pointer to the current MessageBuffers for a message type
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

        /// tag == message type
        int tag;

        /// deleted copy constructor
        explicit MessageTypeInfo(const MessageTypeInfo& other) = delete;

        /// deleted copy assignment operator
        MessageTypeInfo& operator=(const MessageTypeInfo& other) = delete;

      public:
        /// default constructor.  needed by some containers.
        MessageTypeInfo() : finished(false), tag(-1) {}

        /// constructor
        MessageTypeInfo(const int32_t _tag, const std::shared_ptr<MessageBuffersType>& ptr) :
          bufferPtr(ptr), finished(false), tag(_tag) {}

        /// move constructor
        MessageTypeInfo(MessageTypeInfo&& other) :
          bufferPtr(std::move(other.bufferPtr)), tag(other.tag)
        {
          finished.store(other.finished.exchange(true, std::memory_order_relaxed));
          other.tag = -1;
        }

        /// move assignment operator
        MessageTypeInfo& operator=(MessageTypeInfo&& other) {
          bufferPtr = std::move(other.bufferPtr);
          finished.store(other.finished.exchange(true, std::memory_order_relaxed));

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

        /// get the Message Buffers class associated with this message type.
        std::shared_ptr<MessageBuffersType> getBuffers() {
          return std::shared_ptr<MessageBuffersType>(bufferPtr);
        }

    };

    /**
     * @class EpochInfo
     * @brief Used to manage active communication episodes for all message types.
     * @details
     *          This class handles registering an epoch number, store the epoch number
     *          count of control messages with matching tag/epoch for each epoch.
     *
     *          It provides thread safe access to the epoch-to-count mapping,
     *          and provides a mechanism for calling thread to wait for the count to
     *          reach 0 for a TaggedEpoch.  This is done via mutex and condition variables.
     *
     */
    class EpochInfo {
      public:

      protected:
        /// number of processes involved in communication epoch
        mutable int commSize;

        /// map from epoch to count to track active epochs.
        std::unordered_map<uint64_t, int > epoch_capacities;

        /// for locking access to epoch_capacities
        mutable std::mutex mutex;
        /// for use to wait for an epoch to finish
        mutable std::condition_variable condVar;

        /// deleted copy constructor
        explicit EpochInfo(const EpochInfo& other) = delete;

        /// deleted copy assignment operator
        EpochInfo& operator=(const EpochInfo& other) = delete;

        /// mutex locked move constructor.
        EpochInfo(EpochInfo&& other, const std::lock_guard<std::mutex> &) :
          commSize(other.commSize),
          epoch_capacities(std::move(other.epoch_capacities)) {

          other.commSize = 0;
        }


      public:
        /// default constructor.  needed by some containers.
        EpochInfo() : commSize(0) {}

        /// constructor
        EpochInfo(const int epoch_capacity) : commSize(epoch_capacity) {}

        /// move constructor
        explicit EpochInfo(EpochInfo&& other) :
            EpochInfo(std::forward<EpochInfo>(other), std::lock_guard<std::mutex>(other.mutex)) {}

        /// move assignment operator
        EpochInfo& operator=(EpochInfo&& other) {
          std::unique_lock<std::mutex> lock(mutex, std::defer_lock);
          std::unique_lock<std::mutex> otherlock(other.mutex, std::defer_lock);
          std::lock(lock, otherlock);

          commSize = other.commSize; other.commSize = 0;
          epoch_capacities = std::move(other.epoch_capacities);

          return *this;
        }


        /// default constructor
        virtual ~EpochInfo() {
        }

        std::string toString() {
          std::lock_guard<std::mutex> lock(mutex);
          std::stringstream ss;
          for (auto el : epoch_capacities) {
            ss << "(" << el.first << "=" << el.second << "),";
          }
          return ss.str();
        }


        /**
         * @brief  get the next epoch id, and insert the epoch into the epoch-capacity map.
         * @detail   note that the mapping may be inserted before the epoch was created.
         *        this can happen if a remote process sent out a control message with that epoch before local process acquired.
         */
        int registerEpoch(const uint64_t & epoch, std::atomic<int>& totalCount) {

          std::lock_guard<std::mutex> lock(mutex);
          if (epoch_capacities.count(epoch) == 0) {
            epoch_capacities[epoch] = commSize;
            totalCount.fetch_add(1, std::memory_order_release);
          }
          return commSize;
        }
        /**
         * @brief  reduce the capacity associated with an epoch by 1.
         * @detail   note that the mapping may not exist yet when a control MPI message is received,
         *        as the local process not yet acquired the epoch.
         *        In this case, a new mapping is created then decremented.
         *        this can happen if a remote process sent out a control message with that epoch before local process acquired.
         */
        int countdownEpoch(const uint64_t & epoch, std::atomic<int>& totalCount, int decrement = 1) throw (std::out_of_range) {
          std::lock_guard<std::mutex> lock(mutex);
          if (epoch_capacities.count(epoch) == 0) {
            epoch_capacities[epoch] = commSize;
            totalCount.fetch_add(decrement, std::memory_order_release);
            DEBUGF("epoch does not exist to count down. registered %lu", epoch);
          }
          return epoch_capacities.at(epoch) -= decrement;
        }

        /// check if an epoch is finished.
        bool isEpochFinished(const uint64_t & epoch) throw (std::out_of_range) {
          std::lock_guard<std::mutex> lock(mutex);
          if (epoch_capacities.count(epoch) == 0) {
            ERRORF("epoch does not exist to check for Finished: %lu", epoch);
            throw std::out_of_range("ERROR: epoch does not exist.");
          }
          return epoch_capacities.at(epoch) <= 0;
        }

        /// check if an epoch is active.
        bool isEpochPresent(const uint64_t & epoch) {
          std::lock_guard<std::mutex> lock(mutex);
          return epoch_capacities.count(epoch) > 0;
        }

        /// wait for the capacity associated with an epoch to go down to 0. (i.e. epoch's communication is done).  uses condition variable
        void waitForEpochRelease(const uint64_t & epoch) {
          std::unique_lock<std::mutex> lock(mutex);
          while (epoch_capacities.count(epoch) > 0) {
            condVar.wait(lock);
          }
          lock.unlock();
        }

        /// release an epoch.  notifies any waiting thread that epoch was released.
        void releaseEpoch(const uint64_t & epoch, std::atomic<int>& totalCount) {
          std::unique_lock<std::mutex> lock(mutex);
          epoch_capacities.erase(epoch);
          totalCount.fetch_sub(1, std::memory_order_release);
          lock.unlock();

          condVar.notify_all();
        }

    };

  } /* namespace io */
} /* namespace bliss */

#endif /* MESSAGE_TYPE_INFO_HPP_ */
