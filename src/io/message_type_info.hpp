/**
 * @file    message_type_info.hpp
 * @ingroup
 * @author  tpan
 * @brief
 * @details
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

    typedef int64_t TaggedEpoch;
    inline int32_t getTagFromTaggedEpoch(const TaggedEpoch & te) {
      return te >> 32;
    }
    inline uint32_t getEpochFromTaggedEpoch(const TaggedEpoch & te) {
      return te & 0x00000000FFFFFFFF;
    }


    constexpr TaggedEpoch CONTROL_TAGGED_EPOCH = 0;
    constexpr int32_t CONTROL_TAG = 0;

    /**
     * Metadata re. a message type (tag), used to manage whether communication with this message type is completed, the epoch number of the communication,
     * and lock/waits.
     */
    template<typename MessageBuffersType>
    class MessageTypeInfo {
      public:

        // ========= static variables for control messages.
        // unioned data structure is not friendly with both endianness and atomic operations.
        // since Tag+Epoch is used only for control messages, just use methods to extract tag and epoch when needed.
        // this also avoids defining out own hash and equal functions.

       /// constant indicating message type is control (FOC) message.


       /// constant indicating that there is no active tag. (e.g. during flushing)
       //static constexpr int32_t NO_TAG = -1;

       /// terminal epoch id.
       //static constexpr int32_t NO_EPOCH = -1;

       //==== type definitions

      protected:
        std::shared_ptr<MessageBuffersType> bufferPtr;

        std::atomic<bool>         finished;  // do we need this?
        std::atomic<TaggedEpoch>  nextEpoch;

        mutable int maxEpochCapacity;

        /// map from epoch to count
        std::unordered_map<TaggedEpoch, int > epoch_capacities;

        mutable std::mutex mutex;                    // for locking access to epoch capacities
        mutable std::condition_variable condVar;     // for wait.

        explicit MessageTypeInfo(const MessageTypeInfo& other) = delete;
        MessageTypeInfo& operator=(const MessageTypeInfo& other) = delete;

          MessageTypeInfo(MessageTypeInfo&& other, const std::lock_guard<std::mutex> &) :
            bufferPtr(std::move(other.bufferPtr)),
            maxEpochCapacity(other.maxEpochCapacity),
            epoch_capacities(std::move(other.epoch_capacities)) {

            finished.store(other.finished.exchange(true, std::memory_order_relaxed));
            nextEpoch.store(other.nextEpoch.exchange(-1, std::memory_order_relaxed));

            // have to set other to nullptr else destructor will clean up the reassigned buffer.
            other.maxEpochCapacity = 0;
          }

//        MessageTypeInfo(MessageTypeInfo<MessageBuffersType>&& other) = delete;
//        MessageTypeInfo<MessageBuffersType>& operator=(MessageTypeInfo<MessageBuffersType>&& other) = delete;



      public:
        // needed by unordered_map during initial insert, before actual value is assigned.
        MessageTypeInfo() : finished(false), nextEpoch(-1), maxEpochCapacity(0) {}

        MessageTypeInfo(const int32_t tag, const std::shared_ptr<MessageBuffersType>& ptr, const int epoch_capacity) :
          bufferPtr(ptr), finished(false), nextEpoch(static_cast<TaggedEpoch>(tag) << 32), maxEpochCapacity(epoch_capacity) {}

        explicit MessageTypeInfo(MessageTypeInfo&& other) :
            MessageTypeInfo(std::forward<MessageTypeInfo>(other), std::lock_guard<std::mutex>(other.mutex)) {}

        MessageTypeInfo& operator=(MessageTypeInfo&& other) {
          std::unique_lock<std::mutex> lock(mutex, std::defer_lock);
          std::unique_lock<std::mutex> otherlock(other.mutex, std::defer_lock);
          std::lock(lock, otherlock);

          bufferPtr = std::move(other.bufferPtr);

          finished.store(other.finished.exchange(true, std::memory_order_relaxed));
          // -1 in64 is 0xFFFFFFFFFFFFFFFF, which is -1, -1 for 2 ints.
          nextEpoch.store(other.nextEpoch.exchange(-1, std::memory_order_relaxed));
          maxEpochCapacity = other.maxEpochCapacity; other.maxEpochCapacity = 0;
          epoch_capacities = std::move(other.epoch_capacities);

          return *this;
        }



        virtual ~MessageTypeInfo() {}

        //===== accessors
        bool finish() {
          return finished.exchange(true, std::memory_order_acq_rel);
        }
        bool isFinished() {
          return finished.load(std::memory_order_relaxed);
        }


        TaggedEpoch acquireEpoch() {
          TaggedEpoch te = nextEpoch.fetch_add(1, std::memory_order_relaxed);

          std::lock_guard<std::mutex> lock(mutex);
          epoch_capacities[te] = maxEpochCapacity;
          return te;
        }
        int countdownEpoch(const TaggedEpoch & te) throw (std::out_of_range) {
          std::lock_guard<std::mutex> lock(mutex);
          if (epoch_capacities.count(te) == 0) {
            epoch_capacities[te] = maxEpochCapacity;
//            ERRORF("epoch does not exist to count down: %d %d", getTagFromTaggedEpoch(te), getEpochFromTaggedEpoch(te));
//            throw std::out_of_range("ERROR: epoch does not exist.");
          }
          return --epoch_capacities.at(te);
        }
        bool isEpochFinished(const TaggedEpoch & te) throw (std::out_of_range) {
          std::lock_guard<std::mutex> lock(mutex);
          if (epoch_capacities.count(te) == 0) {
            ERRORF("epoch does not exist to check for Finished: %d %d", getTagFromTaggedEpoch(te), getEpochFromTaggedEpoch(te));
            throw std::out_of_range("ERROR: epoch does not exist.");
          }
          return epoch_capacities.at(te) <= 0;
        }
        bool isEpochPresent(const TaggedEpoch & te) {
          std::lock_guard<std::mutex> lock(mutex);
          return epoch_capacities.count(te) > 0;
        }
        void waitForEpochRelease(const TaggedEpoch & te) {
          std::unique_lock<std::mutex> lock(mutex);
          while (epoch_capacities.count(te) > 0) {
            condVar.wait(lock);
          }
          lock.unlock();
        }
        void releaseEpoch(const TaggedEpoch & te) {
          std::unique_lock<std::mutex> lock(mutex);
          epoch_capacities.erase(te);
          lock.unlock();

          condVar.notify_all();
        }


//        //===== conveience methods
//        static inline int32_t getTag(const TaggedEpoch & te) {
//          return te >> 32;
//        }
//
        std::shared_ptr<MessageBuffersType> getBuffer() {
          return std::shared_ptr<MessageBuffersType>(bufferPtr);
        }


//        //===== lock the data structure.
//        // TODO: have to return a new lock object each time, else either deadlock with recursive mutex, or system_error trying to relock already locked object.
//        std::unique_lock<std::mutex> lock() {
//          return std::move(std::unique_lock<std::mutex>(mutex));
//        }
//        void unlock(std::unique_lock<std::mutex>&& l) {
//          l.unlock();
//        }
//
//        //===== condition variable handling.  this can be shared, but then we need to pass in the lock.
//        void wait(std::unique_lock<std::mutex>& l) {
//          condVar.wait(l);
//        }
//        void notifyAll() {
//          condVar.notify_all();
//        }
//        void notifyOne() {
//          condVar.notify_one();
//        }

    };


    //==== static variable definitions.

//    /// CONTROL_TAG definition.  (declaration and initialization inside class).
//    template<typename MessageBuffersType>
//    constexpr typename MessageTypeInfo<MessageBuffersType>::TaggedEpoch MessageTypeInfo<MessageBuffersType>::CONTROL_TAGGED_EPOCH;
//
//    template<typename MessageBuffersType>
//    constexpr int32_t MessageTypeInfo<MessageBuffersType>::CONTROL_TAG;


  } /* namespace io */
} /* namespace bliss */

#endif /* MESSAGE_TYPE_INFO_HPP_ */
