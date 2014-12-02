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



namespace bliss
{
  namespace io
  {

    /**
     * Metadata re. a message type (tag), used to manage whether communication with this message type is completed, the epoch number of the communication,
     * and lock/waits.
     */
    class MessageTypeInfo {
      public:

        using TaggedEpoch = int64_t;

        // ========= static variables for control messages.
        // unioned data structure is not friendly with both endianness and atomic operations.
        // since Tag+Epoch is used only for control messages, just use methods to extract tag and epoch when needed.
        // this also avoids defining out own hash and equal functions.

       /// constant indicating message type is control (FOC) message.
       static constexpr TaggedEpoch CONTROL_TAGGED_EPOCH = 0;
       static constexpr int32_t CONTROL_TAG = 0;

       /// constant indicating that there is no active tag. (e.g. during flushing)
       //static constexpr int32_t NO_TAG = -1;

       /// terminal epoch id.
       //static constexpr int32_t NO_EPOCH = -1;

       //==== type definitions



      protected:
        std::atomic<bool> finished;                    // do we need this?
        std::atomic<TaggedEpoch>  tagged_epoch;

        // TODO: static epoch id?

        std::mutex mutex;                    // do we need this?
        std::condition_variable condVar;           // do we need this?

        explicit MessageTypeInfo(const MessageTypeInfo& other) = delete;
        MessageTypeInfo& operator=(const MessageTypeInfo& other) = delete;

        MessageTypeInfo(MessageTypeInfo&& other, const std::lock_guard<std::mutex> &) {
          finished.exchange(other.finished.load());          other.finished = true;
          tagged_epoch.exchange(other.tagged_epoch.load());  other.tagged_epoch = -1;
          // have to set other to nullptr else destructor will clean up the reassigned buffer.
        }

      public:
        // needed by unordered_map during initial insert, before actual value is assigned.
        MessageTypeInfo() : finished(false), tagged_epoch(-1) {}

        MessageTypeInfo(const int32_t tag) :
          finished(false), tagged_epoch(static_cast<TaggedEpoch>(tag) << 32) {}

        explicit MessageTypeInfo(MessageTypeInfo&& other) :
            MessageTypeInfo(std::forward<MessageTypeInfo>(other), std::lock_guard<std::mutex>(other.mutex)) {}

        MessageTypeInfo& operator=(MessageTypeInfo&& other) {
          std::unique_lock<std::mutex> lock(mutex, std::defer_lock);
          std::unique_lock<std::mutex> otherlock(other.mutex, std::defer_lock);
          std::lock(lock, otherlock);

          finished.exchange(other.finished.load());              other.finished = true;
          // -1 in64 is 0xFFFFFFFFFFFFFFFF, which is -1, -1 for 2 ints.
          tagged_epoch.exchange(other.tagged_epoch.load());      other.tagged_epoch = -1;

          return *this;
        }

        virtual ~MessageTypeInfo() {}

        //===== accessors
        bool finish() {
          return finished.exchange(true, std::memory_order_seq_cst);
        }
        bool isFinished() {
          return finished.load(std::memory_order_seq_cst);
        }
        TaggedEpoch getNextEpoch() {
          return tagged_epoch.fetch_add(1);
        }
        //===== conveience methods
        static inline int32_t getTag(const TaggedEpoch & te) {
          return te >> 32;
        }


        //===== lock the data structure.
        // TODO: have to return a new lock object each time, else either deadlock with recursive mutex, or system_error trying to relock already locked object.
        std::unique_lock<std::mutex> lock() {
          return std::move(std::unique_lock<std::mutex>(mutex));
        }
        void unlock(std::unique_lock<std::mutex>&& l) {
          l.unlock();
        }

        //===== condition variable handling.  this can be shared, but then we need to pass in the lock.
        void wait(std::unique_lock<std::mutex>& l) {
          condVar.wait(l);
        }
        void notifyAll() {
          condVar.notify_all();
        }
        void notifyOne() {
          condVar.notify_one();
        }

    };


    //==== static variable definitions.

    /// CONTROL_TAG definition.  (declaration and initialization inside class).
    constexpr MessageTypeInfo::TaggedEpoch MessageTypeInfo::CONTROL_TAGGED_EPOCH;

    constexpr int32_t MessageTypeInfo::CONTROL_TAG;

//    /// NO_TAG definition.  (declaration and initialization inside class).
//    template <typename MessageBuffersType>
//    constexpr int MessageTypeInfo<MessageBuffersType>::NO_TAG;
//
//    /// NO_EPOCH definition.  (declaration and initialization inside class).
//    template <typename MessageBuffersType>
//    constexpr int MessageTypeInfo<MessageBuffersType>::NO_EPOCH;

  } /* namespace io */
} /* namespace bliss */

#endif /* MESSAGE_TYPE_INFO_HPP_ */
