/**
 * @file    message_types.hpp
 * @ingroup bliss::io
 * @author  Tony Pan <tpan7@gatech.edu>
 * @brief   various types of mpi messages:  sent data, received data, control messages, base class
 * @details  background:
 *          For control message, the tag sent/received under is the control tag.
 *            its payload has 2 parts:  tag to control, and the epoch of the communication.
 *              tag to control is associated to a type of data message.
 *              epoch is a unique number for that tag to indicate the particular communication episode for that message type
 *
 *          For Send data message, the tag is the message type.  a pointer is stored, but space not managed.
 *          For Recv data message, the tag is the message type.  a pointer is stored, and space is managed.
 *
 *        Here we need to rely on non-overtaking behavior of MPI to preserve order between data messages (tag > 0) and control messages (tag == 0).
 *
 *
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef MESSAGE_TYPES_HPP_
#define MESSAGE_TYPES_HPP_

#include <io/mpi_utils.hpp>


namespace bliss
{
  namespace io
  {

    //========= internal message queue data types.


    /**
     * base class for storing mpi tag and src for a received message.
     */
    class MPIMessage
    {
      public:
        /// constant indicating that all MPI ranks are involved.
        static const int ALL_RANKS = -1;
        static const int ALL_TAGS = -1;
        /// static atomic for tracking the current/next epoch in a process.
        static std::atomic<uint64_t> nextEpoch;

      protected:
        /// The message tag and epoch. tag indicates the type of message (how to interpret the message, control vs data)
        /// epoch indicates the phase, which is global.
        int tag;

        /// The message source id
        int rank;

      public:

        /**
         * @brief constructor using a pre-existing memory block
         *
         * @param tag       the MPI message tag, indicating the type of message
         * @param src       the MPI source process rank
         */
        MPIMessage(int _tag, int _rank)
          : tag(_tag), rank(_rank) {}
 
        /// default constructor
        MPIMessage() = default;

        /// needed to create a virtual function table, only then is polymorphism allowed.  (inheritance not sufficient)
        virtual ~MPIMessage() {};

        /// get the message's tag
        inline int getTag() { return tag; };

        /// get the messages's target process rank
        inline int getRank() { return rank; };

        /// get the message's epoch
        virtual uint64_t getEpoch() = 0;

        /// get the payload (Data + metadata)
        virtual uint8_t* getPayload() = 0;

        /// get the payload (data + metadata) size in bytes
        virtual size_t getPayloadSize() = 0;

        /// get the payload
        virtual uint8_t* getData() = 0;

        /// get the payload size in bytes
        virtual size_t getDataSize() = 0;
    };
    std::atomic<uint64_t> MPIMessage::nextEpoch(0);

    /**
     * @brief a MPI message representing a control message.  control messages are NOT associated with a message type.   It is used to mark the start and end of communication.
     * @details tag is used to annotate the message type/format, here it would be CONTROL_TAG.
     *        The payload has a tag component, interpreted as the data message type being controlled.
     *        The epoch component indicates the communication episode this control message is associated with
     *
     *        The control messages are used to flush the receiver's queue of a particular message type.  "flush" is considered
     *        a collective operation, with synchronization barrier satisfied when control messages for the same message type are received from all sources.
     *
     *        if the payload's tag inforation is CONTROL_TAG, then it indicates that the calling application is
     *          completely finished with communication.
     *
     */
    class ControlMessage : public MPIMessage
    {
      protected:
        /// current epoch [0] and next epoch [1].  these are set externally instead of using the atomic MPIMessage::nextEpoch so we can use ControlMessage
        /// to flush multiple tags.
        uint64_t epochs[2];
      public:

      /**
       * @brief Constructs a new instance of this struct and sets all members as given.
       *
       * @param curr_epoch  epoch in which to send the message
       * @param next_epoch  the next epoch.
       * @param _dst  destination of the message
       */
      ControlMessage(uint64_t curr_epoch, uint64_t next_epoch, int _dst)
        : MPIMessage(bliss::io::CONTROL_TAG, _dst) {
        epochs[0] = curr_epoch;
        epochs[1] = next_epoch;
      }


      /// default constructor
      ControlMessage() = default;

      /// default destructor
      virtual ~ControlMessage() {};

      /// get the current epoch
      virtual uint64_t getEpoch() { return epochs[0]; };

      /// get the next epoch
      uint64_t getNextEpoch() { return epochs[1]; };

      /// get pointer to the payload (data + metadata)
      virtual uint8_t* getPayload() {
        return reinterpret_cast<uint8_t*>(epochs);
      }

      /// get size in bytes of the payload (data + metadata)
      virtual size_t getPayloadSize() {
        return sizeof(uint64_t) * 2;
      }


      /// get pointer to the data
      virtual uint8_t* getData() {
        return reinterpret_cast<uint8_t*>(epochs);
      }

      /// get size in bytes of the data.
      virtual size_t getDataSize() {
        return sizeof(uint64_t) * 2;
      }

    };


    /**
     * @brief Structure to hold all associated information of a received MPI message.
     * @details tag is used to annotate the message type/format.  The message is processed according to this tag.
     *
     *        This class only handles received DATA MESSAGES.
     *
     *        On construction, this class is assigned a pointer to a byte array.  From this point on, this class
     *        manages that byte array.
     */
    class DataMessageReceived  : public MPIMessage
    {
      protected:
        /// The received data
        std::unique_ptr<uint8_t[]> data;

        /// The number of bytes received
        std::size_t count;
      
      public:
      /**
       * @brief constructor using a pre-existing memory block
       *
       * @param data      in memory block of data (received from remote proc, to be processed here).  unique_ptr
       * @param count     number of bytes in the data block
       * @param tag       the MPI message tag, indicating the type of message
       * @param src       the MPI source process rank
       */
      DataMessageReceived(std::unique_ptr<uint8_t[]>&& _data, std::size_t count, int tag, int src)
        : MPIMessage(tag, src), data(std::move(_data)), count(count) {}

      /**
       * @brief constructor using a pre-existing memory block
       *
       * @param raw_data  in memory block of data (received from remote proc, to be processed here).  raw pointer
       * @param count     number of bytes in the data block
       * @param tag       the MPI message tag, indicating the type of message
       * @param src       the MPI source process rank
       */
      DataMessageReceived(uint8_t* raw_data, std::size_t count, int tag, int src)
        : MPIMessage(tag, src), data(raw_data), count(count) {}


      /// default constructor
      DataMessageReceived() = default;

      /// default destructor
      virtual ~DataMessageReceived() {};

      /// get the message's epoch embedded as the first sizeof(uint64_t) bytes.
      virtual uint64_t getEpoch() {
        return (reinterpret_cast<uint64_t*>(data.get()))[0];
      }

      /// get pointer to the payload (data + metadata)
      virtual uint8_t* getPayload() {
        return data.get();
      }

      /// get size in bytes of the payload (data + metadata)
      virtual size_t getPayloadSize() {
        return count;
      }


      /// get pointer to data, exclude the metadata portion  (for whole data, access member directly
      virtual uint8_t* getData() { return (data == nullptr) ? nullptr : (data.get() + sizeof(uint64_t)); }

      /// get size of data, exclude the metadata portion
      virtual size_t getDataSize() { return (count == 0) ? 0 : (count - sizeof(uint64_t)); }
    };


    /**
     * @brief Structure to hold all associated information of a message to be sent via MPI.
     * @details MessageToSend is a small metadata stucture to be used for message queuing (to be sent)
     *          MessageToSend does not actually hold any data.  it points to a in-memory Buffer by pointer.
     *          It does NOT manage data.
     *
     *          control message is NOT sent this way, since there is no Buffer associated with a control message.
     *
     * @tparam  BufferPtrType  type of Buffer id.  Obtained from SendMessageBuffer, parameterized by Thread Safety.
     */
    template<typename BufferPtrType>
    struct DataMessageToSend : public MPIMessage
    {
      protected:
        /// pointer to the message buffer
        BufferPtrType data;

      public:
      /**
       * @brief Constructs a new instance of this struct and sets all members as given.
       *
       * @param _id   Id of MessageBuffer that will be sent
       * @param _tag  type of the message being sent
       * @param _dst  destination of the message
       */
      DataMessageToSend(BufferPtrType _ptr, int _tag, int _dst)
        : MPIMessage(_tag, _dst), data(_ptr) {}


      /// default constructor
      DataMessageToSend() = default;

      /// default destructor
      virtual ~DataMessageToSend() {};

      /// get the message's epoch embedded as the first sizeof(uint64_t) bytes.
      virtual uint64_t getEpoch() {
        return ((uint64_t*)data)[0];
      }


      /// get pointer to the payload (data + metadata)
      virtual uint8_t* getPayload() {
        return (data == nullptr) ? nullptr : data->operator uint8_t*();
      }

      /// get size in bytes of the payload (data + metadata)
      virtual size_t getPayloadSize() {
        return (data == nullptr) ? 0 : data->getSize();
      }


      /// get pointer to data, exclude the metadata portion  (for whole data, access member directly
      virtual uint8_t* getData() { return (getPayload() == nullptr) ? nullptr : (getPayload() + sizeof(uint64_t)); }

      /// get size of data, exclude the metadata portion
      virtual size_t getDataSize() { return (data == nullptr) ? 0 : (data->getSize() - sizeof(uint64_t)); }


      BufferPtrType& getBuffer() { return data; }

    };



  } /* namespace io */
} /* namespace bliss */

#endif /* MESSAGE_TYPES_HPP_ */
