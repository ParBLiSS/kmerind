/**
 * @file    message_types.hpp
 * @ingroup
 * @author  tpan
 * @brief   various types of mpi messages:  sent data, received data, control messages, base class
 * @details
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef MESSAGE_TYPES_HPP_
#define MESSAGE_TYPES_HPP_

#include <io/message_type_info.hpp>

namespace bliss
{
  namespace io
  {

    //========= internal message queue data types.

    /**
     * base class for storing mpi tag and sr for a received message.
     */
    struct MPIMessage
    {
        /// The message tag, indicates the type of message (how to interpret the message)
        int tag;

        /// The message source id
        int rank;

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
    };


    /**
     * @brief Structure to hold all associated information of a received MPI message.
     * @details tag is used to annotate the message type/format.  The message is processed according to this tag.
     *
     *        tag == CONTROL_TAG:  control message indicating the last of a type of data messages for a synchronization epoch (period between synchronizations)
     *                    data contains the tag for the data message type that is ending, as well as an identifier for a synchronization epoch.
     *                    if data contains a tag of 0, then the application is ending.
     *        otherwise:  data messages, with data containing payload.
     *
     *        Here we need to rely on non-overtaking behavior of MPI to preserve order between data messages (tag > 0) and control messages (tag == 0).
     *
     *        The control messages are used to flush the receiver's queue of a particular message type.  "flush" is considered
     *        a collective operation, with synchronization barrier satisfied when control messages for the same message type are received from all sources.
     *
     */
    struct ControlMessage : public MPIMessage
    {
       using MessageTypeInfo = bliss::io::MessageTypeInfo;
       using TaggedEpoch = typename MessageTypeInfo::TaggedEpoch;
      /// The received data
       TaggedEpoch tagged_epoch;

      /**
       * @brief constructor using a pre-existing memory block
       *
       * @param data      in memory block of data (received from remote proc, to be processed here)
       * @param count     number of bytes in the data block
       * @param tag       the MPI message tag, indicating the type of message
       * @param src       the MPI source process rank
       */
      ControlMessage(TaggedEpoch _control, int rank)
        : MPIMessage(MessageTypeInfo::CONTROL_TAG, rank), tagged_epoch(_control) {}

      /**
       * @brief Constructs a new instance of this struct and sets all members as given.
       *
       * @param _id   Id of MessageBuffer that will be sent
       * @param _tag  type of the message being sent
       * @param _epoch  epoch in which to send the message (epoch of the tag)
       * @param _dst  destination of the message
       */
      ControlMessage(int _tag, int _epoch, int rank)
        : MPIMessage(MessageTypeInfo::CONTROL_TAG, rank), tagged_epoch(static_cast<TaggedEpoch>(_tag) << 32 | static_cast<TaggedEpoch>(_epoch)) {}


      /// default constructor
      ControlMessage() = default;

      virtual ~ControlMessage() {};

    };


    /**
     * @brief Structure to hold all associated information of a received MPI message.
     * @details tag is used to annotate the message type/format.  The message is processed according to this tag.
     *
     *        tag == CONTROL_TAG:  control message indicating the last of a type of data messages for a synchronization epoch (period between synchronizations)
     *                    data contains the tag for the data message type that is ending, as well as an identifier for a synchronization epoch.
     *                    if data contains a tag of 0, then the application is ending.
     *        otherwise:  data messages, with data containing payload.
     *
     *        Here we need to rely on non-overtaking behavior of MPI to preserve order between data messages (tag > 0) and control messages (tag == 0).
     *
     *        The control messages are used to flush the receiver's queue of a particular message type.  "flush" is considered
     *        a collective operation, with synchronization barrier satisfied when control messages for the same message type are received from all sources.
     *
     *        data is a unique pointer to a byte array.  received data buffer is managed by DataMessageRecevied.
     */
    struct DataMessageReceived  : public MPIMessage
    {
      /// The received data
      std::unique_ptr<uint8_t[]> data;

      /// The number of bytes received
      std::size_t count;

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

      virtual ~DataMessageReceived() {};

    };


    /**
     * @brief Structure to hold all associated information of a message to be sent via MPI.
     * @details MessageToSend is a small metadata stucture to be used for message queuing (to be sent)
     *          MessageToSend does not actually hold any data.  it points to a in-memory Buffer by pointer.
     *
     *          Does NOT manage data.
     *
     *          control message is NOT sent this way, since there is no Buffer associated with a control message.
     *
     * @tparam  BufferPtrType  type of Buffer id.  Obtained from SendMessageBuffer, parameterized by Thread Safety.
     */
    template<typename BufferPtrType>
    struct DataMessageToSend : public MPIMessage
    {
      /// pointer to the message buffer
      BufferPtrType ptr;

      /**
       * @brief Constructs a new instance of this struct and sets all members as given.
       *
       * @param _id   Id of MessageBuffer that will be sent
       * @param _tag  type of the message being sent
       * @param _dst  destination of the message
       */
      DataMessageToSend(BufferPtrType _ptr, int _tag, int _dst)
        : MPIMessage(_tag, _dst), ptr(_ptr) {}


      /// default constructor
      DataMessageToSend() = default;

      virtual ~DataMessageToSend() {};
    };



  } /* namespace io */
} /* namespace bliss */

#endif /* MESSAGE_TYPES_HPP_ */
