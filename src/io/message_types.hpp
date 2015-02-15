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

#include <io/message_type_info.hpp>

namespace bliss
{
  namespace io
  {

    //========= internal message queue data types.

    /**
     * base class for storing mpi tag and src for a received message.
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
     * @brief a MPI message representing a control message associated with a data message type.
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
    struct ControlMessage : public MPIMessage
    {
      /// The received data.  for control message, this is a tag + an epoch number
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
        : MPIMessage(CONTROL_TAG, rank), tagged_epoch(_control) {}

      /**
       * @brief Constructs a new instance of this struct and sets all members as given.
       *
       * @param _id   Id of MessageBuffer that will be sent
       * @param _tag  type of the message being sent
       * @param _epoch  epoch in which to send the message (epoch of the tag)
       * @param _dst  destination of the message
       */
      ControlMessage(int _tag, int _epoch, int rank)
        : MPIMessage(CONTROL_TAG, rank), tagged_epoch(static_cast<TaggedEpoch>(_tag) << 32 | static_cast<TaggedEpoch>(_epoch)) {}


      /// default constructor
      ControlMessage() = default;

      /// default destructor
      virtual ~ControlMessage() {};

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

      /// default destructor
      virtual ~DataMessageReceived() {};

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

      /// default destructor
      virtual ~DataMessageToSend() {};
    };



  } /* namespace io */
} /* namespace bliss */

#endif /* MESSAGE_TYPES_HPP_ */
