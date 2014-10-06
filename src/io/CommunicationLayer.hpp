/**
 * @file    CommunicationLayer.hpp
 * @ingroup bliss::io
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @author  Tony Pan <tpan@gatech.edu>
 * @brief   contains a class that abstracts MPI point-to-point, non-blocking messaging.
 *
 * Copyright (c) TODO
 *
 * TODO add Licence
 */

#ifndef BLISS_COMMUNICATION_LAYER_HPP
#define BLISS_COMMUNICATION_LAYER_HPP

#include <mpi.h>
#include <omp.h>

// C stdlib includes
#include <assert.h>
#include <cstdio>  // fflush

// STL includes
#include <queue>
#include <mutex>
#include <unordered_map>
#include <unordered_set>
#include <thread>
#include <utility>    // std::move

// BLISS includes
#include <utils/logging.h>
#include <concurrent/threadsafe_queue.hpp>
#include <concurrent/concurrent.hpp>
#include <io/message_buffers.hpp>
#include <io/io_exception.hpp>


namespace bliss
{
namespace io
{


/**
 * @brief Abstracts asynchronous and buffered point to point communication between all processes via MPI.
 * @details   This class encapsulates the MPI message handling for point to point asynchronous messaging
 *        for a MPI communicator.  It supports multiple concurrent message producers and consumers
 *
 *        The following requirements are satisfied:
 *        1. support application where there are different message types (e.g. query-response).
 *        2. support collectively, blocking flushing all messages of a particular type.
 *        3. provide flexible way of handling received messages based on message type
 *        4. Thread Safe, non-blocking send (allow overlap of computation with communication)
 *        5. Thread Safe processing of received messages
 *        6. single threaded control of the Communication Layer (flush messages, finish, etc)
 *
 *        Design and Implementation  (MORE DETAILS BELOW)
 *        1. support application where there are different message types (e.g. query-response).
 *          each data message type is assigned an id ranging from 1 to MAX_INT.
 *          FOC messages have id 0, but carries in payload the id of the message type that is being controlled.
 *
 *        2. support collectively, blocking flushing all messages of a particular type.
 *          To indicate to the system that there is no further messages of a particular type, collective flush
 *          operations have been implemented.  There are 2 types:
 *          a. "flush":  flush a single message type
 *          b. "finish": flush all message types (and mark end of application)
 *
 *          All processes need to call these functions together.  When called, the function blocks the calling thread
 *          until all processes have completely processed all matching FOC messages from all senders
 *
 *        3. provide flexible way of handling received messages based on message type
 *          Callbacks can be registered, one per message type.  When processing received messages,
 *          the matching callback is used.
 *
 *        4. Thread Safe, non-blocking send (allow overlap of computation with communication)
 *          A compute thread can call sendMessage on the Communication Layer object.  each invocation
 *          atomically gets a location to insert the message into a batching buffer via CAS operation.
 *          (thread safe insert)
 *
 *          When a buffer is full and sent, it is inserted into a thread-safe send queue for the comm thread to
 *          send via MPI  (non-blocking).
 *
 *        5. Thread Safe processing of received messages
 *          Multiple receive threads can process the received messages in the receive queue, which is thread safe.
 *
 *        6. single threaded control of the Communication Layer (flush messages, finish, etc)
 *          flush and finish functions are collectively called by the control thread (owner of Communication Layer object)
 *
 *
 *        Basic Concepts and Design:
 *        A. Threading
 *          The following are the intended threading pattern:
 *          * single control thread (owner of Communication Layer instance.  can be one of the compute threads)
 *          * one or more compute threads (producer of messages to be sent)
 *          * single MPI communication thread (internal)
 *          * one or more receive threads (consumers of messages received, internal)
 *
 *        The choice of using a single MPI communication thread is based on
 *          1. Single threaded MPI libraries are more readily available on different platforms.
 *          2. Single threaded MPI communication enforces message ordering (especially when sending control messages)
 *          3. Single threaded MPI communication reduce opportunities for deadlock due to thread timing.
 *
 *        The single MPI comm thread and multiple producer/consumer pattern requires thread-safe queues as a mechanism
 *          for inter-thread communication.
 *
 *        To avoid race condition amongst multiple compute threads or receive threads:
 *          1. Compute threads insert messages into Buffers (bliss::io::Buffer), which is thread safe for insertion by using CAS operations to get position for insert.
 *          2. When a buffer is full, compute thread enqueues it for MPI send, and a empty buffer is swapped in from BufferPool.  BufferPool is also thread safe.
 *          3. When multiple receive threads dequeue to process received messages buffers, they do so from a thread safe queue.
 *
 *
 *        B. Data Types
 *          Message: a message has a payload and a type id.
 *
 *          There are 2 families of messages:  data and FOC messages.
 *          Data Messages have actual data as its payload, and the message type id as its MPI tag.
 *          FOC messages have 2 integers as its payload: type id of the message to be controlled, and an epoch id
 *            FOC messages have 0 as its message type.
 *
 *          For why epoch id is needed, see "Rationales"
 *
 *        C. MPI Messaging
 *          Data messages are transmitted using MPI non-blocking Point to Point communication: isend, iprobe + irecv
 *            due to the non-blocking calls, also use queue to manage send and receive in progress.
 *          FOC messages are transmitted using MPI bsend to simplify memory management.
 *          Local messages (where src and destination ranks are the same) are sent via memory (direct insert into receive queue)
 *
 *          Collective operations (flush, finish; for requirement 2) are implemented using Point to Point communication of FOC messages.
 *          The messaging pattern is peer-to-peer (no master-slave relationship).
 *            as a consequence, collective operations require all to send to all (in a non-blocking way)
 *            TODO: perhaps use master-slave, with rotating master based on epoch id.  reduces comm from p^2 for each epoch to 2p.
 *
 *          flush requires:
 *            1. called collectively by all control threads
 *            2. blocks until the process has received AND PROCESSED the FOC message from ALL senders. (to ensure collectiveness)
 *                a. As message ordering is strict, all messages prior to the FOC messages, including all data messages of the specified
 *                  type, are required to be processed as well.
 *          Flush and finish therefore serve as synchronization points in the application logic, and can also be viewed as MPI messaging fences.
 *          The synchronization point defines EPOCHs, which are periods between synchronization calls.  each epoch is assigned an id.
 *
 *          Typical flow of control for sending data is:
 *            1. compute thread adds message to CommunicationLayer, which is added to a buffer.
 *                buffer is enqueued in send-queue when full.)
 *            2. MPI thread pops from send queue and initiates MPI isend
 *              a. isend requests are queued in 'sendInProgress' queue, and checked for completion.
 *            3. MPI thread probes for data.  if found,
 *              a. buffer allocated and irecv initiated, with request pushed into recvInProgress queue.
 *              b. on completion of receive, data is moved to receive queue.
 *            4. receive threads pops from receive queue and processed
 *
 *          Typical flow of control for flushing:
 *            1. control thread calls flush or finish, which enqueues FOC messages, then blocks
 *            2. MPI thread pops from queue and initiates blocking buffered send (bsend)
 *            3. MPI thread probes for data.  if FOC message found,
 *              a. buffer allocated and irecv initiated, with request pushed into recvInProgress queue.
 *              b. on completion of receive, decrement the FOC message count for the specified epoch.
 *                when 0 (received from all senders), an END message is inserted into receive queue.
 *            4. when FOC messages are popped from receive queue, receive thread unblocks control thread's flush call.
 *
 *          Rationale for NOT using MPI collective functions: see Rationales below.
 *
 *          For more detail, see Rationale and Design and Implementation point 2 below
 *
 *        D. Termination signal propagation and comm/receive thread termination
 *          Since we have multiple processes, and multiple threads queuing message for sent, receiving, and processing received messages,
 *          it is important to define how termination signals are propagated through the application.
 *
 *          The CommThread loop terminates when both sending and receiving are complete.
 *            a. sending is completed when finishCommunication has been called and all pending messages are sent,
 *               and all in-progress send requests are completed.
 *            b. receiving is completed when no Application Termination FOC messages (tag == CONTROL_TAG, the last control messages)
 *                are left to receive.  At init, we set the number of Application Termination FOC messages
 *                to the communication size, which are only decremented after finishCommunication is called.
 *
 *          The Receive Thread terminates when receiving is completed and there are no pending received messages to process
 *            a. receiving is completed according to criteria above (b in previous paragraph).
 *                A single Application Termination FOC message would then be queued for processing.
 *            b. Queue for received messages that are pending processing is exhausted.
 *
 *
 *        Rationales
 *          1. Why do we need epoch Ids?
 *              Epoch ids allows the receiver to group all FOC messages with matching tag and epoch from different sources
 *              together.  It also ensure that a receiver gets EXACTLY ONE FOC message from each source during a epoch.
 *              Thus it helps to separate and preserve order of collective calls.
 *
 *              This is particularly important when different MPI processses may be sending FOC messages at different rates.
 *              For example, without epoch id, and using number of FOC messages received == number of senders as flush unblocking
 *              criteria, we would encounter the following problem:  Process A, B, and C generate 0, 1, and 2
 *              flushes for tag X.  B and C, listening for FOC messages with tag X, will receive 3 FOC messages
 *              with tag X, thus unblocked from flush, but they never received messages from A. A, when it finally enters flush, will
 *              have 4 FOC messages with tag X.  If count of FOC message received were reset at this point, subsequent Flush could
 *              be waiting for a FOC message that has already been processed in the last Flush, resulting in deadlock.
 *
 *              An INEFFICIENT alternative is to count the FOC messages per tag from each sender, but then each process
 *              would need O(m*p) memory (m = number of tags, p = number of processes) and O(p^2) time (p FOC messages,
 *              each FOC message with tag X received results in checking an array of size p for possible completion of a flush)
 *
 *              The final design is to use epoch.  For each epoch id and tag id, we count down the number of FOC messages
 *              with that epoch and tag, until the count reaches 0, at which time the flush for that epoch completes.
 *
 *              While this ensures collectiveness of flush call, it does not allow MPI processes to send FOC messages out of order.
 *              If process A flushes tag 1 then tag 2, and process B flushes Tag 2 then Tag 1, because they each locally enter
 *              a blocking wait for the flush operation to complete, the application will deadlock.
 *              example:
 *                proc 1              proc 2
 *                  tag   1 2 1 2         tag   1 2 2 1
 *                  epoch 0 1 2 3         epoch 0 1 2 3
 *              Note that this is NOT different than standard MPI behavior for point-to-point or collective operations.
 *
 *
 *          2. why not use MPI collectives?
 *            1. collective functions called in MPI comm thread needs to occur AFTER
 *                comm thread have cleared send queue, all pending send and receive requests,
 *                and receive threads have processed all received data messages - inter-thread coordination required.
 *            2. a slow process may post send AFTER a fast process already entered collective function call,
 *                resulting in data messages received out of order, likely creating deadlock (as slow process
 *                is waiting for receive before calling collective function) - inter-process coordination required.
 *            The inter-process coordination would require some marker to indicate "last message"
 *              prior synchronization, which is the exact design we have above, without needing MPI collective functions.
 *            In addition, by using Point-To-Point communication, the FOC messages are enqueued like other messages,
 *              so no explicit "wait until queue cleared" is required, thus reducing idle time.
 *
 *        Usage
 *          Control thread instantiates CommunicationLayer with MPI communicator
 *          Control thread addReceiveCallbacks to CommLayer
 *          Control thread calls initCommunication()
 *
 *          in parallel
 *            compute threads call sendMessage(...)
 *
 *          Control thread call flush(tag)
 *
 *          in parallel
 *            compute threads call sendMessage(...)
 *
 *          Control thread call flush(tag)
 *
 *          ...
 *          (internally CommLayer's receive threads call the callbacks.)
 *          ...
 *
 *          Control thread call finish(tag)
 *
 *          control thread call finishCommunication()
 *
 *          Control thread call finalize (TO MAKE SURE ALL THREADS ARE FINISHED)
 *
 * TODO - replace uint8_t by typedef
 *
 */
class CommunicationLayer
{
  protected:
    //==== type definitions


    /// alias MessageBuffersType to Use the thread safe version of the SendMessageBuffers
    using MessageBuffersType = SendMessageBuffers<bliss::concurrent::THREAD_SAFE>;

    /// alias BufferPoolType to MessageBuffersType's
    using BufferPoolType = typename MessageBuffersType::BufferPoolType;

    /// alias BufferIdType to MessageBuffersType's
    using BufferIdType = typename MessageBuffersType::BufferIdType;


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
    struct MessageReceived
    {
      /// The received data
      uint8_t* data;

      /// The number of bytes received
      std::size_t count;

      /// The message tag, indicates the type of message (how to interpret the message)
      int tag;

      /// The message source id
      int src;

      /**
       * @brief constructor using a pre-existing memory block
       *
       * @param data      in memory block of data (received from remote proc, to be processed here)
       * @param count     number of bytes in the data block
       * @param tag       the MPI message tag, indicating the type of message
       * @param src       the MPI source process rank
       */
      MessageReceived(uint8_t* data, std::size_t count, int tag, int src)
        : data(data), count(count), tag(tag), src(src) {}

      /// default constructor
      MessageReceived() = default;

      /**
       * @brief   accessor to get the tag of the data.  returns only data message tag types.
       * @details if tag != CONTROL_TAG, then it's returned (data message type).
       *          if tag == CONTROL_TAG, then this is a control message, the tag is the first integer element in the payload.
       *
       * @return  the (data) message tag associated with this MessageReceived.
       */
      int getTag() const {
        if (tag == CommunicationLayer::CONTROL_TAG)
          return ((int*)data)[0];
        else
          return tag;
      }

      /**
       * @brief   accessor to get the epoch id of a data message flush.
       * @note    See Rationale above for detail about Epoch.
       *
       * @return  the synchronization epoch id.  data messages do not have epoch number so -1 is returned.
       */
      int getEpoch() const {
        if (tag == CommunicationLayer::CONTROL_TAG)
          return ((int*)data)[1];
        else
          return TERMINAL_EPOCH;
      }
    };

    /**
     * @brief Structure to hold all associated information of a message to be sent via MPI.
     * @details MessageToSend is a small metadata stucture to be used for message queuing (to be sent)
     *          MessageToSend does not actually hold any data.  it points to a in-memory Buffer by id.
     *
     *          control message is NOT sent this way, since there is no Buffer associated with a control message.
     *
     * @tparam  BufferIdType  type of Buffer id.  Obtained from SendMessageBuffer, parameterized by Thread Safety.
     */
    template<typename BufferIdType>
    struct MessageToSend
    {
      /// The id of the message buffer
      BufferIdType bufferId;

      /// Tag (type) of the message
      int tag;

      /// Destination rank, i.e., id of the process where the message is to be sent
      int dst;

      /**
       * @brief Constructs a new instance of this struct and sets all members as given.
       *
       * @param _id   Id of MessageBuffer that will be sent
       * @param _tag  type of the message being sent
       * @param _dst  destination of the message
       */
      MessageToSend(BufferIdType _id, int _tag, int _dst)
        : bufferId(_id), tag(_tag), dst(_dst) {}

      /// default constructor
      MessageToSend() = default;
    };


    /// alias SendQueueElementType to MessageToSend, for convenience.
    using SendQueueElementType = MessageToSend<BufferIdType>;


    /// constant indicating message type is control (FOC) message.
    static const int CONTROL_TAG = 0;

    /// constant indicating that there is no active tag. (e.g. during flushing)
    static const int NO_TAG = -1;

    /// terminal epoch id.
    static const int TERMINAL_EPOCH = -1;

public:


  /**
   * @brief Constructor for the CommLayer.
   * @details have to have the comm_size parameter since we can't get it during member variable initialization and recvQueue needs it.
   *
   * @param communicator    The MPI Communicator used in this communicator
   *                        layer.
   * @param comm_size       The size of the MPI Communicator.
   */
  CommunicationLayer (const MPI_Comm& communicator, const int comm_size)
    : sendQueue(3 * comm_size), recvQueue(3 * comm_size), sendAccept(),
      flushing(NO_TAG), finishing(false), sendDone(false), recvDone(false),
      comm(communicator), epoch(0), finalized(false)
  {
    // init communicator rank and size
    MPI_Comm_size(comm, &commSize);
    assert(comm_size == commSize);
    MPI_Comm_rank(comm, &commRank);

    // initialize the mpi buffer for control messaging
    int mpiBufSize = 0;
    // arbitrary, 4X commSize, each has 2 ints, one for tag, one for epoch id,.
    MPI_Pack_size(commSize * 4 * 2, MPI_INT, comm, &mpiBufSize);
    mpiBufSize += commSize * 4 * MPI_BSEND_OVERHEAD;
    char* mpiBuf = (char*)malloc(mpiBufSize);
    MPI_Buffer_attach(mpiBuf, mpiBufSize);
  }

  /**
   * @brief Destructor of the CommmunicationLayer.
   * @details Blocks till all threads are finished.
   * @note    If finish() was not called for all used tags, this will block forever.
   */
  virtual ~CommunicationLayer () {
    finalize();
  }

  /**
   * @brief waits for all threads to finish.  this provides a mechanism to ensure all communications
   *    are complete before subsequent MPI calls are made, such as MPI_FINALIZE.
   */
  void finalize() {
    DEBUGF("FINALIZING COMM LAYER");

    if (finalized) return;

    // TODO: check that we are finished??
	  finishCommunication();

    // wait for both threads to quit
    comm_thread.join();
    callback_thread.join();

    // detach the buffer
    int mpiBufSize;
    char* mpiBuf;
    MPI_Buffer_detach(&mpiBuf, &mpiBufSize);
    free(mpiBuf);

    finalized = true;
  };


  /**
   * @brief Initializes all communication and starts the communication and
   *        callback threads.
   *
   * @note Must be called before any other functions can be called.
   */
  void initCommunication()
  {
    sendDone.store(false);
    recvDone.store(false);
    finishing.store(false);
    flushing.store(NO_TAG);
    recvRemaining[TERMINAL_EPOCH] = commSize;
    comm_thread = std::thread(&CommunicationLayer::commThread, this);
    callback_thread = std::thread(&CommunicationLayer::callbackThread, this);
  }


  /**
   * @brief Sends a message (raw bytes) to the given rank asynchronously.
   * @details
   *   Asynchronously sends the given message to the given rank with the given
   *   message tag. The message is buffered and batched before send. This function
   *   returns immediately after buffering, which does not guarantee that the
   *   message has been sent yet. Only calling `flush()` guarantees that all
   *   messages are sent, received, and processed before the function returns.
   *
   *   internally, if the buffer if full, it will be enqueued for asynchronous MPI send.
   *
   *   This function is THREAD SAFE
   *
   * @param data        A pointer to the data to be sent.
   * @param count       The number of bytes of the message.
   * @param dst_rank    The rank of the destination processor.
   * @param tag         The tag (type) of the message.
   */
  void sendMessage(const void* data, std::size_t count, int dst_rank, int tag)
  {
    //== check to see if the target tag is still accepting
    if (sendAccept.find(tag) == sendAccept.end()) {
      throw std::invalid_argument("invalid tag: tag has been finished already");
    }

    //== don't send control tag.
    assert(tag != CONTROL_TAG && tag >= 0);

    //== if there isn't already a tag listed, add the MessageBuffers for that tag.
    // multiple threads may call this.
    if (buffers.find(tag) == buffers.end()) {
      std::unique_lock<std::mutex> lock(buffers_mutex);
      if (buffers.find(tag) == buffers.end()) {
        // create new message buffer
        buffers[tag] = std::move(MessageBuffersType(commSize, 8192));
      }
    }

    //== try to append the new data - repeat until successful.
    // along the way, if a full buffer's id is returned, queue it for sendQueue.
    BufferIdType fullId = BufferPoolType::ABSENT;
    std::pair<bool, BufferIdType> result;
    do {
      // try append.  append fails if there is no buffer or no room
      result = buffers.at(tag).append(data, count, dst_rank);

      fullId = result.second;

      if (fullId != BufferPoolType::ABSENT) {        // have a full buffer
        // verify that the buffer is not empty (may change because of threading)
        if (!(buffers.at(tag).getBackBuffer(fullId).isEmpty())) {
          // have a non-empty buffer - put in send queue.
          if (!sendQueue.waitAndPush(std::move(SendQueueElementType(fullId, tag, dst_rank)))) {
            throw bliss::io::IOException("ERROR: sendQueue is not accepting new SendQueueElementType due to disablePush");
          }
        }
        // full buffer enqueued, reset it.
        fullId = BufferPoolType::ABSENT;
      } // else don't have a full buffer.

      // repeat until success;
    } while (!result.first);
  }


  /**
   * @brief Registers a callback function for the given message tag.
   * @details
   *      The registered callback function will be called for each received message
   *      of the given type (tag). A callback function has to be registered for each
   *      message type sent. Only one callback function can be registered for each
   *      tag, adding more than one results in a `invalid_argument` exception.
   *
   *      A callback function has the signature:
   *         void(uint8_t* msg, std::size_t count, int fromRank)
   *
   * @note  It is STRONGLY ENCOURAGED that the function be THREAD SAFE, as there may be
   *      more than 1 receive thread calling the functions.
   *
   * @param tag                 The message tag for which the callback function
   *                            is registered. Has to be > 0.
   * @param callbackFunction    The callback-function.
   *
   * @throw std::invalid_argument   If `tag <= 0` or if a callback function has
   *                                already been registered for the given tag.
   */
  void addReceiveCallback(int tag, std::function<void(uint8_t*, std::size_t, int)> callbackFunction) throw (std::invalid_argument)
  {
    // check for valid arguments
    assert(tag != CONTROL_TAG && tag >= 0);

    if (callbackFunctions.find(tag) != callbackFunctions.end()) {
      throw std::invalid_argument("callback function already registered for given tag");
    }

    // add the callback function to a lookup table
    callbackFunctions[tag] = callbackFunction;

    // register the tag so that messages of this type is now accepted by the CommLayer.
    sendAccept.insert(tag);
  }

  /**
   * @brief flushes the message buffers associated with a particular tag.  Blocks until completion
   * @details call by a single thread only. throws exception if another thread is flushing at the same time.
   *        This method sends all buffered messages with a particular tag to all target receivers.
   *        Then it sends a control message to all target receivers and blocks
   *        A process, on receiving control messages from all senders, unblocks the local flush function
   *
   *        The effect is that collective call to flush enforces that all received messages prior
   *        to the last control messages are processed before the flush call is completed.
   *
   *        Each call of this function increments a "epoch" id, which is sent as part of the control
   *        message to all receivers.  Receivers use the epoch id to ensure that one control message
   *        from each sender has been received.
   *
   * @param tag
   * @throw bliss::io::IOException
   */
  void flush(int tag) throw (bliss::io::IOException)
  {
    // cannot flush control tag.
    assert(tag != CONTROL_TAG && tag >= 0);

    if (sendAccept.find(tag) == sendAccept.end()) {
      // already finished, no need for further flushing
      // DEBUGF("NO FLUSHing: already finished tag: %d", tag);
      return;
    }

    DEBUGF("Rank %d FLUSHing tag %d, epoch %d", commRank, tag, epoch);

    // flush all data buffers (put them into send queue)
    flushBuffers(tag);


    // track the tag that is undergoing flush in a thread safe way
    // control thread (calling flush) sets it, while receive thread
    // unsets it in the callback function
    // can only set it if nothing else is flushing.
    int t = NO_TAG;
    if (!flushing.compare_exchange_strong(t, tag)) {
      // TODO throw exception?
      std::stringstream ss;
      ss << "ERROR:  another tag is being flushed: " << t << ". attempt to flush tag " << tag;
      throw bliss::io::IOException(ss.str());
    }

    // send the control message
    sendControlMessages(tag);

    // wait till all control messages have been received (not using MPI_Barrier.  see Rationale.)
    waitForControlMessages(tag);
  }

  /**
   * @brief stops accepting messages of a particular type and flushes existing messages.  Blocks until completion
   * @details call by a single thread only. throws exception if another thread is flushing at the same time.
   *
   *        This function follows the same logic as flush, (see flush documentation),
   *        except that it also sets the CommLayer NOT to accept messages of that type any further.
   *
   * @note  This does NOT send an application termination message
   *
   * @param tag
   * @throw bliss::io::IOException
   */
  void finish(int tag) throw (bliss::io::IOException)
  {
    assert(tag != CONTROL_TAG && tag >= 0);

    if (sendAccept.find(tag) == sendAccept.end()) {
      // already finished.
      //DEBUGF("ERROR:  NO FINISHing: already finished tag: %d", tag);
      return;
    }

    // mark as no more coming in for this tag.  ONE THREAD ONLY modifies this.
    sendAccept.erase(tag);


    DEBUGF("Rank %d FINISHing tag %d, epoch %d", commRank, tag, epoch);

    // flush all buffers (put them into send queue)
    flushBuffers(tag);


    // track the tag that is undergoing flush in a thread safe way
    // control thread (calling flush) sets it, while receive thread
    // unsets it in the callback function
    // can only set it if nothing else is flushing.
    int t = NO_TAG;
    if (!flushing.compare_exchange_strong(t, tag)) {
      // TODO throw exception?
      std::stringstream ss;
      ss << "ERROR:  another tag is being flushed: " << t << ". attempt to flush tag " << tag;
      throw bliss::io::IOException(ss.str());
    }

    // send the control message
    sendControlMessages(tag);


    // wait till all control messages have been received (not using MPI_Barrier.  see Rationale.)
    waitForControlMessages(tag);
  }


  /**
   * @brief stops accepting messages of all types and flushes all existing messages.  Blocks until completion
   * @details call by a single thread only. throws exception if another thread is flushing at the same time.
   *        This method follows the same logic as "finish", but for all message types
   *
   *        At the end, sends an application termination control message.
   *
   * @param tag
   * @throw bliss::io::IOException
   */
  void finishCommunication() throw (bliss::io::IOException)
  {
    DEBUGF("Finish Communication");


	  // already finishing.
	  if (finishing.load() == true) return;


	  std::unordered_set<int> temp(sendAccept);


	  for (auto tag : temp) {

      if (tag != CONTROL_TAG && tag >= 0) {
        finish(tag);
      }
    }

    // track the tag that is undergoing flush in a thread safe way
    // control thread (calling flush) sets it, while receive thread
    // unsets it in the callback function
    // can only set it if nothing else is flushing.
    int t = NO_TAG;
    if (!flushing.compare_exchange_strong(t, CONTROL_TAG)) {
      // TODO throw exception?
      std::stringstream ss;
      ss << "ERROR:  another tag is being flushed: " << t << ". attempt to flush tag " << CONTROL_TAG;
      throw bliss::io::IOException(ss.str());
    }

    // send the control message
    sendControlMessages(CONTROL_TAG);

    // mark as finishing.
    finishing.store(true);

    // wait till all control messages have been received (not using MPI_Barrier.  see Rationale.)
    waitForControlMessages(CONTROL_TAG);

  }



  /**
   * @brief Returns the size of the MPI communicator.
   *
   * @return The number of processes that are part of the MPI communicator.
   */
  int getCommSize() const
  {
    return commSize;
  }

  /**
   * @brief Returns the rank of this process in the MPI communicator.
   *
   * @return The rank of this process.
   */
  int getCommRank() const
  {
    return commRank;
  }

protected:


  /**
   * @brief The function for the dedicated communicator thread.
   * @details  loops until all sending and all receiving are done.
   *    these are atomic variables since we have multiple threads.
   *
   *    completes pending send and recv requests, then receive any new and send any queued messages.
   *
   * @note  There should be only 1 thread calling this function.
   */
  void commThread()
  {
    DEBUGF("THREAD STARTED:  comm-thread on %d", commRank);

    // while there's still work to be done:
    while (!sendDone.load() || !recvDone.load())
    {
      // first clean up finished operations
      finishSends();
      finishReceives();

      // start pending receives
      tryStartReceive();
      // start pending sends
      tryStartSend();
    }
    DEBUGF("THREAD FINISHED:  comm-thread on %d", commRank);
  }

  /**
   * @brief The function for the dedicated callback-handler thread.
   * @details   loops until all received messages in queue are processed, and no more possibility of receiving messages.
   *
   *        Control messages are added in the queue if receiver has received that message from all senders, so it marks the
   *        completion of flush or finish or finishCommunication, which are unblocked.
   *
   *        Data messages are processed directly via callbacks according to the message type (in tag).
   *
   * @note  can have multiple instances of this.
   */
  void callbackThread()
  {
    DEBUGF("THREAD STARTED:  recv-thread on %d", commRank);
    MessageReceived msg;

    // while there are messages to process, or we are still receiving
    while (!recvDone.load() || !recvQueue.isEmpty() )
    {
      // get next element from the queue, wait if none is available.
      // waitAndPop will exit out of wait when termination flag is set on the recvQueue
      auto result = std::move(recvQueue.waitAndPop());
      if (!result.first) {
        // no valid result, we must be done with receiving
        assert(recvDone.load());
        break;
      }

      // get the message
      msg = result.second;

      // if the message is a control message (enqueued ONE into RecvQueue only after all control messages
      //  for a tag have been received.)
      if (msg.tag == CONTROL_TAG) {
        // get control message's target tag.
        int tag = msg.getTag();

        // and unblock the flush/finish/finishCommunication call.
        // getting to here means that all control messages for a tag have been received from all senders.
        if (flushing.compare_exchange_strong(tag, NO_TAG)) {
          // so unset 'flushing', and unblock flush() via condition variable.
          flushBarrier.notify_all();
        } else {
          ERRORF("control message with tag=%d, but flushing contained %d", msg.tag, tag);
        }
      } else {
        // data message.  process it.
        (callbackFunctions[msg.tag])(msg.data, msg.count, msg.src);
        // clean up message data
        delete [] msg.data;
      }
    }
    DEBUGF("THREAD FINISHED: recv-thread on %d", commRank);
  }


  /**
   * @brief Sends FOC message for the given tag to every MPI process
   *
   * @param tag     The MPI tag to end.
   * @throw bliss::io::IOException  If the sendQueue has been disabled.
   */
  void sendControlMessages(int tag) throw (bliss::io::IOException)
  {
    for (int i = 0; i < getCommSize(); ++i) {
      // send end tags in circular fashion, send tag to self last
      int target_rank = (i + getCommRank() + 1) % getCommSize();

      //DEBUGF("Rank %d sendEndTags %d epoch %d target_rank = %d ", commRank, tag, epoch, target_rank );

      // send the end message for this tag.
      if (!sendQueue.waitAndPush(std::move(SendQueueElementType(BufferPoolType::ABSENT, tag, target_rank)))) {
        throw bliss::io::IOException("ERROR: sendQueue is not accepting new SendQueueElementType due to disablePush");
      }
    }
  }


  /**
   * @brief Blocks till the flusing of the given tag is completed (i.e., all
   *        FOC messages have been received).
   *
   * @param tag The tag to wait for.
   */
  void waitForControlMessages(int tag) throw (bliss::io::IOException)
  {
    // waiting for all messages of this tag to flush.
    if (flushing.load() == tag) {  // precheck before locking.
      std::unique_lock<std::mutex> lock(flushMutex);
      do {
        // condition variable waits for notification, which is generated by recv thread
        // when all FOC messages with this tag are received.
        flushBarrier.wait(lock);
      } while (flushing.load() == tag);

    } else {
      throw bliss::io::IOException("ERROR: waitForEndTags called but currently flushing a different tag.");
    }

    DEBUGF("Rank %d received all END message for TAG %d, epoch = %d", commRank, tag, epoch);
    ++epoch;
  }


  /**
   * @brief Flushes all data buffers for message tag `tag`.
   * @details All pending/buffered messages of type `tag` will be send out.
   *    enqueues data buffers to be sent, but enqueuing the buffer's id, if buffer is not empty.
   *
   * @param tag The message tag, whose messages are flushed.
   * @throw bliss::io::IOException
   */
  void flushBuffers(int tag) throw (bliss::io::IOException)
  {
    // no MessageBuffers with the associated tag.  end.
    if (buffers.find(tag) == buffers.end())
    {
      DEBUGF("NO BUFFERS for tag: %d", tag);
      return;
    }
    //DEBUGF("Active Buffer Ids in message buffer: %s", buffers.at(tag).bufferIdsToString().c_str());

    // flush out all the send buffers matching a particular tag.
    int idCount = buffers.at(tag).getBufferIds().size();
    assert(idCount == commSize);

    for (int i = 0; i < idCount; ++i) {
      // flush buffers in a circular fashion, starting with the next neighbor
      int target_rank = (i + getCommRank() + 1) % getCommSize();

      auto id = buffers.at(tag).getBufferId(target_rank);
      // flush/send all remaining non-empty buffers
      if ((id != BufferPoolType::ABSENT) && !(buffers.at(tag).getBackBuffer(id).isEmpty())) {
        if (!sendQueue.waitAndPush(std::move(SendQueueElementType(id, tag, target_rank)))) {
          throw bliss::io::IOException("ERROR: sendQueue is not accepting new SendQueueElementType due to disablePush");
        }
      }
    }
  }





  /**
   * @brief Takes queued messages and initiates MPI_ISends for them.
   *
   * @throw bliss::io::IOException
   */
  void tryStartSend() throw (bliss::io::IOException)
  {
    // quit if no more sends are to be expected
    if (sendDone.load()) return;

    // try to get the send element
    SendQueueElementType se;
    auto el = std::move(sendQueue.tryPop());
    if (!el.first) {  // no valid entry in the queue
      return;
    }
    // else there is a valid entry from sendQueue

    se = el.second;
    if (se.bufferId == BufferPoolType::ABSENT) {
      // Is this a control message?  if control message, SendQueueElement will have the
      // "ABSENT" buffer id, while tag is the message type to control

      //DEBUGF("SEND %d -> %d, termination signal for tag %d", commRank, se.dst, se.tag);

      // termination message for this tag and destination
      if (se.dst == commRank) {
        // local, directly handle by creating an byte array and directly
        // insert into the recvInProgress (needed by the receiver decrement
        // logic in "finishReceive", so can't put into receive queue)
        MPI_Request req = MPI_REQUEST_NULL;
        uint8_t *array = new uint8_t[2 * sizeof(int)];
        memcpy(array, &(se.tag), sizeof(int));
        memcpy(array + sizeof(int), &(epoch), sizeof(int));
        recvInProgress.push(std::move(std::pair<MPI_Request,
              MessageReceived>(std::move(req),
                std::move(MessageReceived(array, 2 * sizeof(int), CONTROL_TAG, commRank)))));
      } else {
        // remote.  send a terminating message with tag CONTROL_TAG
        // payload is the message type and epoch.  Message Ordering is guaranteed by MPI.
        // Bsend to allow use of preallocated buffer.
        int stat[2] = {se.tag, epoch};
        MPI_Bsend(stat, 2, MPI_INT, se.dst, CONTROL_TAG, comm);
      }
    } else {
      // this is an actual data message:

      // get message data and it's size
      void* data = const_cast<void*>(buffers.at(se.tag).getBackBuffer(se.bufferId).getData());
      auto count = buffers.at(se.tag).getBackBuffer(se.bufferId).getSize();

      if (count == 0) {
        WARNINGF("WARNING: NOT SEND %d -> %d: 0 byte message for tag %d, bufferid = %d.", commRank, se.dst, se.tag, se.bufferId);
      } else {

        if (se.dst == commRank) {
          // local, directly handle by creating an output object and directly
          // insert into the recv queue (thread safe)
          uint8_t* array = new uint8_t[count];
          memcpy(array, data, count);

          // out of order receive, but can only be early, so okay.
          if (!recvQueue.waitAndPush(std::move(MessageReceived(array, count, se.tag, commRank)))) {
            throw bliss::io::IOException("ERROR: recvQueue is not accepting new receivedMessage due to disablePush");
          }
          // finished inserting directly to local RecvQueue.  release the buffer
          buffers.at(se.tag).releaseBuffer(se.bufferId);
        } else {
          // remote: initiate async MPI message
          MPI_Request req;
          MPI_Isend(data, count, MPI_UINT8_T, se.dst, se.tag, comm, &req);
          sendInProgress.push(std::move(std::pair<MPI_Request, SendQueueElementType>(std::move(req), std::move(se))));
        }
      }
    }

  }

  /**
   * @brief Finishes pending MPI_Send requests.  stops at first unfinished request, left unfinished requests finish in next iteration's call.
   */
  void finishSends()
  {
    // no more messages to send or being sent, so done.
    if (sendDone.load()) return;

    int finished = 0;

    // whiel there is some in progress MPI_send requests,
    while(!sendInProgress.empty())
    {
      // get the first request to check
      std::pair<MPI_Request, SendQueueElementType>& front = sendInProgress.front();

      // check if MPI request has completed.
      if (front.first == MPI_REQUEST_NULL) {
        // this is a control message since it has no "request".
        finished = 1;
      } else {
        // has request.  check for finished
        MPI_Test(&front.first, &finished, MPI_STATUS_IGNORE);
      }

      if (finished)
      {
        // if there is a buffer, then release it.
        if (front.second.bufferId != BufferPoolType::ABSENT) {
          // cleanup, i.e., release the buffer back into the pool
          buffers.at(front.second.tag).releaseBuffer(front.second.bufferId);
        }
        // remove the entry from the in-progress queue.
        sendInProgress.pop();
      }
      else
      {
        // the head of the queue is not finished, so break and get it later.
        break;
      }
    }

    // all sending is done when "finishing" is set, and no more messages are pending send,
    // and no messages are being sent.
    if (finishing.load() && sendQueue.isEmpty() && sendInProgress.empty()) {
      // not using sendAccept as check - sendAccept is cleared BEFORE the last
      // tag end messages starts to be sent, and definitely before the app end
      // message starts to be sent so sendDone may be set to true too early,
      // and elements for the last tag, or tag 0 (termination) may be added to
      // the sendQueue or sendInProfess queue after the empty check here.
      sendDone.store(true);
      //DEBUGF("send Done!");
      //sendQueue.disablePush();
    }
  }


  /**
   * @brief Initiates async MPI receives.
   */
  void tryStartReceive()
  {
    // if no more are expected, then done
    if (recvDone.load()) return;

    // probe for messages
    int hasMessage = 0;
    MPI_Status status;
    MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, comm, &hasMessage, &status);

    // if have message to receive,
    if (hasMessage > 0) {
      // get some message details
      int src = status.MPI_SOURCE;
      int tag = status.MPI_TAG;
      int received_count;
      MPI_Request req;
      uint8_t* msg_data;


      // receive the message data bytes into a vector of bytes
      if (tag == CONTROL_TAG) {  // control messages
        MPI_Get_count(&status, MPI_INT, &received_count);

        // receive data or FOC messages.  use finishReceive to handle the
        // FOC messages and count decrement

        msg_data = (received_count == 0) ? nullptr : new uint8_t[received_count * sizeof(int)];
        MPI_Irecv(msg_data, received_count, MPI_INT, src, tag, comm, &req);

      } else {  // data messages
        MPI_Get_count(&status, MPI_UNSIGNED_CHAR, &received_count);

        msg_data = (received_count == 0) ? nullptr : new uint8_t[received_count];
        MPI_Irecv(msg_data, received_count, MPI_UNSIGNED_CHAR, src, tag, comm, &req);
      }

      // insert into the received InProgress queue.
      recvInProgress.push(std::pair<MPI_Request, MessageReceived>(req,
            MessageReceived(msg_data, received_count, tag, src)));
    }
  }

  /**
   * @brief Finishes pending receive requests.
   * @details   if control message, count down from the total number of possible senders for an epoch
   *            when reaching 0, then enqueue 1 control message for callback handling to process
   *
   *            if data message, then enqueue for callback handling to process.
   *
   *
   * @throw bliss::io::IOException
   */
  void finishReceives() throw (bliss::io::IOException)
  {
    // no more to receive or pending receive requests, so don.
    if (recvDone.load()) return;

    int finished = 0;


    // if there are pending receive requests
    while(!recvInProgress.empty())
    {
      // get the request
      std::pair<MPI_Request, MessageReceived>& front = recvInProgress.front();

      if (front.first == MPI_REQUEST_NULL) {  // if local messages
        finished = 1;  // then this request is complete
      } else {
        // test if the MPI request for receive has been completed.
        MPI_Test(&front.first, &finished, MPI_STATUS_IGNORE);
      }

      if (finished)
      {
        // finished request, so enqueue the message for processing.


        // if it's a control message, then need to count down
        if (front.second.tag == CONTROL_TAG)  {

          // get the message type being controlled, and the epoch
          int tag = front.second.getTag();
          int epoch = front.second.getEpoch();

          //== if the message type is same as CONTROL_TAG, then this is a Application Termination FOC message.
          if (tag == CONTROL_TAG) {
            epoch = TERMINAL_EPOCH;  // then use special epoch TERMINAL_EPOCH, which has counter set at commlayer init.
          } else {
            // else we are controlling a data message type

            // if we haven't seen this epoch, add a new entry.
            if (recvRemaining.find(epoch) == recvRemaining.end())
              recvRemaining[epoch] = commSize;  // insert new.
          }

          //== now decrement the count of control messages for a epoch
          --recvRemaining.at(epoch);
          DEBUGF("RECV rank %d receiving END signal tag %d epoch %d from %d, num senders remaining is %d",
                 commRank, tag, epoch, front.second.src, recvRemaining.at(epoch));

          //== if after decrement the count is 0 for the epoch,
          if (recvRemaining.at(epoch) == 0) {
            // then we've received all end messages for this epoch.  remove it from tracking.
            recvRemaining.erase(epoch);

            //DEBUGF("ALL END received for tag %d, pushing to recv queue", front.second.tag);

            // and enqueue ONE FOC message for this tag to be handled by callback.  (reduction from commSize to 1)
            if (!recvQueue.waitAndPush(std::move(front.second))) {
              throw bliss::io::IOException("ERROR: recvQueue is not accepting new receivedMessage due to disablePush");
            }

          } else if (recvRemaining.at(epoch) < 0) {
            // count decremented to negative.  should not happen.  throw exception.
            std::stringstream ss;
            ss << "ERROR: number of remaining receivers for tag " << tag << " epoch " << epoch << " is now NEGATIVE";
            throw bliss::io::IOException(ss.str());
          }
        } else {
          // Data message.  add the received messages into the recvQueue
          if (!recvQueue.waitAndPush(std::move(front.second))) {
            throw bliss::io::IOException("ERROR: recvQueue is not accepting new receivedMessage due to disablePush");
          }
        }
        // remove pending receive element
        recvInProgress.pop();
      }
      else
      {
        // the recv request is not done, so stop and wait for the next cycle.
        break;
      }
    }

    // see if we are done with all pending receives and all Application Termination FOC messages are received.
    if ((recvRemaining.size() == 0) && recvInProgress.empty()) {
      // if completely done, epoch=TERMINAL_EPOCH would have been erased from recvRemaining.
      recvDone.store(true);
      //DEBUGF("recv Done!");
      recvQueue.disablePush();
    }
  }


protected:
/******************************************************************************
 *                           Protected member fields                            *
 ******************************************************************************/

  /* Queues of pending MPI operations (MPI_Requests)
   * Used ONLY by the comm-thread, purely internal and purely single thread.
   */
  /// Queue of pending MPI receive operations
  std::queue<std::pair<MPI_Request, MessageReceived> > recvInProgress;
  /// Queue of pending MPI send operations
  std::queue<std::pair<MPI_Request, SendQueueElementType> > sendInProgress;


  /* Thread-safe queues for posting sends and getting received messages,
   * these are used as communication medium between the comm-thread and the
   * (potetially multiple) producer threads for the sends; and the comm-thread
   * and the callback handler thread for the receives
   */

  // Outbound message data structure, Multiple-Producer-Single-Consumer queue
  // consumed by the internal MPI-comm-thread
  /// message queue between sendMessage calling threads (src) and comm-thread
  /// (sink)
  bliss::concurrent::ThreadSafeQueue<SendQueueElementType> sendQueue;

  // Inbound message data structure, Single-Producer-Multiple-Consumer queue
  // produced by the internal MPI-comm-thread
  /// message queue between comm thread (src) and callback thread (sink)
  bliss::concurrent::ThreadSafeQueue<MessageReceived> recvQueue;


  /* information per message tag */

  /// Message buffers per message tag (maps each tag to a set of buffers, one buffer for each target rank)
  std::unordered_map<int, MessageBuffersType> buffers;
  /// Mutex for adding new tags to the buffers map (since inserts can be multi-threaded)
  mutable std::mutex buffers_mutex;

  // sending threads read, comm-thread modifies
  /// set of active message tags (still accepting new messages)
  std::unordered_set<int> sendAccept;

  // compute thread set during initialization and setting callbacks, callback-thread read and modify.
  /// Per tag: number of processes that haven't sent the END-TAG message yet
  std::unordered_map<int, int> recvRemaining;


  /// condition variable to make "flush" and "finish" blocking.
  std::condition_variable flushBarrier;
  /// Associated mutex for the flushBarrier condition variable
  mutable std::mutex flushMutex;

  /// Atomic status variable, -1: no flushing, otherwise holds the value of the
  /// tag which is currently flushed
  std::atomic<int> flushing;

  /// Atomic status variable: whether the commlayer is finishing (after all
  /// tags have been finished).
  std::atomic<bool> finishing;

  /// Atomic status variable: set to true if `finishing` and if all sends are
  /// completed.
  std::atomic<bool> sendDone;
  /// Atomic status variable: set to true if `finishing` and if all receives
  /// are completed.
  std::atomic<bool> recvDone;

  /// std::thread object for the dedicated communication thread
  std::thread comm_thread;
  /// std::thread object for the dedicated callback-handler thread
  std::thread callback_thread;
  /// The MPI Communicator object for this communication layer
  MPI_Comm comm;

  /// Registry of callback functions, mapped to by the associated tags
  std::unordered_map<int,std::function<void(uint8_t*, std::size_t, int)> > callbackFunctions;

  /// The MPI Communicator size
  int commSize;

  /// The MPI Communicator rank
  int commRank;

  /// id of the epochs.  aka the periods between barrier/synchronization.  each epoch is associated with a single flush/finish, thus to 1 tag.
  int epoch;

  bool finalized;
};

/// CONTROL_TAG definition.  (declaration and initialization inside class).
const int CommunicationLayer::CONTROL_TAG;

/// NO_TAG definition.  (declaration and initialization inside class).
const int CommunicationLayer::NO_TAG;

/// TERMINAL_EPOCH definition.  (declaration and initialization inside class).
const int CommunicationLayer::TERMINAL_EPOCH;



} // namespace io
} // namespace bliss

#endif // BLISS_COMMUNICATION_LAYER_HPP
