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
#include <functional>
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

    // ========= static variables for control messages.
   union TaggedEpoch {
       struct {  // order:  lower bits = epoch, higher bits = tag.
           int32_t epoch;
           int32_t tag;
       };
       int64_t tagged_epoch;
   };

   struct TaggedEpochHash {
       inline size_t operator()(const TaggedEpoch& v) const
       {
         return std::hash<int64_t>()(v.tagged_epoch);
       }
   };
   struct TaggedEpochEqual {
       inline bool operator()(const TaggedEpoch& lhs, const TaggedEpoch& rhs) const
       {
         return std::equal_to<int64_t>()(lhs.tagged_epoch, rhs.tagged_epoch);
       }
   };

   /// constant indicating message type is control (FOC) message.
   static constexpr TaggedEpoch CONTROL_TAG{ 0  };

   /// constant indicating that there is no active tag. (e.g. during flushing)
   static const int NO_TAG = -1;

   /// terminal epoch id.
   static const int NO_EPOCH = -1;



   class TagMetadata {
     protected:
       std::atomic<bool> finished;
       std::atomic<int> nextEpoch;
       std::unique_ptr<std::mutex> mutex;
       std::unique_ptr<std::condition_variable> condVar;

     public:

       TagMetadata() : finished(false), nextEpoch(0), mutex(new std::mutex), condVar(new std::condition_variable) {};

       TagMetadata(TagMetadata&& other) {
         finished.exchange(other.finished.load());
         nextEpoch.exchange(other.nextEpoch.load());

       };
       TagMetadata& operator=(TagMetadata&& other) {
         finished.exchange(other.finished.load());
         nextEpoch.exchange(other.nextEpoch.load());
         return *this;
       };

       TagMetadata(const TagMetadata& other) = delete;
       TagMetadata& operator=(const TagMetadata& other) = delete;

       void finish() {
         finished.store(true, std::memory_order_release);
       }
       bool isFinished() {
         return finished.load(std::memory_order_consume);
       }

       int getNextEpoch() {
         return nextEpoch.fetch_add(1);
       }

       std::unique_lock<std::mutex> getUniqueLock() {
         return std::unique_lock<std::mutex>(*mutex);
       }
       // TODO: other lock types

       void wait(std::unique_lock<std::mutex>& lock) {
         condVar->wait(lock);
       }

       void notifyAll() {
         condVar->notify_all();
       }
       void notifyOne() {
         condVar->notify_one();
       }

   };


    //==== type definitions


    /// alias MessageBuffersType to Use the thread safe version of the SendMessageBuffers
    using MessageBuffersType = SendMessageBuffers<bliss::concurrent::THREAD_SAFE>;

    /// alias BufferPoolType to MessageBuffersType's
    using BufferPoolType = typename MessageBuffersType::BufferPoolType;

    /// alias BufferIdType to MessageBuffersType's
    using BufferIdType = typename MessageBuffersType::BufferIdType;


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
      /// The received data
      TaggedEpoch control;

      /**
       * @brief constructor using a pre-existing memory block
       *
       * @param data      in memory block of data (received from remote proc, to be processed here)
       * @param count     number of bytes in the data block
       * @param tag       the MPI message tag, indicating the type of message
       * @param src       the MPI source process rank
       */
      ControlMessage(TaggedEpoch _control, int rank)
        : MPIMessage(CONTROL_TAG.tag, rank), control(_control) {}

      /**
       * @brief Constructs a new instance of this struct and sets all members as given.
       *
       * @param _id   Id of MessageBuffer that will be sent
       * @param _tag  type of the message being sent
       * @param _epoch  epoch in which to send the message (epoch of the tag)
       * @param _dst  destination of the message
       */
      ControlMessage(int _tag, int _epoch, int rank)
        : MPIMessage(CONTROL_TAG.tag, rank), control{_tag, _epoch} {}


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
     *          MessageToSend does not actually hold any data.  it points to a in-memory Buffer by id.
     *
     *          control message is NOT sent this way, since there is no Buffer associated with a control message.
     *
     * @tparam  BufferIdType  type of Buffer id.  Obtained from SendMessageBuffer, parameterized by Thread Safety.
     */
    template<typename BufferIdType>
    struct DataMessageToSend : public MPIMessage
    {
      /// The id of the message buffer
      BufferIdType bufferId;

      /**
       * @brief Constructs a new instance of this struct and sets all members as given.
       *
       * @param _id   Id of MessageBuffer that will be sent
       * @param _tag  type of the message being sent
       * @param _dst  destination of the message
       */
      DataMessageToSend(BufferIdType _id, int _tag, int _dst)
        : MPIMessage(_tag, _dst), bufferId(_id) {}


      /// default constructor
      DataMessageToSend() = default;

      virtual ~DataMessageToSend() {};
    };


    //==== Message type aliases for convenience.

    using SendDataElementType = DataMessageToSend<BufferIdType>;



    // ==== send comm thread class.  for organization.
    // TODO: group the code for comm thread
    class SendCommThread
    {
      public:
        SendCommThread(CommunicationLayer& _comm_layer, const MPI_Comm& communicator) :
          commLayer(_comm_layer),
          comm(communicator),
          sendDone(false),
          sendInProgress(),
          td()
        {
        };

        SendCommThread() = delete;

        virtual ~SendCommThread() {
        }

        void start() {
          int comm_size;
          MPI_Comm_size(comm, &comm_size);

          if (!td.joinable()) {

            //===== initialize the mpi buffer for control messaging, for MPI_Bsend
            int mpiBufSize = 0;
            // arbitrary, 4X commSize, each has 2 ints, one for tag, one for epoch id,.
            MPI_Pack_size(comm_size * 4 * 2, MPI_INT, comm, &mpiBufSize);
            mpiBufSize += comm_size * 4 * MPI_BSEND_OVERHEAD;
            char* mpiBuf = (char*)malloc(mpiBufSize);
            MPI_Buffer_attach(mpiBuf, mpiBufSize);

            td = std::move(std::thread(&bliss::io::CommunicationLayer::SendCommThread::run, this));
          } // else already started the thread.
        }

        void finish() {
          if (td.joinable()) {
            td.join();

            DEBUGF("C sendCommThread finished, thread joined.  about to detach mpi buffer.");

            // detach the buffer for Bsend.
            int mpiBufSize;
            char* mpiBuf;
            MPI_Buffer_detach(&mpiBuf, &mpiBufSize);
            free(mpiBuf);
          } // else already joined.

        }

      protected:
        // parent, in which the thread runs.
        CommunicationLayer& commLayer;

        /// The MPI Communicator object for this communication layer
        MPI_Comm comm;

        /// Atomic status variable: set to true if `finishing` and if all sends are
        /// completed.
        std::atomic<bool> sendDone;

        /// Queue of pending MPI send operations
        std::deque<std::pair<MPI_Request, std::unique_ptr<MPIMessage> > > sendInProgress;

        std::thread td;


        /**
         * @brief The function for the dedicated communicator thread.
         * @details  loops until all sending and all receiving are done.
         *    these are atomic variables since we have multiple threads.
         *
         *    completes pending send and recv requests, then receive any new and send any queued messages.
         *
         * @note  There should be only 1 thread calling this function.
         */
        void run()
        {
          DEBUGF("C THREAD STARTED:  comm-thread on %d", commLayer.commRank);

          // while there's still work to be done:
          while (!sendDone.load())
          {
            // first clean up finished operations
            finishSends();
            // start pending sends
            tryStartSend();
          }
          DEBUGF("C THREAD FINISHED:  comm-thread on %d", commLayer.commRank);
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

          //DEBUGF("rank %d sendDone? %s, finishing %s sendqueue %ld sendInProgress %ld", commRank, (sendDone.load() ? "true":"false"), (finishing.load() ? "true": "false"), sendQueue.getSize(), sendInProgress.size());



          // try to get the send element
          auto el = std::move(commLayer.sendQueue.tryPop());
          if (!el.first) {  // no valid entry in the queue
            return;
          }
          // else there is a valid entry from sendQueue
          assert(el.second.get() != nullptr);

          // second is a pointer to ControlMessage
          if (el.second.get()->tag == CONTROL_TAG.tag) {
            // if control message,

            // ControlMessage ownership is maintained by el.
            ControlMessage* se = dynamic_cast<ControlMessage*>(el.second.get());
            TaggedEpoch te = se->control;

            DEBUGF("C SEND %d -> %d, termination signal for tag %d epoch %d", commLayer.commRank, se->rank, te.tag, te.epoch);

            // termination message for this tag and destination
            if (se->rank == commLayer.commRank) {
              // local, directly handle by creating an byte array and directly
              // insert into the recvInProgress (needed by the receiver decrement
              // logic in "finishReceive", so can't put into receive queue)
              MPI_Request req = MPI_REQUEST_NULL;
//              uint8_t *array = new uint8_t[sizeof(TaggedEpoch)];
//              memcpy(array, &(te), sizeof(TaggedEpoch));


              // now move ControlMessage ownership to recvInProgress
              auto pair = std::make_pair(std::move(req),
                                         std::move(el.second));
              assert(el.second.get() == nullptr);
              assert(pair.second.get() != nullptr);

//              commLayer.recvInProgress.push_back(std::move(pair));
//              assert(pair.second.get() == nullptr);
//              assert(commLayer.recvInProgress.back().second.get() != nullptr);

              if (!commLayer.recvInProgress.waitAndPush(std::move(pair)) ) {
                throw bliss::io::IOException("recvInProgress queue is disabled!");
              }

//              DEBUGF("C rank %d tag %d epoch %d local send control message to rank %d", commLayer.commRank, te.tag, te.epoch, commLayer.commRank);
            } else {
//               remote.  send a terminating message with tag CONTROL_TAG
//               payload is the message type and epoch.  Message Ordering is guaranteed by MPI.
//               Bsend to allow use of preallocated buffer.
              MPI_Bsend(&(te), 2, MPI_INT, se->rank, CONTROL_TAG.tag, comm);

//              std::unique_ptr<ControlMessage> dmr(new ControlMessage(TaggedEpoch(), se->rank));
//
//              MPI_Request req;
//              MPI_Isend(&(te), 2, MPI_INT, se->rank, CONTROL_TAG.tag, comm, &req);
//              recvInProgress.push_back(std::move(std::make_pair(std::move(req),
//                                 std::move(dmr))));
              DEBUGF("C rank %d tag %d epoch %d remote send control message to rank %d", commLayer.commRank, te.tag, te.epoch, se->rank);
            }
          } else {
            // this is an actual data message.


            SendDataElementType* se = dynamic_cast<SendDataElementType*>(el.second.get());

            // get message data and it's size
            void* data = const_cast<void*>(commLayer.buffers.at(se->tag).getBackBuffer(se->bufferId).getData());
            auto count = commLayer.buffers.at(se->tag).getBackBuffer(se->bufferId).getSize();

            if (count > 0) {

              if (se->rank == commLayer.commRank) {
                // local, directly handle by creating an output object and directly
                // insert into the recv queue (thread safe)
                std::unique_ptr<uint8_t[]> array(new uint8_t[count]);
                memcpy(array.get(), data, count);

                std::unique_ptr<DataMessageReceived> dmr(new DataMessageReceived(std::move(array), count, se->tag, commLayer.commRank));


                // put into recvInProgress queue instead of recvQueue to minimize out of order (early) receive.
                MPI_Request req = MPI_REQUEST_NULL;

                auto pair = std::make_pair(std::move(req), std::move(dmr));
                assert(dmr.get() == nullptr);
                assert(pair.second.get() != nullptr);

//                commLayer.recvInProgress.push_back(std::move(pair));
//                assert(pair.second.get() == nullptr);
//
//                assert(commLayer.recvInProgress.back().second.get() != nullptr);

                if (!commLayer.recvInProgress.waitAndPush(std::move(pair)) ) {
                  throw bliss::io::IOException("recvInProgress queue is disabled!");
                }

//                // out of order receive, but can only be early, so okay.
//                if (!commLayer.recvQueue.waitAndPush(std::move(dmr))) {
//                  throw bliss::io::IOException("C ERROR: recvQueue is not accepting new receivedMessage due to disablePush");
//                }
                // finished inserting directly to local RecvQueue.  release the buffer
                commLayer.buffers.at(se->tag).releaseBuffer(se->bufferId);
              } else {
                // remote: initiate async MPI message
                MPI_Request req;
                MPI_Isend(data, count, MPI_UINT8_T, se->rank, se->tag, comm, &req);

                auto pair = std::make_pair(std::move(req), std::move(el.second));
                assert(el.second.get() == nullptr);
                assert(pair.second.get() != nullptr);

                sendInProgress.push_back(std::move(pair));
                assert(pair.second.get() == nullptr);

                assert(sendInProgress.back().second.get() != nullptr);

              }
//            } else {
//              WARNINGF("C WARNING: NOT SEND %d -> %d: 0 byte message for tag %d, bufferid = %d.", commLayer.commRank, se->rank, se->tag, se->bufferId);
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

          // while there is some in progress MPI_send requests,
          while(!sendInProgress.empty())
          {
            // get the first request to check
            assert(sendInProgress.front().second.get() != nullptr);


            // check if MPI request has completed.
            if (sendInProgress.front().first == MPI_REQUEST_NULL) {
              // this is a control message since it has no "request".
              finished = 1;
            } else {
              // has request.  check for finished
              MPI_Test(&(sendInProgress.front().first), &finished, MPI_STATUS_IGNORE);
            }

            if (finished)
            {
              // if there is a buffer, then release it.
              if (sendInProgress.front().second.get()->tag != CONTROL_TAG.tag) {
                SendDataElementType* msg = dynamic_cast<SendDataElementType*>(sendInProgress.front().second.get());
                if (msg->bufferId != BufferPoolType::ABSENT) {
                  // cleanup, i.e., release the buffer back into the pool
                  commLayer.buffers.at(msg->tag).releaseBuffer(msg->bufferId);
                }
              }
              // remove the entry from the in-progress queue.  object destruction will clean up data associated with tthe pointer.
              sendInProgress.pop_front();
            }
            else
            {
              // the head of the queue is not finished, so break and get it later.
              //DEBUGF("not done sending!");
              break;
            }
          }

          // all sending is done when "finishing" is set, and no more messages are pending send,
          // and no messages are being sent.
          //if (finishing.load() && sendQueue.isEmpty() && sendInProgress.empty()) {
          if (!commLayer.sendQueue.canPop() && sendInProgress.empty()) {
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

    };




    // ==== comm Recv thread class.  for organization.
    // TODO: group the code for comm thread
    class RecvCommThread
    {
      public:
        RecvCommThread(CommunicationLayer& _comm_layer, const MPI_Comm& communicator) :
          commLayer(_comm_layer),
          comm(communicator),
          recvDone(false),
          td()
        {
        };

        RecvCommThread() = delete;

        virtual ~RecvCommThread() {
        }

        void start() {

          if (!td.joinable()) {
            td = std::move(std::thread(&bliss::io::CommunicationLayer::RecvCommThread::run, this));
          } // else already started the thread.
        }

        void finish() {
          if (td.joinable()) {
            td.join();

            DEBUGF("C recvCommThread finished, thread joined.");

          } // else already joined.

        }

      protected:
        // parent, in which the thread runs.
        CommunicationLayer& commLayer;

        /// The MPI Communicator object for this communication layer
        MPI_Comm comm;

        /// Atomic status variable: set to true if `finishing` and if all receives
        /// are completed.
        std::atomic<bool> recvDone;


        std::thread td;


        /**
         * @brief The function for the dedicated communicator thread.
         * @details  loops until all sending and all receiving are done.
         *    these are atomic variables since we have multiple threads.
         *
         *    completes pending send and recv requests, then receive any new and send any queued messages.
         *
         * @note  There should be only 1 thread calling this function.
         */
        void run()
        {
          DEBUGF("C THREAD STARTED:  comm-thread on %d", commLayer.commRank);

          // while there's still work to be done:
          while (!recvDone.load())
          {
            // first clean up finished operations
            finishReceives();

            // start pending receives
            tryStartReceive();
          }
          DEBUGF("C THREAD FINISHED:  comm-thread on %d", commLayer.commRank);
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

          if (commLayer.commSize > 1)
            MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, comm, &hasMessage, &status);

          // if have message to receive,
          if (hasMessage > 0) {
            // get some message details
            int src = status.MPI_SOURCE;
            int tag = status.MPI_TAG;
            int received_count;
            MPI_Request req;


            // receive the message data bytes into a vector of bytes
            if (tag == CONTROL_TAG.tag) {  // control messages
              MPI_Get_count(&status, MPI_INT, &received_count);
              assert(received_count == 2);

              std::unique_ptr<ControlMessage> msg(new ControlMessage(TaggedEpoch(), src));

              // receive data or FOC messages.  use finishReceive to handle the
              // FOC messages and count decrement

              MPI_Irecv(&(msg->control), received_count, MPI_INT, src, tag, comm, &req);

              DEBUGF("C RECV %d -> %d, termination signal", src, commLayer.commRank);

              // insert into the received InProgress queue.
              auto pair = std::make_pair(std::move(req), std::move(msg));
              assert(msg.get() == nullptr);
              assert(pair.second.get() != nullptr);

//              commLayer.recvInProgress.push_back(std::move(pair));
//              assert(pair.second.get() == nullptr);
//              assert(commLayer.recvInProgress.back().second.get() != nullptr);

              if (!commLayer.recvInProgress.waitAndPush(std::move(pair)) ) {
                throw bliss::io::IOException("recvInProgress queue is disabled!");
              }


            } else {  // data messages
              MPI_Get_count(&status, MPI_UNSIGNED_CHAR, &received_count);

              uint8_t* msg_data = (received_count == 0) ? nullptr : new uint8_t[received_count];
              std::unique_ptr<DataMessageReceived> msg(new DataMessageReceived(msg_data, received_count, tag, src));
              // owner of msg-data is now msg.


              MPI_Irecv(msg_data, received_count, MPI_UNSIGNED_CHAR, src, tag, comm, &req);

              auto pair = std::make_pair(std::move(req), std::move(msg));
              assert(msg.get() == nullptr);
              assert(pair.second.get() != nullptr);

              // insert into the received InProgress queue.
//              commLayer.recvInProgress.push_back(std::move(pair));
//
//              assert(pair.second.get() == nullptr);
//              assert(commLayer.recvInProgress.back().second.get() != nullptr);

              if (!commLayer.recvInProgress.waitAndPush(std::move(pair)) ) {
                throw bliss::io::IOException("recvInProgress queue is disabled!");
              }


            }

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
          while(!commLayer.recvInProgress.isEmpty())
          {

//            assert(commLayer.recvInProgress.front().second.get() != nullptr);
            auto popped = std::move(commLayer.recvInProgress.waitAndPop());
            auto front = std::move(popped.second);  // this thread owns the element now.
            if (!popped.first) {
              // waitAndPop returned becuase recvInProgress is disabled and empty.  break from loop
              break;
            }

            if (front.first == MPI_REQUEST_NULL) {  // if local messages
              finished = 1;  // then this request is complete
            } else {
              // test if the MPI request for receive has been completed.
              MPI_Test(&(front.first), &finished, MPI_STATUS_IGNORE);
            }

            if (finished)
            {
              // finished request, so enqueue the message for processing.

//              // get the request
//              auto front = std::move(commLayer.recvInProgress.front());
//              // remove pending receive element
//              commLayer.recvInProgress.pop_front();

              // if it's a control message, then need to count down
              if (front.second.get()->tag == CONTROL_TAG.tag)  {

                ControlMessage* msg = dynamic_cast<ControlMessage*>(front.second.get());


                // get the message type being controlled, and the ep
                TaggedEpoch te = msg->control;

//                DEBUGF("C rank %d received control message for tag %d epoch %d.  recvRemaining has %ld tags", commRank, tag, epoch, recvRemaining.size());

                //== if the message type is same as CONTROL_TAG, then this is a Application Termination FOC message.
      //          DEBUGF("C rank %d triaged control message for tag %d epoch %d.  recvRemaining has %ld tags", commRank, tag, epoch, recvRemaining.size());

                // if we are not waiting for the tag, then we need to do something - the message has alreday been recieved vis MPI_Test,
                //  so we need to construct a new one and push back into the queue.
                // then skip over the remaining of this branch and leave the message in recvInProgress.
                if (commLayer.ctrlMsgRemainingForEpoch.find(te) == commLayer.ctrlMsgRemainingForEpoch.end()) {
                  //if (iters %100 == 0) WARNINGF("WARN: RANK %d not yet waiting for control message tag %d epoch %d.  cycle through again.", commLayer.commRank, te.tag, te.epoch);

                  front.first = MPI_REQUEST_NULL;
                  assert(front.second.get() != nullptr);
//                  commLayer.recvInProgress.push_back(std::move(front));
//                  assert(front.second.get() == nullptr);
//                  assert(commLayer.recvInProgress.back().second.get() != nullptr);
                  if (!commLayer.recvInProgress.waitAndPush(std::move(front)) ) {
                    throw bliss::io::IOException("recvInProgress queue is disabled!");
                  }


                  break;  // gives the worker threads a chance to call sendControlMessagesAndWait

                  // this delays the processing, which hopefully does not create deadlock.
                  // goal is for worker threads to be able to call sendControlMessagesAndWait to set up listening for a tag.
                  // and NOT use the comm thread to set up listening for a tag.  worker thread will wait, comm thread does not.
                  // using comm thread could cause consecutive "finish" calles to intermix their control messages at the remote receiver.
                }

                //== now decrement the count of control messages for a ep
//                DEBUGF("C RECV PRE rank %d receiving END signal tag %d epoch %d from %d, recvRemaining size: %ld",
//                     commLayer.commRank, te.tag, te.epoch, front.second->rank, commLayer.ctrlMsgRemainingForEpoch.size());
                --(commLayer.ctrlMsgRemainingForEpoch.at(te));
//                DEBUGF("C RECV rank %d receiving END signal tag %d epoch %d from %d, num senders remaining is %d",
//                     commLayer.commRank, te.tag, te.epoch, front.second->rank, commLayer.ctrlMsgRemainingForEpoch.at(te));


                //== if after decrement the count is 0 for the ep,
                if (commLayer.ctrlMsgRemainingForEpoch.at(te) == 0) {

                  //DEBUGF("ALL END received for tag %d, pushing to recv queue", front.second.tag);

//                  if (commLayer.recvQueue.isFull()) fprintf(stderr, "Rank %d recvQueue is full for control!!!\n", commLayer.commRank);

                  // and enqueue ONE FOC message for this tag to be handled by callback.  (reduction from commSize to 1)
                  if (!commLayer.recvQueue.waitAndPush(std::move(front.second))) {
                    throw bliss::io::IOException("C ERROR: recvQueue is not accepting new receivedMessage due to disablePush");
                  }

                } else if (commLayer.ctrlMsgRemainingForEpoch.at(te) < 0) {
                  // count decremented to negative.  should not happen.  throw exception.

                  std::stringstream ss;
                  ss << "C ERROR: number of remaining receivers for tag " << te.tag << " epoch " << te.epoch << " is now NEGATIVE";
                  throw bliss::io::IOException(ss.str());
                }
              } else {

 //               if (commLayer.recvQueue.isFull()) fprintf(stderr, "Rank %d recvQueue is full for data!!!\n", commLayer.commRank);

                //==== Data message.  add the received messages into the recvQueue
                if (!commLayer.recvQueue.waitAndPush(std::move(front.second))) {
                  throw bliss::io::IOException("C ERROR: recvQueue is not accepting new receivedMessage due to disablePush");
                }
              }
            }
            else
            {
              // not ready yet.  so push back in to front of queue
              if (!commLayer.recvInProgress.waitAndPushFront(std::move(front)) ) {
                throw bliss::io::IOException("recvInProgress queue is disabled!");
              }

              // the recv request is not done, so stop and wait for the next cycle.
              break;
            }
          }

          // see if we are done with all pending receives and all Application Termination FOC messages are received.
          if ((commLayer.ctrlMsgRemainingForEpoch.size() == 0)) { //&& !commLayer.recvInProgress.canPop()) {
            // if completely done, epoch=NO_EPOCH would have been erased from recvRemaining.
            recvDone.store(true);
            commLayer.recvQueue.disablePush();
          }
        }

    };

//TODO:  chagne recvInProgress to a threadsafe queue.


    // callback thread,  for organization.
    // TODO: group the code for callback thread.
    class CallbackThread
    {

      public:
        CallbackThread(CommunicationLayer& _comm_layer) :
          commLayer(_comm_layer), callbackFunctions(),
          td() {

        };
        CallbackThread() = delete;

        virtual ~CallbackThread() {
          callbackFunctions.clear();
        }

        void start() {
          if (!td.joinable()) td = std::move(std::thread(&bliss::io::CommunicationLayer::CallbackThread::run, this));
        }

        void finish() {
          if (td.joinable()) td.join();
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
          assert(tag != CONTROL_TAG.tag && tag >= 0);

          // one thread calls this at a time.
          std::unique_lock<std::mutex> lock(commLayer.mutex);

          if (callbackFunctions.find(tag) != callbackFunctions.end()) {
            lock.unlock();
            throw std::invalid_argument("M callback function already registered for given tag");
          }

          // add the callback function to a lookup table
          callbackFunctions[tag] = callbackFunction;

        }

      protected:
        /// parent (containing) class for this thread.
        CommunicationLayer& commLayer;

        /// Registry of callback functions, mapped to by the associated tags
        std::unordered_map<int,std::function<void(uint8_t*, std::size_t, int)> > callbackFunctions;

        /// std::thread object for the dedicated callback-handler thread
        std::thread td;

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
        void run()
        {
          DEBUGF("R THREAD STARTED:  recv-thread on %d", commLayer.commRank);

          // while there are messages to process, or we are still receiving
          while (commLayer.recvQueue.canPop() )  // same as !recvDone || !recvQueue.isEmpty()
          {
            // get next element from the queue, wait if none is available.
            // waitAndPop will exit out of wait when termination flag is set on the recvQueue
            auto result = std::move(commLayer.recvQueue.waitAndPop());
            if (!result.first) {
              // no valid result, we must be done with receiving
              assert(!commLayer.recvQueue.canPush());
              break;
            }

            // get the message (data or control)

            // if the message is a control message (enqueued ONE into RecvQueue only after all control messages
            //  for a tag have been received.)
            if (result.second.get()->tag == CONTROL_TAG.tag) {   // Control message

              ControlMessage* msg = dynamic_cast<ControlMessage*>(result.second.get());

              // get control message's target tag.
              TaggedEpoch te = msg->control;


              //===== and unblock the flush/finish/finishCommunication call.
              // getting to here means that all control messages for a tag have been received from all senders.
              // this means the local "sendControlMessages" has been called.  now, look at order of calls in waitForControlMessage
              //
              // main thread:  lock; flushing = notag->tag; flushing==tag?; cv_wait;
              // callback thread:  lock; flushing = tag->notag ? cv_notify_all : return message to front of queue
              // to ensure the flushing is modified by 1 thread at a time, and notify_all does not get interleaved to before wait,
              // lock is used to enclose the entire sequence on each thread.
              // also if callback thread's block is called before the main thread, we don't want to discard the control message - requeue.
              std::unique_lock<std::mutex> lock = std::move(commLayer.ctrlMsgProperties.at(te.tag).getUniqueLock());

              // if flushing is not tag, then some other tag is being processed, or flushing == notag.
              // if notag - this block is called before waitForControlMessages.  need to requeue
              // if some other tag (tag1), that means we received all tag2 messages in comm thread before receiving all tag1 messages
              // therefore recv thread is processing tag2 before tag1 control message.
              // if the delayed tag1 message is from remote process, and assume that process called sendControlMessages tag1 before tag2,
              //  then strict ordering of MPI messages is violated.
              // if the delayed tag1 message is from local process, because in mem message movement, tag1 is always before tag2.
              // so tag2 cannot come before tag1, provided that all MPI processes call sendControlMessage and waitForControlMessages in the same order
              //  including requirement that callers of sendControlMessage and waitForControlMessages reside in 1 thread only.

              // to deal with multiple threads that can use control messages, we'd need to have an array of "flushing" variables, array of mutices (sharing would render the threads serial),
              // and array of lock and condition variables (maybe)
              // also, it may be good to set "flushing" when first control message arrives.  perhaps recvRemaining is a good substitute.
              // finally, do we still need to have a UnPop function on the queue? - probably not.
              // even when multiple threads creating control messages, a tag can be handled by 1 thread at a time.

              //assert(recvRemaining.at(tag) == 0);
              commLayer.ctrlMsgRemainingForEpoch.erase(te);
              DEBUGF("R Rank %d tag %d epoch %d start waiting", commLayer.commRank, te.tag, te.epoch);


              // so unset 'flushing', and unblock flush() via condition variable.
              commLayer.ctrlMsgProperties.at(te.tag).notifyAll();

              DEBUGF("R rank %d tag %d epoch %d notified all.", commLayer.commRank, te.tag, te.epoch );
              lock.unlock();
            } else {
              DataMessageReceived* msg = dynamic_cast<DataMessageReceived*>(result.second.get());
              // data message.  process it.
              (callbackFunctions[msg->tag])(msg->data.get(), msg->count, msg->rank);
              // delete [] msg.data;  using unique_ptr, don't need to delete msg data.
            }
            // clean up message data
          }
          DEBUGF("R THREAD FINISHED: recv-thread on %d", commLayer.commRank);
        }
    };


public:
    // TODO: leave only code for CommLayer public api

  /**
   * @brief Constructor for the CommLayer.
   * @details have to have the comm_size parameter since we can't get it during member variable initialization and recvQueue needs it.
   *
   * @param communicator    The MPI Communicator used in this communicator
   *                        layer.
   * @param comm_size       The size of the MPI Communicator.
   */
  CommunicationLayer (const MPI_Comm& communicator, const int comm_size)
    : sendQueue(), recvQueue(), recvInProgress(), buffers(),
      ctrlMsgProperties(), ctrlMsgRemainingForEpoch(),
      sendThread(*this, communicator), recvThread(*this, communicator), callbackThread(*this)
  {
    // init communicator rank and size
    MPI_Comm_size(communicator, &commSize);
    assert(comm_size == commSize);
    MPI_Comm_rank(communicator, &commRank);

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
    // only 1 thread does this.
//    std::unique_lock<std::mutex> lock(mutex);

    DEBUGF("M FINALIZING COMMLAYER");

    // TODO: check that we are finished in finishCommunication (turn off sendQueue asap)
	  finishCommunication();

	  // TODO: properly handle this part.
    // wait for both threads to quit
    sendThread.finish();
	  recvThread.finish();
    callbackThread.finish();
  };


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
  void sendMessage(const void* data, std::size_t nbytes, int dst_rank, int tag)
  {
    //== don't send control tag.
    assert(tag != CONTROL_TAG.tag && tag >= 0);

    //check to see if sendQueue is still accepting (if not, then comm is finished)
    if (!sendQueue.canPush()) {
      throw std::logic_error("W sendMessage called after finishCommunication");
    }

    //== check to see if the target tag is still accepting
    std::unique_lock<std::mutex> lock(mutex);

    if (ctrlMsgProperties.find(tag) == ctrlMsgProperties.end() ||
        buffers.find(tag) == buffers.end()) {
      lock.unlock();
      throw std::invalid_argument("W invalid tag: tag has not been registered");
    }
    if (ctrlMsgProperties.at(tag).isFinished()) {
      lock.unlock();
      throw std::invalid_argument("W invalid tag: tag is already FINISHED");
    }
    lock.unlock();

    //== try to append the new data - repeat until successful.
    // along the way, if a full buffer's id is returned, queue it for sendQueue.
    BufferIdType fullId = BufferPoolType::ABSENT;
    std::pair<bool, BufferIdType> result;
    do {
      // try append.  append fails if there is no buffer or no room
      result = buffers.at(tag).append(data, nbytes, dst_rank);

      fullId = result.second;

      if (fullId != CommunicationLayer::BufferPoolType::ABSENT) {        // have a full buffer.  implies !result.first
        // verify that the buffer is not empty (may change because of threading)
        if (!(buffers.at(tag).getBackBuffer(fullId).isEmpty())) {
          // have a non-empty buffer - put in send queue.
        	DEBUGF("Rank %d has full buffer at id %d.  enqueue for send.", commRank, fullId);
          std::unique_ptr<SendDataElementType> msg(new SendDataElementType(fullId, tag, dst_rank));

//          if (sendQueue.isFull()) fprintf(stderr, "Rank %d sendQueue is full for data msgs\n", commRank);

          if (!sendQueue.waitAndPush(std::move(msg))) {
            throw bliss::io::IOException("W ERROR: sendQueue is not accepting new SendQueueElementType due to disablePush");
          } // else successfully pushed into sendQueue.
        }  // else empty back buffer

        // full buffer enqueued, or empty buffer, reset it.
        fullId = CommunicationLayer::BufferPoolType::ABSENT;
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
    assert(tag != CONTROL_TAG.tag && tag >= 0);

    // single thread calls this at a time.
    std::unique_lock<std::mutex> lock(mutex);

    //=======  reigster the tag.
    if (ctrlMsgProperties.find(tag) == ctrlMsgProperties.end()) {
      // register a new tag.
      ctrlMsgProperties[tag] = std::move(TagMetadata());
    }

    //== if there isn't already a tag listed, add the MessageBuffers for that tag.
    // multiple threads may call this.
    if (buffers.find(tag) == buffers.end()) {
      // create new message buffer
      buffers[tag] = std::move(MessageBuffersType(commSize, 8192, 2 * commSize));
    }
    lock.unlock();

    // now tell the callback thread to register the callbacks.
    callbackThread.addReceiveCallback(tag, callbackFunction);


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
    // multiple threads allowed

    // can only flush data tag.
    assert(tag != CONTROL_TAG.tag && tag >= 0);

    if (!sendQueue.canPush()) {
      throw std::logic_error("W flush called after finishCommunication");
    }

    std::unique_lock<std::mutex> lock(mutex);

    if ((ctrlMsgProperties.find(tag) == ctrlMsgProperties.end()) ||
        ctrlMsgProperties.at(tag).isFinished()) {
      lock.unlock();
      throw std::invalid_argument("W tag not registered, or already finished. cannot FLUSH");
    }
    lock.unlock();

    DEBUGF("M FLUSH Rank %d tag %d, buffer ids: %s", commRank, tag, buffers.at(tag).bufferIdsToString().c_str());

    // flush all data buffers (put them into send queue)
    flushBuffers(tag);

    // send the control message - generates an unique epoch for the tag as well.
    // wait till all control messages have been received (not using MPI_Barrier.  see Rationale.)
    DEBUGF("M FLUSH Rank %d tag %d, waiting for control messages to complete", commRank, tag);
    while (!sendControlMessagesAndWait(tag)) {};
    DEBUGF("M FLUSH DONE Rank %d tag %d, finished wait for control messages\n", commRank, tag);
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
    assert(tag != CONTROL_TAG.tag && tag >= 0);

    if (!sendQueue.canPush()) {
      //throw std::logic_error("W finish called after finishCommunication");
      return;
    }

    std::unique_lock<std::mutex> lock(mutex);

    if (ctrlMsgProperties.find(tag) == ctrlMsgProperties.end()) {
      lock.unlock();
      throw std::invalid_argument("W tag not registered. cannot FINISH");
    }
    if (ctrlMsgProperties.at(tag).isFinished()) {
      lock.unlock();
      //throw std::invalid_argument("W tag not registered, or already finished. cannot FINISH");
      return;
    }
    // mark as no more coming in for this tag.  ONE THREAD ONLY modifies this.
    ctrlMsgProperties.at(tag).finish();

    lock.unlock();

    DEBUGF("M FINISH Rank %d tag %d, buffer ids: %s", commRank, tag, buffers.at(tag).bufferIdsToString().c_str());

    // flush all buffers (put them into send queue)
    flushBuffers(tag);

    // send the control message
    // wait till all control messages have been received (not using MPI_Barrier.  see Rationale.)
    DEBUGF("M FINISH Rank %d tag %d, waiting for control messages to complete", commRank, tag);
    while (!sendControlMessagesAndWait(tag)) {};
    DEBUGF("M FINISH Rank %d tag %d, finished wait for control messages", commRank, tag);

  }


  /**
   * @brief Initializes all communication and starts the communication and
   *        callback threads.
   *
   * @note Must be called before any other functions can be called.
   */
  void initCommunication()
  {
    // one thread calls this
    if (!sendQueue.canPush()) {
      throw std::logic_error("W initCommunication called after finishCommunication");
    }

    std::unique_lock<std::mutex> lock(mutex);

    if (ctrlMsgProperties.find(CONTROL_TAG.tag) != ctrlMsgProperties.end()) {
      lock.unlock();
      throw std::logic_error("W initCommunication called multiple times");
      //return;
    }

    //== init with control_tag only.
    ctrlMsgRemainingForEpoch[CONTROL_TAG] = commSize;
    ctrlMsgProperties[CONTROL_TAG.tag] = std::move(TagMetadata());

    // TODO: deal with where to put the std::threads.
    callbackThread.start();
    recvThread.start();
    sendThread.start();
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
    DEBUGF("M FINISH COMM rank %d", commRank);


	  // already finishing.
    if (!sendQueue.canPush()) {
      return;
      //throw std::logic_error("W finishCommunication called after finishCommunication");
    }

    // one thread does one finish at a time.   not deleting entries from ctrlMsgProperties, so no need to lock whole thing.
    //
	  for (auto ctrlMsgIter = ctrlMsgProperties.begin(); ctrlMsgIter != ctrlMsgProperties.end(); ++ctrlMsgIter) {
      if (ctrlMsgIter->first != CONTROL_TAG.tag && ctrlMsgIter->first >= 0) {
        finish(ctrlMsgIter->first);
      }
    }

    // send the control message
    // wait till all control messages have been received (not using MPI_Barrier.  see Rationale.)
    DEBUGF("M FINISH COMM Rank %d tag %d, waiting for CONTROL messages to complete", commRank, CONTROL_TAG.tag);
    while (!sendControlMessagesAndWait(CONTROL_TAG.tag)) {};
    DEBUGF("M FINISH COMM Rank %d tag %d, finished wait for CONTROL messages", commRank, CONTROL_TAG.tag);

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


  bool isTagActive(int tag) {
    // single threaded
    std::unique_lock<std::mutex> lock(mutex);

    return isTagActive(tag, lock);
  }
  bool isTagActive(int tag, std::unique_lock<std::mutex>& lock) {

    //== check to see if the target tag is still accepting
    if (ctrlMsgProperties.find(tag) == ctrlMsgProperties.end()) {
      throw std::invalid_argument("W invalid tag: tag has not been registered.");
    }
    if (ctrlMsgProperties.at(tag).isFinished()) {
      throw std::invalid_argument("W invalid tag: tag is already FINISHED");
    }
    return true;
  }

  /**
   * @brief Sends FOC message for the given tag to every MPI process
   *        Blocks till the flusing of the given tag is completed (i.e., all
   *        FOC messages have been received).
   *
   * @param tag     The MPI tag to end.
   * @throw bliss::io::IOException  If the sendQueue has been disabled.
   * @return bool indicating if control messages were successfully sent.
   */
  bool sendControlMessagesAndWait(int tag) throw (bliss::io::IOException)
  {
    if (!sendQueue.canPush())
      throw bliss::io::IOException("M ERROR: sendControlMessagesAndWait already called with CONTROL_TAG");



    // since each tag has unique epoch for each call to this function, we don't need global lock as much.
    TaggedEpoch te{ctrlMsgProperties.at(tag).getNextEpoch(), tag};

    std::unique_lock<std::mutex> lock = std::move(ctrlMsgProperties.at(tag).getUniqueLock());
    if (tag != CONTROL_TAG.tag) {
      // not a control tag, so don't need to add to recvRemaining.

      // lock it for the tag - no multiple concurrent calls with the same tag.
      ctrlMsgRemainingForEpoch.emplace(te, commSize);  // insert new, and proceed to sending messages
      DEBUGF("M Rank %d tag %d epoch %d added tag to recvRemaining.", commRank, te.tag, te.epoch);
    } // else there is only one epoch for CONTROL_TAG, and it's already added to the list.


    // should be only 1 thread executing this code below with the provided tag.

    for (int i = 0; i < getCommSize(); ++i) {
      // send end tags in circular fashion, send tag to self last
      int target_rank = (i + getCommRank() + 1) % getCommSize();

      //DEBUGF("M Rank %d sendEndTags %d epoch %d target_rank = %d ", commRank, te.tag, te.epoch, target_rank );

      // send the end message for this tag.
      std::unique_ptr<ControlMessage> msg(new ControlMessage(te, target_rank));


     // if (sendQueue.isFull()) fprintf(stderr, "Rank %d sendQueue is full for control msgs\n", commRank);
      if (!sendQueue.waitAndPush(std::move(msg))) {
        lock.unlock();
        throw bliss::io::IOException("M ERROR: sendQueue is not accepting new SendQueueElementType due to disablePush");
      }
    }


    // if we are sending the control tag, then afterward we are done.
    if (tag == CONTROL_TAG.tag)
      sendQueue.disablePush();
    lock.unlock();


    //============= now wait.

    // CONCERN:  tracking the tag being flushed using a shared "flushing" opens up the possible problem
    // when out-of-order calls to flush, finish, and finishCommunications.


    // waiting for all messages of this tag to flush.
    //DEBUGF("M Rank %d Waiting for tag %d", commRank, te.tag);

    lock.lock();

    // TODO:  does this need to be locked?  should recvRemaining be for epoch or tag?


    while (ctrlMsgRemainingForEpoch.find(te) != ctrlMsgRemainingForEpoch.end()) {
      // condition variable waits for notification, which is generated by recv thread
      // when all FOC messages with this tag are received.
      DEBUGF("M PRE WAIT Rank %d Waiting for tag %d epoch %d, currently flushing", commRank, te.tag, te.epoch);
      ctrlMsgProperties.at(tag).wait(lock);
      DEBUGF("M END WAIT Rank %d Waiting for tag %d epoch %d, currently flushing", commRank, te.tag, te.epoch);
    }
    lock.unlock();

    DEBUGF("M Rank %d received all END message for TAG %d, epoch = %d", commRank, te.tag, te.epoch);
    return true;
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

    if (!sendQueue.canPush()) {
      throw std::logic_error("W cannot flush buffer: cannot push to sendQueue.");
    }

//    // no MessageBuffers with the associated tag.  end.
//    if (buffers.find(tag) == buffers.end())
//    {
//      DEBUGF("NO BUFFERS for tag: %d", tag);
//      return;
//    }
//    //DEBUGF("Active Buffer Ids in message buffer: %s", buffers.at(tag).bufferIdsToString().c_str());

    std::unique_lock<std::mutex> lock = std::move(ctrlMsgProperties.at(tag).getUniqueLock());

    // flush out all the send buffers matching a particular tag.
    int idCount = buffers.at(tag).getBufferIdsForAllRanks().size();
    assert(idCount == commSize);


    for (int i = 0; i < idCount; ++i) {
      // flush buffers in a circular fashion, starting with the next neighbor
      int target_rank = (i + getCommRank() + 1) % getCommSize();

      auto id = buffers.at(tag).flushBufferForRank(target_rank);
      // flush/send all remaining non-empty buffers
      if ((id != BufferPoolType::ABSENT) && !(buffers.at(tag).getBackBuffer(id).isEmpty())) {
        std::unique_ptr<SendDataElementType> msg(new SendDataElementType(id, tag, target_rank));

//        if (sendQueue.isFull()) fprintf(stderr, "Rank %d sendQueue is full for flushBuffers\n", commRank);

        if (!sendQueue.waitAndPush(std::move(msg))) {
          lock.unlock();
          throw bliss::io::IOException("M ERROR: sendQueue is not accepting new SendQueueElementType due to disablePush");
        }
      }
    }

    lock.unlock();
  }




protected:
/******************************************************************************
 *                           Protected member fields                            *
 ******************************************************************************/


  /* Thread-safe queues for posting sends and getting received messages,
   * these are used as communication medium between the comm-thread and the
   * (potetially multiple) producer threads for the sends; and the comm-thread
   * and the callback handler thread for the receives
   */

  /// Mutex for adding new tags to the buffers map (since inserts can be multi-threaded)
  mutable std::mutex mutex;

  /// atomic variable indicating if the CommLayer has been finalized (all communications are completed, ready for destruction.)
//  std::atomic<bool> finalized;


  // Outbound message data structure, Multiple-Producer-Single-Consumer queue
  // consumed by the internal MPI-comm-thread
  /// message queue between sendMessage calling threads (src) and comm-thread
  /// (sink)
  bliss::concurrent::ThreadSafeQueue<std::unique_ptr<MPIMessage> > sendQueue;

  // Inbound message data structure, Single-Producer-Multiple-Consumer queue
  // produced by the internal MPI-comm-thread
  /// message queue between comm thread (src) and callback thread (sink)
  bliss::concurrent::ThreadSafeQueue<std::unique_ptr<MPIMessage> > recvQueue;


  /* Queues of pending MPI operations (MPI_Requests)
   * Used ONLY by the comm-thread, purely internal and purely single thread.
   */
  /// Queue of pending MPI receive operations
  bliss::concurrent::ThreadSafeQueue<std::pair<MPI_Request, std::unique_ptr<MPIMessage> > > recvInProgress;


  /* information per message tag */

  // TODO: comments about when they are used.

  /**
   *  Message buffers per message tag (maps each tag to a set of buffers, one buffer for each target rank)
   *  only destroyed at the end, after all sends are done.
   */
  std::unordered_map<int, MessageBuffersType> buffers;

  /**
   *  map between tag and TagMetadata, which contains info aobut next epoch id,
   *  mutex, and condition variables.  Also has flag to indicate if tag is "finished"
   *
   *  use the flag instead of deletiing the TagMetadata because
   *
   *  Ensures that a single thread during sendControlMessageAndWait is changing the condition for the wait.
   *
   */
  std::unordered_map<int, TagMetadata > ctrlMsgProperties;

  /// Per tag: number of processes that haven't sent the END-TAG message yet
  std::unordered_map<TaggedEpoch, int, TaggedEpochHash, TaggedEpochEqual > ctrlMsgRemainingForEpoch;

  SendCommThread sendThread;
  RecvCommThread recvThread;

  CallbackThread callbackThread;

  /// The MPI Communicator size
  int commSize;

  /// The MPI Communicator rank
  int commRank;

};



//==== static variable definitions.

/// CONTROL_TAG definition.  (declaration and initialization inside class).
constexpr CommunicationLayer::TaggedEpoch CommunicationLayer::CONTROL_TAG;

/// NO_TAG definition.  (declaration and initialization inside class).
const int CommunicationLayer::NO_TAG;

/// NO_EPOCH definition.  (declaration and initialization inside class).
const int CommunicationLayer::NO_EPOCH;


} // namespace io
} // namespace bliss

#endif // BLISS_COMMUNICATION_LAYER_HPP
