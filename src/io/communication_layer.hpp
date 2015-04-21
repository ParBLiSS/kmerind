/**
 * @file    communication_layer.hpp
 * @ingroup bliss::io
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @author  Tony Pan <tpan7@gatech.edu>
 * @brief   contains a class that abstracts MPI point-to-point, non-blocking messaging.
 *
 * Copyright (c) 2015 Georgia Institute of Technology
 *
 * TODO add Licence
 */

#ifndef BLISS_COMMUNICATION_LAYER_HPP
#define BLISS_COMMUNICATION_LAYER_HPP

#include "config/relacy_config.hpp"

#include <mpi.h>
#include <omp.h>

// C stdlib includes

// STL includes
#include <queue>
#include <unordered_map>
#include <thread>
#include <functional>
#include <utility>    // std::move

// BLISS includes
#include "utils/logging.h"
#include "concurrent/mutexlock_queue.hpp"
#include "concurrent/copyable_atomic.hpp"
#include "concurrent/concurrent.hpp"
#include "concurrent/referenced_object_pool.hpp"
#include "io/io_exception.hpp"
#include "io/message_buffers.hpp"
#include "io/message_type_info.hpp"
#include "io/message_types.hpp"

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
 *          FOC messages have tag 0, but carries in payload the tag of the message type that is being controlled.
 *          in addition, FOC messages have an additional, sequentially assigned id for each message type to
 *          disambiguate multiple FOC messages with the same tag
 *
 *        2. support collective, blocking flushing of all messages of a particular type.
 *          To indicate to the system that there is no further messages of a particular type, collective flush
 *          operations have been implemented.  There are 3 types:
 *          a. "flush":  flush a single message type, allow future messages of this type
 *          b. "finish": flush a single message type, disallow future messages of this type
 *          c. "finish communication":  flush all message types, and mark end of application communication)
 *
 *          As these are collective functions, all MPI processes need to make matching calls.
 *          The timing of the calls need not be synchronous, but they need to preserve the order of the calls
 *
 *          The functions form barriers, and blocks the calling thread
 *          until all processes have completely processed all matching FOC messages from all senders
 *
 *        3. provide flexible way of handling received messages based on message type
 *          Callbacks can be registered, one per message type.  When processing received messages,
 *          the matching callback is used.
 *
 *        4. Thread Safe, non-blocking send (allow overlap of computation with communication)
 *          A compute thread can call sendMessage on the Communication Layer object.  each invocation
 *          atomically gets a location to insert the message into a batching buffer via thread safe insert/
 *
 *          When a buffer is full and sent, it is inserted into a thread-safe send queue for the sending comm thread to
 *          send via MPI  (non-blocking).
 *
 *        5. Thread Safe, non-blocking receive
 *          A dedicated thread checks for MPI messages to receive remote messages and populate the recvQueue. Separation
 *          of sender and receiver communication threads avoids deadlock when the receiving process is blocked by a full
 *          recvQueue, thus blocking the sending function.
 *
 *        6. Thread Safe processing of received messages
 *          Multiple callback threads can process the received messages in the receive queue, which is thread safe.
 *
 *        7. single threaded control of the Communication Layer (flush messages, finish, etc)
 *          flush and finish functions are collectively called by the control thread (owner of Communication Layer object)
 *
 *        Basic Concepts and Design:
 *        A. Threading
 *          The following are the intended threading pattern:
 *          * single control thread (owner of Communication Layer instance.  can be one of the compute threads)
 *          * one or more compute threads (producer of messages to be sent)
 *          * single MPI communication thread for send and recv (internal)
 *          * one or more receive callback threads (consumers of messages received, internal)
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
 *          not using MPI AllToAllv because that blocks the compute threads
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
 *          control thread call finishCommunication() (TO MAKE SURE ALL THREADS ARE FINISHED)
 *
 *
 * @note:  COMMLAYER USES THREADLOCAL BUFFERS.  Shared BUFFER requires coordinate buffer ptr swap, which is not really safe.
 *
 */
template<bool ThreadLocal = true>
class CommunicationLayer
{
	 static_assert(ThreadLocal, "Communication Layer only supports ThreadLocal mode. (template param set to true)");

protected:

	 using BufferType = bliss::io::Buffer<bliss::concurrent::LockType::NONE, 1048576, sizeof(uint64_t)>;

	   /// alias BufferPoolType to MessageBuffersType's
	   using BufferPoolType = bliss::concurrent::ObjectPool< BufferType, bliss::concurrent::LockType::MUTEX>;


	 // set pool locktype to mutex - for testing with thread sanitizer (suspect spinlock does not work with threadsanitizer well.)
   /// alias MessageBuffersType to Use the thread safe version of the SendMessageBuffers
   using MessageBuffersType = SendMessageBuffers<bliss::concurrent::LockType::THREADLOCAL, BufferPoolType>;



   /// alias BufferPtrType to MessageBuffersType's
   using BufferPtrType = BufferType*;

   /// send MPI Message type aliases for convenience.
   using SendDataElementType = typename bliss::io::DataMessageToSend<BufferPtrType>;


   /**
    * @brief    communication thread to handle all send the receive.
    * @details  combined send and receive thread.
    *           1. supports single threaded MPI library.
    *           2. need to use unlimited send/recv queue, else either blocking would cascade to all MPI processes.
    *
    * TODO: mixing BSend instead of iSend MAY cause message ordering issues?
    */
   template <typename BufferPtrType>
   class CommThread
   {
     protected:
       /// disabled default constructor
       DELETED_FUNC_DECL(CommThread());

       /// The MPI Communicator object for this communication layer
       MPI_Comm comm;

       /// the parent communication layer object
       CommunicationLayer& commLayer;


       /// set to true after all sends are completed
       bool sendDone;

       /// set to true if `finishing` and if all receives are completed.
       bool recvDone;


       /// Queue of pending MPI send operations
       std::deque< std::pair<MPI_Request, MPIMessage* > > sendInProgress;

       /// Queue of pending MPI receive operations
       std::deque< std::pair<MPI_Request, MPIMessage* > > recvInProgress;

       /// actual thread
       std::thread td;
       /// MPI comm rank
       int commRank;
       /// MPI communicator size
       int commSize;

       int lastSrcRank;

     public:

       /// default constructor
       CommThread(const MPI_Comm& communicator, CommunicationLayer & comm_layer) :
         comm(communicator), commLayer(comm_layer),
         sendDone(false),
         recvDone(false),
         sendInProgress(),
         recvInProgress(),
         td(), commRank(0), commSize(1), lastSrcRank()
       {
         MPI_Comm_size(comm, &commSize);
         MPI_Comm_rank(comm, &commRank);

         lastSrcRank = std::rand() % commSize;
       };

       /// default destructor
       virtual ~CommThread() {
       }

       /// start the comm thread and also allocate some buffer for sending control messages.
       void start() {

         if (!td.joinable()) {

           //===== initialize the mpi buffer for control messaging, for MPI_Bsend
           int mpiBufSize = 0;
           // arbitrary, 4X commSize, each has 2 uint64_t, 1 int
           MPI_Pack_size(commSize * 4 * (sizeof(uint64_t) * 2), MPI_UNSIGNED_CHAR, comm, &mpiBufSize);
           mpiBufSize += commSize * 4 * MPI_BSEND_OVERHEAD;
           char* mpiBuf = (char*)malloc(mpiBufSize);
           MPI_Buffer_attach(mpiBuf, mpiBufSize);

           td = std::move(std::thread(&bliss::io::CommunicationLayer<ThreadLocal>::CommThread<BufferPtrType>::run, this));
         } // else already started the thread.
       }

       /// finished the comm thread by waiting for the thread to join.
       void finish() {
         if (td.joinable()) {
           td.join();

           DEBUGF("C CommThread finished, thread joined.  about to detach mpi buffer.");

           // detach the buffer for Bsend.
           int mpiBufSize;
           char* mpiBuf;
           MPI_Buffer_detach(&mpiBuf, &mpiBufSize);
           free(mpiBuf);
         } // else already joined.
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
       void run()
       {
         DEBUGF("C THREAD STARTED:  send comm-thread on %d", commRank);

         bool worked = false;
         // while there's still work to be done:
         while (!sendDone || !recvDone)
         {
           worked = false;
           // first clean up finished operations
           worked |= finishSends();
           worked |= finishReceives();

           // start pending sends
           worked |= tryStartSend();

           // start pending recvs
           worked |= tryStartReceive();

           if (!worked) _mm_pause();
         }
         DEBUGF("C THREAD FINISHED:  send comm-thread on %d", commRank);
       }




       /**
        * @brief Takes queued messages and initiates MPI_ISends for them.
        *
        * @throw bliss::io::IOException
        */
       bool tryStartSend() throw (bliss::io::IOException)
       {
         // quit if no more sends are to be expected
         if (sendDone) return false;

         if (commLayer.sendQueue.isEmpty()) {
           _mm_pause();
           return false;
         }

         // try to get the send element
         bool suc;
         MPIMessage* ptr;
         std::tie(suc, ptr) = commLayer.sendQueue.tryPop();  // type: std::pair<bool, std::unique_ptr<MPIMessage>>
         if (!suc) {  // no valid entry in the queue
           return false;
         }
         // else there is a valid entry from sendQueue
         assert(ptr);

         bool worked = false;
         std::atomic_thread_fence(std::memory_order_acquire);

         // second is a pointer to ControlMessage
         if (ptr->getTag() == CONTROL_TAG) {
           // ===== if control message,

           // ControlMessage ownership is maintained by el.
           ControlMessage* msg = dynamic_cast<ControlMessage*>(ptr);  // MPIMessage*

           DEBUGF("C R %d -> %d SEND CONTROL for tag %d epoch %lu, next epoch %lu, sendInProgress size %lu",
                 commRank, msg->getRank(), msg->getTag(), msg->getEpoch(), msg->getNextEpoch(), sendInProgress.size());

           // send termination message for this tag and destination
//           if (se->rank == commRank) {
//             // =========== send locally
//
//             // local, directly handle by creating an byte array and directly
//             // insert into the recvInProgress (needed by the receiver decrement
//             // logic in "finishReceive", so can't put into receive queue.  also to maintain message orders)
//             MPI_Request req = MPI_REQUEST_NULL;
//
//             // now move ControlMessage ownership to recvInProgress
//             auto pair = std::make_pair(std::move(req),
//                                        std::move(ptr));  // message from send queue is a ControlMessage, which is compatible with recvInProgress queue.
//             recvInProgress.push_back(std::move(pair));
//
//           } else {


            // remote control message.  send a terminating message with tag CONTROL_TAG
            // payload is the message type and epoch.  Message Ordering is guaranteed by MPI.
            // Bsend to allow use of preallocated buffer, so don't have to clean up allocated messages.

           if (msg->getRank() == MPIMessage::ALL_RANKS) {
             int r = commRank;
             for (int i = 0; i < commSize; ++i) {
               r = (r + 1) % commSize;
               MPI_Bsend(msg->getData(), msg->getDataSize(), MPI_UNSIGNED_CHAR, r, CONTROL_TAG, comm);

               DEBUGF("C R %d -> %d Send Control for epoch %lu, nextEpoch %lu", commRank, r, msg->getEpoch(), msg->getNextEpoch());
             }
           } else {
             MPI_Bsend(msg->getData(), msg->getDataSize(), MPI_UNSIGNED_CHAR, msg->getRank(), CONTROL_TAG, comm);
             DEBUGF("C R %d -> %d Send single Control for epoch %lu, nextEpoch %lu", commRank, msg->getRank(), msg->getEpoch(), msg->getNextEpoch());
           }

//             // instead of Bsend, push into sendInProgress, and then there, do barrier.
//             MPI_Request req = MPI_REQUEST_NULL;
//             sendInProgress.push_back(std::move(std::make_pair(std::move(req), std::move(ptr))));

//           }  //if send to self.
           delete ptr;   // delete the previously allocated ControlMessage instance.
           worked = true;
         } else {
           // data message.

           SendDataElementType* msg = dynamic_cast<SendDataElementType*>(ptr);  // MPIMessage*

           // get message data and it's size
           size_t count = msg->getDataSize();

           if (count > 0) {
             std::atomic_thread_fence(std::memory_order_acquire);
//             if (se->rank == commRank) {
//               // local, directly handle by creating an output object and directly
//               // insert into the recv queue (thread safe)
//
//               // since the send buffer memory is managed by Buffer, and the receiveMessage expects to manage memory as well, need to copy.
//               std::unique_ptr<uint8_t[]> array(new uint8_t[count]);
//               memcpy(array.get(), data, count);
//
//               // pointer here to allow dynamic cast.
//               std::unique_ptr<DataMessageReceived> dmr(new DataMessageReceived(std::move(array), count, se->tag, commRank));
//
//
//               // put into recvInProgress queue instead of recvQueue to minimize out of order (early) receive.
//               MPI_Request req = MPI_REQUEST_NULL;
//
//               auto pair = std::make_pair(std::move(req), std::move(dmr));
//
//               recvInProgress.push_back(std::move(pair));
//
//               // finished inserting directly to local RecvInProgress.  release the buffer (since we have a copy)
//               commLayer.getTagInfo(se->tag).getBuffers()->releaseBuffer(std::move(se->ptr));
//             } else {

               // remote: initiate async MPI message
               MPI_Request req;
               MPI_Isend(msg->getPayload(), msg->getPayloadSize(), MPI_UNSIGNED_CHAR, msg->getRank(), msg->getTag(), comm, &req);


               sendInProgress.emplace_back(req, ptr);   // move the SendDataElement to sendInProgress.

//             } // if send to self.
             worked = true;
           }

         }
         return worked;
       }



       /**
        * @brief Initiates async MPI receives.
        */
       bool tryStartReceive()
       {

         // if no more are expected, then done
         if (recvDone) return false;

         // probe for messages
         int hasMessage = 0;
         MPI_Status status;
         //if (commSize > 1)  // iprobe checks low rank first.
         MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, comm, &hasMessage, &status);
//         int srcRank = lastSrcRank;
//         do {
//           MPI_Iprobe(srcRank, MPI_ANY_TAG, comm, &hasMessage, &status);  // fair probe.
//           srcRank = (srcRank + 1) % commSize;
//         } while ((hasMessage == 0) && (srcRank != lastSrcRank));
//         lastSrcRank = srcRank;

         // if have message to receive,
         if (hasMessage > 0) {
           // get some message details
           int src = status.MPI_SOURCE;
           int tag = status.MPI_TAG;
           int received_count;
           MPI_Request req;


           // receive the message data bytes into a vector of bytes
           if (tag == CONTROL_TAG) {  // control messages
             MPI_Get_count(&status, MPI_UNSIGNED_CHAR, &received_count);
             assert(received_count == 16);

             ControlMessage* msg = new ControlMessage(0, 0, src);  // create a new instance

             // receive data or FOC messages.  use finishReceive to handle the
             // FOC messages and count decrement.
             MPI_Irecv(msg->getData(), received_count, MPI_UNSIGNED_CHAR, src, tag, comm, &req);

             DEBUGF("C R %d <- %d RECEIVING CONTROL\ttag %d, epoch %lu, next epoch %lu. recvInProgress size %lu", 
                      commRank, src, msg->getTag(), msg->getEpoch(),msg->getNextEpoch(), recvInProgress.size());
             fflush(stdout);

             // insert into the received InProgress queue.
             recvInProgress.emplace_back(req, msg);


           } else {  // data messages
             MPI_Get_count(&status, MPI_UNSIGNED_CHAR, &received_count);


             // get a buffer from pool
             BufferType* ptr = commLayer.acquireBuffer();
             // ensure capacity is okay
             if(static_cast<size_t>(received_count) > (ptr->getCapacity() + ptr->getMetadataSize())) {
               ERRORF("C R %d <- %d Recvd count = %d, capacity = %ld", commRank, src, received_count, ptr->getCapacity());
               assert(static_cast<size_t>(received_count) <= (ptr->getCapacity() + ptr->getMetadataSize()));
             }
             ptr->reserve(received_count - ptr->getMetadataSize(), sizeof(uint8_t));  // reserve the data size as well.
             ptr->complete_write(received_count - ptr->getMetadataSize());
             ptr->block_writes();  // and then block future writes (reservation)
             //ptr->finish_writes();

             // create a data message
             SendDataElementType *msg = new SendDataElementType(ptr, tag, src);

             // have to receive even if count is 0.
             MPI_Irecv(ptr->metadata_begin(), received_count, MPI_UNSIGNED_CHAR, src, tag, comm, &req);

             DEBUGF("C R %d <- %d RECEIVING DATA\ttag %d", commRank, src, tag);

             // insert into the received InProgress queue.               // owner of msg_data is now msg
             recvInProgress.emplace_back(req, msg);

           }

           return true;
         }
         else return false;
       }

 
       /**
        * @brief Finishes pending MPI_Send requests.  stops at first unfinished request,
        *        leave unfinished requests to finish in next iteration's call.
        */
       bool finishSends()
       {
         // no more messages to send or being sent, so done.
         if (sendDone) return false;

         int mpi_finished = 0;

         bool worked = false;
         // while there is some in progress MPI_send requests,  don't use loop - delay here propagates to all nodes.
         if(!sendInProgress.empty())
         {

      	   MPI_Request req;
		 MPIMessage* ptr;
		 std::tie(req, ptr) = sendInProgress.front();


           // check if MPI request has completed.
           if (req == MPI_REQUEST_NULL) {
             // this is a control message since it has no "request".
             mpi_finished = 1;
           } else {
             // has request.  check for finished
             MPI_Test(&(req), &mpi_finished, MPI_STATUS_IGNORE);
           }

           // get the first request to check - ONLY 1 THREAD CHECKING sendInProgress.
           assert(ptr != nullptr);

           // if finished, then clean up.
           if (mpi_finished)
           {
        	   // since finished, get the first entry.
        	   sendInProgress.pop_front();


        	   if (ptr) {

             // if there is a buffer, then release it.
             if (ptr->getTag() != CONTROL_TAG) {
               DEBUGF("R %d <- %d FINISHED SEND DATA for tag %d, sendInProgress size %lu", 
                    commRank, ptr->getRank(), ptr->getTag(), sendInProgress.size());

               SendDataElementType* msg = dynamic_cast<SendDataElementType*>(ptr);
               if (msg->getBuffer()) {
                 // cleanup, i.e., release the buffer back into the pool
                 commLayer.releaseBuffer(msg->getBuffer());
               }
//             }  else {
//               ControlMessage* msg = dynamic_cast<ControlMessage*>(sendInProgress.front().second.get());
//               // control tag. introduce a barrier, and push into recvInProgress.
//               DEBUGF("R %d ENTER BARRIER FOR FLUSH for tag %d, epoch %u, sendInProgress size %lu", commRank, getTagFromTaggedEpoch(msg->tagged_epoch), getEpochFromTaggedEpoch(msg->tagged_epoch), sendInProgress.size());
//               MPI_Barrier(comm);
//               recvInProgress.push_back(std::move(sendInProgress.front()));

//             } else { // control tag, and mpi send finished, so nothing to clean up.
//               ERRORF("R %d <- %d SHOULD NOT HAVE CONTROL for tag %d epoch %lu, sendInProgress size %lu", commRank,
//            		   sendInProgress.front().second->getRank()
//                   sendInProgress.front().second->getTag(),
//            		   sendInProgress.front().second->getEpoch(),
//            		   sendInProgress.size());
             }
             delete ptr; // clean up the SendDataElement instance.  buffer has been released back to pool.
        	   }
             worked = true;
           }
           else
           {
             // the head of the queue is not finished, so break and get it later.
             worked = false;
           }
         }

         // all sending is done when "finishing" is set, and no more messages are pending send,
         // and no messages are being sent.
         if (!commLayer.sendQueue.canPop() && sendInProgress.empty()) {
            DEBUGF("sendQueue done");

           sendDone = true;
         }
         return worked;
       }


       /**
        * @brief Finishes pending receive requests from iRecv.
        * @details   if control message, count down from the total number of possible senders for an epoch
        *            when reaching 0, then enqueue 1 control message for callback handling to process
        *
        *            if data message, then enqueue for callback handling to process directly.
        *
        * @throw bliss::io::IOException
        */
       bool finishReceives() throw (bliss::io::IOException)
       {
         // no more to receive or pending receive requests, so done.
         if (recvDone) return false;

         int mpi_finished = 0;

         bool worked = false;

         // if there are pending receive requests
         if(!recvInProgress.empty())
         {
           //== get the next request.
           assert(recvInProgress.front().second != nullptr);

           //== check if MPI request has completed.
           if (recvInProgress.front().first == MPI_REQUEST_NULL) {
             // local message, so request completed.
             mpi_finished = 1;
           } else {
             // has request.  check for finished
             MPI_Test(&(recvInProgress.front().first), &mpi_finished, MPI_STATUS_IGNORE);
           }

           if (mpi_finished)
           {
             //== finished request, so enqueue the message for processing.
        	   MPI_Request req;
               MPIMessage* ptr;
               std::tie(req, ptr) = recvInProgress.front();
        	   recvInProgress.pop_front();


             // if it's a control message, then count down
             if (ptr->getTag() == CONTROL_TAG)  {

            	 DEBUGF("B R %d <- %d CONTROL %d epoch %lu",
            			 commLayer.commRank, ptr->getRank(), ptr->getTag(),
            			 ptr->getEpoch());

               // get the control message
               ControlMessage* msg = dynamic_cast<ControlMessage*>(ptr);
          	 DEBUGF("B R %d <- %d CONTROL cast %d epoch %lu, nextEpoch %lu",
          			 commLayer.commRank, msg->getRank(), msg->getTag(),
          			 msg->getEpoch(), msg->getNextEpoch());

               // counts down.  note that if this is a new epoch, a new entry will be created and activeEpochCount will be incremented.
          	 int v;
          	 if (commRank == 0) {  // head node.  count down by 1
               v = commLayer.epochProperties.countdownEpoch(msg->getEpoch(), std::ref(commLayer.activeEpochCount));
               //   doing it this way means we can process the received request right away (in message recv order)

               // if received all, then send messages back to other nodes.
               if (v == 0) {
                 for (int i = 1; i < commSize; ++i) {
                   MPI_Bsend(msg->getData(), msg->getDataSize(), MPI_UNSIGNED_CHAR, i, CONTROL_TAG, comm);
                 }
               }

          	 } else {
          	   // other nodes:  message from head node, so count down the whole thing.
          	   v = commLayer.epochProperties.countdownEpoch(msg->getEpoch(), std::ref(commLayer.activeEpochCount), commSize);
          	 }

            	 DEBUGF("B R %d <- %d CONTROL countdown %d epoch %lu, nextEpoch %lu, count %d",
            			 commLayer.commRank, msg->getRank(), msg->getTag(),
            			 msg->getEpoch(), msg->getNextEpoch(), v);

               //== if after countdown the count is 0 for the epoch, then the TaggedEpoch is done. enqueue a message for Callback Thread.
               if (v == 0) {

                 DEBUGF("R %d <- %d FINISHING RECV CONTROL for tag %d, epoch %lu, next epoch %lu, recvInProgress size %lu, recvQueue size %ld",
                		 commLayer.commRank, msg->getRank(), msg->getTag(), msg->getEpoch(), msg->getNextEpoch(), recvInProgress.size(), commLayer.recvQueue.getSize());

                 // and enqueue ONE FOC message for this tag to be handled by callback.  (reduction from commSize to 1)
                 if (!commLayer.recvQueue.waitAndPush(ptr) ) {
                   throw bliss::io::IOException("C ERROR: recvQueue is not accepting new receivedMessage due to disablePush");
                 }
//                 DEBUGF("R %d <- %d FINISHED RECV CONTROL for tag %d, epoch %lu, next epoch %lu, recvInProgress size %lu, recvQueue size %ld",
//                		 commLayer.commRank, msg->getRank(), msg->getTag(), msg->getEpoch(), msg->getNextEpoch(), recvInProgress.size(), commLayer.recvQueue.getSize());

               } else if (v < 0) {
                 // count decremented to negative.  should not happen.  throw exception.
            	   delete ptr;
                 std::stringstream ss;
                 ss << "R ERROR: number of remaining receivers for tag " << msg->getTag() << " epoch " << msg->getEpoch() << " next epoch " << msg->getNextEpoch() << " is now NEGATIVE";
                 throw bliss::io::IOException(ss.str());
               } else {
            	   delete ptr;
               }

               worked = true;

             } else {  // data message
               DEBUGF("R %d <- %d FINISHED RECV DATA for tag %d epoch %lu, recvInProgress size %lu, recvQueue size %ld",
            		   commRank, ptr->getRank(), ptr->getTag(), ptr->getEpoch(), recvInProgress.size(), commLayer.recvQueue.getSize());

               //==== Data message.  add the received messages into the recvQueue
               if (!commLayer.recvQueue.waitAndPush(ptr) ) {
                 throw bliss::io::IOException("C ERROR: recvQueue is not accepting new receivedMessage due to disablePush");
               }
             }
           }
           else
           {  //== mpi message is not finished yet.

             // the recv request is not done, so stop and wait for the next cycle.
             worked = false;
           }
         } // end while loop for processing recvInProgress.

         //== see if we are done with all pending receives and all Application Termination FOC messages are received.
         if (commLayer.activeEpochCount.load(std::memory_order_acquire) == 0) {
           // if completely done, epoch=CONTROL_TAGGED_EPOCH would have been erased from recvRemaining, after recvInProgress is empty.
           DEBUGF("ACTIVE EPOCH COUNT == 0.  receive done");

           recvDone = true;
           commLayer.recvQueue.disablePush();
         }

         return worked;
       }

   };




    /**
     * @class CallbackThread
     * @brief thread to handle received data using callbacks that are respond to specific message type (MPI tag)
     * @details     This class dequeues from the receive queue, and calls the callback function to process the data.
     */
    class CallbackThread
    {

      public:
       /// constructor
        CallbackThread(CommunicationLayer& _comm_layer) :
          commLayer(_comm_layer), callbackFunctions(),
          td(), currRecvEpoch(0) {

        };
        /// disabled default constructor
        DELETED_FUNC_DECL(CallbackThread());

        /// default destructor
        virtual ~CallbackThread() {
          callbackFunctions.clear();
        }

        /// starts running th thread.
        void start(uint64_t epoch) {
          currRecvEpoch = epoch;
          if (!td.joinable()) td = std::move(std::thread(&bliss::io::CommunicationLayer<ThreadLocal>::CallbackThread::run, this));
        }

        /// finishes the thread and join if possible.
        void finish() {
          if (td.joinable()) {
            td.join();
            DEBUGF("B CallbackThread finished, thread joined.");
          }
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
        void addReceiveCallback(int tag, std::function<void(uint8_t*, std::size_t, int, uint64_t)> callbackFunction, std::unique_lock<std::mutex>&) throw (std::invalid_argument)
        {

          // check for valid arguments
          assert(tag != CONTROL_TAG && tag >= 0);

          // one thread calls this at a time.

          if (callbackFunctions.count(tag) > 0) {
            throw std::invalid_argument("M callback function already registered for given tag");
          }

          // add the callback function to a lookup table
          callbackFunctions[tag] = callbackFunction;
        }

      protected:
        /// parent (containing) class for this thread.
        CommunicationLayer& commLayer;

        /// Registry of callback functions, mapped to by the associated tags
        std::unordered_map<int, std::function<void(uint8_t*, std::size_t, int, uint64_t)> > callbackFunctions;

        /// std::thread object for the dedicated callback-handler thread
        std::thread td;

        /// current epoch of the callback thread - updated on control message receipt.
        uint64_t currRecvEpoch;

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
          DEBUGF("B THREAD STARTED:  callback-thread on %d", commLayer.commRank);

          // while there are messages to process, or we are still receiving
          while ( commLayer.recvQueue.canPop() )  // same as !recvDone || !recvQueue.isEmpty()
          {
            // get next element from the queue, wait if none is available.
            // waitAndPop will exit out of wait when termination flag is set on the recvQueue


            bool suc = false;
            MPIMessage* ptr = nullptr;
            std::tie(suc, ptr) = commLayer.recvQueue.tryPop();
            if (!suc) {
              // no valid result, we must be done with receiving
              //DEBUGF("B R %d recvQueue has nothing. size = %lu", commLayer.commRank, commLayer.recvQueue.getSize());

              //assert(!commLayer.recvQueue.canPush());
              assert(ptr == nullptr);
              continue;
            }
            DEBUGF("B R %d received. src %d, tag %d, epoch %lu queue size = %lu", commLayer.commRank, ptr->getRank(), ptr->getTag(), ptr->getEpoch(), commLayer.recvQueue.getSize());


            if (!ptr) {
            	ERRORF("B R %d received a message with no payload.", commLayer.commRank);

            	continue;
            }

            // if the message is a control message (enqueued ONE into RecvQueue only after all control messages
            //  for a tag have been received.)
            if (ptr->getTag() == CONTROL_TAG) {   // Control message

              ControlMessage* msg = dynamic_cast<ControlMessage*>(ptr);
              DEBUGF("B R %d tag %d epoch %lu control message, nextEpoch %lu",
                       commLayer.commRank, msg->getTag(), msg->getEpoch(), msg->getNextEpoch());

              // here only when all ctrlMsg for the tagged_epoch are received and the flush/finish is complete.
              assert(currRecvEpoch == msg->getEpoch());

              if (currRecvEpoch < msg->getNextEpoch()) {
	              DEBUGF("B R %d tag %d receiver epoch %lu, epoch %lu next epoch %lu notified all.", commLayer.commRank, msg->getTag(), currRecvEpoch, msg->getEpoch(), msg->getNextEpoch());

				  // move to the next epoch.
				  currRecvEpoch = msg->getNextEpoch();

	              // now release the epoch to unblock the corresponding waitForControlMessage.
	              commLayer.epochProperties.releaseEpoch(msg->getEpoch(), std::ref(commLayer.activeEpochCount));
              }

              delete ptr;

            } else {
              // data message.  process it.  msg owns data now.
              SendDataElementType* msg = dynamic_cast<SendDataElementType*>(ptr);
              DEBUGF("B R %d tag %d epoch %lu data message", commLayer.commRank, msg->getTag(), msg->getEpoch());

				if (msg->getEpoch() <= currRecvEpoch) {  // current epoch message, so process it.

				  DEBUGF("B R %d <- %d, tag %d epoch %lu, currRecvEpoch = %lu data message processing.",
												  commLayer.commRank, msg->getRank(), msg->getTag(), msg->getEpoch(), currRecvEpoch);
				  (callbackFunctions.at(msg->getTag()))(msg->getData(), msg->getDataSize(), msg->getRank(), msg->getEpoch());

				  // finished processing.  release back into pool
				  commLayer.releaseBuffer(msg->getBuffer());
				  delete ptr;

				} else {  // else a future epoch message.  put it back.
				  if (! commLayer.recvQueue.waitAndPush(ptr)) {
					  ERRORF("B R %d <- %d, tag %d epoch %lu, currRecvEpoch = %lu data message push back to queue failed.",
							  commLayer.commRank, msg->getRank(), msg->getTag(), msg->getEpoch(), currRecvEpoch);
					  commLayer.releaseBuffer(msg->getBuffer());
					  delete ptr;
					  std::logic_error("failed to push data message back into recv queue");
				  } else
					  DEBUGF("B R %d <- %d, tag %d epoch %lu data message pushed back to queue.  currRecvEpoch = %lu",
						 commLayer.commRank, msg->getRank(), msg->getTag(), msg->getEpoch(), currRecvEpoch);

				}
              // using unique_ptr in control message, don't need to delete msg data.
            }
            // clean up message data
          }
          DEBUGF("B THREAD FINISHED: callback-thread on %d", commLayer.commRank);
        }
    };


public:
  /**
   * @brief Constructor for the CommLayer.
   * @details have to have the comm_size parameter since we can't get it during member variable initialization and recvQueue needs it.
   *
   * @param communicator    The MPI Communicator used in this communicator
   *                        layer.
   * @param comm_size       The size of the MPI Communicator.
   */
  CommunicationLayer (const MPI_Comm& communicator, const int comm_size, const int num_src_threads)
    : sendQueue(), recvQueue(),
      messageBuffers(1, nullptr), epochProperties(comm_size), activeEpochCount(0),
      commThread(communicator, *this), callbackThread(*this),
      numSrcThreads(num_src_threads), initialized(false), currEpoch(CONTROL_TAG_EPOCH), pool()
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
	  if (initialized.load(std::memory_order_relaxed))
		  finishCommunication();

    // final cleanup
    if (!messageBuffers.empty()) {
      for (int i = messageBuffers.size() - 1; i >= 0; --i) {
        if( messageBuffers[i]) {
          ERRORF("tagInfo %d was not cleared prior to finalize\n", i);
          delete messageBuffers[i];
          messageBuffers[i] = nullptr;
        }
      }
      messageBuffers.clear();
    }
  }

  bool isInitialized() {
    return initialized.load(std::memory_order_relaxed);
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
   *   internally, if a buffer is full, it will be enqueued for asynchronous MPI send.
   *
   *   This function is THREAD SAFE
   *
   * @param data        A pointer to the data to be sent.
   * @param count       The number of elements of the message.
   * @param dst_rank    The rank of the destination processor.
   * @param tag         The tag (type) of the message.
   */
  template <typename T>
  void sendMessage(const T* data, const std::size_t count, const int dst_rank, const int tag)
  {
    //== can't have a data message that has a CONTROL_TAG.
    assert(tag != CONTROL_TAG && tag >= 0);
    if (!initialized.load(std::memory_order_relaxed)) throw std::logic_error("calling sendMessage before initialization");

    //check to see if sendQueue is still accepting (if not, then comm is finished)
    if (!sendQueue.canPush()) {
      throw std::logic_error("W sendMessage called after finishCommunication");
    }

    int tid = omp_get_thread_num();

    MessageBuffersType* tptr = getBuffers(tag);

    //== try to append the new data - repeat until successful.
    // along the way, if a full buffer's id is returned, queue it for sendQueue.
    bool suc = false;
    BufferPtrType ptr = nullptr;

    T* data_to_send = const_cast<T*>(data);
    T* data_remain = nullptr;
    uint32_t count_to_send = count;
    uint32_t count_remain = 0;

    int i = 0;
    do {
      // try append.  append fails if there is no buffer or no room

      std::tie(suc, ptr) = tptr->append(data_to_send, count_to_send, data_remain, count_remain, dst_rank, tid);
      data_to_send = data_remain;
      count_to_send = count_remain;

      ++i;

      if (ptr) {        // have a full buffer.  implies !result.first
        // verify that the buffer is not empty (may change because of threading)
        size_t bufSize = ptr->getSize();
        size_t cap = ptr->getCapacity();

        if (bufSize <= cap && bufSize > 0) {
          // have a non-empty buffer - put in send queue.
          SendDataElementType* msg = new SendDataElementType(ptr, tag, dst_rank);

          msg->setEpoch(currEpoch.load(std::memory_order_relaxed));
          //INFOF("W R %d T %d full buffer %p, append %s.  set msg epoch to %lu, curr %lu", commRank, tid, ptr, (suc ? "success" : "failed"), msg->getEpoch(), currEpoch.load(std::memory_order_relaxed));

          if (!sendQueue.waitAndPush(msg).first) {
            ERRORF("W R %d T %d ERROR: sendQueue is not accepting new SendQueueElementType due to disablePush", commRank, tid);
            delete msg;
            releaseBuffer(ptr);
            throw bliss::io::IOException("W ERROR: sendQueue is not accepting new SendQueueElementType due to disablePush");
          } // else successfully pushed into sendQueue.
        } else if (bufSize > cap) {
        	ERRORF("W R %d T %d buffer not finished writing.  should not get here.  Buffer should not have returned from append", commRank, tid);
        } else {// else empty buffer. shared_ptr handles clean up.
        	ERRORF("W R %d T %d empty buffer.  should not get here.  Buffer should have returned a nullptr and commlayer checked earlier.", commRank, tid);
        }

      } // else don't have a full buffer.

      // repeat until success;

      if (i % 1000 == 0) {  // majority is small.
      //if (tag == 15)
        WARNINGF("W R %d T %d: insert into %p took %d iterations: data %p, size %u/%lu, target %d, tag %d, workers %d", commRank, tid, ptr, i, data, count_remain, count, dst_rank, tag, omp_get_num_threads());
      }
    } while (!suc);

    if (i > 200) {  // majority is small.
      WARNINGF("W R %d T %d: NOTICE: insert took %d iterations: data %p, size %u/%lu, target %d, tag %d, workers %d", commRank, tid, i, data, count_remain, count, dst_rank, tag, omp_get_num_threads());
    }


  }


  /**
   * @brief Registers a callback function for the given message tag.  single threaded.
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
  void addReceiveCallback(int tag, std::function<void(uint8_t*, std::size_t, int, uint64_t)> callbackFunction) throw (std::invalid_argument)
  {
    assert(omp_get_num_threads() == 1);
    // check for valid arguments
    assert(tag != CONTROL_TAG && tag >= 0);

    // single thread calls this at a time.
    std::unique_lock<std::mutex> lock(mutex);

    //== if there isn't already a tag listed, add the MessageTypeInfo (and buffer) for that tag.
    if (static_cast<int>(messageBuffers.size()) <= tag) {
      // called by a single thread. no need to lock
      messageBuffers.resize(2 * tag, nullptr);
    }
    messageBuffers[tag] = new MessageBuffersType(pool, commSize, numSrcThreads);

    // now add the callback to callback thread.
    callbackThread.addReceiveCallback(tag, callbackFunction, lock);
  }


  /**
   * @brief flushes the message buffers associated with a particular tag.  Blocks until completion
   * @details call by a single thread only. throws exception if another thread is flushing at the same time.
   *        This method sends all buffered messages with a particular tag to all target receivers.
   *        Then it sends a control message to all target receivers and blocks
   *        A process. on receiving control messages from all senders, unblocks the local flush function
   *
   *        The effect is that collective call to flush enforces that all received messages prior
   *        to the last control messages are processed before the flush call is completed.
   *
   *        Each call of this function increments a "epoch" id, which is sent as part of the control
   *        message to all receivers.  Receivers use the epoch id to ensure that one control message
   *        from each sender has been received.
   * @note okay for multiple threads to call this
   *
   * @param tag
   * @throw bliss::io::IOException
   */
  uint64_t flush(int tag) throw (bliss::io::IOException)
  {
    // can only flush data tag.
    assert(tag != CONTROL_TAG && tag >= 0);
    if (!initialized.load(std::memory_order_relaxed)) throw std::logic_error("calling flush before initialization");

    if (!sendQueue.canPush()) {
      throw std::logic_error("M ERROR flush called but sendQueue is disabled.");
    }

//    MessageInfo& taginfo = getTagInfo(tag);
//
//    //== check to see if the target tag is still accepting
//    if (taginfo.isFinished()) {
//      throw std::invalid_argument("M ERROR tag not registered, or already finished. cannot FLUSH");
//    }

    DEBUGF("M R %d,\t,\t,\t,\tt %d,\tFLUSH BUFFERS", commRank, tag);

    // flush all data buffers (put them into send queue)
    flushBuffers(tag);

    // send the control message - generates an unique epoch for the tag as well.  (also into the send queue, to maintain ordering)
    // wait till all control messages have been received (not using MPI_Barrier.  see Rationale.)
    DEBUGF("M R %d,\t,\t,\t,\tt %d\tnew e ,\tSEND CTRL MSGS", commRank, tag);
//    sendControlMessagesAndWait(tag, currEpoch);
    uint64_t oe = startNextEpoch();
    DEBUGF("M R %d,\t,\t,\t,\tt %d\toe %lu,\tRECVED CTRL MSGS", commRank, tag, oe);

    return oe;

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
  uint64_t finish(int tag) throw (bliss::io::IOException)
  {

    assert(tag != CONTROL_TAG && tag >= 0);
    if (!initialized.load(std::memory_order_relaxed)) throw std::logic_error("called finish when commLayer is not initialized");

    if (!sendQueue.canPush()) {
    	throw std::logic_error("called finish when commLayer sendQueue is disabled.");
    }

//    // Mark tag as finished in the metadata.
//    MessageInfo& taginfo = getTagInfo(tag);
//
//    //== check to see if the target tag is still accepting
//    if (taginfo.finish()) {
//      // already finished.  return.
//      return;
//    }

    DEBUGF("M R %d tag %d FINISH", commRank, tag);


    // flush all data buffers (put them into send queue)
    flushBuffers(tag);

    // send the control message - generates an unique epoch for the tag as well.  (also into the send queue, to maintain ordering)
    // wait till all control messages have been received (not using MPI_Barrier.  see Rationale.)

    DEBUGF("M R %d,\t,\t,\t,\tt %d\tnew e ,\tSEND CTRL MSGS for FINISH", commRank, tag);
//    sendControlMessagesAndWait(tag, currEpoch);
    uint64_t oe = startNextEpoch();
    DEBUGF("M R %d,\t,\t,\t,\tt %d\toe %lu,\tRECVED CTRL MSGS for FINISH", commRank, tag, oe);


    std::unique_lock<std::mutex> lock(mutex);    // ensure a single thread doing this.

    // remove this tag from further use.
    deleteTagInfo(tag, lock);

    return oe;
  }


  /**
   * @brief Initializes all communication and starts the communication and
   *        callback threads.
   *
   * @note Must be called before any other functions can be called.
   */
  uint64_t initCommunication()
  {
    assert(omp_get_num_threads() == 1);

    if (initialized.load(std::memory_order_relaxed)) throw std::logic_error("called initCommunication when commLayer is already initialized");


    // one thread calls this
    if(!sendQueue.canPush()) sendQueue.enablePush();
    if(!recvQueue.canPush()) recvQueue.enablePush();

    std::unique_lock<std::mutex> lock(mutex);

    // reserve space
    if (messageBuffers.size() <= CONTROL_TAG) {
      std::unique_lock<std::mutex> lock(mutex);
      messageBuffers.reserve(CONTROL_TAG * 2 + 1);
    }

    // store Control Tag.
    if (messageBuffers[CONTROL_TAG] == nullptr) {
      //== init with control_tag only.  only 1 thread uses this.
      messageBuffers[CONTROL_TAG] = new MessageBuffersType(pool, commSize, 1);  // only 1 thread needs to use this set.
    }

    lock.unlock();

    // initialize the epochs
    uint64_t ne = startFirstEpoch();

    DEBUGF("R %d initCommunication epoch %lu", commRank, ne);

    callbackThread.start(ne);  // uses constructor to pass in the current epoch.
    commThread.start();

    initialized.store(true, std::memory_order_relaxed);

    return ne;
  }


  /**
   * @brief stops accepting messages of all types and flushes all existing messages.  Blocks until completion
   * @details call by a single thread only. throws exception if another thread is flushing at the same time.
   *        This method follows the same logic as "finish", but for ALL message types
   *
   *        At the end, sends an application termination control message.
   *
   * @param tag
   * @throw bliss::io::IOException
   */
  uint64_t finishCommunication() throw (bliss::io::IOException)
  {
    assert(omp_get_num_threads() == 1);

    if (!initialized.load(std::memory_order_relaxed)) throw std::logic_error("called finishCommunication when commLayer is not initialized");

    DEBUGF("M R %d FINISH COMM ", commRank);


    // go through all tags and "finish" them, except for the final control tag.
    // recall that MessageTypeInfos are stored in a vector with index == tag.
    int mx = messageBuffers.size();
    for (int i = 0; i < mx; ++i) {
      if (i != CONTROL_TAG && messageBuffers[i]) finish(i);
    }

    // send the control message for the final control tag.
    // wait till all control messages have been received (not using MPI_Barrier.  see Rationale.)
    DEBUGF("M R %d FINISH COMM.  next epoch %lu, waiting for CONTROL messages to complete", commRank, CONTROL_TAG_EPOCH);
    //sendControlMessagesAndWait(CONTROL_TAG, CONTROL_TAG_EPOCH);  // only 1 control message needs to be sent.
    uint64_t oe = finishLastEpoch();
    DEBUGF("M R %d FINISH COMM epoch %lu, finished wait for CONTROL messages", commRank, oe);


    std::unique_lock<std::mutex> lock(mutex);    // ensure a single thread doing this.

    deleteTagInfo(CONTROL_TAG, lock);  // delete and set to null

    DEBUGF("remaining epochs: %s", epochProperties.toString().c_str());

    lock.unlock();



    assert(activeEpochCount.load(std::memory_order_relaxed) == 0);

    // wait for both comm and callback threads to quit
    commThread.finish();
    callbackThread.finish();


    if(sendQueue.canPush()) sendQueue.disablePush();
    if(recvQueue.canPush()) recvQueue.disablePush();

    initialized.store(false, std::memory_order_relaxed);

    return oe;
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

  /// retrieve the MessageTypeInfo object for specified tag.
  MessageBuffersType* getBuffers(int tag) throw (std::out_of_range) {
    //== check to see if the target tag is still accepting
    // modification to messageBuffers only at finish(), and finishCommunications(), both are after sendMessage
    // and sendControlMessageAndWait.  these should be single threaded, so we don't need a lock here.

    if (messageBuffers[tag] == nullptr) {
      ERRORF("W R %d ERROR: invalid tag: tag has not been registered %d", commRank, tag);
      throw std::out_of_range("W invalid tag: tag has not been registered");
    }
    return messageBuffers.at(tag);
  }

  /// delete the MessageTypeInfo object for specified tag.
  void deleteTagInfo(int tag, std::unique_lock<std::mutex> const &) throw (std::invalid_argument){
    if (messageBuffers[tag] == nullptr) {
      throw std::invalid_argument("W invalid tag: tag has not been registered");
    }
    delete messageBuffers[tag];
    messageBuffers[tag] = nullptr;
  }


  void releaseBuffer(BufferPtrType ptr) {
      if (ptr) {
        ptr->block_and_flush();  // if already blocked and flushed, no op.

        pool.releaseObject(ptr);
      }
  }

  BufferPtrType acquireBuffer() {
      // get a buffer from pool
      BufferPtrType ptr = pool.tryAcquireObject();
      while (!ptr) {
        _mm_pause();
        ptr = pool.tryAcquireObject();
      }

      ptr->clear_and_unblock_writes();

      return ptr;
  }


  /// start a new epoch.  this assumes no previous epoch.  note that callback thread is not notified here.
  uint64_t startFirstEpoch() {
    // must be initialized
    if (initialized.load(std::memory_order_relaxed)) throw std::logic_error("calling startFirstEpoch after initialization");

    uint64_t ne = MPIMessage::nextEpoch.fetch_add(1, std::memory_order_relaxed);

    // get a new epoch and exchange it.
//    uint64_t ce = bliss::io::MPIMessage::nextEpoch.fetch_add(1, std::memory_order_relaxed);
    uint64_t oe = currEpoch.exchange(ne, std::memory_order_acq_rel);
    // must not have previous epoch.
    assert(oe == CONTROL_TAG_EPOCH);
    DEBUGF("M R %d obtained first epoch %lu, next %lu", commRank, oe, ne);

    // finally, register it.
    epochProperties.registerEpoch(ne, std::ref(activeEpochCount));

    return ne;
	}

  /// finish the previous epoch and start a new epoch.  this assumes there is a previous epoch
  uint64_t startNextEpoch() {
    //=== create a new epoch
    // must be initialized
    if (!initialized.load(std::memory_order_relaxed)) throw std::logic_error("calling startNextEpoch before initialization");

    // queue must be active
    if (!sendQueue.canPush())  throw std::logic_error("M ERROR: sendControlMessagesAndWait already called");


    // get a new epoch and exchange it.
//    uint64_t ce = bliss::io::MPIMessage::nextEpoch.fetch_add(1, std::memory_order_relaxed); // new
    uint64_t ne = MPIMessage::nextEpoch.fetch_add(1, std::memory_order_relaxed);

    DEBUGF("M R %d obtain next epoch %lu, next %lu", commRank, currEpoch.load(std::memory_order_relaxed), ne);

    uint64_t oe = currEpoch.exchange(ne, std::memory_order_acq_rel);  // old
    // must have previous epoch.
    assert(oe != CONTROL_TAG_EPOCH);

    DEBUGF("M R %d obtained next epoch %lu, next %lu", commRank, oe, ne);

    // register it.
    epochProperties.registerEpoch(ne, std::ref(activeEpochCount));
	  
    DEBUGF("M R %d registered next epoch %lu, next %lu", commRank, oe, ne);

    //=== THEN finish the last epoch
    // create a control message to send to all MPI processes, and queue it.  also notifies callback thread of epoch change.
    ControlMessage* msg = new ControlMessage(oe, ne, 0);
    if (!sendQueue.waitAndPush(msg).first) {
    	DEBUGF("M R %d inserted new control message for epoch %lu, next epoch %lu", commRank, oe, ne);
    	delete msg;
      throw std::logic_error("M ERROR: sendQueue is not accepting new SendQueueElementType due to disablePush");
    }

    //============= now wait for the old epoch to finish.  the next has been set up.
    DEBUGF("M R %d WAIT for epoch %lu, next %lu, currently flushing.  recvQueue size %lu", commRank, oe, ne, recvQueue.getSize());
    epochProperties.waitForEpochRelease(oe);  // TODO:  deadlock caused by this not returning/.
    DEBUGF("M R %d WAITED for epoch %lu, next %lu, currently flushing", commRank, oe, ne);
    // then clean up.

    //============= DEBUGGING
//    if (recvQueue.getSize() > 0) {
//      WARNINGF("R %d finished waiting for epoch %lu.  recvQueue: %lu, sendQueue %lu.\n", commRank, oe, recvQueue.getSize(), sendQueue.getSize());
//    }

    DEBUGF("M R %d <- -1 received all END message for epoch = %lu, next %lu", commRank, oe, ne);

    return oe;
  }

  /// finishes the last epoch.  this assumes previous epoch is set.
  uint64_t finishLastEpoch() {
    // must be initialized
	if (!initialized.load(std::memory_order_relaxed)) throw std::logic_error("calling finishLastEpoch before initialization");

    // queue must be active
    if (!sendQueue.canPush())  throw std::logic_error("M ERROR: sendControlMessagesAndWait already called");
    
    // get a new epoch and exchange it.
    uint64_t oe = currEpoch.exchange(CONTROL_TAG_EPOCH, std::memory_order_relaxed);
    assert(oe != CONTROL_TAG_EPOCH);

    // create a control message to send to all MPI processes, and queue it.  notifies the callback thread of epoch change.
    // note that use of CONTROL_TAG_EPOCH = FFFFFFFF causes all previous messages to flush.
    ControlMessage* msg = new ControlMessage(oe, CONTROL_TAG_EPOCH, 0);
    if (!sendQueue.waitAndPush(msg).first) {
    	delete msg;
      throw std::logic_error("M ERROR: sendQueue is not accepting new SendQueueElementType due to disablePush");
    }

    // last message sent, stop the send queue.
    sendQueue.disablePush();

    //============= now wait for the current epoch to finish.  the next has been set up.
    DEBUGF("M R %d WAIT for epoch %lu, currently flushing", commRank, oe);
    epochProperties.waitForEpochRelease(oe);
    DEBUGF("M R %d WAITED for epoch %lu, currently flushing", commRank, oe);


    //============= DEBUGGING
//    if (recvQueue.getSize() > 0) {
//        WARNINGF("R %d finished waiting for epoch %lu.  recvQueue: %lu, sendQueue %lu.\n", commRank, oe, recvQueue.getSize(), sendQueue.getSize());
//    }

    DEBUGF("M R %d <- -1 received all END message for epoch = %lu", commRank, oe);

    return oe;
  }



  /**
   * @brief Flushes all data buffers for message tag `tag`.
   * @details All pending/buffered messages of type `tag` will be send out.
   *    enqueues data buffers to be sent, if buffer is not empty.
   *
   * @param tag           The message tag, whose messages are flushed.
   * @param next_epoch    The next epoch - new bufers that are swapped in will have this value.  to just flush, pass in the same value
   * @throw bliss::io::IOException
   */
  uint64_t flushBuffers(int tag) throw (bliss::io::IOException)
  {
	if (!initialized.load(std::memory_order_relaxed)) throw std::logic_error("calling flushBuffers before initialization");

    if (!sendQueue.canPush()) {
      throw std::logic_error("W cannot flush buffer: cannot push to sendQueue.");
    }

    MessageBuffersType* tptr = getBuffers(tag);
    uint64_t ce = currEpoch.load(std::memory_order_acquire);

    int target_rank = commRank;
    for (int i = 0; i < commSize; ++i) {
      // flush buffers in a circular fashion, starting with the next neighbor
      target_rank = (target_rank + 1) % commSize;

      auto ptrs = tptr->flushBufferForRank(target_rank);  // returns vector of buffer pointers, one per thread
      // flush/send all remaining non-empty buffers
      for (auto bufPtr : ptrs) {
        if (bufPtr) {
          if (!(bufPtr->isEmpty())) {
            // put all buffers onto sendQueue
            SendDataElementType* msg = new SendDataElementType(bufPtr, tag, target_rank);

            msg->setEpoch(ce);
            //DEBUGF("W Flush tag %d.  epoch set to %lu actual %lu", tag, ce, msg->getEpoch());

            if (!sendQueue.waitAndPush(msg).first) {
              releaseBuffer(bufPtr);
              delete msg;

              throw bliss::io::IOException("M ERROR: sendQueue is not accepting new SendQueueElementType");
            }
          } else {
            releaseBuffer(bufPtr);
          }

        }

      }
    }
    return ce;
  }



protected:
/******************************************************************************
 *                           Protected member fields                            *
 ******************************************************************************/

  /// Mutex for adding/removing from local arrays and maps.
  mutable std::mutex mutex;

  /// message queue between sendMessage calling threads (src) and comm-thread
  /// (sink).  multiple producer, single consumer thread-safe queue.
  bliss::concurrent::ThreadSafeQueue< MPIMessage* , bliss::concurrent::LockType::MUTEX> sendQueue;

  /// message queue between comm thread (src) and callback thread (sink)
  /// single producer potentially multiple consumer queue.
  bliss::concurrent::ThreadSafeQueue< MPIMessage* , bliss::concurrent::LockType::MUTEX> recvQueue;


  // Message buffers per tag (message type) is stored in MessageTypeInfo.

  /// DEFINE message type info class based on MessageBuffersType
  //using MessageInfo = bliss::io::MessageTypeInfo<MessageBuffersType>;
  using EpochInfo = bliss::io::EpochInfo;

  /**
   *  map between tag (message type) and MessageTypeInfo<MessageBuffersType>
   *  which contains active epochs, MessageBuffers Pointer, and mutexes.
   *
   *  using vector because unordered map requires locking.
   *  vector is allocated so that the tag number directly map to vector index.
   */
  //std::vector< MessageInfo* > messageBuffers;
  // TODO:  replace with vector of MessageBuffers
  std::vector< MessageBuffersType* > messageBuffers;
  bliss::io::EpochInfo epochProperties;

  /// total count of active epochs across all message types
  std::atomic<int> activeEpochCount;
  // per message type count is stored in MessageTypeInfo

  /// communication thread for send and recv messages.
  CommThread<BufferPtrType> commThread;

  /// call back thread to process the received objects.
  CallbackThread callbackThread;

  /// The MPI Communicator size
  mutable int commSize;

  /// The MPI Communicator rank
  mutable int commRank;

  /// the max number of source threads that Comm Layer can serve/use
  mutable int numSrcThreads;

  /// flag indicating that the comm layer has been initialized
  std::atomic<bool> initialized;

  /// current epoch.
  std::atomic<uint64_t> currEpoch;   // TODO: does it need to be atomic?

  /// object pool for buffers
  BufferPoolType pool;
};


} // namespace io
} // namespace bliss

#endif // BLISS_COMMUNICATION_LAYER_HPP
