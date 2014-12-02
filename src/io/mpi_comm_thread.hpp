/**
 * @file    mpi_comm_thread.hpp
 * @ingroup
 * @author  tpan
 * @brief
 * @details
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef MPI_COMM_THREAD_HPP_
#define MPI_COMM_THREAD_HPP_

#include <memory>
#include <mpi.h>
#include <xmmintrin.h>

#include <concurrent/threadsafe_queue.hpp>
#include <io/message_types.hpp>

namespace bliss
{
  namespace io
  {


    // ==== send comm thread class.  for organization.
    // TODO: group the code for comm thread
    template <typename BufferPtrType>
    class CommThread
    {
      protected:


        using TaggedEpoch = typename MessageTypeInfo::TaggedEpoch;
        using SendDataElementType = typename bliss::io::DataMessageToSend<BufferPtrType>;


        CommThread() = delete;


        /// The MPI Communicator object for this communication layer
        MPI_Comm comm;

        /// Atomic status variable: set to true if `finishing` and if all sends are
        /// completed.
        bool sendDone;

        /// Atomic status variable: set to true if `finishing` and if all receives
        /// are completed.
        bool recvDone;


        /// Queue of pending MPI send operations
        std::deque<std::pair<MPI_Request, std::unique_ptr<MPIMessage> > > sendInProgress;

        // pending mpi messages.  also holds local messages (send to self) between send thread and recv thread.
        //   (2 producers, 1 consumer) thread.
        /// Queue of pending MPI receive operations
        //bliss::concurrent::ThreadSafeQueue<std::pair<MPI_Request, std::unique_ptr<MPIMessage> > > recvInProgress;
        std::deque<std::pair<MPI_Request, std::unique_ptr<MPIMessage> > > recvInProgress;

        std::thread td;
        int commRank;
        int commSize;



      public:


        // Outbound message data structure, Multiple-Producer-Single-Consumer queue
        // consumed by the internal MPI-comm-thread, produced from worker threads
        /// message queue between sendMessage calling threads (src) and comm-thread
        /// (sink)
        bliss::concurrent::ThreadSafeQueue<std::unique_ptr<MPIMessage> > sendQueue;

        // Inbound message data structure, Single-Producer-Multiple-Consumer queue
        // produced by the internal MPI-comm-thread, consumed by one or more callback threads
        /// message queue between comm thread (src) and callback thread (sink)
        bliss::concurrent::ThreadSafeQueue<std::unique_ptr<MPIMessage> > recvQueue;



        CommThread(const MPI_Comm& communicator) :
          comm(communicator),
          sendDone(false),
          recvDone(false),
          sendQueue(),
          recvQueue(),
          sendInProgress(),
          recvInProgress(),
          td()
        {
        };

        virtual ~CommThread() {}

        void start() {
          MPI_Comm_size(comm, &commSize);
          MPI_Comm_rank(comm, &commRank);

          if (!td.joinable()) {

            //===== initialize the mpi buffer for control messaging, for MPI_Bsend
            int mpiBufSize = 0;
            // arbitrary, 4X commSize, each has 2 ints, one for tag, one for epoch id,.
            MPI_Pack_size(commSize * 4, MPI_LONG, comm, &mpiBufSize);
            mpiBufSize += commSize * 4 * MPI_BSEND_OVERHEAD;
            char* mpiBuf = (char*)malloc(mpiBufSize);
            MPI_Buffer_attach(mpiBuf, mpiBufSize);

            td = std::move(std::thread(&bliss::io::CommThread::run, this));
          } // else already started the thread.
        }

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

          // TODO: use thread lock mutex and cond var to help reduce load.

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

          //DEBUGF("rank %d sendDone? %s, finishing %s sendqueue %ld sendInProgress %ld", commRank, (sendDone ? "true":"false"), (finishing.load() ? "true": "false"), sendQueue.getSize(), sendInProgress.size());



          // try to get the send element
          auto el = std::move(sendQueue.tryPop());  // type: std::pair<bool, std::unique_ptr<MPIMessage>>
          if (!el.first) {  // no valid entry in the queue
            return false;
          }
          // else there is a valid entry from sendQueue
          assert(el.second);

          bool worked = false;

          // second is a pointer to ControlMessage
          if (el.second->tag == MessageTypeInfo::CONTROL_TAG) {
            // if control message,

            // ControlMessage ownership is maintained by el.
            ControlMessage* se = dynamic_cast<ControlMessage*>(el.second.get());  // MPIMessage*
            TaggedEpoch te = se->tagged_epoch;

            DEBUGF("C SEND %d -> %d, termination signal for tag %d epoch %ld", commRank, se->rank, MessageTypeInfo::getTag(te), te & 0x00000000FFFFFFFF);

            // termination message for this tag and destination
            if (se->rank == commRank) {
              // local, directly handle by creating an byte array and directly
              // insert into the recvInProgress (needed by the receiver decrement
              // logic in "finishReceive", so can't put into receive queue)
              MPI_Request req = MPI_REQUEST_NULL;
//              uint8_t *array = new uint8_t[sizeof(TaggedEpoch)];
//              memcpy(array, &(te), sizeof(TaggedEpoch));



              // now move ControlMessage ownership to recvInProgress
              auto pair = std::make_pair(std::move(req),
                                         std::move(el.second));
              assert(el.second == nullptr);
              assert(pair.second != nullptr);

              recvInProgress.push_back(std::move(pair));
              assert(pair.second.get() == nullptr);
              assert(recvInProgress.back().second.get() != nullptr);


              DEBUGF("C rank %d tag %d epoch %ld local send control message to rank %d", commRank, MessageTypeInfo::getTag(te), te & 0x00000000FFFFFFFF, commRank);
            } else {


//               remote control message.  send a terminating message with tag CONTROL_TAG
//               payload is the message type and epoch.  Message Ordering is guaranteed by MPI.
//               Bsend to allow use of preallocated buffer, so don't have to clean up allocated messages.
              MPI_Bsend(&(te), 1, MPI_LONG, se->rank, MessageTypeInfo::CONTROL_TAG, comm);


              DEBUGF("C rank %d tag %d epoch %ld remote send control message to rank %d", commRank, MessageTypeInfo::getTag(te), te & 0x00000000FFFFFFFF, se->rank);
            }

            worked = true;
          } else {
            // data message.

            SendDataElementType* se = dynamic_cast<SendDataElementType*>(el.second.get());  // MPIMessage*

            // get message data and it's size
            uint8_t* data = se->ptr->operator uint8_t*();  // BufferPtrType then getData
            auto count = se->ptr->getSize();

            if (count > 0) {

              if (se->rank == commRank) {
                // local, directly handle by creating an output object and directly
                // insert into the recv queue (thread safe)

                // since the send buffer memory is managed by Buffer, and the receiveMessage expects to manage memory as well, need to copy.
                std::unique_ptr<uint8_t[]> array(new uint8_t[count]);
                memcpy(array.get(), data, count);

                // pointer here to allow dynamic cast.
                std::unique_ptr<DataMessageReceived> dmr(new DataMessageReceived(std::move(array), count, se->tag, commRank));


                // put into recvInProgress queue instead of recvQueue to minimize out of order (early) receive.
                MPI_Request req = MPI_REQUEST_NULL;

                auto pair = std::make_pair(std::move(req), std::move(dmr));
                assert(dmr == nullptr);
                assert(pair.second != nullptr);

                recvInProgress.push_back(std::move(pair));
                assert(pair.second.get() == nullptr);

                assert(recvInProgress.back().second.get() != nullptr);

                // finished inserting directly to local RecvInProgress.  release the buffer

                // TODO:  don't release the buffer for local transmit.  do it after callback thread is done with it.
                commLayer.ctrlMsgProperties.at(se->tag).getBuffers()->releaseBuffer(std::move(se->ptr));
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
//              WARNINGF("C WARNING: NOT SEND %d -> %d: 0 byte message for tag %d, ptr = %d.", commRank, se->rank, se->tag, se->bufferId);
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

          if (commSize > 1)
            MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, comm, &hasMessage, &status);

          // if have message to receive,
          if (hasMessage > 0) {
            // get some message details
            int src = status.MPI_SOURCE;
            int tag = status.MPI_TAG;
            int received_count;
            MPI_Request req;


            // receive the message data bytes into a vector of bytes
            if (tag == MessageTypeInfo::CONTROL_TAG) {  // control messages
              MPI_Get_count(&status, MPI_LONG, &received_count);
              assert(received_count == 1);

              std::unique_ptr<ControlMessage> msg(new ControlMessage(MessageTypeInfo::CONTROL_TAGGED_EPOCH, src));

              // receive data or FOC messages.  use finishReceive to handle the
              // FOC messages and count decrement

              MPI_Irecv(&(msg->tagged_epoch), received_count, MPI_LONG, src, tag, comm, &req);

              DEBUGF("R RECV %d -> %d, termination signal", src, commRank);

              // insert into the received InProgress queue.
              auto pair = std::make_pair(std::move(req), std::move(msg));
              assert(msg.get() == nullptr);
              assert(pair.second.get() != nullptr);

              recvInProgress.push_back(std::move(pair));
              assert(pair.second.get() == nullptr);
              assert(recvInProgress.back().second.get() != nullptr);


            } else {  // data messages
              MPI_Get_count(&status, MPI_UNSIGNED_CHAR, &received_count);

              // have to receive even if count is 0.
              uint8_t* msg_data = (received_count == 0) ? nullptr : new uint8_t[received_count];
              std::unique_ptr<DataMessageReceived> msg(new DataMessageReceived(msg_data, received_count, tag, src));
              // owner of msg-data is now msg.


              MPI_Irecv(msg_data, received_count, MPI_UNSIGNED_CHAR, src, tag, comm, &req);

              auto pair = std::make_pair(std::move(req), std::move(msg));
              assert(msg.get() == nullptr);
              assert(pair.second.get() != nullptr);

              // insert into the received InProgress queue.
              recvInProgress.push_back(std::move(pair));

              assert(pair.second.get() == nullptr);
              assert(recvInProgress.back().second.get() != nullptr);

            }
            return true;
          }
          else return false;
        }


        /**
         * @brief Finishes pending MPI_Send requests.  stops at first unfinished request, left unfinished requests finish in next iteration's call.
         */
        bool finishSends()
        {
          // no more messages to send or being sent, so done.
          if (sendDone) return false;

          int mpi_finished = 0;

          bool worked = false;
          // while there is some in progress MPI_send requests,
          while(!sendInProgress.empty())
          {
            // get the first request to check - ONLY 1 THREAD CHECKING sendInProgress.
            assert(sendInProgress.front().second != nullptr);


            // check if MPI request has completed.
            if (sendInProgress.front().first == MPI_REQUEST_NULL) {
              // this is a control message since it has no "request".
              mpi_finished = 1;
            } else {
              // has request.  check for finished
              MPI_Test(&(sendInProgress.front().first), &mpi_finished, MPI_STATUS_IGNORE);
            }

            // if finished, then clean up.
            if (mpi_finished)
            {
              // if there is a buffer, then release it.
              if (sendInProgress.front().second->tag != MessageTypeInfo::CONTROL_TAG) {
                SendDataElementType* msg = dynamic_cast<SendDataElementType*>(sendInProgress.front().second.get());
                if (msg->ptr) {
                  // cleanup, i.e., release the buffer back into the pool
                  commLayer.ctrlMsgProperties.at(msg->tag).getBuffers()->releaseBuffer(std::move(msg->ptr));
                }
              }  // control tag, and mpi send finished, so nothing to clean up.


              // remove the entry from the in-progress queue.  object destruction will clean up data associated with tthe pointer.
              sendInProgress.pop_front();

              worked = true;
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
          if (!sendQueue.canPop() && sendInProgress.empty()) {
            // not using sendAccept as check - sendAccept is cleared BEFORE the last
            // tag end messages starts to be sent, and definitely before the app end
            // message starts to be sent so sendDone may be set to true too early,
            // and elements for the last tag, or tag 0 (termination) may be added to
            // the sendQueue or sendInProfess queue after the empty check here.
            sendDone = true;
            //
            //sendQueue.disablePush();
          }
          return worked;
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
        bool finishReceives() throw (bliss::io::IOException)
        {
          // no more to receive or pending receive requests, so don.
          if (recvDone) return false;

          int mpi_finished = 0;

          bool worked = false;

          // if there are pending receive requests
          while(!recvInProgress.empty())
          {

            assert(recvInProgress.front().second.get() != nullptr);


            // check if MPI request has completed.
            if (recvInProgress.front().first == MPI_REQUEST_NULL) {  // locla message
              // this is a control message since it has no "request".
              mpi_finished = 1;
            } else {
              // has request.  check for finished
              MPI_Test(&(recvInProgress.front().first), &mpi_finished, MPI_STATUS_IGNORE);
            }


            // if finished, then clean up.
            if (mpi_finished)
            {
              // finished request, so enqueue the message for processing.

              // get the request
              auto front = std::move(recvInProgress.front());
              // remove pending receive element
              recvInProgress.pop_front();


              // if it's a control message, then need to count down
              if (front.second.get()->tag == MessageTypeInfo::CONTROL_TAG)  {

                // get the control message
                ControlMessage* msg = dynamic_cast<ControlMessage*>(front.second.get());
                TaggedEpoch te = msg->tagged_epoch;

//                DEBUGF("C rank %d received control message for tag %d epoch %ld.  recvRemaining has %ld tags", commRank, tag, epoch, recvRemaining.size());

                //==== handle case where local process is not yet listening for tagged epoch.
                // if this process is not yet waiting for this epoch, that means it has not flushed or finished that tag,
                // then we need to push this message back into recvInProgress.
                // need to construct a new message with a NULL request.
                // then skip over the remaining of this branch and leave the message in recvInProgress.
                if (commLayer.ctrlMsgRemainingForEpoch.count(te) == 0) {
                  //if (iters %100 == 0) WARNINGF("WARN: RANK %d not yet waiting for control message tag %d epoch %ld.  cycle through again.", commRank, te.tag, te.epoch);

                  front.first = MPI_REQUEST_NULL;
                  assert(front.second.get() != nullptr);
                  recvInProgress.push_back(std::move(front));
                  assert(front.second.get() == nullptr);
                  assert(recvInProgress.back().second.get() != nullptr);

                  break;  // gives the worker threads a chance to call sendControlMessagesAndWait

                  // this delays the processing, which hopefully does not create deadlock.
                  // Better for worker thread to wait via  sendControlMessagesAndWait, instead
                  // of comm thread waiting for the taggedepoch to become present.
                }


                //== now decrement the count of control messages for a epoch
//                DEBUGF("C RECV PRE rank %d receiving END signal tag %d epoch %ld from %d, recvRemaining size: %ld",
//                     commRank, te.tag, te.epoch, front.second->rank, commLayer.ctrlMsgRemainingForEpoch.size());
                auto v = --(commLayer.ctrlMsgRemainingForEpoch.at(te).val);
//                DEBUGF("C RECV rank %d receiving END signal tag %d epoch %ld from %d, num senders remaining is %d",
//                     commRank, te.tag, te.epoch, front.second->rank, commLayer.ctrlMsgRemainingForEpoch.at(te));


                //== if after decrement the count is 0 for the ep,
                if (v == 0) {

                  //DEBUGF("ALL END received for tag %d, pushing to recv queue", front.second.tag);


                  // and enqueue ONE FOC message for this tag to be handled by callback.  (reduction from commSize to 1)
                  if (!recvQueue.waitAndPush(std::move(front.second)).first) {
                    throw bliss::io::IOException("C ERROR: recvQueue is not accepting new receivedMessage due to disablePush");
                  }

                } else if (v < 0) {
                  // count decremented to negative.  should not happen.  throw exception.

                  std::stringstream ss;
                  ss << "R ERROR: number of remaining receivers for tag " << MessageTypeInfo::getTag(te) << " epoch " << (te & 0x00000000FFFFFFFF) << " is now NEGATIVE";
                  throw bliss::io::IOException(ss.str());
                }

                worked = true;


              } else {  // data message


                //==== Data message.  add the received messages into the recvQueue
                if (!recvQueue.waitAndPush(std::move(front.second)).first) {
                  throw bliss::io::IOException("C ERROR: recvQueue is not accepting new receivedMessage due to disablePush");
                }
              }
            }
            else
            {  // mpi message is not finished yet.

              // the recv request is not done, so stop and wait for the next cycle.
              break;
            }
          } // end while loop for processing recvInProgress.

          //== see if we are done with all pending receives and all Application Termination FOC messages are received.
          if (ctrlMsgRemainingForEpoch.empty()) { //&& !recvInProgress.canPop()) {
            // if completely done, epoch=CONTROL_TAGGED_EPOCH would have been erased from recvRemaining.
            recvDone = true;
            recvQueue.disablePush();
          }

          return worked;
        }

    };


  } /* namespace io */
} /* namespace bliss */

#endif /* MPI_COMM_THREAD_HPP_ */
