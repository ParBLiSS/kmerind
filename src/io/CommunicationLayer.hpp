/**
 * @file    CommunicationLayer.hpp
 * @ingroup index
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @author  Tony Pan <tpan@gatech.edu>
 * @brief
 *
 * Copyright (c) TODO
 *
 * TODO add Licence
 *
 *
 * TODO:  stopping the communication.
 *    need 1. mapping between message types.  assumption, each incoming message type corresponds to one or more outgoing message types.
 *          2. method to mark a message type producer as being done. - single thread.
 *              also notify remote receiver with done (for that message type/tag)
 *          3. receiver decrements sender count down to 0, at which point, mark the corresponding response send message as done.
 *          4. when all send types are done, and all receive types are done, the comm thread terminates.
 *          5. when all receive type are done, the recv thread terminates.
 *
 */

#ifndef BLISS_COMMUNICATION_LAYER_HPP
#define BLISS_COMMUNICATION_LAYER_HPP

#include <mpi.h>
#include <omp.h>

// C stdlib includes
#include <assert.h>
#include <cstdio>  // fflush

// STL includes
#include <vector>
#include <queue>
#include <mutex>
#include <unordered_map>
#include <unordered_set>
#include <thread>

#include <utility>


// system includes
// TODO: make this system indepenedent!?
#include <unistd.h> // for usleep()

// BLISS includes
#include <concurrent/threadsafe_queue.hpp>
#include <concurrent/concurrent.hpp>
#include <io/message_buffers.hpp>

typedef bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_SAFE> BuffersType;


struct ReceivedMessage
{
  uint8_t* data;
  std::size_t count;
  int tag;
  int src;

  ReceivedMessage(uint8_t* data, std::size_t count, int tag, int src)
    : data(data), count(count), tag(tag), src(src) {}
  ReceivedMessage() = default;
};

// TODO: rename (either this or the ReceivedMessage) for identical naming scheme
struct SendQueueElement
{
  typename BuffersType::BufferIdType bufferId;
  int tag;
  int dst;

  SendQueueElement(BuffersType::BufferIdType _id, int _tag, int _dst)
    : bufferId(_id), tag(_tag), dst(_dst) {}
  SendQueueElement() = default;
};


class CommunicationLayer
{
public:

  constexpr static int default_tag = 0;

protected:

  /// USED BY comm thread, purely internal and purely single thread.
  // request, data pointer, data size
  std::queue<std::pair<MPI_Request, ReceivedMessage> > recvInProgress;
  // request, data pointer, tag
  std::queue<std::pair<MPI_Request, SendQueueElement> > sendInProgress;

  /// message queue between sendMessage calling threads (src) and comm thread (sink)
  // Outbound message structure, Multiple-Producer-Single-Consumer queue
  // consumed by the internal MPI-comm thread
  bliss::concurrent::ThreadSafeQueue<SendQueueElement> sendQueue;

  /// message queue between comm thread (src) and callback thread (sink)
  // Inbound message structure, Single-Producer-Multiple-Consumer queue
  // produced by the internal MPI-comm thread
  bliss::concurrent::ThreadSafeQueue<ReceivedMessage> recvQueue;

  /// buffering multiple message from sendMessage threads.
  // outbound temporary data buffer type for the producer threads.  ThreadSafe version for now.
  //typedef SendMessageBuffers<> BuffersType;
  // outbound temporary data buffers.  one BuffersType per tag value.
  std::unordered_map<int, BuffersType> buffers;

  /// list of tags that are accepting messages.
  // sending threads read, comm thread modify.
  std::unordered_set<int> sendAccept;      // tags


  /// number of senders per tag
  // recv thread modify, callback thread read.
  std::unordered_map<int, int> recvRemaining;   // tag to number of mpi processes


  typedef typename BuffersType::BufferIdType BufferIdType;


  /// condition variable to wake up callbackThread
  //std::condition_variable waitForRecvQueue;


  /// condition variable to make "flush" and "finish" blocking.
  mutable std::mutex mutex;
  mutable std::mutex flushMutex;
  std::condition_variable flushBarrier;
  /// additional check after flushBarrier wakes up, to make sure we are flushing the target tag.
  std::atomic<int> flushing;

  std::atomic<bool> finishing;
  std::atomic<bool> sendDone;
  std::atomic<bool> recvDone;

  std::thread comm_thread;
  std::thread callback_thread;

public:

  CommunicationLayer (const MPI_Comm& communicator, const int comm_size)
    : sendQueue(2 * omp_get_num_threads()), recvQueue(2 * comm_size), flushing(-1), finishing(false), sendDone(false), recvDone(false),
      comm(communicator)
  {
    // init communicator rank and size
    MPI_Comm_size(comm, &commSize);
    assert(comm_size == commSize);
    MPI_Comm_rank(comm, &commRank);
  }

  virtual ~CommunicationLayer () {
    callback_thread.join();
    comm_thread.join();
  };

  // FIXME: rename me
  void startThreads()
  {
    sendDone.store(false);
    recvDone.store(false);
    finishing.store(false);
    flushing.store(-1);
    recvRemaining[0] = commSize;
    comm_thread = std::thread(&CommunicationLayer::commThread, this);
    callback_thread = std::thread(&CommunicationLayer::callbackThread, this);
  }

  // adding the callback function with signature:
  // void(uint8_t* msg, std::size_t count, int fromRank)
  void addReceiveCallback(int tag, std::function<void(uint8_t*, std::size_t, int)> callbackFunction)
  {
    assert(tag != 0);

    if (callbackFunctions.find(tag) != callbackFunctions.end()) {
      printf("function already registered for tag %d\n", tag);
      return;
    }

    /*
    if (recvRemaining[tag] == 0) {
      printf("function already had finished processing for tag %d.\n", tag);
      return;
    }
    */


    if (callbackFunctions.empty())
    {
      // this is the first registered callback, thus spawn the callback
      // executer thread
      // TODO
    }
    // add the callback function to a lookup table
    callbackFunctions[tag] = callbackFunction;

    // also set the number of potential senders
    recvRemaining[tag] = commSize;
    sendAccept.insert(tag);
  }



  /**
   * sends a single message.  the message is buffered and batched before send.
   * may be called by multiple threads concurrently.
   *
   * @param data
   * @param count
   * @param dst_rank
   * @param tag
   */
  void sendMessage(const void* data, std::size_t count, int dst_rank, int tag=default_tag) throw (bliss::io::IOException)
  {


    // check to see if the target tag is still accepting
    if (sendAccept.find(tag) == sendAccept.end()) {
      // TODO:  change to using FATAL on logger.
      printf("ERROR: calling CommunicationLayer::sendMessage with a tag that has been flushed already.  tag=%d\n", tag);
      return;
    }

    /// if there isn't already a tag listed, add the MessageBuffers for that tag.
    // multiple threads may call this.
    if (buffers.find(tag) == buffers.end()) {
      std::unique_lock<std::mutex> lock(mutex);
      if (buffers.find(tag) == buffers.end()) {
        std::cout << std::this_thread::get_id() << " create buffers for tag " << tag << std::endl;
        buffers[tag] = std::move(BuffersType(commSize, 8192));
//        int t2 = tag;
//        buffers.insert(std::move(std::make_pair<int, BuffersType>(std::move(t2), std::move(BuffersType(commSize, 8192)))));
      }
    }
 //   printf("%lu ", count);

    /// try to append the new data - repeat until successful.
    /// along the way, if a full buffer's id is returned, queue it for sendQueue.
    BufferIdType fullId = -1;
    std::pair<bool, BufferIdType> result;
    //bool success = false;
    do {
      result = buffers.at(tag).append(data, count, dst_rank);

      fullId = result.second;

      // repeat until success;
      if (fullId != -1) {
        if (!(buffers.at(tag).getBackBuffer(fullId).isEmpty())) {
        // have a full buffer - put in send queue.
//        SendQueueElement v(fullId, tag, dst_rank);
//        while (!sendQueue.tryPush(std::move(v))) {
//          usleep(20);
//        }
          //printf("full buffer %d\n", fullId );
          if (!sendQueue.waitAndPush(std::move(SendQueueElement(fullId, tag, dst_rank)))) {
            throw bliss::io::IOException("ERROR: sendQueue is not accepting new SendQueueElement due to disablePush");
          }
        }
      }
      fullId = -1;
      usleep(10);
    } while (!result.first);

  }

  /**
   * flushes the message buffers asscoiated with a particular tag.  should be called by a single thread only.
   * @param tag
   */
  void flush(int tag)
  {
    if (sendAccept.find(tag) == sendAccept.end()) {
      // already finished, no need for further flushing
      printf("NO FLUSHing: already finished tag: %d\n", tag);
      return;
    }

    if (flushing.load() != -1) {
      printf("ERROR:  Another tag is being flushed: %d.  target is %d\n", flushing.load(), tag);
      return;
    }

    printf("FLUSHing tag %d.\n", tag);

    // set the tag that is being flushed, this is being reseted in the
    // callback function

    // flush all buffers
    flushBuffers(tag);

    int t = -1;
    if (!flushing.compare_exchange_strong(t, tag)) {
      printf("ERROR:  another tag is being flushed. %d.  target is %d\n", t, tag);
    } else {
      // send the END tags
      sendEndTags(tag);
      // wait till all END tags have been received
      waitForEndTags(tag);

      // reset the count of number of processes active on this tag  ONLY AFTER  waitForEngTags returns
      recvRemaining[tag] = commSize;
    }
  }

  /**
   * flushes the message buffers asscoiated with a particular tag.  should be called by a single thread only.
   * @param tag
   */
  void finishTag(int tag)
  {
    if (sendAccept.find(tag) == sendAccept.end()) {
      printf("ERROR:  NO FINISHing: already finished tag: %d\n", tag);
      // already finished.
      return;
    }

    if (flushing.load() != -1) {
      printf("ERROR:  Another tag is being flushed: %d.  target is %d\n", flushing.load(), tag);
      return;
    }

    printf("FINISHing tag %d.\n", tag);

    /// mark as no more coming in for this tag.  ONE THREAD ONLY modifies this.
    sendAccept.erase(tag);


    // flush all buffers
    flushBuffers(tag);

    int t = -1;
    // set the tag that is being flushed, this is being reseted in the
    // callback function
    if (!flushing.compare_exchange_strong(t, tag)) {
      printf("ERROR:  another tag is being flushed. %d.  target is %d\n", t, tag);
    } else {
      // send the END tags
      sendEndTags(tag);
      // wait till all END tags have been received
      waitForEndTags(tag);
    }


    /// if no more, then send the application termination signal.
    if (sendAccept.empty()) {

      t = -1;
      if (!flushing.compare_exchange_strong(t, 0)) {
        printf("ERROR:  another tag is being flushed. %d.  target is %d\n", t, 0);
      } else {
        sendEndTags(0);
        finishing.store(true);

        waitForEndTags(0);
      }
    }
  }


  // TODO: make this private
  void sendEndTags(int tag) throw (bliss::io::IOException)
  {
    for (int i = 0; i < getCommSize(); ++i) {
      // send end tags in circular fashion, send tag to self last
      int target_rank = (i + getCommRank() + 1) % getCommSize();

      printf("sendEndTags target_rank = %d \n", target_rank );

      // send the end message for this tag.
      if (!sendQueue.waitAndPush(std::move(SendQueueElement(-1, tag, target_rank)))) {
        throw bliss::io::IOException("ERROR: sendQueue is not accepting new SendQueueElement due to disablePush");
      }
    }
  }

  // TODO: make private
  void waitForEndTags(int tag)
  {
    printf("wait for end tag.  flushing= %d, tag= %d\n", flushing.load(), tag);
    assert(flushing.load() == tag);
    // waiting for all messages of this tag to flush.
    std::unique_lock<std::mutex> lock(flushMutex);
    while (flushing.load() == tag) {
      flushBarrier.wait(lock);
    }
  }

  // TODO: make private
  void flushBuffers(int tag) throw (bliss::io::IOException)
  {
    // no MessageBuffers with the associated tag.  end.
    if (buffers.find(tag) == buffers.end())
    {
      printf("NO BUFFERS for tag: %d\n", tag);
      return;
    }
    // flush out all the send buffers matching a particular tag.
    printf("Active Buffer Ids in message buffer: %s\n", buffers.at(tag).activeIdsToString().c_str());

    int idCount = buffers.at(tag).getActiveIds().size();
    assert(idCount == commSize);
    for (int i = 0; i < idCount; ++i) {
      int target_rank = (i + getCommRank()) % getCommSize();
      printf("flushBuffers : target_rank = %d \n", target_rank );

      auto id = buffers.at(tag).getActiveId(target_rank);
      // flush/send all remaining non-empty buffers
      if ((id != -1) && !(buffers.at(tag).getBackBuffer(id).isEmpty())) {
        if (!sendQueue.waitAndPush(std::move(SendQueueElement(id, tag, target_rank)))) {
          throw bliss::io::IOException("ERROR: sendQueue is not accepting new SendQueueElement due to disablePush");
        }
      }
    }
  }

  int getCommSize() const
  {
    return commSize;
  }

  int getCommRank() const
  {
    return commRank;
  }


      /*
      // active receiving (by polling) for when there is no callback set
      // these must be thread-safe!
      Message receiveAnyOne();
      std::vector<message> receiveAnyAll();

      Message receiveOne(tag);
      std::vector<message> receiveAll(tag);
    */



    //
    void commThread()
    {
      // TODO: implement termination condition
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

      printf("comm thread finished on %d\n", commRank);
    }

    void callbackThread()
    {
      // TODO: add termination condition
      ReceivedMessage msg;
      while (!recvDone.load() || !recvQueue.empty() )
      {


        // get next element from the queue, wait if none is available
        // TODO: check if the tag exists as callback function
        auto result = std::move(recvQueue.waitAndPop());
//        auto result = std::move(recvQueue.tryPop());
        if (! result.first) {
          // no valid result.
          continue;
        }

        // else have a valid result.
        msg = result.second;
        if (msg.data == nullptr && msg.count == 0) {   // handles all end messages, including app end message with tag 0 (for which there is no callback)
          // end message.  - once this message is reached, barrier in producer thread (one that called flush or finish) can be breached.
          int tag = msg.tag;
          if (flushing.compare_exchange_strong(tag, -1)) {
            flushBarrier.notify_all();
          } else {
            printf("WARNING: empty message with tag=%d, but flushing contained %d\n", msg.tag, tag); fflush(stdout);
          }
        } else {  // only call callback Function if not an end message.
          // call the matching callback function
          (callbackFunctions[msg.tag])(msg.data, msg.count, msg.src);
          delete [] msg.data;
        }

      }
      printf("recv thread finished on %d\n", commRank);

    }

  protected:

    //// Order of message matters so make sure end message will be sent/recv'd after data messages.
    //// same source/destination messages goes between queues directly, bypassing MPI, so need to enqueue local termination message into "inprogress" (data can go into recv queue)


    void tryStartSend() throw (bliss::io::IOException)
    {
      if (sendDone.load()) return;


      // try to get the send element
      SendQueueElement se;
      auto el = std::move(sendQueue.tryPop());
//      auto el = std::move(sendQueue.waitAndPop());
      if (!el.first) {
        return;   // sendQueue disabled.
      }

      // else there is a valid entry from sendQueue
      se = el.second;
      if (se.bufferId == -1) {
        // termination message for this tag and destination


        if (se.dst == commRank) {
          // local, directly handle by creating an output object and directly
          // insert into the recvInProgress (need to use the receiver decrement logic in "finishReceive")
          MPI_Request req = MPI_REQUEST_NULL;
          recvInProgress.push(std::move(std::pair<MPI_Request, ReceivedMessage>(std::move(req), std::move(ReceivedMessage(nullptr, 0, se.tag, commRank)))));

        } else {
          // remote.  send a terminating message. with the same tag to ensure message ordering.
          MPI_Request req;
          MPI_Isend(nullptr, 0, MPI_UINT8_T, se.dst, se.tag, comm, &req);
          sendInProgress.push(std::move(std::pair<MPI_Request, SendQueueElement>(std::move(req), std::move(se))));
        }

      } else {  // real data.
        void* data = const_cast<void*>(buffers.at(se.tag).getBackBuffer(se.bufferId).getData());
        auto count = buffers.at(se.tag).getBackBuffer(se.bufferId).getSize();

//          printf ("sending %d\n", se.bufferId);
        if (se.dst == commRank) {
          // local, directly handle by creating an output object and directly
          // insert into the recv queue (thread safe)
          uint8_t* array = new uint8_t[count];
          memcpy(array, data, count);

          if (!recvQueue.waitAndPush(std::move(ReceivedMessage(array, count, se.tag, commRank)))) {
            throw bliss::io::IOException("ERROR: recvQueue is not accepting new receivedMessage due to disablePush");
          }

          // finished inserting directly to local RecvQueue.  release the buffer
          buffers.at(se.tag).releaseBuffer(se.bufferId);
//            printf("released %d.  recvQueue size is %lu\n", se.bufferId, recvQueue.size());

        } else {

          MPI_Request req;
          MPI_Isend(data, count, MPI_UINT8_T, se.dst, se.tag, comm, &req);
          sendInProgress.push(std::move(std::pair<MPI_Request, SendQueueElement>(std::move(req), std::move(se))));
        }
      }

    }

    void finishSends()
    {
      if (sendDone.load()) return;

      int finished = 0;
      while(!sendInProgress.empty())
      {
        std::pair<MPI_Request, SendQueueElement>& front = sendInProgress.front();

        // check if MPI request has completed.
        if (front.first == MPI_REQUEST_NULL) finished = 1;
        else
          MPI_Test(&front.first, &finished, MPI_STATUS_IGNORE);

        if (finished)
        {
          if (front.second.bufferId != -1) {
            // cleanup, i.e., release the buffer back into the pool
            buffers.at(front.second.tag).releaseBuffer(front.second.bufferId);
          }
          sendInProgress.pop();
        }
        else
        {
          break;
        }
      }

      if (finishing.load() && sendQueue.empty() && sendInProgress.empty()) {
        // not using sendAccept as check - sendAccept is cleared BEFORE the last tag end messages starts to be sent, and definitely before the app end message starts to be sent
        // so sendDone may be set to true too early, and elements for the last tag, or tag 0 (termination) may be added to the sendQueue or sendInProfess queue after the empty check here.
        sendDone.store(true);
        printf("send Done!\n");
//        sendQueue.disablePush();
      }
    }

    void tryStartReceive()
    {
      if (recvDone.load()) return;

      /// probe for messages
      int hasMessage = 0;
      MPI_Status status;
      MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, comm, &hasMessage, &status);

      if (hasMessage > 0) {
        int src = status.MPI_SOURCE;
        int tag = status.MPI_TAG;
        int received_count;
        MPI_Get_count(&status, MPI_UNSIGNED_CHAR, &received_count);

        // receive whether it's empty or not.  use finishReceive to handle the termination message and count decrement

        // receive the message data bytes into a vector of bytes
        uint8_t* msg_data = (received_count == 0) ? nullptr : new uint8_t[received_count];
        MPI_Request req;
        MPI_Irecv(msg_data, received_count, MPI_UNSIGNED_CHAR, src, tag, comm, &req);

        // insert into the in-progress queue
        recvInProgress.push(std::pair<MPI_Request, ReceivedMessage>(req,
                                                                 ReceivedMessage(msg_data, received_count, tag, src)));

      }
    }


    // Not thread safe
    void finishReceives() throw (bliss::io::IOException)
    {
      if (recvDone.load()) return;

      int finished = 0;
      //MPI_Status status;
      while(!recvInProgress.empty())
      {
        std::pair<MPI_Request, ReceivedMessage>& front = recvInProgress.front();

        // test if the MPI request for receive has been completed.
        if (front.first == MPI_REQUEST_NULL) finished = 1;
        else
          MPI_Test(&front.first, &finished, MPI_STATUS_IGNORE);

        if (finished)
        {

          if (front.second.count == 0) {  // terminating message.
            // end of messaging.
            --recvRemaining[front.second.tag];
            //printf("RECV rank %d receiving END signal %d from %d, num senders remaining is %d\n", commRank, front.second.tag, front.second.src, recvRemaining[front.second.tag]);

            if (recvRemaining[front.second.tag] == 0) {
              // received all end messages.  there may still be message in progress and in recvQueue from this and other sources.
              recvRemaining.erase(front.second.tag);

              printf("ALL END received for tag %d, pushing to recv queue\n", front.second.tag);
              fflush(stdout);
              if (!recvQueue.waitAndPush(std::move(front.second))) {
                throw bliss::io::IOException("ERROR: recvQueue is not accepting new receivedMessage due to disablePush");
              }

            } else if (recvRemaining[front.second.tag] < 0) {
              printf("ERROR: number of remaining receivers for tag %d is now NEGATIVE\n", front.second.tag);
            }

          } else {

          // add the received messages into the recvQueue
    //        ReceivedMessage msg(std::get<1>(front), std::get<2>(front),
    //                            status.MPI_TAG, status.MPI_SOURCE);
          // TODO: (maybe) one queue per tag?? Where does this
          //       sorting/categorizing happen?
//          while (!recvQueue.tryPush(std::move(front.second))) {
//            usleep(50);
//          }
            if (!recvQueue.waitAndPush(std::move(front.second))) {
              throw bliss::io::IOException("ERROR: recvQueue is not accepting new receivedMessage due to disablePush");
            }
          }
          // remove moved element
          recvInProgress.pop();
        }
        else
        {
          // stop receiving for now
          break;
        }
      }

      // see if we are done with receiving and everything is in the recvQueue
      if ((recvRemaining.find(0) == recvRemaining.end()) && recvInProgress.empty()) {
        recvDone.store(true);
        printf("recv Done!\n");
        recvQueue.disablePush();
      }
    }


private:
  /* data */

  /// The MPI Communicator object for this communication layer
  MPI_Comm comm;

  /// Registry of callback functions, mapped to by the associated tags
  std::map<int,std::function<void(uint8_t*, std::size_t, int)> > callbackFunctions;

  /// The MPI Communicator size
  int commSize;

  /// The MPI Communicator rank
  int commRank;
};

#endif // BLISS_COMMUNICATION_LAYER_HPP
