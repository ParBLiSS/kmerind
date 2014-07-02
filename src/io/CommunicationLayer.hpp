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
  // request, data pointer, data size
  std::queue<std::pair<MPI_Request, ReceivedMessage> > recvInProgress;
  // request, data pointer, tag
  std::queue<std::pair<MPI_Request, SendQueueElement> > sendInProgress;

  // Outbound message structure, Multiple-Producer-Single-Consumer queue
  // consumed by the internal MPI-comm thread
  bliss::concurrent::ThreadSafeQueue<SendQueueElement> sendQueue;

  // Inbound message structure, Single-Producer-Multiple-Consumer queue
  // produced by the internal MPI-comm thread
  bliss::concurrent::ThreadSafeQueue<ReceivedMessage> recvQueue;

  // outbound temporary data buffer type for the producer threads.  ThreadSafe version for now.
  //typedef SendMessageBuffers<> BuffersType;
  // outbound temporary data buffers.  one BuffersType per tag value.
  std::unordered_map<int, BuffersType> buffers;

  std::unordered_set<int> sendAccept;      // tags
  std::unordered_map<int, int> recvRemaining;   // tag to number of mpi processes


  typedef typename BuffersType::BufferIdType BufferIdType;

  mutable std::mutex mutex;
  std::condition_variable cond_var;
  std::atomic<int> flushing;

  std::thread comm_thread;
  std::thread callback_thread;

public:

  CommunicationLayer (const MPI_Comm& communicator, const int comm_size)
    : sendQueue(2 * omp_get_num_threads()), recvQueue(2 * comm_size), flushing(-1),
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
    comm_thread = std::thread(&CommunicationLayer::commThread, this);
    callback_thread = std::thread(&CommunicationLayer::callbackThread, this);
  }

  // adding the callback function with signature:
  // void(uint8_t* msg, std::size_t count, int fromRank)
  void addReceiveCallback(int tag, std::function<void(uint8_t*, std::size_t, int)> callbackFunction)
  {
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
  void sendMessage(const void* data, std::size_t count, int dst_rank, int tag=default_tag)
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
        int t2 = tag;
        buffers[tag] = std::move(BuffersType(commSize, 8192));
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
          sendQueue.waitAndPush(std::move(SendQueueElement(fullId, tag, dst_rank)));
        }
      }
      fullId = -1;
      usleep(20);
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

    printf("FLUSHing tag %d.\n", tag);

    // set the tag that is being flushed, this is being reseted in the
    // callback function
    flushing.store(tag);

    // flush all buffers
    flushBuffers(tag);
    // send the END tags
    sendEndTags(tag);
    // wait till all END tags have been received
    waitForEndTags(tag);

    // reset the count of number of processes active on this tag
    recvRemaining[tag] = commSize;
  }

  /**
   * flushes the message buffers asscoiated with a particular tag.  should be called by a single thread only.
   * @param tag
   */
  void finishTag(int tag)
  {
    if (sendAccept.find(tag) == sendAccept.end()) {
      printf("NO FINISHing: already finished tag: %d\n", tag);
      // already finished.
      return;
    }

    printf("FINISHing tag %d.\n", tag);

    /// mark as no more coming in for this tag.
    sendAccept.erase(tag);

    // set the tag that is being flushed, this is being reseted in the
    // callback function
    flushing.store(tag);

    // flush all buffers
    flushBuffers(tag);
    // send the END tags
    sendEndTags(tag);
    // wait till all END tags have been received
    waitForEndTags(tag);
  }


  // TODO: make this private
  void sendEndTags(int tag)
  {
    for (int i = 0; i < getCommSize(); ++i) {
      // send end tags in circular fashion, send tag to self last
      int target_rank = (i + getCommRank() + 1) % getCommSize();
      // send the end message for this tag.
      sendQueue.waitAndPush(std::move(SendQueueElement(-1, tag, target_rank)));
    }
  }

  // TODO: make private
  void waitForEndTags(int tag)
  {
    assert(flushing == tag);
    // waiting for all messages of this tag to flush.
    std::unique_lock<std::mutex> lock(mutex);
    while (flushing == tag) {
      cond_var.wait(lock);
    }
  }

  // TODO: make private
  void flushBuffers(int tag)
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
      auto id = buffers.at(tag).getActiveId(target_rank);
      // flush/send all remaining non-empty buffers
      if ((id != -1) && !(buffers.at(tag).getBackBuffer(id).isEmpty())) {
        sendQueue.waitAndPush(std::move(SendQueueElement(id, tag, target_rank)));
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
      while ((sendAccept.size() > 0) || (sendQueue.size() > 0) || (sendInProgress.size() > 0)
              || (recvRemaining.size() > 0) || (recvInProgress.size() > 0) || (recvQueue.size() > 0))
      {
        // first clean up finished operations
        finishReceives();
        finishSends();

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
      while ((recvRemaining.size() > 0) || (recvInProgress.size() > 0) || (recvQueue.size() > 0))
      {
        // get next element from the queue, wait if none is available
        // TODO: check if the tag exists as callback function
        auto result = std::move(recvQueue.tryPop());
        if (result.first) {
          msg = result.second;
          if (msg.data == nullptr && msg.count == 0) {
            if (flushing == msg.tag) {
              std::unique_lock<std::mutex> lock(mutex);
              flushing = -1;
              lock.unlock();
              cond_var.notify_one();
            } else {
              printf("WARNING: empty message with flushing=%d, but looking for %d\n", flushing.load(), msg.tag);
            }
          }

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


    void tryStartSend()
    {
      // try to get the send element
      SendQueueElement se;
      auto el = std::move(sendQueue.tryPop());
      if (el.first) {
        se = el.second;
        if (se.bufferId == -1) {
          // termination message for this tag and destination


          if (se.dst == commRank) {
            // local, directly handle by creating an output object and directly
            // insert into the recvInProgress (need to use the receiver decrement logic in "finishReceive")
            MPI_Request req = MPI_REQUEST_NULL;
            recvInProgress.push(std::pair<MPI_Request, ReceivedMessage>(req, ReceivedMessage(nullptr, 0, se.tag, commRank)));

          } else {
            // remote.  send a terminating message. with the same tag to ensure message ordering.
            MPI_Request req;
            MPI_Isend(nullptr, 0, MPI_UINT8_T, se.dst, se.tag, comm, &req);
            sendInProgress.push(std::move(std::pair<MPI_Request, SendQueueElement>(req, se)));
          }

        } else {  // real data.
          void* data = const_cast<void*>(buffers.at(se.tag).getBackBuffer(se.bufferId).getData());
          auto count = buffers.at(se.tag).getBackBuffer(se.bufferId).getSize();

//          printf ("sending %d\n", se.bufferId);
          if (se.dst == commRank) {
            // local, directly handle by creating an output object and directly
            // insert into the recv queue
            uint8_t* array = new uint8_t[count];
            memcpy(array, data, count);

            recvQueue.waitAndPush(std::move(ReceivedMessage(array, count, se.tag, commRank)));

            // finished inserting directly to local RecvQueue.  release the buffer
            buffers.at(se.tag).releaseBuffer(se.bufferId);
//            printf("released %d.  recvQueue size is %lu\n", se.bufferId, recvQueue.size());

          } else {

            MPI_Request req;
            MPI_Isend(data, count, MPI_UINT8_T, se.dst, se.tag, comm, &req);
            sendInProgress.push(std::move(std::pair<MPI_Request, SendQueueElement>(req, se)));
          }
        }
      }
    }

    void finishSends()
    {
      int finished = 0;
      while(!sendInProgress.empty())
      {
        std::pair<MPI_Request, SendQueueElement>& front = sendInProgress.front();
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
    }

    void tryStartReceive()
    {
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
    void finishReceives()
    {
      int finished = 0;
      //MPI_Status status;
      while(!recvInProgress.empty())
      {
        std::pair<MPI_Request, ReceivedMessage>& front = recvInProgress.front();

        if (front.first == MPI_REQUEST_NULL) finished = 1;
        else
          MPI_Test(&front.first, &finished, MPI_STATUS_IGNORE);

        if (finished)
        {

          if (front.second.count == 0) {  // terminating message.
            // end of messaging.
            --recvRemaining[front.second.tag];
            printf("RECV rank %d receiving END signal %d from %d, num senders remaining is %d\n", commRank, front.second.tag, front.second.src, recvRemaining[front.second.tag]);

            if (recvRemaining[front.second.tag] == 0) {
              // received all end messages.  there may still be message in progress and in recvQueue from this and other sources.
              recvRemaining.erase(front.second.tag);

              printf("ALL END received, pushing to recv queue");
              fflush(stdout);
              recvQueue.waitAndPush(std::move(front.second));

            } else if (recvRemaining[front.second.tag] == 0) {
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
            recvQueue.waitAndPush(std::move(front.second));
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
