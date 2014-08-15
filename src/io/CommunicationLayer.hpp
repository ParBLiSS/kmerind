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
 * TODO: cleanup
 *  - replace uint8_t by typedef
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
#include <unistd.h> // for usleep() FIXME: replace by std::this_thread::sleep_for

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


// Use the threadsafe version of the SendMessageBuffers as message buffer
typedef SendMessageBuffers<bliss::concurrent::THREAD_SAFE> BuffersType;
typedef typename BuffersType::BufferPoolType                BufferPoolType;


/**
 * @brief Structure to hold all associated information of a received MPI
 *        message.
 */
struct ReceivedMessage
{
  /// The received data
  uint8_t* data;
  /// The number of bytes received
  std::size_t count;
  /// The message tag
  int tag;
  /// The message source id
  int src;

  /**
   * @brief Constructs a new instance of this struct.
   *
   * Sets all data members to the given values.
   */
  ReceivedMessage(uint8_t* data, std::size_t count, int tag, int src)
    : data(data), count(count), tag(tag), src(src) {}

  // Enable default construction
  ReceivedMessage() = default;

  int getTag() const {
    if (tag == 0)
      return ((int*)data)[0];
    else
      return tag;
  }

  int getEpoch() const {
    if (tag == 0)
      return ((int*)data)[1];
    else
      return -1;
  }
};

// TODO: rename (either this or the ReceivedMessage) for identical naming scheme
/**
 * @brief Structure to hold all associated information of a message to be sent.
 */
struct SendQueueElement
{
  /// The BufferID for the message
  typename BuffersType::BufferIdType bufferId;
  /// Tag of the message
  int tag;
  /// Destination rank, i.e., id of the process where the message is to be sent
  int dst;

  /**
   * @brief Constructs a new instance of this struct and sets all members as
   * given.
   *
   * @param _id
   * @param _tag
   * @param _dst
   */
  SendQueueElement(BuffersType::BufferIdType _id, int _tag, int _dst)
    : bufferId(_id), tag(_tag), dst(_dst) {}
  // enable default construction
  SendQueueElement() = default;
};


/**
 * @brief Abstracts asynchronous and buffered communication between all
 *        processes via MPI.
 *
 * Prior to using the sendMessage(), flush(), and finishTag() functions,
 * the CommunicationLayer first needs to be initialized by calling the
 *      `initCommunication()`
 * function.
 *
 *
 * Design requirements:
 *  asynchronous MPI messaging using isend/iprobe/irecv
 *  message queues for sent and recved message handling
 *  tag for each message type
 *  tag = 0 means application termination
 *  flushTag() inserts zero byte messages with specific tag to indicate that processes need to wait until all messages of that tag type are done.
 *  finishTag() also sends application termination tag.
 *
 *  each process handles the flush and finish independently without barrier synchronization.
 *  thread calling flush and finish essentially "request" that "barriers" be erected.
 *    Recv Thread handles the actual within-process barrier to ensure all messages of a particular tag are processed.
 *    flushTag() resets the counter
 *
 *  potential issue with using just counters:  if multiple phases for a particular tag's messages, and one process is particularly slow,
 *    it may receive flush signal from multiple phases of other processes, e.g. process 0 in phase 1,
 *    and receives flush signal from process 1 from phase 1 and 2.  Using just counter could cause process 0 to not properly
 *    barrier each phase.
 */
class CommunicationLayer
{

public:

    static const int CONTROL_TAG = 0;


  /**
   * @brief Constructor for the CommLayer.
   *
   * @param communicator    The MPI Communicator used in this communicator
   *                        layer.
   * @param comm_size       The size of the MPI Communicator.
   */
  CommunicationLayer (const MPI_Comm& communicator, const int comm_size)
    : sendQueue(2 * omp_get_num_threads()), recvQueue(2 * comm_size),
      flushing(-1), finishing(false), sendDone(false), recvDone(false),
      comm(communicator), epoch(0)
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
   *
   * Blocks till all threads are finished. In case finish() has not been called
   * for all used tags, this will block forever.
   */
  virtual ~CommunicationLayer () {
    // TODO: check that we are finished??
    // wait for both threads to quit
    callback_thread.join();
    comm_thread.join();

    // detach the buffer
    int mpiBufSize;
    char* mpiBuf;
    MPI_Buffer_detach(&mpiBuf, &mpiBufSize);
    free(mpiBuf);
  };


  /**
   * @brief Initializes all communication and starts the communication and
   *        callback threads.
   *
   * Must be called before any other functions can be called.
   */
  void initCommunication()
  {
    sendDone.store(false);
    recvDone.store(false);
    finishing.store(false);
    flushing.store(-1);
    recvRemaining[-1] = commSize;
    comm_thread = std::thread(&CommunicationLayer::commThread, this);
    callback_thread = std::thread(&CommunicationLayer::callbackThread, this);
  }


  /**
   * @brief Sends a message (raw bytes) to the given rank asynchronously.
   *
   * Asynchronously sends the given message to the given rank with the given
   * message tag. The message is buffered and batched before send. This function
   * returns immediately after buffering, which does not guarantuee that the
   * message has been sent yet. Only calling `flush()` guarantuees that all
   * messages are sent, received, and processed before the function returns.
   *
   * This function is thread-safe.
   *
   * @param data        A pointer to the data to be sent.
   * @param count       The number of bytes of the message.
   * @param dst_rank    The rank of the destination processor.
   * @param tag         The tag for the message.
   */
  void sendMessage(const void* data, std::size_t count, int dst_rank, int tag)
  {
    // check to see if the target tag is still accepting
    if (sendAccept.find(tag) == sendAccept.end()) {
      throw std::invalid_argument("invalid tag: tag has been finished already");
    }

    assert(tag != CONTROL_TAG && tag >= 0);

    /// if there isn't already a tag listed, add the MessageBuffers for that tag.
    // multiple threads may call this.
    if (buffers.find(tag) == buffers.end()) {
      std::unique_lock<std::mutex> lock(buffers_mutex);
      if (buffers.find(tag) == buffers.end()) {
        // create new message buffer
        buffers[tag] = std::move(BuffersType(commSize, 8192));
      }
    }

    // try to append the new data - repeat until successful.
    // along the way, if a full buffer's id is returned, queue it for sendQueue.
    BufferIdType fullId = BufferPoolType::ABSENT;
    std::pair<bool, BufferIdType> result;
    do {
      result = buffers.at(tag).append(data, count, dst_rank);

      fullId = result.second;

      // repeat until success;
      if (fullId != BufferPoolType::ABSENT) {
        if (!(buffers.at(tag).getBackBuffer(fullId).isEmpty())) {
          // have a full buffer - put in send queue.
          if (!sendQueue.waitAndPush(std::move(SendQueueElement(fullId, tag, dst_rank)))) {
            throw bliss::io::IOException("ERROR: sendQueue is not accepting new SendQueueElement due to disablePush");
          }
        }
      }
      fullId = BufferPoolType::ABSENT;
    } while (!result.first);
  }

  /**
   * @brief Registers a callback function for the given message tag.
   *
   * The registered callback function will be called for each received message
   * of the given type. A callback function has to be registered for each
   * message type sent. Only one callback function can be registered for each
   * tag, adding more than one results in a `invalid_argument` expception begin
   * thrown.
   *
   * A callback function has the signature:
   *    void(uint8_t* msg, std::size_t count, int fromRank)
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

    // also set the number of potential senders
    // now accepting messages of this type
    sendAccept.insert(tag);
  }

  /**
   * flushes the message buffers asscoiated with a particular tag.  should be called by a single thread only.
   * @param tag
   */
  void flush(int tag)
  {
    assert(tag != CONTROL_TAG && tag >= 0);

    if (sendAccept.find(tag) == sendAccept.end()) {
      // already finished, no need for further flushing
      DEBUGF("NO FLUSHing: already finished tag: %d", tag);
      return;
    }

    INFOF("FLUSHing tag %d, epoch %d", tag, epoch);


    // flush all buffers (put them into send queue)
    flushBuffers(tag);

    // set the tag that is being flushed, this is being reseted in the
    // callback function
    int t = -1;
    if (!flushing.compare_exchange_strong(t, tag)) {
      // TODO throw exception?
      ERRORF("ERROR:  another tag is being flushed. %d.  target is %d", t, tag);
      return;
    }

    // send the END tags
    sendEndTags(tag);
    // wait till all END tags have been received (not using Barrier because comm layer still needs to buffer messages locally)
    waitForEndTags(tag);

    // reset the count of number of processes active on this tag  ONLY AFTER  waitForEngTags returns
//    DEBUGF("FLUSHED %d tag.  reseting RecvRemaining to commsize.", commRank);
//    recvRemaining[tag] = commSize;
  }

  /**
   * flushes the message buffers asscoiated with a particular tag.  should be called by a single thread only.
   * @param tag
   */
  void finishTag(int tag)
  {
    assert(tag != CONTROL_TAG && tag >= 0);

    if (sendAccept.find(tag) == sendAccept.end()) {
      //DEBUGF("ERROR:  NO FINISHing: already finished tag: %d", tag);
      // already finished.
      return;
    }


    INFOF("FINISHing tag %d, epoch %d", tag, epoch);

    /// mark as no more coming in for this tag.  ONE THREAD ONLY modifies this.
    sendAccept.erase(tag);

    // flush all buffers
    flushBuffers(tag);

    // set the tag that is being flushed, this is being reseted in the
    // callback function
    int t = -1;
    if (!flushing.compare_exchange_strong(t, tag)) {
      ERRORF("another tag is being flushed. %d.  target is %d", t, tag);
      return;
    }
    // send the END tags
    sendEndTags(tag);
    // wait till all END tags have been received
    waitForEndTags(tag);

    /// if no more, then send the application termination signal.
    if (sendAccept.empty()) {
      t = -1;
      if (!flushing.compare_exchange_strong(t, CONTROL_TAG)) {
        ERRORF("ERROR:  another tag is being flushed. %d.  target is %d", t, CONTROL_TAG);
        return;
      }
      sendEndTags(CONTROL_TAG);
      finishing.store(true);

      waitForEndTags(CONTROL_TAG);

    }
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


  /**
   * @brief The function for the dedicated communicator thread.
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
   */
  void callbackThread()
  {
    INFOF("THREAD STARTED:  recv-thread on %d", commRank);
    ReceivedMessage msg;
    // while there is still work to be done:
    while (!recvDone.load() || !recvQueue.isEmpty() )
    {
      // get next element from the queue, wait if none is available
      auto result = std::move(recvQueue.waitAndPop());
      if (! result.first) {
        // no valid result, we must be done with receiving
        assert(recvDone.load());
        break;
      }

      // get the message
      msg = result.second;
      // if the message is emtpy -> END-TAG
      if (msg.tag == CONTROL_TAG) {
        // handles all end messages, including app end message with tag 0 (for
        // which there is no callback)
        int tag = msg.getTag();
        // expectation: the received message tag is the current tag that is
        // being flushed or finished
        if (flushing.compare_exchange_strong(tag, -1)) {
          flushBarrier.notify_all();
        } else {
          ERRORF("empty message with tag=%d, but flushing contained %d", msg.tag, tag);
        }
      } else {
        // only call callback Function if not an end message.
        (callbackFunctions[msg.tag])(msg.data, msg.count, msg.src);
        // clean up message data
        delete [] msg.data;
      }
    }
    INFOF("THREAD FINISHED: recv-thread on %d", commRank);
  }

protected:

  /**
   * @brief Sends END tags for the given tag to every MPI process
   *
   * @param tag     The MPI tag to end.
   *
   * @throw bliss::io::IOException  If the sendQueue has been disabled.
   */
  void sendEndTags(int tag) throw (bliss::io::IOException)
  {
    for (int i = 0; i < getCommSize(); ++i) {
      // send end tags in circular fashion, send tag to self last
      int target_rank = (i + getCommRank() + 1) % getCommSize();

      DEBUGF("Rank %d sendEndTags %d epoch %d target_rank = %d ", commRank, tag, epoch, target_rank );

      // send the end message for this tag.
      if (!sendQueue.waitAndPush(std::move(SendQueueElement(BufferPoolType::ABSENT, tag, target_rank)))) {
        throw bliss::io::IOException("ERROR: sendQueue is not accepting new SendQueueElement due to disablePush");
      }
    }
  }


  /**
   * @brief Blocks till the flusing of the given tag is completed (i.e., all
   * END tags have been received).
   *
   * @param tag The tag to wait for.
   */
  void waitForEndTags(int tag) throw (bliss::io::IOException)
  {

    // waiting for all messages of this tag to flush.
    if (flushing.load() == tag) {  // precheck before locking.
      std::unique_lock<std::mutex> lock(flushMutex);
      do {
        flushBarrier.wait(lock);
      } while (flushing.load() == tag);
    } else {
      throw bliss::io::IOException("ERROR: waitForEndTags called but currently flushing a different tag.");
    }
    DEBUGF("Rank %d WAITED FOR END TAG.  flushing= %d, tag= %d, epoch = %d", commRank, flushing.load(), tag, epoch);
    ++epoch;
  }


  /**
   * @brief Flushes all buffers for message tag `tag`.
   *
   * All pending/buffered messages for tag `tag` will be send out.
   *
   * @param tag The message tag, whose messages are flushed.
   *
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
    // flush out all the send buffers matching a particular tag.
    //DEBUGF("Active Buffer Ids in message buffer: %s", buffers.at(tag).bufferIdsToString().c_str());

    int idCount = buffers.at(tag).getBufferIds().size();
    assert(idCount == commSize);
    for (int i = 0; i < idCount; ++i) {
      // flush buffers in a circular fashion, starting with the next neighbor
      int target_rank = (i + getCommRank()) % getCommSize();
      //DEBUGF("flushBuffers : target_rank = %d ", target_rank );

      auto id = buffers.at(tag).getBufferId(target_rank);
      // flush/send all remaining non-empty buffers
      if ((id != BufferPoolType::ABSENT) && !(buffers.at(tag).getBackBuffer(id).isEmpty())) {
        if (!sendQueue.waitAndPush(std::move(SendQueueElement(id, tag, target_rank)))) {
          throw bliss::io::IOException("ERROR: sendQueue is not accepting new SendQueueElement due to disablePush");
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
    SendQueueElement se;
    auto el = std::move(sendQueue.tryPop());
    //auto el = std::move(sendQueue.waitAndPop());
    if (!el.first) {
      return;   // sendQueue disabled.
    }

    // else there is a valid entry from sendQueue
    se = el.second;
    // Is this an END tag/ empty message?
    if (se.bufferId == BufferPoolType::ABSENT) {
      //DEBUGF("SEND %d -> %d, termination signal for tag %d", commRank, se.dst, se.tag);
      // termination message for this tag and destination
      if (se.dst == commRank) {
        // local, directly handle by creating an output object and directly
        // insert into the recvInProgress (need to use the receiver decrement
        // logic in "finishReceive")
        MPI_Request req = MPI_REQUEST_NULL;
        uint8_t *array = new uint8_t[2 * sizeof(int)];
        memcpy(array, &(se.tag), sizeof(int));
        memcpy(array + sizeof(int), &(epoch), sizeof(int));
        recvInProgress.push(std::move(std::pair<MPI_Request,
              ReceivedMessage>(std::move(req),
                std::move(ReceivedMessage(array, 2 * sizeof(int), CONTROL_TAG, commRank)))));
      } else {
        // remote.  send a terminating message. with the same tag to ensure
        // message ordering.
        int stat[2] = {se.tag, epoch};
        MPI_Bsend(stat, 2, MPI_INT, se.dst, CONTROL_TAG, comm);
      }
    // this is an actual message:
    } else {

      // get message data and it's size
      void* data = const_cast<void*>(buffers.at(se.tag).getBackBuffer(se.bufferId).getData());
      auto count = buffers.at(se.tag).getBackBuffer(se.bufferId).getSize();

      if (count == 0) {
        ERRORF("ERROR: SEND %d -> %d, trying to send 0 byte message as data for tag %d.", commRank, se.dst, se.tag);
      }

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
      } else {
        // remote: initiate async MPI message
        MPI_Request req;
        MPI_Isend(data, count, MPI_UINT8_T, se.dst, se.tag, comm, &req);
        sendInProgress.push(std::move(std::pair<MPI_Request, SendQueueElement>(std::move(req), std::move(se))));
      }
    }

  }

  /**
   * @brief Finishes all completed, pending MPI_Send requests.
   */
  void finishSends()
  {
    if (sendDone.load()) return;

    int finished = 0;
    while(!sendInProgress.empty())
    {
      std::pair<MPI_Request, SendQueueElement>& front = sendInProgress.front();

      // check if MPI request has completed.
      if (front.first == MPI_REQUEST_NULL) {
        finished = 1;
      } else {
        MPI_Test(&front.first, &finished, MPI_STATUS_IGNORE);
      }

      if (finished)
      {
        if (front.second.bufferId != BufferPoolType::ABSENT) {
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
    if (recvDone.load()) return;

    /// probe for messages
    int hasMessage = 0;
    MPI_Status status;
    MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, comm, &hasMessage, &status);

    if (hasMessage > 0) {
      // get message details
      int src = status.MPI_SOURCE;
      int tag = status.MPI_TAG;
      int received_count;
      MPI_Get_count(&status, MPI_UNSIGNED_CHAR, &received_count);

      // receive whether it's empty or not.  use finishReceive to handle the
      // termination message and count decrement

      uint8_t* msg_data = (received_count == 0) ? nullptr : new uint8_t[received_count];
      MPI_Request req;

      // receive the message data bytes into a vector of bytes
      if (tag == CONTROL_TAG) {
        MPI_Irecv(msg_data, received_count / sizeof(int), MPI_INT, src, tag, comm, &req);

      } else {
        MPI_Irecv(msg_data, received_count, MPI_UNSIGNED_CHAR, src, tag, comm, &req);
      }
      // insert into the in-progress queue
      recvInProgress.push(std::pair<MPI_Request, ReceivedMessage>(req,
            ReceivedMessage(msg_data, received_count, tag, src)));
    }
  }


  /**
   * @brief Finishes pending but completed receives.
   *
   * @throw bliss::io::IOException
   */
  void finishReceives() throw (bliss::io::IOException)
  {
    if (recvDone.load()) return;

    int finished = 0;
    while(!recvInProgress.empty())
    {
      std::pair<MPI_Request, ReceivedMessage>& front = recvInProgress.front();

      // test if the MPI request for receive has been completed.
      if (front.first == MPI_REQUEST_NULL) {
        finished = 1;   // this happens with local request.
      } else {
        MPI_Test(&front.first, &finished, MPI_STATUS_IGNORE);
      }

      if (finished)
      {
        if (front.second.tag == CONTROL_TAG)  {
          // control message

          int tag = front.second.getTag();
          int epoch = front.second.getEpoch();

          if (tag == CONTROL_TAG) {
            epoch = -1;  // already set at construction
          } else {
            // if we haven't seen this epoch, add a new entry.
            if (recvRemaining.find(epoch) == recvRemaining.end())
              recvRemaining[epoch] = commSize;  // insert new.
          }

          // termination message
          --recvRemaining.at(epoch);
          DEBUGF("RECV rank %d receiving END signal tag %d epoch %d from %d, num senders remaining is %d",
                 commRank, tag, epoch, front.second.src, recvRemaining.at(epoch));

          if (recvRemaining.at(epoch) == 0) {
            // received all end messages.  there may still be messages in
            // progress and in recvQueue from this and other sources.

            recvRemaining.erase(epoch);

            //DEBUGF("ALL END received for tag %d, pushing to recv queue", front.second.tag);
            fflush(stdout);
            if (!recvQueue.waitAndPush(std::move(front.second))) {
              throw bliss::io::IOException("ERROR: recvQueue is not accepting new receivedMessage due to disablePush");
            }

          } else if (recvRemaining.at(epoch) < 0) {
            ERRORF("ERROR: number of remaining receivers for tag %d epoch %d is now NEGATIVE", tag, epoch);
          }
        } else {
        // add the received messages into the recvQueue
          if (!recvQueue.waitAndPush(std::move(front.second))) {
            throw bliss::io::IOException("ERROR: recvQueue is not accepting new receivedMessage due to disablePush");
          }
        }
        // remove moved element
        recvInProgress.pop();
      }
      else
      {
        // stop receiving for now - wait for requests to finish
        break;
      }
    }

    // see if we are done with receiving and everything is in the recvQueue
    if ((recvRemaining.size() == 0) && recvInProgress.empty()) {
      // if completely done, epoch=-1 would have been erased from recvRemaining.
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
  std::queue<std::pair<MPI_Request, ReceivedMessage> > recvInProgress;
  /// Queue of pending MPI send operations
  std::queue<std::pair<MPI_Request, SendQueueElement> > sendInProgress;


  /* Thread-safe queues for posting sends and getting received messages,
   * these are used as communication medium between the comm-thread and the
   * (potetially multiple) producer threads for the sends; and the comm-thread
   * and the callback handler thread for the receives
   */

  /// message queue between sendMessage calling threads (src) and comm-thread
  /// (sink)
  // Outbound message structure, Multiple-Producer-Single-Consumer queue
  // consumed by the internal MPI-comm-thread
  bliss::concurrent::ThreadSafeQueue<SendQueueElement> sendQueue;

  /// message queue between comm thread (src) and callback thread (sink)
  // Inbound message structure, Single-Producer-Multiple-Consumer queue
  // produced by the internal MPI-comm-thread
  bliss::concurrent::ThreadSafeQueue<ReceivedMessage> recvQueue;


  /* information per message tag */

  /// Message buffers per message tag (maps each tag to a buffer)
  std::unordered_map<int, BuffersType> buffers;
  /// Mutex for adding new tags to the buffers map (since inserts can be
  /// multi-threaded)
  mutable std::mutex buffers_mutex;

  /// set of active message tags (still accepting new messages)
  // sending threads read, comm-thread modifies
  std::unordered_set<int> sendAccept;

  /// Per tag: number of processes that haven't sent the END-TAG message yet
  // compute thread set during initialization and setting callbacks, callback-thread read and modify.
  std::unordered_map<int, int> recvRemaining;

  /// The type of the BufferIDs
  typedef typename BuffersType::BufferIdType BufferIdType;

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
  std::map<int,std::function<void(uint8_t*, std::size_t, int)> > callbackFunctions;

  /// The MPI Communicator size
  int commSize;

  /// The MPI Communicator rank
  int commRank;

  /// id of the epochs.  aka the periods between barrier/synchronization.  each epoch is associated with a single flush/finish, thus to 1 tag.
  int epoch;
  //int finalEpoch;
};

const int CommunicationLayer::CONTROL_TAG;

} // namespace io
} // namespace bliss

#endif // BLISS_COMMUNICATION_LAYER_HPP
