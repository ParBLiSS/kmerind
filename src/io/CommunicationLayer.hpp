/**
 * @file    CommunicationLayer.hpp
 * @ingroup index
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @author  Tony Pan <tpan@gatech.edu>
 * @brief   descr
 *
 * Copyright (c) TODO
 *
 * TODO add Licence
 */

#ifndef BLISS_COMMUNICATION_LAYER_HPP
#define BLISS_COMMUNICATION_LAYER_HPP

#include <mpi.h>
//#include <omp.h>

// C stdlib includes
#include <assert.h>

// STL includes
#include <vector>
#include <queue>
#include <mutex>

// system includes
// TODO: make this system indepenedent!?
#include <unistd.h> // for usleep()

// BLISS includes
#include <concurrent/threadsafe_queue.hpp>
#include <concurrent/concurrent.hpp>
#include <io/message_buffers.hpp>

typedef bliss::io::MessageBuffers<bliss::concurrent::THREAD_SAFE> BuffersType;


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
  typename bliss::io::MessageBuffers<bliss::concurrent::THREAD_SAFE>::BufferIdType bufferId;
  int tag;
  int dst;

  SendQueueElement(bliss::io::MessageBuffers<bliss::concurrent::THREAD_SAFE>::BufferIdType _id, int _tag, int _dst)
    : bufferId(_id), tag(_tag), dst(_dst) {}
  SendQueueElement() = delete;
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

  /// Outbound message structure, Multiple-Producer-Single-Consumer queue
  /// consumed by the internal MPI-comm thread
  bliss::concurrent::ThreadSafeQueue<SendQueueElement> sendQueue;

  // Inbound message structure, Single-Producer-Multiple-Consumer queue
  // produced by the internal MPI-comm thread
  bliss::concurrent::ThreadSafeQueue<ReceivedMessage> recvQueue;

  // outbound temporary data buffer type for the producer threads.  ThreadSafe version for now.
  typedef bliss::io::MessageBuffers<bliss::concurrent::THREAD_SAFE> BuffersType;
  // outbound temporary data buffers.  one per tag.
  std::unordered_map<int, BuffersType> buffers;
  typedef typename bliss::io::MessageBuffers<bliss::concurrent::THREAD_SAFE>::BufferIdType BufferIdType;

  mutable std::mutex mutex;

public:

  CommunicationLayer (const MPI_Comm& communicator, const int comm_size)
    : sendQueue(2 * omp_get_num_threads()), recvQueue(2 * comm_size),
      comm(communicator)
  {
    // init communicator rank and size
    MPI_Comm_size(comm, &commSize);
    assert(comm_size == commSize);
    MPI_Comm_rank(comm, &commRank);
  }

  virtual ~CommunicationLayer ();


  void sendMessage(const void* data, std::size_t count, int dst_rank, int tag=default_tag)
  {
    /// if there isn't already a tag listed, add the MessageBuffers for that tag.

    // multiple threads may call this.
    std::unique_lock<std::mutex> lock(mutex);
    if (buffers.find(tag) == buffers.end())
      buffers[tag] = std::move(BuffersType(commSize, 8192));
    lock.unlock();

    /// try to append the new data - repeat until successful.
    /// along the way, if a full buffer's id is returned, queue it for send.
    BufferIdType fullId = -1;
    while (!buffers[tag].append(data, count, dst_rank, fullId)) {
      // repeat until success;
      if (fullId != -1) {
        // have a full buffer - put in send queue.
        SendQueueElement v(fullId, tag, dst_rank);
        while (!sendQueue.tryPush(std::move(v))) {
          usleep(20);
        }
      }

      usleep(20);
    }

    /// that's it.
  }

  void flush(int tag)
  {
    // flush out all the send buffers matching a particular tag.
    int i = 0;
    for (auto id : buffers[tag].getActiveIds()) {
      if (id != -1) {
        SendQueueElement v(id, tag, i);
        while (!sendQueue.tryPush(std::move(v))) {
          usleep(20);
        }
      }
    }
  }


  //
  void commThread()
  {
    // TODO: implement termination condition
    while(true)
    {
      // first clean up finished operations
      finishReceives();
      finishSends();

      // start pending receives
      tryStartReceive();
      // start pending sends
      tryStartSend();
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

      /*
      if (tag == MPIBufferType::END_TAG) {
        // end of messaging.
        MPI_Recv(nullptr, received_count, MPI_UNSIGNED_CHAR, src, tag, comm, MPI_STATUS_IGNORE);
#pragma omp atomic
        --mpi_senders;
        printf("RECV rank %d receiving END signal %d from %d, num sanders remaining is %d\n", rank, received, src, mpi_senders);
#pragma omp flush(mpi_senders)
      } else {
    */
        // receive the message data bytes into a vector of bytes
        uint8_t* msg_data = new uint8_t[received_count];
        MPI_Request req;
        MPI_Irecv(msg_data, received_count, MPI_UNSIGNED_CHAR, src, tag, comm, &req);

        // insert into the in-progress queue
        std::tuple<MPI_Request, ReceivedMessage> msg(req, ReceivedMessage(msg_data, received_count, tag, src));
        recvInProgress.push(std::move(msg));

     // }
    }
  }

  void tryStartSend()
  {
    // try to get the send element
    SendQueueElement se;
    if (sendQueue.tryPop(se)) {
      auto data = buffers[se.tag].getBackBuffer(se.bufferId).getData();
      auto count = buffers[se.tag].getBackBuffer(se.bufferId).getSize();


      if (se.dst == commRank) {
        // local, directly handle by creating an output object and directly
        // insert into the recv queue
        uint8_t* array = new uint8_t[count];
        memcpy(array, data, count);
        ReceivedMessage msg(array, count, se.tag, se.dst);
        while(!recvQueue.tryPush(std::move(msg))) {
          usleep(50);
        }

        // finished inserting.  release the buffer
        buffers[se.tag].releaseBuffer(se.bufferId);
      } else {
        MPI_Request req;
        MPI_Isend(data, count, MPI_UINT8_T, se.dst, se.tag, comm, &req);
        sendInProgress.push(std::move(std::pair<MPI_Request, SendQueueElement>(req, se)));
      }
    }
  }

  void finishSends()
  {
    int finished = 0;
    while(!sendInProgress.empty())
    {
      std::pair<MPI_Request, SendQueueElement>& front = sendInProgress.front();
      MPI_Test(&front.first, &finished, MPI_STATUS_IGNORE);

      if (finished)
      {
        // cleanup, i.e., release the buffer back into the pool
        buffers[front.second.tag].releaseBuffer(front.second.bufferId);
        sendInProgress.pop();
      }
      else
      {
        break;
      }
    }
  }

  // Not thread safe
  void finishReceives()
  {
    int finished = 0;
    MPI_Status status;
    while(!recvInProgress.empty())
    {
      std::pair<MPI_Request, ReceivedMessage>& front = recvInProgress.front();
      MPI_Test(&std::get<0>(front), &finished, &status);

      if (finished)
      {
        assert(front.second.tag == status.MPI_TAG);
        assert(front.second.src == status.MPI_SOURCE);

        // add the received messages into the recvQueue
//        ReceivedMessage msg(std::get<1>(front), std::get<2>(front),
//                            status.MPI_TAG, status.MPI_SOURCE);
        // TODO: (maybe) one queue per tag?? Where does this
        //       sorting/categorizing happen?
        while (!recvQueue.tryPush(std::move(front.second))) {
          usleep(50);
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

  int getCommSize() const
  {
    return commSize;
  }

  int getCommRank() const
  {
    return commRank;
  }



  void callbackThread()
  {
    // TODO: add termination condition
    while(true)
    {
      // get next element from the queue, wait if none is available
      ReceivedMessage msg;
      recvQueue.waitAndPop(msg);

      // TODO: check if the tag exists as callback function

      // call the matching callback function
      (callbackFunctions[msg.tag])(msg.data, msg.count, msg.src);
      delete [] msg.data;
    }
  }

  // adding the callback function with signature:
  // void(uint8_t* msg, std::size_t count, int fromRank)
  void addReceiveCallback(int tag, std::function<void(uint8_t*, std::size_t, int)> callbackFunction)
  {
    if (callbackFunctions.empty())
    {
      // this is the first registered callback, thus spawn the callback
      // executer thread
      // TODO
    }
    // add the callback function to a lookup table
    callbackFunctions[tag] = callbackFunction;
  }

  /*
  // active receiving (by polling) for when there is no callback set
  // these must be thread-safe!
  Message receiveAnyOne();
  std::vector<message> receiveAnyAll();

  Message receiveOne(tag);
  std::vector<message> receiveAll(tag);
*/

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
