/**
 * @file    CommunicationLayer.hpp
 * @ingroup index
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   descr
 *
 * Copyright (c) TODO
 *
 * TODO add Licence
 */

#ifndef BLISS_COMMUNICATION_LAYER_HPP
#define BLISS_COMMUNICATION_LAYER_HPP

#include <mpi.h>

class CommunicationLayer
{
public:

  constexpr static int default_tag = 0;

  CommunicationLayer (const MPI_Comm& communicator)
    : comm(communicator)
  {
    // init communicator rank and size
    MPI_Comm_size(comm, &commSize);
    MPI_Comm_rank(comm, &commRank);
  }

  virtual ~CommunicationLayer ();

  /// potentially buffered communication
  /// TODO: how to handle different message types?
  ///    also for lookups and answers (two way comm)
  ///    and end tag messages (should also be definable from the outside)
  /// TODO: what is the `message` per se?
  void sendMessage(const void* data, int count, int dst_rank, int tag=default_tag)
  {
    MPI_Request req;
    MPI_ISend(data, count, MPI_UINT8_T, dst_rank, tag, comm, &req);
  }


  int getCommSize() const
  {
    return commSize;
  }

  int getCommRank() const
  {
    return commRank;
  }

  void flush();

  // This _requires_ some kind of async threading system
  void addReceiveCallback(int sometag, std::function something);

  // active receiving (by polling) for when there is no callback set
  // these must be thread-safe!
  Message receiveAnyOne();
  std::vector<message> receiveAnyAll();

  Message receiveOne(tag);
  std::vector<message> receiveAll(tag);

private:
  /* data */

  /// The MPI Communicator object for this communication layer
  MPI_Comm comm;

  /// The MPI Communicator size
  int commSize;

  /// The MPI Communicator rank
  int commRank;
};

#endif // BLISS_COMMUNICATION_LAYER_HPP
