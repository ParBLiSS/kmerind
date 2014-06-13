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

class CommunicationLayer
{
public:

  constexpr static int default_tag = 0;

  CommunicationLayer (const MPI_Comm& communicator);
  virtual ~CommunicationLayer ();

  /// potentially buffered communication
  /// TODO: how to handle different message types?
  ///    also for lookups and answers (two way comm)
  ///    and end tag messages (should also be definable from the outside)
  /// TODO: what is the `message` per se?
  template<typename T>
  void sendMessage(int dst_rank, const T* data, std::size_t data_size, int tag=default_tag)
  {
    
  }

  template<typename T>
  void sendMessage(int dst_rank, const T& data, int tag=default_tag)
  {
    int count = sizeof(T);
    // how to deserialize the data?
  }



  void flush();

  // This _requires_ some kind of async threading system
  void addReceiveCallback(sometag, std::function something);

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
};

#endif // BLISS_COMMUNICATION_LAYER_HPP
