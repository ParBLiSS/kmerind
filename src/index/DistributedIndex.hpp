/**
 * @file    DistributedIndex.hpp
 * @ingroup index
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   descr
 *
 * Copyright (c) TODO
 *
 * TODO add Licence
 */

#ifndef BLISS_DISTRIBUTED_INDEX_HPP
#define BLISS_DISTRIBUTED_INDEX_HPP

#include <utility> // for std::pair

template<typename T, typename K, typename CommunicationLayer, typename HashFunction>
class DistributedIndex
{
public:
  DistributedIndex (arguments);
  virtual ~DistributedIndex ();

  void remoteInsert(const K& key, const T& element)
  {
    commLayer.sendMessage(/* TODO */);
  }

  /// Flushes all buffered elements to be inserted at the target processor.
  /// Blocks till all elements have been received at their destination.
  void flush();

  // TODO: Actually could be much more, if there are multiple occurances
  std::pair<K, T> lookup(const K& key);


protected:
  void localInsert(const K& key, const T& element);

  void receivedElementCallback(const K& key, const T& element);

  CommunicationLayer commLayer;

private:
  /* data */
};

#endif // BLISS_DISTRIBUTED_INDEX_HPP
