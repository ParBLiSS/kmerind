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
#include <unordered_map> // local storage hash table

template<typename K, typename T, typename CommunicationLayer, typename DistrFunction, typename LocalContainer=std::unordered_multimap<K, T> >
class DistributedIndex
{
public:

  static constexpr int INSERT_MPI_TAG = 13;
  static constexpr int LOOKUP_MPI_TAG = 14;
  static constexpr int LOOKUP_ANSWER_MPI_TAG = 15;

  DistributedIndex (CommunicationLayer& commLayer)
      : commLayer(commLayer)
  {
    // TODO: add callback function for commLayer receive
  }

  virtual ~DistributedIndex () {}

  void remoteInsert(const K& key, const T& value)
  {
    int targetRank = getTargetRank(key);
    sendPair(key, value, targetRank, INSERT_MPI_TAG);
  }

  /// Flushes all buffered elements to be inserted at the target processor.
  /// Blocks till all elements have been received at their destination.
  void flush()
  {
    commLayer.flush();
  }

  void asyncLookup(const K& key)
  {
    const int targetRank = getTargetRank(key);
    sendKey(key, targetRank, LOOKUP_MPI_TAG);
  }

  void syncLookup(const K& key)
  {
    const int targetRank = getTargetRank(key);
    sendKey(key, targetRank, LOOKUP_MPI_TAG);
    // TODO: wait for answer somehow
  }

protected:
  // for positional index
  void receivedInsertCallback(uint8_t* msg, std::size_t count, int fromRank)
  {
    // deserialize
    std::pair<K, T>* elements = reinterpret_cast<std::pair<K, T>*>(msg);
    int element_count = count / sizeof(std::pair<K, T>);

    // insert all elements into the hash table
    for (int i = 0; i < element_count; ++i)
    {
      hashTable.insert(elements[i]);
    }
  }

  // for counting index
  void receivedCountCallback(uint8_t* msg, std::size_t count, int fromRank)
  {
    // deserialize
    K* keys = reinterpret_cast<K*>(msg);
    int key_count = count / sizeof(K);

    // insert all elements into the hash table
    for (int i = 0; i < key_count; ++i)
    {
      hashTable[keys[i]]++;
    }
  }

  void receivedLookupCallback(uint8_t* msg, std::size_t count, int fromRank)
  {
    // deserialize
    K* keys = reinterpret_cast<K*>(msg);
    int key_count = count / sizeof(K);

    // for all received requests, send the value from the lookup
    for (int i = 0; i < key_count; ++i)
    {
      // check if exists and then send
      auto range = hashTable.equal_range(keys[i]);
      for (auto it = range.first; it != range.second; ++it)
      {
        // send the results to the requesting processor
        sendPair(*it, fromRank, LOOKUP_ANSWER_MPI_TAG);
      }
    }
  }

  void receivedLookupAnswerCallback(uint8_t* msg, std::size_t count, int fromRank)
  {
    // deserialize
    std::pair<K, T>* elements = reinterpret_cast<std::pair<K, T>*>(msg);
    int element_count = count / sizeof(std::pair<K, T>);

    // insert all elements into the hash table
    for (int i = 0; i < element_count; ++i)
    {
      // TODO: call some external callback function, or queue this
    }
  }

  // sends key only
  void sendKey(const K& key, const int dstRank, const int tag)
  {
    // cast key into pointer and get byte size
    const uint8_t* msg = reinterpret_cast<uint8_t*>(&key);
    const std::size_t count = sizeof(key);

    // send the key as a message with the approriate tag
    commLayer.sendMessage(msg, count, dstRank, tag);
  }

  // sends key and value
  void sendPair(const K& key, const T& value, const int dstRank, const int tag)
  {
    // create the pair and call the overloaded function
    std::pair<K, T> keyElement(key, value);
    sendPair(keyElement, dstRank, tag);
  }

  // sends key-value pair
  void sendPair(const std::pair<K, T>& keyValue, const int dstRank, const int tag)
  {
    // create key-value pair and serialize as pointer
    const uint8_t* msg = reinterpret_cast<uint8_t*>(&keyValue);
    const std::size_t count = sizeof(keyValue);

    // send the message
    commLayer.sendMessage(msg, count, dstRank, tag);
  }

  // returns the target rank for a given key (uses the distribution function)
  int getTargetRank(const K& key)
  {
    // get the target rank for the processor
    int size = commLayer.getCommSize();
    return hashFunct(key) % size;
  }

  CommunicationLayer commLayer;

  DistrFunction hashFunct;

  LocalContainer hashTable;

private:
  /* data */
};

#endif // BLISS_DISTRIBUTED_INDEX_HPP
