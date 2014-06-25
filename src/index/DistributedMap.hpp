/**
 * @file    DistributedMap.hpp
 * @ingroup index
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   descr
 *
 * Copyright (c) TODO
 *
 * TODO add Licence
 */

#ifndef BLISS_DISTRIBUTED_MAP_HPP
#define BLISS_DISTRIBUTED_MAP_HPP

#include <utility> // for std::pair
#include <unordered_map> // local storage hash table
#include <vector>
#include <functional> // for std::function and std::hash
#include <limits>
#include <stdexcept>
#include <algorithm> // for std::max

// include MPI
#include <mpi.h>

// TODO LIST:
//  - [ ] split up distributed index into multimap and counting map
//  - [ ] populate(Iterator)
//  - [ ] WAIT FOR COMM_LAYER: flush()
//  - [ ] expose local iterators
//  - [X] finish the count Histrogram
//  - [x] filter() function -> (maybe with broadcast??)
//  - [ ] take apart commlayer into one-way and two-way comm with proper destruction

template<typename K, typename T, typename CommunicationLayer, typename LocalContainer=std::unordered_multimap<K, T> >
class DistributedMap
{
public:

  static constexpr int INSERT_MPI_TAG = 13;
  static constexpr int LOOKUP_MPI_TAG = 14;
  static constexpr int LOOKUP_ANSWER_MPI_TAG = 15;

  DistributedMap (CommunicationLayer& commLayer, MPI_Comm mpi_comm, std::function<std::size_t(K)> hashFunction = std::hash<K>())
      : commLayer(commLayer), comm(mpi_comm), hashFunct(hashFunction)
  {
    // TODO: add callback function for commLayer receive
  }

  virtual ~DistributedMap () {}

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


  /**
   * @brief Removes all (key,value) pairs with a key count of less than `count`.
   *
   * This function has to be called on each MPI process with the same parameter
   * value. This function only operates locally, there is no communication
   * involved.
   *
   * @param count   The key count threshold. Everything lower than this will be
   *                removed.
   */
  void filter(const std::size_t count)
  {
    // iterate through all keys
    for (auto iter = hashTable.begin(); iter!=hashTable.end();)
    {
      // get end of the range of identical key
      auto cur_end_iter = hashTable.equal_range(iter->first)->second;
      std::size_t cur_count = getLocalCount(iter->first);
      if (cur_count < count)
      {
        // remove all entries with this key, this will invalidate the `iter`
        // iterator
        hashTable.erase(iter, cur_end_iter);
      }
      // advance loop iterator
      iter = cur_end_iter;
    }
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

  // implementation depends on wheather we use multimap or map saving counts
  // explicitly
  std::size_t getLocalCount(const K& key)
  {
    return hashTable.count(key);
  }

  // TODO: - granularity/resolution of histogram (i.e., bin size)
  //       - sampling rather than full histogram (estimation)
  std::vector<int> countHistrogram()
  {
    // determine some granuarity?
    // TODO

    // first determine the maximum count
    uint64_t local_max_count = 0; // use uint64_t for all systems!
    for (auto iter=hashTable.begin(); iter!=hashTable.end();
         iter=hashTable.equal_range(iter->first)->second)
    {
      std::size_t count = getLocalCount(iter->first);
      local_max_count = std::max<uint64_t>(local_max_count, count);
    }

    // cast max count to int and check that it doesn't overflow
    if (local_max_count >= std::numeric_limits<int>::max())
    {
      throw std::range_error("Histrogram of counts: maximum count exceeds integer range");
    }
    int max_count = static_cast<int>(local_max_count);

    // get max accross all processors
    int all_max_count;
    MPI_Allreduce(&max_count, &all_max_count, 1, MPI_INT, MPI_MAX, comm);

    // count the counts to create local histogram
    std::vector<int> local_count_hist(all_max_count+1, 0);
    for (auto iter=hashTable.begin(); iter!=hashTable.end();
         iter=hashTable.equal_range(iter->first)->second)
    {
      std::size_t count = getLocalCount(iter->first);
      local_count_hist[count]++;
    }

    // then accumulate accross all processors
    std::vector<int> count_hist(all_max_count+1, 0);
    MPI_Allreduce(&local_count_hist[0], &count_hist[0], all_max_count+1, MPI_UINT64_T, MPI_SUM, comm);

    return count_hist;
  }

  // returns the target rank for a given key (uses the distribution function)
  int getTargetRank(const K& key)
  {
    // get the target rank for the processor
    int size = commLayer.getCommSize();
    return hashFunct(key) % size;
  }

  CommunicationLayer commLayer;

  MPI_Comm comm;

  std::function<std::size_t(K)> hashFunct;

  LocalContainer hashTable;

private:
  /* data */
};

#endif // BLISS_DISTRIBUTED_MAP_HPP
