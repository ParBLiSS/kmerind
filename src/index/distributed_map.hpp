/**
 * @file    distributed_map.hpp
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


template<typename K, typename T, typename CommunicationLayer, typename LocalContainer>
class _distributed_map_base
{
public:
  /// The iterator type of the local container type
  typedef typename LocalContainer::iterator local_iterator;
  /// The constant iterator type of the local container type
  typedef typename LocalContainer::const_iterator local_const_iterator;

  /**
   * @brief Returns an iterator to the first element of the local container.
   *
   * @return Iterator to the first local element.
   */
  local_iterator begin()
  {
    return local_map.begin();
  }

  /**
   * @brief Returns an iterator to the element following the last element of
   * the local container.
   *
   * @return Iterator the the element following the last local element.
   */
  local_iterator end()
  {
    return local_map.end();
  }

  /**
   * @brief Returns an iterator to the first element of the local container.
   *
   * @return Iterator to the first local element.
   */
  local_const_iterator begin() const
  {
    return local_map.begin();
  }

  /**
   * @brief Returns an iterator to the element following the last element of
   * the local container.
   *
   * @return Iterator the the element following the last local element.
   */
  local_const_iterator end() const
  {
    return local_map.end();
  }


protected:
  _distributed_map_base(MPI_Comm mpi_comm, int comm_size, std::function<std::size_t(K)> hashFunction = std::hash<K>())
      : commLayer(mpi_comm, comm_size), comm(mpi_comm), hashFunct(hashFunction)
  {
  }

  virtual ~_distributed_map_base() {}

  // sends key only
  void sendKey(const K& key, const int dstRank, const int tag)
  {
    // cast key into pointer and get byte size
    const uint8_t* msg = reinterpret_cast<const uint8_t*>(&key);
    const std::size_t count = sizeof(key);

    // send the key as a message with the approriate tag
    commLayer.sendMessage(msg, count, dstRank, tag);
  }

  // sends key and value
  void sendPair(const K& key, const T& value, const int dstRank, const int tag)
  {
    // create the pair and call the overloaded function
    std::pair<K, T> keyElement(key, value);
    this->sendPair(keyElement, dstRank, tag);
  }

  // sends key-value pair
  void sendPair(const std::pair<K, T>& keyValue, const int dstRank, const int tag)
  {
    // create key-value pair and serialize as pointer
    const uint8_t* msg = reinterpret_cast<const uint8_t*>(&keyValue);
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

protected:
  /******************
   *  Data members  *
   ******************/

  CommunicationLayer commLayer;
  MPI_Comm comm;
  std::function<std::size_t(K)> hashFunct;

  LocalContainer local_map;

};


template<typename K, typename T, typename CommunicationLayer>
class distributed_multimap : public _distributed_map_base<K,T,CommunicationLayer,std::unordered_multimap<K, T> >
{
public:
  /// The baseclass type
  typedef _distributed_map_base<K, T, CommunicationLayer, std::unordered_multimap<K,T> > _base_class;
  /// The iterator type of the local container type
  typedef typename _base_class::local_iterator local_iterator;
  /// The constant iterator type of the local container type
  typedef typename _base_class::local_const_iterator local_const_iterator;


  static constexpr int INSERT_MPI_TAG = 13;
  static constexpr int LOOKUP_MPI_TAG = 14;
  static constexpr int LOOKUP_ANSWER_MPI_TAG = 15;

  distributed_multimap (MPI_Comm mpi_comm, int comm_size,
        std::function<std::size_t(K)> hashFunction = std::hash<K>())
      : _base_class(mpi_comm, comm_size, hashFunction)
  {
    // add comm layer receive callbacks
    using namespace std::placeholders;
    this->commLayer.addReceiveCallback(INSERT_MPI_TAG, std::bind(&distributed_multimap::receivedCountCallback, this, _1, _2, _3));
    this->commLayer.addReceiveCallback(LOOKUP_MPI_TAG, std::bind(&distributed_multimap::receivedLookupCallback, this, _1, _2, _3));
    this->commLayer.addReceiveCallback(LOOKUP_ANSWER_MPI_TAG, std::bind(&distributed_multimap::receivedLookupAnswerCallback, this, _1, _2, _3));

    // start the threads in the comm layer (if not already running)
    this->commLayer.startThreads();
  }

  virtual ~distributed_multimap() {
    this->commLayer.finishTag(INSERT_MPI_TAG);
    this->commLayer.finishTag(LOOKUP_MPI_TAG);
    this->commLayer.finishTag(LOOKUP_ANSWER_MPI_TAG);
  }

  void remoteInsert(const std::pair<K, T>& keyvalue)
  {
    int targetRank = this->getTargetRank(keyvalue.first);
    this->sendPair(keyvalue, targetRank, INSERT_MPI_TAG);
    has_pending_inserts = true;
  }

  void remoteInsert(const K& key, const T& value)
  {
    int targetRank = this->getTargetRank(key);
    this->sendPair(key, value, targetRank, INSERT_MPI_TAG);
    has_pending_inserts = true;
  }

  template<typename Iterator>
  void populate(const Iterator& begin, const Iterator& end)
  {
    // get the iterator traits
    typedef typename std::iterator_traits<Iterator> traits;
    typedef typename traits::value_type value_type;
    // check for the correct iterator traits
    static_assert(std::is_same<value_type, std::pair<K, T> >::value, "Iterator value_type must be a std::pair of (key,value)");

    // iterate through all elements and insert them
    for (Iterator it = begin; it != end; ++it)
    {
      int targetRank = this->getTargetRank(it->first);
      this->sendPair(*it, targetRank, INSERT_MPI_TAG);
    }
    // flush the insert mesages
    this->commLayer.flush(INSERT_MPI_TAG);
  }

  /// Flushes all buffered elements to be inserted at the target processor.
  /// Blocks till all elements have been received at their destination.
  void flush()
  {
    if (has_pending_inserts)
    {
      this->commLayer.flush(INSERT_MPI_TAG);
      has_pending_inserts = false;
    }
    if (has_pending_lookups)
    {
      this->commLayer.flush(LOOKUP_MPI_TAG);
      this->commLayer.flush(LOOKUP_ANSWER_MPI_TAG);
      has_pending_lookups = false;
    }
  }

  void asyncLookup(const K& key)
  {
    // check that there is a valid callback function
    if (lookupAnswerCallbackFunc == nullptr)
    {
      throw std::runtime_error("ERROR: Callback function not set!");
    }
    const int targetRank = this->getTargetRank(key);
    this->sendKey(key, targetRank, LOOKUP_MPI_TAG);
    has_pending_lookups = true;
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
    for (auto iter = this->local_map.begin(); iter!=this->local_map.end();)
    {
      // get end of the range of identical key
      auto cur_end_iter = this->local_map.equal_range(iter->first)->second;
      std::size_t cur_count = getLocalCount(iter->first);
      if (cur_count < count)
      {
        // remove all entries with this key, this will invalidate the `iter`
        // iterator
        this->local_map.erase(iter, cur_end_iter);
      }
      // advance loop iterator
      iter = cur_end_iter;
    }
  }

  void setLookupAnswerCallback(const std::function<void(std::pair<K, T>&)>& callbackFunction)
  {
    lookupAnswerCallbackFunc = callbackFunction;
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
      this->local_map.insert(elements[i]);
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
      auto range = this->local_map.equal_range(keys[i]);
      for (auto it = range.first; it != range.second; ++it)
      {
        // send the results to the requesting processor
        this->sendPair(*it, fromRank, LOOKUP_ANSWER_MPI_TAG);
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
      lookupAnswerCallbackFunc(elements[i]);
    }
  }


  // implementation depends on wheather we use multimap or map saving counts
  // explicitly
  std::size_t getLocalCount(const K& key)
  {
    return this->local_map.count(key);
  }

  // TODO: - granularity/resolution of histogram (i.e., bin size)
  //       - sampling rather than full histogram (estimation)
  std::vector<int> countHistrogram()
  {
    // determine some granuarity?
    // TODO

    // first determine the maximum count
    uint64_t local_max_count = 0; // use uint64_t for all systems!
    for (auto iter=this->local_map.begin(); iter!=this->local_map.end();
         iter=this->local_map.equal_range(iter->first)->second)
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
    MPI_Allreduce(&max_count, &all_max_count, 1, MPI_INT, MPI_MAX, this->comm);

    // count the counts to create local histogram
    std::vector<int> local_count_hist(all_max_count+1, 0);
    for (auto iter=this->local_map.begin(); iter!=this->local_map.end();
         iter=this->local_map.equal_range(iter->first)->second)
    {
      std::size_t count = getLocalCount(iter->first);
      local_count_hist[count]++;
    }

    // then accumulate accross all processors
    std::vector<int> count_hist(all_max_count+1, 0);
    MPI_Allreduce(&local_count_hist[0], &count_hist[0], all_max_count+1,
                  MPI_INT, MPI_SUM, this->comm);

    return count_hist;
  }

private:
  /* data */
  volatile bool has_pending_inserts = false;
  volatile bool has_pending_lookups = false;
  std::function<void(std::pair<K, T>&)> lookupAnswerCallbackFunc;
};


typedef uint32_t count_t;

template<typename K, typename CommunicationLayer>
class distributed_counting_map
 : public _distributed_map_base<K, count_t, CommunicationLayer,
                                std::unordered_map<K, count_t> >
{
public:
  /// The baseclass type
  typedef _distributed_map_base<K, count_t, CommunicationLayer, std::unordered_map<K,count_t> > _base_class;
  /// The iterator type of the local container type
  typedef typename _base_class::local_iterator local_iterator;
  /// The constant iterator type of the local container type
  typedef typename _base_class::local_const_iterator local_const_iterator;

  // the value type of the (key,value) pairs in the hash table
  typedef count_t T;

  static constexpr int INSERT_MPI_TAG = 13;
  static constexpr int LOOKUP_MPI_TAG = 14;
  static constexpr int LOOKUP_ANSWER_MPI_TAG = 15;

  distributed_counting_map (MPI_Comm mpi_comm, int comm_size,
          std::function<std::size_t(K)> hashFunction = std::hash<K>())
      : _base_class(mpi_comm, comm_size, hashFunction)
  {
    // add comm layer receive callbacks
    using namespace std::placeholders;
    this->commLayer.addReceiveCallback(INSERT_MPI_TAG, std::bind(&distributed_counting_map::receivedCountCallback, this, _1, _2, _3));
    this->commLayer.addReceiveCallback(LOOKUP_MPI_TAG, std::bind(&distributed_counting_map::receivedLookupCallback, this, _1, _2, _3));
    this->commLayer.addReceiveCallback(LOOKUP_ANSWER_MPI_TAG, std::bind(&distributed_counting_map::receivedLookupAnswerCallback, this, _1, _2, _3));

    // start the threads in the comm layer (if not already running)
    this->commLayer.initCommunication();
  }

  virtual ~distributed_counting_map()
  {
    this->commLayer.finishTag(INSERT_MPI_TAG);
    this->commLayer.finishTag(LOOKUP_MPI_TAG);
    this->commLayer.finishTag(LOOKUP_ANSWER_MPI_TAG);
  }

  void remoteInsert(const K& key)
  {
    int targetRank = this->getTargetRank(key);
    this->sendKey(key, targetRank, INSERT_MPI_TAG);
    has_pending_inserts = true;
  }

  template<typename Iterator>
  void populate(const Iterator& begin, const Iterator& end)
  {
    // get the iterator traits
    typedef typename std::iterator_traits<Iterator> traits;
    typedef typename traits::value_type value_type;
    // check for the correct iterator traits
    static_assert(std::is_same<value_type, K>::value,
                  "Iterator value_type must be the same as the key type `K`");

    // iterate through all elements and insert them
    for (Iterator it = begin; it != end; ++it)
    {
      int targetRank = this->getTargetRank(*it);
      this->sendKey(*it, targetRank, INSERT_MPI_TAG);
    }
    this->commLayer.flush(INSERT_MPI_TAG);
  }

  /// Flushes all buffered elements to be inserted at the target processor.
  /// Blocks till all elements have been received at their destination.
  void flush()
  {
    if (has_pending_inserts)
    {
      this->commLayer.flush(INSERT_MPI_TAG);
      has_pending_inserts = false;
    }
    if (has_pending_lookups)
    {
      this->commLayer.flush(LOOKUP_MPI_TAG);
      this->commLayer.flush(LOOKUP_ANSWER_MPI_TAG);
      has_pending_lookups = false;
    }
  }

  void asyncLookup(const K& key)
  {
    // check that there is a valid callback function
    if (lookupAnswerCallbackFunc == nullptr)
    {
      throw std::runtime_error("ERROR: Callback function not set!");
    }

    const int targetRank = this->getTargetRank(key);
    this->sendKey(key, targetRank, LOOKUP_MPI_TAG);
    has_pending_lookups = true;
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
    for (auto iter = this->local_map.begin(); iter!=this->local_map.end();)
    {
      // get end of the range of identical key
      auto cur_end_iter = this->local_map.equal_range(iter->first)->second;
      std::size_t cur_count = iter->second;
      if (cur_count < count)
      {
        // remove all entries with this key, this will invalidate the `iter`
        // iterator
        this->local_map.erase(iter, cur_end_iter);
      }
      // advance loop iterator
      iter = cur_end_iter;
    }
  }

  void setLookupAnswerCallback(const std::function<void(std::pair<K, T>&)>& callbackFunction)
  {
    lookupAnswerCallbackFunc = callbackFunction;
  }



protected:

  // for counting index
  void receivedCountCallback(uint8_t* msg, std::size_t count, int fromRank)
  {
    // deserialize
    K* keys = reinterpret_cast<K*>(msg);
    int key_count = count / sizeof(K);

    // insert all elements into the hash table
    for (int i = 0; i < key_count; ++i)
    {
      this->local_map[keys[i]]++;
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
      auto range = this->local_map.equal_range(keys[i]);
      for (auto it = range.first; it != range.second; ++it)
      {
        // send the results to the requesting processor
        this->sendPair(*it, fromRank, LOOKUP_ANSWER_MPI_TAG);
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
      // call the external callback function
      lookupAnswerCallbackFunc(elements[i]);
    }
  }

  // TODO: - granularity/resolution of histogram (i.e., bin size)
  //       - sampling rather than full histogram (estimation)
  std::vector<int> countHistrogram()
  {
    // determine some granuarity?
    // TODO

    // first determine the maximum count
    int max_count = 0; // use uint64_t for all systems!
    for (auto iter=this->local_map.begin(); iter!=this->local_map.end(); ++iter)
    {
      max_count = std::max<int>(max_count, iter->second);
    }

    // get max accross all processors
    int all_max_count;
    MPI_Allreduce(&max_count, &all_max_count, 1, MPI_INT, MPI_MAX, this->comm);

    // count the counts to create local histogram
    std::vector<int> local_count_hist(all_max_count+1, 0);
    for (auto iter=this->local_map.begin(); iter!=this->local_map.end(); ++iter)
    {
      local_count_hist[iter->second]++;
    }

    // then accumulate accross all processors
    std::vector<int> count_hist(all_max_count+1, 0);
    MPI_Allreduce(&local_count_hist[0], &count_hist[0], all_max_count+1,
                  MPI_INT, MPI_SUM, this->comm);

    return count_hist;
  }

private:
  /* data */
  volatile bool has_pending_inserts = false;
  volatile bool has_pending_lookups = false;
  std::function<void(std::pair<K, T>&)> lookupAnswerCallbackFunc;
};

#endif // BLISS_DISTRIBUTED_MAP_HPP
