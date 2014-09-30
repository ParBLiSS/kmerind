/**
 * @file    distributed_map.hpp
 * @ingroup index
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   Implements the distributed_multimap and distributed_counting_map
 *          data structures.
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
#include <atomic>

// include MPI
#include <mpi.h>

namespace bliss
{
namespace index
{

/// Type for the counting map
typedef uint32_t count_t;

/**
 * @brief A shared base class for the distributed_multimap and
 *        distributed_counting_map implementations.
 *
 * @tparam K                    The key type.
 * @tparam T                    The value type.
 * @tparam CommunicationLayer   The CommunicationLayer class type.
 * @tparam LocalContainer       The local container type (unordered_map or
 *                              unodered_multimap)
 */
template<typename K, typename T, typename CommunicationLayer,
         typename LocalContainer>
class _distributed_map_base
{
public:
  /// The iterator type of the local container type
  typedef typename LocalContainer::iterator local_iterator;
  /// The constant iterator type of the local container type
  typedef typename LocalContainer::const_iterator local_const_iterator;

  typedef typename LocalContainer::hasher HasherType;

  /// MPI message tag for inserts
  static constexpr int INSERT_MPI_TAG = 13;
  /// MPI message tag for lookup queries
  static constexpr int LOOKUP_MPI_TAG = 14;
  /// MPI message tag for answers to lookup queries
  static constexpr int LOOKUP_ANSWER_MPI_TAG = 15;

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

  /**
   * @brief   Returns the local map's size.  primarily for debugging.
   * @return  size of the local map.
   */
  const size_t local_size() const
  {
    return local_map.size();
  }

  /**
   * @brief Flushes all pending operations.
   *
   * Since all insert and lookup operations are executed asynchronosly,
   * calling this function ensures that all pending operations have been
   * executed, including that all pending lookups have returned an answer.
   */
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

  /**
   * @brief Posts an asynchronous lookup for the given key.
   *
   * This function returns immediately, the answer to the lookup will be
   * returned by the callback function set by `setLookupAnswerCallback()`.
   *
   * @param key The key to look-up.
   */
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
   * @brief Sets the callback function for received answers to asynchronously
   *        posted lookups.
   *
   * @param callbackFunction    The function to call with all received answers
   *                            to lookups.
   */
  void setLookupAnswerCallback(const std::function<void(std::pair<K, T>&)>& callbackFunction)
  {
    lookupAnswerCallbackFunc = callbackFunction;
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
      std::size_t cur_count = getLocalCount(local_map, iter);
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

  /**
   * @brief Returns a histrogram of counts for all elements in the distributed
   *        hash table.
   *
   * The valid histrogram is returned on ALL MPI processes. This function has
   * to be called by ALL MPI processes of the given MPI communicator.
   *
   * TODO: [ ] granularity/resolution of histogram (i.e., bin sizes as param)
   *       [ ] sampling rather than full histogram (estimation)
   *
   * @return The histrogram of counts as std::vector<int>.
   */
  std::vector<int> countHistrogram()
  {
    // first determine the maximum count
    uint64_t local_max_count = 0; // use uint64_t for all systems!
    for (auto iter=this->local_map.begin(); iter!=this->local_map.end();
         iter=this->local_map.equal_range(iter->first)->second)
    {
      std::size_t count = getLocalCount(local_map, *iter);
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
      std::size_t count = getLocalCount(local_map, *iter);
      local_count_hist[count]++;
    }

    // then accumulate accross all processors
    std::vector<int> count_hist(all_max_count+1, 0);
    MPI_Allreduce(&local_count_hist[0], &count_hist[0], all_max_count+1,
                  MPI_INT, MPI_SUM, this->comm);

    return count_hist;
  }

protected:
  /**
   * @brief Constructs the shared base class.
   *
   * @param mpi_comm        The MPI communicator to pass onto the communication
   *                        layer.
   * @param comm_size       The size of the MPI communicator, needed for
   *                        initialization.
   * @param hashFunction    The hash function to use for distributing elements
   *                        accross processors. Defaults to std::hash<K>().
   */
  _distributed_map_base(MPI_Comm mpi_comm, int comm_size, HasherType hashFunction = HasherType())
      : commLayer(mpi_comm, comm_size), comm(mpi_comm), hashFunct(hashFunction),
        has_pending_inserts(false), has_pending_lookups(false)
  {
  }

  /**
   * @brief
   */
  virtual ~_distributed_map_base() {}

  /**
   * @brief Helper function to send the given key to the given rank.
   *
   * @param key     The key to send.
   * @param dstRank The destination rank.
   * @param tag     The message tag used for sending.
   */
  void sendKey(const K& key, const int dstRank, const int tag)
  {
    // cast key into pointer and get byte size
    const uint8_t* msg = reinterpret_cast<const uint8_t*>(&key);
    const std::size_t count = sizeof(key);

    // send the key as a message with the approriate tag
    commLayer.sendMessage(msg, count, dstRank, tag);
  }

  /**
   * @brief Helper function to send the given (key,value) pair to the given rank.
   *
   * @param key     The key to send.
   * @param value   The value to send.
   * @param dstRank The destination rank.
   * @param tag     The message tag used for sending.
   */
  void sendPair(const K& key, const T& value, const int dstRank, const int tag)
  {
    // create the pair and call the overloaded function
    std::pair<K, T> keyElement(key, value);
    this->sendPair(keyElement, dstRank, tag);
  }

  /**
   * @brief Helper function to send the given (key,value) pair to the given rank.
   *
   * @param keyValue    The std::pair of (key,value) to send.
   * @param dstRank     The destination rank.
   * @param tag         The message tag used for sending.
   */
  void sendPair(const std::pair<K, T>& keyValue, const int dstRank, const int tag)
  {
    // create key-value pair and serialize as pointer
    const uint8_t* msg = reinterpret_cast<const uint8_t*>(&keyValue);
    const std::size_t count = sizeof(keyValue);

    // send the message
    commLayer.sendMessage(msg, count, dstRank, tag);
  }

  // returns the target rank for a given key (uses the distribution function)
  /**
   * @brief Returns the target rank for a given key.
   *
   * Uses the given hash function to determine where an element with the given
   * key has to be send to or retrieved from.
   *
   * @param key     The key for which to determine the target rank.
   *
   * @return        The target rank of the processor to which the given key
   *                belongs.
   */
  int getTargetRank(const K& key)
  {
    // get the target rank for the processor
    int size = commLayer.getCommSize();
    return hashFunct(key) % size;
  }

  /**
   * @brief Returns the local count of a given key.
   *
   * @param localMap    The local map data structure.
   * @param item        The item to lookup the count for.
   *
   * @return    The count of elements with the given item's key.
   */
  inline std::size_t getLocalCount(const std::unordered_multimap<K, T>& localMap, const std::pair<K, T>& item)
  {
    // use the maps count() function
    return localMap.count(item.first);
  }

  /**
   * @brief Returns the local count of a given key.
   *
   * @param localMap    The local map data structure.
   * @param item        The item to lookup the count for.
   *
   * @return    The count of elements with the given item's key.
   */
  inline std::size_t getLocalCount(const std::unordered_map<K, T>& localMap, const std::pair<K, T>& item)
  {
    // this is a counting map: i.e. the value of (key,value) is the count
    assert(localMap.find(item.first) != localMap.end());
    return static_cast<std::size_t>(item.second);
  }


  /**
   * @brief Callback function for received lookup messages.
   *
   * Performs the actual lookup in the local data structure and sends a message
   * as reply.
   *
   * @param msg         The message data received.
   * @param count       The number of bytes received.
   * @param fromRank    The source rank of the message.
   */
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

  /**
   * @brief Callback function for received lookup answers.
   *
   * Deserializes the answer to lookup queries and calls the
   * given callback function.
   *
   * @param msg         The message data received.
   * @param count       The number of bytes received.
   * @param fromRank    The source rank of the message.
   */
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

protected:
  /******************
   *  Data members  *
   ******************/

  /// The async communication layer abstraction, used for sending and receiving
  /// messages
  CommunicationLayer commLayer;

  /// The MPI communicator used for distributing the input
  MPI_Comm comm;

  /// The hash function used for distributing input accross processors
  HasherType hashFunct;

  /// The local container data-structure. For most cases this is either
  /// a std::unordered_map or std::unordered_multimap
  LocalContainer local_map;

  /// Whether there are pending insert operations
  std::atomic<bool> has_pending_inserts;
  /// Whether there are pending lookup operations
  std::atomic<bool> has_pending_lookups;

  /// The async callback function, which is called when an answer to a
  /// lookup query is received
  std::function<void(std::pair<K, T>&)> lookupAnswerCallbackFunc;
};


/**
 * @brief   A distributed, asynchronous multimap.
 *
 * This can hold multiple elements with identical key.
 *
 * @tparam K                    The key type.
 * @tparam T                    The value type.
 * @tparam CommunicationLayer   The CommunicationLayer class type.
 */
template<typename K, typename T, typename CommunicationLayer, typename HASHER = std::hash<K> >
class distributed_multimap
  : public _distributed_map_base<K, T,
              CommunicationLayer, std::unordered_multimap<K, T, HASHER> >
{
public:
  /// The baseclass type
  typedef _distributed_map_base<K, T, CommunicationLayer,
                                std::unordered_multimap<K, T, HASHER> > _base_class;
  /// The iterator type of the local container type
  typedef typename _base_class::local_iterator local_iterator;
  /// The constant iterator type of the local container type
  typedef typename _base_class::local_const_iterator local_const_iterator;
  /// The value type of the map
  typedef std::pair<K, T> value_type;


  /**
   * @brief Constructs the distributed multimap
   *
   * @param mpi_comm        The MPI communicator to pass onto the communication
   *                        layer.
   * @param comm_size       The size of the MPI communicator, needed for
   *                        initialization.
   * @param hashFunction    The hash function to use for distributing elements
   *                        accross processors. Defaults to std::hash<K>().
   */
  distributed_multimap (MPI_Comm mpi_comm, int comm_size,
        HASHER hashFunction = HASHER())
      : _base_class(mpi_comm, comm_size, hashFunction)
  {
    // add comm layer receive callbacks
    using namespace std::placeholders;
    this->commLayer.addReceiveCallback(_base_class::INSERT_MPI_TAG,
        std::bind(&distributed_multimap::receivedInsertCallback,
                  this, _1, _2, _3));
    this->commLayer.addReceiveCallback(_base_class::LOOKUP_MPI_TAG,
        std::bind(&distributed_multimap::receivedLookupCallback,
                  this, _1, _2, _3));
    this->commLayer.addReceiveCallback(_base_class::LOOKUP_ANSWER_MPI_TAG,
        std::bind(&distributed_multimap::receivedLookupAnswerCallback,
                  this, _1, _2, _3));

    // start the threads in the comm layer (if not already running)
    this->commLayer.initCommunication();
  }

  /**
   * @brief Destructor.
   */
  virtual ~distributed_multimap() {
    // finish all three tags in order
//    this->commLayer.finishTag(_base_class::INSERT_MPI_TAG);
//    this->commLayer.finishTag(_base_class::LOOKUP_MPI_TAG);
//    this->commLayer.finishTag(_base_class::LOOKUP_ANSWER_MPI_TAG);
	  this->commLayer.finishCommunication();
  }

  /**
   * @brief Inserts the given (key,value) pair into this distributed map.
   *
   * The pair might be inserted on a remote processor based on the hashing
   * function.
   *
   * @param keyvalue    The (key,value) pair (std::pair) to insert.
   */
  void insert(const std::pair<K, T>& keyvalue)
  {
    int targetRank = this->getTargetRank(keyvalue.first);
    this->sendPair(keyvalue, targetRank, _base_class::INSERT_MPI_TAG);
    this->has_pending_inserts = true;
  }

  /**
   * @brief Inserts the given (key,value) pair into this distributed map.
   *
   * The pair might be inserted on a remote processor based on the hashing
   * function.
   *
   * @param key     The key to insert.
   * @param value   The value to insert.
   */
  void insert(const K& key, const T& value)
  {
    int targetRank = this->getTargetRank(key);
    this->sendPair(key, value, targetRank, _base_class::INSERT_MPI_TAG);
    this->has_pending_inserts = true;
  }

  /**
   * @brief   construct a new element in place, with same semantic as std::unordered_map::emplace.
   *
   * construct new element to insert in place.  using variadic template parameters allows
   * both multimap and counting map to provide the same interface.
   * Does not return iterator and bool like unordered_map.
   *
   * @param args  variable number of arguments.
   */
  template <class... Args>
  void emplace(Args&&... args) {
      insert(args...);
  }

  /**
   * @brief Populates the distributed map by bulk-inserting elements from the
   *        given iterators.
   *
   * This function blocks until all elements haven been inserted from all
   * processors.
   *
   * @tparam Iterator   An input iterator with value_type std::pair<K,T>.
   * @param begin       Iterator to the first element which will be inserted.
   * @param end         Iterator to one element past the last element to insert.
   */
  template<typename Iterator>
  void populate(const Iterator& begin, const Iterator& end)
  {
    // get the iterator traits
    typedef typename std::iterator_traits<Iterator> traits;
    typedef typename traits::value_type value_type;
    // check for the correct iterator traits
    static_assert(std::is_same<value_type, std::pair<K, T> >::value,
                  "Iterator value_type must be a std::pair of (key,value)");

    // iterate through all elements and insert them
    for (Iterator it = begin; it != end; ++it)
    {
      int targetRank = this->getTargetRank(it->first);
      this->sendPair(*it, targetRank, _base_class::INSERT_MPI_TAG);
    }
    // flush the insert mesages
    this->commLayer.flush(_base_class::INSERT_MPI_TAG);
  }

protected:
  /**
   * @brief Callback function for received inserts.
   *
   * Deserializes the message data and performs the local insertions.
   *
   * @param msg         The message data received.
   * @param count       The number of bytes received.
   * @param fromRank    The source rank of the message.
   */
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

};


/**
 * @brief   A distributed, asynchronous map which counts how many times an
 *          element with identical key is inserted.
 *
 * This is a map with (key->count). Each time an element with key `k` is
 * inserted, the count of that element is increased: `map[k]++`.
 *
 * @tparam K                    The key type.
 * @tparam CommunicationLayer   The CommunicationLayer class type.
 */
template<typename K, typename CommunicationLayer, typename HASHER = std::hash<K>>
class distributed_counting_map
 : public _distributed_map_base<K, count_t, CommunicationLayer,
                                std::unordered_map<K, count_t, HASHER> >
{
public:
  /// The baseclass type
  typedef _distributed_map_base<K, count_t, CommunicationLayer,
                                std::unordered_map<K,count_t, HASHER> > _base_class;
  /// The iterator type of the local container type
  typedef typename _base_class::local_iterator local_iterator;
  /// The constant iterator type of the local container type
  typedef typename _base_class::local_const_iterator local_const_iterator;
  /// The value type of the map
  typedef std::pair<K, count_t> value_type;

  // the value type of the (key,value) pairs in the hash table
  typedef count_t T;

  /**
   * @brief Constructs the distributed counting map.
   *
   * @param mpi_comm        The MPI communicator to pass onto the communication
   *                        layer.
   * @param comm_size       The size of the MPI communicator, needed for
   *                        initialization.
   * @param hashFunction    The hash function to use for distributing elements
   *                        accross processors. Defaults to std::hash<K>().
   */
  distributed_counting_map (MPI_Comm mpi_comm, int comm_size,
          HASHER hashFunction = HASHER())
      : _base_class(mpi_comm, comm_size, hashFunction)
  {
    // add comm layer receive callbacks
    using namespace std::placeholders;
    this->commLayer.addReceiveCallback(_base_class::INSERT_MPI_TAG,
        std::bind(&distributed_counting_map::receivedCountCallback,
                  this, _1, _2, _3));
    this->commLayer.addReceiveCallback(_base_class::LOOKUP_MPI_TAG,
        std::bind(&distributed_counting_map::receivedLookupCallback,
                  this, _1, _2, _3));
    this->commLayer.addReceiveCallback(_base_class::LOOKUP_ANSWER_MPI_TAG,
        std::bind(&distributed_counting_map::receivedLookupAnswerCallback,
                  this, _1, _2, _3));

    // start the threads in the comm layer (if not already running)
    this->commLayer.initCommunication();
  }

  /**
   * @brief Destructor.
   */
  virtual ~distributed_counting_map()
  {
//    this->commLayer.finishTag(_base_class::INSERT_MPI_TAG);
//    this->commLayer.finishTag(_base_class::LOOKUP_MPI_TAG);
//    this->commLayer.finishTag(_base_class::LOOKUP_ANSWER_MPI_TAG);
	  this->commLayer.finishCommunication();
  }

  /**
   * @brief Inserts the given key into the distributed counting map.
   *
   * If this key has been previously inserted, this increases its count by 1.
   * Otherwise the key is newly inserted and its count set to 1.
   *
   * @param key     The key to insert.
   */
  void insert(const K& key)
  {
    int targetRank = this->getTargetRank(key);
    this->sendKey(key, targetRank, _base_class::INSERT_MPI_TAG);
    this->has_pending_inserts = true;
  }

  /**
   * @brief   construct a new element in place, with same semantic as std::unordered_map::emplace.
   *
   * construct new element to insert in place.  using variadic template parameters allows
   * both multimap and counting map to provide the same interface.
   * Does not return iterator and bool like unordered_map.
   *
   * @param args  variable number of arguments.
   */
  template <class... Args>
  void emplace(Args&&... args) {
      insert(args...);
  }

  /**
   * @brief Populates the distributed map by bulk-inserting elements from the
   *        given iterators.
   *
   * This function blocks until all elements haven been inserted from all
   * processors.
   *
   * @tparam Iterator   An input iterator with value_type std::pair<K,T>.
   * @param begin       Iterator to the first element which will be inserted.
   * @param end         Iterator to one element past the last element to insert.
   */
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
      this->sendKey(*it, targetRank, _base_class::INSERT_MPI_TAG);
    }
    this->commLayer.flush(_base_class::INSERT_MPI_TAG);
  }

protected:
  /**
   * @brief Callback function for received inserts.
   *
   * Deserializes the message data and sets or increases the local count of the
   * received keys.
   *
   * @param msg         The message data received.
   * @param count       The number of bytes received.
   * @param fromRank    The source rank of the message.
   */
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
};

} // namespace bliss
} // namespace index

#endif // BLISS_DISTRIBUTED_MAP_HPP
