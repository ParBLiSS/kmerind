/**
 * @file    distributed_map.hpp
 * @ingroup index
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @author  Tony Pan <tpan7@gatech.edu>
 * @brief   Implements the distributed_multimap and distributed_counting_map
 *          data structures.
 *
 * Copyright (c) 2014 Georgia Institute of Technology. All Rights Reserved.
 *
 * TODO add Licence
 */

#ifndef BLISS_DISTRIBUTED_MAP_HPP
#define BLISS_DISTRIBUTED_MAP_HPP

#include <utility> 			// for std::pair
#include <unordered_map> 	// local storage hash table  // for multimap
//#include <sparsehash/dense_hash_map>  // local storage hash table, google dense hash map.
#include <vector>
#include <functional> 		// for std::function and std::hash
#include <limits>
#include <stdexcept>
#include <algorithm> 		// for std::max
#include <atomic>
#include <tuple>

#include "iterators/concatenating_iterator.hpp"

// include MPI
#include <mpi.h>

namespace bliss
{
namespace index
{

/// Type for the counting map
typedef uint32_t count_t;

// TODO: are the "have_pending_insert and have_pending_lookup necessary?

/**
 * @brief A shared base class for the distributed_multimap and
 *        distributed_counting_map implementations.
 * @detail  distributed map supports Insert and Query functions, but not delete.
 * 		each of these is supported in an asynchronous, distributed fashion.
 * 		insert and query functions essentially initiates communication with remote
 * 		processes.
 *
 * 		remote processes perform actual insert and query via callback
 *		functions.  internally, we use thread-local unordered_map
 *		to store the actual key-value pairs.
 *
 *		3 levels of hashing exists:  MPI process, thread, and local container.
 *		The user can specify 3 separate hash functions, e.g. use prefix, infix, and suffix
 *		of a kmer.  MPI processes therefore have disjoint subspace of the keys.  Threads
 *		also have disjoint subspace of the keys.  Using different hash functions
 *		avoid overloading particular hashtable buckets.
 *
 * @note  supplied insert and lookup callbacks should be sufficient, and hidden from user
 * 		lookup result callback should be supplied by user, and should be able to operator
 * 		on bulk data.
 *
 * @tparam K                    The key type.
 * @tparam T                    The value type.
 * @tparam CommunicationLayer   The CommunicationLayer class type.
 * @tparam LocalContainer       The local container type (unordered_map or
 *                              unodered_multimap)
 * @tparam ThreadHasher		    Hash function for hashing to key to thread id
 * @tparam MPIHasher		    Hash function for hashing to key to MPI rank.
 */
template<typename CommunicationLayer,
         typename LocalContainer, typename ThreadHasher, typename MPIHasher>
class _distributed_map_base
{
protected:
	// iterator for a the local container.
  typedef typename LocalContainer::iterator local_thread_iterator;
    // const iterator for a local container.
  typedef typename LocalContainer::const_iterator local_thread_const_iterator;

  typedef typename LocalContainer::key_type K;
  typedef typename LocalContainer::mapped_type T;


  /// MPI message tag for reserve size.  this is meant to be sent only to self.
  static constexpr int RESERVE_TAG = 12;

  typedef std::unordered_map<K,
//  typedef dense_hash_map<K,
                             T,
                             typename LocalContainer::hasher,
                             typename LocalContainer::key_equal,
                             typename LocalContainer::allocator_type> MapType;
  typedef std::unordered_multimap<K,
                             T,
                             typename LocalContainer::hasher,
                             typename LocalContainer::key_equal,
                             typename LocalContainer::allocator_type> MultimapType;

public:
  /// The iterator type of the entire set of thread-local containers on an MPI process
  typedef typename bliss::iterator::ConcatenatingIterator<local_thread_iterator>  local_iterator;
  /// The constant iterator type of the entire set of thread-local containers on an MPI process
  typedef typename bliss::iterator::ConcatenatingIterator<local_thread_const_iterator>  local_const_iterator;

  /// MPI message tag for inserts
  static constexpr int INSERT_MPI_TAG = 13;
  /// MPI message tag for lookup queries
  static constexpr int LOOKUP_MPI_TAG = 14;
  /// MPI message tag for answers to lookup queries
  static constexpr int LOOKUP_ANSWER_MPI_TAG = 15;

  /**
   * @brief Returns an iterator to the first element of the local container(s).
   *
   * @return Iterator to the first local element.
   */
  local_iterator begin()
  {
    local_iterator it;
    for (int i = 0; i < nThreads; ++i)
	it.addRange(local_map[i].begin(), local_map[i].end());
    return it;
  }

  /**
   * @brief Returns an iterator to the element following the last element of
   * the local container(s).
   *
   * @return Iterator to the element following the last local element.
   */
  local_iterator end()
  {
    return local_iterator();
  }

  /**
   * @brief Returns a const iterator to the first element of the local container(s).
   *
   * @return Iterator to the first local element.
   */
  local_const_iterator cbegin() const
  {
    local_const_iterator it;
    for (int i = 0; i < nThreads; ++i)
	it.addRange(local_map[i].cbegin(), local_map[i].cend());
    return it;
  }

  /**
   * @brief Returns a const iterator to the element following the last element of
   * the local container(s).
   *
   * @return Iterator to the element following the last local element.
   */
  local_const_iterator cend() const
  {
    return local_const_iterator();
  }

  /**
   * @brief initiates a local map size increment request
   * @details  this generates a local message that INCREASES the size of the local_map by size.
   *        this is sent as a message to self because only Callback Thread can access the local_map
   *
   * @param size
   */
  void reserve(const size_t& size) {
    this->commLayer.sendMessage(&size, sizeof(size), commRank, RESERVE_TAG);
    this->commLayer.flush(RESERVE_TAG);
  }

  /**
   * @brief Reserves size
   * @param size_hint
   */
  void mapReserveCallback(uint8_t* msg, std::size_t count, int fromRank) {
    assert(count == sizeof(size_t));
    assert(fromRank == commRank);

    size_t addl = *(reinterpret_cast<size_t*>(msg)) / nThreads;

    // resize the hashtables.
#pragma omp parallel num_threads(nThreads) shared(addl)
    {
      int tid = omp_get_thread_num();
      size_t size = this->local_map[tid].size();
      this->local_map[tid].reserve(size+addl);   // causes a rehash
    }
  }

  /**
   * @brief   Returns the local map's size.  primarily for debugging.
   * @return  size of the local map.
   */
  const size_t local_size() const
  {
    size_t result = 0;
    for (int i = 0; i < nThreads; ++i) {
      result += local_map[i].size();
    }
    return result;
  }

  /**
   * @brief   Returns the local map's size.  primarily for debugging.
   * @return  size of the local map.
   */
  std::vector<size_t> local_sizes() const
  {
    std::vector<size_t> result;
    for (int i = 0; i < nThreads; ++i) {
      result.push_back(local_map[i].size());
    }
    return result;
  }

  /**
   * @brief flush all the insertion messages (across all MPI processes)
   * @details  exits this call when all incoming messages from other MPI processes
   * 	have arrived.
   */
  void flushInsert()
  {
    DEBUGF("FLUSH DISTR MAP Insert");
    if (has_pending_inserts.exchange(false, std::memory_order_acquire))
    {
      this->commLayer.flush(INSERT_MPI_TAG);
     }
  }
  /**
   * @brief flush all the lookup messages (across all MPI processes)
   * @details  exits this call when all incoming messages from other MPI processes
   * 	have arrived.
   */
  void flushLookup()
  {
   DEBUGF("FLUSH DISTR MAP Lookup");
   if (has_pending_lookups.exchange(false, std::memory_order_acquire))
   {
    this->commLayer.flush(LOOKUP_MPI_TAG);
    this->commLayer.flush(LOOKUP_ANSWER_MPI_TAG);
   }
  }

  /**
   * @brief Flushes all pending operations.
   * @note  SHOULD BE CALLED BY A SINGLE THREAD
   * Since all insert and lookup operations are executed asynchronosly,
   * calling this function ensures that all pending operations have been
   * executed, including that all pending lookups have returned an answer.
   */
  void flush()
  {
    flushInsert();
    flushLookup();
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

    DEBUG("Rank " << commRank << " LOOKUP " << key << " at rank " << targetRank);

    this->sendKey(key, targetRank, LOOKUP_MPI_TAG);
    has_pending_lookups.store(true);
  }

  /**
   * @brief Sets the callback function for received answers to asynchronously
   *        posted lookups.
   * @note  previously set callback function will be replaced.
   *
   * @param callbackFunction    The function to call with all received answers
   *                            to lookups.
   */
  void setLookupAnswerCallback(const std::function<void(std::pair<K, T>*, std::size_t)>& callbackFunction)
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
   * performed using openmp.
   *
   * @param count   The key count threshold. Everything lower than this will be
   *                removed.
   */
  void filter(const std::size_t lower_threshold)
  {

    int nt = this->nThreads;

#pragma omp parallel default(none) num_threads(nt)
    {
    // iterate through all keys
      int tid = omp_get_thread_num();

      for (auto iter = this->local_map[tid].begin(); iter!=this->local_map[tid].end();)
      {
        // get end of the range of identical key
        auto cur_end_iter = this->local_map[tid].equal_range(iter->first)->second;
        std::size_t cur_count = getLocalCount(local_map[tid], iter);  // note we're only operation on 1 map at a time.
        if (cur_count < lower_threshold)
        {
          // remove all entries with this key, this will invalidate the `iter`
          // iterator
          this->local_map[tid].erase(iter, cur_end_iter);
        }
        // advance loop iterator
        iter = cur_end_iter;
      }
    }
  }
//  void sequentialFilter(const std::size_t lower_threshold)
//  {
//    for (int tid = 0; tid < nThreads; ++tid) {
//      // iterate through all keys
//      for (auto iter = this->local_map[tid].begin(); iter!=this->local_map[tid].end();)
//      {
//        // get end of the range of identical key
//        auto cur_end_iter = this->local_map[tid].equal_range(iter->first)->second;
//        std::size_t cur_count = getLocalCount(local_map[tid], iter);  // note we're only operation on 1 map at a time.
//        if (cur_count < lower_threshold)
//        {
//          // remove all entries with this key, this will invalidate the `iter`
//          // iterator
//          this->local_map[tid].erase(iter, cur_end_iter);
//        }
//        // advance loop iterator
//        iter = cur_end_iter;
//      }
//    }
//  }

  /**
   * @brief Returns a histrogram of counts for all elements in the distributed
   *        hash table.
   *
   * The valid histrogram is returned on ALL MPI processes. This function has
   * to be called by ALL MPI processes of the given MPI communicator.
   *
   * performed using openmp
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
    int nt = this->nThreads;
#pragma omp parallel default(none) num_threads(nt) reduction(max: local_max_count)
    {
      int tid = omp_get_thread_num();
      for (auto iter=this->local_map[tid].begin(); iter!=this->local_map[tid].end();
           iter=this->local_map[tid].equal_range(iter->first)->second)
      {
        std::size_t count = getLocalCount(local_map[tid], *iter);
        local_max_count = std::max<uint64_t>(local_max_count, count);
      }
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
#pragma omp parallel default(none) num_threads(nt)
    {
      int tid = omp_get_thread_num();
      for (auto iter=this->local_map[tid].begin(); iter!=this->local_map[tid].end();
           iter=this->local_map[tid].equal_range(iter->first)->second)
      {
        std::size_t count = getLocalCount(local_map[tid], *iter);
        local_count_hist[count]++;
      }
    }
    // then accumulate across all processors
    std::vector<int> count_hist(all_max_count+1, 0);
    MPI_Allreduce(&local_count_hist[0], &count_hist[0], all_max_count+1,
                  MPI_INT, MPI_SUM, this->comm);

    return count_hist;
  }
//  std::vector<int> sequentialCountHistrogram()
//  {
//    // first determine the maximum count
//    uint64_t local_max_count = 0; // use uint64_t for all systems!
//    for (int i = 0; i < nThreads; ++i) {
//      for (auto iter=this->local_map[i].begin(); iter!=this->local_map[i].end();
//           iter=this->local_map[i].equal_range(iter->first)->second)
//      {
//        std::size_t count = getLocalCount(local_map[i], *iter);
//        local_max_count = std::max<uint64_t>(local_max_count, count);
//      }
//    }
//
//    // cast max count to int and check that it doesn't overflow
//    if (local_max_count >= std::numeric_limits<int>::max())
//    {
//      throw std::range_error("Histrogram of counts: maximum count exceeds integer range");
//    }
//    int max_count = static_cast<int>(local_max_count);
//
//    // get max accross all processors
//    int all_max_count;
//    MPI_Allreduce(&max_count, &all_max_count, 1, MPI_INT, MPI_MAX, this->comm);
//
//    // count the counts to create local histogram
//    std::vector<int> local_count_hist(all_max_count+1, 0);
//    for (int i = 0; i < nThreads; ++i) {
//      for (auto iter=this->local_map[i].begin(); iter!=this->local_map[i].end();
//           iter=this->local_map[i].equal_range(iter->first)->second)
//      {
//        std::size_t count = getLocalCount(local_map[i], *iter);
//        local_count_hist[count]++;
//      }
//    }
//
//    // then accumulate across all processors
//    std::vector<int> count_hist(all_max_count+1, 0);
//    MPI_Allreduce(&local_count_hist[0], &count_hist[0], all_max_count+1,
//                  MPI_INT, MPI_SUM, this->comm);
//
//    return count_hist;
//  }

  /**
   * initializes the communication for a distributed map.
   */
  void init() {
    DEBUGF("Initialize distributed multimap's communication");
    // start the threads in the comm layer (if not already running)
    this->commLayer.initCommunication();
  }


  /**
   * waits for commLayer to finish, which guarantees that all communications are done (to avoid timing issue with MPI_FINALIZE and distributed map object clean up.
   */
  void finish() {
    this->commLayer.finish(INSERT_MPI_TAG);
    this->commLayer.finish(LOOKUP_MPI_TAG);
    this->commLayer.finish(LOOKUP_ANSWER_MPI_TAG);
    this->commLayer.finish(RESERVE_TAG);


    this->commLayer.finishCommunication();
  }

protected:
  /**
   * @brief Constructs the shared base class.  protected:  use derived class
   *
   * @param mpi_comm        The MPI communicator to pass onto the communication
   *                        layer.
   * @param comm_size       The size of the MPI communicator, needed for
   *                        initialization.
   */
  _distributed_map_base(MPI_Comm mpi_comm, int comm_size, int num_threads = 1, ThreadHasher thread_hash_func = ThreadHasher(), MPIHasher mpi_hash_func = MPIHasher())
      : commLayer(mpi_comm, comm_size, num_threads), comm(mpi_comm), nThreads(num_threads), threadHashFunc(thread_hash_func), mpiHashFunc(mpi_hash_func),
        local_map(num_threads), threadKeys(num_threads),
        has_pending_inserts(false), has_pending_lookups(false)
  {
	  MPI_Comm_rank(comm, &commRank);

    using namespace std::placeholders;
    this->commLayer.addReceiveCallback(RESERVE_TAG,
        std::bind(&_distributed_map_base::mapReserveCallback,
                  this, _1, _2, _3));


  }

  /**
   * @brief  virtual destructor.  waits for underlying commLayer to finish
   */
  virtual ~_distributed_map_base() {
    this->finish();
  }


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
    this->sendPair(std::make_pair(key, value), dstRank, tag);
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
    return mpiHashFunc(key) % size;
  }

  /**
   * @brief Returns the target thread id for a given key.
   *
   * Uses the given hash function to determine the thread responsible for storing
   * the key and its metadata.
   *
   * @param key     The key for which to determine the target thread.
   *
   * @return        The target thread to which the given key
   *                belongs.
   *
   */
  int getTargetThread(const K& key)
  {
    // get the target rank for the processor
    return threadHashFunc(key) % nThreads;
  }


  /**
   * @brief Returns the local count of a given key in one unordered multimap
   *
   * @param localMap    The local map data structure.
   * @param item        The item to lookup the count for.
   *
   * @return    The count of elements with the given item's key.
   */
  template <typename Container = LocalContainer>
  inline typename std::enable_if<std::is_same<Container, MultimapType>::value, std::size_t>::type
  getLocalCount(const Container& localMap, const std::pair<K, T>& item)
  {
    // use the maps count() function
    return localMap.count(item.first);
  }

  /**
   * @brief Returns the local count of a given key in one unordered multimap.  this one is for count map.
   *
   * @param localMap    The local map data structure.
   * @param item        The item to lookup the count for.
   *
   * @return    The count of elements with the given item's key.
   */
  template <typename Container = LocalContainer>
  inline typename std::enable_if<std::is_same<Container, MapType>::value, std::size_t>::type
  getLocalCount(const Container& localMap, const std::pair<K, T>& item)
  {
    // this is a counting map: i.e. the value of (key,value) is the count
    assert(localMap.find(item.first) != localMap.end());
    return static_cast<std::size_t>(item.second);
  }


  /**
   * @brief Callback function for received lookup messages.
   *
   * Performs the actual lookup in the local data structure and sends a message
   * as reply.  OMP parallelized
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

    int nt = this->nThreads;

    // for all received requests, send the value from the lookup
#pragma omp parallel for default(none) num_threads(nt) shared(key_count, keys, fromRank)
    for (int i = 0; i < key_count; ++i)
    {
      // check if exists and then send
      auto range = this->local_map[this->getTargetThread(keys[i])].equal_range(keys[i]);
      for (auto it = range.first; it != range.second; ++it)
      {
        // send the results to the requesting processor
        this->sendPair(*it, fromRank, LOOKUP_ANSWER_MPI_TAG);
      }
    }
  }

//  /*
//   * @brief Callback function for received lookup messages.
//   *
//   * Performs the actual lookup in the local data structure and sends a message
//   * as reply.
//   *
//   * @param msg         The message data received.
//   * @param count       The number of bytes received.
//   * @param fromRank    The source rank of the message.
//   */
//  void sequentialReceivedLookupCallback(uint8_t* msg, std::size_t count, int fromRank)
//  {
//    // deserialize
//    K* keys = reinterpret_cast<K*>(msg);
//    int key_count = count / sizeof(K);
//
//    // for all received requests, send the value from the lookup
//    for (int i = 0; i < key_count; ++i)
//    {
//      // check if exists and then send
//      auto range = this->local_map[this->getTargetThread(keys[i])].equal_range(keys[i]);
//      for (auto it = range.first; it != range.second; ++it)
//      {
//        // send the results to the requesting processor
//        this->sendPair(*it, fromRank, LOOKUP_ANSWER_MPI_TAG);
//      }
//    }
//  }


  /**
   * @brief Callback function for received lookup answers.
   *
   * Deserializes the answer to lookup queries and calls the
   * given callback function.
   *
   * @note	user supplied callback function should operate on these in bulk
   * 	so that they can better express parallelism.
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

    lookupAnswerCallbackFunc(elements, element_count);
  }

protected:
  /******************
   *  Data members  *
   ******************/

  /// The async communication layer abstraction, used for sending and receiving messages
  CommunicationLayer commLayer;

  /// The MPI communicator used for distributing the input
  MPI_Comm comm;

  /// number of threads participating in callback functions.
  int nThreads;

  /// The hash functijon used for distributing elements across threads
  ThreadHasher threadHashFunc;

  /// The hash function used for distributing input across processors
  MPIHasher mpiHashFunc;

  /// The local container data-structure. For most cases this is either
  /// a std::unordered_map or std::unordered_multimap
  std::vector<LocalContainer> local_map;

  /// storage for temporarily generated key to thread mapping.  reused for each received buffer to process.
  std::vector<int> threadKeys;

  /// Whether there are pending insert operations
  std::atomic<bool> has_pending_inserts;
  /// Whether there are pending lookup operations
  std::atomic<bool> has_pending_lookups;

  /// The async callback function, which is called when an answer to a
  /// lookup query is received
  std::function<void(std::pair<K, T>*, std::size_t)> lookupAnswerCallbackFunc;

  std::mutex mutex;

  int commRank;

};


/**
 * @brief   A distributed, asynchronous multimap.
 *
 * This can hold multiple elements with identical key.
 * @note	this is useful for position kmer index, and position+quality score kmer index.
 *
 * @tparam K                    The key type.
 * @tparam T                    The value type.
 * @tparam CommunicationLayer   The CommunicationLayer class type.
 * @tparam LocalHasher		    Hash function used by local unordered_map backend store.
 * @tparam ThreadHasher		    Hash function for hashing to key to thread id
 * @tparam MPIHasher		    Hash function for hashing to key to MPI rank.
 */
template<typename K, typename T, typename CommunicationLayer, typename LocalHasher = std::hash<K>, typename ThreadHasher = std::hash<K>, typename MPIHasher = std::hash<K> >
class distributed_multimap
  : public _distributed_map_base<
              CommunicationLayer, std::unordered_multimap<K, T, LocalHasher>, ThreadHasher, MPIHasher >
{
public:
  /// The baseclass type
  typedef _distributed_map_base<CommunicationLayer,
                                std::unordered_multimap<K, T, LocalHasher>, ThreadHasher, MPIHasher > _base_class;
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
  distributed_multimap (MPI_Comm mpi_comm, int comm_size, int num_threads = 1,
        ThreadHasher thread_hash_func = ThreadHasher(),
        MPIHasher mpi_hash_func = MPIHasher())
      : _base_class(mpi_comm, comm_size, num_threads, thread_hash_func, mpi_hash_func)
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
  }

  /**
   * @brief Destructor.
   */
  virtual ~distributed_multimap() {
    // finishCommunication called in base class dtor
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
    this->has_pending_inserts.store(true, std::memory_order_release);
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
    this->has_pending_inserts.store(true, std::memory_order_release);
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
   * @note internally, the function is parallel.
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

    //====  now use parallel for to compute the hash,
    int nt = this->nThreads;
    // clear and make room.  callbacks are invoked in series so only 1 thread uses threadKeys.
    this->threadKeys.clear();
    if (this->threadKeys.capacity() < element_count) this->threadKeys.reserve(element_count);

    // compute threadhash and store in vector in parallel, (perfectly parallelizable, compute heavy)
#pragma omp parallel for num_threads(nt) shared(element_count, elements) default(none)
    for (int i = 0; i < element_count; ++i) {
      this->threadKeys[i] = this->getTargetThread(elements[i].first);  // insert hash value
    }

    // in parallel, each thread walks through entire vector (read only) and if hash matches insert into per thread map  (read vec is not compute heavy.)
#pragma omp parallel num_threads(nt) shared(element_count, elements) default(none)
    {
      int tid = omp_get_thread_num();
      for (int i = 0; i < element_count; ++i) {
        // insert value if hash matches thread id
        if (this->threadKeys[i] == tid) this->local_map[tid].insert(elements[i]);
      }
    }

    // also alternatively
    // T^2 vectors of pointers.
    // first parallel generate hash, each thread insert vector into vector[hash][tid]
    // then parallel insert  each thread process all vectors in vector[tid]
    // potentially more cost.
      // performance is kind of similar for small number of thread.  there is more variability here.


    // alternatively.
    // store entry in per thread concurrentqueues
    // use omp parallel to pop from queue and insert into the per thread map.
    // this approach is likely more costly, as contended insertion into concurrent queue (up to nThreads conflicts) is serialized.
    //    and popping may be more costly because of internal state tracking.

  }
//  void sequentialReceivedInsertCallback(uint8_t* msg, std::size_t count, int fromRank)
//  {
//    // deserialize
//    std::pair<K, T>* elements = reinterpret_cast<std::pair<K, T>*>(msg);
//    int element_count = count / sizeof(std::pair<K, T>);
//
//    for (int i = 0; i < element_count; ++i) {
//
//      this->local_map[this->getTargetThread(elements[i].first)].insert(elements[i]);
//    }
//  }
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
 * @tparam LocalHasher		    Hash function used by local unordered_map backend store.
 * @tparam ThreadHasher		    Hash function for hashing to key to thread id
 * @tparam MPIHasher		    Hash function for hashing to key to MPI rank.
 */
template<typename K, typename CommunicationLayer,
	typename LocalHasher = std::hash<K>, typename ThreadHasher = std::hash<K>, typename MPIHasher = std::hash<K>>
class distributed_counting_map
 : public _distributed_map_base<CommunicationLayer,
                               std::unordered_map<K, count_t, LocalHasher>, ThreadHasher, MPIHasher >
//                               google::dense_hash_map<K, count_t, LocalHasher>, ThreadHasher, MPIHasher >
{
public:
  /// The baseclass type
  typedef _distributed_map_base<CommunicationLayer,
                                std::unordered_map<K,count_t, LocalHasher>, ThreadHasher, MPIHasher > _base_class;
                                //google::dense_hash_map<K,count_t, LocalHasher>, ThreadHasher, MPIHasher > _base_class;
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
  distributed_counting_map (MPI_Comm mpi_comm, int comm_size, int num_threads = 1,
          ThreadHasher thread_hash_func = ThreadHasher(),
          MPIHasher mpi_hash_func = MPIHasher())
      : _base_class(mpi_comm, comm_size, num_threads, thread_hash_func, mpi_hash_func)
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

  }


  /**
   * @brief Destructor.
   */
  virtual ~distributed_counting_map()
  {
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
    this->has_pending_inserts.store(true, std::memory_order_release);
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
   * @note parallel internally vis OMP.
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
    size_t key_count = count / sizeof(K);

    //====  now use parallel for to compute the hash,
    int nt = this->nThreads;
    // clear and make room.  callbacks are invoked in series so only 1 thread uses threadKeys.
    this->threadKeys.clear();
    if (this->threadKeys.capacity() < key_count) this->threadKeys.reserve(key_count);

    // compute threadhash and store in vector in parallel, (perfectly parallelizable, compute heavy)
#pragma omp parallel for num_threads(nt) shared(key_count, keys) default(none)
    for (size_t i = 0; i < key_count; ++i) {
      this->threadKeys[i] = this->getTargetThread(keys[i]);  // insert hash value
    }

    // in parallel, each thread walks through entire vector (read only) and if hash matches insert into per thread map  (read vec is not compute heavy.)
#pragma omp parallel num_threads(nt) shared(key_count, keys) default(none)
    {
      int tid = omp_get_thread_num();
      for (size_t i = 0; i < key_count; ++i) {
        // insert value if hash matches thread id
        if (this->threadKeys[i] == tid) this->local_map[tid][keys[i]]++;
      }
    }

    // also alternatively
    // T^2 vectors of pointers.
    // first parallel generate hash, each thread insert vector into vector[hash][tid]
    // then parallel insert  each thread process all vectors in vector[tid]
    // potentially more cost.
      // performance is kind of similar for small number of thread.  there is more variability here.


    // alternatively.
    // store entry in per thread concurrentqueues
    // use omp parallel to pop from queue and insert into the per thread map.
    // this approach is likely more costly, as contended insertion into concurrent queue (up to nThreads conflicts) is serialized.
    //    and popping may be more costly because of internal state tracking.

  }
//  void sequentialReceivedCountCallback(uint8_t* msg, std::size_t count, int fromRank)
//  {
//    // deserialize
//    K* keys = reinterpret_cast<K*>(msg);
//    int key_count = count / sizeof(K);
//
//    // insert all elements into the hash table
//    for (int i = 0; i < key_count; ++i)
//    {
//      this->local_map[this->getTargetThread(keys[i])][keys[i]]++;
//    }
//  }

};

} // namespace bliss
} // namespace index

#endif // BLISS_DISTRIBUTED_MAP_HPP
