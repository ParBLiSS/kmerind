/**
 * @file    MPISendBuffer.hpp
 * @ingroup
 * @author  tpan
 * @brief
 * @details
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef MPISENDBUFFER_HPP_
#define MPISENDBUFFER_HPP_

#include "mpi.h"

#include <cassert>
#include <unistd.h>  // for usleep

#include <vector>

#include "config.hpp"

namespace bliss
{
namespace io
{

/**
 * @class     bliss::io::MPISendBuffer
 * @brief     a data structure to buffer and send MPI messages
 * @details   THREAD_SAFE enables omp_lock for OMP thread safety.  this likely
 *            introduces a good amount of contention.
 *
 * TODO: 1. abstract the MPI requirement.  in fact, abstract the send/flush calls - make it a generic (thread-safe) buffer.
 * TODO: replace all the printf debug output with proper logging
 */
template<typename T, bool THREAD_SAFE=true>
class MPISendBuffer
{
public:
  /// The MPI Tag for signaling that no more messages will be sent.
  static constexpr int END_TAG = 0;

  /// The value type of the elements
  typedef T ValueType;

  /**
   * @brief Creates a new MPISendBuffer instance.
   *
   * @param _comm       The MPI Communicator to be used for buffered
   *                    communication.
   * @param _targetRank The rank of the target process to which the messages are
   *                    sent.
   * @param nbytes      The size of the buffer in bytes. Once this many bytes
   *                    have been added to the buffer, the buffered data will be
   *                    sent to the target process.
   */
  MPISendBuffer(MPI_Comm _comm, const int _targetRank, const size_t nbytes)
     : accepting(true), targetRank(_targetRank), comm(_comm),
       send_request(MPI_REQUEST_NULL), nthreads(0), rank(0), total(0)
  {
    // initialize internal buffers (double buffering)
    capacity = nbytes / sizeof(T);
    active_buf.reserve(capacity);
    inactive_buf.reserve(capacity);

   // printf("targetrank is %d \n", targetRank);

    // get MPI Communicator rank and size
    MPI_Comm_rank(comm, &rank);
    int nprocs;
    MPI_Comm_size(comm, &nprocs);

    // initialize locks in case OpenMP is used
#ifdef USE_OPENMP
    if (THREAD_SAFE)
      omp_init_lock(&lock);
    nthreads = omp_get_thread_num();
#endif
  }


  /**
   * @brief Destroys the MPIBuffer instance.
   * @details
   *    All buffers should be emtpy (i.e. all elements sent) by the time
   *    of calling this function.
   */
  virtual ~MPISendBuffer()
  {
    // both buffers should be empty at this point.
    assert(active_buf.size() == 0);
    assert(inactive_buf.size() == 0);

    // destroy the OpenMP lock
#ifdef USE_OPENMP
    if (THREAD_SAFE)
      omp_destroy_lock(&lock);
#endif
  }

  /**
   * @brief Adds a new element to the send buffer.
   *
   * @param val  The new element that is added to the send buffer.
   */
  void buffer(T val)
  {
    // if this assert fails, the calling thread's logic is faulty.
    assert(accepting);

    // locking.
#ifdef USE_OPENMP
    if (THREAD_SAFE)
      omp_set_lock(&lock);
#endif
    // store value
    //printf("val: %lu %lu %f\n", val.id.composite, val.kmer, val.qual);

    // careful.  when growing, the content is automatically copied to new and
    // destroyed in old.  the destruction could cause double free error or segv.

    active_buf.push_back(val);
    // if full, call send (block if other buffer is sending)
    if (active_buf.size() == capacity) {
      send();
    }

#ifdef USE_OPENMP
    if (THREAD_SAFE)
      omp_unset_lock(&lock);
#endif
  }

  /**
   * @brief Flushes the buffer, i.e. sends all current elements to the target
   *        processes.
   */
  void flush()
  {
    //printf("flushing %d target %d.\n", rank, targetRank);

    accepting = false;

    // no more entries to buffer. send all remaining.
    send();

    // one more time to clear double buffer and wait for previous send to finish.
    send();

    // send termination signal (send empty message).
    //printf("%d done sending to %d\n", rank, targetRank); fflush(stdout);
    MPI_Send(nullptr, 0, MPI_INT, targetRank, END_TAG, comm);


    printf("%d.%d  %ld flushed to %d.\n", rank, nthreads, total, targetRank);
    fflush(stdout);
  }

  /**
   * @brief Returns the size of the active buffer.
   *
   * @return The number of elements in the active buffer.
   */
  size_t size()
  {
    return active_buf.size();
  }

protected:

  // 2 buffers per target (double buffering)
  /// The front buffer (this one is used for adding new elements)
  std::vector<T> active_buf;
  /// The back buffer (this one is used in sending)
  std::vector<T> inactive_buf;

  /// Whether the buffer is still accepting new elements (still free capacity)
  bool accepting;

  /// The rank of the target processor
  const int targetRank;

  /// The capacity of each buffer (in number of elements)
  size_t capacity;

  /// The MPI Communicator for sending
  MPI_Comm comm;

  /// Send status of currently inactive buffer
  MPI_Request send_request;

  /// The number of active OpenMP threads
  int nthreads;

  /// The MPI rank of this process
  int rank;

  /// The total number of elements sent (for stats/logging purposes)
  size_t total;

#ifdef USE_OPENMP
    omp_lock_t lock;
#endif


  void send()
  {
    //printf("SEND %d -> %d,  inactive size: %lu, active size: %lu\n",
    //       rank, targetRank, inactive_buf.size(), active_buf.size());

    // wait for previous send operations to complete
    wait_for_finish(send_request);

    // clear the inactive buf
    inactive_buf.clear();

    // swap active and inactive buffer  (this swaps the contents)
    int s = active_buf.size();
    inactive_buf.swap(active_buf);
    assert(active_buf.size() == 0);
    assert(inactive_buf.size() == s);

    // async send inactive buffer if it's not empty.  Requires that remote side
    // uses probe to get the size of data first.
    if (inactive_buf.size() > 0) {
      total += inactive_buf.size();
      // TODO: use proper MPI datatypes to send the buffer content rather than
      //       raw bytes
      MPI_Isend(inactive_buf.data(), inactive_buf.size() * sizeof(T),
                MPI_UNSIGNED_CHAR, targetRank, nthreads + 1, comm, &send_request);
      //printf("SEND %d -> %d, number of entries %ld, total %ld\n",
      //       rank, targetRank, inactive_buf.size(), total); fflush(stdout);
    }
  }

private:

  /**
   * @brief Blocks till the status returns a `finished`.
   *
   * @param status_request  The MPI_Request object.
   */
  void wait_for_finish(MPI_Request& status_request)
  {
    /*
     * MPI_Wait results in busy CPU waiting for some MPI implementations.
     * This happend with OpenMPI.
     *
     * Thus we implement our own "busy" waiting routing that continuously
     * calls MPI_Test, but yields the CPU by sleeping in between calls to
     * MPI_Test.
     *
     * The sleep duration is set in the class constructor, depending on the
     * number of processors.
     */
    int finished = 0;
    do {
      // repeatedly check to see if request was completed.
      // okay to test against a NULL request.

      // if being sent, wait for that to complete
      MPI_Test(&status_request, &finished, MPI_STATUS_IGNORE);
      // sleep for 1 ms (this yields the CPU)
      usleep(1000);
    } while(finished == 0);
  }
};

} /* namespace io */
} /* namespace bliss */

#endif /* MPISENDBUFFER_HPP_ */
