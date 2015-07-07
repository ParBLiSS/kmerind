/**
 * @file    mxx_support.hpp
 * @ingroup
 * @author  tpan
 * @brief
 * @details
 *
 * Copyright (c) 2015 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef SRC_IO_MXX_SUPPORT_HPP_
#define SRC_IO_MXX_SUPPORT_HPP_

#include "mxx/collective.hpp"
#include "mxx/reduction.hpp"

#include "common/kmer.hpp"
#include "io/fastq_loader.hpp"

#include "utils/logging.h"

//#include "utils/system_utils.hpp"

#include <algorithm>  // for std::min


namespace mxx {

  template<unsigned int size, typename A, typename WT>
  class datatype<typename bliss::common::Kmer<size, A, WT> > :
    public datatype_contiguous<typename bliss::common::Kmer<size, A, WT>::KmerWordType,
      bliss::common::Kmer<size, A, WT>::nWords> {};

  template<unsigned int size, typename A, typename WT>
  class datatype<const typename bliss::common::Kmer<size, A, WT> > :
    public datatype_contiguous<typename bliss::common::Kmer<size, A, WT>::KmerWordType,
      bliss::common::Kmer<size, A, WT>::nWords> {};

  template<>
  class datatype<bliss::io::FASTQ::SequenceId > :
    public datatype_contiguous<decltype(bliss::io::FASTQ::SequenceId::file_pos), 1> {};

}





namespace mxx2 {


  /// in place bucketing.  uses an extra vector of same size as msgs.  hops around msgs vector until no movement is possible.
  /// complexity is O(b) * O(n).  scales badly, same as bucketing_copy, but with a factor of 2 for large data.
  /// when fixed n, slight increase with b..  else increases with n.
  template<typename count_t, typename T, typename _TargetP>
  std::vector<count_t> bucketing(std::vector<T>& buffer, _TargetP target_p_fun, int p) {

    //== track the send count
    std::vector<count_t> send_counts(p, 0);

    if (buffer.size() == 0) return send_counts;

    //== track target process assignment
    std::vector<int> pids(buffer.size());
    for (size_t i = 0; i < buffer.size(); ++i)
    {
      pids[i] = target_p_fun(buffer[i]);
      ++(send_counts[pids[i]]);
    }
    // at this point, have target assignment for each data element, and also count for each process bucket.

    // compute the offsets within the buffer
    std::vector<count_t> offset = ::mxx::get_displacements(send_counts);
    std::vector<count_t> maxes = offset;
    for (int i = 0; i < p; ++i) {
      maxes[i] += send_counts[i];
    }


    //== swap elements around.
    T val;
    count_t tar_pos, start_pos;

    int tar_rank;

    // while loop will stop under 2 conditions: 
    //      1. returned to starting position (looped), or
    //      2, tar_pos is the current pos.
    // either way, we need a new starting point.  instead of searching through buffer O(N), search
    // for incomplete buckets via offset O(p).

    for (int i = 0; i < p;) {
      // determine the starting position.
      if (offset[i] == maxes[i]) {
        ++i;  // skip all completed buckets
        continue;  // have the loop check value.
      }
      // get the start pos.
      start_pos = offset[i];

      // set up the variable with the current entry.
      tar_rank = pids[start_pos];
      if (tar_rank > -1) {
        val = ::std::move(buffer[start_pos]);  // value to move
        pids[start_pos] = -2;                // special value to indicate where we started from.

        while (tar_rank > -1) {  // if -1 or -2, then either visited or beginning of chain.
          tar_pos = offset[tar_rank]++;  // compute new position.  earlier offset values for the same pid are should have final values already.

          tar_rank = pids[tar_pos];

//          if (tar_rank == -1) throw ::std::logic_error("target position target rank indicates that it's been visited already.");

          // save the info at tar_pos;
          ::std::swap(val, buffer[tar_pos]);  // put what's in src into buffer at tar_pos, and save what's at buffer[tar_pos]
          pids[tar_pos] = -1;               // mark as visited.

        }  // else already visited, so done.
      }
    }

    // clean up
    std::vector<int>().swap(pids);  // clear pids reallocating.
    std::vector<count_t>().swap(offset);
    std::vector<count_t>().swap(maxes);

    return send_counts;
  }

  /// bucketing using temporary storage.  uses an extra vector of same size as buffer.  hops around buffer vector until no movement is possible.
  template<typename count_t, typename T, typename _TargetP>
  std::vector<count_t> bucketing_copy(std::vector<T> const & buffer, std::vector<T> &bucketed, _TargetP target_p_fun, int p) {

    //== track the send count
    std::vector<count_t> send_counts(p, 0);

    if (buffer.size() == 0) return send_counts;

    //== track target process assignment
    std::vector<int> pids(buffer.size());
    for (size_t i = 0; i < buffer.size(); ++i)
    {
      pids[i] = target_p_fun(buffer[i]);
      ++(send_counts[pids[i]]);
    }
    // at this point, have target assignment for each data element, and also count for each process bucket.

    // compute the offsets within the buffer
    std::vector<count_t> offset = ::mxx::get_displacements(send_counts);

    //== copy.  need to be able to track current position within each block.
    bucketed.resize(buffer.size());
    for (size_t i = 0; i < buffer.size(); ++i)
    {
      bucketed[(offset[pids[i]])++] = ::std::move(buffer[i]);
    }

    // clean up
    std::vector<int>().swap(pids);  // clear pids reallocating.
    std::vector<count_t>().swap(offset);

    return send_counts;
  }

  /// keep the unique entries.  complexity is b * O(N/b), where b is the bucket size, and O(N/b) is complexity of inserting into and copying from set.
  /// when used within bucket, scales with O(N/b), not with b.  this is as good as it gets wrt complexity.
  template<typename SET, typename EQUAL, typename T = typename SET::key_type, typename count_t = int>
  void retain_unique(std::vector<T>& input, std::vector<count_t> &send_counts, bool sorted = false) {
  
    auto newstart = input.begin();
    auto newend = newstart;
    auto start = input.begin();
    auto end = start;
    
    count_t max = *(::std::max_element(send_counts.begin(), send_counts.end()));
    SET set(max);

    for (size_t i = 0; i < send_counts.size(); ++i) {
      end = start + send_counts[i];
        
      if (sorted) {  // already sorted, then just get the unique stuff and remove rest.
        if (i == 0)
          newend = ::std::unique(start, end, EQUAL());
        else 
          newend = ::std::unique_copy(start, end,
                                      newstart, EQUAL());
      } else {  // not sorted, so use an unordered_set to keep the first occurence.

        // sorting is SLOW and not scalable.  use unordered set instead.
        // unordered_set for large data is memory intensive.  depending on use, bucket per processor first.
        set.clear();
        set.insert(start, end);
        newend = ::std::copy(set.begin(), set.end(), newstart);
      }
      
      send_counts[i] = ::std::distance(newstart, newend);  

      start = end;
      newstart = newend;
    }

    // compact.
    input.erase(newend, input.end());
  }



  template<typename MAP, typename REDUCER, typename KEY = typename MAP::key_type, typename T = typename MAP::mapped_type, typename count_t = int>
  void bucket_reduce(std::vector<::std::pair<KEY, T> >& input, std::vector<count_t> &send_counts, REDUCER const& r = REDUCER(), bool sorted = false) {
    auto newstart = input.begin();
    auto newend = newstart;
    auto start = input.begin();
    auto end = start;
    
    count_t max = *(::std::max_element(send_counts.begin(), send_counts.end()));
    MAP map(max);

    KEY key;
    T val;
  
    for (int i = 0; i < send_counts.size(); ++i) {
      end = start + send_counts[i];
  
      map.clear();
      for (auto it = start; it != end; ++it) {
        key = it->first;
        val = it->second;
        if (map.find(key) == map.end()) map.emplace(key, val);
        else map.at(key) = r(map.at(key), val);
      }
  
      newend = ::std::copy(map.begin(), map.end(), newstart);

      send_counts[i] = map.size();

      start = end;
      newstart = newend;
    }

    // compact.
    input.erase(newend, input.end());
  }



  /// modified version of mxx all2all, to replace content of buffer, and to return the recv counts.
  /// also, determine whether to use the large data version of mxx all2all at runtime.
  template<typename T, typename count_t = int>
  std::vector<count_t> all2all(std::vector<T>& buffer, const std::vector<count_t>& send_counts, MPI_Comm comm = MPI_COMM_WORLD)
  {
    //== first determine the max chunk to be sent around.  it should be smaller than max_int / p.
    //== otherwise there may be intermediate or final buffer larger than max_int.
    count_t max_count = *(::std::max_element(send_counts.begin(), send_counts.end()));
    int p;
    MPI_Comm_size(comm, &p);
    bool smaller = (max_count < (::std::numeric_limits<int>::max() / p));
    smaller = mxx::test_all(smaller, comm);

    if (!smaller && ::std::is_same<count_t, int>::value)
      // more data than datatype or api can SAFELY handle.  throw error.
      throw std::logic_error("Specified send count using int, but data size is greater than (max_int / p).  please convert to using size_t.");


    if (smaller && ::std::is_same<count_t, size_t>::value) {
      // convert count_t to int for send_count.  uses a2a.

      std::vector<int> sc(send_counts.size());
      ::std::transform(send_counts.begin(), send_counts.end(), sc.begin(), [](count_t const &x) { return x; });

      // get counts and displacements
      std::vector<int> rc = mxx::all2all(sc, 1, comm);
      // get total size
      std::size_t recv_size = ::std::accumulate(rc.begin(), rc.end(), 0);
      std::vector<T> recv_buffer(recv_size);
      // all-2-all  - automatically chooses between isend/irecv version or collective version based on type of count_t
      ::mxx::all2all(buffer.begin(), recv_buffer.begin(), sc, rc, comm);

      // return the buffer
      buffer.swap(recv_buffer);

      // convert recv_count from int to count_t
      std::vector<count_t> recv_counts(rc.size());
      ::std::transform(rc.begin(), rc.end(), recv_counts.begin(), [](int const &x) { return x; });

      // return the receive buffer
      return recv_counts;

    } else {  // else smaller & using int, or larger & using size_t.  use mxx api as is.


      // get counts and displacements
      std::vector<count_t> recv_counts = mxx::all2all(send_counts, 1, comm);
      // get total size
      std::size_t recv_size = ::std::accumulate(recv_counts.begin(), recv_counts.end(), 0);
      std::vector<T> recv_buffer(recv_size);
      // all-2-all  - automatically chooses between isend/irecv version or collective version based on type of count_t
      ::mxx::all2all(buffer.begin(), recv_buffer.begin(), send_counts, recv_counts, comm);

      // return the buffer
      buffer.swap(recv_buffer);
      // return the receive buffer
      return recv_counts;
    }
  }

//
//  /// send only unique elements.  bucketing via in-place bucketing, and unique is computed per bucket.
//  template<typename T, typename _TargetP, typename count_t = int>
//  std::vector<count_t> all2all_unique(std::vector<T>& buffer, _TargetP target_p_fun, MPI_Comm comm)
//  {
//      // get comm parameters
//      int p, rank;
//      MPI_Comm_size(comm, &p);
//      MPI_Comm_rank(comm, &rank);
//
//      // bucket input by their target processor
//      ::std::vector<count_t> send_counts = ::mxx2::bucketing(msgs, target_p_fun, p);
//
//      // retain unique and change send_counts.  note that inplace bucketing does not allow
//      retain_unique<>(msgs, send_counts);
//
//      // send and return the recv count
//      return mxx2::all2all(msgs, send_counts, comm);
//  }
//
//
//  /// send only unique elements.  bucketing via in-place bucketing, and unique is computed per bucket.
//  template<typename T, typename _TargetP, typename count_t = int>
//  std::vector<count_t> all2all_reduce(std::vector<T>& msgs, _TargetP target_p_fun, MPI_Comm comm)
//  {
//      // get comm parameters
//      int p, rank;
//      MPI_Comm_size(comm, &p);
//      MPI_Comm_rank(comm, &rank);
//
//      // bucket input by their target processor
//      ::std::vector<count_t> send_counts = ::mxx2::bucketing(msgs, target_p_fun, p);
//
//      // retain unique and change send_counts.  note that inplace bucketing does not allow
//      retain_unique(msgs, send_counts);
//
//      // send and return the recv count
//      return mxx2::all2all(msgs, send_counts, comm);
//  }


/* deprecated.  use bucketing + all2all in this file instead - more flexible and allows large data (size > int)
  // ~ 5 to 15% faster compared to standard version, but requires more memory.
  template<typename T, typename _TargetP>
  std::vector<int> msgs_all2all(std::vector<T>& msgs, _TargetP target_p_fun, MPI_Comm comm)
  {
      // get comm parameters
      int p, rank;
      MPI_Comm_size(comm, &p);
      MPI_Comm_rank(comm, &rank);


      // bucket input by their target processor
      // TODO: in-place bucketing??
      std::vector<int> send_counts(p, 0);
      std::vector<int> pids(msgs.size());
      if (msgs.size() > 0) {
//        pids.resize(msgs.size());
        for (int i = 0; i < msgs.size(); ++i)
        {
          pids[i] = target_p_fun(msgs[i]);
          ++(send_counts[pids[i]]);
        }
      }

      // get all2all params
      std::vector<int> recv_counts = mxx::all2all(send_counts, 1, comm);
      std::vector<int> send_displs = mxx::get_displacements(send_counts);
      std::vector<int> recv_displs = mxx::get_displacements(recv_counts);

      std::vector<int> offset = send_displs;

      // copy.  need to be able to track current position within each block.
      std::vector<T> send_buffer(msgs.size());
      if (msgs.size() > 0) {  // still need to participate even if msg size is 0.
//        send_buffer.resize(msgs.size());
        for (int i = 0; i < msgs.size(); ++i)
        {
          send_buffer[(offset[pids[i]])++] = msgs[i];
        }
      }
      std::vector<int>().swap(pids);  // clear pids reallocating.

      // resize messages to fit recv
      std::size_t recv_size = recv_displs[p-1] + recv_counts[p-1];
      //msgs.clear();
      //msgs.shrink_to_fit();
      msgs.resize(recv_size);   // create a vector with exact expected size.  content will be overwritten during all2all.
          // if we clear first, then resize, then we are creating and inserting the entire recv_size.
          // calling resize only: if realloc, then create and insert recv_size.  else only the difference.
      //msgs = std::vector<T>(recv_size);

      // get MPI type
      mxx::datatype<T> dt;
      MPI_Datatype mpi_dt = dt.type();
      MPI_Aint s;
      MPI_Type_extent(mpi_dt, &s);
      if (s != sizeof(T)) printf("ERROR: local rank %d data size for type %lu, mpitype size %ld\n", rank, sizeof(T), s);

      // all2all
      MPI_Alltoallv(&send_buffer[0], &send_counts[0], &send_displs[0], mpi_dt,
                    &msgs[0], &recv_counts[0], &recv_displs[0], mpi_dt, comm);
      // done, result is returned in vector of input messages

      return recv_counts;
  }
  */

  template <typename T, typename Func>
  T reduce(T& x, Func func, MPI_Comm comm = MPI_COMM_WORLD, int root = 0)
  {
    // get user op
    MPI_Op op = ::mxx::create_user_op<T, Func>(func);
      // get type
      ::mxx::datatype<T> dt;
      // perform reduction
      T result;
      MPI_Reduce(&x, &result, 1, dt.type(), op, root, comm);
      // clean up op
      ::mxx::free_user_op<T>(op);
      // return result
      return result;
  }

  template <typename T, typename Func>
  ::std::vector<T> reduce(::std::vector<T> & v, Func func, MPI_Comm comm = MPI_COMM_WORLD, int root = 0)
  {
    // get user op
    MPI_Op op = ::mxx::create_user_op<T, Func>(func);
    // get type
    ::mxx::datatype<T> dt;
    // get the max size.
    size_t s = v.size();
    size_t min_s = ::mxx::allreduce(s, [](size_t const &x, size_t const &y) { return ::std::min(x, y); }, comm);

    if (s > min_s) {
      int rank;
      MPI_Comm_rank(comm, &rank);
      WARNINGF("WARNING: process %d has larger vector than the min size during vector reduction.", rank);
    }

    // perform reduction
    ::std::vector<T> results(v.size());
    MPI_Reduce(&(v[0]), &(results[0]), min_s, dt.type(), op, root, comm);
    // clean up op
    ::mxx::free_user_op<T>(op);
    // return result
    return results;
  }


  template <typename T, typename Func>
  T reduce_scatter(::std::vector<T> & v, Func func, MPI_Comm comm = MPI_COMM_WORLD)
  {
    // get user op
    MPI_Op op = ::mxx::create_user_op<T, Func>(func);
    // get type
    ::mxx::datatype<T> dt;

    int p;
    MPI_Comm_size(comm, &p);

    if (v.size() != p) {
      int rank;
      MPI_Comm_rank(comm, &rank);
      WARNINGF("WARNING: process %d has vector size %lu that is not same as communicator size %d.", rank, v.size(), p);
    }

    // perform reduction
    T result;
    MPI_Reduce_scatter_block(&(v[0]), &result, 1, dt.type(), op, comm);
    // clean up op
    ::mxx::free_user_op<T>(op);
    // return result
    return result;
  }

}

#endif /* SRC_IO_MXX_SUPPORT_HPP_ */
