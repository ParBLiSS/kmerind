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
#include "mxx/shift.hpp"

#include "common/kmer.hpp"
#include "io/fastq_loader.hpp"

#include "utils/logging.h"

//#include "utils/system_utils.hpp"

#include <algorithm>  // for std::min

#include <unistd.h>   // for gethostname

#include "utils/container_traits.hpp"
#include "containers/container_utils.hpp"

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
    public datatype<decltype(bliss::io::FASTQ::SequenceId::file_pos)> {};

  template <typename T>
  class datatype<bliss::partition::range<T> > :
    public datatype_contiguous<T, sizeof(bliss::partition::range<T>) / sizeof(T)> {};

  template <typename T>
  class datatype<const bliss::partition::range<T> > :
    public datatype_contiguous<T, sizeof(bliss::partition::range<T>) / sizeof(T)> {};


}

template <typename T1, typename T2>
std::ostream &operator<<(std::ostream &os, std::pair<T1, T2> const &t) {
    return os << t.first << ":" << static_cast<uint64_t>(t.second);
}



namespace mxx2 {

  template <typename T>
  std::vector<T> right_shift(std::vector<T> & t, MPI_Comm comm = MPI_COMM_WORLD)
  {
      // get datatype
      mxx::datatype<T> dt;
      MPI_Datatype mpi_dt = dt.type();

      // get communication parameters
      // TODO: mxx::comm or boost::mpi::comm
      int p, rank;
      MPI_Comm_size(comm, &p);
      MPI_Comm_rank(comm, &rank);

      // TODO: handle tags with MXX (get unique tag function)
      int tag = 13;

      std::vector<T> left_value(t.size()); // result is the value that lies on the previous processor
      MPI_Request recv_req;
      if (rank > 0) // if not last processor
      {
          MPI_Irecv(&(left_value[0]), t.size(), mpi_dt, rank-1, tag,
                    comm, &recv_req);
      }
      if (rank < p-1) // if not first processor
      {
          // send my most right element to the right
          MPI_Send(&(t[0]), t.size(), mpi_dt, rank+1, tag, comm);
      }
      if (rank > 0)
      {
          // wait for the async receive to finish
          MPI_Wait(&recv_req, MPI_STATUS_IGNORE);
      }
      return left_value;
  }

  template <typename T>
  std::vector<T> left_shift(std::vector<T>& t, MPI_Comm comm = MPI_COMM_WORLD)
  {
      // get datatype
      mxx::datatype<T> dt;
      MPI_Datatype mpi_dt = dt.type();

      // get communication parameters
      // TODO: mxx comm or boost::mpi::comm
      int p, rank;
      MPI_Comm_size(comm, &p);
      MPI_Comm_rank(comm, &rank);

      // TODO: handle tags with MXX (get unique tag function)
      int tag = 15;

      std::vector<T> right_value(t.size()); // result is the value that lies on the previous processor
      MPI_Request recv_req;
      if (rank < p-1) // if not last processor
      {
          MPI_Irecv(&(right_value[0]), t.size(), mpi_dt, rank+1, tag,
                    comm, &recv_req);
      }
      if (rank > 0) // if not first processor
      {
          // send my most right element to the right
          MPI_Send(&(t[0]), t.size(), mpi_dt, rank-1, tag, comm);
      }
      if (rank < p-1)
      {
          // wait for the async receive to finish
          MPI_Wait(&recv_req, MPI_STATUS_IGNORE);
      }
      return right_value;
  }


#if !defined(HOST_NAME_MAX)
#define HOST_NAME_MAX 256
#endif

  /**
   * @brief     create a parent communicator and a child communicator, grouped by host names
   * @details   first we gather all the host names, and then 1 process will sort and pick unique host/rank pairs.
   *            the master than scatter, for each node, 2 colors, the first for the parent communicator, and second for the child communicator
   *
   *            all processes then participate to create the communicators.  for the non-selected processes (for the parent communicator), a null communicator is created.
   *
   * @param comm
   * @return    pair.  the first is a communicator with 1 proc per host machine.  the second is a communicator for all procs for that host machine.
   */
  ::std::pair<MPI_Comm, MPI_Comm> split_communicator_by_host(MPI_Comm comm = MPI_COMM_WORLD) {

    ::std::vector<char> host(HOST_NAME_MAX, 0);
    gethostname(&(host[0]), HOST_NAME_MAX);

    int rank;
    int p;
    MPI_Comm_size(comm, &p);
    MPI_Comm_rank(comm, &rank);


    ::std::vector<char> hosts = ::mxx::gather_vectors(host, comm);
    using ElemType = ::std::tuple<int, ::std::vector<char>, int > ;
    ::std::vector<ElemType> host_strs;

    ::std::vector<int > colors;
    int color;

    if (rank == 0) {
      // convert to strings
      auto b = hosts.begin();
      auto e = b;
      for (int i = 0; i < p; ++i) {
        e += HOST_NAME_MAX;
        host_strs.emplace_back(i, ::std::vector<char>(b, e), MPI_UNDEFINED);
        b = e;
      }
      // sort the rank-hostname tuples.
      ::std::sort(host_strs.begin(), host_strs.end(), [](ElemType const &x,
                                                         ElemType const &y) {
        int same = strncmp(&((::std::get<1>(x))[0]), &((::std::get<1>(y))[0]), HOST_NAME_MAX);
        if (same < 0) return true;
        else if ((same == 0) && (::std::get<0>(x) < ::std::get<0>(y))) return true;
        else return false;
      });


      // identify the one for each host name with the lowest rank and mark them
      color = 0;
      ::std::get<2>(host_strs[0]) = color;
      int same;
      for (int i = 1; i < p; ++i) {
        same = strncmp(&((::std::get<1>(host_strs[i-1]))[0]), &((::std::get<1>(host_strs[i]))[0]), HOST_NAME_MAX);
        if (same < 0) {
          ++color;
        }
        ::std::get<2>(host_strs[i]) = color;
      }

      // for each host name, color the result.
      // print them for debug
      DEBUGF("sorted hostnames:\n");
      for (int i = 0; i < p; ++i) {
        DEBUGF("%d, %s, %d\n", ::std::get<0>(host_strs[i]), &((::std::get<1>(host_strs[i]))[0]), ::std::get<2>(host_strs[i]));
      }


      // unsort
      ::std::sort(host_strs.begin(), host_strs.end(), [](ElemType const &x,
                                                         ElemType const &y) {
        return ::std::get<0>(x) < ::std::get<0>(y);
      });

      // copy to output vector.
      for (int i = 0; i < p; ++i) {
        colors.emplace_back(::std::get<2>(host_strs[i]));
      }
    }

    // scatter colors.
    color = MPI_UNDEFINED;
    MPI_Scatter(&(colors[0]), 1, MPI_INT, &color, 1, MPI_INT, 0, comm);

    // now create the sub communicators for each host's processes
    MPI_Comm group;
    MPI_Comm_split(comm, color, rank, &group);

    int group_rank;
    int group_size;
    MPI_Comm_rank(group, &group_rank);
    MPI_Comm_size(group, &group_size);

    // for the lowest rank process in the subcommunicator, create a spannng communicator.
    if (group_rank == 0) {
      color = 0;
    } else color = MPI_UNDEFINED;


    MPI_Comm group_leaders;
    MPI_Comm_split(comm, color, rank, &group_leaders);

    int group_leader_rank = MPI_UNDEFINED;
    int group_leader_size = MPI_UNDEFINED;
    if (group_leaders != MPI_COMM_NULL) {
      MPI_Comm_rank(group_leaders, &group_leader_rank);
      MPI_Comm_size(group_leaders, &group_leader_size);
    }

    INFOF("world %d/%d, host group %d/%d, host leaders %d/%d\n", rank, p, group_rank, group_size, group_leader_rank, group_leader_size);


    return ::std::make_pair(group_leaders, group);
  }

  // function returns color to split by.
  template <typename Func>
  void split_communicator_by_function(MPI_Comm comm, Func f, MPI_Comm & new_comm) {
    int rank;
    MPI_Comm_rank(comm, &rank);
    int color = f();
    MPI_Comm_split(comm, color, rank, &new_comm);
  }


  /// in place bucketing.  expect buffer and map both to be sorted by key.
  /// (in fact, no need to rebucket. just compute the send_counts.)
  template<typename count_t, typename T, typename Key, typename Comp>
  std::vector<count_t> bucketing(std::vector<T>& buffer, ::std::vector<::std::pair<Key, int> > map, Comp const & comp) {

    //== track the send count
    std::vector<count_t> send_counts(map.size() + 1, 0);

    if (buffer.size() == 0) return send_counts;

    //== since both sorted, search map entry in buffer - mlog(n).  the otherway around is nlog(m).
    // note that map contains splitters.  range defined as [map[i], map[i+1]), so when searching in buffer using entry in map, search via lower_bound
    auto b = buffer.begin();
    auto e = b;

    // if mlog(n) is larger than n, then linear would be faster
    bool linear = (map.size() * std::log(buffer.size()) ) > buffer.size();

    if (linear)
      for (int i = 0; i < map.size(); ++i) {
        e = ::fsc::lower_bound<true>(b, buffer.end(), map[i].first, comp);
        send_counts[i] = ::std::distance(b, e);
        b = e;
      }
    else
      for (int i = 0; i < map.size(); ++i) {
        e = ::fsc::lower_bound<false>(b, buffer.end(), map[i].first, comp);
        send_counts[i] = ::std::distance(b, e);
        b = e;
      }

    // last 1
    send_counts[map.size()] = ::std::distance(b, buffer.end());

    return send_counts;
  }


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
  std::vector<count_t> all2all_nocast(std::vector<T>& buffer, const std::vector<count_t>& send_counts, MPI_Comm comm = MPI_COMM_WORLD)
  {
      // get counts and displacements
      std::vector<count_t> recv_counts = mxx::all2all(send_counts, 1, comm);
      // get total size
      std::size_t recv_size = 0;
      for (int i = 0; i < recv_counts.size(); ++i) {
        recv_size += recv_counts[i];
      }
      std::vector<T> recv_buffer(recv_size);
      // all-2-all  - automatically chooses between isend/irecv version or collective version based on type of count_t
      ::mxx::all2all(buffer.begin(), recv_buffer.begin(), send_counts, recv_counts, comm);

      // return the buffer
      buffer.swap(recv_buffer);
      // return the receive buffer
      return recv_counts;
  }


  /// modified version of mxx all2all, to replace content of buffer, and to return the recv counts.
  /// also, determine whether to use the large data version of mxx all2all at runtime.
  template<typename internal_count_t, typename T, typename count_t = int>
  std::vector<count_t> all2all_cast(std::vector<T>& buffer, const std::vector<count_t>& send_counts, MPI_Comm comm = MPI_COMM_WORLD)
  {
      std::vector<internal_count_t> sc(send_counts.size());
      ::std::transform(send_counts.begin(), send_counts.end(), sc.begin(), [](count_t const &x) { return static_cast<internal_count_t>(x); });

      // get counts and displacements
      std::vector<internal_count_t> rc = ::mxx2::all2all_nocast(buffer, sc, comm);

      // convert recv_count from int to count_t
      std::vector<count_t> recv_counts(rc.size());
      ::std::transform(rc.begin(), rc.end(), recv_counts.begin(), [](int const &x) { return static_cast<count_t>(x); });

      // return the receive buffer
      return recv_counts;
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

    if (smaller && ::std::is_same<count_t, size_t>::value)
      return ::mxx2::all2all_cast<int>(buffer, send_counts, comm);
    else if (smaller && ::std::is_same<count_t, int>::value)
      return ::mxx2::all2all_nocast(buffer, send_counts, comm);
    else if (!smaller && ::std::is_same<count_t, size_t>::value)
      return ::mxx2::all2all_nocast(buffer, send_counts, comm);
    else
      return ::mxx2::all2all_cast<size_t>(buffer, send_counts, comm);
  }




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



  /// gathers n elements from each process to the root
  template <typename T>
  ::std::vector<T> gathern(::std::vector<T> & v, MPI_Comm comm = MPI_COMM_WORLD, int root = 0) {
    // get data type
    ::mxx::datatype<T> dt;

    // allocate result
    int p;
    int rank;
    MPI_Comm_size(comm, &p);
    MPI_Comm_rank(comm, &rank);

    // get the max size.
    size_t count = v.size();
    size_t min_s = ::mxx::allreduce(count, [](size_t const &x, size_t const &y) { return ::std::min(x, y); }, comm);

    if (count > min_s) {
      WARNINGF("WARNING: gathern process %d trying to send %lu, more than min of %lu.", rank, count, min_s);
    }
    if (min_s > (::std::numeric_limits<int>::max() / p)) {
      WARNINGF("WARNING: gathern process %d trying to send min of %lu, more than max_int/p.", rank, min_s);
    }

    ::std::vector<T> results;

    if (rank == root) {
      results.resize(p * min_s);
    }

    // perform gather.  does not explicitly support large count.
    MPI_Gather(&(v[0]), min_s, dt.type(), &(results[0]), min_s, dt.type(), root, comm);
    // return result
    return results;
  }

  /// scatters n elements from root to each process
  template <typename T>
  ::std::vector<T> scattern(::std::vector<T> & v, MPI_Comm comm = MPI_COMM_WORLD, int root = 0) {
    // get data type
    ::mxx::datatype<T> dt;

    // allocate result
    int p;
    MPI_Comm_size(comm, &p);

    size_t count = v.size() / p;
    ::mxx::datatype<size_t> dt2;
    MPI_Bcast(&count, 1, dt2.type(), root, comm);

    if (count == 0)
      WARNINGF("ERROR: mxx2::scattern called with %lu elements, not enough elements to sent to communicator size %d", v.size(), p);
    if ((v.size() %p) != 0)
      WARNINGF("WARNING: mxx2::scattern called with %lu elements, some elements not sent %lu", v.size(), v.size() % p);

    ::std::vector<T> results(count);

    // perform gather  - does not explicitly support large count.
    MPI_Scatter(&(v[0]), count, dt.type(), &(results[0]), count, dt.type(), root, comm);
    // return result
    return results;
  }

  /// default gather (1)
  template <typename T>
  ::std::vector<T> gather(T &v, MPI_Comm comm = MPI_COMM_WORLD, int root = 0) {
    // get data type
    ::mxx::datatype<T> dt;

    // allocate result
    int p;
    int rank;
    MPI_Comm_size(comm, &p);
    MPI_Comm_rank(comm, &rank);

    ::std::vector<T> results;

    if (rank == root) {
      results.resize(p);
    }

    // perform gather.  does not explicitly support large count.
    MPI_Gather(&v, 1, dt.type(), &(results[0]), 1, dt.type(), root, comm);
    // return result
    return results;
  }

  /// default scatter (1)
  template <typename T>
  T scatter(::std::vector<T>& v, MPI_Comm comm = MPI_COMM_WORLD, int root = 0) {
    // get data type
    ::mxx::datatype<T> dt;

    // allocate result
    int p;
    int rank;
    MPI_Comm_size(comm, &p);
    MPI_Comm_rank(comm, &rank);

    if (rank == root && v.size() < p)
      WARNINGF("ERROR: mxx2::scatter called with %lu elements, less than p (%d)", v.size(), p);

    T results;

    // perform gather  - does not explicitly support large count.
    MPI_Scatter(&(v[0]), 1, dt.type(), &results, 1, dt.type(), root, comm);
    // return result
    return results;
  }

  /// scatter with variable number of elements.  integer send_counts
  template <typename Iterator, typename T = typename ::std::iterator_traits<Iterator>::value_type>
  ::std::vector<T> 
  scatterv(Iterator begin, Iterator end, ::std::vector<int>& send_counts, MPI_Comm comm = MPI_COMM_WORLD, int root = 0) {

    // now determine if we should send via scatterv or isend/irecv
    int rank;
    int p;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &p);

    // get the amount to receive.
    int count;
    ::mxx::datatype<int> count_dt;
    MPI_Scatter(&(send_counts[0]), 1, count_dt.type(), &count, 1, count_dt.type(), root, comm);

    // get the data's type
    ::mxx::datatype<T> dt;

    // copy the input if not contiguous
    T* src = nullptr;
    ::std::vector<T> source;
    if (::bliss::utils::is_contiguous(begin))
      src = &(*begin);
    else {
      if (rank == root) source.insert(source.end(), begin, end);
      src = &(source[0]);
    }

    // allocate the output
    ::std::vector<T> results;
    results.resize(count);

    // displacements in source.
    auto send_displ = ::mxx::get_displacements(send_counts);
   

    size_t total = std::accumulate(send_counts.begin(), send_counts.end(), 0UL, std::plus<size_t>());
    mxx::datatype<size_t> size_dt;
    MPI_Bcast(&total, 1, size_dt.type(), root, comm);

    bool smaller = (total <= static_cast<size_t>(::std::numeric_limits<int>::max()));

    if (smaller) {  // total amount to send is smaller than max int, so can use scatterv.
        // using int.  no conversion needed.

      MPI_Scatterv(src, &(send_counts[0]), &(send_displ[0]), dt.type(), &(results[0]), count, dt.type(), root, comm);
    } else {  // not smaller.  need isend/irecv.

      // check to see if each isend has fewer than max_int number of elements.
      size_t max_count = *(std::max_element(send_counts.begin(), send_counts.end()));

      MPI_Bcast(&max_count, 1, size_dt.type(), root, comm);

      if (max_count > static_cast<size_t>(::std::numeric_limits<int>::max()))
        throw ::std::logic_error("trying to send more than max_int number of elements.  implementation to handle this is not yet in place");


      // set up the requests
      ::std::vector<MPI_Request> reqs((rank == root) ? (p+1) : 1);

      // everyone receives
      MPI_Irecv(&(results[0]), count, dt.type(), root, rank, comm, &(reqs[0]));

      // root send
      if (rank == root) {
        for (int i = 0; i < p; ++i) {
          MPI_Isend(src + send_displ[i], send_counts[i], dt.type(), i, i, comm, &(reqs[i + 1]));
        }
      }

      MPI_Waitall((rank == root) ? (p+1) : 1, &(reqs[0]), MPI_STATUSES_IGNORE);
    }


    return results;
  }

  /// scatter with variable number of elements.  size_t send_counts
  template <typename Iterator, typename T = typename ::std::iterator_traits<Iterator>::value_type>
  ::std::vector<T>
  scatterv(Iterator begin, Iterator end, ::std::vector<size_t>& send_counts, MPI_Comm comm = MPI_COMM_WORLD, int root = 0) {
    // now determine if we should send via scatterv or isend/irecv
    int rank;
    int p;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &p);

    // get the amount to receive.
    size_t count;
    ::mxx::datatype<size_t> count_dt;
    MPI_Scatter(&(send_counts[0]), 1, count_dt.type(), &count, 1, count_dt.type(), root, comm);

    // get the data's type
    ::mxx::datatype<T> dt;

    // allocate the output
    ::std::vector<T> results;
    results.resize(count);

    // copy the input if root
    T* src = nullptr;
    ::std::vector<T> source;
    if (::bliss::utils::is_contiguous(begin))
      src = &(*begin);
    else {
      if (rank == root) source.insert(source.end(), begin, end);
      src = &(source[0]);
    }

    // displacements in source.
    std::vector<size_t> send_displ = ::mxx::get_displacements(send_counts);

    size_t total = std::accumulate(send_counts.begin(), send_counts.end(), 0UL, std::plus<size_t>());

    MPI_Bcast(&total, 1, count_dt.type(), root, comm);

    bool smaller = (total <= static_cast<size_t>(::std::numeric_limits<int>::max()));

    if (smaller) {  // total amount to send is smaller than max int, so can use scatterv.
        // convert to int
        ::std::vector<int> sc(send_counts.size());
        ::std::transform(send_counts.begin(), send_counts.end(), sc.begin(), [](size_t const& x) { return static_cast<int>(x); });
      
        // get displacements
        ::std::vector<int> sd(send_counts.size());
        ::std::transform(send_displ.begin(), send_displ.end(), sd.begin(), [](size_t const& x) { return static_cast<int>(x); });

        MPI_Scatterv(src, &(sc[0]), &(sd[0]), dt.type(), &(results[0]), count, dt.type(), root, comm);
    } else {  // not smaller.  need isend/irecv.
      INFOF("NOTE: scatterv using isend version.");
      // check to see if each isend has fewer than max_int number of elements.
      size_t max_count = *(std::max_element(send_counts.begin(), send_counts.end()));
      MPI_Bcast(&max_count, 1, count_dt.type(), root, comm);

      if (max_count > static_cast<size_t>(::std::numeric_limits<int>::max()))
        throw ::std::logic_error("trying to send more than max_int number of elements.  implementation to handle this is not yet in place");


      // set up the requests
      ::std::vector<MPI_Request> reqs((rank == root) ? (p+1) : 1);

      // everyone receives
      MPI_Irecv(&(results[0]), count, dt.type(), root, rank, comm, &(reqs[0]));

      // root sends
      if (rank == root) {
        for (int i = 0; i < p; ++i) {
          MPI_Isend(src + send_displ[i], send_counts[i], dt.type(), i, i, comm, &(reqs[i + 1]));
        }
      }

      MPI_Waitall((rank == root) ? (p+1) : 1, &(reqs[0]), MPI_STATUSES_IGNORE);
    }
  
  
    return results;
  }



  struct scan_op {

      /// MPI Scan
      template <typename T, typename Func = std::plus<T>, bool exclusive = false >
      static T scan(T & x, Func f, MPI_Comm comm = MPI_COMM_WORLD) {
        MPI_Op op = mxx::create_user_op<T, Func >(f);
        mxx::datatype<T> dt;
        T result;
        if (exclusive)
          MPI_Exscan(&x, &result, 1, dt.type(), op, comm );
        else
          MPI_Scan(&x, &result, 1, dt.type(), op, comm );
        mxx::free_user_op<T>(op);
        return result;
      }

      template <typename T, typename Func = std::plus<T>, bool exclusive = false  >
      static T scan_reverse(T & x, Func f, MPI_Comm comm = MPI_COMM_WORLD) {
        MPI_Comm rev_comm;
        mxx::rev_comm(comm, rev_comm);
        auto result = scan<T, Func, exclusive>(x, f, rev_comm);
        MPI_Comm_free(&rev_comm);
        return result;
      }

      template <typename T, typename Func = std::plus<T> >
      static T exscan(T & x, Func f, MPI_Comm comm = MPI_COMM_WORLD) {
        return scan<T, Func, true>(x, f, comm);
      }

      template <typename T, typename Func = std::plus<T> >
      static T exscan_reverse(T & x, Func f, MPI_Comm comm = MPI_COMM_WORLD) {
        return scan_reverse<T, Func, true>(x, f, comm);
      }

      template <typename T, typename Func = std::plus<T>, bool exclusive = false >
      static std::vector<T> scan_elementwise(std::vector<T> & x, Func f, MPI_Comm comm = MPI_COMM_WORLD) {
        size_t s = x.size();
        s = mxx::allreduce(s, [](size_t &x , size_t & y){ return ::std::min(x, y); }, comm);
        if (s < x.size()) {
          int rank;
          MPI_Comm_rank(comm, &rank);
          WARNINGF("WARNING: Rank %d  %lu elements out of %lu used in scan.", rank, s, x.size());
        }

        MPI_Op op = mxx::create_user_op<T, Func >(f);
        mxx::datatype<T> dt;
        std::vector<T> results(s);
        if (exclusive)
          MPI_Exscan(&(x[0]), &(results[0]), s, dt.type(), op, comm );
        else
          MPI_Scan(&(x[0]), &(results[0]), s, dt.type(), op, comm );
        mxx::free_user_op<T>(op);
        return results;
      }

      template <typename T, typename Func = std::plus<T>, bool exclusive = false >
      static std::vector<T> scan_reverse_elementwise(std::vector<T> & x, Func f, MPI_Comm comm = MPI_COMM_WORLD) {
        MPI_Comm rev_comm;
        mxx::rev_comm(comm, rev_comm);
        auto results = scan_elementwise<T, Func, exclusive>(x, f, rev_comm);
        MPI_Comm_free(&rev_comm);
        return results;
      }

      template <typename T, typename Func = std::plus<T> >
      static std::vector<T> exscan_elementwise(std::vector<T> & x, Func f, MPI_Comm comm = MPI_COMM_WORLD) {
        return scan_elementwise<T, Func, true>(x, f, comm);
      }

      template <typename T, typename Func = std::plus<T> >
      static std::vector<T> exscan_reverse_elementwise(std::vector<T> & x, Func f, MPI_Comm comm = MPI_COMM_WORLD) {
        return scan_reverse_elementwise<T, Func, true>(x, f, comm);
      }


      template <typename T, typename Func = std::plus<T> >
      static std::vector<T> scanv(std::vector<T> & x, Func f, MPI_Comm comm = MPI_COMM_WORLD) {

        std::vector<T> results(x.size());

        // create subcommuinicator with processors that have non-zero values.  this avoids special cases.
        MPI_Comm in_comm;
        ::mxx2::split_communicator_by_function(comm, [&x](){ return (x.size() > 0) ? 1 : MPI_UNDEFINED; }, in_comm);

        if (in_comm != MPI_COMM_NULL) {  // only do this for valid comm (i.e. x.size > 0)

          int rank, p;
          MPI_Comm_rank(in_comm, &rank);
          MPI_Comm_size(in_comm, &p);

          // first scan locally.
          std::partial_sum(x.begin(), x.end(), results.begin(), f);

          // then do global exscan
          T last = results.back();
          last = exscan(last, f, in_comm);

          // finally update local values
          if (rank > 0) {  // rank 0 does not participate
            // note that the order of argument when calling f matters for segmented scan - not communitive
            std::for_each(results.begin(), results.end(), [&last, &f](T & x) { x = f(last, x); });
          }

          MPI_Comm_free(&in_comm);
        }

        return results;
      }

      template <typename T, typename Func = std::plus<T> >
      static std::vector<T> scanv_reverse(std::vector<T> & x, Func f, MPI_Comm comm = MPI_COMM_WORLD) {
        std::vector<T> results(x.size());

        // create subcommuinicator with processors that have non-zero values.  this avoids special cases.
        MPI_Comm in_comm;
        ::mxx2::split_communicator_by_function(comm, [&x](){ return (x.size() > 0) ? 1 : MPI_UNDEFINED; }, in_comm);

        if (in_comm != MPI_COMM_NULL) {  // only do this for valid comm (i.e. x.size > 0)

          // first scan locally.
          std::partial_sum(x.rbegin(), x.rend(), results.rbegin(), f);

          // then do global exscan
          T first = results.front();
          first = exscan_reverse(first, f, in_comm);

          // finally update local values
          int rank, p;
          MPI_Comm_rank(in_comm, &rank);
          MPI_Comm_size(in_comm, &p);
          if (rank < p - 1) {  // last rank does not participate
            // note that the order of argument when calling f matters for segmented scan - not communitive
            std::for_each(results.begin(), results.end(), [&first, &f](T & x) { x = f(first, x); });
          }

          MPI_Comm_free(&in_comm);
        }

        return results;
      }


      template <typename T, typename Func = std::plus<T> >
      static std::vector<T> exscanv(std::vector<T> & x, Func f, MPI_Comm comm = MPI_COMM_WORLD) {
        std::vector<T> results(x.size());

        // create subcommuinicator with processors that have non-zero values.  this avoids special cases.
        MPI_Comm in_comm;
        ::mxx2::split_communicator_by_function(comm, [&x](){ return (x.size() > 0) ? 1 : MPI_UNDEFINED; }, in_comm);

        if (in_comm != MPI_COMM_NULL) {  // only do this for valid comm (i.e. x.size > 0)
          int rank, p;
          MPI_Comm_rank(in_comm, &rank);
          MPI_Comm_size(in_comm, &p);


//          for (int i = 0; i < p; ++i) {
//            std::cout << "R " << rank << "," << i << " scan input " << x[i] << std::endl;
//          }

          // first scan locally.  scan 0 to m-1 and store at 1 to m
          std::partial_sum(x.begin(), x.end() - 1, results.begin() + 1, f);

//          for (int i = 0; i < p; ++i) {
//            std::cout << "R " << rank << "," << i << " local prefix " << results[i] << std::endl;
//          }


          // then do global exscan with result of local _SCAN_ (not exscan)
          // note that the order of argument when calling f matters for segmented scan - not communitive
          T last = (x.size() == 1) ? x.back() : f(results.back(), x.back());  // special case if length == 1
          last = exscan(last, f, in_comm);
//          std::cout << "R " << rank << "," << " last " << last << std::endl;

          // finally update local values
          if (rank > 0) {  // rank 0 does not participate
            results.front() = last;  // exclusive scan takes value from last proc
            // note that the order of argument when calling f matters for segmented scan - not communitive
            std::for_each(results.begin() + 1, results.end(), [&last, &f](T & x) { x = f(last, x); });
          }

//          for (int i = 0; i < p; ++i) {
//            std::cout << "R " << rank << "," << i << " global prefix " << results[i] << std::endl;
//          }

          MPI_Comm_free(&in_comm);
        }

        return results;
      }

      template <typename T, typename Func = std::plus<T> >
      static std::vector<T> exscanv_reverse(std::vector<T> & x, Func f, MPI_Comm comm = MPI_COMM_WORLD) {
        std::vector<T> results(x.size());

        // create subcommuinicator with processors that have non-zero values.  this avoids special cases.
        MPI_Comm in_comm;
        ::mxx2::split_communicator_by_function(comm, [&x](){ return (x.size() > 0) ? 1 : MPI_UNDEFINED; }, in_comm);

        if (in_comm != MPI_COMM_NULL) {  // only do this for valid comm (i.e. x.size > 0)

          // first scan locally.  scan 0 to m-1 and store at 1 to m
          std::partial_sum(x.rbegin(), x.rend() - 1, results.rbegin() + 1, f);

          // then do global exscan with result of local _SCAN_ (not exscan)
          // note that the order of argument when calling f matters for segmented scan - not communitive
          T first = (x.size() == 1) ? x.front() : f(results.front(), x.front());  // special case if length == 1
          first = exscan_reverse(first, f, in_comm);

          // finally update local values
          int rank, p;
          MPI_Comm_rank(in_comm, &rank);
          MPI_Comm_size(in_comm, &p);
          if (rank < p - 1) {  // last rank does not participate
            results.back() = first;  // exclusive scan takes value from last proc
            // note that the order of argument when calling f matters for segmented scan - not communitive
            std::for_each(results.begin(), results.end() - 1, [&first, &f](T & x) { x = f(first, x); });
          }
          MPI_Comm_free(&in_comm);
        }

        return results;
      }

  };

  struct reduce_op {
      template <typename T, typename Func = std::plus<T> >
      static T reduce(T & x, Func f, MPI_Comm comm = MPI_COMM_WORLD, int root = 0) {
        MPI_Op op = mxx::create_user_op<T, Func >(f);
        mxx::datatype<T> dt;
        T result;
        MPI_Reduce(&x, &result, 1, dt.type(), op, root, comm );
        mxx::free_user_op<T>(op);
        return result;
      }
      template <typename T, typename Func = std::plus<T> >
      static std::vector<T> reduce_elementwise(std::vector<T> & x, Func f, MPI_Comm comm = MPI_COMM_WORLD, int root = 0) {
        size_t s = x.size();
        s = mxx::allreduce(s, [](size_t &x , size_t & y){ return ::std::min(x, y); }, comm);
        if (s < x.size()) {
          int rank;
          MPI_Comm_rank(comm, &rank);
          WARNINGF("WARNING: Rank %d  %lu elements out of %lu used in scan.", rank, s, x.size());
        }

        MPI_Op op = mxx::create_user_op<T, Func >(f);
        mxx::datatype<T> dt;
        std::vector<T> results(s);
        MPI_Reduce(&(x[0]), &(results[0]), s, dt.type(), op, root, comm );
        mxx::free_user_op<T>(op);
        return results;
      }
      template <typename T, typename Func = std::plus<T> >
      static T reducev(std::vector<T> & x, Func f, MPI_Comm comm = MPI_COMM_WORLD, int root = 0) {

        T result;

        // create subcommuinicator with processors that have non-zero values.  this avoids special cases.
        MPI_Comm in_comm;
        ::mxx2::split_communicator_by_function(comm, [&x](){ return (x.size() > 0) ? 1 : MPI_UNDEFINED; }, in_comm);

        int in_rank = MPI_UNDEFINED;
        if (in_comm != MPI_COMM_NULL) {  // only do this for valid comm (i.e. x.size > 0)

          // first reduce locally.
          result = std::accumulate(x.begin() + 1, x.end(), x.front(), f);

          // then do global reduce, directly to the root if possible.
          result = reduce_op::reduce(result, f, in_comm, 0);

          // assume rank 0 in in_comm is not root rank.
          MPI_Comm_rank(in_comm, &in_rank);
          //printf("in_rank == %d\n", in_rank);

          MPI_Comm_free(&in_comm);

        }

        int rank;
        MPI_Comm_rank(comm, &rank);

        MPI_Request reqs[2];

        ::mxx::datatype<T> dt;
        if (rank == root) {
          //printf("receiving on %d\n", root);
          MPI_Irecv(&result, 1, dt.type(), MPI_ANY_SOURCE, 770, comm,  reqs);
        }
        if (in_rank == 0) {
          //printf("sending from %d (in_rank == 0) \n", rank);
          MPI_Isend(&result, 1, dt.type(), root, 770, comm, reqs + 1);
        }

        if (in_rank == 0)
          MPI_Wait(reqs + 1, MPI_STATUS_IGNORE);

        if (rank == root)
          MPI_Wait(reqs, MPI_STATUS_IGNORE);

        return result;
      }


      /// reduce operation with minloc, elementwise
      template <typename T>
      static ::std::pair<T, int > reduce_loc(T & v, MPI_Op op, MPI_Comm comm = MPI_COMM_WORLD, int root = 0)
      {
        if ((op != MPI_MINLOC) && (op != MPI_MAXLOC)) throw std::logic_error("reduce_loc only support MPI_MINLOC and MPI_MAXLOC");


        // get the mpi type
        MPI_Datatype dt;
        if (::std::is_same<T, float>::value) dt = MPI_FLOAT_INT;
        else if (::std::is_same<T, double>::value) dt = MPI_DOUBLE_INT;
        else if (::std::is_same<T, long>::value) dt = MPI_LONG_INT;
        else if (::std::is_same<T, int>::value) dt = MPI_2INT;
        else if (::std::is_same<T, short>::value) dt = MPI_SHORT_INT;
        else if (::std::is_same<T, long double>::value) dt = MPI_LONG_DOUBLE_INT;
        else
            throw std::logic_error("type param T is not one of float, double, long, int, short, or long double.");


        // struct because pair may not have the right ordering in memory.
        struct vi {
            T val;
            int rank;
        };

        int rank;
        MPI_Comm_rank(comm, &rank);

        // set up the data structure
        vi temp = {v, rank};

        // perform reduction
        vi temp2;
        MPI_Reduce(&temp, &temp2, 1, dt, op, root, comm);

        // return result
        return ::std::make_pair(temp2.val, temp2.rank);
      }

      /// reduce operation with minloc, elementwise
      template <typename T>
      static ::std::pair<::std::vector<T>, ::std::vector<int> > reduce_loc_elementwise(::std::vector<T> & v, MPI_Op op, MPI_Comm comm = MPI_COMM_WORLD, int root = 0)
      {
        if ((op != MPI_MINLOC) && (op != MPI_MAXLOC)) throw std::logic_error("reduce_loc only support MPI_MINLOC and MPI_MAXLOC");


        // get the mpi type
        MPI_Datatype dt;
        if (::std::is_same<T, float>::value) dt = MPI_FLOAT_INT;
        else if (::std::is_same<T, double>::value) dt = MPI_DOUBLE_INT;
        else if (::std::is_same<T, long>::value) dt = MPI_LONG_INT;
        else if (::std::is_same<T, int>::value) dt = MPI_2INT;
        else if (::std::is_same<T, short>::value) dt = MPI_SHORT_INT;
        else if (::std::is_same<T, long double>::value) dt = MPI_LONG_DOUBLE_INT;
        else
            throw std::logic_error("type param T is not one of float, double, long, int, short, or long double.");


        // struct because pair may not have the right ordering in memory.
        struct vi {
            T val;
            int rank;
        };

        int rank;
        MPI_Comm_rank(comm, &rank);

        // get the max size.
        size_t s = v.size();
        size_t min_s = ::mxx::allreduce(s, [](size_t const &x, size_t const &y) { return ::std::min(x, y); }, comm);

        if (s > min_s) {
          WARNINGF("WARNING: process %d has larger vector than the min size during vector reduction.", rank);
        }


        // set up the data structure
        ::std::vector< vi > temp(s);
        for (int i = 0; i < s; ++i) {
          temp[i].val = v[i];
          temp[i].rank = rank;
        }

        // perform reduction
        ::std::vector< vi > temp2(s);
        MPI_Reduce(&(temp[0]), &(temp2[0]), min_s, dt, op, root, comm);

        // now extract the results
        ::std::vector<T> outVal(s);
        ::std::vector<int> outIdx(s);
        for (int i = 0; i < s; ++i) {
          outVal[i] = temp2[i].val;
          outIdx[i] = temp2[i].rank;
        }

        // return result
        return ::std::make_pair(outVal, outIdx);
      }

      template <typename T, typename Func = std::plus<T> >
      static T allreduce(T & x, Func f, MPI_Comm comm = MPI_COMM_WORLD) {
        MPI_Op op = mxx::create_user_op<T, Func >(f);
        mxx::datatype<T> dt;
        T result;
        MPI_Allreduce(&x, &result, 1, dt.type(), op, comm );
        mxx::free_user_op<T>(op);
        return result;
      }
      template <typename T, typename Func = std::plus<T> >
      static std::vector<T> allreduce_elementwise(std::vector<T> & x, Func f, MPI_Comm comm = MPI_COMM_WORLD) {
        size_t s = x.size();
        s = mxx::allreduce(s, [](size_t &x , size_t & y){ return ::std::min(x, y); }, comm);
        if (s < x.size()) {
          int rank;
          MPI_Comm_rank(comm, &rank);
          WARNINGF("WARNING: Rank %d  %lu elements out of %lu used in scan.", rank, s, x.size());
        }

        MPI_Op op = mxx::create_user_op<T, Func >(f);
        mxx::datatype<T> dt;
        std::vector<T> results(s);
        MPI_Allreduce(&(x[0]), &(results[0]), s, dt.type(), op, comm );
        mxx::free_user_op<T>(op);
        return results;
      }
      template <typename T, typename Func = std::plus<T> >
      static T allreducev(std::vector<T> & x, Func f, MPI_Comm comm = MPI_COMM_WORLD) {

        T result;

        // create subcommuinicator with processors that have non-zero values.  this avoids special cases.
        MPI_Comm in_comm;
        ::mxx2::split_communicator_by_function(comm, [&x](){ return (x.size() > 0) ? 1 : MPI_UNDEFINED; }, in_comm);

        int root = MPI_UNDEFINED;
        if (in_comm != MPI_COMM_NULL) {  // only do this for valid comm (i.e. x.size > 0)

          // first reduce locally.
          result = std::accumulate(x.begin() + 1, x.end(), x.front(), f);

          // then do global reduce.
          result = reduce_op::reduce(result, f, in_comm, 0);

          // mark the root
          int rank;
          MPI_Comm_rank(in_comm, &rank);
          if (rank == 0)
            MPI_Comm_rank(comm, &root);  // only 1 will have defined rank (as root)

          MPI_Comm_free(&in_comm);
        }

        // let every one know the rank of the root.
        if (root == MPI_UNDEFINED) root = -1;  // reset all others to -1.  root has a rank so >= 0
        root = reduce_op::allreduce(root, [](int & x, int & y){ return ::std::max(x, y); }, comm);

        // broadcast the results
        mxx::datatype<T> dt;
        //printf("broadcast from %d\n", root);
        MPI_Bcast(&result, 1, dt.type(), root, comm);

        return result;
      }

  };



  template <typename T>
  struct segment {

    static uint8_t is_start(T const & seg_id, MPI_Comm comm = MPI_COMM_WORLD) {
      // right shift
      T prev = mxx::right_shift(seg_id, comm);

      // compare equal
      int rank;
      MPI_Comm_rank(comm, &rank);
      return (rank == 0 || prev != seg_id) ? 1 : 0;
    }

    static std::vector<uint8_t> is_start_elementwise(std::vector<T> & seg_ids, MPI_Comm comm = MPI_COMM_WORLD) {
      // right shift
      std::vector<T> prev = mxx2::right_shift(seg_ids, comm);

      // compare equal
      int rank;
      MPI_Comm_rank(comm, &rank);

      std::vector<uint8_t> results(seg_ids.size());
      std::transform(seg_ids.begin(), seg_ids.end(), prev.begin(), results.begin(),
                      [&rank](T const &x, T const &y) { return (rank == 0 || x != y) ? 1 : 0; });

      return results;
    }


    static std::vector<uint8_t> is_start_v(std::vector<T> & seg_ids, MPI_Comm comm = MPI_COMM_WORLD) {
      // exclude empty procs
	  MPI_Comm work_comm;
	  mxx2::split_communicator_by_function(comm, [&seg_ids](){ return (seg_ids.size() == 0) ? MPI_UNDEFINED : 1; }, work_comm);

	  // compute all except first
	  std::vector<uint8_t> results(seg_ids.size());

	  if (seg_ids.size() > 0)
	    std::transform(seg_ids.begin(), seg_ids.end() - 1, seg_ids.begin() + 1, results.begin() + 1,
			  	  	  [](T const &x, T const &y) { return (x == y) ? 0 : 1; });

	  // now compute the first.

	  // right shift last
	  if (work_comm != MPI_COMM_NULL) {  // also means non-zero vector
		  T prev = mxx::right_shift(seg_ids.back(), work_comm);

		  // local compare equal first element
      int rank;
      MPI_Comm_rank(comm, &rank);
		  results[0] = (rank == 0 || prev != seg_ids.front()) ? 1 : 0;

      MPI_Comm_free(&work_comm);

	  }

      return results;
    }


    static uint8_t is_end(T const & seg_id, MPI_Comm comm = MPI_COMM_WORLD) {
      // (more complicated to use reverse communicator)
      // left shift
      T next = mxx::left_shift(seg_id, comm);

      // compare equal?
      int rank;
      MPI_Comm_rank(comm, &rank);
      int p;
      MPI_Comm_size(comm, &p);
      return (rank == p - 1 || next != seg_id) ? 1 : 0;
    }

    static std::vector<uint8_t> is_end_elementwise(std::vector<T> & seg_ids, MPI_Comm comm = MPI_COMM_WORLD) {
      // right shift
      std::vector<T> next = mxx2::left_shift(seg_ids, comm);

      // compare equal
      int rank;
      MPI_Comm_rank(comm, &rank);
      int p;
      MPI_Comm_size(comm, &p);

      std::vector<uint8_t> results(seg_ids.size());
      std::transform(seg_ids.begin(), seg_ids.end(), next.begin(), results.begin(),
                      [&rank, &p](T const &x, T const &y) { return (rank == p - 1 || x != y) ? 1 : 0; });

      return results;
    }


    static std::vector<uint8_t> is_end_v(std::vector<T> & seg_ids, MPI_Comm comm = MPI_COMM_WORLD) {
      // exclude empty procs
  	  MPI_Comm work_comm;
  	  mxx2::split_communicator_by_function(comm, [&seg_ids](){ return (seg_ids.size() == 0) ? MPI_UNDEFINED : 1; }, work_comm);

  	  // compute all except first
  	  std::vector<uint8_t> results(seg_ids.size());

  	  if (seg_ids.size() > 0)
  	    std::transform(seg_ids.begin()+1, seg_ids.end(), seg_ids.begin(), results.begin(),
  			  	  	  [](T const &x, T const &y) { return (x == y) ? 0 : 1; });

  	  // now compute the last.

  	  // left shift last
  	  if (work_comm != MPI_COMM_NULL) {  // also means non-zero vector
  		  T next = mxx::left_shift(seg_ids.front(), work_comm);

  		  // local compare equal first element
        int rank;
        MPI_Comm_rank(comm, &rank);
        int p;
        MPI_Comm_size(comm, &p);

  		  results.back() = (rank == p - 1 || next != seg_ids.back()) ? 1 : 0;

  		  MPI_Comm_free(&work_comm);
  	  }

        return results;
    }


    static T to_unique_segment_id_from_start(uint8_t const & start, uint8_t const & non_start_value = 0, MPI_Comm comm = MPI_COMM_WORLD) {
      // convert to 1 and 0.
    	T seg_id = (start == non_start_value) ? 0 : 1;

      // scan
    	return scan_op::scan(seg_id, std::plus<T>(), comm);
    }

    static std::vector<T> to_unique_segment_id_from_start_elementwise(std::vector<uint8_t> & starts, uint8_t const & non_start_value = 0, MPI_Comm comm = MPI_COMM_WORLD) {
      // convert to 1 and 0
      std::vector<T> seg_ids(starts.size());
      std::transform(starts.begin(), starts.end(), seg_ids.begin(), [&non_start_value](uint8_t & x){ return (x == non_start_value) ? 0 : 1; });

      // scan
      return scan_op::scan_elementwise(seg_ids, std::plus<T>(), comm);
    }


    static std::vector<T> to_unique_segment_id_from_start_v(std::vector<uint8_t> & starts, uint8_t const & non_start_value = 0, MPI_Comm comm = MPI_COMM_WORLD) {
      // convert to 1 and 0
      std::vector<T> seg_ids(starts.size());
      if (starts.size() > 0)
        std::transform(starts.begin(), starts.end(), seg_ids.begin(), [&non_start_value](uint8_t & x){ return (x == non_start_value) ? 0 : 1; });
      
      // scan
      return scan_op::scanv(seg_ids, std::plus<T>(), comm);
    }



    static T to_unique_segment_id_from_end(uint8_t const & end, uint8_t const & non_end_value = 0, MPI_Comm comm = MPI_COMM_WORLD) {
      // convert to 1 and 0.
    	T seg_id = (end == non_end_value) ? 0 : 1;

      // reverse scan
    	return scan_op::scan_reverse(seg_id, std::plus<T>(), comm);
    }

    static std::vector<T> to_unique_segment_id_from_end_elementwise(std::vector<uint8_t> & ends, uint8_t const & non_end_value = 0, MPI_Comm comm = MPI_COMM_WORLD) {
      // convert to 1 and 0
      std::vector<T> results(ends.size());
      std::transform(ends.begin(), ends.end(), results.begin(), [&non_end_value](uint8_t & x) { return (x == non_end_value) ? 0 : 1; });

      // reverse scan
      return scan_op::scan_reverse_elementwise(results, std::plus<T>(), comm);
    }

    static std::vector<T> to_unique_segment_id_from_end_v(std::vector<uint8_t> & ends, uint8_t const & non_end_value = 0, MPI_Comm comm = MPI_COMM_WORLD) {
      // convert to 1 and 0
      std::vector<T> results(ends.size());
      if (ends.size() > 0)
        std::transform(ends.begin(), ends.end(), results.begin(), [&non_end_value](uint8_t & x) { return (x == non_end_value) ? 0 : 1; });
      
      // reverse scan
      return scan_op::scanv_reverse(results, std::plus<T>(), comm);
    }

    static T to_unique_segment_id(T const & non_unique_seg_id, MPI_Comm comm = MPI_COMM_WORLD) {
      auto start = is_start(non_unique_seg_id, comm);
      return to_unique_segment_id_from_start(start, 0, comm);
    }

    static std::vector<T> to_unique_segment_id_elementwise(std::vector<T> & non_unique_seg_ids, MPI_Comm comm = MPI_COMM_WORLD) {
      auto start = is_start_elementwise(non_unique_seg_ids, comm);
      return to_unique_segment_id_from_start_elementwise(start, 0, comm);
    }

    static std::vector<T> to_unique_segment_id_v(std::vector<T> & non_unique_seg_ids, MPI_Comm comm = MPI_COMM_WORLD) {
      auto start = is_start_v(non_unique_seg_ids, comm);
      return to_unique_segment_id_from_start_v(start, 0, comm);
    }
    

  };



  /// segmented scan operations.  requires that the segments have unique ids and all elements in the same segment id
    struct seg_scan {

      /**
       * @brief segmented scan helper function - converts the user specified operation into a segmented scan operation
       */
      template <typename T, typename SegT, typename Func>
      struct seg_scan_op {
          Func f;
          const uint8_t non_terminus_val;

          seg_scan_op(Func _f, uint8_t middle_val = 0) : f(_f), non_terminus_val(middle_val) {}

          /// performs the reduction operator, using segment value as guide.  Note that there are no "segment start marker" version of this since we'd need a corresponding "segment end marker" version, whcih complicates things.
          std::pair<T, SegT> operator()(std::pair<T, SegT> const & x, std::pair<T, SegT> const & y) {
            return (x.second == y.second) ?
              std::pair<T, SegT>(this->f(x.first, y.first), y.second) : y;
          }
      };



      /// MPI Scan  with uniquely marked segments.  we can then split the comm and do normal scan
      template <typename T, typename SegT, typename Func = std::plus<T>, bool exclusive = false >
      static T scan(T & x, SegT seg, Func f, MPI_Comm comm = MPI_COMM_WORLD) {
        int rank;
        MPI_Comm_rank(comm, &rank);

        MPI_Comm seg_comm;
        MPI_Comm_split(comm, seg, rank, &seg_comm);

        // calling normal scan function.
        T result = scan_op::scan<T, Func, exclusive>(x, f, seg_comm);

        MPI_Comm_free(&seg_comm);
        return result;
      }

      template <typename T, typename SegT, typename Func = std::plus<T>, bool exclusive = false  >
      static T scan_reverse(T & x, SegT seg, Func f, MPI_Comm comm = MPI_COMM_WORLD) {
        // reverse the comm, then split it using the same coloring, results in individual reversed comms.
        int rank, p;
        MPI_Comm_rank(comm, &rank);
        MPI_Comm_size(comm, &p);

        MPI_Comm rev_seg_comm;
        MPI_Comm_split(comm, seg, p - rank, &rev_seg_comm);

        // now call local scan
        T result = scan_op::scan<T, Func, exclusive>(x, f, rev_seg_comm);

        MPI_Comm_free(&rev_seg_comm);
        return result;
      }

      /// exclusive scan
      template <typename T, typename SegT, typename Func = std::plus<T> >
      static T exscan(T & x, SegT seg, Func f, MPI_Comm comm = MPI_COMM_WORLD) {
        return seg_scan::scan<T, SegT, Func, true>(x, seg, f, comm);
      }

      template <typename T, typename SegT, typename Func = std::plus<T> >
      static T exscan_reverse(T & x, SegT seg, Func f, MPI_Comm comm = MPI_COMM_WORLD) {
        return seg_scan::scan_reverse<T, SegT, Func, true>(x, seg, f, comm);
      }



      /// scan of a vector, element-wise across all processors.
      template <typename T, typename SegT, typename Func = std::plus<T>, bool exclusive = false >
      static std::vector<T> scan_elementwise(std::vector<T> & x, std::vector<SegT> &seg, Func f, MPI_Comm comm = MPI_COMM_WORLD) {
        // can't do this via comm split because segments may not line up at different element positions.

        if (seg.size() != x.size())
          throw std::logic_error("ERROR: segmented scan input and segment label are not the same size");

        // minimum size:
        size_t s = x.size();
        s = ::mxx::allreduce(s, [](size_t &x , size_t & y){ return ::std::min(x, y); }, comm);
        if (s < x.size())
          WARNINGF("WARNING: %lu elements out of %lu used in scan.", s, x.size());

        // converted data structure
        std::vector<std::pair<T, SegT> > marked(s);
        std::transform(x.begin(), x.begin() + s, seg.begin(), marked.begin(), [](T & x, SegT & y) { return std::pair<T, SegT>(x, y); } );

        // converted op.
        seg_scan_op<T, SegT, Func> func(f);

        std::vector<std::pair<T, SegT> > out = scan_op::scan_elementwise<std::pair<T, SegT>, seg_scan_op<T, SegT, Func>, exclusive>(marked, func, comm);

        // convert data back to original form
        std::vector<T> results(out.size());
        std::transform(out.begin(), out.end(), results.begin(), [](std::pair<T, SegT> & x) { return x.first; } );

        return results;

      }

      template <typename T, typename SegT, typename Func = std::plus<T>, bool exclusive = false >
      static std::vector<T> scan_reverse_elementwise(std::vector<T> & x, std::vector<SegT> &seg, Func f, MPI_Comm comm = MPI_COMM_WORLD) {
        // since we cannot do this via comm splitting (segments at different position in vector may not line up
        // we can just call the scan_elementwise function with a reverse comm.

        MPI_Comm rev_comm;
        mxx::rev_comm(comm, rev_comm);

        // now call local scan
        auto results = seg_scan::scan_elementwise<T, SegT, Func, exclusive>(x, seg, f, rev_comm);

        MPI_Comm_free(&rev_comm);
        return results;

      }


      /// exclusive segmented scan.  first element of a segment is not valid.
      template <typename T, typename SegT, typename Func = std::plus<T> >
      static std::vector<T> exscan_elementwise(std::vector<T> & x, std::vector<SegT> &seg, Func f, MPI_Comm comm = MPI_COMM_WORLD) {
        return seg_scan::scan_elementwise<T, SegT, Func, true>(x, seg, f, comm);
      }

      /// exclusive segmented scan in reverse.  last element of a segment is not valid.
      template <typename T, typename SegT, typename Func = std::plus<T> >
      static std::vector<T> exscan_reverse_elementwise(std::vector<T> & x, std::vector<SegT> &seg, Func f, MPI_Comm comm = MPI_COMM_WORLD) {
        return seg_scan::scan_reverse_elementwise<T, SegT, Func, true>(x, seg, f, comm);
      }


      /**
       * @brief Scan of partitioned vector.
       * @note  cannot split communicator.
       * @details if the segment ids are not unique for each segment, the last (or first) element in the vector
       *   from 2 adjacent processes may have the same segment value when they are not in same segment due to intervening vector elements.
       *   simplest solution is to convert to Unique segment ids first.
       */
      template <typename T, typename SegT, typename Func = std::plus<T>, bool exclusive = false  >
      static std::vector<T> scanv(std::vector<T> & x, std::vector<SegT> &seg, Func f, MPI_Comm comm = MPI_COMM_WORLD) {
        // convert make sure seg is converted to Unique first.

        // can't do this via comm split because segments may not line up at different element positions.

        if (seg.size() != x.size())
          throw std::logic_error("ERROR: segmented scan input and segment label are not the same size");

        // converted data structure
        std::vector<std::pair<T, SegT> > marked(x.size());
        std::transform(x.begin(), x.end(), seg.begin(), marked.begin(), [](T & x, SegT & y) { return std::pair<T, SegT>(x, y); } );

        int rank, p;
        MPI_Comm_rank(comm, &rank);
        MPI_Comm_size(comm, &p);

        // converted op.
        seg_scan_op<T, SegT, Func> func(f);

        std::vector<std::pair<T, SegT> > out;
        if (exclusive)
          scan_op::exscanv(marked, func, comm).swap(out);
        else
          scan_op::scanv(marked, func, comm).swap(out);

        // convert data back to original form
        std::vector<T> results(out.size());
        std::transform(out.begin(), out.end(), results.begin(), [](std::pair<T, SegT> & x) { return x.first; } );

        return results;

      }



      template <typename T, typename SegT, typename Func = std::plus<T>, bool exclusive = false  >
      static std::vector<T> scanv_reverse(std::vector<T> & x, std::vector<SegT> &seg, Func f, MPI_Comm comm = MPI_COMM_WORLD) {
        // convert make sure seg is converted to Unique first.

        // can't do this via comm split because segments may not line up at different element positions.

        if (seg.size() != x.size())
          throw std::logic_error("ERROR: segmented scan input and segment label are not the same size");

        // converted data structure
        std::vector<std::pair<T, SegT> > marked(x.size());
        std::transform(x.begin(), x.end(), seg.begin(), marked.begin(), [](T & x, SegT & y) { return std::pair<T, SegT>(x, y); } );

        // converted op.
        seg_scan_op<T, SegT, Func> func(f);


        std::vector<std::pair<T, SegT> > out;
        if (exclusive)
          scan_op::exscanv_reverse(marked, func, comm).swap(out);
        else
          scan_op::scanv_reverse(marked, func, comm).swap(out);

        // convert data back to original form
        std::vector<T> results(out.size());
        std::transform(out.begin(), out.end(), results.begin(), [](std::pair<T, SegT> & x) { return x.first; } );

        return results;

      }


      // NOTE for segmented exclusive scan, the "first" element of a segment is NOT valid.
      template <typename T, typename SegT, typename Func = std::plus<T> >
      static std::vector<T> exscanv(std::vector<T> & x, std::vector<SegT> &seg, Func f, MPI_Comm comm = MPI_COMM_WORLD) {
        return seg_scan::scanv<T, SegT, Func, true>(x, seg, f, comm);
      }

      template <typename T, typename SegT, typename Func = std::plus<T> >
      static std::vector<T> exscanv_reverse(std::vector<T> & x, std::vector<SegT> &seg, Func f, MPI_Comm comm = MPI_COMM_WORLD) {
        return seg_scan::scanv_reverse<T, SegT, Func, true>(x, seg, f, comm);

      }

  };



  /// perform segmented reduction.  results has same distribution as original data.  segment ids should be unique, with each segment marked with same id.
  struct seg_reduce {

        /**
         * @brief segmented scan helper function - converts the user specified operation into a segmented scan operation
         */
        template <typename T, typename SegT>
        struct replace_op {
            const uint8_t non_terminus_val;

            replace_op(uint8_t middle_val = 0) : non_terminus_val(middle_val) {}

            std::pair<T, SegT> operator()(std::pair<T, SegT> const & x, std::pair<T, SegT> const & y) {
              return (x.second == y.second) ?
                x : y;
            }
        };



        template <typename T, typename SegT, typename Func = std::plus<T> >
        static T reduce(T & x, SegT & seg, Func f, MPI_Comm comm = MPI_COMM_WORLD) {

          int rank;
          MPI_Comm_rank(comm, &rank);

          // split the comm
          MPI_Comm seg_comm;
          MPI_Comm_split(comm, seg, rank, &seg_comm);


          // then do allreduce for each comm

          // get user op note that this is NOT communitative
          MPI_Op op = mxx::create_user_op<T, Func >(f);
          // get type
          mxx::datatype<T > dt;
          T result;
          // perform reduction
          MPI_Allreduce(&x, &result, 1, dt.type(), op, seg_comm);

          // clean up op
          mxx::free_user_op<T >(op);
          MPI_Comm_free(&seg_comm);

          return result;
        }


        template <typename T, typename SegT, typename Func = std::plus<T> >
        static std::vector<T> reducen(std::vector<T> & x, std::vector<SegT> & seg, Func f, MPI_Comm comm = MPI_COMM_WORLD) {

          if (seg.size() != x.size())
            throw std::logic_error("ERROR: segmented scan input and segment label are not the same size");

          // minimum size:
          size_t s = x.size();
          s = ::mxx::allreduce(s, [](size_t &x , size_t & y){ return ::std::min(x, y); }, comm);
          if (s < x.size())
            WARNINGF("WARNING: %lu elements out of %lu used in scan.", s, x.size());

          // converted data structure
          std::vector<std::pair<T, SegT> > marked(s);
          std::transform(x.begin(), x.begin() + s, seg.begin(), marked.begin(), [](T & x, SegT & y) { return std::pair<T, SegT>(x, y); } );


          // first do a forward segmented scan with the user function.
          seg_scan::seg_scan_op<T, SegT, Func> func(f);
          std::vector<std::pair<T, SegT> > out = scan_op::scan_elementwise(marked, func, comm);


          // then do a reverse segmented scan with a "replace" function
          MPI_Comm rev_comm;
          mxx::rev_comm(comm, rev_comm);
          replace_op<T, SegT> func2;
          out = scan_op::scan_reverse_elementwise(out, func2, comm);
          MPI_Comm_free(&rev_comm);

          // convert data back to original form
          std::vector<T> results(out.size());
          std::transform(out.begin(), out.end(), results.begin(), [](std::pair<T, SegT> & x) { return x.first; } );

          return results;
        }


        /**
         * @brief  perform allreduce for each segment.
         * @details if the segment ids are not unique for each segment, the last (or first) element in the vector
         *   from 2 adjacent processes may have the same segment value when they are not in same segment due to intervening vector elements.
         *   simplest solution is to convert to Unique segment ids first.
         *
         * @param x
         * @param seg
         * @param f
         * @param comm
         * @return
         */
        template <typename T, typename SegT, typename Func = std::plus<T> >
        static std::vector<T> reducev(std::vector<T> & x, std::vector<SegT> & seg, Func f, MPI_Comm comm = MPI_COMM_WORLD) {
          if (seg.size() != x.size())
            throw std::logic_error("ERROR: segmented scan input and segment label are not the same size");


          // converted data structure
          std::vector<std::pair<T, SegT> > marked(x.size());
          std::transform(x.begin(), x.end(), seg.begin(), marked.begin(), [](T & x, SegT & y) { return std::pair<T, SegT>(x, y); } );


          // first do a forward segmented scan with the user function.
          seg_scan::seg_scan_op<T, SegT, Func> func(f);
          std::vector<std::pair<T, SegT> > out = scan_op::scanv(marked, func, comm);


          // then do a reverse segmented scan with a "replace" function
          MPI_Comm rev_comm;
          mxx::rev_comm(comm, rev_comm);
          replace_op<T, SegT> func2;
          out = scan_op::scanv_reverse(out, func2, comm);
          MPI_Comm_free(&rev_comm);

          // convert data back to original form
          std::vector<T> results(out.size());
          std::transform(out.begin(), out.end(), results.begin(), [](std::pair<T, SegT> & x) { return x.first; } );

          return results;

        }

  };
}

#endif /* SRC_IO_MXX_SUPPORT_HPP_ */
