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

#include <unistd.h>   // for gethostname

#include "utils/container_traits.hpp"

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





namespace mxx2 {
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
    for (int i = 0; i < map.size(); ++i) {
      e = ::std::lower_bound(b, buffer.end(), map[i].first, comp);
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

  /// reduce operation with minloc
  template <typename T>
  ::std::pair<::std::vector<T>, ::std::vector<int> > reduce_loc(::std::vector<T> & v, MPI_Op op, MPI_Comm comm = MPI_COMM_WORLD, int root = 0)
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


    // struct
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
      printf("NOTE: scatterv using isend version.\n");
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


}

#endif /* SRC_IO_MXX_SUPPORT_HPP_ */
