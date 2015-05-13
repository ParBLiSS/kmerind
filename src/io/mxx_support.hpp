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
    public datatype_contiguous<decltype(bliss::io::FASTQ::SequenceId::file_pos),1> {};

}





namespace mxx2 {
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
      int pid = 0;
      for (int i = 0; i < msgs.size(); ++i)
      {
          pids[i] = target_p_fun(msgs[i]);
          send_counts[pids[i]]++;
      }

      // get all2all params
      std::vector<int> recv_counts = mxx::all2all(send_counts, 1, comm);
      std::vector<int> send_displs = mxx::get_displacements(send_counts);
      std::vector<int> recv_displs = mxx::get_displacements(recv_counts);

      // copy.  need to be able to track current position within each block.
      std::vector<int> offset = send_displs;
      std::vector<T> send_buffer;
      if (msgs.size() > 0)
          send_buffer.resize(msgs.size());
      for (int i = 0; i < msgs.size(); ++i)
      {
          send_buffer[offset[pids[i]]++] = msgs[i];
      }


      // resize messages to fit recv
      std::size_t recv_size = recv_displs[p-1] + recv_counts[p-1];
      msgs.clear();
      //msgs.shrink_to_fit();
      msgs.resize(recv_size);
      //msgs = std::vector<T>(recv_size);

      // get MPI type
      mxx::datatype<T> dt;
      MPI_Datatype mpi_dt = dt.type();

      // all2all
      MPI_Alltoallv(&send_buffer[0], &send_counts[0], &send_displs[0], mpi_dt,
                    &msgs[0], &recv_counts[0], &recv_displs[0], mpi_dt, comm);
      // done, result is returned in vector of input messages

      return recv_counts;
  }


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


}

#endif /* SRC_IO_MXX_SUPPORT_HPP_ */
