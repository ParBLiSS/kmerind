/**
 * @file    mxx_comm_filter.hpp
 * @ingroup
 * @author  tpan
 * @brief
 * @details
 *
 * Copyright (c) 2016 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef SRC_UTILS_MXX_FAST_COMM_HPP_
#define SRC_UTILS_MXX_FAST_COMM_HPP_

#include "mxx/env.hpp"
#include "mxx/comm.hpp"
#include "mxx/benchmark.hpp"

namespace bliss {

  namespace mxx {

   bool get_fast_nodes(::mxx::comm const & comm, ::mxx::env const & e, int target_node_count) {
     // create shared-mem MPI+MPI hybrid communicator
     ::mxx::hybrid_comm hc(comm);

     // assert same number processors per node
     int proc_per_node = hc.local.size();
     if (!::mxx::all_same(proc_per_node, comm)) {
         throw std::logic_error("Error: this benchmark assumes the same number of processors per node");
     }

     // assert we have an even number of nodes
     int num_nodes = hc.num_nodes();
     if (num_nodes % 2 != 0) {
         throw std::logic_error("Error: this benchmark assumes an even number of nodes");
     }

     // number to vote off.
     int n_vote_off = num_nodes - target_node_count;

     // split by pairwise bandwidth
     std::vector<double> bw_row = ::mxx::pairwise_bw_matrix(hc);
     ::mxx::print_bw_matrix_stats(hc, bw_row);
     bool part = ::mxx::vote_off(hc, n_vote_off, bw_row);

     // do some reporting.
     if (hc.global.rank() == 0)
         std::cout << "Before vote off: " << std::endl;
     ::mxx::bw_all2all(hc.global, hc.local);

     if (hc.global.rank() == 0)
         std::cout << "After vote off: " << std::endl;
     hc.with_nodes(part, [&](const ::mxx::hybrid_comm& subhc) {
             ::mxx::bw_all2all(subhc.global, subhc.local);
     });

     // create a subcommunicator and send it back
     return part;

   }

  }

}


#endif /* SRC_UTILS_MXX_FAST_COMM_HPP_ */
