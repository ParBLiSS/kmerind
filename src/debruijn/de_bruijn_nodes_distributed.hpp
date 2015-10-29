/*
 * de_bruijn_nodes_distributed.hpp
 *
 *  Created on: Aug 6, 2015
 *      Author: yongchao
 */

#ifndef DE_BRUIJN_NODES_DISTRIBUTED_HPP_
#define DE_BRUIJN_NODES_DISTRIBUTED_HPP_

#include <unordered_map>  // local storage hash table  // for multimap
#include <unordered_set>  // local storage hash table  // for multimap
#include <utility> 			  // for std::pair

//#include <sparsehash/dense_hash_map>  // not a multimap, where we need it most.
#include <functional> 		// for std::function and std::hash
#include <algorithm> 		// for sort, stable_sort, unique, is_sorted
#include <iterator>  // advance, distance

#include <cstdint>  // for uint8, etc.

#include <type_traits>
#include "debruijn/de_bruijn_node_trait.hpp"	//node trait data structure storing the linkage information to the node
#include "containers/distributed_map_base.hpp"
#include "containers/distributed_unordered_map.hpp"
#include "containers/unordered_vecmap.hpp"

#include "utils/timer.hpp"  // for timing.
#include "utils/logging.h"


namespace bliss{
	namespace de_bruijn{

		template<typename Key, typename T,
		  class Comm,
		  template <typename> class KeyTransform,
		  template <typename, bool> class Hash,
		  class Equal = ::std::equal_to<Key>,
		  class Alloc = ::std::allocator< ::std::pair<const Key, T> >
		  >
		  class de_bruijn_nodes_distributed : public ::dsc::unordered_map<Key, T, Comm, KeyTransform, Hash, Equal, Alloc> {
			  using Base = ::dsc::unordered_map<Key, T, Comm, KeyTransform, Hash, Equal, Alloc>;

			public:
			  using local_container_type = typename Base::local_container_type;

			  // std::unordered_multimap public members.
			  using key_type              = typename local_container_type::key_type;
			  using mapped_type           = typename local_container_type::mapped_type;
			  using value_type            = typename local_container_type::value_type;
			  using hasher                = typename local_container_type::hasher;
			  using key_equal             = typename local_container_type::key_equal;
			  using allocator_type        = typename local_container_type::allocator_type;
			  using reference             = typename local_container_type::reference;
			  using const_reference       = typename local_container_type::const_reference;
			  using pointer               = typename local_container_type::pointer;
			  using const_pointer         = typename local_container_type::const_pointer;
			  using iterator              = typename local_container_type::iterator;
			  using const_iterator        = typename local_container_type::const_iterator;
			  using size_type             = typename local_container_type::size_type;
			  using difference_type       = typename local_container_type::difference_type;

			protected:
			  // defined Communicator as a friend
			  friend Comm;

			  /**
			   * @brief insert new elements in the distributed unordered_multimap.
			   * @param first
			   * @param last
			   */
			  template <class InputIterator>
			  size_t local_insert(InputIterator first, InputIterator last) {
				  int32_t relative_strand;
				  size_t before = this->c.size();

				  /*reserve space*/
				  this->local_reserve(before + ::std::distance(first, last));

				  // allocate a node
          auto node = this->c.find(first->first);

          using EdgeType = decltype(node->second);
          using Alphabet = typename EdgeType::Alphabet;

				  /*iterate each input tuple*/
				  for (auto it = first; it != last; ++it) {

				    node = this->c.find(it->first);

				    /*tranform from <key, int> to <node, node_info>*/
					  if(node == this->c.end()){
						  /*create a new node*/
						  auto ret = this->c.emplace(::std::make_pair(it->first, T()));
						  if(ret.second == false){
							  cerr << "Insertion failed at line " << __LINE__ << " in file " << __FILE__ << endl;
							  exit(-1);
						  }

#if 0
              ::std::cerr << "Not exist in the hash table" << ::std::endl;
              ::std::cerr << bliss::utils::KmerUtils::toASCIIString(ret.first->first)  << ::std::endl;
              ::std::cerr << "0x" << std::hex << it->second << ::std::endl;
              ::std::cerr << bliss::utils::KmerUtils::toASCIIString(it->first) << ::std::endl << ::std::endl;
#endif
						  node = ret.first;  // now use it.

	            // determine if the node in the graph has the same orientation as the one being inserted.
						  relative_strand = bliss::de_bruijn::node::SENSE;
					  } else {
#if 0
             /*update the node*/
              ::std::cerr << "Exist in the hash table" << ::std::endl;
              ::std::cerr << bliss::utils::KmerUtils::toASCIIString(node->first)  << ::std::endl;
              ::std::cerr << "0x" << std::hex << it->second << ::std::endl;
              ::std::cerr << bliss::utils::KmerUtils::toASCIIString(it->first) << ::std::endl << ::std::endl;
#endif

              // determine if the node in the graph has the same orientation as the one being inserted.
              relative_strand = (node->first == it->first) ? bliss::de_bruijn::node::SENSE : bliss::de_bruijn::node::ANTI_SENSE;

					  }


					  // if different, swap and reverse complement the edges.  else, use as is.
					  if (relative_strand == bliss::de_bruijn::node::ANTI_SENSE) {
              node->second.update(bliss::de_bruijn::node::input_edge_utils::reverse_complement_edges<Alphabet>(it->second));

					  } else {
              node->second.update(it->second);
					  }

				  }
				  return this->c.size() - before;
			  }

			  /**
			   * @brief insert new elements in the distributed unordered_multimap.
			   * @param first
			   * @param last
			   */
			  template <class InputIterator, class Predicate>
			  size_t local_insert(InputIterator first, InputIterator last, Predicate const & pred) {
				  int32_t relative_strand;
				  size_t before = this->c.size();

				  this->local_reserve(before + ::std::distance(first, last));
				  auto node = this->c.find(first->first);

          using EdgeType = decltype(node->second);
          using Alphabet = typename EdgeType::Alphabet;

				  for (auto it = first; it != last; ++it) {
					if (pred(*it)) {

					  node = this->c.find(it->first);

					  /*tranform from <key, int> to <node, node_info*/
					  if(node == this->c.end()){
						  /*create a new node*/
						  auto ret = this->c.emplace(::std::make_pair(it->first, T()));
						  if(ret.second == false){
							  cerr << "Insertion failed at line " << __LINE__ << " in file " << __FILE__ << endl;
							  exit(-1);
						  }

						  node = ret.first;  // now use it.

						  /*update the node*/
						  relative_strand = bliss::de_bruijn::node::SENSE;
					  }else{
						 /*update the node*/
						  relative_strand = (node->first == it->first) ? bliss::de_bruijn::node::SENSE : bliss::de_bruijn::node::ANTI_SENSE;
					  }

            // if different, swap and reverse complement the edges.  else, use as is.
            if (relative_strand == bliss::de_bruijn::node::ANTI_SENSE) {
              node->second.update(bliss::de_bruijn::node::input_edge_utils::reverse_complement_edges<Alphabet>(it->second));

            } else {
              node->second.update(it->second);
            }


					}
				  }
				  return this->c.size() - before;

			  }

			public:
			  de_bruijn_nodes_distributed(const mxx::comm& _comm) : Base(_comm) {/*do nothing*/}

			  virtual ~de_bruijn_nodes_distributed() {/*do nothing*/};

			  /*transform function*/

			  /**
				* @brief insert new elements in the distributed unordered_multimap.
				* @param first
				* @param last
				*/
			   template <typename InputEdgeType, typename Predicate = ::dsc::Identity>
			   size_t insert(std::vector<::std::pair<Key, InputEdgeType> >& input, bool sorted_input = false, Predicate const & pred = Predicate()) {
				 // even if count is 0, still need to participate in mpi calls.  if (input.size() == 0) return;
				 TIMER_INIT(insert);

				 TIMER_START(insert);
				 TIMER_END(insert, "start", input.size());


				 // communication part
				 if (this->comm_size > 1) {
				   TIMER_START(insert);
				   // get mapping to proc
				   ::std::vector<size_t> send_counts = mxx::bucketing(input, this->key_to_rank, this->comm_size);
				   TIMER_END(insert, "bucket", input.size());

		 //          TIMER_START(insert);
		 //          // keep unique only.  may not be needed - comm speed may be faster than we can compute unique.
		 //          mxx2::retain_unique<local_container_type, typename Base::TransformedEqual>(input, send_counts, sorted_input);
		 //          TIMER_END(insert, "uniq1", input.size());

				   TIMER_COLLECTIVE_START(insert, "a2a", this->comm);
				   input = mxx::all2allv(input, send_counts, this->comm);
				   TIMER_END(insert, "a2a", input.size());
				 }

				 TIMER_START(insert);
				 // local compute part.  called by the communicator.
				 size_t count = 0;
				 if (!::std::is_same<Predicate, ::dsc::Identity>::value)
				   count = this->local_insert(input.begin(), input.end(), pred);
				 else
				   count = this->local_insert(input.begin(), input.end());
				 TIMER_END(insert, "insert", this->c.size());

				 TIMER_REPORT_MPI(insert, this->comm.rank(), this->comm);

				 return count;
			   }
		};
	}/*de_bruijn*/
}/*bliss*/

#endif /* DE_BRUIJN_NODES_DISTRIBUTED_HPP_ */
