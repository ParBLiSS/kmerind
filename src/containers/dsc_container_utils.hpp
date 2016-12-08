/*
 * Copyright 2015 Georgia Institute of Technology
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/**
 * @file    fsc_container_utils.hpp
 * @ingroup
 * @author  tpan
 * @brief
 * @details
 *
 */
#ifndef SRC_CONTAINERS_DSC_CONTAINER_UTILS_HPP_
#define SRC_CONTAINERS_DSC_CONTAINER_UTILS_HPP_

#include <iterator>  // iterator_traits
#include <unordered_set>
#include <random>
#include <algorithm>  // upper bound, unique, sort, etc.

#include "containers/fsc_container_utils.hpp"

#include "utils/benchmark_utils.hpp"


namespace dsc {

	template <typename CONTAINER>
	bool empty(CONTAINER const & c, mxx::comm const & comm) {
		bool local_empty = (c.size() == 0);
	  if (comm.size() == 1) {
			if (local_empty) printf("rank 0/1 input is EMPTY.\n");
		return local_empty;
	} else { // all reduce
		local_empty =  mxx::all_of(local_empty, comm);
		if (local_empty && comm.rank() == 0) printf("rank ALL/%d inputs are all EMPTY.\n", comm.size());
		return local_empty;
	  }
	}


  // =============== convenience functions for distribution of vector via all to all and a rank mapping function
	// TODO: make this cleaner...


  /**
   * @brief       distribute.  speed version.  no guarantee of output ordering, but actually is same.
   * @param[IN/OUT] vals    vals to distribute.  sortedness is NOT kept because of inplace bucketing.
   * @param[IN/OUT] sorted_input    indicates if input is sorted.  and whether each bucket is sorted.
   *
   * @return received counts.
   */
  template <typename V, typename ToRank>
  ::std::vector<size_t> distribute(::std::vector<V>& vals, ToRank const & to_rank, bool & sorted_input,
                                   mxx::comm const &_comm) {
    BL_BENCH_INIT(distribute);

      // speed over mem use.  mxx all2allv already has to double memory usage. same as stable distribute.

      BL_BENCH_START(distribute);
      // distribute (communication part)
      std::vector<size_t> send_counts = mxx::bucketing(vals, to_rank, _comm.size());
      BL_BENCH_END(distribute, "bucket", vals.size());

      // distribute (communication part)
      BL_BENCH_COLLECTIVE_START(distribute, "a2a", _comm);
      mxx::all2allv(vals, send_counts, _comm).swap(vals);
      BL_BENCH_END(distribute, "a2a", vals.size());

      BL_BENCH_START(distribute);
      std::vector<size_t> recv_counts= mxx::all2all(send_counts, _comm);
      BL_BENCH_END(distribute, "a2a_counts", vals.size());

      BL_BENCH_REPORT_MPI_NAMED(distribute, "map_base:distribute", _comm);

      return recv_counts;
  }


    /**
     * @brief       distribute and ensure uniqueness within each bucket.  speed version.  no guarantee of output ordering.
     *				to_rank uses distribution hash and transform.  unique should use the storage hash and transform
     * @param[IN/OUT] keys    keys to distribute. sortedness is NOT kept because of inplace bucketing.,
     * @param[IN/OUT] sorted_input    indicates if input is sorted.  and whether each bucket is sorted.
     * @return received counts
     */
    template <typename V, typename ToRank, typename Hash, typename Eq>
    ::std::vector<size_t> distribute_unique(::std::vector<V>& vals, ToRank const & to_rank, bool & sorted_input,
                                            mxx::comm const &_comm, const Hash & hash = Hash(), const Eq & equal = Eq()) {

      BL_BENCH_INIT(distribute);
        // go for speed.   mxx::bucketing uses extra copy.  okay, since all2all also does.

        BL_BENCH_START(distribute);
        // distribute (communication part)
        std::vector<size_t> send_counts = mxx::bucketing(vals, to_rank, _comm.size());
        BL_BENCH_END(distribute, "bucket", vals.size());

        // using set is okay.
        BL_BENCH_START(distribute);
        // distribute (communication part)
        ::fsc::bucket_unique(vals, send_counts, sorted_input, hash, equal);
        BL_BENCH_END(distribute, "unique", vals.size());

        // distribute (communication part)
        BL_BENCH_COLLECTIVE_START(distribute, "a2a", _comm);
        mxx::all2allv(vals, send_counts, _comm).swap(vals);
        BL_BENCH_END(distribute, "a2a", vals.size());

        BL_BENCH_START(distribute);
        std::vector<size_t> recv_counts= mxx::all2all(send_counts, _comm);
        BL_BENCH_END(distribute, "a2a_counts", vals.size());

        BL_BENCH_REPORT_MPI_NAMED(distribute, "map_base:distribute_unique", _comm);


        return recv_counts;
    }



    /**
     * @brief       distribute.  order preserving version.  guarantees output ordering.
     *				to_rank uses distribution hash and transform.  unique should use the storage hash and transform
     * @param[IN/OUT] vals    vals to distribute.  sortedness is KEPT within each bucket
     * @param[IN/OUT] sorted_input    indicates if input is sorted.  and whether each bucket is sorted.
     *
     * @return received counts.
     */
    template <typename V, typename ToRank, typename Less>
    ::std::vector<size_t> distribute_sorted(::std::vector<V>& vals, ToRank const & to_rank, bool & sorted_input,
                                            mxx::comm const &_comm, const Less & less = Less()) {
      BL_BENCH_INIT(distribute);

      // ordering over speed/mem use.

        BL_BENCH_START(distribute);
        // distribute (communication part)
        std::vector<size_t> send_counts = mxx::bucketing(vals, to_rank, _comm.size());
        BL_BENCH_END(distribute, "bucket", vals.size());

        // using set is okay.
        BL_BENCH_START(distribute);
        // distribute (communication part)
        ::fsc::bucket_sort(vals, send_counts, sorted_input, less);
        BL_BENCH_END(distribute, "sort", vals.size());

        // distribute (communication part)
        BL_BENCH_COLLECTIVE_START(distribute, "a2a", _comm);
        mxx::all2allv(vals, send_counts, _comm).swap(vals);
        BL_BENCH_END(distribute, "a2a", vals.size());

        BL_BENCH_START(distribute);
        std::vector<size_t> recv_counts= mxx::all2all(send_counts, _comm);
        BL_BENCH_END(distribute, "a2a_counts", vals.size());

        BL_BENCH_REPORT_MPI_NAMED(distribute, "map_base:distribute_sorted", _comm);


        return recv_counts;
    }

    /**
     * @brief       distribute.  order preserving version.  guarantees output ordering.
     *				to_rank uses distribution hash and transform.  unique should use the storage hash and transform
     * @param[IN/OUT] keys    keys to distribute. sortedness is KEPT.,
     * @param[IN/OUT] sorted_input    indicates if input is sorted.  and whether each bucket is sorted.
     * @return received counts
     */
    template <typename V, typename ToRank, typename Less, typename Eq>
    ::std::vector<size_t> distribute_sorted_unique(::std::vector<V>& vals, ToRank const & to_rank, bool & sorted_input,
                                                   mxx::comm const &_comm, const Less & less = Less(), const Eq & equal = Eq()) {
      BL_BENCH_INIT(distribute);
      // ordering over speed/mem use.

        BL_BENCH_START(distribute);
        // distribute (communication part)
        std::vector<size_t> send_counts = mxx::bucketing(vals, to_rank, _comm.size());
        BL_BENCH_END(distribute, "bucket", vals.size());

        BL_BENCH_START(distribute);
        // distribute (communication part)
        ::fsc::bucket_sorted_unique(vals, send_counts, sorted_input, less, equal);
        BL_BENCH_END(distribute, "sort_unique", vals.size());


        // distribute (communication part)
        BL_BENCH_COLLECTIVE_START(distribute, "a2a", _comm);
        mxx::all2allv(vals, send_counts, _comm).swap(vals);
        BL_BENCH_END(distribute, "a2a", vals.size());

        BL_BENCH_START(distribute);
        std::vector<size_t> recv_counts= mxx::all2all(send_counts, _comm);
        BL_BENCH_END(distribute, "a2a_counts", vals.size());

        BL_BENCH_REPORT_MPI_NAMED(distribute, "map_base:distribute_sorted_unique", _comm);


        return recv_counts;
	}


    /**
     * @brief       distribute and ensure uniqueness within each bucket.  speed version.  no guarantee of output ordering.
     *				to_rank uses distribution hash and transform.  unique should use the storage hash and transform
     * @param[IN/OUT] keys    keys to distribute. sortedness is NOT kept because of inplace bucketing.,
     * @param[IN/OUT] sorted_input    indicates if input is sorted.  and whether each bucket is sorted.
     * @return received counts
     */
    template <typename V, typename ToRank, typename Reduc>
    ::std::vector<size_t> distribute_reduce(::std::vector<V>& vals, ToRank const & to_rank, bool & sorted_input,
                                            mxx::comm const &_comm, const Reduc & reducer = Reduc()) {

      BL_BENCH_INIT(distribute);
        // go for speed.   mxx::bucketing uses extra copy.  okay, since all2all also does.

  	if (_comm.rank() == 0) printf("start\n");  fflush(stdout);


        BL_BENCH_START(distribute);
        // distribute (communication part)
        std::vector<size_t> send_counts = mxx::bucketing_inplace(vals, to_rank, _comm.size());
        BL_BENCH_END(distribute, "bucket", vals.size());

    	if (_comm.rank() == 0) printf("bucket\n");  fflush(stdout);

    	double mean, stdev;
    	for (auto it = send_counts.begin(), max = send_counts.end(); it != max; ++it) {
    		mean += *it;
    		stdev += (*it) * (*it);
    	}
    	mean /= _comm.size();
    	stdev -= sqrt((stdev / _comm.size()) - (mean * mean));
    	if (_comm.rank() == 0) printf("mean: %f, stdev %f\n", mean, stdev);


        // using set is okay.
        BL_BENCH_START(distribute);
        // distribute (communication part)
        ::fsc::bucket_reduce(vals, send_counts, sorted_input, reducer);
        BL_BENCH_END(distribute, "reduce", vals.size());

    	if (_comm.rank() == 0) printf("reduce\n");  fflush(stdout);


        // distribute (communication part)
        BL_BENCH_COLLECTIVE_START(distribute, "a2a", _comm);
        mxx::all2allv(vals, send_counts, _comm).swap(vals);
        BL_BENCH_END(distribute, "a2a", vals.size());


    	if (_comm.rank() == 0) printf("a2a\n");  fflush(stdout);

        BL_BENCH_START(distribute);
        std::vector<size_t> recv_counts= mxx::all2all(send_counts, _comm);
        BL_BENCH_END(distribute, "a2a_counts", vals.size());

    	if (_comm.rank() == 0) printf("a2acounts\n");  fflush(stdout);

        BL_BENCH_REPORT_MPI_NAMED(distribute, "map_base:distribute_reduce", _comm);


        return recv_counts;
    }

    /**
     * @brief samples a vector to be evenly split, sorting the samples to get splitters.
     * @note  goal is to make the sample sorting process parallel.
     * 		NOTE that duplicates are removed.
     * 		NOTE sampling uses random number generator.
     *
     *      random sample p.
	 *		local sort the p samples.
	 *		use all2allv with overlapping range to send 2 to each proc.
	 *			proc 1:  1,2
	 *			proc 2:  2,3
	 *			...
	 *			proc p:  p
	 *		local sort the 2p samples and reduce.
	 *		pick the middle, first p-1 procs.
	 *
	 *		allgatherv to form splitters.
	 *
	 *		sort splitters, and reduce
	 *
	 *		use splitters to split (get send_counts
	 *
	 *		all2allv
	 *
	 *		sort and reduce.
	 *
	 *		block decomp.
	 *		get 1st of each as final splitter.
	 *		allgather
	 *  NOT WELL TESTED.
	 * @param vals   vector (distributed on multiple procs. that needs to have its splitters generated.
     */
    template <typename V, typename Less, typename Eq>
    ::std::vector<V> sample_for_splitters(::std::vector<V>& vals, bool & sorted_input,
                    mxx::comm const &_comm, const Less & less = Less(), const Eq & equal = Eq()) {

    	// random number generator setup.
		std::default_random_engine generator;
		std::uniform_int_distribution<size_t> distribution(0, (vals.size() - 1));
		auto long_rand = std::bind ( distribution, generator );

		// random sample p.
		std::vector< V > samples;
		for (int i = 0; i < _comm.size(); ++i) {
			samples.emplace_back(vals[long_rand()]);
		}

//		// local sort p
//		std::sort(samples.begin(), samples.end(), less);
//
//		// shuffle, 2 items each
//		std::vector<size_t> send_counts(_comm.size(), 2);
//		send_counts[_comm.size() - 1] = 1;
//		std::vector<size_t> send_displs(_comm.size(), 0);
//		std::iota(send_displs.begin(), send_displs.end(), 0);
//
//		std::vector<size_t> recv_counts = ::mxx::all2all(send_counts, _comm);
//		std::vector<size_t> recv_displs = ::mxx::impl::get_displacements(send_counts);
//		size_t recv_total = recv_dspls.back() + recv_counts.back();
//		::std::vector<Key> moved_samples(recv_total);
//		mxx::all2allv(&(samples[0]), send_counts, send_displs,
//				&(moved_samples[0]), recv_counts, recv_displs,
//				, _comm);
		::std::vector<V> moved_samples = mxx::all2all(samples, _comm);

		// local sort again, with reduction this time.
		std::sort(moved_samples.begin(), moved_samples.end(), less);
		auto mend = std::unique(moved_samples.begin(), moved_samples.end(), equal);
		moved_samples.erase(mend, moved_samples.end());

//		// at this point, we expect overlaps, with worst case
//		// being p(p-1) for the split.
		std::vector< V > splitters;
//		if (_comm.rank() < (_comm.size() - 1)) splitters.emplace_back(moved_samples.back());
		if (_comm.rank() > 0) splitters.emplace_back(moved_samples.front());
		::mxx::allgatherv(splitters, _comm).swap(splitters);

		// allgather splitters and sort it
		std::sort(splitters.begin(), splitters.end(), less);
		auto s_end = std::unique(splitters.begin(), splitters.end(), equal);
		splitters.erase(s_end, splitters.end());

		// compute the send counts for samples - no identical splitters, so won't see things crossing boundaries.
		std::vector<size_t> send_counts =
				::mxx::impl::split(samples.begin(), samples.end(), less, splitters, _comm);

		// and all2all the samples (so there are no repeated values between nodes)
		mxx::all2allv(samples, send_counts, _comm).swap(samples);

		// sort and reduce the split samples
		std::sort(samples.begin(), samples.end(), less);
		mend = std::unique(samples.begin(), samples.end(), equal);
		samples.erase(mend, samples.end());

		// block decompose again
		mxx::stable_distribute(samples, _comm).swap(moved_samples);

		splitters.clear();
		if (_comm.rank() > 0) splitters.emplace_back(moved_samples.front());
		mxx::allgatherv(splitters, _comm).swap(splitters);

		return splitters;

    }

}  // namespace dsc


#endif /* SRC_CONTAINERS_DSC_CONTAINER_UTILS_HPP_ */
