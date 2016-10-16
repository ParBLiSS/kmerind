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
#include <algorithm>  // upper bound, unique, sort, etc.
#include <random>

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
	 * @brief   compute the element index mapping between input and bucketed output.
	 *
	 * this is an "deconstructed" version of mxx bucketing function
	 *
	 * It returns an array with the bucketed position of the input entries, i.e.
	 *   output to input mapping
	 * 		an array x, such that input[i] = output[x[i]]; and the mapping and input share the same coordinates.
	 * 		then to reconstruct input, we can walk through the mapping and random access output array,
	 * 		which should be faster as the writes are consecutive in cache while the reads are random.
	 * 		also, this is more friendly with emplacing into a vector.
	 * the mapping is useful to avoid recomputing key_func.
	 *	 *
	 * this is faster than mxx by about 2x (prob due to compiting key_func just once.), but uses another vector of size 8N.
	 *
	 * @tparam T            Input type
	 * @tparam Func         Type of the key function.
	 * @param input[in|out] Contains the values to be bucketed.
	 *                      This vector is both input and output.
	 * @param key_func      A function taking a type T and returning the bucket index
	 *                      in the range [0, num_buckets).
	 * @param num_buckets   The total number of buckets.
	 *
	 * @return              The number of elements in each bucket.
	 *                      The size of this vector is `num_buckets`.
	 *                      The source position for each entry in the bucketed input.
	 *                      the size of this vector is size of input.
	 */
	template <typename T, typename Func>
	std::pair<std::vector<size_t>, std::vector<size_t> >
	assign_to_buckets(std::vector<T> const & input, Func const & key_func, size_t num_buckets) {
	    if (num_buckets == 0)
	    	return std::make_pair(std::vector<size_t>(), std::vector<size_t>());

		// initialize number of elements per bucket
	    std::vector<size_t> bucket_counts(num_buckets, 0);

	    // if input is empty, simply return
	    if (input.size() == 0)
	        return std::make_pair(std::vector<size_t>(), bucket_counts);

	    std::vector<size_t> i2o;
	    i2o.reserve(input.size());

	    // [1st pass]: compute bucket counts and input to bucket assignment.
	    size_t p;
	    for (auto it = input.begin(); it != input.end(); ++it) {
	    	p = key_func(*it);

	        assert((0 <= p) && ((size_t)p < num_buckets));

	        i2o.emplace_back(p);
	        ++bucket_counts[p];
	    }

	    // get offsets of where buckets start (= exclusive prefix sum)
	    // iterator once.
	    std::vector<std::size_t> offset(bucket_counts.size(), 0);
	    for (size_t i = 1, j = 0; i < num_buckets; ++i, ++j) {
	    	offset[i] = offset[j] + bucket_counts[j];
	    }

	    // [2nd pass]: saving elements into correct position, and save the final position.
	    for (size_t i = 0; i < input.size(); ++i) {
	    	i2o[i] = offset[i2o[i]]++; 	// output position from bucket id (i2o[i])
	    }

	    return std::make_pair(i2o, bucket_counts);
	}


	/**
	 * @brief   compute the element index mapping between input and bucketed output.
	 *
	 * this is an "deconstructed" version of mxx bucketing function
	 *
	 * It returns an array with the bucketed position of the input entries, i.e.
	 *   output to input mapping
	 * 		an array x, such that input[i] = output[x[i]]; and the mapping and input share the same coordinates.
	 * 		then to reconstruct input, we can walk through the mapping and random access output array,
	 * 		which should be faster as the writes are consecutive in cache while the reads are random.
	 * 		also, this is more friendly with emplacing into a vector.
	 * the mapping is useful to avoid recomputing key_func.
	 *	 *
	 * this is faster than mxx by about 2x (prob due to compiting key_func just once.), but uses another vector of size 8N.
	 *
	 * @tparam T            Input type
	 * @tparam Func         Type of the key function.
	 * @param input[in|out] Contains the values to be bucketed.
	 *                      This vector is both input and output.
	 * @param key_func      A function taking a type T and returning the bucket index
	 *                      in the range [0, num_buckets).
	 * @param num_buckets   The total number of buckets.
	 *
	 * @return              The number of elements in each bucket.
	 *                      The size of this vector is `num_buckets`.
	 *                      The source position for each entry in the bucketed input.
	 *                      the size of this vector is size of input.
	 */
	template <typename T, typename Func>
	std::vector<size_t>
	assign_and_bucket(std::vector<T>& input, Func const & key_func, size_t num_buckets) {
	    if (num_buckets == 0)
	    	return std::vector<size_t>();

		// initialize number of elements per bucket
	    std::vector<size_t> bucket_counts(num_buckets, 0);

	    // if input is empty, simply return
	    if (input.size() == 0)
	        return bucket_counts;

	    if ((num_buckets >> 32) > 0) {  // more than 32 bit of buckets.
			std::vector<size_t> i2o;
			i2o.reserve(input.size());

			// [1st pass]: compute bucket counts and input to bucket assignment.
			size_t p;
			for (auto it = input.begin(); it != input.end(); ++it) {
				p = key_func(*it);

				assert((0 <= p) && ((size_t)p < num_buckets));

				i2o.emplace_back(p);
				++bucket_counts[p];
			}

			// get offsets of where buckets start (= exclusive prefix sum)
			// iterator once.
			std::vector<std::size_t> offset(bucket_counts.size(), 0);
			for (size_t i = 1, j = 0; i < num_buckets; ++i, ++j) {
				offset[i] = offset[j] + bucket_counts[j];
			}


			// [2nd pass]: saving elements into correct position, and save the final position.
		    std::vector<T> tmp_result(input.size());
			for (size_t i = 0; i < input.size(); ++i) {
					tmp_result[offset[i2o[i]]++] = input[i];
			}

			std::cout << "bucket64 SIZES i2o " << sizeof(i2o) << " tmp results " << sizeof(tmp_result) << std::endl;

			input.swap(tmp_result);
			return bucket_counts;
	    } else if ((num_buckets >> 16) > 0) {
			std::vector<uint32_t> i2o;
			i2o.reserve(input.size());

			// [1st pass]: compute bucket counts and input to bucket assignment.
			uint32_t p;
			for (auto it = input.begin(); it != input.end(); ++it) {
				p = key_func(*it);

				assert((0 <= p) && ((size_t)p < num_buckets));

				i2o.emplace_back(p);
				++bucket_counts[p];
			}

			// get offsets of where buckets start (= exclusive prefix sum)
			// iterator once.
			std::vector<std::size_t> offset(bucket_counts.size(), 0);
			for (size_t i = 1, j = 0; i < num_buckets; ++i, ++j) {
				offset[i] = offset[j] + bucket_counts[j];
			}

			// [2nd pass]: saving elements into correct position, and save the final position.
		    std::vector<T> tmp_result(input.size());
			for (size_t i = 0; i < input.size(); ++i) {
		        tmp_result[offset[i2o[i]]++] = input[i];
			}
			std::cout << "bucket32 SIZES i2o " << sizeof(i2o) << " tmp results " << sizeof(tmp_result) << std::endl;
			input.swap(tmp_result);

			return bucket_counts;

	    } else if ((num_buckets >> 8) > 0) {
			std::vector<uint16_t> i2o;
			i2o.reserve(input.size());

			// [1st pass]: compute bucket counts and input to bucket assignment.
			uint16_t p;
			for (auto it = input.begin(); it != input.end(); ++it) {
				p = key_func(*it);

				assert((0 <= p) && ((size_t)p < num_buckets));

				i2o.emplace_back(p);
				++bucket_counts[p];
			}

			// get offsets of where buckets start (= exclusive prefix sum)
			// iterator once.
			std::vector<std::size_t> offset(bucket_counts.size(), 0);
			for (size_t i = 1, j = 0; i < num_buckets; ++i, ++j) {
				offset[i] = offset[j] + bucket_counts[j];
			}

			// [2nd pass]: saving elements into correct position, and save the final position.
		    std::vector<T> tmp_result(input.size());
			for (size_t i = 0; i < input.size(); ++i) {
		        tmp_result[offset[i2o[i]]++] = input[i];
			}
			std::cout << "bucket16 SIZES i2o " << sizeof(i2o) << " tmp results " << sizeof(tmp_result) << std::endl;

			input.swap(tmp_result);

			return bucket_counts;

	    } else {
			std::vector<uint8_t> i2o;
			i2o.reserve(input.size());

			// [1st pass]: compute bucket counts and input to bucket assignment.
			uint8_t p;
			for (auto it = input.begin(); it != input.end(); ++it) {
				p = key_func(*it);

				assert((0 <= p) && ((size_t)p < num_buckets));

				i2o.emplace_back(p);
				++bucket_counts[p];
			}

			// get offsets of where buckets start (= exclusive prefix sum)
			// iterator once.
			std::vector<std::size_t> offset(bucket_counts.size(), 0);
			for (size_t i = 1, j = 0; i < num_buckets; ++i, ++j) {
				offset[i] = offset[j] + bucket_counts[j];
			}

			// [2nd pass]: saving elements into correct position, and save the final position.
		    std::vector<T> tmp_result(input.size());
			for (size_t i = 0; i < input.size(); ++i) {
		        tmp_result[offset[i2o[i]]++] = input[i];
			}
			std::cout << "bucket8 SIZES i2o " << sizeof(i2o) << " tmp results " << sizeof(tmp_result) << std::endl;

			input.swap(tmp_result);

			return bucket_counts;
	    }
	}

//	/**
//	 * @brief   bucketing of values into `num_buckets` buckets.  Customized from mxx version.
//	 *
//	 * This particular implementation uses an internal temporary buffer
//	 * of the same size as the input. Thus requiring that amount of additional
//	 * memory space. For a version that doesn't need O(n) additional memory,
//	 * use the (somewhat slower) `bucketing_inplace()` function below.
//	 *
//	 * in addition, it returns an array with the bucketed position of the input entries, i.e.
//	 *   output to input mapping
//	 * 		an array x, such that input[i] = output[x[i]]; and the mapping and input share the same coordinates.
//	 * 		then to reconstruct input, we can walk through the mapping and random access output array,
//	 * 		which should be faster as the writes are consecutive in cache while the reads are random.
//	 * 		also, this is more friendly with emplacing into a vector.
//	 * the mapping is useful to avoid recomputing key_func.
//	 *
//	 * THIS VERSION ADDRESSES BELOW.  It's about 50% slower than the version above.
//	 * the question is whether that'd be faster for bucketing as well.  we can keep an array of first entries
//	 *      and a separate array of offsets (for input) from last bucket entry.  together with the bucket size,
//	 *      can then be used to sequentially create output, with 1 random load for offsets, and 1 random load for input to read,
//	 *      and 1 seq write per input entry.
//	 *
//	 *  this is about 50% faster than mxx, but uses 8N more memory, and the ordered write is not helping.
//	 *
//	 * @tparam T            Input type
//	 * @tparam Func         Type of the key function.
//	 * @param input[in|out] Contains the values to be bucketed.
//	 *                      This vector is both input and output.
//	 * @param key_func      A function taking a type T and returning the bucket index
//	 *                      in the range [0, num_buckets).
//	 * @param num_buckets   The total number of buckets.
//	 *
//	 * @return              The number of elements in each bucket.
//	 *                      The size of this vector is `num_buckets`.
//	 *                      The source position for each entry in the bucketed input.
//	 *                      the size of this vector is size of input.
//	 */
//	template <typename T, typename Func>
//	std::pair<std::vector<size_t>, std::vector<size_t> >
//	bucketing_ordered_write(std::vector<T>& input, Func key_func, size_t num_buckets) {
//	    if (num_buckets == 0)
//	    	return std::make_pair(std::vector<size_t>(), std::vector<size_t>());
//
//		// initialize number of elements per bucket
//	    std::vector<size_t> bucket_counts(num_buckets, 0);
//
//	    // if input is empty, simply return
//	    if (input.size() == 0)
//	        return std::make_pair(bucket_counts, std::vector<size_t>());
//
//	    std::vector<size_t> i2o;   // initially used as a linked list (offset to next entry in same bucket.
//	    i2o.reserve(input.size());
//
//	    // [1st pass]: compute bucket counts, input to bucket assignments, and interleaved linked lists for each bucket.
//	    size_t p;
//	    std::vector<std::size_t> heads(bucket_counts.size(), 0);
//	    std::vector<std::size_t> tails(bucket_counts.size(), 0);
//	    size_t i = 0;
//	    for (auto it = input.begin(); it != input.end(); ++it, ++i) {
//	    	p = key_func(*it);   // bucket mapping
//	        assert((0 <= p) && ((size_t)p < num_buckets));
//
//	        if (bucket_counts[p] == 0) { // first time for this bucket
//	        	heads[p] = i;	// point to first input element for this bucket.
//	        } else {
//	        	i2o[tails[p]] = i;	// record in last bucket entry the next bucket entry's position.
//	        						// since it's a prev position, traversal is generally moving forward.
//	        }
//        	tails[p] = i;
//
//        	++bucket_counts[p];
//	    }
//
//
//	    // [2nd pass]: iterate over buckets, populate output, and generate i2o list.
//	    std::vector<T> tmp_result;
//	    tmp_result.reserve(input.size());
//	    size_t in_pos = 0, tmp;
//	    size_t out_pos = 0;
//	    for (size_t i = 0; i < num_buckets; ++i) {
//
//	    	if (bucket_counts[i] == 0) continue;  // skip empty buckets.
//
//	    	in_pos = heads[i];   // get the head of linked list.
//
//	    	// iterate over rest of bucket.
//	    	for (size_t j = 0; j < bucket_counts[i]; ++j, ++out_pos) {
//				tmp_result.emplace_back(input[in_pos]);
//				tmp = i2o[in_pos];   // get the next position.
//				i2o[in_pos] = out_pos;  // convert to from input poition linked list to input to output mapping.
//				in_pos = tmp;
//	    	}
//
//	    }
//
//	    // replacing input with tmp result buffer and return the number of elements
//	    // in each bucket
//	    input.swap(tmp_result);
//	    return std::make_pair(bucket_counts, i2o);
//	}


	/**
	 * @brief  given a bucketed array and the mapping from unbucketed to bucketed, undo the bucketing.
	 * @details  note that output and i2o are traversed in order (write and read respectively),
	 * 			 while input is randomly read.
	 */
	template <typename T>
	std::vector<T> permute(std::vector<T> const & unbucketed, std::vector<size_t> const & i2o) {
	    // if input is empty, simply return
	    if ((unbucketed.size() == 0) || (i2o.size() == 0))
	        return std::vector<T>();

	    assert(unbucketed.size() == i2o.size());

	    // saving elements into correct position
	    std::vector<T> tmp_result(unbucketed.size());
	    for (size_t i = 0; i < unbucketed.size(); ++i) {
	        tmp_result[i2o[i]] = unbucketed[i];
	    }

	    return tmp_result;
	}

	/**
	 * @brief  given a bucketed array and the mapping from unbucketed to bucketed, undo the bucketing.
	 * @details  note that output and i2o are traversed in order (write and read respectively),
	 * 			 while input is randomly read.
	 */
	template <typename T>
	std::vector<T> unpermute(std::vector<T> const & bucketed, std::vector<size_t> const & i2o) {
	    // if input is empty, simply return
	    if ((bucketed.size() == 0) || (i2o.size() == 0))
	        return std::vector<T>();

	    assert(bucketed.size() == i2o.size());

	    // saving elements into correct position
	    std::vector<T> tmp_result;
		tmp_result.reserve(bucketed.size());
	    for (size_t i = 0; i < bucketed.size(); ++i) {
	        tmp_result.emplace_back(bucketed[i2o[i]]);
	    }

	    return tmp_result;
	}


	  /**
	   * @brief       distribute.  speed version.  no guarantee of output ordering, but actually is same.
	   * @param[IN] vals    vals to distribute.  should already be bucketed.
	   * @param[IN/OUT] send_counts    count of each target bucket of input, and each source bucket after call.
	   *
	   * @return distributed values.
	   */
	  template <typename V>
	  ::std::vector<V> distribute_bucketed(::std::vector<V> const & vals,
			  	  ::std::vector<size_t> & send_counts,
	                                   mxx::comm const &_comm) {
		  ::std::vector<size_t> temp(send_counts);
		  mxx::all2all(send_counts, _comm).swap(send_counts);

	      return mxx::all2allv(vals, temp, _comm);
	  }


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

  	if (_comm.rank() == 0) { printf("start\n");  fflush(stdout); }


        BL_BENCH_START(distribute);
        // distribute (communication part)
        std::vector<size_t> send_counts = mxx::bucketing_inplace(vals, to_rank, _comm.size());
        BL_BENCH_END(distribute, "bucket", vals.size());

    	if (_comm.rank() == 0) { printf("bucket\n");  fflush(stdout); }

    	double mean, stdev;
    	for (auto it = send_counts.begin(), max = send_counts.end(); it != max; ++it) {
    		mean += *it;
    		stdev += (*it) * (*it);
    	}
    	mean /= _comm.size();
    	stdev -= sqrt((stdev / _comm.size()) - (mean * mean));
    	if (_comm.rank() == 0) { printf("mean: %f, stdev %f\n", mean, stdev); fflush(stdout); }


        // using set is okay.
        BL_BENCH_START(distribute);
        // distribute (communication part)
        ::fsc::bucket_reduce(vals, send_counts, sorted_input, reducer);
        BL_BENCH_END(distribute, "reduce", vals.size());

    	if (_comm.rank() == 0) { printf("reduce\n");  fflush(stdout); }


        // distribute (communication part)
        BL_BENCH_COLLECTIVE_START(distribute, "a2a", _comm);
        mxx::all2allv(vals, send_counts, _comm).swap(vals);
        BL_BENCH_END(distribute, "a2a", vals.size());


    	if (_comm.rank() == 0) { printf("a2a\n");  fflush(stdout); }

        BL_BENCH_START(distribute);
        std::vector<size_t> recv_counts= mxx::all2all(send_counts, _comm);
        BL_BENCH_END(distribute, "a2a_counts", vals.size());

    	if (_comm.rank() == 0) { printf("a2acounts\n");  fflush(stdout); }

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
