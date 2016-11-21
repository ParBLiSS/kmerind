/*
 * Copyright 2016 Georgia Institute of Technology
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
 * @file    incremental mxx.hpp
 * @ingroup
 * @author  tpan
 * @brief   extends mxx to send/recv incrementally, so as to minimize memory use and allocation
 * @details
 *
 */

#include <algorithm>
#include <mxx/datatypes.hpp>
#include <mxx/comm.hpp>
#include <mxx/collective.hpp>

#include "utils/benchmark_utils.hpp"
#include "utils/function_traits.hpp"

#include "containers/fsc_container_utils.hpp"

namespace imxx
{

  // local version of MPI
  namespace local {



    /* ====
     *
     * options:
     * 	3 possibilities
     *
     * 	1x hashing	|	stable	|	in place	|	algo			|	comment
     * 	======================================================================================================
     * 	n			      |	n		|	n			|					|	we can be stable for free with non-in-place, so this does not make sense.  also 2x hash is hard to offset in computation
     * 	n			      |	n		|	y			|	mxx inplace		|	not stable
     * 	n			      |	y		|	n			|	mxx 			|
     * 	n			      |	y		|	y			|		 			|	no possible strictly.  to be stable means we need ordering information, especially with multiple overlapping ranges (for each bucket entries.)
     *
     * 		1x hash means we need temp storage O(n) for bucket ids, so can be stable automatically
     *
     * 	y						|	y		|	n			|	tony's       	|	in place means to not have a copy of the data, just index.  strictly speaking it's not in place.
     * 	y						|	y		|	y   	|	tony's  			|	possible, if one-to-one mapping between in and out is available.  (extra O(n) memory)
     *
     *  TODO: minimize memory by comparing temp storage element sizes vs index array size - then choose to do "in place" or "not".
     *
     *
     */



    /**
     * @brief   implementation function for use by assign_and_bucket
     * @details uses the smallest data type (ASSIGN_TYPE) given the number of buckets.
     *          this version requires extra O(n) for the mapping, and extra O(n) for data movement.
     *
     *          for version that calls Func 2x (no additional mapping space) and O(n) for data movement, use mxx::bucket
     *          basically a counting sort impl.
     */
    template <typename T, typename Func, typename ASSIGN_TYPE, typename SIZE>
    void
    bucketing_impl(std::vector<T>& input,
                           Func const & key_func,
                           ASSIGN_TYPE const num_buckets,
                           std::vector<SIZE> & bucket_sizes,
                           size_t first = 0,
                           size_t last = std::numeric_limits<size_t>::max()) {

      static_assert(::std::is_integral<ASSIGN_TYPE>::value, "ASSIGN_TYPE should be integral, preferably unsigned");

      bucket_sizes.clear();

      // no bucket.
      if (num_buckets == 0) return;

      // initialize number of elements per bucket
      bucket_sizes.resize(num_buckets, 0);

      // ensure valid range
      size_t f = std::min(first, input.size());
      size_t l = std::min(last, input.size());
      assert((f <= l) && "first should not exceed last" );

      if (f == l) return;  // no data in question.

      size_t len = l - f;

      // single bucket.
      if (num_buckets == 1) {
        bucket_sizes[0] = len;

        return;
      }

      // output to input mapping
      std::vector<ASSIGN_TYPE> i2o;
      i2o.reserve(len);

      // [1st pass]: compute bucket counts and input to bucket assignment.
      ASSIGN_TYPE p;
      for (size_t i = f; i < l; ++i) {
          p = key_func(input[i]);

          assert(((0 <= p) && ((size_t)p < num_buckets)) && "assigned bucket id is not valid");

          i2o.emplace_back(p);
          ++bucket_sizes[p];
      }

      // get offsets of where buckets start (= exclusive prefix sum)
      // use bucket_sizes temporarily.
      bucket_sizes.back() = len - bucket_sizes.back();
      // checking for >= 0 with unsigned and decrement is a bad idea.  need to check for > 0
      for (size_t i = num_buckets - 1; i > 0; --i) {
        bucket_sizes[i-1] = bucket_sizes[i] - bucket_sizes[i-1];
      }
      assert(bucket_sizes.front() == 0);  // first one should be 0 at this point.

      // [2nd pass]: saving elements into correct position, and save the final position.
      std::vector<T> tmp_result(len);
      for (size_t i = f; i < l; ++i) {
          tmp_result[bucket_sizes[i2o[i-f]]++] = input[i];
      }
      //std::cout << "bucket64 SIZES i2o " << sizeof(i2o) << " tmp results " << sizeof(tmp_result) << std::endl;
      if (len == input.size()) input.swap(tmp_result);   // if can swap, swap
      else memcpy(input.data() + f, tmp_result.data(), len * sizeof(T));   // else memcpy.

      // this process should have turned bucket_sizes to an inclusive prefix sum
      assert(bucket_sizes.back() == len);
      // convert inclusive prefix sum back to counts.
      // hopefully this is a fast process when compared to allocating memory.
      for (size_t i = num_buckets - 1; i > 0; --i) {
        bucket_sizes[i] -= bucket_sizes[i-1];
      }
    }

    template <typename T, typename Func, typename ASSIGN_TYPE, typename SIZE>
    void
    bucketing_impl(std::vector<T>const & input,
                           Func const & key_func,
                           ASSIGN_TYPE const num_buckets,
                           std::vector<SIZE> & bucket_sizes,
						   std::vector<T> & results,
                           size_t first = 0,
                           size_t last = std::numeric_limits<size_t>::max()) {

      static_assert(::std::is_integral<ASSIGN_TYPE>::value, "ASSIGN_TYPE should be integral, preferably unsigned");
  	assert(((input.size() == 0) || (input.data() != results.data())) &&
  			"input and output should not be the same.");

      bucket_sizes.clear();

      // no bucket.
      if (num_buckets == 0) return;

      // initialize number of elements per bucket
      bucket_sizes.resize(num_buckets, 0);

      // ensure valid range
      size_t f = std::min(first, input.size());
      size_t l = std::min(last, input.size());
      assert((f <= l) && "first should not exceed last" );

      if (f == l) return;  // no data in question.

      size_t len = l - f;

      // single bucket.
      if (num_buckets == 1) {
    	// set output buckets sizes
        bucket_sizes[0] = len;

        // set output values
        memcpy(results.data() + f, input.data() + f, len * sizeof(T));

        return;
      }

      // output to input mapping
      std::vector<ASSIGN_TYPE> i2o;
      i2o.reserve(len);

      // [1st pass]: compute bucket counts and input to bucket assignment.
      ASSIGN_TYPE p;
      for (size_t i = f; i < l; ++i) {
          p = key_func(input[i]);

          assert(((0 <= p) && ((size_t)p < num_buckets)) && "assigned bucket id is not valid");

          i2o.emplace_back(p);
          ++bucket_sizes[p];
      }

      // get offsets of where buckets start (= exclusive prefix sum)
      // use bucket_sizes temporarily.
      bucket_sizes.back() = l - bucket_sizes.back();
      // checking for >= 0 with unsigned and decrement is a bad idea.  need to check for > 0
      for (size_t i = num_buckets - 1; i > 0; --i) {
        bucket_sizes[i-1] = bucket_sizes[i] - bucket_sizes[i-1];
      }
      assert(bucket_sizes.front() == f);  // first one should be 0 at this point.

      // [2nd pass]: saving elements into correct position, and save the final position.
      results.resize(input.size());
      for (size_t i = f; i < l; ++i) {
          results[bucket_sizes[i2o[i-f]]++] = input[i];
      }

      // this process should have turned bucket_sizes to an inclusive prefix sum
      assert(bucket_sizes.back() == l);
      // convert inclusive prefix sum back to counts.
      // hopefully this is a fast process when compared to allocating memory.
      for (size_t i = num_buckets - 1; i > 0; --i) {
        bucket_sizes[i] -= bucket_sizes[i-1];
      }
      bucket_sizes[0] -= f;

    }



    /**
     * @brief   compute the element index mapping between input and bucketed output.
     *
     *  good for a2av, all at once.
     * this is an "deconstructed" version of assign_to_buckets
     *
     * It returns an array with the output position for the input entries, i.e.
     *    an array x, such that input[i] = output[x[i]]; and the mapping and input share the same coordinates.
     *
     *    To reconstruct input, we can walk through the mapping and random access output array,
     *    which should be no slower than before as the writes are consecutive in cache while the reads are random.
     *
     *    also, this is more friendly with emplacing into a vector.
     *    the mapping is useful to avoid recomputing key_func.
     *
     * this is faster than mxx by about 2x (for hashing k-mers 1x.), but uses another vector of size 8N.
     *
     *    ** operate only within the input range [first, last), and assign to same output range.
     *
     * @tparam T            Input type
     * @tparam Func         Type of the key function.
     * @param input         Contains the values to be bucketed.
     * @param key_func      A function taking a type T and returning the bucket index
     *                      in the range [0, num_buckets).
     * @param num_buckets     number of buckets.  cannot be more than num procs in MPI
     * @param bucket_sizes[out]  The number of elements in each bucket.  reset during each call.
     * @param i2o[out]      input to output position mapping.  updated between first and last.  value is between first and last
     * @param first         first position to bucket in input
     * @param last          last position to bucket in input
     */
    template <typename T, typename Func, typename SIZE = size_t>
    void
    assign_to_buckets(std::vector<T> const & input,
                      Func const & key_func,
                      size_t const & num_buckets,
                      std::vector<SIZE> & bucket_sizes,
                      std::vector<SIZE> & i2o,
                      size_t first = 0,
                      size_t last = std::numeric_limits<size_t>::max()) {
        bucket_sizes.clear();

        // no bucket.
        if (num_buckets == 0) throw std::invalid_argument("ERROR: number of buckets is 0");

        // initialize number of elements per bucket
        bucket_sizes.resize(num_buckets, 0);

        // ensure valid range
        size_t f = std::min(first, input.size());
        size_t l = std::min(last, input.size());
        assert((f <= l) && "first should not exceed last" );

        if (f == l) return;  // no data in question.

        i2o.resize(input.size(), 0);  // need exactly input size.

        // single bucket.
        if (num_buckets == 1) {
          bucket_sizes[0] = l - f;

          memset(&(i2o[f]), 0, bucket_sizes[0] * sizeof(SIZE));

          return;
        }


        // [1st pass]: compute bucket counts and input2bucket assignment.
        // store input2bucket assignment in i2o temporarily.
        size_t p;
        for (size_t i = f; i < l; ++i) {
            p = key_func(input[i]);

            assert(((0 <= p) && ((size_t)p < num_buckets)) && "assigned bucket id is not valid");

            i2o[i] = p;
            ++bucket_sizes[p];
        }

    }

    // operate only within the input range [first, last), and assign to same output range.
    template <typename SIZE = size_t>
    void
    bucket_to_permutation(std::vector<SIZE> & bucket_sizes,
                          std::vector<SIZE> & i2o,
                          size_t first = 0,
                          size_t last = std::numeric_limits<size_t>::max()) {

      // no bucket.
      if (bucket_sizes.size() == 0) throw std::invalid_argument("bucket_sizes has 0 buckets.");

      // ensure valid range
      size_t f = std::min(first, i2o.size());
      size_t l = std::min(last, i2o.size());
      assert((f <= l) && "first should not exceed last" );

      if (f == l) return;  // no data in question.


      // get offsets of where buckets start (= exclusive prefix sum), offset by the range.
      // use bucket_sizes temporarily.
      bucket_sizes.back() = l - bucket_sizes.back();  // offset from f, or alternatively, l.
      // checking for >= 0 with unsigned and decrement is a bad idea.  need to check for > 0
      for (size_t i = bucket_sizes.size() - 1; i > 0; --i) {
        bucket_sizes[i-1] = bucket_sizes[i] - bucket_sizes[i-1];
      }
      assert((bucket_sizes.front() == f) && "prefix sum resulted in incorrect starting position");  // first one should be 0 at this point.

      // [2nd pass]: saving elements into correct position, and save the final position.
      for (size_t i = f; i < l; ++i) {
        i2o[i] = bucket_sizes[i2o[i]]++;  // output position from bucket id (i2o[i]).  post increment.  already offset by f.
      }

      // this process should have turned bucket_sizes to an inclusive prefix sum
      assert((bucket_sizes.back() == l) && "mapping assignment resulted in incorrect ending position");
      // convert inclusive prefix sum back to counts.
      // hopefully this is a fast process when compared to allocating memory.
      for (size_t i = bucket_sizes.size() - 1; i > 0; --i) {
        bucket_sizes[i] -= bucket_sizes[i-1];
      }
      bucket_sizes[0] -= f;
    }

//    template <typename SIZE = size_t>
//    void
//    permutation_to_bucket(std::vector<SIZE> & bucket_sizes,
//                          std::vector<SIZE> & i2o,
//                          size_t first = 0,
//                          size_t last = std::numeric_limits<size_t>::max()) {
//
//      // TO IMPLEMENT?
//    }




    //===  permute and unpermute are good for communication bucketing for partitions of INPUT data, using A2Av for communication (not for A2A)

    /**
     * @brief  given a bucketed array and the mapping from unbucketed to bucketed, undo the bucketing.
     * @details  note that output and i2o are traversed in order (write and read respectively),
     *       while input is randomly read.
     *
     *       permute a range in input.  (entire output range may be affected.)
     *
     *		 if OT is not random access iterator, this will be very slow.
     *
     * @param unbucketed	input to be bucketed.  only processed between first and last.
     * @param i2o			mapping.  only the part between first and last is used.  values should be between first and last.
     * @param results		bucketed output, only processed between first and last.
     */
    template <typename IT, typename OT, typename MT>
    void permute_for_input_range(IT unbucketed, IT unbucketed_end,
    		MT i2o, OT bucketed, OT bucketed_end,
    		size_t const & bucketed_pos_offset) {

    	static_assert(std::is_same<typename std::iterator_traits<IT>::value_type,
    			typename std::iterator_traits<OT>::value_type>::value,
				"ERROR: IT and OT should be iterators with same value types");
    	static_assert(std::is_integral<typename std::iterator_traits<MT>::value_type>::value,
				"ERROR: MT should be an iterator of integral type value");

    	size_t in_len = std::distance(unbucketed, unbucketed_end);  // input range.
    	size_t out_len = std::distance(bucketed, bucketed_end);   // number of output to look for.

        // if input is empty, simply return
        if ((in_len == 0) || (out_len == 0)) return;

        assert((out_len >= in_len) && "input range is larger than output range");

        size_t bucketed_max = bucketed_pos_offset + out_len - 1;
        assert((*(std::min_element(i2o, i2o + in_len)) >= bucketed_pos_offset) &&
        		(*(std::max_element(i2o, i2o + in_len)) <= bucketed_max) &&
				"ERROR, i2o [0, len) does not map to itself");

        // saving elements into correct position
        for (; unbucketed != unbucketed_end; ++i2o, ++unbucketed) {
            *(bucketed + (*i2o - bucketed_pos_offset)) = *unbucketed;
        }
    }


    /**
     * @brief  given a bucketed array and the mapping from unbucketed to bucketed, undo the bucketing.
     * @details  note that output and i2o are traversed in order (write and read respectively),
     *       while input is randomly read.
     *
     *       permute a range in output.  (entire input range may be read..)
     *
     * @param unbucketed	input to be bucketed.  only processed between first and last.
     * @param i2o			mapping.  only the part between first and last is used.  values should be between first and last.
     * @param results		bucketed output, only processed between first and last.
     */
    template <typename IT, typename OT, typename MT>
    void permute_for_output_range(IT unbucketed, IT unbucketed_end,
    		MT i2o, OT bucketed, OT bucketed_end,
    		size_t const & bucketed_pos_offset) {

    	static_assert(std::is_same<typename std::iterator_traits<IT>::value_type,
    			typename std::iterator_traits<OT>::value_type>::value,
				"ERROR: IT and OT should be iterators with same value types");
    	static_assert(std::is_integral<typename std::iterator_traits<MT>::value_type>::value,
				"ERROR: MT should be an iterator of integral type value");

    	size_t in_len = std::distance(unbucketed, unbucketed_end);  // input range.
    	size_t out_len = std::distance(bucketed, bucketed_end);   // number of output to look for.

        // if input is empty, simply return
        if ((in_len == 0) || (out_len == 0)) return;

        assert((in_len >= out_len) && "output range is larger than input range");

        size_t bucketed_max = bucketed_pos_offset + out_len - 1;
        assert((*(std::min_element(i2o, i2o + in_len)) <= bucketed_pos_offset) &&
        		(*(std::max_element(i2o, i2o + in_len)) >= bucketed_max) &&
				"ERROR, i2o [0, len) does not map to itself");

        // saving elements into correct position
        size_t pos = 0;
        size_t count = 0;
        for (; unbucketed != unbucketed_end; ++i2o, ++unbucketed) {
        	pos = *i2o;
        	if ((pos < bucketed_pos_offset) || (pos > bucketed_max)) continue;

			*(bucketed + (pos - bucketed_pos_offset)) = *unbucketed;
			++count;

        	if (count == out_len) break;  // filled.  done.
        }
    }

    /**
     * @brief  given a bucketed array and the mapping from unbucketed to bucketed, undo the bucketing.
     * @details  note that output and i2o are traversed in order (write and read respectively),
     *       while input is randomly read.
     *
     *       permute a range in input.  (the range maps to itself.)
     *
     * @param unbucketed	input to be bucketed.  only processed between first and last.
     * @param i2o			mapping.  only the part between first and last is used.  values should be between first and last.
     * @param results		bucketed output, only processed between first and last.
     */
    template <typename IT, typename OT, typename MT>
    void permute(IT unbucketed, IT unbucketed_end,
    		MT i2o, OT bucketed,
    		size_t const & bucketed_pos_offset) {

    	static_assert(std::is_same<typename std::iterator_traits<IT>::value_type,
    			typename std::iterator_traits<OT>::value_type>::value,
				"ERROR: IT and OT should be iterators with same value types");
    	static_assert(std::is_integral<typename std::iterator_traits<MT>::value_type>::value,
				"ERROR: MT should be an iterator of integral type value");

    	size_t len = std::distance(unbucketed, unbucketed_end);
        // if input is empty, simply return
        if (len == 0) return;

        assert((*(std::min_element(i2o, i2o + len)) == bucketed_pos_offset) &&
        		(*(std::max_element(i2o, i2o + len)) == ( bucketed_pos_offset + len - 1)) &&
				"ERROR, i2o [first, last) does not map to itself");

        // saving elements into correct position
        for (; unbucketed != unbucketed_end; ++unbucketed, ++i2o) {
            *(bucketed + (*i2o - bucketed_pos_offset)) = *unbucketed;
        }

    }



    /**
     * @brief inplace permute.  at most 2n steps, as we need to check each entry for completion.
     * @details	also need i2o to be a signed integer, so we can use the sign bit to check
     * 			that an entry has been visited.
     *
     *      permute_inplace uses cycle following, so entire range is affected, unless it is known a priori
     *      that [first, last) maps to itself.
     *
     * 			profiling:  heaptrack shows that there are unlikely unnecessary allocation
     * 						cachegrind and gprof show that bulk of time is in inplace_permute, and there are no further details
     * 						perf stats show that cache miss rate is about the same for permute as is for inplace permute
     * 						perf record shows that the majority of the cost is in swap (or in the case of inplace_unpermute, in copy assignment)
     * 						checking counter shows that the array is visited approximately 2x. reducing this is unlikely to improve performance
     * 							due to move/copy cost.
     * 						move/copy cost is associated with std::pair.
     *				loop unrolling did not help much.  trying to change read/write dependence does not seem to help much either.
     *				issue is random read and write, probably.  perf profiling identifies the problem as such.
     *
     */
    template <typename T>
    void permute_inplace(std::vector<T> & unbucketed, std::vector<size_t> & i2o,
    		size_t first = 0,
			size_t last = std::numeric_limits<size_t>::max()) {

    	if (unbucketed.size() > static_cast<size_t>(std::numeric_limits<int64_t>::max()))
    		throw std::invalid_argument("input is limited to max signed long in size");
    	// TODO: separate logic using a bit vector.

    	// ensure valid range
        size_t f = std::min(first, unbucketed.size());
        size_t l = std::min(last, unbucketed.size());
        assert((f <= l) && "first should not exceed last" );
        size_t len = l-f;

        // if input is empty, simply return
        if ((len == 0) || (unbucketed.size() == 0) || (i2o.size() == 0)) return;

        assert((*(std::min_element(i2o.begin() + f, i2o.begin() + l)) == f) &&
        		(*(std::max_element(i2o.begin() + f, i2o.begin() + l)) == (l - 1)) &&
				"ERROR, i2o [first, last) does not map to itself");

        assert(unbucketed.size() == i2o.size());

        // saving elements into correct position
        // modeled after http://stackoverflow.com/questions/7365814/in-place-array-reordering
        // and https://blog.merovius.de/2014/08/12/applying-permutation-in-constant.html
        size_t j, k;
        T v;
        constexpr size_t mask = ~(0UL) << (sizeof(size_t) * 8 - 1);
//        size_t counter = 0;
        for (size_t i = f; i < l; ++i) {
        	k = i2o[i];
//        	++counter;

        	if (k >= mask) {
        		// previously visited.  unmark and move on.
        		i2o[i] = ~k;  // negate
        		continue;
        	}

        	// get the current position and value
        	j = i;
        	v = unbucketed[j];
        	while (k != i) {  // stop when cycle completes (return to original position)
        		std::swap(unbucketed[k], v);
        		j = k;
        		k = i2o[j];
//        		++counter;

        		i2o[j] = ~k;  			  // mark the last position as visited.
        	}
        	unbucketed[i] = v;				  // cycle complete, swap back the current value (from some other location)
        }

//        std::cout << "i2o read " << counter << " times.  size= " << i2o.size() << std::endl;

    }


    /**
     * @brief  given a bucketed array and the mapping from unbucketed to bucketed, undo the bucketing.
     * @details  note that output and i2o are traversed in order (write and read respectively),
     *       while input is randomly read.
     *
     *       unpermute input range (entire output range is involved, last - first entries are written.)
     *       i2o and unbucketed should be same size.
     */
    template <typename IT, typename OT, typename MT>
    void unpermute_by_input_range(IT bucketed, IT bucketed_end,
    		MT i2o, OT unbucketed, OT unbucketed_end,
    		size_t const & bucketed_pos_offset) {

    	static_assert(std::is_same<typename std::iterator_traits<IT>::value_type,
    			typename std::iterator_traits<OT>::value_type>::value,
				"ERROR: IT and OT should be iterators with same value types");
    	static_assert(std::is_integral<typename std::iterator_traits<MT>::value_type>::value,
				"ERROR: MT should be an iterator of integral type value");

    	size_t in_len = std::distance(bucketed, bucketed_end);   // number of output to look for.
    	size_t out_len = std::distance(unbucketed, unbucketed_end);  // input range.

        // if input is empty, simply return
        if ((in_len == 0) || (out_len == 0)) return;

        assert((out_len >= in_len) && "input range is larger than output range");

        size_t bucketed_max = bucketed_pos_offset + in_len - 1;
        assert((*(std::min_element(i2o, i2o + out_len)) <= bucketed_pos_offset) &&
        		(*(std::max_element(i2o, i2o + out_len)) >= bucketed_max) &&
				"ERROR, i2o [0, len) does not map to itself");

        MT i2o_end = i2o + out_len;

        size_t pos = 0;
        size_t count = 0;
        for (; i2o != i2o_end; ++i2o, ++unbucketed) {
        	pos = *i2o;
        	if ((pos < bucketed_pos_offset) || (pos > bucketed_max)) continue;

        	*unbucketed = *(bucketed + (pos - bucketed_pos_offset));
			++count;

        	if (count == in_len) break;  // filled.  done.
        }

    }


    /**
     * @brief  given a bucketed array and the mapping from unbucketed to bucketed, undo the bucketing.
     * @details  note that output and i2o are traversed in order (write and read respectively),
     *       while input is randomly read.
     *
     *       unpermute output range (entire input range is involved, last - first entries are read.)
     */
    template <typename IT, typename OT, typename MT>
    void unpermute_by_output_range(IT bucketed, IT bucketed_end,
    		MT i2o, OT unbucketed, OT unbucketed_end,
    		size_t const & bucketed_pos_offset) {

    	static_assert(std::is_same<typename std::iterator_traits<IT>::value_type,
    			typename std::iterator_traits<OT>::value_type>::value,
				"ERROR: IT and OT should be iterators with same value types");
    	static_assert(std::is_integral<typename std::iterator_traits<MT>::value_type>::value,
				"ERROR: MT should be an iterator of integral type value");

    	size_t out_len = std::distance(unbucketed, unbucketed_end);  // input range.
    	size_t in_len = std::distance(bucketed, bucketed_end);   // number of output to look for.

        // if input is empty, simply return
        if ((in_len == 0) || (out_len == 0)) return;

        assert((out_len <= in_len) && "input range is larger than output range");

        size_t bucketed_max = bucketed_pos_offset + in_len - 1;
        assert((*(std::min_element(i2o, i2o + out_len)) >= bucketed_pos_offset) &&
        		(*(std::max_element(i2o, i2o + out_len)) <= bucketed_max) &&
				"ERROR, i2o [0, len) does not map to itself");

        // saving elements into correct position
        for (; unbucketed < unbucketed_end; ++i2o, ++unbucketed) {
            *unbucketed = *(bucketed + (*i2o - bucketed_pos_offset));
        }
    }

    /**
     * @brief  given a bucketed array and the mapping from unbucketed to bucketed, undo the bucketing.
     * @details  note that output and i2o are traversed in order (write and read respectively),
     *       while input is randomly read.
     *
     *       unpermute output range (requires [first, last) map to itself.)
     */
    template <typename IT, typename OT, typename MT>
    void unpermute(IT bucketed, IT bucketed_end,
    		MT i2o, OT unbucketed,
    		size_t const & bucketed_pos_offset) {

    	static_assert(std::is_same<typename std::iterator_traits<IT>::value_type,
    			typename std::iterator_traits<OT>::value_type>::value,
				"ERROR: IT and OT should be iterators with same value types");
    	static_assert(std::is_integral<typename std::iterator_traits<MT>::value_type>::value,
				"ERROR: MT should be an iterator of integral type value");

    	size_t len = std::distance(bucketed, bucketed_end);
        // if input is empty, simply return
        if (len == 0) return;

        assert((*(std::min_element(i2o, i2o + len)) == bucketed_pos_offset) &&
        		(*(std::max_element(i2o, i2o + len)) == ( bucketed_pos_offset + len - 1)) &&
				"ERROR, i2o [first, last) does not map to itself");

        // saving elements into correct position
        OT unbucketed_end = unbucketed + len;
        for (; unbucketed != unbucketed_end; ++unbucketed, ++i2o) {
            *unbucketed = *(bucketed + (*i2o - bucketed_pos_offset));
        }
    }


    /**
     *
    * 			profiling:  heaptrack shows that there are unlikely unnecessary allocation
    * 						cachegrind and gprof show that bulk of time is in inplace_permute, and there are no further details
    * 						perf stats show that cache miss rate is about the same for permute as is for inplace permute
    * 						perf record shows that the majority of the cost is in swap (or in the case of inplace_unpermute, in copy assignment)
    * 						checking counter shows that the array is visited approximately 2x. reducing this is unlikely to improve performance
    * 							due to move/copy cost.
    * 						move/copy cost is associated with std::pair.
    *
    * 						unpermute input range (entire bucketed range is used.)
    *
     */
    template <typename T>
    void unpermute_inplace(std::vector<T> & bucketed, std::vector<size_t> & i2o,
    		size_t first = 0,
    					size_t last = std::numeric_limits<size_t>::max()) {

    	if (bucketed.size() > static_cast<size_t>(std::numeric_limits<int64_t>::max()))
    		throw std::invalid_argument("input is limited to max signed long in size");
    	// TODO: separate logic using a bit vector.

    	// ensure valid range
        size_t f = std::min(first, bucketed.size());
        size_t l = std::min(last, bucketed.size());
        assert((f <= l) && "first should not exceed last" );
        size_t len = l-f;

        // if input is empty, simply return
        if ((len == 0) || (bucketed.size() == 0) || (i2o.size() == 0)) return;

        assert(bucketed.size() == i2o.size());

        // saving elements into correct position
        // modeled after http://stackoverflow.com/questions/7365814/in-place-array-reordering
        // and https://blog.merovius.de/2014/08/12/applying-permutation-in-constant.html
        size_t j, k;
        T v;
//        size_t counter = 0;
        constexpr size_t mask = ~(0UL) << (sizeof(size_t) * 8 - 1);
        for (size_t i = f; i < l; ++i) {
        	k = i2o[i];				// position that the unbucketed was moved to in the bucketed list; source data
//        	++counter;

        	if (k >= mask) {
        		// previously visited.  unmark and move on.
        		i2o[i] = ~k;  // negate
        		continue;
        	}

        	// get the current position and value
        	v = bucketed[i];		// curr value - to be moved forward until its correct position is found
        	j = i;  				// curr position in bucketed list == , where the value is to be replaced.
        	while (k != i) {  // stop when cycle completes (return to original position)
//        	while ((k & mask) == 0) {  // avoid previously visited, including i.  each entry visited just once.  THIS DOES NOT WORK - k may not be set to negative yet.
        		// copy assign operator of the data type.  could be costly.
        		bucketed[j] = bucketed[k];   // move value back from the bucketed list to the unbucketed.
        		i2o[j] = ~k;  			  // mark the curr position as visited.

        		j = k;  				  // advance the curr position
        		k = i2o[j];				  // get the next position
//        		++counter;
        	}
        	bucketed[j] = v;		// cycle complete, swap back the current value (from some other location)
        	i2o[j] = ~k;		// cycle complete

        	i2o[i] = ~i2o[i];   // unmark the i2o map entry- never visiting prev i again.
        }

//        std::cout << "i2o read " << counter << " times.  size= " << i2o.size() << std::endl;
    }


    /**
     * @brief perform permutation on i2o to produce x number of blocks, each with p number of buckets with specified per-bucket size s.  remainder are placed at end
     *        and output bucket sizes are the remainders for each bucket.
     * @details  kind of a division on the buckets.
     *
     *        operates within range [first, last) only.
     *
     * @param in_block_bucket_size   size of a bucket inside a block - all buckets are same size.
     * @param bucket_sizes[in/out]   entire size of local buckets.
     * @param i2o[in/out]            from bucket assignment, to permutation index.
     * @return number of FULL blocks in the permutation.
     */
    template <typename SIZE = size_t>
    void
    bucket_to_block_permutation(SIZE const & block_bucket_size, SIZE const & nblocks,
                                std::vector<SIZE> & bucket_sizes,
                                std::vector<SIZE> & i2o,
                                size_t first = 0,
                                size_t last = std::numeric_limits<size_t>::max()) {

      if (bucket_sizes.size() == 0)
        throw std::invalid_argument("bucket_sizes is 0");


      // if block_bucket_size is 0, use normal stuff
      if ((block_bucket_size * nblocks * bucket_sizes.size()) > i2o.size()) {
        throw std::invalid_argument("block is larger than total size");
      }

        // no block.  so do it without blocks
      if ((block_bucket_size == 0) || (nblocks == 0)) {
        // no blocks. use standard permutation
        bucket_to_permutation(bucket_sizes, i2o, first, last);

        return;
      }

    	if (bucket_sizes.size() == 1)
    	{
    	  // all go into the block.  no left overs.
    		bucket_to_permutation(bucket_sizes, i2o, first, last);
    		bucket_sizes[0] = i2o.size() - block_bucket_size * nblocks;
    		return;
    	}

      SIZE min_bucket_size = *(std::min_element(bucket_sizes.begin(), bucket_sizes.end()));  // find the minimum
      if (nblocks > (min_bucket_size / block_bucket_size)) {
        // no blocks. use standard permutation
        throw std::invalid_argument("min bucket is smaller than the number of requested blocks. ");
      }

      // else we do blocked permutation

      size_t block_size = block_bucket_size * bucket_sizes.size();
      size_t front_size = block_size * nblocks;

      // ensure valid range
      size_t f = std::min(first, i2o.size());
      size_t l = std::min(last, i2o.size());
      assert((f <= l) && "first should not exceed last" );

      if (f == l) return;  // no data in question.


      // initialize the offset, and reduce the bucket_sizes
      std::vector<SIZE> offsets(bucket_sizes.size(), 0);
      for (size_t i = 0; i < bucket_sizes.size(); ++i) {
        offsets[i] = i * block_bucket_size;
        bucket_sizes[i] -= block_bucket_size * nblocks;
      }

      //== convert bucket sizes to offsets as well.
      // get offsets of where buckets start (= exclusive prefix sum), offset by the range.
      // use bucket_sizes temporarily.
      bucket_sizes.back() = l - bucket_sizes.back();  // offset from f, or alternatively, l.
      // checking for >= 0 with unsigned and decrement is a bad idea.  need to check for > 0
      for (size_t i = bucket_sizes.size() - 1; i > 0; --i) {
        bucket_sizes[i-1] = bucket_sizes[i] - bucket_sizes[i-1];
      }
      assert((bucket_sizes.front() == (f + front_size)) && "prefix sum resulted in incorrect starting position");  // first one should be 0 at this point.


      // walk through all.
      size_t bucket;
      for (size_t i = f; i < l; ++i) {
        bucket = i2o[i];
        i2o[i] = offsets[bucket]++;  // output position from bucket id (i2o[i]).  post increment.  already offset by f.

        // if the offset entry indicates full, move the offset to next block, unless we already have max number of blocks.
        if (offsets[bucket] > front_size) {
          continue;  // last part with NONBLOCKS, so let it continue;
//        } else if (offsets[bucket] == front_size) {  // last element of last bucket in last block in the front portion
//          offsets[bucket] = bucket_sizes[bucket];
        } else if (offsets[bucket] % block_bucket_size == 0) {  // full bucket.
          if (((offsets[bucket] + block_size - 1) / block_size) == nblocks) { // in last block
            offsets[bucket] = bucket_sizes[bucket];
          } else {
            offsets[bucket] += block_size - block_bucket_size;
          }
        }
      }

      // convert inclusive prefix sum back to counts.
      // hopefully this is a fast process when compared to allocating memory.
      for (size_t i = 0; i < bucket_sizes.size() - 1; ++i) {
        bucket_sizes[i] = bucket_sizes[i+1] - bucket_sizes[i];
      }
      bucket_sizes.back() = l - bucket_sizes.back();

    }

  } // local namespace


  //== get bucket assignment. - complete grouping by processors (for a2av), or grouping by processors in communication blocks (for a2a, followed by 1 a2av)
  // parameter: in block per bucket send count (each block is for an all2all operation.  max block size is capped by max int.
  // 2 sets of send counts


  // need to 1. get assignment + bucket count
  //         2. globally compute block size
  //         3. do a blocks at a time ()
  // once the block permutation matrix is available, we need to permute the input completely (e.g. using the inplace permute and unpermute algorithm),
  //     or gather the part we need into local buffers (inplace gather and scatter does NOT make sense.  also, local buffers allows inplace MPI_alltoall)
  //

  /*
     *        from this data can be processed in a few ways:
     *          1. permute wholly in place, then communicate (wholly or in blocks) not in place.  up to 2x memory.
     *              communicate wholly is wasteful because of complete intermediate buffer is needed but may be faster.
     *              comm in blocks saves memory but may be slower.
     *              origin not preserved, permuted is preserved.  okay for count and exists.
     *              permute slower, comm faster
     *          2. permute wholly in place, then communicate (wholly or in blocks) in place for the first blocks.
     *              original not preserved. permuted not preserved.  okay for erase, update, insert, even find.
     *              permute slower, comm probably slower
     *          3. permute wholly not in place, then communicate not in place - might be faster.  3x memory now.
     *              original and permuted preserved.  output may be in permuted order.
     *              permute fast, comm fast
     *          4. permute wholly not in place, then communicate in place for first blocks.
     *              original preserved, permuted not preserved.  no ordering for count and exists
     *              permute fast, comm probably slower.
     *          5. gather locally then communicate .  still need local buffer.  probably slower than permute wholly in place since we need to search.  orig order preserved.
   *
   */




  /**
   * @brief       distribute.  speed version.  no guarantee of output ordering, but actually is same.
   * @param[IN] vals    vals to distribute.  should already be bucketed.
   * @param[IN/OUT] send_counts    count of each target bucket of input, and each source bucket after call.
   *
   * @return distributed values.
   */
//      template <typename V>
//      ::std::vector<V> distribute_bucketed(::std::vector<V> const & vals,
//              ::std::vector<size_t> & send_counts,
//                                       mxx::comm const &_comm) {
//        ::std::vector<size_t> temp(send_counts);
//        mxx::all2all(send_counts, _comm).swap(send_counts);
//
//          return mxx::all2allv(vals, temp, _comm);
//      }



// specialized for lower memory, incremental operations
// also  later, non-blocking mxx.
// requirements:  in place in original or fixed-size buffer(s)
// iterative processing




//== bucket actually.

//== all2allv using in place all2all as first step, then all2allv + buffer as last part.  pure communication.
// to use all2all, the data has to be contiguous.  2 ways to do this:  1. in place.  2. separate buffer
  /**
   * @brief use all2all to send 1 block of data, block = send_count * comm.size()
   *
   * @param input should be already bucketed and arranged as blocks with p * b, where b is size of each bucket.
   * @param send_count   number to send to each processor
   * @param output    buffer to store the result
   * @param send_offset    offset in input from which to start sending.  end of range is offset + send_count * comm.size()
   * @param recv_offset    offset in output from which to start writing.  end of range is offset + send_count * comm.size()
   * @param comm    communicator
   */
  template <typename T>
  void block_all2all(std::vector<T> const & input, size_t send_count, std::vector<T> & output,
                     size_t send_offset = 0, size_t recv_offset = 0, mxx::comm const & comm = mxx::comm()) {

	    bool empty = ((input.size() == 0) || (send_count == 0));
	    empty = mxx::all_of(empty);
	    if (empty) {
	      return;
	    }

    // ensure enough space
    assert((input.size() >= (send_offset + send_count * comm.size())) && "input for block_all2all not big enough");
    assert((output.size() >= (recv_offset + send_count * comm.size())) && "output for block_all2all not big enough");

    // send via mxx all2all - leverage any large message support from mxx.  (which uses datatype.contiguous() to increase element size and reduce element count to 1)
    ::mxx::all2all(&(input[send_offset]), send_count, &(output[recv_offset]), comm);
  }

  /**
   * @brief use all2all and all2allv to send data in blocks
   * @details   first part of the vector is sent via alltoall.  second part sent by all2allv.
   *            sends 1 block.
   *
   * @param input[in|out] should be already bucketed and arranged as blocks with p * b, where b is size of each bucket.
   * @param send_count   number to send to each processor
   * @param send_offset    offset in input from which to start sending.  end of range is offset + send_count * comm.size()
   * @param comm    communicator
   */
  template <typename T>
  void block_all2all_inplace(std::vector<T> & input, size_t send_count,
                     size_t offset = 0, mxx::comm const & comm = mxx::comm()) {

	    bool empty = ((input.size() == 0) || (send_count == 0));
	    empty = mxx::all_of(empty);
	    if (empty) {
	      return;
	    }


    // ensure enough space
    assert((input.size() >= (offset + send_count * comm.size())) && "input for block_all2all not big enough");
    assert((send_count < mxx::max_int) && "send count too large for mpi ");
//    assert((send_count * comm.size() < mxx::max_int) && "send count too large for mpi ");

    // get datatype based on size.
    mxx::datatype dt;
    if ((send_count * comm.size()) < mxx::max_int) {
      dt = mxx::get_datatype<T>();
    } else {
      // create a contiguous data type to support large messages (send_count * comm.size() > mxx::max_int).
      dt = mxx::get_datatype<T>().contiguous(send_count);
      send_count = 1;
    }

    // send using special mpi keyword MPI_IN_PLACE. should work for MPI_Alltoall.
    MPI_Alltoall(MPI_IN_PLACE, send_count, dt.type(), const_cast<T*>(&(input[offset])), send_count, dt.type(), comm);
  }


  /**
   * @brief distribute function.  input is transformed, but remains the original input with original order.  buffer is used for output.
   * @details
   * @tparam SIZE     type for the i2o mapping and recv counts.  should be large enough to represent max of input.size() and output.size()
   */
  template <typename V, typename ToRank, typename SIZE>
  void distribute(::std::vector<V>& input, ToRank const & to_rank,
                  ::std::vector<SIZE> & recv_counts,
                  ::std::vector<SIZE> & i2o,
                  ::std::vector<V>& output,
                  ::mxx::comm const &_comm, bool const & preserve_input = true) {
    BL_BENCH_INIT(distribute);

    BL_BENCH_COLLECTIVE_START(distribute, "empty", _comm);
    bool empty = input.size() == 0;
    empty = mxx::all_of(empty);
    BL_BENCH_END(distribute, "empty", input.size());

    if (empty) {
      BL_BENCH_REPORT_MPI_NAMED(distribute, "map_base:distribute", _comm);
      return;
    }
    // speed over mem use.  mxx all2allv already has to double memory usage. same as stable distribute.

    BL_BENCH_START(distribute);
    std::vector<SIZE> send_counts(_comm.size(), 0);
    i2o.resize(input.size());
    BL_BENCH_END(distribute, "alloc_map", input.size());

    // bucketing
    BL_BENCH_START(distribute);
    imxx::local::assign_to_buckets(input, to_rank, _comm.size(), send_counts, i2o, 0, input.size());
    BL_BENCH_END(distribute, "bucket", input.size());

    BL_BENCH_START(distribute);
    // distribute (communication part)
    imxx::local::bucket_to_permutation(send_counts, i2o, 0, input.size());
    BL_BENCH_END(distribute, "to_pos", input.size());

    BL_BENCH_START(distribute);
    if (output.capacity() < input.size()) output.clear();
    output.resize(input.size());
    output.swap(input);  // swap the 2.
    BL_BENCH_END(distribute, "alloc_permute", output.size());

    BL_BENCH_START(distribute);
    // distribute (communication part)
    imxx::local::permute(output.begin(), output.end(), i2o.begin(), input.begin(), 0);  // input now holds permuted entries.
    BL_BENCH_END(distribute, "permute", input.size());

    // distribute (communication part)
    BL_BENCH_START(distribute);
    recv_counts.resize(_comm.size());
    mxx::all2all(send_counts.data(), 1, recv_counts.data(), _comm);
    size_t total = std::accumulate(recv_counts.begin(), recv_counts.end(), static_cast<size_t>(0));
    BL_BENCH_END(distribute, "a2a_count", recv_counts.size());

    BL_BENCH_START(distribute);
    // now resize output
    if (output.capacity() < total) output.clear();
    output.resize(total);
    BL_BENCH_END(distribute, "realloc_out", output.size());

    BL_BENCH_START(distribute);
    mxx::all2allv(input.data(), send_counts, output.data(), recv_counts, _comm);
    BL_BENCH_END(distribute, "a2a", output.size());

    if (preserve_input) {
      BL_BENCH_START(distribute);
      // unpermute.  may be able to work around this so leave it as "_inplace"
      imxx::local::unpermute_inplace(input, i2o, 0, input.size());
      BL_BENCH_END(distribute, "unpermute_inplace", input.size());
    }
    BL_BENCH_REPORT_MPI_NAMED(distribute, "map_base:distribute", _comm);

  }

  template <typename V, typename SIZE>
  void undistribute(::std::vector<V> const & input,
                  ::std::vector<SIZE> const & recv_counts,
                  ::std::vector<SIZE> & i2o,
                  ::std::vector<V>& output,
                  ::mxx::comm const &_comm, bool const & restore_order = true) {
    BL_BENCH_INIT(undistribute);

    BL_BENCH_COLLECTIVE_START(undistribute, "empty", _comm);
    bool empty = input.size() == 0;
    empty = mxx::all_of(empty);
    BL_BENCH_END(undistribute, "empty", input.size());

    if (empty) {
      BL_BENCH_REPORT_MPI_NAMED(undistribute, "map_base:undistribute", _comm);
      return;
    }


    BL_BENCH_START(undistribute);
    std::vector<size_t> send_counts(recv_counts.size());
    mxx::all2all(recv_counts.data(), 1, send_counts.data(), _comm);
    size_t total = std::accumulate(send_counts.begin(), send_counts.end(), static_cast<size_t>(0));
    BL_BENCH_END(undistribute, "recv_counts", input.size());

    BL_BENCH_START(undistribute);
    if (output.capacity() < total) output.clear();
    output.resize(total);
    BL_BENCH_END(undistribute, "realloc_out", output.size());

    BL_BENCH_START(undistribute);
    mxx::all2allv(input.data(), recv_counts, output.data(), send_counts, _comm);
    BL_BENCH_END(undistribute, "a2av", input.size());

    if (restore_order) {
      BL_BENCH_START(undistribute);
      // distribute (communication part)
      imxx::local::unpermute_inplace(output, i2o, 0, output.size());
      BL_BENCH_END(undistribute, "unpermute_inplace", output.size());

    }
    BL_BENCH_REPORT_MPI_NAMED(undistribute, "map_base:undistribute", _comm);

  }

  /**
   * @brief distribute function.  input is transformed, but remains the original input with original order.  buffer is used for output.
   *
   */
  template <typename V, typename ToRank, typename SIZE>
  void distribute_2part(::std::vector<V>& input, ToRank const & to_rank,
                        ::std::vector<SIZE> & recv_counts,
                        ::std::vector<SIZE> & i2o,
                        ::std::vector<V>& output,
                        ::mxx::comm const &_comm, bool const & preserve_input = true) {
    BL_BENCH_INIT(distribute);

      // speed over mem use.  mxx all2allv already has to double memory usage. same as stable distribute.
    BL_BENCH_COLLECTIVE_START(distribute, "empty", _comm);
    bool empty = input.size() == 0;
    empty = mxx::all_of(empty);
    BL_BENCH_END(distribute, "empty", input.size());

    if (empty) {
      BL_BENCH_REPORT_MPI_NAMED(distribute, "map_base:distribute", _comm);
      return;
    }

      // do assignment.
      BL_BENCH_START(distribute);
      std::vector<SIZE> send_counts(_comm.size(), 0);
      i2o.resize(input.size());
      BL_BENCH_END(distribute, "alloc_map", input.size());

      // bucketing
      BL_BENCH_START(distribute);
      imxx::local::assign_to_buckets(input, to_rank, _comm.size(), send_counts, i2o, 0, input.size());
      BL_BENCH_END(distribute, "bucket", input.size());

      // compute minimum block size.
      BL_BENCH_START(distribute);
      SIZE min_bucket_size = *(::std::min_element(send_counts.begin(), send_counts.end()));
      min_bucket_size = ::mxx::allreduce(min_bucket_size, mxx::min<SIZE>(), _comm);
      BL_BENCH_END(distribute, "min_bucket_size", min_bucket_size);

      // compute the permutations from block size and processor mapping.  send_counts modified to the remainders.
      BL_BENCH_START(distribute);
      ::imxx::local::bucket_to_block_permutation(min_bucket_size, 1UL, send_counts, i2o, 0, input.size());
      SIZE first_part = _comm.size() * min_bucket_size;
      BL_BENCH_END(distribute, "to_pos", input.size());

      // compute receive counts and total
      BL_BENCH_START(distribute);
      recv_counts.resize(_comm.size());
      mxx::all2all(send_counts.data(), 1, recv_counts.data(), _comm);
      SIZE total = std::accumulate(recv_counts.begin(), recv_counts.end(), static_cast<size_t>(0));
      total += first_part;
      BL_BENCH_END(distribute, "a2av_count", total);

      BL_BENCH_START(distribute);
      if (output.capacity() < input.size()) output.clear();
      output.resize(input.size());
      output.swap(input);
      BL_BENCH_END(distribute, "alloc_out", output.size());

      // permute
      BL_BENCH_START(distribute);
      imxx::local::permute(output.begin(), output.end(), i2o.begin(), input.begin(), 0);
      BL_BENCH_END(distribute, "permute", input.size());

      BL_BENCH_START(distribute);
      if (output.capacity() < total) output.clear();
      output.resize(total);
      BL_BENCH_END(distribute, "alloc_out", output.size());

      BL_BENCH_START(distribute);
      block_all2all(input, min_bucket_size, output, 0, 0, _comm);
      BL_BENCH_END(distribute, "a2a", first_part);


      BL_BENCH_START(distribute);
      mxx::all2allv(input.data() + first_part, send_counts,
                    output.data() + first_part, recv_counts, _comm);
      BL_BENCH_END(distribute, "a2av", total - first_part);

      // permute
      if (preserve_input) {
        BL_BENCH_START(distribute);
        ::imxx::local::unpermute_inplace(input, i2o, 0, input.size());
        BL_BENCH_END(distribute, "unpermute_inplace", input.size());
      }

      BL_BENCH_REPORT_MPI_NAMED(distribute, "map_base:distribute", _comm);
  }


  /**
   * @param recv_counts  counts for each bucket that is NOT PART OF FIRST BLOCK.
   */
  template <typename V, typename SIZE>
  void undistribute_2part(::std::vector<V> const & input,
                  ::std::vector<SIZE> const & recv_counts,
                  ::std::vector<SIZE> & i2o,
                  ::std::vector<V>& output,
                  ::mxx::comm const &_comm, bool const & restore_order = true) {
    BL_BENCH_INIT(undistribute);

    BL_BENCH_COLLECTIVE_START(undistribute, "empty", _comm);
    bool empty = input.size() == 0;
    empty = mxx::all_of(empty);
    BL_BENCH_END(undistribute, "empty", input.size());

    if (empty) {
      BL_BENCH_REPORT_MPI_NAMED(undistribute, "map_base:undistribute", _comm);
      return;
    }


    BL_BENCH_START(undistribute);
    size_t second_part = std::accumulate(recv_counts.begin(), recv_counts.end(), static_cast<size_t>(0));
    size_t first_part = input.size() - second_part;
    assert((first_part % _comm.size() == 0) && "the first block should be evenly distributed to buckets.");

    std::vector<size_t> send_counts(recv_counts.size());
    mxx::all2all(recv_counts.data(), 1, send_counts.data(), _comm);
    second_part = std::accumulate(send_counts.begin(), send_counts.end(), static_cast<size_t>(0));
    BL_BENCH_END(undistribute, "recv_counts", second_part);

    BL_BENCH_START(undistribute);
    if (output.capacity() < (first_part + second_part)) output.clear();
    output.resize(first_part + second_part);
    BL_BENCH_END(undistribute, "realloc_out", output.size());

    BL_BENCH_START(undistribute);
    mxx::all2all(input.data(), first_part / _comm.size(), output.data(), _comm);
    BL_BENCH_END(undistribute, "a2a", first_part);

    BL_BENCH_START(undistribute);
    mxx::all2allv(input.data() + first_part, recv_counts, output.data() + first_part, send_counts, _comm);
    BL_BENCH_END(undistribute, "a2av", second_part);

    if (restore_order) {
      BL_BENCH_START(undistribute);
      // distribute (communication part)
      imxx::local::unpermute_inplace(output, i2o, 0, output.size());
      BL_BENCH_END(undistribute, "unpermute_inplace", output.size());

    }
    BL_BENCH_REPORT_MPI_NAMED(undistribute, "map_base:undistribute", _comm);

  }


  /**
   * @brief distribute, compute, send back.  one to one.  result matching input in order at then end.
   * @detail   this is the memory inefficient version
   *
   *
   */
  template <typename V, typename ToRank, typename Operation, typename SIZE = size_t,
      typename T = typename bliss::functional::function_traits<Operation, V>::return_type>
  void scatter_compute_gather(::std::vector<V>& input, ToRank const & to_rank,
                              Operation const & op,
                              ::std::vector<SIZE> & i2o,
                              ::std::vector<T>& output,
                              ::std::vector<V>& in_buffer, std::vector<T>& out_buffer,
                              ::mxx::comm const &_comm,
                              bool const & preserve_input = true) {
      BL_BENCH_INIT(scat_comp_gath);

      // speed over mem use.  mxx all2allv already has to double memory usage. same as stable distribute.
      BL_BENCH_COLLECTIVE_START(scat_comp_gath, "empty", _comm);
      bool empty = input.size() == 0;
      empty = mxx::all_of(empty);
      BL_BENCH_END(scat_comp_gath, "empty", input.size());

      if (empty) {
        BL_BENCH_REPORT_MPI_NAMED(scat_comp_gath, "map_base:scat_comp_gath", _comm);
        return;
      }

      // do assignment.
      BL_BENCH_START(scat_comp_gath);
      std::vector<SIZE> recv_counts(_comm.size(), 0);
      i2o.resize(input.size());
      BL_BENCH_END(scat_comp_gath, "alloc_map", input.size());

      // distribute
      BL_BENCH_START(scat_comp_gath);
      distribute(input, to_rank, recv_counts, i2o, in_buffer, _comm, false);
      BL_BENCH_END(scat_comp_gath, "distribute", in_buffer.size());

      // allocate out_buffer - output is same size as input
      BL_BENCH_START(scat_comp_gath);
      if (out_buffer.capacity() < (in_buffer.size())) out_buffer.clear();
      out_buffer.resize(in_buffer.size());
      BL_BENCH_END(scat_comp_gath, "alloc_outbuf", out_buffer.size());

      // process
      BL_BENCH_START(scat_comp_gath);
      op(in_buffer.begin(), in_buffer.end(), out_buffer.begin());
      BL_BENCH_END(scat_comp_gath, "compute", out_buffer.size());

      // allocate output - output is same size as input
      BL_BENCH_START(scat_comp_gath);
      if (output.capacity() < (input.size())) output.clear();
      output.resize(input.size());
      BL_BENCH_END(scat_comp_gath, "alloc_out", output.size());

      // distribute data back to source
      BL_BENCH_START(scat_comp_gath);
      undistribute(out_buffer, recv_counts, i2o, output, _comm, false);
      BL_BENCH_END(scat_comp_gath, "undistribute", output.size());


      // permute
      if (preserve_input) {
        BL_BENCH_START(scat_comp_gath);
        ::imxx::local::unpermute(input.begin(), input.end(), i2o.begin(), in_buffer.begin(), 0);
        in_buffer.swap(input);
        ::imxx::local::unpermute(output.begin(), output.end(), i2o.begin(), out_buffer.begin(), 0);
        out_buffer.swap(output);
        BL_BENCH_END(scat_comp_gath, "unpermute_inplace", output.size());
      }

      BL_BENCH_REPORT_MPI_NAMED(scat_comp_gath, "map_base:scat_comp_gath", _comm);
  }


  /**
   * @brief distribute, compute, send back.  one to one.  result matching input in order at then end.
   * @details  this is the memory efficient version.  this has to be incremental.
   *            this version uses a permute buffer. (in_buffer)
   */
  template <typename V, typename ToRank, typename Operation, typename SIZE = size_t,
      typename T = typename bliss::functional::function_traits<Operation, V>::return_type>
  void scatter_compute_gather_2part(::std::vector<V>& input, ToRank const & to_rank,
                              Operation const & op,
                              ::std::vector<SIZE> & i2o,
                              ::std::vector<T>& output,
                              ::std::vector<V>& in_buffer, std::vector<T>& out_buffer,
                              ::mxx::comm const &_comm,
                              bool const & preserve_input = true) {
      BL_BENCH_INIT(scat_comp_gath_2);

      // speed over mem use.  mxx all2allv already has to double memory usage. same as stable scat_comp_gath_2.
      BL_BENCH_COLLECTIVE_START(scat_comp_gath_2, "empty", _comm);
      bool empty = input.size() == 0;
      empty = mxx::all_of(empty);
      BL_BENCH_END(scat_comp_gath_2, "empty", input.size());

      if (empty) {
        BL_BENCH_REPORT_MPI_NAMED(scat_comp_gath_2, "map_base:scat_comp_gath_2", _comm);
        return;
      }

      // do assignment.
      BL_BENCH_START(scat_comp_gath_2);
      std::vector<SIZE> send_counts(_comm.size(), 0);
      std::vector<SIZE> recv_counts(_comm.size(), 0);
      i2o.resize(input.size());
      BL_BENCH_END(scat_comp_gath_2, "alloc_map", input.size());

      // first bucketing
      BL_BENCH_START(scat_comp_gath_2);
      imxx::local::assign_to_buckets(input, to_rank, _comm.size(), send_counts, i2o, 0, input.size());
      BL_BENCH_END(scat_comp_gath_2, "bucket", input.size());

      // then compute minimum block size.
      BL_BENCH_START(scat_comp_gath_2);
      SIZE min_bucket_size = *(::std::min_element(send_counts.begin(), send_counts.end()));
      min_bucket_size = ::mxx::allreduce(min_bucket_size, mxx::min<SIZE>(), _comm);
      SIZE first_part = _comm.size() * min_bucket_size;
      BL_BENCH_END(scat_comp_gath_2, "min_bucket_size", first_part);

      // compute the permutations from block size and processor mapping.  send_counts modified to the remainders.
      BL_BENCH_START(scat_comp_gath_2);
      ::imxx::local::bucket_to_block_permutation(min_bucket_size, 1UL, send_counts, i2o, 0, input.size());
      BL_BENCH_END(scat_comp_gath_2, "to_pos", input.size());

      // compute receive counts and total
      BL_BENCH_START(scat_comp_gath_2);
      recv_counts.resize(_comm.size());
      mxx::all2all(send_counts.data(), 1, recv_counts.data(), _comm);
      SIZE second_part = std::accumulate(recv_counts.begin(), recv_counts.end(), static_cast<size_t>(0));
      BL_BENCH_END(scat_comp_gath_2, "a2av_count", second_part);

      // allocate input buffer (for permuted data.)
      BL_BENCH_START(scat_comp_gath_2);
      if (in_buffer.capacity() < input.size()) in_buffer.clear();
      in_buffer.resize(input.size());
      BL_BENCH_END(scat_comp_gath_2, "alloc_inbuf", in_buffer.size());

      // permute
      BL_BENCH_START(scat_comp_gath_2);
      imxx::local::permute(input.begin(), input.end(), i2o.begin(), in_buffer.begin(), 0);
      in_buffer.swap(input);       // input is now permuted.
//      in_buffer.resize(second_part);
      BL_BENCH_END(scat_comp_gath_2, "permute", input.size());


      // allocate output - output is same size as input
      BL_BENCH_START(scat_comp_gath_2);
      if (output.capacity() < (input.size())) output.clear();
      output.resize(input.size());
      BL_BENCH_END(scat_comp_gath_2, "alloc_out", output.size());

      //== process first part.  communicate in place
      BL_BENCH_START(scat_comp_gath_2);
      block_all2all(input, min_bucket_size, in_buffer, 0, 0, _comm);
      BL_BENCH_END(scat_comp_gath_2, "a2a_inplace", first_part);

      // process
      BL_BENCH_START(scat_comp_gath_2);
      op(in_buffer.begin(), in_buffer.begin() + first_part, output.begin());
      BL_BENCH_END(scat_comp_gath_2, "compute1", first_part);

      // send the results back.  and reverse the input
      // undo a2a, so that result data matches.
      BL_BENCH_START(scat_comp_gath_2);
      block_all2all_inplace(output, min_bucket_size, 0, _comm);
      BL_BENCH_END(scat_comp_gath_2, "inverse_a2a_inplace", first_part);

      //======= process the second part
      // allocate output - output is same size as input
      BL_BENCH_START(scat_comp_gath_2);
      if (out_buffer.capacity() < second_part) out_buffer.clear();
      out_buffer.resize(second_part);
      in_buffer.resize(second_part);
      BL_BENCH_END(scat_comp_gath_2, "alloc_outbuf", out_buffer.size());

      // send second part.  reuse entire in_buffer
      BL_BENCH_START(scat_comp_gath_2);
	  mxx::all2allv(input.data() + first_part, send_counts,
                    in_buffer.data(), recv_counts, _comm);
      BL_BENCH_END(scat_comp_gath_2, "a2av", in_buffer.size());

      // process the second part.
      BL_BENCH_START(scat_comp_gath_2);
      op(in_buffer.begin(), in_buffer.begin() + second_part, out_buffer.begin());
      BL_BENCH_END(scat_comp_gath_2, "compute2", out_buffer.size());

      // send the results back
      BL_BENCH_START(scat_comp_gath_2);
	  mxx::all2allv(out_buffer.data(), recv_counts,
                    output.data() + first_part, send_counts, _comm);
      BL_BENCH_END(scat_comp_gath_2, "inverse_a2av", output.size());


      // permute
      if (preserve_input) {
        BL_BENCH_START(scat_comp_gath_2);
        // in_buffer was already allocated to be same size as input.
        ::imxx::local::unpermute(input.begin(), input.end(), i2o.begin(), in_buffer.begin(), 0);
        in_buffer.swap(input);
        // out_buffer is small, so should do this inplace.
        ::imxx::local::unpermute_inplace(output, i2o, 0, output.size());
        BL_BENCH_END(scat_comp_gath_2, "unpermute_inplace", output.size());
      }

      BL_BENCH_REPORT_MPI_NAMED(scat_comp_gath_2, "map_base:scat_comp_gath_2", _comm);
  }

  /**
   * @brief distribute, compute, send back.  one to one.  result matching input in order at then end.
   * @details  this is the memory efficient version.  this has to be incremental.
   *            this version uses a permute buffer. (in_buffer)
   *
   *            low mem version.
   *
   *            the second half (a2av) require the space it requires.  for in_buffer and out_buffer.
   *            the first half can be reduced in space usage by incrementally populate the in_buffer,
   *            	compute the output, then populate the output.  then 1 all2all inplace for the output
   *            	at the requester, then unpermute somehow.
   *
      // compute the size of the buffers to use.  set a maximum limit.
       * let the final block_bucket_size be bbs. then the first_part x = bbs * blocks * p., p = comm_size.
       * let the sum of UNEVEN parts in orig input be y1, and after a2av be y2.  the minimum space needed here is y1+y2 for the a2av part.
       *   the uneven part is defined by the difference between a bucket and the global min bucket.
       *
       * APPROACH:  use either y1+y2 or 1/2 min_bucket, whichever is larger, as buffer
       * 	if y1+y2 > 1/2 (input + y2), use traditional
       * 	else
       * 		use larger of y1+y2, and 1/2 min bucket.
       *
       *
       * however, bbs may not perfectly divide the min_bucket..  let the remainder be r.  0 <= r < bbs
       * the UNEVEN part then needs to be appended with r in each bucket.  Y1 = y1 + r * p, and Y2 = y2 + r * p.  post communication, the remainder size remains same.
       * minimum space needed is revised to Y1 + Y2 = y1 + y2 + 2 rp.
       *
       * if y1 + y2 > input_size, then we should just the normal scatter_compute_gather, or the 2 part version, since buffer requirement is high.
       *
       * case 1
       * bbs * p > Y1 + Y2.  each block still needs to got through a2a, so makes no sense to have more than 1 blocks in the buffer and we should reduce the number of iterations
       * to have some form of savings, we should use bbs * p < input.  then input > Y2 + Y2.  if this does not hold, then use full memory version.
       *
       * for all processors, the following must hold to have buffer savings.
       * y1 + rp + bbs * blocks * p > bbs * p > y1 + y2 + 2rp.   r + bbs * blocks = min_bucket_size.
       * 	p * min_bucket_size > y2 + 2rp,  so r < min_bucket_size - y2 / 2p,
       * 	mod(min_bucket_size, bbs) = r < bbs.
       * so let bbs = min_bucket_size - y2 / 2p.  (bbs can be at most this quantity, else r may exceed bound)
       * 	now this can be pretty big, or really small.
       *
       * what is the lower bound of bbs?
       *
       *
       *
       *
       *
       * for there to be space savings, y1+y2 is the minimum buffer size and the maximum buffer size should not exceed 1/2 of input
       * buffer size should be std::max(y1+y2, comm_size * block_bucket_size)   1 block...
       * where block_bucket_size = std::min(min_bucket_size, 1/2 max_bucket_size)    // at most min_bucket_size, at least 1/2 max_bucket_size
       * this will affect y1 and y2, of course.
       * 	y1 has size input % (comm_size * block_bucket_size), with max value of (comm_size * block_bucket_size - 1)
       * 	y2 has max comm_size * (comm_size  * block_bucket_size - 1 ) - 1
       *
       * note at the end, we'd still need to permute either the input or the output to make them consistent with each other.
      // max of "second_part" is the minimum required space for the a2av portion.  let this be y.
      //    reducing this requires O(p) type communication rather than O(log(p))
      // to effect some space savings, we should use max(1/2 input_size, y).
      //	we can reduce the 1/2 input size, at the expense of more iterations, each requiring a complete scan of input during permuting.
      //	if y > 1/2 input, then on one processor the min bucket is less than 1/2 input,  we'd be using this min bucket anyways.
      //	if y < 1/2 input, then we can potentially do more in a2a phase, but that would use more mem.
      // so configure the buffer size to be max(1/2 input_size, y) = Y.
      //    for bucket size, this means min(min_bucket, (Y + comm_size - 1) / comm_size))
      // a simpler way to look at bucket size is to do min(max_bucket/2, min_bucket), and the buffer size is max(1/2 input_size, y)
       *
       * the buffer is set to largest of first part, local second part, remote second part.
       *
       * break up the input into smaller chunks imply iterative process.  each iteration requires full input scan.  if we do logarithmic number of iterations, we can't scale very large.
       * instead we break up into 2 pieces only, (or maybe break it up to 3 or 4), a low number, and accept the overhead incurred.
   */
  template <typename V, typename ToRank, typename Operation, typename SIZE = size_t,
      typename T = typename bliss::functional::function_traits<Operation, V>::return_type>
  void scatter_compute_gather_lowmem(::std::vector<V>& input, ToRank const & to_rank,
                              Operation const & op,
                              ::std::vector<SIZE> & i2o,
                              ::std::vector<T>& output,
                              ::std::vector<V>& in_buffer, std::vector<T>& out_buffer,
                              ::mxx::comm const &_comm,
                              bool const & preserve_input = true) {
      BL_BENCH_INIT(scat_comp_gath_lm);

      BL_BENCH_COLLECTIVE_START(scat_comp_gath_lm, "empty", _comm);
      bool empty = input.size() == 0;
      empty = mxx::all_of(empty);
      BL_BENCH_END(scat_comp_gath_lm, "empty", input.size());

      if (empty) {
        BL_BENCH_REPORT_MPI_NAMED(scat_comp_gath_lm, "map_base:scat_comp_gath_lm", _comm);
        return;
      }



      // do assignment.
      BL_BENCH_START(scat_comp_gath_lm);
      std::vector<SIZE> send_counts(_comm.size(), 0);
      std::vector<SIZE> recv_counts(_comm.size(), 0);
      i2o.resize(input.size());
      BL_BENCH_END(scat_comp_gath_lm, "alloc_map", input.size());

      // first bucketing
      BL_BENCH_START(scat_comp_gath_lm);
      imxx::local::assign_to_buckets(input, to_rank, _comm.size(), send_counts, i2o, 0, input.size());
      BL_BENCH_END(scat_comp_gath_lm, "bucket", input.size());

      // then compute minimum block size.
      BL_BENCH_START(scat_comp_gath_lm);
      SIZE min_bucket_size = *(::std::min_element(send_counts.begin(), send_counts.end()));
      min_bucket_size = ::mxx::allreduce(min_bucket_size, mxx::min<SIZE>(), _comm);
      SIZE block_bucket_size = min_bucket_size / 2;  // block_bucket_size is at least 1/2 as large as the largest bucket.
      SIZE block_size = _comm.size() * block_bucket_size;
      SIZE first_part = 2 * block_size;   // this is at least 1/2 of the largest input.
      SIZE second_part_local = input.size() - first_part;
      recv_counts.resize(_comm.size());
      ::mxx::all2all(send_counts.data(), 1, recv_counts.data(), _comm);
      SIZE second_part_remote = std::accumulate(recv_counts.begin(), recv_counts.end(), static_cast<size_t>(0)) - first_part;   // second part size.

      SIZE second = second_part_local + second_part_remote;
      SIZE input_plus = input.size() + second_part_remote;
      bool traditional = (second > (input_plus / 2));
      traditional = mxx::any_of(traditional, _comm);
      BL_BENCH_END(scat_comp_gath_lm, "a2av_count", first_part);

      if (traditional) {
    	  BL_BENCH_START(scat_comp_gath_lm);

    	  scatter_compute_gather_2part(input, to_rank, op, i2o, output, in_buffer, out_buffer, _comm, preserve_input);
          BL_BENCH_END(scat_comp_gath_lm, "switch_to_trad", output.size());


        BL_BENCH_REPORT_MPI_NAMED(scat_comp_gath_lm, "map_base:scat_comp_gath_lm", _comm);
        return;
      }

      // compute the permutations from block size and processor mapping.  send_counts modified to the remainders.
      BL_BENCH_START(scat_comp_gath_lm);
      ::imxx::local::bucket_to_block_permutation(block_bucket_size, 2UL, send_counts, i2o, 0, input.size());
      BL_BENCH_END(scat_comp_gath_lm, "to_pos", input.size());

      // allocate input buffer (for permuted data.)
      BL_BENCH_START(scat_comp_gath_lm);
      SIZE buffer_size = std::max(block_size, second_part_local + second_part_remote);
      if (in_buffer.capacity() < buffer_size) in_buffer.clear();
      in_buffer.resize(buffer_size);
      BL_BENCH_END(scat_comp_gath_lm, "alloc_inbuf", in_buffer.size());

      // allocate output - output is same size as input
      BL_BENCH_START(scat_comp_gath_lm);
      if (output.capacity() < (input.size())) output.clear();
      output.resize(input.size());
      BL_BENCH_END(scat_comp_gath_lm, "alloc_out", output.size());

      //== process first part - 2 iterations..  communicate in place
      for (size_t i = 0; i < 2; ++i) {
		  // permute
		  BL_BENCH_START(scat_comp_gath_lm);
		  ::imxx::local::permute_for_output_range(input.begin(), input.end(), i2o.begin(), in_buffer.begin(), in_buffer.begin() + block_size, i * block_size);
		  BL_BENCH_END(scat_comp_gath_lm, "permute_block", block_size);

		  BL_BENCH_START(scat_comp_gath_lm);
		  block_all2all_inplace(in_buffer, block_bucket_size, 0, _comm);
		  BL_BENCH_END(scat_comp_gath_lm, "a2a_inplace", block_size);

		  // process
		  BL_BENCH_START(scat_comp_gath_lm);
		  op(in_buffer.begin(), in_buffer.begin() + block_size, output.begin() + i * block_size);
		  BL_BENCH_END(scat_comp_gath_lm, "compute1", block_size);

		  // send the results back.  and reverse the input
		  // undo a2a, so that result data matches.
		  BL_BENCH_START(scat_comp_gath_lm);
		  block_all2all_inplace(output, block_bucket_size, i * block_size, _comm);
		  BL_BENCH_END(scat_comp_gath_lm, "inverse_a2a_inplace", block_size);

      }

      //======= process the second part
      // allocate output - output is same size as input
      BL_BENCH_START(scat_comp_gath_lm);
      if (out_buffer.capacity() < second_part_remote) out_buffer.clear();
      out_buffer.resize(second_part_remote);
      BL_BENCH_END(scat_comp_gath_lm, "alloc_outbuf", out_buffer.size());

	  // permute
	  BL_BENCH_START(scat_comp_gath_lm);
	  ::imxx::local::permute_for_output_range(input.begin(), input.end(), i2o.begin(),
			  in_buffer.begin(), in_buffer.begin() + second_part_local, first_part);
	  BL_BENCH_END(scat_comp_gath_lm, "permute_block", second_part_local);

      // send second part.  reuse entire in_buffer
      BL_BENCH_START(scat_comp_gath_lm);
      ::mxx::all2all(send_counts.data(), 1, recv_counts.data(), _comm);
	  mxx::all2allv(in_buffer.data(), send_counts,
                    in_buffer.data() + second_part_local, recv_counts, _comm);
      BL_BENCH_END(scat_comp_gath_lm, "a2av", second_part_local);

      // process the second part.
      BL_BENCH_START(scat_comp_gath_lm);
      op(in_buffer.begin() + second_part_local, in_buffer.begin() + second_part_local + second_part_remote, out_buffer.begin());
      BL_BENCH_END(scat_comp_gath_lm, "compute2", second_part_remote);

      // send the results back
      BL_BENCH_START(scat_comp_gath_lm);
	  mxx::all2allv(out_buffer.data(), recv_counts,
                    output.data() + first_part, send_counts, _comm);
      BL_BENCH_END(scat_comp_gath_lm, "inverse_a2av", second_part_remote);

      // permute
      if (preserve_input) {
        // in_buffer was already allocated to be same size as input.
        ::imxx::local::unpermute_inplace(input, i2o, 0, input.size());
        // out_buffer is small, so should do this inplace.
        ::imxx::local::unpermute_inplace(output, i2o, 0, output.size());
        BL_BENCH_END(scat_comp_gath_lm, "unpermute_inplace", output.size());
      }

      BL_BENCH_REPORT_MPI_NAMED(scat_comp_gath_lm, "map_base:scat_comp_gath_lm", _comm);
  }

  // TODO: non-one-to-one version.

  /**
   * @brief distribute, compute, send back.  one to one.  result matching input in order at then end.
   * @detail   this is the memory inefficient version
   *
   *
   */
  template <typename V, typename ToRank, typename Operation, typename SIZE = size_t,
      typename T = typename bliss::functional::function_traits<Operation, V>::return_type>
  void scatter_compute_gather_v(::std::vector<V>& input, ToRank const & to_rank,
                              Operation const & op,
                              ::std::vector<SIZE> & i2o,
                              ::std::vector<T>& output,
                              ::std::vector<V>& in_buffer, std::vector<T>& out_buffer,
                              ::mxx::comm const &_comm,
                              bool const & preserve_input = true) {
      BL_BENCH_INIT(scat_comp_gath_v);

      // speed over mem use.  mxx all2allv already has to double memory usage. same as stable distribute.
      BL_BENCH_COLLECTIVE_START(scat_comp_gath_v, "empty", _comm);
      bool empty = input.size() == 0;
      empty = mxx::all_of(empty);
      BL_BENCH_END(scat_comp_gath_v, "empty", input.size());

      if (empty) {
        BL_BENCH_REPORT_MPI_NAMED(scat_comp_gath_v, "map_base:scat_comp_gath_v", _comm);
        return;
      }

      // do assignment.
      BL_BENCH_START(scat_comp_gath_v);
      std::vector<SIZE> recv_counts(_comm.size(), 0);
      i2o.resize(input.size());
      BL_BENCH_END(scat_comp_gath_v, "alloc_map", input.size());

      // distribute
      BL_BENCH_START(scat_comp_gath_v);
      distribute(input, to_rank, recv_counts, i2o, in_buffer, _comm, false);
      BL_BENCH_END(scat_comp_gath_v, "distribute", in_buffer.size());

      // allocate out_buffer - output is same size as input
      BL_BENCH_START(scat_comp_gath_v);
      if (out_buffer.capacity() < (in_buffer.size())) out_buffer.clear();
      out_buffer.reserve(in_buffer.size());
      ::fsc::back_emplace_iterator<std::vector<T> > emplacer(out_buffer);
      BL_BENCH_END(scat_comp_gath_v, "alloc_outbuf", out_buffer.size());

      // process
      BL_BENCH_START(scat_comp_gath_v);
      size_t s;
      auto it = in_buffer.begin();
      for (size_t i = 0; i < _comm.size(); ++i) {
    	  s = recv_counts[i];

          recv_counts[i] = op(it, it + s, emplacer);
          std::advance(it, s);
      }
      BL_BENCH_END(scat_comp_gath_v, "compute", out_buffer.size());

      // allocate output - output is same size as input
      BL_BENCH_START(scat_comp_gath_v);
      std::vector<SIZE> send_counts(_comm.size(), 0);
      ::mxx::all2all(recv_counts.data(), 1, send_counts.data(), _comm);
      size_t total = ::std::accumulate(send_counts.begin(), send_counts.end(), static_cast<size_t>(0));
      if (output.capacity() < total) output.clear();
      output.resize(total);
      BL_BENCH_END(scat_comp_gath_v, "alloc_out", output.size());

      // distribute data back to source
      BL_BENCH_START(scat_comp_gath_v);
      undistribute(out_buffer, recv_counts, i2o, output, _comm, false);
      BL_BENCH_END(scat_comp_gath_v, "undistribute", output.size());

      // permute
      if (preserve_input) {
        BL_BENCH_START(scat_comp_gath_v);
        ::imxx::local::unpermute(input.begin(), input.end(), i2o.begin(), in_buffer.begin(), 0);
        in_buffer.swap(input);
        ::imxx::local::unpermute(output.begin(), output.end(), i2o.begin(), out_buffer.begin(), 0);
        out_buffer.swap(output);
        BL_BENCH_END(scat_comp_gath_v, "unpermute_inplace", output.size());
      }

      BL_BENCH_REPORT_MPI_NAMED(scat_comp_gath_v, "map_base:scat_comp_gath_v", _comm);
  }

  //TODO:
//
//  /**
//   * @brief distribute, compute, send back.  one to one.  result matching input in order at then end.
//   * @details  this is the memory efficient version.  this has to be incremental.
//   *            this version uses a permute buffer. (in_buffer)
//   *
//   *            low mem version.
//   *
//   *            the second half (a2av) require the space it requires.  for in_buffer and out_buffer.
//   *            the first half can be reduced in space usage by incrementally populate the in_buffer,
//   *            	compute the output, then populate the output.  then 1 all2all inplace for the output
//   *            	at the requester, then unpermute somehow.
//   *
//      // compute the size of the buffers to use.  set a maximum limit.
//       * let the final block_bucket_size be bbs. then the first_part x = bbs * blocks * p., p = comm_size.
//       * let the sum of UNEVEN parts in orig input be y1, and after a2av be y2.  the minimum space needed here is y1+y2 for the a2av part.
//       *   the uneven part is defined by the difference between a bucket and the global min bucket.
//       *
//       * APPROACH:  use either y1+y2 or 1/2 min_bucket, whichever is larger, as buffer
//       * 	if y1+y2 > 1/2 (input + y2), use traditional
//       * 	else
//       * 		use larger of y1+y2, and 1/2 min bucket.
//       *
//       *
//       * however, bbs may not perfectly divide the min_bucket..  let the remainder be r.  0 <= r < bbs
//       * the UNEVEN part then needs to be appended with r in each bucket.  Y1 = y1 + r * p, and Y2 = y2 + r * p.  post communication, the remainder size remains same.
//       * minimum space needed is revised to Y1 + Y2 = y1 + y2 + 2 rp.
//       *
//       * if y1 + y2 > input_size, then we should just the normal scatter_compute_gather, or the 2 part version, since buffer requirement is high.
//       *
//       * case 1
//       * bbs * p > Y1 + Y2.  each block still needs to got through a2a, so makes no sense to have more than 1 blocks in the buffer and we should reduce the number of iterations
//       * to have some form of savings, we should use bbs * p < input.  then input > Y2 + Y2.  if this does not hold, then use full memory version.
//       *
//       * for all processors, the following must hold to have buffer savings.
//       * y1 + rp + bbs * blocks * p > bbs * p > y1 + y2 + 2rp.   r + bbs * blocks = min_bucket_size.
//       * 	p * min_bucket_size > y2 + 2rp,  so r < min_bucket_size - y2 / 2p,
//       * 	mod(min_bucket_size, bbs) = r < bbs.
//       * so let bbs = min_bucket_size - y2 / 2p.  (bbs can be at most this quantity, else r may exceed bound)
//       * 	now this can be pretty big, or really small.
//       *
//       * what is the lower bound of bbs?
//       *
//       *
//       *
//       *
//       *
//       * for there to be space savings, y1+y2 is the minimum buffer size and the maximum buffer size should not exceed 1/2 of input
//       * buffer size should be std::max(y1+y2, comm_size * block_bucket_size)   1 block...
//       * where block_bucket_size = std::min(min_bucket_size, 1/2 max_bucket_size)    // at most min_bucket_size, at least 1/2 max_bucket_size
//       * this will affect y1 and y2, of course.
//       * 	y1 has size input % (comm_size * block_bucket_size), with max value of (comm_size * block_bucket_size - 1)
//       * 	y2 has max comm_size * (comm_size  * block_bucket_size - 1 ) - 1
//       *
//       * note at the end, we'd still need to permute either the input or the output to make them consistent with each other.
//      // max of "second_part" is the minimum required space for the a2av portion.  let this be y.
//      //    reducing this requires O(p) type communication rather than O(log(p))
//      // to effect some space savings, we should use max(1/2 input_size, y).
//      //	we can reduce the 1/2 input size, at the expense of more iterations, each requiring a complete scan of input during permuting.
//      //	if y > 1/2 input, then on one processor the min bucket is less than 1/2 input,  we'd be using this min bucket anyways.
//      //	if y < 1/2 input, then we can potentially do more in a2a phase, but that would use more mem.
//      // so configure the buffer size to be max(1/2 input_size, y) = Y.
//      //    for bucket size, this means min(min_bucket, (Y + comm_size - 1) / comm_size))
//      // a simpler way to look at bucket size is to do min(max_bucket/2, min_bucket), and the buffer size is max(1/2 input_size, y)
//       *
//       * the buffer is set to largest of first part, local second part, remote second part.
//       *
//       * break up the input into smaller chunks imply iterative process.  each iteration requires full input scan.  if we do logarithmic number of iterations, we can't scale very large.
//       * instead we break up into 2 pieces only, (or maybe break it up to 3 or 4), a low number, and accept the overhead incurred.
//   */
//  template <typename V, typename ToRank, typename Operation, typename SIZE = size_t,
//      typename T = typename bliss::functional::function_traits<Operation, V>::return_type>
//  void scatter_compute_gather_v_lowmem(::std::vector<V>& input, ToRank const & to_rank,
//                              Operation const & op,
//                              ::std::vector<SIZE> & i2o,
//                              ::std::vector<T>& output,
//                              ::std::vector<V>& in_buffer, std::vector<T>& out_buffer,
//                              ::mxx::comm const &_comm,
//                              bool const & preserve_input = true) {
//      BL_BENCH_INIT(scat_comp_gath_v_lm);
//
//      bool empty = input.size() == 0;
//      empty = mxx::all_of(empty);
//      if (empty) {
//        BL_BENCH_REPORT_MPI_NAMED(scat_comp_gath_v_lm, "map_base:scat_comp_gath_v_lm", _comm);
//        return;
//      }
//
//
//
//      // do assignment.
//      BL_BENCH_START(scat_comp_gath_v_lm);
//      std::vector<SIZE> send_counts(_comm.size(), 0);
//      std::vector<SIZE> recv_counts(_comm.size(), 0);
//      i2o.resize(input.size());
//      BL_BENCH_END(scat_comp_gath_v_lm, "alloc_map", input.size());
//
//      // first bucketing
//      BL_BENCH_START(scat_comp_gath_v_lm);
//      imxx::local::assign_to_buckets(input, to_rank, _comm.size(), send_counts, i2o, 0, input.size());
//      BL_BENCH_END(scat_comp_gath_v_lm, "bucket", input.size());
//
//      // then compute minimum block size.
//      BL_BENCH_COLLECTIVE_START(scat_comp_gath_v_lm, "a2av_count", _comm);
//      SIZE min_bucket_size = *(::std::min_element(send_counts.begin(), send_counts.end()));
//      min_bucket_size = ::mxx::allreduce(min_bucket_size, mxx::min<SIZE>(), _comm);
//      SIZE block_bucket_size = min_bucket_size / 2;  // block_bucket_size is at least 1/2 as large as the largest bucket.
//      SIZE block_size = _comm.size() * block_bucket_size;
//      SIZE first_part = 2 * block_size;   // this is at least 1/2 of the largest input.
//      SIZE second_part_local = input.size() - first_part;
//      recv_counts.resize(_comm.size());
//      ::mxx::all2all(send_counts.data(), 1, recv_counts.data(), _comm);
//      SIZE second_part_remote = std::accumulate(recv_counts.begin(), recv_counts.end(), static_cast<size_t>(0)) - first_part;   // second part size.
//
//      SIZE second = second_part_local + second_part_remote;
//      SIZE input_plus = input.size() + second_part_remote;
//      bool traditional = (second > (input_plus / 2));
//      traditional = mxx::any_of(traditional, _comm);
//      BL_BENCH_END(scat_comp_gath_v_lm, "a2av_count", first_part);
//
//      if (traditional) {
//    	  BL_BENCH_START(scat_comp_gath_v_lm);
//
//    	  scatter_compute_gather_2part(input, to_rank, op, i2o, output, in_buffer, out_buffer, _comm, preserve_input);
//          BL_BENCH_END(scat_comp_gath_v_lm, "switch_to_trad", output.size());
//
//
//        BL_BENCH_REPORT_MPI_NAMED(scat_comp_gath_v_lm, "map_base:scat_comp_gath_v_lm", _comm);
//        return;
//      }
//
//      // compute the permutations from block size and processor mapping.  send_counts modified to the remainders.
//      BL_BENCH_START(scat_comp_gath_v_lm);
//      ::imxx::local::bucket_to_block_permutation(block_bucket_size, 2UL, send_counts, i2o, 0, input.size());
//      BL_BENCH_END(scat_comp_gath_v_lm, "to_pos", input.size());
//
//      // allocate input buffer (for permuted data.)
//      BL_BENCH_START(scat_comp_gath_v_lm);
//      SIZE buffer_size = std::max(block_size, second_part_local + second_part_remote);
//      if (in_buffer.capacity() < buffer_size) in_buffer.clear();
//      in_buffer.resize(buffer_size);
//      BL_BENCH_END(scat_comp_gath_v_lm, "alloc_inbuf", in_buffer.size());
//
//      // allocate output - output is same size as input
//      BL_BENCH_START(scat_comp_gath_v_lm);
//      if (output.capacity() < (input.size())) output.clear();
//      output.resize(input.size());
//      BL_BENCH_END(scat_comp_gath_v_lm, "alloc_out", output.size());
//
//      //== process first part - 2 iterations..  communicate in place
//      for (size_t i = 0; i < 2; ++i) {
//		  // permute
//		  BL_BENCH_START(scat_comp_gath_v_lm);
//		  ::imxx::local::permute_for_output_range(input.begin(), input.end(), i2o.begin(), in_buffer.begin(), in_buffer.begin() + block_size, i * block_size);
//		  BL_BENCH_END(scat_comp_gath_v_lm, "permute_block", block_size);
//
//		  BL_BENCH_COLLECTIVE_START(scat_comp_gath_v_lm, "a2a_inplace", _comm);
//		  block_all2all_inplace(in_buffer, block_bucket_size, 0, _comm);
//		  BL_BENCH_END(scat_comp_gath_v_lm, "a2a_inplace", block_size);
//
//		  // process
//		  BL_BENCH_START(scat_comp_gath_v_lm);
//		  op(in_buffer.begin(), in_buffer.begin() + block_size, output.begin() + i * block_size);
//		  BL_BENCH_END(scat_comp_gath_v_lm, "compute1", block_size);
//
//		  // send the results back.  and reverse the input
//		  // undo a2a, so that result data matches.
//		  BL_BENCH_COLLECTIVE_START(scat_comp_gath_v_lm, "inverse_a2a_inplace", _comm);
//		  block_all2all_inplace(output, block_bucket_size, i * block_size, _comm);
//		  BL_BENCH_END(scat_comp_gath_v_lm, "inverse_a2a_inplace", block_size);
//
//      }
//
//
//      //======= process the second part
//      // allocate output - output is same size as input
//      BL_BENCH_START(scat_comp_gath_v_lm);
//      if (out_buffer.capacity() < second_part_remote) out_buffer.clear();
//      out_buffer.resize(second_part_remote);
//      BL_BENCH_END(scat_comp_gath_v_lm, "alloc_outbuf", out_buffer.size());
//
//	  // permute
//	  BL_BENCH_START(scat_comp_gath_v_lm);
//	  ::imxx::local::permute_for_output_range(input.begin(), input.end(), i2o.begin(),
//			  in_buffer.begin(), in_buffer.begin() + second_part_local, first_part);
//	  BL_BENCH_END(scat_comp_gath_v_lm, "permute_block", second_part_local);
//
//      // send second part.  reuse entire in_buffer
//      BL_BENCH_COLLECTIVE_START(scat_comp_gath_v_lm, "a2av", _comm);
//      ::mxx::all2all(send_counts.data(), 1, recv_counts.data(), _comm);
//	  mxx::all2allv(in_buffer.data(), send_counts,
//                    in_buffer.data() + second_part_local, recv_counts, _comm);
//      BL_BENCH_END(scat_comp_gath_v_lm, "a2av", second_part_local);
//
//      // process the second part.
//      BL_BENCH_START(scat_comp_gath_v_lm);
//      op(in_buffer.begin() + second_part_local, in_buffer.begin() + second_part_local + second_part_remote, out_buffer.begin());
//      BL_BENCH_END(scat_comp_gath_v_lm, "compute2", second_part_remote);
//
//      // send the results back
//      BL_BENCH_COLLECTIVE_START(scat_comp_gath_v_lm, "inverse_a2av", _comm);
//	  mxx::all2allv(out_buffer.data(), recv_counts,
//                    output.data() + first_part, send_counts, _comm);
//      BL_BENCH_END(scat_comp_gath_v_lm, "inverse_a2av", second_part_remote);
//
//      // permute
//      if (preserve_input) {
//        // in_buffer was already allocated to be same size as input.
//        ::imxx::local::unpermute_inplace(input, i2o, 0, input.size());
//        // out_buffer is small, so should do this inplace.
//        ::imxx::local::unpermute_inplace(output, i2o, 0, output.size());
//        BL_BENCH_END(scat_comp_gath_v_lm, "unpermute_inplace", output.size());
//      }
//
//      BL_BENCH_REPORT_MPI_NAMED(scat_comp_gath_v_lm, "map_base:scat_comp_gath_v_lm", _comm);
//  }
//

//== TODO: transform before or after communication


//== communicate, process (not one to one), return communication  - order does not matter


// TODO mpi3 versions?

} // namespace imxx
