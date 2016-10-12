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
 * @file    algorithm.hpp
 * @ingroup iterators
 * @author  Tony Pan <tpan7@gatech.edu>
 * @brief   set of algorithms that are used with iterators.
 *
 */

#ifndef ALGORITHM_HPP
#define ALGORITHM_HPP

#include <iterator>
#include <algorithm>

namespace bliss
{

namespace iterator
{

namespace algorithm
{

/**
 * @brief conditional transform - predicate determines which transform to apply and which output iterator to populate.
 * @return  number marked as true.
 */
template <typename InputIterator, typename Predicate,
	typename TrueTransform, typename FalseTransform,
	typename TrueOutputIterator, typename FalseOutputIterator>
std::pair<TrueOutputIterator, FalseOutputIterator>
conditional_transform(InputIterator start, InputIterator end,
		TrueOutputIterator true_output, TrueTransform const & t_trans,
		FalseOutputIterator false_output, FalseTransform const & f_trans,
		Predicate const & pred) {

	for (; start != end; ++start) {
		auto val = *start;
		if (pred(val)) {
			*true_output = t_trans(val);
			++true_output;
		} else {
			*false_output = f_trans(val);
			++false_output;
		}
	}
	return ::std::make_pair(true_output, false_output);
}


/**
 * @brief conditional transform - transform is applied only if predicate returns true, else entry discarded..
 * @note	predicate is applied to input, not output.
 * @return  output iterator one past the last inserted position
 */
template <typename InputIterator, typename OutputIterator,
	typename Transform, typename Predicate>
OutputIterator
transform_if(InputIterator start, InputIterator end,
		OutputIterator output, Transform const & trans,
		Predicate const & pred) {

	for (; start != end; ++start) {
		if (pred(*start)) {
			*output = trans(*start);
			++output;
		}
	}
	return output;
}


}  // namespace algorithm
} // iterator
} // bliss
#endif /* ALGORITHM_HPP */
