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
 * @file    filter_utils.hpp
 * @ingroup bliss::utils
 * @author  tpan
 * @brief   functors related to filter operaitons.
 *
 */
#ifndef FILTER_UTILS_HPP_
#define FILTER_UTILS_HPP_

namespace bliss {

  namespace filter
  {

	  struct TruePredicate {
		  template <typename T>
		  bool operator()(T const & x) const { return true; }
		  template <typename Iter>
		  bool operator()(Iter b, Iter e) const { return true; }
	  };

  } // namespace filter

} // namespace bliss

#endif /* FILTER_UTILS_HPP_ */
