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
 * @file    transform_utils.hpp
 * @ingroup bliss::utils
 * @author  tpan
 * @brief   functors related to filter operaitons.
 *
 */
#ifndef TRANSFORM_UTILS_HPP_
#define TRANSFORM_UTILS_HPP_

namespace bliss {

  namespace transform
  {

  // identity transform
  template <typename Key>
  struct identity {
      inline Key operator()(Key const & x) const {
        return std::forward<const Key>(x);
      }
      template <typename VAL>
      inline ::std::pair<Key, VAL> operator()(std::pair<Key, VAL> const & x) const {
        return std::forward<const std::pair<Key, VAL> >(x);
      }
      template <typename VAL>
      inline ::std::pair<const Key, VAL> operator()(std::pair<const Key, VAL> const & x) const {
        return std::forward<const std::pair<const Key, VAL> >(x);
      }

  };

  } // namespace filter

} // namespace bliss

#endif /* TRANSFORM_UTILS_HPP_ */
