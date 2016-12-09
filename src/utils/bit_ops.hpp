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
 * @file    bit_ops.hpp
 * @ingroup
 * @author  tpan
 * @brief   collection of bit operations.
 */
#ifndef SRC_UTILS_BIT_OPS_HPP_
#define SRC_UTILS_BIT_OPS_HPP_

#include <array>

namespace bliss {

  namespace utils {

    namespace bit_ops {


      template <typename DUMMY = void>
      struct pop_cnt {
          static constexpr std::array<uint8_t, 16> LUT =
          {{
            0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4
          }};

          template <typename T>
          uint8_t operator()(T const & x) const {
            uint8_t cnt = 0;
            T v = x;
            for (; v; v >>= 4) {
              // while v is not 0
              cnt += LUT[v & static_cast<T>(0xF)];  // get lowest 4 bits, then look up
            }
            return cnt;
          }
      };
      template <typename DUMMY>
      constexpr std::array<uint8_t, 16> pop_cnt<DUMMY>::LUT;


    } // namespace bit_ops


  } // namespace utils

} // namespace bliss








#endif /* SRC_UTILS_BIT_OPS_HPP_ */
