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
 * @file    file_utils.hpp
 * @ingroup
 * @author  tpan
 * @brief
 * @details
 *
 */
#ifndef SRC_UTILS_FILE_UTILS_HPP_
#define SRC_UTILS_FILE_UTILS_HPP_

#include <string>

namespace bliss {
  namespace utils {

    namespace file {

      std::string get_file_extension(std::string const & filename) {
        // find the last part.
        size_t pos = filename.find_last_of('.');
        if (pos == std::string::npos)  return std::string();
        else
          return filename.substr(pos + 1);  // from next char to end.
      }

      struct NotEOL {
        template <typename CharType>
        bool operator()(CharType const & x) {
        return (x != '\n') && (x != '\r' );
        }

        template <typename CharType, typename MDType>
        bool operator()(std::pair<CharType, MDType> const & x) {
          return (x.first != '\n') && (x.first != '\r');
        }
      };

    }

  }
}

#endif /* SRC_UTILS_FILE_UTILS_HPP_ */
