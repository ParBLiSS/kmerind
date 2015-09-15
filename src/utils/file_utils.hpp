/**
 * @file    file_utils.hpp
 * @ingroup
 * @author  tpan
 * @brief
 * @details
 *
 * Copyright (c) 2015 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
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

    }
  }
}

#endif /* SRC_UTILS_FILE_UTILS_HPP_ */
