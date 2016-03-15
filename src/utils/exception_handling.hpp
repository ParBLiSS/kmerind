/**
 * @file    exception_handling.hpp
 * @ingroup
 * @author  tpan
 * @brief
 * @details
 *
 * Copyright (c) 2016 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef SRC_UTILS_EXCEPTION_HANDLING_HPP_
#define SRC_UTILS_EXCEPTION_HANDLING_HPP_

#include <string>
#include <iostream>  // for cerr

namespace bliss {

  namespace utils {

    template <typename EXCEPTION>
    inline EXCEPTION make_exception(::std::string const & msg) {
      ::std::cerr << msg << ::std::endl;
      return EXCEPTION(msg);
    }

  }
}

#endif /* SRC_UTILS_EXCEPTION_HANDLING_HPP_ */
