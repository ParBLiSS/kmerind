/**
 * @file		concurrent.hpp
 * @ingroup
 * @author	tpan
 * @brief
 * @details
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef CONCURRENT_HPP_
#define CONCURRENT_HPP_

namespace bliss {
  namespace concurrent {

    enum ThreadSafety : uint8_t {
      THREAD_UNSAFE = 0,
      THREAD_SAFE = 1
    };

  }
}



#endif /* CONCURRENT_HPP_ */
