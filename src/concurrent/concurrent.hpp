/**
 * @file		concurrent.hpp
 * @ingroup bliss::concurrent
 * @author	tpan
 * @brief   this file contains predefined constants for concurrency and thread safety
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

    /// Type for indicating Thread Safety.
    typedef bool ThreadSafety;

    /// Constant indicating Thread Safe
    constexpr bool THREAD_SAFE = true;

    /// Constant indicating Thread Unsafe
    constexpr bool THREAD_UNSAFE = false;


    // types of threading locks.
    enum class LockType { NONE = 0, MUTEX = 1, SPINLOCK = 2, LOCK_FREE = 4 };

  }
}



#endif /* CONCURRENT_HPP_ */
