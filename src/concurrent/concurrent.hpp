/**
 * @file		concurrent.hpp
 * @ingroup concurrent
 * @author	Tony Pan <tpan7@gatech.edu>
 * @brief   this file contains predefined constants for concurrency and thread safety
 * @details
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef CONCURRENT_HPP_
#define CONCURRENT_HPP_

#include <atomic>

namespace bliss {
  namespace concurrent {

    /// Type for indicating Thread Safety.
    typedef bool ThreadSafety;

    /// Constant indicating Thread Safe
    constexpr bool THREAD_SAFE = true;

    /// Constant indicating Thread Unsafe
    constexpr bool THREAD_UNSAFE = false;


    /// types of thread safety mechanisms
    enum class LockType { NONE = 0, MUTEX = 1, SPINLOCK = 2, LOCKFREE = 4, THREADLOCAL = 8 };


    constexpr std::memory_order MO_RELAXED = std::memory_order_relaxed;
    constexpr std::memory_order MO_CONSUME = std::memory_order_consume;
    constexpr std::memory_order MO_ACQUIRE = std::memory_order_acquire;
    constexpr std::memory_order MO_RELEASE = std::memory_order_release;
    constexpr std::memory_order MO_ACQ_REL = std::memory_order_acq_rel;
    constexpr std::memory_order MO_SEQ_CST = std::memory_order_seq_cst;

  }
}



#endif /* CONCURRENT_HPP_ */
