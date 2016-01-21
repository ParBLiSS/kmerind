/**
 * @file    copyable_atomic.hpp
 * @ingroup concurrent
 * @author  Tony Pan <tpan7@gatech.edu>
 * @brief   atomic type that supports copy consturctor, suitable for map container.
 * @details std::atomic types are compatible with vector, but not map container.
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef COPYABLE_ATOMIC_HPP_
#define COPYABLE_ATOMIC_HPP_

#include "config/relacy_config.hpp"

namespace bliss
{
  namespace concurrent
  {

    /**
     * @class    bliss::concurrent::copyable_atomic
     * @brief     subclass of std::atomic that is CopyContructible and CopyAssignable.
     * @details   all remaining operations use std::atomic's
     *
     */
    template <typename T>
    struct copyable_atomic : public std::atomic<T>
    {
      protected:
        using base = std::atomic<T>;

      public:
        copyable_atomic() noexcept : base() {}

        /// constructor. initialization is not atomic
        constexpr copyable_atomic(T desired) noexcept : base(desired) {}

        /// copy constructor (from std::atomic).  initialization is not atomic
        copyable_atomic(const base& other) noexcept : base(other.load(std::memory_order_relaxed)) {}

        /// copy constructor.  initialization is not atomic.  required copy constructor
        copyable_atomic(const copyable_atomic<T>& other) noexcept : base(other.load(std::memory_order_relaxed)) {}

        /// copy constructor.  initialization is not atomic.  required copy constructor
        copyable_atomic(copyable_atomic<T>&& other) noexcept : base(other.exchange(T(), std::memory_order_relaxed)) {
        }


        // nonvirtual destructor, since atomic has trivial destructors (does nothing, not virtual, members have trivial destructor)
        // in fact, make it trivial (implicitly defined.)
        //~copyable_atomic() noexcept {}

        /// assignemnt operator
        T operator=(T desired) noexcept {
          base::store(desired, std::memory_order_relaxed);
          return desired;
        }
        /// assignment operator
        T operator=(T desired) volatile noexcept {
          base::store(desired, std::memory_order_relaxed);
          return desired;
        }

        /// copy assign from std::atomic
        copyable_atomic& operator=(const base& other) noexcept {
          base::store(other.load(std::memory_order_relaxed), std::memory_order_relaxed);
          return *this;
        }
        /// copy assign from std::atomic
        copyable_atomic& operator=(const base& other) volatile noexcept {
          base::store(other.load(std::memory_order_relaxed), std::memory_order_relaxed);
          return *this;
        }

        /// copy assign
        copyable_atomic& operator=(const copyable_atomic<T>& other) noexcept {
          base::store(other.load(std::memory_order_relaxed), std::memory_order_relaxed);
          return *this;
        }
        /// copy assign
        copyable_atomic& operator=(const copyable_atomic<T>& other) volatile noexcept {
          base::store(other.load(std::memory_order_relaxed), std::memory_order_relaxed);
          return *this;
        }

        /// move assign
        copyable_atomic& operator=(copyable_atomic<T>&& other) noexcept {
          this->store(other.exchange(T(), std::memory_order_relaxed), std::memory_order_relaxed);
          return *this;
        }
        /// move assign
        copyable_atomic& operator=(copyable_atomic<T>&& other) volatile noexcept {
          this->store(other.exchange(T(), std::memory_order_relaxed), std::memory_order_relaxed);
          return *this;
        }


    };

  } /* namespace concurrent */
} /* namespace bliss */

#endif /* COPYABLE_ATOMIC_HPP_ */
