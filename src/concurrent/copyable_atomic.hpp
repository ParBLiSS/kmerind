/**
 * @file    copyable_atomic.hpp
 * @ingroup
 * @author  tpan
 * @brief
 * @details
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef COPYABLE_ATOMIC_HPP_
#define COPYABLE_ATOMIC_HPP_

#include <atomic>

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

        // initialization is not atomic
        constexpr copyable_atomic(T desired) noexcept : base(desired) {}

        // initialization is not atomic
        copyable_atomic(const base& other) noexcept : base(other.load(std::memory_order_relaxed)) {}

        // initialization is not atomic.  required copy constructor
        copyable_atomic(const copyable_atomic<T>& other) noexcept : base(other.load(std::memory_order_relaxed)) {}

        // nonvirtual destructor, since atomic does not have a virtual destructor.
        ~copyable_atomic() noexcept {
          base::~atomic();
        }

        T operator=(T desired) noexcept {
          base::store(desired, std::memory_order_relaxed);
          return desired;
        }
        T operator=(T desired) volatile noexcept {
          base::store(desired, std::memory_order_relaxed);
          return desired;
        }

        copyable_atomic& operator=(const base& other) noexcept {
          base::store(other.load(std::memory_order_relaxed), std::memory_order_relaxed);
          return *this;
        }
        copyable_atomic& operator=(const base& other) volatile noexcept {
          base::store(other.load(std::memory_order_relaxed), std::memory_order_relaxed);
          return *this;
        }

        copyable_atomic& operator=(const copyable_atomic<T>& other) noexcept {
          base::store(other.load(std::memory_order_relaxed), std::memory_order_relaxed);
          return *this;
        }
        copyable_atomic& operator=(const copyable_atomic<T>& other) volatile noexcept {
          base::store(other.load(std::memory_order_relaxed), std::memory_order_relaxed);
          return *this;
        }

    };

  } /* namespace concurrent */
} /* namespace bliss */

#endif /* COPYABLE_ATOMIC_HPP_ */
