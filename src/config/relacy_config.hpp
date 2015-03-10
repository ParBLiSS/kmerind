/**
 * @file    relacy_config.hpp
 * @ingroup
 * @author  tpan
 * @brief
 * @details
 *
 * Copyright (c) 2015 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef RELACY_CONFIG_HPP_
#define RELACY_CONFIG_HPP_

// NOTE: relacy simulates threads, but not actually spawning threads.
#if defined(USE_RELACY)

// std::unique_lock and lock_guard does not work with relacy's mutex API. need to redefine those.
// #include <mutex>  // for lock_guard and unique_lock, but needed to be before relacy so that mutex is redefined.



#include "relacy/relacy_std.hpp"
// VAR_T, VAR, and TLS_T all have been defined in relacy.hpp.  just use throughout the code.

// for non-atomic, non thread_local, VAR_T(T) defines a proxy type to T.  calls to T then requires x($), which in effect accesses the data through the proxy.
// for non-atomic, thread_local, TLS_T(T) defines a proxy type to T.  calls to T then requires x($), which in effect accesses the data through the proxy.
// for atomic (in other words the rl::atomic_proxy), x($) then access the relacy atomic class, which does requires a extra debuginfo param on all atomic operators.
// lession:  for atomic, use std::atomic as before
//            for non-atomic, non-thread_local variables, use VAR(T) x, then use VAR(x)
//            for non-atomic, thread_local variables, use TLS_T(T) x, then use VAR(x)

// NOTE:  relacy remaps it's atomic_proxy type and memory_order to std::'s atomic and memory_order, respectively.
//        this means that calls to std::atomic can remain the same.

// needed to make concurrentqueue not load std atomic, and load relacy instead
#define MCDBGQ_USE_RELACY

// relacy redefines "delete".  so for "CTOR(const type&) = delete;", need to use the form below.
// same for "type& operator=(const type&) = delete;".  form  "CTOR(const type&);"
#define DELETED_FUNC_DECL(F) F

// relacy's use of condition_variable also does not follow standard convention
#define CV_WAIT(CV, LOCK)  do { CV.wait(LOCK, $); } while (false)
#define CV_NOTIFY_ONE(CV)  do { CV.notify_one($); } while (false)
#define CV_NOTIFY_ALL(CV)  do { CV.notify_all($); } while (false)


// includes support for atomic_flags
#define INIT_ATOMIC_FLAG(f) std::atomic_flag f     // initializer is kind of special and we can't use copy or move constructor or assignment operators.
// concurrent queue's definition uses a atomic<bool>.  changed and included here for use.
// also, RelacyThreadExitNotifier is critical  - use in "void thread(unsigned thread_index)" method.
// else we will get error "long jump causes uninitialized stack frame"
#include <relacy_shims.h>

// std::unique_lock and lock_guard from <mutex> does not work with relacy's mutex API. need to redefine those.
//#include <system_error>  // can't use system_error.  there is some pure virtual function in here or deleted functions.  delete clashes with RL_DELETE_PROXY
#include <stdexcept>
namespace std
{

  /// @brief  Scoped lock idiom.
  // Acquire the mutex here with a constructor call, then release with
  // the destructor call in accordance with RAII style.
  template<typename _Mutex>
    class lock_guard
    {
    public:
      typedef _Mutex mutex_type;

      explicit lock_guard(mutex_type& __m) : _M_device(__m)
      { _M_device.lock($); }

      ~lock_guard()
      { _M_device.unlock($); }

      lock_guard(const lock_guard&);
      lock_guard& operator=(const lock_guard&);

    private:
      mutex_type&  _M_device;
    };

/// Do not acquire ownership of the mutex.
  struct defer_lock_t { };
  /// Try to acquire ownership of the mutex without blocking.
  struct try_to_lock_t { };
  constexpr defer_lock_t        defer_lock { };
  constexpr try_to_lock_t       try_to_lock { };



  /// unique_lock
  template<typename _Mutex>
    class unique_lock
    {
    public:
      typedef _Mutex mutex_type;

      unique_lock() noexcept
      : _M_device(0), _M_owns(false)
      { }

      explicit unique_lock(mutex_type& __m)
      : _M_device(&__m), _M_owns(false)
      {
        lock();
        _M_owns = true;
      }

      unique_lock(mutex_type& __m, defer_lock_t) noexcept
      : _M_device(&__m), _M_owns(false)
      { }

      unique_lock(mutex_type& __m, try_to_lock_t)
      : _M_device(&__m), _M_owns(_M_device->try_lock())
      { }

      ~unique_lock()
      {
        if (_M_owns)
          unlock();
      }

      unique_lock(const unique_lock&);
      unique_lock& operator=(const unique_lock&);

      unique_lock(unique_lock&& __u) noexcept
      : _M_device(__u._M_device), _M_owns(__u._M_owns)
      {
        __u._M_device = 0;
        __u._M_owns = false;
      }

      unique_lock& operator=(unique_lock&& __u) noexcept
      {
        if(_M_owns)
          unlock();

        unique_lock(std::move(__u)).swap(*this);

        __u._M_device = 0;
        __u._M_owns = false;

        return *this;
      }

      void
      lock()
      {
        if (!_M_device)
//          __throw_system_error(int(errc::operation_not_permitted));
          throw std::logic_error("operation not permitted");
        else if (_M_owns)
//          __throw_system_error(int(errc::resource_deadlock_would_occur));
          throw std::logic_error("resource deadlock would occur");
        else
          {
            _M_device->lock($);
            _M_owns = true;
          }
      }
      bool
      try_lock()
      {
        if (!_M_device)
//          __throw_system_error(int(errc::operation_not_permitted));
          throw std::logic_error("operation not permitted");
        else if (_M_owns)
//          __throw_system_error(int(errc::resource_deadlock_would_occur));
          throw std::logic_error("resource deadlock would occur");
        else
          {
            _M_owns = _M_device->try_lock($);
            return _M_owns;
          }
      }

      void
      unlock()
      {
        if (!_M_owns)
//          __throw_system_error(int(errc::operation_not_permitted));
          throw std::logic_error("operation not permitted");
        else if (_M_device)
          {
            _M_device->unlock($);
            _M_owns = false;
          }
      }
      mutex_type*
      release() noexcept
      {
        mutex_type* __ret = _M_device;
        _M_device = 0;
        _M_owns = false;
        return __ret;
      }


      void
      lock(rl::debug_info_param d)
      {
        if (!_M_device)
//          __throw_system_error(int(errc::operation_not_permitted));
          throw std::logic_error("operation not permitted");
        else if (_M_owns)
//          __throw_system_error(int(errc::resource_deadlock_would_occur));
          throw std::logic_error("resource deadlock would occur");
        else
          {
            _M_device->lock(d);
            _M_owns = true;
          }
      }
      bool
      try_lock(rl::debug_info_param d)
      {
        if (!_M_device)
//          __throw_system_error(int(errc::operation_not_permitted));
          throw std::logic_error("operation not permitted");
        else if (_M_owns)
//          __throw_system_error(int(errc::resource_deadlock_would_occur));
          throw std::logic_error("resource deadlock would occur");
        else
          {
            _M_owns = _M_device->try_lock(d);
            return _M_owns;
          }
      }

      void
      unlock(rl::debug_info_param d)
      {
        if (!_M_owns)
//          __throw_system_error(int(errc::operation_not_permitted));
          throw std::logic_error("operation not permitted");
        else if (_M_device)
          {
            _M_device->unlock(d);
            _M_owns = false;
          }
      }



    private:
      mutex_type*       _M_device;
      bool              _M_owns; // XXX use atomic_bool
    };

  template<typename _Lock>
    unique_lock<_Lock>
    __try_to_lock(_Lock& __l)
    { return unique_lock<_Lock>(__l, try_to_lock); }

  template<int _Idx, bool _Continue = true>
    struct __try_lock_impl
    {
      template<typename... _Lock>
        static void
        __do_try_lock(tuple<_Lock&...>& __locks, int& __idx)
        {
          __idx = _Idx;
          auto __lock = __try_to_lock(std::get<_Idx>(__locks));
          if (__lock.owns_lock())
            {
              __try_lock_impl<_Idx + 1, _Idx + 2 < sizeof...(_Lock)>::
                __do_try_lock(__locks, __idx);
              if (__idx == -1)
                __lock.release();
            }
        }
    };



  template<int _Idx>
    struct __try_lock_impl<_Idx, false>
    {
      template<typename... _Lock>
        static void
        __do_try_lock(tuple<_Lock&...>& __locks, int& __idx)
        {
          __idx = _Idx;
          auto __lock = __try_to_lock(std::get<_Idx>(__locks));
          if (__lock.owns_lock())
            {
              __idx = -1;
              __lock.release();
            }
        }
    };

  /** @brief Generic lock.
   *  @param __l1 Meets Mutex requirements (try_lock() may throw).
   *  @param __l2 Meets Mutex requirements (try_lock() may throw).
   *  @param __l3 Meets Mutex requirements (try_lock() may throw).
   *  @throw An exception thrown by an argument's lock() or try_lock() member.
   *  @post All arguments are locked.
   *
   *  All arguments are locked via a sequence of calls to lock(), try_lock()
   *  and unlock().  If the call exits via an exception any locks that were
   *  obtained will be released.
   */
  template<typename _L1, typename _L2, typename ..._L3>
    void
    lock(_L1& __l1, _L2& __l2, _L3&... __l3)
    {
      while (true)
        {
          unique_lock<_L1> __first(__l1);
          int __idx;
          auto __locks = std::tie(__l2, __l3...);
          __try_lock_impl<0, sizeof...(_L3)>::__do_try_lock(__locks, __idx);
          if (__idx == -1)
            {
              __first.release();
              return;
            }
        }
    }


}



#else

// #<new>               // new and delete
#include <cstdlib>      // malloc etc.
#include <cassert>      // assert
#include <atomic>       // atomic and memory_order
#include <mutex>
#include <condition_variable>

#define VAR_T(T) T      // type of normal variable, non-atomic
#define VAR(x) x        // use of normal variable, non-atomic
#define TLS_T(T) T      // type of thread local variable

// needed because relacy redefines "delete" so can't just delete constructor and assignment operators.
#define DELETED_FUNC_DECL(F) F = delete

// initializer is kind of special and we can't use copy or move constructor or assignment operators.
#define INIT_ATOMIC_FLAG(f) std::atomic_flag f = ATOMIC_FLAG_INIT


// relacy's use of condition_variable also does not follow standard convention
#define CV_WAIT(CV, LOCK)  do { CV.wait(LOCK); } while (false)
#define CV_NOTIFY_ONE(CV)  do { CV.notify_one(); } while (false)
#define CV_NOTIFY_ALL(CV)  do { CV.notify_all(); } while (false)

#endif

#endif /* RELACY_CONFIG_HPP_ */
