/**
 * @file		Buffer.hpp
 * @ingroup
 * @author	tpan
 * @brief   a thread safe memory buffer.
 * @details provides a reusable, thread safe memory buffer.  uses atomic structure if "THREAD_SAFE" is turned on.
 *          using atomic instead of mutex (which are release acquire operations
 *
 *          default memory model is seq_cst.  we can avoid that sync.
 *
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef BUFFER_HPP_
#define BUFFER_HPP_

#include <cstring>
#include <cassert>

#include <memory>
#include <atomic>
#include <cstdlib>
#include <utility>
#include <mutex>

#include "concurrent/concurrent.hpp"

namespace bliss
{
  namespace io
  {


    template<bliss::concurrent::ThreadSafety THREAD_SAFE>
    class Buffer;


    /**
     * Thread_safe version of buffer
     *
     * uses capacity to enforce that we can write.
     */
    template<>
    class Buffer<bliss::concurrent::THREAD_SAFE>
    {
        /*
         * Friend the Buffer<bliss::concurrent::THREAD_UNSAFE> class, because we cannot reference the member function of the Buffer<bliss::concurrent::THREAD_UNSAFE> class
         * without first declarining it - would create circular reference.
         *
         * alternative is to use inheritance, which would require virtual functions that are potentially expensive.
         */
      friend class Buffer<bliss::concurrent::THREAD_UNSAFE>;

      protected:
        typedef std::unique_ptr<uint8_t[], std::default_delete<uint8_t[]> >  DataType;

        size_t capacity;

        mutable std::atomic<size_t> size;     // mutable so can call append from a const Buffer object.  conceptually size and data are constant members of the Buffer.
        DataType data;

        mutable std::mutex mutex;  // mutable so can lock
      private:
        /**
         * called after acquiring mutex on "other"
         * @param other
         * @param
         */
        Buffer(Buffer<bliss::concurrent::THREAD_SAFE>&& other, const std::lock_guard<std::mutex> &) : capacity(other.capacity), data(std::move(other.data)) {
          size = other.size.load(std::memory_order_relaxed);
          other.size = 0;
          other.capacity = 0;
          other.data = nullptr;
        };

      public:
        /**
         * Constructor.  allocation of memory array is automatic.
         * @param _capacity   The number of bytes to store.
         */
        explicit Buffer(const size_t _capacity) : capacity(_capacity), size(0), data(new uint8_t[_capacity]) {
          memset(data.get(), 0, _capacity * sizeof(uint8_t));
          assert(_capacity > 0);
        };

        Buffer(void* _data, const size_t count) : capacity(count), size(count), data(static_cast<uint8_t*>(_data)) {
          assert(count > 0);
        }

        /**
         * deallocation is automatic.
         */
        virtual ~Buffer() {};

        // move semantics
        /**
         * Move constructor.  internal data memory moved by std::unique_ptr semantics.
         * need to be thread safe.  using approach from http://www.justsoftwaresolutions.co.uk/threading/thread-safe-copy-constructors.html
         *
         * forward to another constructor.
         *
         * @param other
         */
        explicit Buffer(Buffer<bliss::concurrent::THREAD_SAFE>&& other) : Buffer<bliss::concurrent::THREAD_SAFE>(std::move(other), std::lock_guard<std::mutex>(other.mutex) ) {};

        /**
         * Move assignment operator.  internal data memory moved by std::unique_ptr semantics.
         * need to be thread safe.  use lock_guard.
         * @param other
         * @return
         */
        Buffer<bliss::concurrent::THREAD_SAFE>& operator=(Buffer<bliss::concurrent::THREAD_SAFE>&& other) {
          if (this->data != other.data) {

            std::unique_lock<std::mutex> myLock(mutex, std::defer_lock),
                                          otherLock(other.mutex, std::defer_lock);
            std::lock(myLock, otherLock);

            /// move the internal memory.
            capacity = other.capacity;
            size = other.size.load(std::memory_order_relaxed);
            data = std::move(other.data);

            other.capacity = 0;
            other.size = 0;
            other.data = nullptr;
          }
          return *this;
        }

        /**
         * move constructor from none-threadsafe buffer
         * @param other
         */
        explicit Buffer(Buffer<bliss::concurrent::THREAD_UNSAFE>&& other);

        /**
         * move assignment operator
         * @param other
         * @return
         */
        Buffer<bliss::concurrent::THREAD_SAFE>& operator=(Buffer<bliss::concurrent::THREAD_UNSAFE>&& other);

        explicit Buffer(const Buffer<bliss::concurrent::THREAD_SAFE>& other) = delete;
        Buffer<bliss::concurrent::THREAD_SAFE>& operator=(const Buffer<bliss::concurrent::THREAD_SAFE>& other) = delete;


        /// copy ctor/assign are automatically NOT generated because of available destructor and move ctor/assign.
        /**
         * get the current size of the buffer.
         * @return
         */
        size_t getSize() const {
          return size.load(std::memory_order_consume);
        }

        /**
         * get the capacity of the buffer.
         * @return
         */
        const size_t & getCapacity() const {
          return capacity;
        }

        /**
         * note that the data should not be deleted by something else.  also, read it.
         * read only.
         * no reason to pass back the data unique_ptr - no functions to be called by user.
         * can't template this since can't overload function by return type only.
         */
        const void* getData() const {
          return data.get();
        }


        /**
         * clear the buffer. (just set the size to 0)
         */
        void clear() const {
          size.store(0, std::memory_order_release);
        }

        /**
         * check if a buffer is full.  if return true, guaranteed to be full.  false does not guarantee - size may be modified before returning.
         * @return
         */
        bool isFull() const {
          size_t s = size.load(std::memory_order_consume);
          return s >= capacity;
        }
        bool isEmpty() const {
          size_t s = size.load(std::memory_order_consume);
          return s == 0;
        }


        /**
         * append data to the buffer.  semantically, a "false" return value does not mean
         *  the buffer if full.  it means that there is not enough room for the new data.
         * @param[in] add_data
         * @param[in] add_size
         * @return  bool indicated whether operation succeeded.
         */
        bool append(const void* typed_data, const size_t count) const {     // const method because the user will have const reference to Buffers.
                                                                         // because of the const-ness, we have mutable data and size.
          if (capacity == 0) return false;

          // can't use memory ordering alone
          // we need to check if we have room to add.  if not, then don't add.
          // the check and add have to happen atomically, so fetch_add does not work.
          std::unique_lock<std::mutex> lock(mutex);
          size_t s = size.load(std::memory_order_relaxed);  // no memory ordering within mutex lock
          size_t newS = s + count;
          if (newS > capacity) return false;
          size.store(newS, std::memory_order_relaxed);
          lock.unlock();

          std::memcpy(data.get() + s, typed_data, count);
          return true;
        }
        /**
         * lockfree version.  sync version is faster on 4 threads.  using compare_exchange_weak vs strong - no major difference.
         *    using relaxed and release vs acquire and acq_rel - no major difference.
         *
         *    DO NO USE.
         * @param[in] typed_data
         * @param[in] count
         * @return   success/fail
         */
        bool append_lockfree(const void* typed_data, const size_t count) const {    // const method because the user will have const reference to Buffers.
                                                                                 // because of the const-ness, we have mutable data and size.
          if (capacity == 0) return false;


          // can't use memory ordering alone
          // we need to check if we have room to add.  if not, then don't add.
          // the check and add have to happen atomically, so fetch_add does not work

          // try with compare and exchange weak.
          size_t s = size.load(std::memory_order_consume);
          size_t newS = s;
          do {
            newS = s + count;
            if (newS > capacity) return false;
          } while (!size.compare_exchange_strong(s, newS, std::memory_order_acq_rel, std::memory_order_consume));  // choosing strong, since we return.


          std::memcpy(data.get() + s, typed_data, count);
          return true;
        }

    };


    /**
     *  non-thread safe version.
     */
    template<>
    class Buffer<bliss::concurrent::THREAD_UNSAFE>
    {
        /*
         * Friend the Buffer<bliss::concurrent::THREAD_SAFE> class, because we cannot reference the member function of the Buffer<bliss::concurrent::THREAD_SAFE> class
         * without first declarining it - would create circular reference.
         *
         * alternative is to use inheritance, which would require virtual functions that are potentially expensive.
         */
        friend class Buffer<bliss::concurrent::THREAD_SAFE>;


      protected:
        typedef std::unique_ptr<uint8_t[], std::default_delete<uint8_t[]> >  DataType;


        size_t capacity;

        mutable size_t size;        // mutable so can call append on a const Buffer object. conceptually size and data are constant members of the Buffer.
        DataType data;      // mutable so can call append on a const Buffer object. conceptually size and data are constant members of the Buffer.


      private:
        /**
         * called after acquiring mutex on "other"
         * @param other
         * @param
         */
        Buffer(Buffer<bliss::concurrent::THREAD_SAFE>&& other, const std::lock_guard<std::mutex> &);

      public:
        /**
         * Constructor.  allocation of memory array is automatic.
         * @param _capacity   The number of bytes to store.
         */
        explicit Buffer(const size_t _capacity) : capacity(_capacity), size(0), data(new uint8_t[_capacity]) {
          memset(data.get(), 0, _capacity * sizeof(uint8_t));
          assert(_capacity > 0);
        };

        Buffer(void* _data, const size_t count) : capacity(count), size(count), data(static_cast<uint8_t*>(_data)) {
          assert(count > 0);
        }

        /**
         * deallocation is automatic because of unique_ptr.
         */
        virtual ~Buffer() {};

        // move semantics
        /**
         * Move constructor.  internal data memory moved by std::unique_ptr semantics.
         * need to be thread safe.  using approach from http://www.justsoftwaresolutions.co.uk/threading/thread-safe-copy-constructors.html
         *
         * forward to another constructor.
         *
         * @param other
         */
        explicit Buffer(Buffer<bliss::concurrent::THREAD_SAFE>&& other) : Buffer<bliss::concurrent::THREAD_UNSAFE>(std::move(other), std::lock_guard<std::mutex>(other.mutex)) {};
        /**
         * Move assignment operator.  internal data memory moved by std::unique_ptr semantics.
         * need to be thread safe.  use lock_guard.
         * @param other
         * @return
         */
        Buffer<bliss::concurrent::THREAD_UNSAFE>& operator=(Buffer<bliss::concurrent::THREAD_SAFE>&& other);

        /**
         * move constructor from none-threadsafe buffer
         * @param other
         */
        explicit Buffer(Buffer<bliss::concurrent::THREAD_UNSAFE>&& other) : capacity(other.capacity), size(other.size), data(std::move(other.data)) {
          other.size = 0;
          other.capacity = 0;
          other.data = nullptr;
        };

        /**
         * move assignment operator
         * @param other
         * @return
         */
        Buffer<bliss::concurrent::THREAD_UNSAFE>& operator=(Buffer<bliss::concurrent::THREAD_UNSAFE>&& other) {
          if (this->data != other.data) {

            /// move the internal memory.
            capacity = other.capacity;
            size = other.size;
            data = std::move(other.data);

            other.capacity = 0;
            other.size = 0;
            other.data = nullptr;
          }

          return *this;
        }

        explicit Buffer(const Buffer<bliss::concurrent::THREAD_UNSAFE>& other) = delete;
        Buffer<bliss::concurrent::THREAD_UNSAFE>& operator=(const Buffer<bliss::concurrent::THREAD_UNSAFE>& other) = delete;

        /// copy ctor/assign are automatically NOT generated because of available destructor and move ctor/assign.
        /**
         * get the current size of the buffer.
         * @return
         */
        size_t getSize() const {
          return size;
        }

        /**
         * get the capacity of the buffer.
         * @return
         */
        const size_t& getCapacity() const {
          return capacity;
        }

        /**
         * note that the data should not be deleted by something else.  also, read it.
         * read only.
         * no reason to pass back the data unique_ptr - no functions to be called by user.
         * can't template this since can't overload function by return type only.
         */
        const void* getData() const {
          return data.get();
        }


        /**
         * clear the buffer. (just set the size to 0)
         */
        void clear() const {
          size = 0;
        }

        /**
         * check if buffer if full.
         * @return
         */
        bool isFull() const {
          return size >= capacity;
        }

        bool isEmpty() const {
          return size == 0;
        }


        /**
         * append data to the buffer.  semantically, a "false" return value does not mean
         *  the buffer if full.  it means that there is not enough room for the new data.
         * @param[in] add_data
         * @param[in] add_size
         * @return success/failure
         */
        bool append(const void* typed_data, const size_t count) const {   // const method because the user will have const reference to Buffers.
                                                                       // because of the const-ness, we have mutable data and size.
          if (capacity == 0) return false;

          if ((size + count) > capacity) return false;

          std::memcpy(data.get() + size, typed_data, count);
          size += count;
          return true;
        }
    };

    ///////////// following functions are defined here because they use each other's definitions.



    //////////////////////// Thread_safe version of buffer
    /*
     * move constructor from none-threadsafe buffer
     * @param other
     */
    Buffer<bliss::concurrent::THREAD_SAFE>::Buffer(Buffer<bliss::concurrent::THREAD_UNSAFE>&& other) {
      capacity = other.capacity;
      size = other.size;
      data = std::move(other.data);

      other.capacity = 0;
      other.size = 0;
      other.data = nullptr;
    };

    /*
     * move assignment operator
     * @param other
     * @return
     */
    Buffer<bliss::concurrent::THREAD_SAFE>& Buffer<bliss::concurrent::THREAD_SAFE>::operator=(Buffer<bliss::concurrent::THREAD_UNSAFE>&& other) {
      if (this->data != other.data) {

        std::unique_lock<std::mutex> myLock(mutex);

        /// move the internal memory.
        capacity = other.capacity;
        size = other.size;
        data = std::move(other.data);

        other.capacity = 0;
        other.size = 0;
        other.data = nullptr;
      }

      return *this;
    };



    /////////////// non thread safe
    /*
     * move constructor from none-threadsafe buffer
     * @param other
     */
    Buffer<bliss::concurrent::THREAD_UNSAFE>::Buffer(Buffer<bliss::concurrent::THREAD_SAFE>&& other, const std::lock_guard<std::mutex> &)  :
      capacity(other.capacity), size(other.size.load(std::memory_order_relaxed)), data(std::move(other.data)) {
      other.size = 0;
      other.capacity = 0;
      other.data = nullptr;
    };

    /*
     * Move assignment operator.  internal data memory moved by std::unique_ptr semantics.
     * need to be thread safe.  use lock_guard.
     * @param other
     * @return
     */
    Buffer<bliss::concurrent::THREAD_UNSAFE>& Buffer<bliss::concurrent::THREAD_UNSAFE>::operator=(Buffer<bliss::concurrent::THREAD_SAFE>&& other) {
      if (this->data != other.data) {

        std::unique_lock<std::mutex> otherLock(other.mutex);

        /// move the internal memory.
        capacity = other.capacity;
        size = other.size.load(std::memory_order_relaxed);
        data = std::move(other.data);

        other.capacity = 0;
        other.size = 0;
        other.data = nullptr;
      }

      return *this;
    }




  } /* namespace io */
} /* namespace bliss */

#endif /* BUFFERPOOL_HPP_ */
