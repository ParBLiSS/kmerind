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
#ifndef BUFFER_HPP_A
#define BUFFER_HPP_

#include <string.h>

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
        size_t capacity;
        std::atomic<size_t> size;

        std::unique_ptr<unsigned char[], std::default_delete<unsigned char[]> > data;

        mutable std::mutex mutex;
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
        Buffer(const size_t _capacity) : capacity(_capacity), size(0), data(new unsigned char[_capacity]) {
          assert(_capacity > 0);
        };
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
        Buffer(Buffer<bliss::concurrent::THREAD_SAFE>&& other) : Buffer<bliss::concurrent::THREAD_SAFE>(std::move(other), std::lock_guard<std::mutex>(other.mutex) ) {};

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
        Buffer(Buffer<bliss::concurrent::THREAD_UNSAFE>&& other);

        /**
         * move assignment operator
         * @param other
         * @return
         */
        Buffer<bliss::concurrent::THREAD_SAFE>& operator=(Buffer<bliss::concurrent::THREAD_UNSAFE>&& other);


        /// copy ctor/assign are automatically NOT generated because of available destructor and move ctor/assign.
        /**
         * get the current size of the buffer.
         * @return
         */
        size_t getSize() {
          return size.load(std::memory_order_consume);
        }

        /**
         * get the capacity of the buffer.
         * @return
         */
        size_t getCapacity() {
          return capacity;
        }

        /**
         * note that the data should not be deleted by something else.  also, read it.
         * read only.
         */
        const unsigned char* getData() const {
          return data.get();
        }



        /**
         * clear the buffer. (just set the size to 0)
         */
        void clear() {
          size.store(0, std::memory_order_release);
        }

        /**
         * check if a buffer is full.
         * @return
         */
        bool isFull() {
          size_t s = size.load(std::memory_order_consume);
          return s >= capacity;
        }


        /**
         * append data to the buffer.
         * @param add_data
         * @param add_size
         */
        bool append(const void* add_data, const size_t add_size) {
          // can't use memory ordering alone, because we need
          // to check if we have room to add.  if not, then don't add.
          // the check and add have to happen atomically, so fetch_add does not work.
          if (capacity == 0) return false;

          std::unique_lock<std::mutex> lock(mutex);
          size_t s = size.load(std::memory_order_relaxed);  // no memory ordering within mutex lock
          size_t newS = s + add_size;
          if (newS > capacity) {
            return false;
          }

          size.store(newS, std::memory_order_relaxed);
          lock.unlock();

          memcpy(data.get() + s, add_data, add_size);
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
        size_t capacity;
        size_t size;

        std::unique_ptr<unsigned char[], std::default_delete<unsigned char[]> > data;

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
        Buffer(const size_t _capacity) : capacity(_capacity), size(0), data(new unsigned char[_capacity]) {
          assert(_capacity > 0);
        };
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
        Buffer(Buffer<bliss::concurrent::THREAD_SAFE>&& other) : Buffer<bliss::concurrent::THREAD_UNSAFE>(std::move(other), std::lock_guard<std::mutex>(other.mutex)) {};
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
        Buffer(Buffer<bliss::concurrent::THREAD_UNSAFE>&& other) : capacity(other.capacity), size(other.size), data(std::move(other.data)) {
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


        /// copy ctor/assign are automatically NOT generated because of available destructor and move ctor/assign.
        /**
         * get the current size of the buffer.
         * @return
         */
        size_t getSize() {
          return size;
        }

        /**
         * get the capacity of the buffer.
         * @return
         */
        size_t getCapacity() {
          return capacity;
        }

        /**
         * note that the data should not be deleted by something else.  also, read it.
         * read only.
         */
        const unsigned char* getData() const {
          return data.get();
        }


        /**
         * clear the buffer. (just set the size to 0)
         */
        void clear() {
          size = 0;
        }

        /**
         * check if buffer if full.
         * @return
         */
        bool isFull() {
          return size >= capacity;
        }

        /**
         * append data to the buffer.
         * @param add_data
         * @param add_size
         * @return
         */
        bool append(const void* add_data, const size_t add_size) {
          if (capacity == 0) return false;

          if ((size + add_size) > capacity) return false;

          memcpy(data.get() + size, add_data, add_size);
          size += add_size;
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
