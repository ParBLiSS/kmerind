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

#include <memory>
#include <atomic>
#include <cstdlib>
#include <string.h>
#include <utility>

namespace bliss
{
  namespace io
  {


    template<bool THREAD_SAFE>
    class Buffer;


    /**
     * Thread_safe version of buffer
     */
    template<>
    class Buffer<true>
    {
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
        Buffer(Buffer<true>&& other, const std::lock_guard<std::mutex> &) : capacity(other.capacity), size(other.size.load()), data(other.data) {
          other.size = 0;
          other.data = nullptr;
        };

      public:
        /**
         * Constructor.  allocation of memory array is automatic.
         * @param _capacity   The number of bytes to store.
         */
        Buffer(const size_t _capacity) : capacity(_capacity), size(0), data(new unsigned char[_capacity]) {};
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
        Buffer(Buffer<true>&& other) : Buffer<true>(other, std::lock_guard<std::mutex>(other.mutex) ) {};

        /**
         * move constructor from none-threadsafe buffer
         * @param other
         */
        Buffer(Buffer<false>&& other) : capacity(other.capacity), size(other.size), data(other.data) {
          other.size = 0;
          other.data = nullptr;
        };

        /**
         * Move assignment operator.  internal data memory moved by std::unique_ptr semantics.
         * need to be thread safe.  use lock_guard.
         * @param other
         * @return
         */
        Buffer<true>& operator=(Buffer<true>&& other) {
          if (this != other) {

            std::unique_lock<std::mutex> myLock(mutex, std::defer_lock),
                                          otherLock(other.mutex, std::defer_lock);
            std::lock(myLock, otherLock);

            /// move the internal memory.
            capacity = other.capacity;
            size = other.size;
            data = other.data;

            other.size = 0;
            other.data = nullptr;
          }
        }

        /**
         * move assignment operator
         * @param other
         * @return
         */
        Buffer<true>& operator=(Buffer<false>&& other) {
          if (this != other) {

            std::unique_lock<std::mutex> myLock(mutex);

            /// move the internal memory.
            capacity = other.capacity;
            size = other.size;
            data = other.data;

            other.size = 0;
            other.data = nullptr;
          }
        }


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
         * note that the data should not be deleted by something else.  just use it.
         */
        unsiqned char* getData() {
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
    class Buffer<false>
    {
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
        Buffer(Buffer<true>&& other, const std::lock_guard<std::mutex> &) : capacity(other.capacity), size(other.size.load()), data(other.data) {
          other.size = 0;
          other.data = nullptr;
        };

      public:
        /**
         * Constructor.  allocation of memory array is automatic.
         * @param _capacity   The number of bytes to store.
         */
        Buffer(const size_t _capacity) : capacity(_capacity), size(0), data(new unsigned char[_capacity]) {};
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
        Buffer(Buffer<true>&& other) : Buffer<false>(other, std::lock_guard<std::mutex>(other.mutex) ) {};

        /**
         * move constructor from none-threadsafe buffer
         * @param other
         */
        Buffer(Buffer<false>&& other) : capacity(other.capacity), size(other.size), data(other.data) {
          other.size = 0;
          other.data = nullptr;
        };

        /**
         * Move assignment operator.  internal data memory moved by std::unique_ptr semantics.
         * need to be thread safe.  use lock_guard.
         * @param other
         * @return
         */
        Buffer<true>& operator=(Buffer<true>&& other) {
          if (this != other) {

            std::unique_lock<std::mutex> otherLock(other.mutex);

            /// move the internal memory.
            capacity = other.capacity;
            size = other.size;
            data = other.data;

            other.size = 0;
            other.data = nullptr;
          }
        }

        /**
         * move assignment operator
         * @param other
         * @return
         */
        Buffer<true>& operator=(Buffer<false>&& other) {
          if (this != other) {

            /// move the internal memory.
            capacity = other.capacity;
            size = other.size;
            data = other.data;

            other.size = 0;
            other.data = nullptr;
          }
        }


        /// copy ctor/assign are automatically NOT generated because of available destructor and move ctor/assign.


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
          if (size + add_size >= capacity) return false;

          memcpy(data.get() + size, add_data, add_size);
          size += add_size;
          return true;
        }
    };


  } /* namespace io */
} /* namespace bliss */

#endif /* BUFFERPOOL_HPP_ */
