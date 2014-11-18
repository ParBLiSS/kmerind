/**
 * @file    Buffer.hpp
 * @ingroup bliss::io
 * @author  tpan
 * @brief   fixed sized memory buffer declaration and implementations.
 * @details templated memory buffer classes with fixed sized byte array for storage, including thread safe and thread unsafe declarations and implementations.
 *
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef BUFFER_HPP_
#define BUFFER_HPP_

#include <cassert>
#include <cstring>   // memset
//#include <cstdlib>

#include <atomic>
#include <mutex>

#include <memory>     // unique_ptr
#include <utility>    // move, forward, swap, make_pair
#include <stdexcept>
#include <iostream>   //for std::cout

#include "concurrent/concurrent.hpp"   // LockType boolean constants
#include "utils/logging.h"

#include <xmmintrin.h>  // _mm_pause, instead of usleep.

namespace bliss
{
  namespace io
  {

    /**
     * @brief Thread safe/unsafe Memory Buffer, which is a fixed size allocated block of memory where data can be appended to, and read from via pointer.
     * @details this class uses unique_ptr internally to manage the memory, and supports only MOVE semantics.
     *
     * Thread Safety is enforced using atomic variables - writing to offsets that are atomically computed.
     * Memory ordering uses seq_cst at the moment
     *
     * Only when necessary, mutex lock is used. (in move constructor, move assignment operator, destructor, and clear function)
     *
     * thread-safe and thread-unsafe versions can be move constructed/assigned from each other.
     *
     * life cycle:  pool acquire():  buffer unblock().  -> allow rw.
     *              application:  buffer block() -> allow ro
     *              pool release():  buffer clear().
     *
     * There are 3 types of threads:  target write area is
     *        1. completely within buffer:  these threads should proceed to memcpy with the local ptr var
     *        2. completely outside buffer: these threads will not memcpy.  they all disabled buffer, and
     *           only 1 thread from 2) or 3) should swap in a new buffer atomically.
     *            - if curr buffer is not disabled (already swapped), then don't swap further.
     *        3. crossing buffer boundary: this thread will not memcpy.  disabled.  and it should retract the pointer advance.
     *
     *
     * @note  Buffer swap and free is expected to be conducted by a single thread.  other threads may access the buffer via pointer during free or swap.
     *         It is therefore important to swap via atomic operation and to free outside of parallel region (symptom - heap corruption, write after free to reserved)
     *         reserve function has been updated to prevent write if reserved is null, but this is not a guarantee.
     *            2 simple alternatives:
     *             1. have each thread hold its own buffer
     *             2. save the buffer into a list/vector and clean out later.
     *         also not a guarantee - using compare-exchange in the presence of allocator that can reuse memory - could result in ABA problem.
     *            Solution to this is to use cmpxchg16b instruction with custom reference counted pointers. (gcc intrinsic, enabled with -march=native)
     *
     * @tparam  LockType  controls whether this class is thread safe or not
     */
    template<bliss::concurrent::LockType LockType>
    class Buffer
    {
        /*
         * Declare Buffer<!LockType> class as a friend class, so we can reference its member functions and variables
         * in move constructor and assignment operators that moves data between instances with different thread safety.
         *
         * alternative is to use inheritance, which would require virtual functions that are potentially expensive.
         */
        friend class Buffer<bliss::concurrent::LockType::NONE> ;
        friend class Buffer<bliss::concurrent::LockType::MUTEX> ;
        friend class Buffer<bliss::concurrent::LockType::SPINLOCK> ;
        friend class Buffer<bliss::concurrent::LockType::LOCK_FREE> ;

        template<bliss::concurrent::LockType TS>
        friend std::ostream & operator<<(std::ostream &os, const Buffer<TS>& p);

      protected:

        /// internal data storage
        mutable uint8_t* start_ptr; // const, does not change

        /// start pointer shifted by maximum capacity.  unrealistic to use size_t - can't possibly allocate.
        mutable int64_t capacity;   // const, does not change

        /// pointer to current head of reservation
        volatile typename std::conditional<
            LockType == bliss::concurrent::LockType::LOCK_FREE,
            std::atomic<int64_t>, int64_t>::type reserved;

        /// pointer to FINAL  end of data.   only updated when buffer is blocked or when full (from reserved)
        mutable typename std::conditional<
            LockType == bliss::concurrent::LockType::LOCK_FREE,
            std::atomic<int64_t>, int64_t>::type size;

        /// represent amount of data written.  will not update beyond the FINAL size
        volatile typename std::conditional<
            LockType == bliss::concurrent::LockType::LOCK_FREE,
            std::atomic<int64_t>, int64_t>::type written;

        /// mutex for locking access to the buffer.  available in both thread safe and unsafe versions so we on't need to extensively enable_if or inherit
        mutable std::mutex mutex;
        mutable std::atomic_flag spinlock = ATOMIC_FLAG_INIT;

      private:
        /**
         * @brief Move constructor with a mutex lock on the source object.
         * @details This constructor is private and only use as constructor delegation target.
         *  the source object is locked before this function is called and data moved to the newly constructed object.
         * This version copies between 2 objects with the SAME thread safety.
         *
         * @param other   the source Buffer
         * @param l       the mutex lock on the source Buffer.
         */
        Buffer(Buffer<LockType> && other, const std::lock_guard<std::mutex> &l)
            : start_ptr(other.start_ptr), capacity(other.capacity),
              reserved((int64_t)(other.reserved)), size((int64_t)(other.size)),
              written((int64_t)(other.written))
        {

          other.start_ptr = nullptr;
          other.size = 0;
          other.capacity = 0;
          other.written = 0;
          other.reserved = capacity + 1;

        }
        ;

      public:
        /**
         * @brief Normal constructor.  Allocate and initialize memory with size specified as parameter.
         * @param _capacity   The maximum capacity of the Buffer in bytes
         */
        explicit Buffer(const uint32_t _capacity)
        {
          if (_capacity == 0)
            throw std::invalid_argument(
                "Buffer constructor parameter capacity is given as 0");

          std::lock_guard<std::mutex> lock(mutex);
          start_ptr = new uint8_t[_capacity](); // parenthesis initializes the memory to 0

          // buffer empty
          written = 0;
          capacity = _capacity;
          // buffer blocked.
          size = 0;
          reserved = capacity + 1;   // block from insertion.

        }
        ;

        /**
         * @brief Destructor.  waits for all writes and then deallocate memory manually.
         */
        virtual ~Buffer()
        {

          block_and_flush();

          std::lock_guard<std::mutex> lock(mutex);
          // blocked.
          reserved = capacity + 1;
          size = 0;

          capacity = 0;
          written = 0;

          if (start_ptr != nullptr)
          {
            delete[] start_ptr;
            start_ptr = nullptr;
          }
        }
        ;

        /**
         * @brief Move constructs from a Buffer with the SAME LockType property
         * @details  internal data memory moved by std::unique_ptr semantics.
         * Thread Safe always, using approach from http://www.justsoftwaresolutions.co.uk/threading/thread-safe-copy-constructors.html
         *
         * Constructs a mutex lock and then delegates to another constructor.
         *
         *
         * @param other   Source object to move
         */
        explicit Buffer(Buffer<LockType> && other)
            : Buffer<LockType>(std::move(other),
                               std::lock_guard<std::mutex>(other.mutex))
        {
        }
        ;

        /**
         * @brief Move assignment operator, between Buffers of the SAME LockType property.
         * @details  Internal data memory moved by std::unique_ptr semantics.
         * The move is done in a thread safe way always.
         *
         * @param other     Source Buffer to move
         * @return          target Buffer reference
         */
        Buffer<LockType>& operator=(Buffer<LockType> && other)
        {

          if (this->start_ptr != other.start_ptr)
          {

            std::unique_lock<std::mutex> myLock(mutex, std::defer_lock),
                otherLock(other.mutex, std::defer_lock);
            std::lock(myLock, otherLock);

            if (start_ptr != nullptr)
            {
              delete[] start_ptr;
              start_ptr = nullptr;
            }

            /// move the internal memory.

            start_ptr = other.start_ptr;
            written = (int64_t)(other.written);
            size = (int64_t)(other.size);
            capacity = other.capacity;
            reserved = (int64_t)(other.reserved);

            other.size = 0;
            other.start_ptr = nullptr;
            other.capacity = 0;
            other.written = 0;
            other.reserved = capacity + 1;
          }
          return *this;
        }

        /// remove copy constructor and copy assignement operators.
        explicit Buffer(const Buffer<LockType>& other) = delete;
        /// remove copy constructor and copy assignement operators.
        Buffer<LockType>& operator=(const Buffer<LockType>& other) = delete;

        /// remove default constructor
        Buffer() = delete;

        /**
         * @brief get the current size of the Buffer.
         * @detail  implemented as pointer difference because this is less frequently called compared to append.
         * @return    current size
         */
        const int64_t getSize() const
        {
          return (int64_t)size;
        }

      protected:

        /**
         * @brief get the current approximate size of the Buffer.
         * @return    current size
         */
        const int64_t getApproximateSize() const
        {
          return (int64_t)reserved;
        }

        /**
         * @brief get the current written data size of the Buffer.
         * @note  written data may be scattered.
         * @return    current size
         */
        const int64_t getWrittenSize() const
        {
          return (int64_t)written;
        }

      public:

        /**
         * @brief Get a pointer to the buffer data memory block.
         *
         * @note the data should not be deleted by a calling function/thread.
         * The access is read only.  There is no reason to return the unique_ptr.
         *
         * const because the caller will have a const reference to the buffer
         */
        template<typename T>
        operator T*() const
        {
          return reinterpret_cast<T*>(start_ptr);
        }

        /**
         * @brief get the capacity of the buffer.
         * @return    maximum capacity of buffer.
         */
        const int64_t getCapacity() const
        {
          return capacity;
        }

        /**
         * @brief Checks if a buffer is empty.
         * @note The return value is approximate due to threading
         *
         * @return    true if the buffer is empty, false otherwise.
         */
        const bool isEmpty() const
        {
          return (getWrittenSize() == 0);
        }

        // TODO:  make this more clear.
        // need 1. flag to indicate no more new writers
        //      2. wait for current writers to finish
        //      3. 3 phases:  writing, flushing, reading.  -> reserve, wirte_unlock, flush_begin, block, unblock
        //          reserve (counter ++)
        //          wirte_unlock (counter --)
        //          flush (MSB = 1)
        //          block ( wait for MSB==1 && counter == 0)
        //          unblock ( MSB = 0)
        // state checks:  is_writing (MSB==0|1 && counter > 0)
        //                is_flushing (MSB==1 && counter> 0)
        //                is_reading ( MSB==1 && counter == 0)

        template<bliss::concurrent::LockType TS = LockType>
        typename std::enable_if<TS == bliss::concurrent::LockType::MUTEX, bool>::type is_writing() const
        {
          std::lock_guard<std::mutex> lock(mutex);
          return is_writing<bliss::concurrent::LockType::NONE>();
        }

        template<bliss::concurrent::LockType TS = LockType>
        typename std::enable_if<TS == bliss::concurrent::LockType::MUTEX, bool>::type is_full() const
        {
          std::lock_guard<std::mutex> lock(mutex);
          return is_full<bliss::concurrent::LockType::NONE>();
        }

        template<bliss::concurrent::LockType TS = LockType>
        typename std::enable_if<TS == bliss::concurrent::LockType::MUTEX, bool>::type is_flushing() const
        {
          std::lock_guard<std::mutex> lock(mutex);
          return is_flushing<bliss::concurrent::LockType::NONE>();
        }

        template<bliss::concurrent::LockType TS = LockType>
        typename std::enable_if<TS == bliss::concurrent::LockType::MUTEX, bool>::type is_reading() const
        {
          std::lock_guard<std::mutex> lock(mutex);
          return is_reading<bliss::concurrent::LockType::NONE>();
        }

        template<bliss::concurrent::LockType TS = LockType>
        typename std::enable_if<TS == bliss::concurrent::LockType::SPINLOCK,
            bool>::type is_writing() const
        {
          while (spinlock.test_and_set())
            ;
          bool res = is_writing<bliss::concurrent::LockType::NONE>();
          spinlock.clear();
          return res;
        }

        template<bliss::concurrent::LockType TS = LockType>
        typename std::enable_if<TS == bliss::concurrent::LockType::SPINLOCK,
            bool>::type is_full() const
        {
          while (spinlock.test_and_set())
            ;
          bool res = is_full<bliss::concurrent::LockType::NONE>();
          spinlock.clear();
          return res;
        }

        template<bliss::concurrent::LockType TS = LockType>
        typename std::enable_if<TS == bliss::concurrent::LockType::SPINLOCK,
            bool>::type is_flushing() const
        {
          while (spinlock.test_and_set())
            ;
          bool res = is_flushing<bliss::concurrent::LockType::NONE>();
          spinlock.clear();
          return res;
        }

        template<bliss::concurrent::LockType TS = LockType>
        typename std::enable_if<TS == bliss::concurrent::LockType::SPINLOCK,
            bool>::type is_reading() const
        {
          while (spinlock.test_and_set())
            ;
          bool res = is_reading<bliss::concurrent::LockType::NONE>();
          spinlock.clear();
          return res;
        }

        template<bliss::concurrent::LockType TS = LockType>
        typename std::enable_if<TS == bliss::concurrent::LockType::NONE, bool>::type is_writing() const
        {
          return written < size;
        }

        template<bliss::concurrent::LockType TS = LockType>
        typename std::enable_if<TS == bliss::concurrent::LockType::NONE, bool>::type is_full() const
        {
          return reserved > capacity;
        }

        template<bliss::concurrent::LockType TS = LockType>
        typename std::enable_if<TS == bliss::concurrent::LockType::NONE, bool>::type is_flushing() const
        {
          return (reserved > capacity) && (written < size);
        }

        template<bliss::concurrent::LockType TS = LockType>
        typename std::enable_if<TS == bliss::concurrent::LockType::NONE, bool>::type is_reading() const
        {
          return (reserved > capacity) && (written >= size);
        }

        // TODO:  make this more clear.
        // need 1. flag to indicate no more new writers
        //      2. wait for current writers to finish
        //      3. 3 phases:  writing, flushing, reading.  -> reserve, wirte_unlock, flush_begin, block, unblock
        //          reserve (counter ++)
        //          wirte_unlock (counter --)
        //          flush (MSB = 1)
        //          block ( wait for MSB==1 && counter == 0)
        //          unblock ( MSB = 0)
        // state checks:  is_writing (MSB==0|1 && counter > 0)
        //                is_flushing (MSB==1 && counter> 0)
        //                is_reading ( MSB==1 && counter == 0)

        template<bliss::concurrent::LockType TS = LockType>
        typename std::enable_if<TS == bliss::concurrent::LockType::LOCK_FREE,
            bool>::type is_writing() const
        {
          std::atomic_thread_fence(std::memory_order_seq_cst);
          return (int64_t)written < (int64_t)size;
        }

        template<bliss::concurrent::LockType TS = LockType>
        typename std::enable_if<TS == bliss::concurrent::LockType::LOCK_FREE,
            bool>::type is_full() const
        {
          std::atomic_thread_fence(std::memory_order_seq_cst);
          return (int64_t)reserved > capacity;
        }

        template<bliss::concurrent::LockType TS = LockType>
        typename std::enable_if<TS == bliss::concurrent::LockType::LOCK_FREE,
            bool>::type is_flushing() const
        {
          return (is_full<TS>()) && (is_writing<TS>());
        }

        template<bliss::concurrent::LockType TS = LockType>
        typename std::enable_if<TS == bliss::concurrent::LockType::LOCK_FREE,
            bool>::type is_reading() const
        {
          return (is_full<TS>()) && (!is_writing<TS>());
        }

      protected:

        /**
         * @brief       reserves a position in the buffer at which to insert the count number of bytes.
         * @details     increment only if we have room (so done via atomic CAS.).  if full, return nullptr. if becoming full, set the size.  else return the insertion point
         * @param count
         * @return      pointer at which to insert.
         */
        template<bliss::concurrent::LockType TS = LockType>
        typename std::enable_if<TS == bliss::concurrent::LockType::MUTEX,
            int64_t>::type reserve(const uint32_t count)
        {

          std::unique_lock<std::mutex> lock(mutex);

          // this part reduces likelihood of accessing freed memory (a different thread tries to reserve to a freed region, i.e. )
          // does not prevent ABA problem (if the memory is reallocated.)
          int64_t curr = reserved;

          if (curr > capacity)
          {  // full.  another thread saved to size
            lock.unlock();
            return -1;
          }
          else
          { //if curr == capacity, buffer must be FILLED in this cycle.  need to process it as if not full yet.
            reserved += count;
            lock.unlock();
            if (curr + count > capacity)
            {  // just filled.  set size.  reserved is already set.
              // multiple threads could reach here if they have different counts? NO. because of CAS is atomic.

              // since only 1 thread reaching here, just set size, no need to lock
              size = curr;

            }  // else normal append
          }
          return curr;
        }
        /**
         * @brief       reserves a position in the buffer at which to insert the count number of bytes.
         * @details     increment only if we have room (so done via atomic CAS.).  if full, return nullptr. if becoming full, set the size.  else return the insertion point
         * @param count
         * @return      pointer at which to insert.
         */
        template<bliss::concurrent::LockType TS = LockType>
        typename std::enable_if<TS == bliss::concurrent::LockType::SPINLOCK,
            int64_t>::type reserve(const uint32_t count)
        {

          while (spinlock.test_and_set())
            ;

          // this part reduces likelihood of accessing freed memory (a different thread tries to reserve to a freed region, i.e. )
          // does not prevent ABA problem (if the memory is reallocated.)
          int64_t curr = reserved;

          if (curr > capacity)
          {  // full.  another thread saved to size
            spinlock.clear();
            return -1;
          }
          else
          { //if curr == capacity, buffer must be FILLED in this cycle.  need to process it as if not full yet.
            reserved += count;
            spinlock.clear();
            if (curr + count > capacity)
            {  // just filled.  set size.  reserved is already set.
              // multiple threads could reach here if they have different counts? NO. because of CAS is atomic.

              // since only 1 thread reaching here, just set size, no need to lock
              size = curr;

            }  // else normal append
          }
          return curr;
        }

        /**
         * @brief       reserves a position in the buffer at which to insert the count number of bytes.
         * @details     increment only if we have room (so done via atomic CAS.).  if full, return nullptr. if becoming full, set the size.  else return the insertion point
         * @param count
         * @return      pointer at which to insert.
         */
        template<bliss::concurrent::LockType TS = LockType>
        typename std::enable_if<TS == bliss::concurrent::LockType::LOCK_FREE,
            int64_t>::type reserve(const uint32_t count)
        {

          // this part reduces likelihood of accessing freed memory (a different thread tries to reserve to a freed region, i.e. )
          // does not prevent ABA problem (if the memory is reallocated.)
          std::atomic_thread_fence(std::memory_order_seq_cst);
          int64_t curr = reserved.load();

          // increment
          bool exchanged = false;
          while (curr <= capacity && !exchanged)
          {
            exchanged = reserved.compare_exchange_weak(curr, curr + count);
          }
          std::atomic_thread_fence(std::memory_order_seq_cst);

          if (curr > capacity)
          {  // full.  another thread saved to size
            return -1;
          }
          else
          { //if curr == capacity, buffer must be FILLED in this cycle.  need to process it as if not full yet.
            if (curr + count > capacity)
            {  // just filled.  set size.  reserved is already set.
              // multiple threads could reach here if they have different counts? NO. because of CAS is atomic.

              // since only 1 thread reaching here, just set size, no need to check to make sure we store minimum.
              size.store(curr);

            }  // else normal append
          }
          return curr;
        }

        template<bliss::concurrent::LockType TS = LockType>
        typename std::enable_if<TS == bliss::concurrent::LockType::NONE, int64_t>::type reserve(
            const uint32_t count)
        {
          int64_t curr = reserved;
          if (curr > capacity)
          {
            return -1;
          }
          else
          { //if ptr == capacity, buffer must be FILLED in this cycle.  need to process it as if not full yet.
            reserved += count;
            if (curr + count > capacity)
            {
              // again, only 1 thread eaches here, so

              size = curr;
            }  // else normal append

          }
          return curr;
        }

        /// marked write as completed.
        template<bliss::concurrent::LockType TS = LockType>
        typename std::enable_if<TS == bliss::concurrent::LockType::MUTEX, void>::type complete_write(
            const int count)
        {
          std::lock_guard<std::mutex> lock(mutex);

          complete_write<bliss::concurrent::LockType::NONE>(count);
        }

        /// marked write as completed.
        template<bliss::concurrent::LockType TS = LockType>
        typename std::enable_if<TS == bliss::concurrent::LockType::SPINLOCK,
            void>::type complete_write(const int count)
        {
          while (spinlock.test_and_set())
            ;

          complete_write<bliss::concurrent::LockType::NONE>(count);
          spinlock.clear();
        }

        /// marked write as completed.
        template<bliss::concurrent::LockType TS = LockType>
        typename std::enable_if<TS == bliss::concurrent::LockType::LOCK_FREE,
            void>::type complete_write(const int count)
        {
          std::atomic_thread_fence(std::memory_order_seq_cst);
          written += count;
        }

        template<bliss::concurrent::LockType TS = LockType>
        typename std::enable_if<TS == bliss::concurrent::LockType::NONE, void>::type complete_write(
            const int count)
        {
          written += count;
        }

      public:

        //===== read lock trumps write.
        /**
         * @brief  read lock prevents future writes.  read lock can be turned on while there are writes in progress.  once read lock is on, reserve will returns nullptr.
         *
         * @details purpose of this method is to swap the size and reserved, so that reserved is at max+1, and size is the smallest of all locking thread's pointers
         */
        template<bliss::concurrent::LockType TS = LockType>
        typename std::enable_if<TS == bliss::concurrent::LockType::MUTEX, void>::type block(
            int64_t _curr = -1)
        { // not using xor since it toggles and complete_write is not going to change the sign bit (so we don't need xor)
          // disable reserved and set the size
          std::lock_guard<std::mutex> lock(mutex);

          block<bliss::concurrent::LockType::NONE>(_curr);
        }
        //===== read lock trumps write.
        /**
         * @brief  read lock prevents future writes.  read lock can be turned on while there are writes in progress.  once read lock is on, reserve will returns nullptr.
         *
         * @details purpose of this method is to swap the size and reserved, so that reserved is at max+1, and size is the smallest of all locking thread's pointers
         */
        template<bliss::concurrent::LockType TS = LockType>
        typename std::enable_if<TS == bliss::concurrent::LockType::SPINLOCK,
            void>::type block(int64_t _curr = -1)
        { // not using xor since it toggles and complete_write is not going to change the sign bit (so we don't need xor)
          // disable reserved and set the size
          while (spinlock.test_and_set())
            ;
          block<bliss::concurrent::LockType::NONE>(_curr);
          spinlock.clear();
        }

        //===== read lock trumps write.
        /**
         * @brief  read lock prevents future writes.  read lock can be turned on while there are writes in progress.  once read lock is on, reserve will returns nullptr.
         *
         * @details purpose of this method is to swap the size and reserved, so that reserved is at max+1, and size is the smallest of all locking thread's pointers
         */
        template<bliss::concurrent::LockType TS = LockType>
        typename std::enable_if<TS == bliss::concurrent::LockType::LOCK_FREE,
            void>::type block(int64_t _curr = -1)
        { // not using xor since it toggles and complete_write is not going to change the sign bit (so we don't need xor)
          // disable reserved and put the content in size
          int64_t end = (int64_t)size;
          int64_t curr = reserved.exchange(capacity + 1);

          std::atomic_thread_fence(std::memory_order_seq_cst); // make sure all writes are done.

          curr = (_curr >= 0 ? _curr : curr);
          bool stop = false;

          while (end > curr && !stop)
          {
            stop = size.compare_exchange_weak(end, curr); // get smallest, indicating the first one to reach block
          }
          std::atomic_thread_fence(std::memory_order_seq_cst); // make sure all updates are done.
        }

        template<bliss::concurrent::LockType TS = LockType>
        typename std::enable_if<TS == bliss::concurrent::LockType::NONE, void>::type block(
            int64_t _curr = -1)
        {
          // disable reserved and set the size
          int64_t curr = (_curr >= 0 ? _curr : reserved);
          reserved = capacity + 1;

          if (size > curr)
          {
            size = curr;
          }
        }

        /**
         * @brief  read lock prevents future writes.  read lock can be turned on while there are writes in progress.  once read lock is on, reserve will returns nullptr.
         *
         * @details purpose of this method is to swap the size and reserved, so that reserved is at max+1, and size is the smallest of all locking thread's pointers
         */
        template<bliss::concurrent::LockType TS = LockType>
        typename std::enable_if<TS == bliss::concurrent::LockType::MUTEX, void>::type unblock()
        {
          std::lock_guard<std::mutex> lock(mutex);

          unblock<bliss::concurrent::LockType::NONE>();

        }

        /**
         * @brief  read lock prevents future writes.  read lock can be turned on while there are writes in progress.  once read lock is on, reserve will returns nullptr.
         *
         * @details purpose of this method is to swap the size and reserved, so that reserved is at max+1, and size is the smallest of all locking thread's pointers
         */
        template<bliss::concurrent::LockType TS = LockType>
        typename std::enable_if<TS == bliss::concurrent::LockType::SPINLOCK,
            void>::type unblock()
        {
          while (spinlock.test_and_set())
            ;

          unblock<bliss::concurrent::LockType::NONE>();
          spinlock.clear();

        }

        /**
         * @brief  read lock prevents future writes.  read lock can be turned on while there are writes in progress.  once read lock is on, reserve will returns nullptr.
         *
         * @details purpose of this method is to swap the size and reserved, so that reserved is at max+1, and size is the smallest of all locking thread's pointers
         */
        template<bliss::concurrent::LockType TS = LockType>
        typename std::enable_if<TS == bliss::concurrent::LockType::LOCK_FREE,
            void>::type unblock()
        {
          // put the size into curr if it hasn't been done.
          int64_t curr = (int64_t)reserved;
          int64_t end = size.exchange(capacity + 1);
          std::atomic_thread_fence(std::memory_order_seq_cst); // make sure all writes are done.

          bool stop = false;
          while (curr > capacity && !stop)
          {
            stop = reserved.compare_exchange_weak(curr, end);
          }
          std::atomic_thread_fence(std::memory_order_seq_cst); // make sure all writes are done.

        }

        template<bliss::concurrent::LockType TS = LockType>
        typename std::enable_if<TS == bliss::concurrent::LockType::NONE, void>::type unblock()
        {
          // put the size into curr if it hasn't been done.
          if (reserved > capacity)
            reserved = size;
          // then reset size to max
          size = capacity + 1;
        }

        /**
         * @brief   Sets the buffer as read only and wait for all pending writes.
         */
        void block_and_flush()
        {
          block<LockType>();
          while (is_writing<LockType>())
            _mm_pause();
        }

        /**
         * @brief Clears the buffer. (set the size to 0, leaving the capacity and the memory allocation intact)
         * @note  const because the caller will have a const reference to the buffer.
         */
        template<bliss::concurrent::LockType TS = LockType>
        typename std::enable_if<TS == bliss::concurrent::LockType::SPINLOCK,
            void>::type clear()
        {
          while (spinlock.test_and_set())
            ;
          clear<bliss::concurrent::LockType::NONE>();
          spinlock.clear();
        }

        template<bliss::concurrent::LockType TS = LockType>
        typename std::enable_if<TS != bliss::concurrent::LockType::SPINLOCK,
            void>::type clear()
        {
          std::lock_guard<std::mutex> lock(mutex);
          // blocked
          reserved = capacity + 1;
          size = 0;
          written = 0;
        }

        /**
         * @brief Append data to the buffer, THREAD SAFE, LOCK FREE.
         * @details The function updates the current occupied size of the Buffer
         *  and memcopies the supplied data into the internal memory block.
         *
         *  This is the THREAD SAFE version using MUTEX LOCK.
         *
         * Semantically, a "false" return value means there is not enough room for the new data, not that the buffer is flushing.
         *
         * method is const because the caller will have const reference to the Buffer.
         *
         * This method relies on the current reserved incrementing monotonically (by count) until it exceeds capacity, and the written
         * incrementing monotonically (by count) as threads finish writes, chasing the reserved thread.
         *
         * The first thread to exceed the capacity will set the size.  written will never exceed that, and reserved will not
         * be smaller than size.
         *
         *  single primitive type is better as it can be returned in 1 op.
         *
         * @tparam TS Choose thread safe vs unsafe implementation. defaults to same as the parent class.
         * @param[in] _data   pointer to data to be copied into the Buffer..  this SHOULD NOT be shared between threads
         * @param[in] count   number of bytes to be copied into the Buffer
         * @return            unsigned char, bit 0 indicated whether operation succeeded, bit 1 indicating whether buffer swap is needed.
         */
        unsigned int append(const void* _data, const uint32_t count)
        {

          if (count == 0 || _data == nullptr)
            throw std::invalid_argument("_data is nullptr");

          // reserve a spot.
          int64_t pos = reserve<LockType>(count);

          // if fails, then already full
          if (pos == -1)
          { // write lock returns nullptr if ptr points to outside of the range.
            // already read locked.

            return 0x0;  // full and not swapping.
          }
          else
          {  // pos starts inside allocated buffer range.

            if ((pos + count) > this->capacity)
            { // thread that filled the buffer

              while (is_writing<LockType>())
                _mm_pause();

              return 0x2;  // was not full, now full.

            }
            else
            { // can insert.  may make buffer full

              // write
              std::memcpy(start_ptr + pos, _data, count);
              complete_write<LockType>(count); // all full buffers lock the read and unlock the writer

              return 0x1;   // not full, successfully inserted.
            }

          }
        }

    };

    /**
     * @brief << operator to write out DataBlock object's actual data.
     * @tparam Iterator   Source data iterator type.
     * @tparam Range      Range data type
     * @tparam Container  container type for buffer.  defaults to std::vector.
     * @param[in/out] ost   output stream to which the content is directed.
     * @param[in]     db    BufferedDataBlock object to write out
     * @return              output stream object
     */
    template<bliss::concurrent::LockType LockType>
    std::ostream& operator<<(std::ostream& ost, const Buffer<LockType> & buffer)
    {
      ost << "LockType=" << static_cast<int>(LockType)
      << " BUFFER: data_ptr/size="
      << static_cast<const void *>(buffer.start_ptr) << "/"
      << (int64_t)(buffer.size) << " currptr/maxptr="
      << (int64_t)(buffer.reserved) << "/" << buffer.capacity << " written="
      << (int64_t)(buffer.written) << " approx,size/cap="
      << buffer.getApproximateSize() << "," << buffer.getSize() << "/"
      << buffer.getCapacity() << " R? " << (buffer.is_reading() ? "y" : "n")
      << " F? " << (buffer.is_flushing() ? "y" : "n") << " W? "
      << (buffer.is_writing() ? "y" : "n") << std::flush;

      return ost;
    }

  } /* namespace io */
} /* namespace bliss */

#endif /* BUFFER_HPP_ */
