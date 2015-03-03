/**
 * @file    buffer.hpp
 * @ingroup concurrent
 * @author  Tony Pan <tpan7@gatech.edu>
 * @brief   fixed sized memory buffer declaration and implementations.
 * @details templated memory buffer classes with fixed sized byte array for storage, including thread safe and thread unsafe versions.
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
#include <xmmintrin.h>  // _mm_pause, instead of usleep.

#include <atomic>
#include <mutex>

#include <memory>     // unique_ptr
#include <utility>    // move, forward, swap, make_pair
#include <stdexcept>
#include <iostream>   //for std::cout

#include "concurrent/concurrent.hpp"   // LockType boolean constants
#include "utils/logging.h"


namespace bliss
{
  namespace io
  {


    /**
     * @brief Base Memory Buffer.  provides base storage and some common functionalities. All thread safety related implementations
     *        are handled in subclasses.
     * @note  does NOT use virtual functions for append, for performance reason.
     */
    template<int64_t Capacity, size_t MetadataSize>
    class BufferBase
    {
        static_assert(Capacity > 0, "Buffer Capacity is given as 0");

      protected:
        /// internal data storage.  includes metadata
        mutable uint8_t* start_ptr; // const, does not change
        
        /// pointer to where the actual data starts
        mutable uint8_t* data_ptr;  

        /// mutex for locking access to the buffer.  available in both thread safe and unsafe versions so we on't need to extensively enable_if or inherit
        mutable std::mutex mutex;

        /// pointer to FINAL  end of data.   only updated when buffer is blocked or when full (from reserved)
        /// capacity + 1 indicates buffer is unblocked.  other values indicate that buffer is blocked or full.
        ///  exchanges with reserved
        mutable int64_t size;


        /**
         * @brief Move constructor with a mutex lock on the source object.
         * @details This constructor is private and only use as constructor delegation target.
         *  the source object is locked before this function is called and data moved to the newly constructed object.
         *
         * @param other   the source Buffer
         * @param l       the mutex lock on the source Buffer.
         */
        BufferBase(BufferBase<Capacity, MetadataSize> && other, const std::lock_guard<std::mutex> &)
            : start_ptr(other.start_ptr), data_ptr(other.data_ptr), size(other.size)
        {
          other.start_ptr = nullptr;
          other.data_ptr = nullptr;
          other.size = 0;
          std::atomic_thread_fence(std::memory_order_release);
        };


        /**
         * @brief Move assignment operator, requiring mutex lock. For between Buffers of the SAME LockType property.
         * @details  Internal data memory moved.
         * The move is done in a thread safe way always.
         *
         * @param other     Source Buffer to move
         */
        void move_assigner(BufferBase<Capacity, MetadataSize> && other,
            const std::unique_lock<std::mutex> &, const std::unique_lock<std::mutex> &)
        {

          if (this->start_ptr != other.start_ptr)
          {
            // move the internal memory.
            start_ptr = other.start_ptr;
            data_ptr = other.data_ptr;

            std::atomic_thread_fence(std::memory_order_acquire);
            size = other.size;

            other.start_ptr = nullptr;
            other.data_ptr = nullptr;
            other.size = 0;
            std::atomic_thread_fence(std::memory_order_release);
          }
        }


        /// remove copy constructor and copy assignement operators.
        explicit BufferBase(const BufferBase<Capacity, MetadataSize>& other) = delete;
        /// remove copy constructor and copy assignement operators.
        BufferBase<Capacity, MetadataSize>& operator=(const BufferBase<Capacity, MetadataSize>& other) = delete;


        /**
         * @brief Constructor.  Allocate and initialize memory with size specified as parameter.
         */
        BufferBase() : start_ptr(new uint8_t[Capacity + MetadataSize]()), data_ptr(start_ptr + MetadataSize), size(0)
        {}
        //start_ptr(new uint8_t[Capacity](), []( uint8_t *p ) { delete [] p; })  if using shared ptr

        /**
         * @brief Destructor.  waits for all writes and then deallocate memory manually.
         */
        virtual ~BufferBase()
        {
          delete [] start_ptr;
        }

        /**
         * @brief Move constructs from a Buffer with the SAME LockType property
         * @details
         * Thread Safe always, using approach from http://www.justsoftwaresolutions.co.uk/threading/thread-safe-copy-constructors.html
         *
         * Constructs a mutex lock and then delegates to another constructor.
         *
         * @param other   Source object to move
         */
        explicit BufferBase(BufferBase<Capacity, MetadataSize> && other)
            : BufferBase<Capacity, MetadataSize>(std::move(other),
                               std::lock_guard<std::mutex>(other.mutex))
        {}





        /**
         * @brief Move assignment operator, between Buffers of the SAME LockType property.
         * @details  Internal data memory moved.
         * The move is done in a thread safe way always.
         *
         * @param other     Source Buffer to move
         * @return          target Buffer reference
         */
        BufferBase<Capacity, MetadataSize>& operator=(BufferBase<Capacity, MetadataSize> && other)
        {

          if (this->start_ptr != other.start_ptr)
          {

            std::unique_lock<std::mutex> myLock(mutex, std::defer_lock),
                otherLock(other.mutex, std::defer_lock);
            std::lock(myLock, otherLock);

            // delegate the assignment operator.
            this->move_assigner(std::forward<BufferBase<Capacity, MetadataSize> >(other), myLock, otherLock);
            
          }
          return *this;
        }

      public:

        /// Get size of buffer.  only meaningful after buffer has been "blocked"
        int64_t getSize() const
        {
          std::atomic_thread_fence(std::memory_order_acquire);
          return size;
        }

        /**
         * @brief pointer to the start of the metadata block in buffer
         * @note  DO NOT DEALLOCATE THE MEMORY SPACE RETURNED.
         *
         * const because the caller will have a const reference to the buffer
         */
        uint8_t* metadata_begin() const 
        {
          return start_ptr;
        }
  
        /// pointer to the end of the metadata block in buffer
        uint8_t* metadata_end() const
        {
          return begin();
        }


        /**
         * @brief Get a pointer to the buffer data memory block.
         *
         * @note  DO NOT DEALLOCATE THE MEMORY SPACE RETURNED.
         *
         * const because the caller will have a const reference to the buffer
         *
         */
        uint8_t* begin() const
        {
          return data_ptr;
        }
        
        /// pointer to the end of the data block in buffer.
        uint8_t* end() const
        {
          std::atomic_thread_fence(std::memory_order_acquire);
          return data_ptr + size;
        }

        /**
         * @brief cast the internal data buffer to type T and return pointer to user.
         * @note  this breaks the encapsulation a little.
         * @tparam  T   desired data type.
         */
        template<typename T>
        operator T*() const {
          std::atomic_thread_fence(std::memory_order_acquire);
          return reinterpret_cast<T*>(data_ptr);
        }

        /**
         * @brief   get the capacity of the buffer.
         * @return  maximum capacity of buffer.
         */
        const int64_t getCapacity() const
        {
          return Capacity;
        }


      	/// get the metadata size
      	const size_t getMetadataSize() const
       	{
      	  return MetadataSize;
       	}


    };


    /**
     * @brief   Mutex or Spinlocked Memory Buffer, which is a fixed size allocated block of memory where data can be appended to, and read from via pointer.
     * @details this class supports only MOVE semantics.
     *
     * Thread Safety is enforced using mutex or spinlock.  when write, reservation is first made and then memory is written.
     *
     * life cycle:  application acquires a buffer from objectpool.
     *              application threads write into buffer until full
     *              a single application thread swaps in a new, empty buffer from objectpool, and
     *                  dispatches the full buffer for communication
     *              communication thread completes and releases buffer back to object pool
     *
     * When application thread writes, there are 4 possible scenarios:
     *        1. reservation succeeds and is completely within the buffer.  write succeeds in this case.
     *        2. reservation fails as the target range is completely outside buffer.  write fails in this case.
     *        3. reservation succeeds and exactly fills the buffers.  write succeeds in this case, but an
     *            empty buffer is also swapped in and the full buffer dispatched for communication.
     *        4. reservation would span buffer maximum size limit.  write fails in this case, and an
     *            empty buffer is swapped in and the full buffer dispatched for communication.
     *
     * @note  because of the last 2 cases, the replacement of a buffer is done by a single thread only.
     *        however, because other threads may still have pointer to the full buffer, shared_ptr is a better choice.
     *
     *        symptoms of data races from swap include heap corruption, write after free errors as reported by ASAN.
     *        also, such swaps are subject to the ABA problem.
     *
     *        recommendation is to have each thread hold its own non-thread-safe buffer, or use CMPXCHG16B with reference counted pointers.
     *
     * @tparam  LockType  controls the type of thread safely.  See concurrent/concurrent.h for allowed values.
     * @tparam  Capacity  size of the allocated byte array
     */
    template<bliss::concurrent::LockType LockType, int64_t Capacity = 8192, size_t MetadataSize = 0>
    class Buffer : public BufferBase<Capacity, MetadataSize>
    {
        template<bliss::concurrent::LockType LT, int64_t C, size_t MDS>
        friend std::ostream & operator<<(std::ostream &os, const Buffer<LT, C, MDS>& p);

      protected:
        using BaseType = BufferBase<Capacity, MetadataSize>;

        // instead of a single "size" variable and a "blocked" variable (to prevent further writes to the buffer)
        //   (which require concurrent updates to both)
        // we use Capacity+1 to indicate a blocked state: if reserved is set to Capacity+1, it is blocked..
        // also, because reservation and actual write may occur separately, we need a "reserved" and a "written"
        //   for separate incrementing.
        // having the "written" variable allows us to wait for all writes to complete.


        /// position of current start of free space.  for reservation.
        /// capacity + 1 indicates buffer is blocked/full.  other values indicate that buffer is unblocked.
        /// exchanges with size.
        volatile int64_t reserved;

        /// amount of data written so far.  will not update beyond the FINAL size
        volatile int64_t written;


        /// spinlock for locking access to the buffer.  available in both thread safe and unsafe versions so we on't need to extensively enable_if or inherit
        mutable std::atomic_flag spinlock = ATOMIC_FLAG_INIT;

      private:
        /**
         * @brief Move constructor with a mutex lock on the source object.
         * @details This constructor is private and only use as constructor delegation target.
         *  the source object is locked before this function is called and data moved to the newly constructed object.
         *
         * @param other   the source Buffer
         * @param l       the mutex lock on the source Buffer.
         */
        Buffer(Buffer<LockType, Capacity, MetadataSize> && other, const std::lock_guard<std::mutex> &l)
            : BaseType(other, l), reserved(other.reserved), written(other.written)
        {
          other.written = 0;
          other.reserved = Capacity + 1;
          std::atomic_thread_fence(std::memory_order_release);

        };

        /// remove copy constructor and copy assignement operators.
        explicit Buffer(const Buffer<LockType, Capacity, MetadataSize>& other) = delete;
        /// remove copy constructor and copy assignement operators.
        Buffer<LockType, Capacity, MetadataSize>& operator=(const Buffer<LockType, Capacity, MetadataSize>& other) = delete;


      public:
        /**
         * @brief Constructor.  Allocate and initialize memory with size specified as parameter.
         */
        Buffer() : BaseType(),
         reserved(Capacity + 1), written(0L)
        {}
        //start_ptr(new uint8_t[Capacity](), []( uint8_t *p ) { delete [] p; })  if using shared ptr

        /**
         * @brief Destructor.  waits for all writes and then deallocate memory manually.
         */
        virtual ~Buffer()
        {
          block_and_flush();
        }

        /**
         * @brief Move constructs from a Buffer with the SAME LockType property
         * @details
         * Thread Safe always, using approach from http://www.justsoftwaresolutions.co.uk/threading/thread-safe-copy-constructors.html
         *
         * Constructs a mutex lock and then delegates to another constructor.
         *
         * @param other   Source object to move
         */
        explicit Buffer(Buffer<LockType, Capacity, MetadataSize> && other)
            : Buffer<LockType, Capacity, MetadataSize>(std::move(other),
                               std::lock_guard<std::mutex>(other.mutex))
        {}

        /**
         * @brief Move assignment operator, between Buffers of the SAME LockType property.
         * @details  Internal data memory moved.
         * The move is done in a thread safe way always.
         *
         * @param other     Source Buffer to move
         * @return          target Buffer reference
         */
        Buffer<LockType, Capacity, MetadataSize>& operator=(Buffer<LockType, Capacity, MetadataSize> && other)
        {

          if (this->start_ptr != other.start_ptr)
          {

            std::unique_lock<std::mutex> myLock(this->mutex, std::defer_lock),
                otherLock(other.mutex, std::defer_lock);
            std::lock(myLock, otherLock);

            // move the internal memory.
            std::atomic_thread_fence(std::memory_order_acquire);

            written = (int64_t)(other.written);
            reserved = (int64_t)(other.reserved);

            this->move_assigner(std::move(other), myLock, otherLock);

            other.written = 0;
            other.reserved = Capacity + 1;
            std::atomic_thread_fence(std::memory_order_release);
          }
          return *this;
        }


      protected:

        /**
         * @brief get the current approximate size of the Buffer.
         * @return    current size
         */
        int64_t getApproximateSize() const
        {
            std::atomic_thread_fence(std::memory_order_acquire);
          return reserved;
        }

        /**
         * @brief get the current written data size of the Buffer.
         * @note  written data may be scattered.
         * @return    current size
         */
        int64_t getWrittenSize() const
        {
            std::atomic_thread_fence(std::memory_order_acquire);
          return written;
        }

        /// Check if buffer has been blocked from further reservation but still waiting for all writes to finish. via mutex lock
        template<bliss::concurrent::LockType LT = LockType>
        typename std::enable_if<LT == bliss::concurrent::LockType::MUTEX, bool>::type is_flushing() const
        {
          std::lock_guard<std::mutex> lock(this->mutex);
          std::atomic_thread_fence(std::memory_order_acquire);
          return (reserved > Capacity) && (written < this->size);
        }

        /// Check if buffer has been blocked from further reservation but still waiting for all writes to finish. via spin lock
        template<bliss::concurrent::LockType LT = LockType>
        typename std::enable_if<LT == bliss::concurrent::LockType::SPINLOCK,
            bool>::type is_flushing() const
        {
          while (spinlock.test_and_set(std::memory_order_acquire));
          bool res = (reserved > Capacity) && (written < this->size);
          spinlock.clear(std::memory_order_release);
          return res;
        }


      public:


        /**
         * @brief   Checks if a buffer is empty.
         * @note    The return value is approximate due to threading
         *
         * @return  true if the buffer is empty, false otherwise.
         */
        bool isEmpty() const
        {
          return (getWrittenSize() == 0);
        }

        // TODO:  make this more clear.
        // need 1. flag to indicate no more new writers
        //      2. wait for current writers to finish
        //      3. 3 phases:  writing, flushing, reading.  -> reserve, wirte_unlock, flush_begin, block_writes, unblock_writes
        //          reserve (counter ++)
        //          wirte_unlock (counter --)
        //          flush (MSB = 1)
        //          block_writes ( wait for MSB==1 && counter == 0)
        //          unblock_writes ( MSB = 0)
        // state checks:  is_writing (MSB==0|1 && counter > 0)
        //                is_flushing (MSB==1 && counter> 0)
        //                is_read_only ( MSB==1 && counter == 0)

        /// check if there are threads writing to the buffer.  uses mutex lock
        template<bliss::concurrent::LockType LT = LockType>
        typename std::enable_if<LT == bliss::concurrent::LockType::MUTEX, bool>::type is_writing() const
        {
          std::lock_guard<std::mutex> lock(this->mutex);
          std::atomic_thread_fence(std::memory_order_acquire);

          return written < this->size;  // reserved > written if is writing. but reserve can increase while buffer is unblocked.
        }

        /// check if the buffer is blocked from further writes and no threads need to write to it.  uses mutex
        template<bliss::concurrent::LockType LT = LockType>
        typename std::enable_if<LT == bliss::concurrent::LockType::MUTEX, bool>::type is_read_only() const
        {
          std::lock_guard<std::mutex> lock(this->mutex);
          std::atomic_thread_fence(std::memory_order_acquire);

          return (reserved > Capacity) && (written >= this->size);
        }

        /// check if there are threads writing to the buffer.  uses mutex lock
        template<bliss::concurrent::LockType LT = LockType>
        typename std::enable_if<LT == bliss::concurrent::LockType::SPINLOCK,
            bool>::type is_writing() const
        {
          while (spinlock.test_and_set(std::memory_order_acquire));
          bool res = written < this->size;
          spinlock.clear(std::memory_order_release);
          return res;
        }

        /// check if the buffer is blocked from further writes and no threads need to write to it.  uses spinlock
        template<bliss::concurrent::LockType LT = LockType>
        typename std::enable_if<LT == bliss::concurrent::LockType::SPINLOCK,
            bool>::type is_read_only() const
        {
          while (spinlock.test_and_set(std::memory_order_acquire));
          bool res = (reserved > Capacity) && (written >= this->size);
          spinlock.clear(std::memory_order_release);
          return res;
        }


      protected:

        /**
         * @brief       reserves a position in the buffer at which to insert the count number of bytes.
         * @details     increment only if we have room (protected by mutex or spinlock).
         *                if full, return -1.
         *                if becoming over-capacity, mark as full and return -2 to indicate it just became full.
         *                if becoming exactly full, mark as full and return the reservation
         *                else return the reservation.
         * @param       count of bytes to reserve
         * @return      the position in buffer where the reservation starts.
         */
        int64_t internal_reserve(const uint32_t count)
        {
          std::atomic_thread_fence(std::memory_order_acquire);
          int64_t curr = reserved;

          // this part reduces likelihood of accessing freed memory (a different thread tries to reserve to a freed region, i.e. )
          // does not prevent ABA problem (if the memory is reallocated.)
          if (curr >= Capacity)
          {
            curr = -1;
          }
          else
          {
            if (curr + count > Capacity)  // filled to past capacity, so no valid position is returned.
            {
              reserved = Capacity + 1 ;
              this->size = curr;

              curr = -2;

            }
            else if (curr + count == Capacity)  // filled to capacity, so valid position is returned.
            {
              reserved = Capacity + 1 ;
              this->size = Capacity;

            }  // else normal append
            else {
              reserved += count;

            }
            std::atomic_thread_fence(std::memory_order_release);
          }
          return curr;
        }




        /// marked write as completed. mutex version
        template<bliss::concurrent::LockType LT = LockType>
        typename std::enable_if<LT == bliss::concurrent::LockType::MUTEX, void>::type complete_write(
            const int count)
        {
          std::lock_guard<std::mutex> lock(this->mutex);
          std::atomic_thread_fence(std::memory_order_acquire);

          written += count;
          std::atomic_thread_fence(std::memory_order_release);

        }

        /// marked write as completed. spinlock version
        template<bliss::concurrent::LockType LT = LockType>
        typename std::enable_if<LT == bliss::concurrent::LockType::SPINLOCK,
            void>::type complete_write(const int count)
        {
          while (spinlock.test_and_set(std::memory_order_acquire));

          written += count;

          spinlock.clear(std::memory_order_release);
        }


      public:


        /**
         * @brief       reserves a position in the buffer at which to insert count number of bytes.
         * @details     increment only if we have room (so done via mutex).
         *              see documentation of "internal_reserve"
         * @param count
         * @return      pointer at which to insert.
         */
        template<bliss::concurrent::LockType LT = LockType>
        typename std::enable_if<LT == bliss::concurrent::LockType::MUTEX,
            int64_t>::type reserve(const uint32_t count)
        {
          std::lock_guard<std::mutex> lock(this->mutex);

          return internal_reserve(count);
        }
        /**
         * @brief       reserves a position in the buffer at which to insert count number of bytes.
         * @details     increment only if we have room (so done via spinlock).
         *              see documentation of "internal_reserve"
         * @param count
         * @return      pointer at which to insert.
         */
        template<bliss::concurrent::LockType LT = LockType>
        typename std::enable_if<LT == bliss::concurrent::LockType::SPINLOCK,
            int64_t>::type reserve(const uint32_t count)
        {

          while (spinlock.test_and_set(std::memory_order_acquire));

          int64_t curr = internal_reserve(count);

          spinlock.clear(std::memory_order_release);
          return curr;
        }


        /**
         * @brief  prevents future writes.  block can be turned on while there are writes in progress.  once block is on, reserve will fail.
         *
         * @details purpose of this method is to swap the size and reserved, so that reserved is at max+1, and size is the smallest of all locking thread's pointers
         *         user can specify a size at which which the block occurs.
         */
        template<bliss::concurrent::LockType LT = LockType>
        typename std::enable_if<LT == bliss::concurrent::LockType::MUTEX, void>::type block_writes(
            int64_t _curr = -1)
        {
          std::lock_guard<std::mutex> lock(this->mutex);

          // if supplied a desired size at which to lock, use it.
          std::atomic_thread_fence(std::memory_order_acquire);
          int64_t curr = (_curr >= 0 ? _curr : reserved);

          // set reserve so no further reservations occur.
          reserved = Capacity + 1;

          // choose smallest curr value amount all threads.
          if (this->size > curr)
          {
        	  this->size = curr;
          }
          std::atomic_thread_fence(std::memory_order_release);

        }
        /**
         * @brief  prevents future writes.  block can be turned on while there are writes in progress.  once block is on, reserve will fail.
         *
         * @details purpose of this method is to swap the size and reserved, so that reserved is at max+1, and size is the smallest of all locking thread's pointers
         *         user can specify a size at which which the block occurs.
         */
        template<bliss::concurrent::LockType LT = LockType>
        typename std::enable_if<LT == bliss::concurrent::LockType::SPINLOCK,
            void>::type block_writes(int64_t _curr = -1)
        {
          while (spinlock.test_and_set(std::memory_order_acquire));

          // if supplied a desired size at which to lock, use it.

          int64_t curr = (_curr >= 0 ? _curr : reserved);

          // set reserve so no further reservations occur.
          reserved = Capacity + 1;

          // choose smallest curr value amount all threads.
          if (this->size > curr)
          {
        	  this->size = curr;
          }

          spinlock.clear(std::memory_order_release);
        }


        /**
         * @brief  allow future writes.  block can be turned off while there are writes in progress.  once block is off, reserve may proceed
         *
         * @details purpose of this method is to swap the size and reserved, so that reserved is again below max, and size is the at max+1
         */
        template<bliss::concurrent::LockType LT = LockType>
        typename std::enable_if<LT == bliss::concurrent::LockType::MUTEX, void>::type unblock_writes()
        {
          std::lock_guard<std::mutex> lock(this->mutex);
          std::atomic_thread_fence(std::memory_order_acquire);

          if (reserved > Capacity)
            reserved = this->size;

          // then reset size to max
          this->size = Capacity + 1;
          std::atomic_thread_fence(std::memory_order_release);
        }

        /**
         * @brief  allow future writes.  block can be turned off while there are writes in progress.  once block is off, reserve may proceed
         *
         * @details purpose of this method is to swap the size and reserved, so that reserved is again below max, and size is the at max+1
         */
        template<bliss::concurrent::LockType LT = LockType>
        typename std::enable_if<LT == bliss::concurrent::LockType::SPINLOCK, void>::type unblock_writes()
        {
          while (spinlock.test_and_set(std::memory_order_acquire));

          if (reserved > Capacity)
            reserved = this->size;

          // then reset size to max
          this->size = Capacity + 1;

          spinlock.clear(std::memory_order_release);

        }



        /**
         * @brief   Sets the buffer as blocked and wait for all pending writes.
         */
        void block_and_flush()
        {
          block_writes<LockType>();
          while (is_writing<LockType>())  _mm_pause();
        }



        /**
         * @brief Clears the buffer. (set the size to 0, leaving the Capacity and the memory allocation intact).  also blocks future writes
         * @note  const because the caller will have a const reference to the buffer.
         */
        template<bliss::concurrent::LockType LT = LockType>
        typename std::enable_if<LT == bliss::concurrent::LockType::SPINLOCK,
            void>::type clear_and_block_writes()
        {
          while (spinlock.test_and_set(std::memory_order_acquire)) ;
          reserved = Capacity + 1;
          this->size = 0;
          written = 0;
          spinlock.clear(std::memory_order_release);
        }

        /**
         * @brief Clears the buffer. (set the size to 0, leaving the Capacity and the memory allocation intact).  also blocks future writes
         * @note  const because the caller will have a const reference to the buffer.
         */
        template<bliss::concurrent::LockType LT = LockType>
        typename std::enable_if<LT == bliss::concurrent::LockType::MUTEX,
            void>::type clear_and_block_writes()
        {
          std::lock_guard<std::mutex> lock(this->mutex);
          reserved = Capacity + 1;
          this->size = 0;
          written = 0;
          std::atomic_thread_fence(std::memory_order_release);

        }

        /**
         * @brief Clears the buffer. (set the size to 0, leaving the Capacity and the memory allocation intact).  also unblocks the buffer for writes
         * @note  const because the caller will have a const reference to the buffer.
         */
        template<bliss::concurrent::LockType LT = LockType>
        typename std::enable_if<LT == bliss::concurrent::LockType::SPINLOCK,
            void>::type clear_and_unblock_writes()
        {
          while (spinlock.test_and_set(std::memory_order_acquire));
          reserved = 0;
          this->size = Capacity + 1;
          written = 0;

          spinlock.clear(std::memory_order_release);

        }


        /**
         * @brief Clears the buffer. (set the size to 0, leaving the Capacity and the memory allocation intact).  also unblocks the buffer for writes
         * @note  const because the caller will have a const reference to the buffer.
         */
        template<bliss::concurrent::LockType LT = LockType>
        typename std::enable_if<LT == bliss::concurrent::LockType::MUTEX,
            void>::type clear_and_unblock_writes()
        {
          std::lock_guard<std::mutex> lock(this->mutex);
          reserved = 0;
          this->size = Capacity + 1;
          written = 0;
          std::atomic_thread_fence(std::memory_order_release);

        }


        /**
         * @brief   Append data to the buffer.  DEBUGGING VERSION
         * @details
         *          The function first reserves a block of memory,
         *          if successful, then memcpy the data into that block.
         *
         *          reserve is the part that is thread aware.
         *
         * method is const because the caller will have const reference to the Buffer.
         *
         * This method relies on the current reserved incrementing monotonically (by count) until it exceeds Capacity, and the written
         * incrementing monotonically (by count) as threads finish writes, chasing the reserved thread.
         *
         * The first thread to exceed the Capacity will set the size.  written will never exceed that, and reserved will not
         * be smaller than size.
         *
         *
         * @param[in] _data   pointer to data to be copied into the Buffer.  this SHOULD NOT be shared between threads
         * @param[in] count   number of bytes to be copied into the Buffer
         * @param[out] _inserted  parameter for debugging.  points to location of insertion.
         *
         * @return            2 bits:  LSB indicates if data was written.  MSB indicates if buffer is full and require swapping.
         */
        unsigned int append(const void* _data, const uint32_t count, void* &_inserted)  // _inserted is for DEBUGGING only.
        {

          if (count == 0 || _data == nullptr)
            throw std::invalid_argument("_data is nullptr");

          // reserve a spot.
          int64_t pos = reserve<LockType>(count);
          _inserted = nullptr;

          // if fails, then already full
          if (pos == -1)
          {
            // already read locked at this point

            return 0x0;  // full and not swapping.
          }
          else if (pos == -2)
          { // filled to beyond capacity

            while (is_writing<LockType>())  // wait for other writes to complete
              _mm_pause();
            return 0x2;  // full, swapping.
          }
          else
          { // valid position returned, so can write.

            // write
            std::memcpy(this->data_ptr + pos, _data, count);
            complete_write<LockType>(count); // all full buffers lock the read and unlock the writer

            // for DEBUGGING.
            _inserted = this->data_ptr + pos;

            if ((pos + count) == Capacity)
            { // thread that JUST filled the buffer

              // wait for other threads to finish
              while (is_writing<LockType>())
                _mm_pause();

              return 0x3;  // write to full, swapping, and success
            }
            else
            {
              return 0x1;   // not full, not swapping, successfully inserted.
            } // no case for pos + count > Capacity.

          }
        }

        /**
         * @brief   Append data to the buffer.  Wrapper for DEBUGGING version to present NON_DEBUG API to user.
         *
         * @param[in] _data   pointer to data to be copied into the Buffer.  this SHOULD NOT be shared between threads
         * @param[in] count   number of bytes to be copied into the Buffer
         *
         * @return            2 bits:  LSB indicates if data was written.  MSB indicates if buffer is full and require swapping.
         */
        unsigned int append(const void* _data, const uint32_t count)
        {
          void* output = nullptr;
          unsigned int result = append(_data, count, output);

//          // DEBUG ONLY.
//          if (result & 0x1) {
//            if (output == nullptr) {
//              ERRORF("ERROR: successful insert but result pointer is null.");
//            } else if (std::memcmp(_data, output, count) != 0) {
//              ERRORF("ERROR: A successful insert but input not same as what was memcpy'd.");
//            }
//          }

          return result;
        }
    };



    /**
     * @brief   NON-THREAD-SAFE Memory Buffer, which is a fixed size allocated block of memory where data can be appended to, and read from via pointer.
     * @details this class supports only MOVE semantics.
     *
     * Thread Safety is not enforced so it should be used only locally.  when write, reservation is first made and then memory is written.
     *
     * life cycle:  application acquires a buffer from objectpool.
     *              application threads write into buffer until full
     *              a single application thread swaps in a new, empty buffer from objectpool, and
     *                  dispatches the full buffer for communication
     *              communication thread completes and releases buffer back to object pool
     *
     * When application thread writes, there are 4 possible scenarios:
     *        1. reservation succeeds and is completely within the buffer.  write succeeds in this case.
     *        2. reservation fails as the target range is completely outside buffer.  write fails in this case.
     *        3. reservation succeeds and exactly fills the buffers.  write succeeds in this case, but an
     *            empty buffer is also swapped in and the full buffer dispatched for communication.
     *        4. reservation would span buffer maximum size limit.  write fails in this case, and an
     *            empty buffer is swapped in and the full buffer dispatched for communication.
     *
     * @note  because of the last 2 cases, the replacement of a buffer is done by a single thread only.
     *        however, because other threads may still have pointer to the full buffer, shared_ptr is a better choice.
     *
     *        symptoms of data races from swap include heap corruption, write after free errors as reported by ASAN.
     *        also, such swaps are subject to the ABA problem.
     *
     *        recommendation is to have each thread hold its own non-thread-safe buffer, or use CMPXCHG16B with reference counted pointers.
     *
     * @tparam  LockType  controls the type of thread safely.  See concurrent/concurrent.h for allowed values.
     * @tparam  Capacity  size of the allocated byte array
     */
    template<int64_t Capacity, size_t MetadataSize>
    class Buffer<bliss::concurrent::LockType::NONE, Capacity, MetadataSize> : public BufferBase<Capacity, MetadataSize>
    {

        template<bliss::concurrent::LockType LT, int64_t C, size_t MDS>
        friend std::ostream & operator<<(std::ostream &os, const Buffer<LT, C, MDS>& p);


      protected:
        using BaseType = BufferBase<Capacity, MetadataSize>;
      
        using BufferType = Buffer<bliss::concurrent::LockType::NONE, Capacity, MetadataSize>;
        // since not threaded, uses a single size member variable  and blocked flag.

        /// inidicate if the buffer is accepting further append
        mutable bool blocked;


      private:
        /**
         * @brief Move constructor with a mutex lock on the source object.
         * @details This constructor is private and only use as constructor delegation target.
         *  the source object is locked before this function is called and data moved to the newly constructed object.
         *
         * @param other   the source Buffer
         * @param l       the mutex lock on the source Buffer.
         */
        Buffer(BufferType && other, const std::lock_guard<std::mutex> &l)
            : BaseType(other, l),
              blocked(other.blocked)
        {
          other.blocked = Capacity + 1;
        };

        /// remove copy constructor and copy assignement operators.
        explicit Buffer(const BufferType& other) = delete;
        /// remove copy constructor and copy assignement operators.
        BufferType& operator=(const BufferType& other) = delete;


      public:
        /**
         * @brief Normal constructor.  Allocate and initialize memory with size specified as parameter.
         */
        Buffer() : BaseType(), blocked(true)
        {}
        //start_ptr(new uint8_t[Capacity](), []( uint8_t *p ) { delete [] p; })  for shared ptr

        /**
         * @brief Destructor.  waits for all writes and then deallocate memory manually.
         */
        virtual ~Buffer()
        {
          block_and_flush();
        }

        /**
         * @brief Move constructs from a Buffer with the SAME LockType property
         * @details
         * Thread Safe always, using approach from http://www.justsoftwaresolutions.co.uk/threading/thread-safe-copy-constructors.html
         *
         * Constructs a mutex lock and then delegates to another constructor.
         *
         * @param other   Source object to move
         */
        explicit Buffer(BufferType && other)
            : BufferType(std::move(other),
                               std::lock_guard<std::mutex>(other.mutex))
        {}

        /**
         * @brief Move assignment operator, between Buffers of the SAME LockType property.
         * @details  Internal data memory moved.
         * The move is done in a thread safe way always.
         *
         * @param other     Source Buffer to move
         * @return          target Buffer reference
         */
        BufferType& operator=(BufferType && other)
        {

          if (this->start_ptr != other.start_ptr)
          {

            std::unique_lock<std::mutex> myLock(this->mutex, std::defer_lock),
                otherLock(other.mutex, std::defer_lock);
            std::lock(myLock, otherLock);

            /// move the internal memory.
            blocked = other.blocked;

            this->move_assigner(std::move(other), myLock, otherLock);

            other.blocked = true;
          }
          return *this;
        }

      protected:

        /**
         * @brief get the current approximate size of the Buffer.
         * @return    current size
         */
        int64_t getApproximateSize() const
        {
          return this->size;
        }
        /**
         * @brief get the current written data size of the Buffer.
         * @note  written data may be scattered.
         * @return    current size
         */
        int64_t getWrittenSize() const
        {
          return this->size;
        }

        /// Check if buffer has been blocked from further reservation but still waiting for all writes to finish.
        bool is_flushing() const
        {
          return false;
        }



      public:

        /**
         * @brief Checks if a buffer is empty.
         *
         * @return    true if the buffer is empty, false otherwise.
         */
        bool isEmpty() const
        {
          return this->size == 0;
        }

        // TODO:  make this more clear.
        // need 1. flag to indicate no more new writers
        //      2. wait for current writers to finish
        //      3. 3 phases:  writing, flushing, reading.  -> reserve, wirte_unlock, flush_begin, block_writes, unblock_writes
        //          reserve (counter ++)
        //          wirte_unlock (counter --)
        //          flush (MSB = 1)
        //          block_writes ( wait for MSB==1 && counter == 0)
        //          unblock_writes ( MSB = 0)
        // state checks:  is_writing (MSB==0|1 && counter > 0)
        //                is_flushing (MSB==1 && counter> 0)
        //                is_read_only ( MSB==1 && counter == 0)

        /// check if there are threads writing to the buffer.
        bool is_writing() const
        {
          return false;
        }

        // check if the buffer is blocked from further writes and no threads need to write to it.
        bool is_read_only() const
        {
          return blocked;
        }


        /**
          * @brief       reserves a position in the buffer at which to insert count number of bytes.
          * @details     increment only if we have room
         *                if full, return -1.
         *                if becoming over-capacity, mark as full and return -2 to indicate it just became full.
         *                if becoming exactly full, mark as full and return the reservation
         *                else return the reservation.
          * @param count
          * @return      pointer at which to insert.
          */
        int64_t reserve(const uint32_t count)
        {
          int64_t curr = this->size;

          if (blocked)
          {
            curr = -1;
          }
          else
          {
            if (curr + count > Capacity)  // filled to past capacity, so no valid position is returned.
            {
              blocked = true;
              curr = -2;

            } else if (curr + count == Capacity)  // filled to capacity, so valid position is returned.
            {
              blocked = true ;
              this->size = Capacity;
            }  // else normal append
            else
              this->size += count;

          }
          return curr;
        }


        /**
         * @brief  prevents future writes.  block can be turned on while there are writes in progress.  once block is on, reserve will fail.
         *
         * @details purpose of this method is to swap the size and reserved, so that reserved is at max+1, and size is the smallest of all locking thread's pointers
         *         user can specify a size at which which the block occurs.
         */
        void block_writes(int64_t _curr = -1)
        {
          this->size = (_curr >= 0 ? _curr : this->size);
          blocked = true;
        }

        /**
         * @brief  allow future writes.  block can be turned off while there are writes in progress.  once block is off, reserve may proceed
         *
         * @details purpose of this method is to swap the size and reserved, so that reserved is again below max, and size is the at max+1
         */
        void unblock_writes()
        {
          blocked = false;
        }

        /**
         * @brief   Sets the buffer as read only and wait for all pending writes.
         */
        void block_and_flush()
        {
          block_writes();
        }



        /**
         * @brief Clears the buffer. (set the size to 0, leaving the Capacity and the memory allocation intact).  also blocks future writes
         * @note  const because the caller will have a const reference to the buffer.
         */
        void clear_and_block_writes()
        {
          // blocked
          blocked = true;
          this->size = 0;
        }

        /**
         * @brief Clears the buffer. (set the size to 0, leaving the Capacity and the memory allocation intact).  also unblocks the buffer for writes
         * @note  const because the caller will have a const reference to the buffer.
         */
        void clear_and_unblock_writes()
        {
          // blocked
          blocked = false;
          this->size = 0;
        }

        /**
         * @brief   Append data to the buffer.  DEBUGGING VERSION
         * @details
         *          The function first reserves a block of memory,
         *          if successful, then memcpy the data into that block.
         *
         *          reserve is the part that is thread aware.
         *
         * method is const because the caller will have const reference to the Buffer.
         *
         * This method relies on the current reserved incrementing monotonically (by count) until it exceeds Capacity, and the written
         * incrementing monotonically (by count) as threads finish writes, chasing the reserved thread.
         *
         * The first thread to exceed the Capacity will set the size.  written will never exceed that, and reserved will not
         * be smaller than size.
         *
         *
         * @param[in] _data   pointer to data to be copied into the Buffer.  this SHOULD NOT be shared between threads
         * @param[in] count   number of bytes to be copied into the Buffer
         * @param[out] _inserted  parameter for debugging.  points to location of insertion.
         *
         * @return            2 bits:  LSB indicates if data was written.  MSB indicates if buffer is full and require swapping.
         */
        unsigned int append(const void* _data, const uint32_t count, void* &_inserted)  // _inserted is for DEBUGGING only.
        {

          if (count == 0 || _data == nullptr)
            throw std::invalid_argument("_data is nullptr");

          // reserve a spot.
          int64_t pos = reserve(count);
          _inserted = nullptr;

          // if fails, then already full
          if (pos == -1)
          {
            // already read locked.
            return 0x0;  // full and not swapping.
          }
          else if (pos == -2)
          { // filled to beyond capacity

            // no other threads, no wait for other thread to finish writing.
            return 0x2;  // full, swapping.
          }
          else
          { // valid position returned, so can write.

            // write
            std::memcpy(this->data_ptr + pos, _data, count);

            // for DEBUGGING.
            _inserted = this->data_ptr + pos;

            if ((pos + count) == Capacity)
            { // thread that JUST filled the buffer

              // no other threads, so no wait for other threads to finish
              return 0x3;  // write to full, swapping, and success
            }
            else
            {
              return 0x1;   // not full, not swapping, successfully inserted.
            }  // no case for pos + count > Capacity.

          }
        }

        /**
         * @brief   Append data to the buffer.  Wrapper for DEBUGGING version to present NON_DEBUG API to user.
         *
         * @param[in] _data   pointer to data to be copied into the Buffer.  this SHOULD NOT be shared between threads
         * @param[in] count   number of bytes to be copied into the Buffer
         *
         * @return            2 bits:  LSB indicates if data was written.  MSB indicates if buffer is full and require swapping.
         */
        unsigned int append(const void* _data, const uint32_t count)  // _inserted is for DEBUGGING only.
        {
          void* output = nullptr;
          unsigned int result = append(_data, count, output);

//          // DEBUG ONLY
//          if (result & 0x1) {
//            if (output == nullptr) {
//              ERRORF("ERROR: successful insert but result pointer is null.");
//            } else if (std::memcmp(_data, output, count) != 0) {
//              ERRORF("ERROR: B successful insert but input not same as what was memcpy'd.");
//            }
//          }

          return result;
        }
    };




    /**
     * @brief   LockFree Memory Buffer, which is a fixed size allocated block of memory where data can be appended to, and read from via pointer.
     * @details this class supports only MOVE semantics.
     *
     * Thread Safety is enforced using atomic variable.  when write, reservation is first made and then memory is written.
     *
     * life cycle:  application acquires a buffer from objectpool.
     *              application threads write into buffer until full
     *              a single application thread swaps in a new, empty buffer from objectpool, and
     *                  dispatches the full buffer for communication
     *              communication thread completes and releases buffer back to object pool
     *
     * When application thread writes, there are 4 possible scenarios:
     *        1. reservation succeeds and is completely within the buffer.  write succeeds in this case.
     *        2. reservation fails as the target range is completely outside buffer.  write fails in this case.
     *        3. reservation succeeds and exactly fills the buffers.  write succeeds in this case, but an
     *            empty buffer is also swapped in and the full buffer dispatched for communication.
     *        4. reservation would span buffer maximum size limit.  write fails in this case, and an
     *            empty buffer is swapped in and the full buffer dispatched for communication.
     *
     * @note  because of the last 2 cases, the replacement of a buffer is done by a single thread only.
     *        however, because other threads may still have pointer to the full buffer, shared_ptr is a better choice.
     *
     *        symptoms of data races from swap include heap corruption, write after free errors as reported by ASAN.
     *        also, such swaps are subject to the ABA problem.
     *
     *        recommendation is to have each thread hold its own non-thread-safe buffer, or use CMPXCHG16B with reference counted pointers.
     *
     * @tparam  LockType  controls the type of thread safely.  See concurrent/concurrent.h for allowed values.
     * @tparam  Capacity  size of the allocated byte array
     */
    template<int64_t Capacity, size_t MetadataSize>
    class Buffer<bliss::concurrent::LockType::LOCKFREE, Capacity, MetadataSize> : public BufferBase<Capacity, MetadataSize>
    {
        static_assert(Capacity > 0, "Buffer Capacity is given as 0");

        template<bliss::concurrent::LockType LT, int64_t C, size_t MDS>
        friend std::ostream & operator<<(std::ostream &os, const Buffer<LT, C, MDS>& p);


      protected:
        using BufferType = Buffer<bliss::concurrent::LockType::LOCKFREE, Capacity, MetadataSize>;
        
        using BaseType = BufferBase<Capacity, MetadataSize>;

        // instead of a single "size" variable and a "blocked" variable (to prevent further writes to the buffer)
        //   (which require concurrent updates to both)
        // we use the sign bit to indicate "blocked".
        // also, because reservation and actual write may occur separately, we need a "reserved" and a "written"
        //   for separate incrementing.
        // having the "written" variable allows us to wait for all writes to complete.

        /// pointer to current head of reservation.
        /// capacity + 1 indicates buffer is blocked/full.  other values indicate that buffer is unblocked.
        /// exchanges with size.
        volatile std::atomic<int64_t> reserved;


        /// represent amount of data written.  will not update beyond the FINAL size
        volatile std::atomic<int64_t> written;

      private:
        /**
         * @brief Move constructor with a mutex lock on the source object.
         * @details This constructor is private and only use as constructor delegation target.
         *  the source object is locked before this function is called and data moved to the newly constructed object.
         *
         * @param other   the source Buffer
         * @param l       the mutex lock on the source Buffer.
         */
        Buffer(BufferType && other, const std::lock_guard<std::mutex> &l) :
              BaseType(other, l),
        {
			reserved.store(other.reserved.exchange(Capacity+1, std::memory_order_acquire), std::memory_order_release);
			written.store(other.written.exchange(0, std::memory_order_acquire), std::memory_order_release);
        };

        /// remove copy constructor and copy assignement operators.
        explicit Buffer(const BufferType& other) = delete;
        /// remove copy constructor and copy assignement operators.
        BufferType& operator=(const BufferType& other) = delete;


      public:
        /**
         * @brief Normal constructor.  Allocate and initialize memory with size specified as parameter.
         */
        Buffer() : BaseType(),
         reserved(Capacity + 1), written(0L)
        {}
        //start_ptr(new uint8_t[Capacity](), []( uint8_t *p ) { delete [] p; })  for shared ptr

        /**
         * @brief Destructor.  waits for all writes and then deallocate memory manually.
         */
        virtual ~Buffer()
        {
          block_and_flush();
        }

        /**
         * @brief Move constructs from a Buffer with the SAME LockType property
         * @details
         * Thread Safe always, using approach from http://www.justsoftwaresolutions.co.uk/threading/thread-safe-copy-constructors.html
         *
         * Constructs a mutex lock and then delegates to another constructor.
         *
         * @param other   Source object to move
         */
        explicit Buffer(BufferType && other)
            : BufferType(std::move(other),
                               std::lock_guard<std::mutex>(other.mutex))
        {}

        /**
         * @brief Move assignment operator, between Buffers of the SAME LockType property.
         * @details  Internal data memory moved by std::unique_ptr semantics.
         * The move is done in a thread safe way always.
         *
         * @param other     Source Buffer to move
         * @return          target Buffer reference
         */
        BufferType& operator=(BufferType && other)
        {

          if (this->start_ptr != other.start_ptr)
          {

            std::unique_lock<std::mutex> myLock(this->mutex, std::defer_lock),
                otherLock(other.mutex, std::defer_lock);
            std::lock(myLock, otherLock);

            this->move_assigner(std::move(other), myLock, otherLock);

            reserved.store(other.reserved.exchange(Capacity+1, std::memory_order_acquire), std::memory_order_release);
            written.store(other.written.exchange(0, std::memory_order_acquire), std::memory_order_release);
          }
          return *this;
        }

      protected:

        /**
         * @brief get the current approximate size of the Buffer.
         * @return    current size
         */
        int64_t getApproximateSize() const
        {
          return reserved.load(std::memory_order_relaxed);
        }

        /**
         * @brief get the current written data size of the Buffer.
         * @note  written data may be scattered.
         * @return    current size
         */
        int64_t getWrittenSize() const
        {
          return written.load(std::memory_order_relaxed);
        }

        /// Check if buffer has been blocked from further reservation but still waiting for all writes to finish.
        bool is_flushing() const
        {

        	std::atomic_thread_fence(std::memory_order_acquire);
        	auto s = this->size;
          return (reserved.load(std::memory_order_acquire) > Capacity) && (s > written.load(std::memory_order_relaxed));
        }



      public:

        /**
         * @brief Checks if a buffer is empty.
         * @note The return value is approximate due to threading
         *
         * @return    true if the buffer is empty, false otherwise.
         */
        bool isEmpty() const
        {
          return written.load(std::memory_order_relaxed) == 0;
        }

        // TODO:  make this more clear.
        // need 1. flag to indicate no more new writers
        //      2. wait for current writers to finish
        //      3. 3 phases:  writing, flushing, reading.  -> reserve, wirte_unlock, flush_begin, block_writes, unblock_writes
        //          reserve (counter ++)
        //          wirte_unlock (counter --)
        //          flush (MSB = 1)
        //          block_writes ( wait for MSB==1 && counter == 0)
        //          unblock_writes ( MSB = 0)
        // state checks:  is_writing (MSB==0|1 && counter > 0)
        //                is_flushing (MSB==1 && counter> 0)
        //                is_read_only ( MSB==1 && counter == 0)

        /// check if there are threads writing to the buffer.
        bool is_writing() const
        {
          // is writing - really only for debugging external to this class
          //int64_t s = size.load(std::memory_order_relaxed);
          return written.load(std::memory_order_relaxed) < this->size;
        }

        /// check if the buffer is blocked from further writes and no threads need to write to it.
        bool is_read_only() const
        {
          // is read only - really only for debugging external to this class
        	std::atomic_thread_fence(std::memory_order_acquire);
        	auto s = this->size;
          return (reserved.load(std::memory_order_acquire) > Capacity) && (s <= written.load(std::memory_order_relaxed));
        }


        /**
         * @brief       reserves a position in the buffer at which to insert the count number of bytes.
         * @details     increment only if we have room (protected by mutex or spinlock).
         *                if full, return -1.
         *                if becoming over-capacity, mark as full and return -2 to indicate it just became full.
         *                if becoming exactly full, mark as full and return the reservation
         *                else return the reservation.
         * @param       count of bytes to reserve
         * @return      the position in buffer where the reservation starts.
         */
        int64_t reserve(const uint32_t count)
        {

          // blocked buffer, from full or explicit block call.
          if (reserved.load(std::memory_order_acquire) > Capacity) {
            // have this here to minimize fetch_add calls.

            return -1;
          } else {
            // not yet full, so reserve
            int64_t curr = reserved.fetch_add(count, std::memory_order_acq_rel);

            if (curr >= Capacity) {
              return -1;
            }
            else if (curr + count > Capacity) {
              reserved.store(Capacity + 1, std::memory_order_relaxed);
              // just filled. curr position not a valid insertion point.

              // since only 1 thread reaching here, just set size, no need to check to make sure we store minimum.
              //size.store(curr, std::memory_order_release);
              this->size = curr;
              std::atomic_thread_fence(std::memory_order_release);
              return -2;
            } else if (curr + count == Capacity) {
              reserved.store(Capacity + 1, std::memory_order_relaxed);
              // just filled.  curr position is a valid insertion point.

              // since only 1 thread reaching here, just set size, no need to check to make sure we store minimum.
              //size.store(Capacity, std::memory_order_release);
              this->size = Capacity;
              std::atomic_thread_fence(std::memory_order_release);
              return curr;
            } else return curr;
          }
        }

      protected:

        /// marked write as completed.
        void complete_write(const int count)
        {
          written.fetch_add(count, std::memory_order_acq_rel);
        }

      public:

        /**
         * @brief  prevents future writes.  block can be turned on while there are writes in progress.  once block is on, reserve will fail.
         *
         * @details purpose of this method is to swap the size and reserved, so that reserved is at max+1, and size is the smallest of all locking thread's pointers
         *         user can specify a size at which which the block occurs.
         *
         *         threading cases:  reserve <= capacity;  reserve = capacity+1;  reserver > capacity
         *         					size = capacity + 1; size <= capacity.
         *
         *         					size > reserve:  reserve <= cap, size = cap + 1:  size = reserve
         *         									 reserve < size <= cap:	  		  size = reserve.
         *         					size <= reserve: size < reserve <= cap:			  size = size
         *         		reserve.exchange guarantees that only 1 thread will be able to match condition above.
         *         		want minimum of reserve and size to set as new size.
         */
        void block_writes(int64_t _curr = -1)
        {
          int64_t curr = reserved.exchange(Capacity + 1, std::memory_order_acq_rel);
          // curr could be greater than capacity, or not.

          //int64_t end = size.load(std::memory_order_relaxed);  // size initially at Capacity + 1

          curr = (_curr >= 0 ? _curr : curr);
          //bool stop = false;

          // multiple threads:  smallest wins.
//          while (end > curr && !stop)
//          {
//            stop = size.compare_exchange_weak(end, curr, std::memory_order_relaxed);
//          }
          this->size = std::min(curr, this->size);
          std::atomic_thread_fence(std::memory_order_release);

        }


        /**
         * @brief  allow future writes.  block can be turned off while there are writes in progress.  once block is off, reserve may proceed
         *
         * @details purpose of this method is to swap the size and reserved, so that reserved is again below max, and size is the at max+1
         * 			compare_exchange_strong ensures only 1 thread succeeds with this function.
         *
         * 			thread safety is not guaranteed if unblock_writes and block_writes are called concurrently.
         */
        void unblock_writes()
        {
//          // put the size into curr if it hasn't been done.
//          int64_t curr = reserved.load(std::memory_order_relaxed);
//          int64_t end = size.exchange(Capacity + 1, std::memory_order_relaxed);
//
//          bool stop = false;
//          while (curr > Capacity && !stop)
//          {
//            stop = reserved.compare_exchange_weak(curr, end, std::memory_order_relaxed);
//          }

        	int64_t curr = Capacity + 1;
        	if (reserved.compare_exchange_strong(curr, this->size, std::memory_order_acq_rel)) {
        		this->size = Capacity + 1;
                std::atomic_thread_fence(std::memory_order_release);
        	} else {
        		ERRORF("unblock_writes called, but reserved is not at Capacity+1. reserved: %ld, size: %ld", curr, this->size);
        	}
        }


        /**
         * @brief   Sets the buffer as read only and wait for all pending writes.
         */
        void block_and_flush()
        {
          block_writes();
          while (is_writing()) {
            _mm_pause();
          }
        }

        /**
         * @brief Clears the buffer. (set the size to 0, leaving the Capacity and the memory allocation intact).  also blocks future writes
         * @note  const because the caller will have a const reference to the buffer.
         */
        void clear_and_block_writes()
        {
          std::atomic_thread_fence(std::memory_order_acquire);
          // blocked
          reserved.store(Capacity + 1, std::memory_order_relaxed);
          written.store(0, std::memory_order_relaxed);
          this->size = 0;
          std::atomic_thread_fence(std::memory_order_release);

        }


        /**
         * @brief Clears the buffer. (set the size to 0, leaving the Capacity and the memory allocation intact).  also unblocks the buffer for writes
         * @note  const because the caller will have a const reference to the buffer.
         */
        void clear_and_unblock_writes()
        {
          std::atomic_thread_fence(std::memory_order_acquire);
          written.store(0, std::memory_order_relaxed);
          reserved.store(0, std::memory_order_relaxed);
          this->size = Capacity + 1;
          std::atomic_thread_fence(std::memory_order_release);
        }

        /**
         * @brief   Append data to the buffer.  DEBUGGING VERSION
         * @details
         *          The function first reserves a block of memory,
         *          if successful, then memcpy the data into that block.
         *
         *          reserve is the part that is thread aware.
         *
         * method is const because the caller will have const reference to the Buffer.
         *
         * This method relies on the current reserved incrementing monotonically (by count) until it exceeds Capacity, and the written
         * incrementing monotonically (by count) as threads finish writes, chasing the reserved thread.
         *
         * The first thread to exceed the Capacity will set the size.  written will never exceed that, and reserved will not
         * be smaller than size.
         *
         *
         * @param[in] _data   pointer to data to be copied into the Buffer.  this SHOULD NOT be shared between threads
         * @param[in] count   number of bytes to be copied into the Buffer
         * @param[out] _inserted  parameter for debugging.  points to location of insertion.
         *
         * @return            2 bits:  LSB indicates if data was written.  MSB indicates if buffer is full and require swapping.
         */
        unsigned int append(const void* _data, const uint32_t count, void* &_inserted)  // _inserted is for DEBUGGING only.
        {

          if (count == 0 || _data == nullptr)
            throw std::invalid_argument("_data is nullptr");

          // reserve a spot.
          int64_t pos = reserve(count);
          _inserted = nullptr;

          // if fails, then already full
          if (pos == -1)
          {
            // already read locked.

            return 0x0;  // full and not swapping.
          }
          else if (pos == -2)
          { // filled to beyond capacity

            while (is_writing())  // wait for other writes to complete
              _mm_pause();
            return 0x2;  // full, swapping.
          }
          else
          { // valid position returned, so can write.

            // write
            std::memcpy(this->data_ptr + pos, _data, count);
            complete_write(count); // all full buffers lock the read and unlock the writer

            // for DEBUGGING.
            _inserted = this->data_ptr + pos;

            if ((pos + count) == Capacity)
            { // thread that JUST filled the buffer

              // wait for other threads to finish
              while (is_writing())
                _mm_pause();

              return 0x3;  // write to full, swapping, and success
            }
            else
            { // not filled.

              return 0x1;   // not full, not swapping, successfully inserted.
            } // no case for pos + count > Capacity.

          }
        }

        /**
         * @brief   Append data to the buffer.  Wrapper for DEBUGGING version to present NON_DEBUG API to user.
         *
         * @param[in] _data   pointer to data to be copied into the Buffer.  this SHOULD NOT be shared between threads
         * @param[in] count   number of bytes to be copied into the Buffer
         *
         * @return            2 bits:  LSB indicates if data was written.  MSB indicates if buffer is full and require swapping.
         */
        unsigned int append(const void* _data, const uint32_t count)  // _inserted is for DEBUGGING only.
        {
          void* output = nullptr;
          //uint8_t* ptr = this->data_ptr;
          unsigned int result = append(_data, count, output);

          // DEBUG ONLY
//          if (result & 0x1) {
//            if (output == nullptr) {
//              ERRORF("ERROR: successful insert but result pointer is null.");
//            } else if (std::memcmp(_data, output, count) != 0) {
//              ERRORF("ERROR: C successful insert but input not same as what was memcpy'd.");
//            }
//          }

          return result;
        }
    };





    /**
     * @brief << operator to write out DataBlock object's actual data to ostream.
     * @tparam Iterator   Source data iterator type.
     * @tparam Range      Range data type
     * @tparam Container  container type for buffer.  defaults to std::vector.
     * @param[in/out] ost   output stream to which the content is directed.
     * @param[in]     db    BufferedDataBlock object to write out
     * @return              output stream object
     */
    template<bliss::concurrent::LockType LT, int64_t C = 8192, size_t MDS = 0>
    std::ostream& operator<<(std::ostream& ost, const Buffer<LT, C, MDS> & buffer)
    {
      ost << "LockType=" << static_cast<int>(LT)
      << " BUFFER: this->data_ptr/size="
      << static_cast<const void *>(buffer.start_ptr) << "/"
      << buffer.getSize() << " written="
      << buffer.getWrittenSize() << " approx,size/cap="
      << buffer.getApproximateSize() << "," << buffer.getSize() << "/" << buffer.getMetadataSize() << "+"
      << buffer.getCapacity() << " R? " << (buffer.is_read_only() ? "y" : "n")
      << " F? " << (buffer.is_flushing() ? "y" : "n") << " W? "
      << (buffer.is_writing() ? "y" : "n") << std::flush;

      return ost;
    }


  } /* namespace io */
} /* namespace bliss */

#endif /* BUFFER_HPP_ */
