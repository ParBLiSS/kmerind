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
     * life cycle:  pool acquire():  buffer unblock_writes().  -> allow rw.
     *              application:  buffer block_writes() -> allow ro
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
     * data array access via shared_ptr
     *
     * @tparam  LockType  controls whether this class is thread safe or not
     * @tparam  Capacity  size of the allocated byte array
     */
    template<bliss::concurrent::LockType LockType, int64_t Capacity = 8192>
    class Buffer
    {

        static_assert(Capacity > 0, "Buffer Capacity is given as 0");

        /*
         * Declare Buffer<!LockType> class as a friend class, so we can reference its member functions and variables
         * in move constructor and assignment operators that moves data between instances with different thread safety.
         *
         * alternative is to use inheritance, which would require virtual functions that are potentially expensive.
         */
        template<bliss::concurrent::LockType LT, int64_t C>
        friend std::ostream & operator<<(std::ostream &os, const Buffer<LT, C>& p);



      protected:

        /// internal data storage
        mutable uint8_t* start_ptr; // const, does not change


        /// pointer to current head of reservation.
        /// capacity + 1 indicates buffer is blocked/full.  other values indicate that buffer is unblocked.
        /// exchanges with size.
        volatile int64_t reserved;

        /// pointer to FINAL  end of data.   only updated when buffer is blocked or when full (from reserved)
        /// capacity + 1 indicates buffer is unblocked.  other values indicate that buffer is blocked or full.
        ///  exchanges with reserved
        mutable int64_t size;

        /// represent amount of data written.  will not update beyond the FINAL size
        volatile int64_t written;

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
        Buffer(Buffer<LockType, Capacity> && other, const std::lock_guard<std::mutex> &l)
            : start_ptr(other.start_ptr),
              reserved(other.reserved), size(other.size),
              written(other.written)
        {

          other.start_ptr = nullptr;
          other.size = 0;
          other.written = 0;
          other.reserved = Capacity + 1;

        };

        /// remove copy constructor and copy assignement operators.
        explicit Buffer(const Buffer<LockType, Capacity>& other) = delete;
        /// remove copy constructor and copy assignement operators.
        Buffer<LockType, Capacity>& operator=(const Buffer<LockType, Capacity>& other) = delete;


      public:
        /**
         * @brief Normal constructor.  Allocate and initialize memory with size specified as parameter.
         */
        //start_ptr(new uint8_t[Capacity](), []( uint8_t *p ) { delete [] p; })  for shared ptr
        Buffer() : start_ptr(new uint8_t[Capacity]()),
         reserved(Capacity + 1), size(0L), written(0L)
        {}

        /**
         * @brief Destructor.  waits for all writes and then deallocate memory manually.
         */
        virtual ~Buffer()
        {
          block_and_flush();
        }

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
        explicit Buffer(Buffer<LockType, Capacity> && other)
            : Buffer<LockType, Capacity>(std::move(other),
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
        Buffer<LockType, Capacity>& operator=(Buffer<LockType, Capacity> && other)
        {

          if (this->start_ptr != other.start_ptr)
          {

            std::unique_lock<std::mutex> myLock(mutex, std::defer_lock),
                otherLock(other.mutex, std::defer_lock);
            std::lock(myLock, otherLock);

//            if (start_ptr != nullptr)
//            {
//              delete[] start_ptr;
//              start_ptr = nullptr;
//            }

            /// move the internal memory.

            start_ptr = other.start_ptr;
            written = (int64_t)(other.written);
            size = (int64_t)(other.size);
            reserved = (int64_t)(other.reserved);

            other.size = 0;
            other.start_ptr = nullptr;
            other.written = 0;
            other.reserved = Capacity + 1;
          }
          return *this;
        }

        int64_t getSize() const
        {
          return size;
        }

      protected:

        /**
         * @brief get the current approximate size of the Buffer.
         * @return    current size
         */
// need to fine tune the memory order so this method may not be useful
        int64_t getApproximateSize() const
        {
          return reserved;
        }
        /**
         * @brief get the current written data size of the Buffer.
         * @note  written data may be scattered.
         * @return    current size
         */
// need to fine tune the memory order so this method may not be useful
        int64_t getWrittenSize() const
        {
          return written;
        }


        template<bliss::concurrent::LockType LT = LockType>
        typename std::enable_if<LT == bliss::concurrent::LockType::MUTEX, bool>::type is_flushing() const
        {
          std::lock_guard<std::mutex> lock(mutex);
          return is_flushing<bliss::concurrent::LockType::NONE>();
        }

        template<bliss::concurrent::LockType LT = LockType>
        typename std::enable_if<LT == bliss::concurrent::LockType::SPINLOCK,
            bool>::type is_flushing() const
        {
          while (spinlock.test_and_set())
            ;
          bool res = is_flushing<bliss::concurrent::LockType::NONE>();
          spinlock.clear();
          return res;
        }

        template<bliss::concurrent::LockType LT = LockType>
        typename std::enable_if<LT == bliss::concurrent::LockType::NONE, bool>::type is_flushing() const
        {
          return (reserved > Capacity) && (written < size);
        }



      public:

        /**
         * @brief Get a pointer to the buffer data memory block.
         *
         * @note the data should not be deleted by a calling function/thread.
         * The access is read only.  There is no reason to return the unique_ptr.
         *
         * const because the caller will have a const reference to the buffer
         *
         * return weak pointer for use in multithreaded environment.  This is useful when
         * Buffer is not changing, but the data array is directly accessed by multiple threads.
         */
        uint8_t* getData() const
        {
          return start_ptr;
        }

        /**
         * @brief return a casted pointer to the internal data array.
         * @note  This is better suited when the Buffer itself is shared via shared_ptr and each thread has read access only.
         * @note  this breaks the encapsulation a little.
         */
        template<typename T>
        operator T*() const {
          return reinterpret_cast<T*>(start_ptr);
        }

        /**
         * @brief get the capacity of the buffer.
         * @return    maximum capacity of buffer.
         */
        const int64_t getCapacity() const
        {
          return Capacity;
        }

        /**
         * @brief Checks if a buffer is empty.
         * @note The return value is approximate due to threading
         *
         * @return    true if the buffer is empty, false otherwise.
         */
        // is this used?
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

        template<bliss::concurrent::LockType LT = LockType>
        typename std::enable_if<LT == bliss::concurrent::LockType::MUTEX, bool>::type is_writing() const
        {
          std::lock_guard<std::mutex> lock(mutex);
          return is_writing<bliss::concurrent::LockType::NONE>();
        }

        template<bliss::concurrent::LockType LT = LockType>
        typename std::enable_if<LT == bliss::concurrent::LockType::MUTEX, bool>::type is_read_only() const
        {
          std::lock_guard<std::mutex> lock(mutex);
          return is_read_only<bliss::concurrent::LockType::NONE>();
        }

        template<bliss::concurrent::LockType LT = LockType>
        typename std::enable_if<LT == bliss::concurrent::LockType::SPINLOCK,
            bool>::type is_writing() const
        {
          while (spinlock.test_and_set())
            ;
          bool res = is_writing<bliss::concurrent::LockType::NONE>();
          spinlock.clear();
          return res;
        }

        template<bliss::concurrent::LockType LT = LockType>
        typename std::enable_if<LT == bliss::concurrent::LockType::SPINLOCK,
            bool>::type is_read_only() const
        {
          while (spinlock.test_and_set())
            ;
          bool res = is_read_only<bliss::concurrent::LockType::NONE>();
          spinlock.clear();
          return res;
        }

        template<bliss::concurrent::LockType LT = LockType>
        typename std::enable_if<LT == bliss::concurrent::LockType::NONE, bool>::type is_writing() const
        {
          return written < size;
        }

        template<bliss::concurrent::LockType LT = LockType>
        typename std::enable_if<LT == bliss::concurrent::LockType::NONE, bool>::type is_read_only() const
        {
          return (reserved > Capacity) && (written >= size);
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



      protected:

        /**
         * @brief       reserves a position in the buffer at which to insert the count number of bytes.
         * @details     increment only if we have room (so done via atomic CAS.).  if full, return nullptr. if becoming full, set the size.  else return the insertion point
         * @param count
         * @return      pointer at which to insert.
         */
        template<bliss::concurrent::LockType LT = LockType>
        typename std::enable_if<LT == bliss::concurrent::LockType::MUTEX,
            int64_t>::type reserve(const uint32_t count)
        {

          std::unique_lock<std::mutex> lock(mutex);

          int64_t curr = reserve<bliss::concurrent::LockType::NONE>(count);

          lock.unlock();
          return curr;

        }
        /**
         * @brief       reserves a position in the buffer at which to insert the count number of bytes.
         * @details     increment only if we have room (so done via atomic CAS.).  if full, return nullptr. if becoming full, set the size.  else return the insertion point
         * @param count
         * @return      pointer at which to insert.
         */
        template<bliss::concurrent::LockType LT = LockType>
        typename std::enable_if<LT == bliss::concurrent::LockType::SPINLOCK,
            int64_t>::type reserve(const uint32_t count)
        {

          while (spinlock.test_and_set());

          int64_t curr = reserve<bliss::concurrent::LockType::NONE>(count);

          spinlock.clear();
          return curr;
        }


        /**
         * @brief       reserves a position in the buffer at which to insert the count number of bytes.
         * @details     increment only if we have room (so done via atomic CAS.).  if full, return nullptr. if becoming full, set the size.  else return the insertion point
         * @param count
         * @return      pointer at which to insert.
         */
        template<bliss::concurrent::LockType LT = LockType>
        typename std::enable_if<LT == bliss::concurrent::LockType::NONE, int64_t>::type reserve(
            const uint32_t count)
        {
          int64_t curr = reserved;

          // this part reduces likelihood of accessing freed memory (a different thread tries to reserve to a freed region, i.e. )
          // does not prevent ABA problem (if the memory is reallocated.)
          if (curr >= Capacity)
          {
            curr = -1;
          }
          else
          { //if ptr == Capacity, buffer must be FILLED in this cycle.  need to process it as if not full yet.

            if (curr + count > Capacity)  // filled to past capacity, so no valid position is returned.
            {
              reserved = Capacity + 1 ;
              size = curr;
              curr = -2;

            } else if (curr + count == Capacity)  // filled to capacity, so valid position is returned.
            {
              reserved = Capacity + 1 ;
              size = Capacity;
            }  // else normal append
            else
              reserved += count;

          }
          return curr;
        }

        /// marked write as completed.
        template<bliss::concurrent::LockType LT = LockType>
        typename std::enable_if<LT == bliss::concurrent::LockType::MUTEX, void>::type complete_write(
            const int count)
        {
          std::lock_guard<std::mutex> lock(mutex);

          complete_write<bliss::concurrent::LockType::NONE>(count);
        }

        /// marked write as completed.
        template<bliss::concurrent::LockType LT = LockType>
        typename std::enable_if<LT == bliss::concurrent::LockType::SPINLOCK,
            void>::type complete_write(const int count)
        {
          while (spinlock.test_and_set())
            ;

          complete_write<bliss::concurrent::LockType::NONE>(count);
          spinlock.clear();
        }


        template<bliss::concurrent::LockType LT = LockType>
        typename std::enable_if<LT == bliss::concurrent::LockType::NONE, void>::type complete_write(
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
        template<bliss::concurrent::LockType LT = LockType>
        typename std::enable_if<LT == bliss::concurrent::LockType::MUTEX, void>::type block_writes(
            int64_t _curr = -1)
        { // not using xor since it toggles and complete_write is not going to change the sign bit (so we don't need xor)
          // disable reserved and set the size
          std::lock_guard<std::mutex> lock(mutex);

          block_writes<bliss::concurrent::LockType::NONE>(_curr);
        }
        //===== read lock trumps write.
        /**
         * @brief  read lock prevents future writes.  read lock can be turned on while there are writes in progress.  once read lock is on, reserve will returns nullptr.
         *
         * @details purpose of this method is to swap the size and reserved, so that reserved is at max+1, and size is the smallest of all locking thread's pointers
         */
        template<bliss::concurrent::LockType LT = LockType>
        typename std::enable_if<LT == bliss::concurrent::LockType::SPINLOCK,
            void>::type block_writes(int64_t _curr = -1)
        { // not using xor since it toggles and complete_write is not going to change the sign bit (so we don't need xor)
          // disable reserved and set the size
          while (spinlock.test_and_set())
            ;
          block_writes<bliss::concurrent::LockType::NONE>(_curr);
          spinlock.clear();
        }

//        //===== read lock trumps write.
        template<bliss::concurrent::LockType LT = LockType>
        typename std::enable_if<LT == bliss::concurrent::LockType::NONE, void>::type block_writes(
            int64_t _curr = -1)
        {
          // disable reserved and set the size
          int64_t curr = (_curr >= 0 ? _curr : reserved);
          reserved = Capacity + 1;

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
        template<bliss::concurrent::LockType LT = LockType>
        typename std::enable_if<LT == bliss::concurrent::LockType::MUTEX, void>::type unblock_writes()
        {
          std::lock_guard<std::mutex> lock(mutex);

          unblock_writes<bliss::concurrent::LockType::NONE>();

        }

        /**
         * @brief  read lock prevents future writes.  read lock can be turned on while there are writes in progress.  once read lock is on, reserve will returns nullptr.
         *
         * @details purpose of this method is to swap the size and reserved, so that reserved is at max+1, and size is the smallest of all locking thread's pointers
         */
        template<bliss::concurrent::LockType LT = LockType>
        typename std::enable_if<LT == bliss::concurrent::LockType::SPINLOCK,
            void>::type unblock_writes()
        {
          while (spinlock.test_and_set())
            ;

          unblock_writes<bliss::concurrent::LockType::NONE>();
          spinlock.clear();

        }



        template<bliss::concurrent::LockType LT = LockType>
        typename std::enable_if<LT == bliss::concurrent::LockType::NONE, void>::type unblock_writes()
        {
          // put the size into curr if it hasn't been done.
          if (reserved > Capacity)
            reserved = size;
          // then reset size to max
          size = Capacity + 1;
        }

        /**
         * @brief   Sets the buffer as read only and wait for all pending writes.
         */
        void block_and_flush()
        {
          //printf("INFO:  set block start:  reserved %ld, size %ld, written %ld \n", (int64_t)reserved, (int64_t)size, (int64_t)written);
          block_writes<LockType>();
          //printf("INFO:  set block end:  reserved %ld, size %ld, written %ld \n", (int64_t)reserved, (int64_t)size, (int64_t)written);
          while (is_writing<LockType>()) {
            //printf("INFO:  wait for writes:  reserved %ld, size %ld, written %ld \n", (int64_t)reserved, (int64_t)size, (int64_t)written);
            _mm_pause();
          }
        }



        /**
         * @brief Clears the buffer. (set the size to 0, leaving the Capacity and the memory allocation intact)
         * @note  const because the caller will have a const reference to the buffer.
         */
        template<bliss::concurrent::LockType LT = LockType>
        typename std::enable_if<LT == bliss::concurrent::LockType::SPINLOCK,
            void>::type clear_and_block_writes()
        {
          while (spinlock.test_and_set())
            ;
          clear_and_block_writes<bliss::concurrent::LockType::NONE>();
          spinlock.clear();
        }

        template<bliss::concurrent::LockType LT = LockType>
        typename std::enable_if<LT == bliss::concurrent::LockType::MUTEX,
            void>::type clear_and_block_writes()
        {
          std::lock_guard<std::mutex> lock(mutex);
          // blocked
          reserved = Capacity + 1;
          size = 0;
          written = 0;
        }
        template<bliss::concurrent::LockType LT = LockType>
        typename std::enable_if<LT == bliss::concurrent::LockType::NONE,
            void>::type clear_and_block_writes()
        {
          // blocked
          reserved = Capacity + 1;
          size = 0;
          written = 0;
        }

        /**
         * @brief Clears the buffer. (set the size to 0, leaving the Capacity and the memory allocation intact)
         * @note  const because the caller will have a const reference to the buffer.
         */
        template<bliss::concurrent::LockType LT = LockType>
        typename std::enable_if<LT == bliss::concurrent::LockType::SPINLOCK,
            void>::type clear_and_unblock_writes()
        {
          while (spinlock.test_and_set())
            ;
          clear_and_unblock_writes<bliss::concurrent::LockType::NONE>();
          spinlock.clear();
        }


        template<bliss::concurrent::LockType LT = LockType>
        typename std::enable_if<LT == bliss::concurrent::LockType::MUTEX,
            void>::type clear_and_unblock_writes()
        {
          std::lock_guard<std::mutex> lock(mutex);
          // blocked
          reserved = 0;
          size = Capacity + 1;
          written = 0;
        }

        template<bliss::concurrent::LockType LT = LockType>
        typename std::enable_if<LT == bliss::concurrent::LockType::NONE,
            void>::type clear_and_unblock_writes()
        {
          // blocked
          reserved = 0;
          size = Capacity + 1;
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
         * This method relies on the current reserved incrementing monotonically (by count) until it exceeds Capacity, and the written
         * incrementing monotonically (by count) as threads finish writes, chasing the reserved thread.
         *
         * The first thread to exceed the Capacity will set the size.  written will never exceed that, and reserved will not
         * be smaller than size.
         *
         *  single primitive type is better as it can be returned in 1 op.
         *
         * @tparam LT Choose thread safe vs unsafe implementation. defaults to same as the parent class.
         * @param[in] _data   pointer to data to be copied into the Buffer..  this SHOULD NOT be shared between threads
         * @param[in] count   number of bytes to be copied into the Buffer
         * @return            unsigned char, bit 0 indicated whether operation succeeded, bit 1 indicating whether buffer swap is needed.
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
          { // write lock returns nullptr if ptr points to outside of the range.
            // already read locked.

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
            std::memcpy(start_ptr + pos, _data, count);
            complete_write<LockType>(count); // all full buffers lock the read and unlock the writer

            // for DEBUGGING.
            _inserted = start_ptr + pos;

            if ((pos + count) == Capacity)
            { // thread that JUST filled the buffer

              // wait for other threads to finish
              while (is_writing<LockType>())
                _mm_pause();

              return 0x3;  // write to full, swapping, and success
            }
            else  // no case for pos + count > Capacity.
            { // not filled.

              return 0x1;   // not full, not swapping, successfully inserted.
            }

          }
        }

/// DEBUG WRAPPER
        unsigned int append(const void* _data, const uint32_t count)  // _inserted is for DEBUGGING only.
        {
          void* output = nullptr;
          unsigned int result = append(_data, count, output);

          if (result & 0x1) {
            if (output == nullptr) {
              ERRORF("ERROR: successful insert but result pointer is null.");
            } else if (std::memcmp(_data, output, count) != 0) {
              ERRORF("ERROR: successful insert but input not same as what was memcpy'd.");
            }
          }

          return result;
        }
    };



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
     * life cycle:  pool acquire():  buffer unblock_writes().  -> allow rw.
     *              application:  buffer block_writes() -> allow ro
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
     * data array access via shared_ptr
     *
     * @tparam  LockType  controls whether this class is thread safe or not
     * @tparam  Capacity  size of the allocated byte array
     */
    template<int64_t Capacity>
    class Buffer<bliss::concurrent::LockType::NONE, Capacity>
    {

        static_assert(Capacity > 0, "Buffer Capacity is given as 0");

        /*
         * Declare Buffer<!LockType> class as a friend class, so we can reference its member functions and variables
         * in move constructor and assignment operators that moves data between instances with different thread safety.
         *
         * alternative is to use inheritance, which would require virtual functions that are potentially expensive.
         */
        template<bliss::concurrent::LockType LT, int64_t C>
        friend std::ostream & operator<<(std::ostream &os, const Buffer<LT, C>& p);

        using BufferType = Buffer<bliss::concurrent::LockType::NONE, Capacity>;


      protected:

        /// internal data storage
        mutable uint8_t* start_ptr; // const, does not change

        mutable bool blocked;

        /// pointer to FINAL  end of data.   only updated when buffer is blocked or when full (from reserved)
        /// capacity + 1 indicates buffer is unblocked.  other values indicate that buffer is blocked or full.
        ///  exchanges with reserved
        mutable int64_t size;

        /// mutex for locking access to the buffer.  available in both thread safe and unsafe versions so we on't need to extensively enable_if or inherit
        mutable std::mutex mutex;

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
        Buffer(BufferType && other, const std::lock_guard<std::mutex> &l)
            : start_ptr(other.start_ptr),
              blocked(other.blocked), size(other.size)
        {

          other.start_ptr = nullptr;
          other.size = 0;
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
        //start_ptr(new uint8_t[Capacity](), []( uint8_t *p ) { delete [] p; })  for shared ptr
        Buffer() : start_ptr(new uint8_t[Capacity]()),
         blocked(true), size(0L)
        {}

        /**
         * @brief Destructor.  waits for all writes and then deallocate memory manually.
         */
        virtual ~Buffer()
        {
          block_and_flush();
        }

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

            std::unique_lock<std::mutex> myLock(mutex, std::defer_lock),
                otherLock(other.mutex, std::defer_lock);
            std::lock(myLock, otherLock);

//            if (start_ptr != nullptr)
//            {
//              delete[] start_ptr;
//              start_ptr = nullptr;
//            }

            /// move the internal memory.

            start_ptr = other.start_ptr;
            size = other.size;
            blocked = other.blocked;

            other.size = 0;
            other.start_ptr = nullptr;
            other.blocked = true;
          }
          return *this;
        }

        int64_t getSize() const
        {
          return size;
        }

      protected:

        /**
         * @brief get the current approximate size of the Buffer.
         * @return    current size
         */
// need to fine tune the memory order so this method may not be useful
        int64_t getApproximateSize() const
        {
          return size;
        }
        /**
         * @brief get the current written data size of the Buffer.
         * @note  written data may be scattered.
         * @return    current size
         */
// need to fine tune the memory order so this method may not be useful
        int64_t getWrittenSize() const
        {
          return size;
        }


        bool is_flushing() const
        {
          return false;
        }



      public:

        /**
         * @brief Get a pointer to the buffer data memory block.
         *
         * @note the data should not be deleted by a calling function/thread.
         * The access is read only.  There is no reason to return the unique_ptr.
         *
         * const because the caller will have a const reference to the buffer
         *
         * return weak pointer for use in multithreaded environment.  This is useful when
         * Buffer is not changing, but the data array is directly accessed by multiple threads.
         */
        uint8_t* getData() const
        {
          return start_ptr;
        }

        /**
         * @brief return a casted pointer to the internal data array.
         * @note  This is better suited when the Buffer itself is shared via shared_ptr and each thread has read access only.
         * @note  this breaks the encapsulation a little.
         */
        template<typename T>
        operator T*() const {
          return reinterpret_cast<T*>(start_ptr);
        }

        /**
         * @brief get the capacity of the buffer.
         * @return    maximum capacity of buffer.
         */
        const int64_t getCapacity() const
        {
          return Capacity;
        }

        /**
         * @brief Checks if a buffer is empty.
         * @note The return value is approximate due to threading
         *
         * @return    true if the buffer is empty, false otherwise.
         */
        // is this used?
        bool isEmpty() const
        {
          return size == 0;
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


        bool is_writing() const
        {
          return false;
        }

        bool is_read_only() const
        {
          return blocked;
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



      protected:

        /**
         * @brief       reserves a position in the buffer at which to insert the count number of bytes.
         * @details     increment only if we have room (so done via atomic CAS.).  if full, return nullptr. if becoming full, set the size.  else return the insertion point
         * @param count
         * @return      pointer at which to insert.
         */


        /**
         * @brief       reserves a position in the buffer at which to insert the count number of bytes.
         * @details     increment only if we have room (so done via atomic CAS.).  if full, return nullptr. if becoming full, set the size.  else return the insertion point
         * @param count
         * @return      pointer at which to insert.
         */
        int64_t reserve(const uint32_t count)
        {
          int64_t curr = size;

          // this part reduces likelihood of accessing freed memory (a different thread tries to reserve to a freed region, i.e. )
          // does not prevent ABA problem (if the memory is reallocated.)
          if (blocked)
          {
            curr = -1;
          }
          else
          { //if ptr == Capacity, buffer must be FILLED in this cycle.  need to process it as if not full yet.

            if (curr + count > Capacity)  // filled to past capacity, so no valid position is returned.
            {
              blocked = true;
              curr = -2;

            } else if (curr + count == Capacity)  // filled to capacity, so valid position is returned.
            {
              blocked = true ;
              size = Capacity;
            }  // else normal append
            else
              size += count;

          }
          return curr;
        }


      public:


//        //===== read lock trumps write.
        void block_writes(int64_t _curr = -1)
        {
          // disable blocked and set the size
          size = (_curr >= 0 ? _curr : size);
          blocked = true;
        }


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
          // no other threads, so no wait for other threads to finish writing.
        }



        /**
         * @brief Clears the buffer. (set the size to 0, leaving the Capacity and the memory allocation intact)
         * @note  const because the caller will have a const reference to the buffer.
         */
        void clear_and_block_writes()
        {
          // blocked
          blocked = true;
          size = 0;
        }

        /**
         * @brief Clears the buffer. (set the size to 0, leaving the Capacity and the memory allocation intact)
         * @note  const because the caller will have a const reference to the buffer.
         */
        void clear_and_unblock_writes()
        {
          // blocked
          blocked = false;
          size = 0;
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
         * This method relies on the current reserved incrementing monotonically (by count) until it exceeds Capacity, and the written
         * incrementing monotonically (by count) as threads finish writes, chasing the reserved thread.
         *
         * The first thread to exceed the Capacity will set the size.  written will never exceed that, and reserved will not
         * be smaller than size.
         *
         *  single primitive type is better as it can be returned in 1 op.
         *
         * @tparam LT Choose thread safe vs unsafe implementation. defaults to same as the parent class.
         * @param[in] _data   pointer to data to be copied into the Buffer..  this SHOULD NOT be shared between threads
         * @param[in] count   number of bytes to be copied into the Buffer
         * @return            unsigned char, bit 0 indicated whether operation succeeded, bit 1 indicating whether buffer swap is needed.
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
          { // write lock returns nullptr if ptr points to outside of the range.
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
            std::memcpy(start_ptr + pos, _data, count);
            //complete_write(count); // all full buffers lock the read and unlock the writer

            // for DEBUGGING.
            _inserted = start_ptr + pos;

            if ((pos + count) == Capacity)
            { // thread that JUST filled the buffer

              // no other threads, so no wait for other threads to finish
              return 0x3;  // write to full, swapping, and success
            }
            else  // no case for pos + count > Capacity.
            { // not filled.

              return 0x1;   // not full, not swapping, successfully inserted.
            }

          }
        }

/// DEBUG WRAPPER
        unsigned int append(const void* _data, const uint32_t count)  // _inserted is for DEBUGGING only.
        {
          void* output = nullptr;
          unsigned int result = append(_data, count, output);

          if (result & 0x1) {
            if (output == nullptr) {
              ERRORF("ERROR: successful insert but result pointer is null.");
            } else if (std::memcmp(_data, output, count) != 0) {
              ERRORF("ERROR: successful insert but input not same as what was memcpy'd.");
            }
          }

          return result;
        }
    };



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
     * life cycle:  pool acquire():  buffer unblock_writes().  -> allow rw.
     *              application:  buffer block_writes() -> allow ro
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
     * data array access via shared_ptr
     *
     * @tparam  LockType  controls whether this class is thread safe or not
     * @tparam  Capacity  size of the allocated byte array
     */
    template<int64_t Capacity>
    class Buffer<bliss::concurrent::LockType::LOCKFREE, Capacity>
    {

        static_assert(Capacity > 0, "Buffer Capacity is given as 0");

        /*
         * Declare Buffer<!LockType> class as a friend class, so we can reference its member functions and variables
         * in move constructor and assignment operators that moves data between instances with different thread safety.
         *
         * alternative is to use inheritance, which would require virtual functions that are potentially expensive.
         */
        template<bliss::concurrent::LockType LT, int64_t C>
        friend std::ostream & operator<<(std::ostream &os, const Buffer<LT, C>& p);

        using BufferType = Buffer<bliss::concurrent::LockType::LOCKFREE, Capacity>;

      protected:

        /// internal data storage
        mutable uint8_t* start_ptr; // const, does not change


        /// pointer to current head of reservation.
        /// capacity + 1 indicates buffer is blocked/full.  other values indicate that buffer is unblocked.
        /// exchanges with size.
        volatile std::atomic<int64_t> reserved;

        /// pointer to FINAL  end of data.   only updated when buffer is blocked or when full (from reserved)
        /// capacity + 1 indicates buffer is unblocked.  other values indicate that buffer is blocked or full.
        ///  exchanges with reserved
        mutable std::atomic<int64_t> size;

        /// represent amount of data written.  will not update beyond the FINAL size
        volatile std::atomic<int64_t> written;

        /// mutex for locking access to the buffer.  available in both thread safe and unsafe versions so we on't need to extensively enable_if or inherit
        mutable std::mutex mutex;

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
        Buffer(BufferType && other, const std::lock_guard<std::mutex> &l)
            : start_ptr(std::move(other.start_ptr)),
              reserved(other.reserved.load()), size(other.size.load()),
              written(other.written.load())
        {

          other.start_ptr = nullptr;
          other.size = 0;
          other.written = 0;
          other.reserved = Capacity + 1;

        };

        /// remove copy constructor and copy assignement operators.
        explicit Buffer(const BufferType& other) = delete;
        /// remove copy constructor and copy assignement operators.
        BufferType& operator=(const BufferType& other) = delete;


      public:
        /**
         * @brief Normal constructor.  Allocate and initialize memory with size specified as parameter.
         */
        //start_ptr(new uint8_t[Capacity](), []( uint8_t *p ) { delete [] p; })  for shared ptr
        Buffer() : start_ptr(new uint8_t[Capacity]()),
         reserved(Capacity + 1), size(0L), written(0L)
        {}

        /**
         * @brief Destructor.  waits for all writes and then deallocate memory manually.
         */
        virtual ~Buffer()
        {
          block_and_flush();
        }

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

            std::unique_lock<std::mutex> myLock(mutex, std::defer_lock),
                otherLock(other.mutex, std::defer_lock);
            std::lock(myLock, otherLock);

//            if (start_ptr != nullptr)
//            {
//              delete[] start_ptr;
//              start_ptr = nullptr;
//            }

            /// move the internal memory.

            //start_ptr.reset(); std::swap(start_ptr, other.start_ptr);
            start_ptr = other.start_ptr;
            written = (int64_t)(other.written);
            size = (int64_t)(other.size);
            reserved = (int64_t)(other.reserved);

            other.size = 0;
            other.start_ptr = nullptr;
            other.written = 0;
            other.reserved = Capacity + 1;
          }
          return *this;
        }

        /**
         * @brief get the current size of the Buffer.
         * @detail  implemented as pointer difference because this is less frequently called compared to append.
         * @return    current size
         */
        int64_t getSize() const
        {
          return size.load(std::memory_order_relaxed);
        }

      protected:

        /**
         * @brief get the current approximate size of the Buffer.
         * @return    current size
         */
// need to fine tune the memory order so this method may not be useful
        int64_t getApproximateSize() const
        {
          return reserved.load(std::memory_order_relaxed);
        }

        /**
         * @brief get the current written data size of the Buffer.
         * @note  written data may be scattered.
         * @return    current size
         */
// need to fine tune the memory order so this method may not be useful
        int64_t getWrittenSize() const
        {
          return written.load(std::memory_order_relaxed);
        }

        bool is_flushing() const
        {
          // is read only - really only for debugging external to this class
          int64_t s = size.load(std::memory_order_relaxed);
          return (reserved.load(std::memory_order_relaxed) > Capacity) && (written.load(std::memory_order_relaxed) < s);
        }



      public:

        /**
         * @brief Get a pointer to the buffer data memory block.
         *
         * @note the data should not be deleted by a calling function/thread.
         * The access is read only.  There is no reason to return the unique_ptr.
         *
         * const because the caller will have a const reference to the buffer
         *
         * return weak pointer for use in multithreaded environment.  This is useful when
         * Buffer is not changing, but the data array is directly accessed by multiple threads.
         */
        uint8_t* getData() const
        {
          return start_ptr;
        }

        /**
         * @brief return a casted pointer to the internal data array.
         * @note  This is better suited when the Buffer itself is shared via shared_ptr and each thread has read access only.
         * @note  this breaks the encapsulation a little.
         */
        template<typename T>
        operator T*() const {
          std::atomic_thread_fence(std::memory_order_release);
          return reinterpret_cast<T*>(start_ptr);
        }

        /**
         * @brief get the capacity of the buffer.
         * @return    maximum capacity of buffer.
         */
        const int64_t getCapacity() const
        {
          return Capacity;
        }

        /**
         * @brief Checks if a buffer is empty.
         * @note The return value is approximate due to threading
         *
         * @return    true if the buffer is empty, false otherwise.
         */
        // is this used?
        bool isEmpty() const
        {
          // if not empty, then need to do some work on the data, e.g. write to stream.  so want side effect visibility
          return written.load(std::memory_order_acquire) == 0;
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

        bool is_writing() const
        {
          // is writing - really only for debugging external to this class
          int64_t s = size.load(std::memory_order_relaxed);
          return written.load(std::memory_order_acquire) < s;
        }


        bool is_read_only() const
        {
          // is read only - really only for debugging external to this class
          int64_t s = size.load(std::memory_order_relaxed);
          return (reserved.load(std::memory_order_relaxed) > Capacity) && (written.load(std::memory_order_acquire) >= s);
        }


      protected:

        /**
         * @brief       reserves a position in the buffer at which to insert the count number of bytes.
         * @details     increment only if we have room (so done via atomic CAS.).  if full, return nullptr. if becoming full, set the size.  else return the insertion point
         * @param count
         * @return      pointer at which to insert.
         */


        /**
         * @brief       reserves a position in the buffer at which to insert the count number of bytes.
         * @details     increment only if we have room (so done via atomic CAS.).  if full, return nullptr. if becoming full, set the size.  else return the insertion point
         * @param count
         * @return      pointer at which to insert.
         */
        int64_t reserve(const uint32_t count)
        {

          // blocked buffer, from full or explicit block call.
          if (reserved.load(std::memory_order_relaxed) > Capacity) {
            // have this here to minimize fetch_add calls.

            return -1;
          } else {
            // not yet full, so reserve

            int64_t curr = reserved.fetch_add(count, std::memory_order_acquire);

            if (curr >= Capacity) {
              return -1;
            }
            else if (curr + count > Capacity) {
              reserved.store(Capacity + 1, std::memory_order_relaxed);
                // just filled. curr position not valid.

              // multiple threads could reach here if they have different counts? NO. because of CAS is atomic.
              // since only 1 thread reaching here, just set size, no need to check to make sure we store minimum.
              size.store(curr, std::memory_order_release);
              return -2;
            } else if (curr + count == Capacity) {
              reserved.store(Capacity + 1, std::memory_order_relaxed);
              // just filled.  curr position valid.

              // multiple threads could reach here if they have different counts? NO. because of CAS is atomic.
              // since only 1 thread reaching here, just set size, no need to check to make sure we store minimum.
              size.store(Capacity, std::memory_order_release);
              return curr;
            } else return curr;
          }
        }


        /// marked write as completed.
        void complete_write(const int count)
        {
          written.fetch_add(count, std::memory_order_release);
        }

      public:

//        //===== read lock trumps write.
        /**
         * @brief  read lock prevents future writes.  read lock can be turned on while there are writes in progress.  once read lock is on, reserve will returns nullptr.
         *
         * @details purpose of this method is to swap the size and reserved, so that reserved is at max+1, and size is the smallest of all locking thread's pointers
         */
        void block_writes(int64_t _curr = -1)
        {
          int64_t curr = reserved.exchange(Capacity + 1, std::memory_order_relaxed);
          // curr could be greater than capacity, or not.

          int64_t end = size.load(std::memory_order_relaxed);  // size initially at Capacity + 1

          curr = (_curr >= 0 ? _curr : curr);
          bool stop = false;

          // multiple threads:  smallest wins.
          while (end > curr && !stop)
          {
            stop = size.compare_exchange_weak(end, curr, std::memory_order_relaxed);
          }

        }


        /**
         * @brief  read lock prevents future writes.  read lock can be turned on while there are writes in progress.  once read lock is on, reserve will returns nullptr.
         *
         * @details purpose of this method is to swap the size and reserved, so that reserved is at max+1, and size is the smallest of all locking thread's pointers
         */
        void unblock_writes()
        {
          // put the size into curr if it hasn't been done.
          int64_t curr = reserved.load(std::memory_order_relaxed);
          int64_t end = size.exchange(Capacity + 1, std::memory_order_relaxed);

          bool stop = false;
          while (curr > Capacity && !stop)
          {
            stop = reserved.compare_exchange_weak(curr, end, std::memory_order_relaxed);
          }

        }


        /**
         * @brief   Sets the buffer as read only and wait for all pending writes.
         */
        void block_and_flush()
        {
          //printf("INFO:  set block start:  reserved %ld, size %ld, written %ld \n", (int64_t)reserved, (int64_t)size, (int64_t)written);
          block_writes();
          //printf("INFO:  set block end:  reserved %ld, size %ld, written %ld \n", (int64_t)reserved, (int64_t)size, (int64_t)written);
          while (is_writing()) {
            //printf("INFO:  wait for writes:  reserved %ld, size %ld, written %ld \n", (int64_t)reserved, (int64_t)size, (int64_t)written);
            _mm_pause();
          }
        }


        void clear_and_block_writes()
        {
          // blocked
          reserved.store(Capacity + 1, std::memory_order_relaxed);
          written.store(0, std::memory_order_relaxed);
          size.store(0, std::memory_order_relaxed);
        }


        /**
         * @brief Clears the buffer. (set the size to 0, leaving the Capacity and the memory allocation intact)
         * @note  const because the caller will have a const reference to the buffer.
         */

        void clear_and_unblock_writes()
        {
          written.store(0, std::memory_order_relaxed);
          reserved.store(0, std::memory_order_relaxed);
          size.store(Capacity + 1, std::memory_order_relaxed);
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
         * This method relies on the current reserved incrementing monotonically (by count) until it exceeds Capacity, and the written
         * incrementing monotonically (by count) as threads finish writes, chasing the reserved thread.
         *
         * The first thread to exceed the Capacity will set the size.  written will never exceed that, and reserved will not
         * be smaller than size.
         *
         *  single primitive type is better as it can be returned in 1 op.
         *
         * @tparam LT Choose thread safe vs unsafe implementation. defaults to same as the parent class.
         * @param[in] _data   pointer to data to be copied into the Buffer..  this SHOULD NOT be shared between threads
         * @param[in] count   number of bytes to be copied into the Buffer
         * @return            unsigned char, bit 0 indicated whether operation succeeded, bit 1 indicating whether buffer swap is needed.
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
          { // write lock returns nullptr if ptr points to outside of the range.
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
            std::memcpy(start_ptr + pos, _data, count);
            complete_write(count); // all full buffers lock the read and unlock the writer

            // for DEBUGGING.
            _inserted = start_ptr + pos;

            if ((pos + count) == Capacity)
            { // thread that JUST filled the buffer

              // wait for other threads to finish
              while (is_writing())
                _mm_pause();

              return 0x3;  // write to full, swapping, and success
            }
            else  // no case for pos + count > Capacity.
            { // not filled.

              return 0x1;   // not full, not swapping, successfully inserted.
            }

          }
        }

/// DEBUG WRAPPER
        unsigned int append(const void* _data, const uint32_t count)  // _inserted is for DEBUGGING only.
        {
          void* output = nullptr;
          unsigned int result = append(_data, count, output);

          if (result & 0x1) {
            if (output == nullptr) {
              ERRORF("ERROR: successful insert but result pointer is null.");
            } else if (std::memcmp(_data, output, count) != 0) {
              ERRORF("ERROR: successful insert but input not same as what was memcpy'd.");
            }
          }

          return result;
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
    template<bliss::concurrent::LockType LT, int64_t C = 8192>
    std::ostream& operator<<(std::ostream& ost, const Buffer<LT, C> & buffer)
    {
      ost << "LockType=" << static_cast<int>(LT)
      << " BUFFER: data_ptr/size="
      << static_cast<const void *>(buffer.start_ptr) << "/"
      << buffer.getSize() << " written="
      << buffer.getWrittenSize() << " approx,size/cap="
      << buffer.getApproximateSize() << "," << buffer.getSize() << "/"
      << buffer.getCapacity() << " R? " << (buffer.is_read_only() ? "y" : "n")
      << " F? " << (buffer.is_flushing() ? "y" : "n") << " W? "
      << (buffer.is_writing() ? "y" : "n") << std::flush;

      return ost;
    }


  } /* namespace io */
} /* namespace bliss */

#endif /* BUFFER_HPP_ */
