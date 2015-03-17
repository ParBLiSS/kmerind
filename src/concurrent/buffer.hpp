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

#include <cstring>   // memset
#include <xmmintrin.h>  // _mm_pause, instead of usleep.


#include <memory>     // unique_ptr
#include <utility>    // move, forward, swap, make_pair
#include <stdexcept>
#include <iostream>   //for std::cout
#include <type_traits>

#include "concurrent/concurrent.hpp"   // LockType boolean constants
#include "config/relacy_config.hpp"
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
    template<bliss::concurrent::LockType LockType, size_t Capacity, size_t MetadataSize>
    class BufferBase
    {
        static_assert(Capacity > 0, "Buffer Capacity is given as 0");

      protected:
        static constexpr size_t FULL_FLAG = std::numeric_limits<size_t>::max();
        static constexpr size_t ALMOST_FULL_FLAG = std::numeric_limits<size_t>::max() - 1;


        static constexpr size_t MAX_SIZE =
            std::numeric_limits<size_t>::max() >> 1;
        static constexpr size_t DISABLED =
            std::numeric_limits<size_t>::max() ^ MAX_SIZE;


        /// internal data storage.  includes metadata
        mutable VAR_T(uint8_t*) start_ptr; // const, does not change
        
        /// pointer to where the actual data starts
        mutable VAR_T(uint8_t*) data_ptr;

        /// mutex for locking access to the buffer.  available in both thread safe and unsafe versions so we on't need to extensively enable_if or inherit
        mutable std::mutex mutex;

        /// pointer to FINAL  end of data.   only updated when buffer is blocked or when full (from reserved)
        /// capacity + 1 indicates buffer is unblocked.  other values indicate that buffer is blocked or full.
        ///  exchanges with reserved
        mutable typename std::conditional<LockType == bliss::concurrent::LockType::LOCKFREE,
            std::atomic<size_t>,
            VAR_T(size_t)>::type size;


        /**
         * @brief Move constructor with a mutex lock on the source object.
         * @details This constructor is private and only use as constructor delegation target.
         *  the source object is locked before this function is called and data moved to the newly constructed object.
         *
         * @param other   the source Buffer
         * @param l       the mutex lock on the source Buffer.
         */
        template<bliss::concurrent::LockType LT = LockType,
            typename std::enable_if<LT == bliss::concurrent::LockType::LOCKFREE, int>::type = 0>
        BufferBase(BufferBase && other, const std::lock_guard<std::mutex> &)
            : start_ptr(other.VAR(start_ptr)), data_ptr(other.VAR(data_ptr))
        {
          std::atomic_thread_fence(std::memory_order_acquire);
          size.store(other.size.exchange(0, std::memory_order_relaxed), std::memory_order_relaxed);
          other.VAR(start_ptr) = nullptr;
          other.VAR(data_ptr) = nullptr;
          std::atomic_thread_fence(std::memory_order_release);
        };

        /**
         * @brief Move constructor with a mutex lock on the source object.
         * @details This constructor is private and only use as constructor delegation target.
         *  the source object is locked before this function is called and data moved to the newly constructed object.
         *
         * @param other   the source Buffer
         * @param l       the mutex lock on the source Buffer.
         */
        template<bliss::concurrent::LockType LT = LockType,
            typename std::enable_if<LT != bliss::concurrent::LockType::LOCKFREE, int>::type = 0>
        BufferBase(BufferBase && other, const std::lock_guard<std::mutex> &)
            : start_ptr(other.VAR(start_ptr)), data_ptr(other.VAR(data_ptr))
        {
          std::atomic_thread_fence(std::memory_order_acquire);
          VAR(size) = other.VAR(size); other.VAR(size) = 0;
          other.VAR(start_ptr) = nullptr;
          other.VAR(data_ptr) = nullptr;
          std::atomic_thread_fence(std::memory_order_release);
        };

        /**
         * @brief Move assignment operator, requiring mutex lock. For between Buffers of the SAME LockType property.
         * @details  Internal data memory moved.
         * The move is done in a thread safe way always.
         *
         * @param other     Source Buffer to move
         */
        template<bliss::concurrent::LockType LT = LockType,
            typename std::enable_if<LT == bliss::concurrent::LockType::LOCKFREE, int>::type = 0>
        void move_assign(BufferBase && other,
            const std::unique_lock<std::mutex> &, const std::unique_lock<std::mutex> &)
        {

          if (this->VAR(start_ptr) != other.VAR(start_ptr))
          {
            // move the internal memory.
            VAR(start_ptr) = other.VAR(start_ptr);
            VAR(data_ptr) = other.VAR(data_ptr);

            std::atomic_thread_fence(std::memory_order_acquire);
            size.store(other.size.exchange(0, std::memory_order_relaxed), std::memory_order_relaxed);

            other.VAR(start_ptr) = nullptr;
            other.VAR(data_ptr) = nullptr;
            std::atomic_thread_fence(std::memory_order_release);
          }
        }

        /**
         * @brief Move assignment operator, requiring mutex lock. For between Buffers of the SAME LockType property.
         * @details  Internal data memory moved.
         * The move is done in a thread safe way always.
         *
         * @param other     Source Buffer to move
         */
        template<bliss::concurrent::LockType LT = LockType,
            typename std::enable_if<LT != bliss::concurrent::LockType::LOCKFREE, int>::type = 0>
        void move_assign(BufferBase && other,
            const std::unique_lock<std::mutex> &, const std::unique_lock<std::mutex> &)
        {

          if (this->VAR(start_ptr) != other.VAR(start_ptr))
          {
            // move the internal memory.
            VAR(start_ptr) = other.VAR(start_ptr);
            VAR(data_ptr) = other.VAR(data_ptr);

            std::atomic_thread_fence(std::memory_order_acquire);
            VAR(size) = other.VAR(size);

            other.VAR(start_ptr) = nullptr;
            other.VAR(data_ptr) = nullptr;
            other.VAR(size) = 0;
            std::atomic_thread_fence(std::memory_order_release);
          }
        }


        /// remove copy constructor and copy assignement operators.
        explicit DELETED_FUNC_DECL(BufferBase(const BufferBase& other));
        /// remove copy constructor and copy assignement operators.
        BufferBase& DELETED_FUNC_DECL(operator=(const BufferBase& other));


        /**
         * @brief Constructor.  Allocate and initialize memory with size specified as parameter.
         */
        BufferBase() : start_ptr(new uint8_t[Capacity + MetadataSize]()), data_ptr(VAR(start_ptr) + MetadataSize), size(0)
        {
          std::atomic_thread_fence(std::memory_order_release);
        }
        //start_ptr(new uint8_t[Capacity](), []( uint8_t *p ) { delete [] p; })  if using shared ptr

        /**
         * @brief Destructor.  waits for all writes and then deallocate memory manually.
         */
        virtual ~BufferBase()
        {
//          DEBUGF("buffer destructor base called");
          delete [] VAR(start_ptr);
          std::atomic_thread_fence(std::memory_order_release);
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
        explicit BufferBase(BufferBase && other)
            : BufferBase(std::move(other),
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
        BufferBase& operator=(BufferBase && other)
        {

          if (this->VAR(start_ptr) != other.VAR(start_ptr))
          {

            std::unique_lock<std::mutex> myLock(mutex, std::defer_lock),
                otherLock(other.mutex, std::defer_lock);
            std::lock(myLock, otherLock);

            // delegate the assignment operator.
            this->move_assign(std::forward<BufferBase >(other), myLock, otherLock);
            
          }
          return *this;
        }

      public:

        /// Get size of buffer.  only meaningful after buffer has been "blocked"
        template<bliss::concurrent::LockType LT = LockType,
            typename std::enable_if<LT == bliss::concurrent::LockType::LOCKFREE, int>::type = 0>
        size_t getSize() const
        {
          return size.load(std::memory_order_relaxed);
        }
        /// Get size of buffer.  only meaningful after buffer has been "blocked"
        template<bliss::concurrent::LockType LT = LockType,
            typename std::enable_if<LT != bliss::concurrent::LockType::LOCKFREE, int>::type = 0>
        size_t getSize() const
        {
          std::atomic_thread_fence(std::memory_order_acquire);
          return VAR(size);
        }

        /**
         * @brief pointer to the start of the metadata block in buffer
         * @note  DO NOT DEALLOCATE THE MEMORY SPACE RETURNED.
         *
         * const because the caller will have a const reference to the buffer
         */
        uint8_t* metadata_begin() const 
        {
          return VAR(start_ptr);
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
          return VAR(data_ptr);
        }
        
        /// pointer to the end of the data block in buffer.
        uint8_t* end() const
        {
          std::atomic_thread_fence(std::memory_order_acquire);
          return VAR(data_ptr) + getSize<LockType>();
        }

        /**
         * @brief cast the internal data buffer to type T and return pointer to user.
         * @note  this breaks the encapsulation a little.
         * @tparam  T   desired data type.
         */
        template<typename T>
        operator T*() const {
          uint8_t* d = VAR(data_ptr);
          return reinterpret_cast<T*>(d);
        }

        /**
         * @brief   get the capacity of the buffer.
         * @return  maximum capacity of buffer.
         */
        const size_t getCapacity() const
        {
          return Capacity;
        }


      	/// get the metadata size
      	const size_t getMetadataSize() const
       	{
      	  return MetadataSize;
       	}


    };


    template<bliss::concurrent::LockType LockType, size_t Capacity = 8192L, size_t MetadataSize = 0UL>
        class Buffer;

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
     *        mutex locked version  relies on mutex to provide memory visibility.
     *
     * @tparam  LockType  controls the type of thread safely.  See concurrent/concurrent.h for allowed values.
     * @tparam  Capacity  size of the allocated byte array
     */
    template<size_t Capacity, size_t MetadataSize>
    class Buffer<bliss::concurrent::LockType::MUTEX, Capacity, MetadataSize> :
    public BufferBase<bliss::concurrent::LockType::MUTEX, Capacity, MetadataSize>
    {
        template<bliss::concurrent::LockType LT, size_t C, size_t MDS>
        friend std::ostream & operator<<(std::ostream &os, const Buffer<LT, C, MDS>& p);

      public:
        static constexpr bliss::concurrent::LockType LockType = bliss::concurrent::LockType::MUTEX;

      protected:

        using BufferType = Buffer<bliss::concurrent::LockType::MUTEX, Capacity, MetadataSize>;
        using BaseType = BufferBase<bliss::concurrent::LockType::MUTEX, Capacity, MetadataSize>;

        // instead of a single "size" variable and a "blocked" variable (to prevent further writes to the buffer)
        //   (which require concurrent updates to both)
        // we use Capacity+1 to indicate a blocked state: if reserved is set to Capacity+1, it is blocked..
        // also, because reservation and actual write may occur separately, we need a "reserved" and a "written"
        //   for separate incrementing.
        // having the "written" variable allows us to wait for all writes to complete.


        /// position of current start of free space.  for reservation.
        /// capacity + 1 indicates buffer is blocked/full.  other values indicate that buffer is unblocked.
        /// exchanges with size.
        VAR_T(size_t) reserved;

        /// amount of data written so far.  will not update beyond the FINAL size
        VAR_T(size_t) written;

        /// condition variable for waiting for waiting for write to catch up
        std::condition_variable condvar;


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
            : BaseType(std::move(other), l)
        {
          VAR(written) = other.VAR(written);   other.VAR(written) = 0;
          VAR(reserved) = other.VAR(reserved); other.VAR(reserved) = Capacity + 1;
        };

        /// remove copy constructor and copy assignement operators.
        explicit DELETED_FUNC_DECL(Buffer(const BufferType& other));
        /// remove copy constructor and copy assignement operators.
        BufferType& DELETED_FUNC_DECL(operator=(const BufferType& other));


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
        explicit Buffer(BufferType && other)
            : BufferType(std::move(other), std::lock_guard<std::mutex>(other.mutex))
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

          if (this->VAR(start_ptr) != other.VAR(start_ptr))
          {

            std::unique_lock<std::mutex> myLock(this->mutex, std::defer_lock),
                otherLock(other.mutex, std::defer_lock);
            std::lock(myLock, otherLock);

            // move the internal memory.
            this->move_assign(std::move(other), myLock, otherLock);

            VAR(written) = (size_t)(other.VAR(written));    other.VAR(written) = 0;
            VAR(reserved) = (size_t)(other.VAR(reserved));  other.VAR(reserved) = Capacity + 1;

          }
          return *this;
        }


        /**
         * @brief   Checks if a buffer is empty.
         * @note    The return value is approximate due to threading
         *
         * @return  true if the buffer is empty, false otherwise.
         */
        bool isEmpty() const
        {
          std::lock_guard<std::mutex> lock(this->mutex);
          return (VAR(written) == 0);
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
        bool is_writing() const
        {
          std::lock_guard<std::mutex> lock(this->mutex);
          return VAR(written) < this->VAR(size);  // reserved > written if is writing. but reserve can increase while buffer is unblocked.
        }

        /// check if the buffer is blocked from further writes and no threads need to write to it.  uses mutex
        bool is_read_only() const
        {
          std::lock_guard<std::mutex> lock(this->mutex);
          return (VAR(reserved) > Capacity) && (VAR(written) >= this->VAR(size));
        }


        /// marked write as completed. internal use only
        void complete_write(const int count)
        {
          std::unique_lock<std::mutex> lock(this->mutex);
          VAR(written) += count;
          lock.unlock();

          CV_NOTIFY_ALL(condvar);
        }

      protected:




        /**
         * @brief       reserves a position in the buffer at which to insert the count number of bytes.
         * @details     increment only if we have room (protected by mutex or spinlock).
         *                if full, return FULL_FLAG.
         *                if becoming over-capacity, mark as full and return ALMOST_FULL_FLAG to indicate it just became full.
         *                if becoming exactly full, mark as full and return the reservation
         *                else return the reservation.
         * @param       count of bytes to reserve
         * @return      the position in buffer where the reservation starts.
         */

        size_t reserve(const uint32_t count)
        {
          std::lock_guard<std::mutex> lock(this->mutex);
          size_t curr = VAR(reserved);

          // this part reduces likelihood of accessing freed memory (a different thread tries to reserve to a freed region, i.e. )
          // does not prevent ABA problem (if the memory is reallocated.)
          if (curr >= Capacity)
          {
            curr = BaseType::FULL_FLAG;
          }
          else
          {
            if (curr + count > Capacity)  // filled to past capacity, so no valid position is returned.
            {
              VAR(reserved) = Capacity + 1 ;
              this->VAR(size) = curr;

              curr = BaseType::ALMOST_FULL_FLAG;

            }
            else if (curr + count == Capacity)  // filled to capacity, so valid position is returned.
            {
              VAR(reserved) = Capacity + 1 ;
              this->VAR(size) = Capacity;

            }  // else normal append
            else {
              VAR(reserved) += count;

            }
          }
          return curr;
        }

      public:

        /**
         * @brief  prevents future writes.  block can be turned on while there are writes in progress.  once block is on, reserve will fail.
         *
         * @details purpose of this method is to swap the size and reserved, so that reserved is at max+1, and size is the smallest of all locking thread's pointers
         *         user can specify a size at which which the block occurs.
         */
        void block_writes()
        {
          std::lock_guard<std::mutex> lock(this->mutex);

          // if supplied a desired size at which to lock, use it.
          size_t curr = VAR(reserved);

          // set reserve so no further reservations occur.
          VAR(reserved) = Capacity + 1;

          // choose smallest curr value amount all threads.
          if (this->VAR(size) > curr)
          {
        	  this->VAR(size) = curr;
          }

        }


        /**
         * @brief  allow future writes.  block can be turned off while there are writes in progress.  once block is off, reserve may proceed
         *
         * @details purpose of this method is to swap the size and reserved, so that reserved is again below max, and size is the at max+1
         */
        void unblock_writes()
        {
          std::lock_guard<std::mutex> lock(this->mutex);

          if (VAR(reserved) > Capacity)
            VAR(reserved) = this->VAR(size);

          // then reset size to max
          this->VAR(size) = Capacity + 1;
        }

        /**
         * @brief   Sets the buffer as blocked and wait for all pending writes.
         */
        void block_and_flush()
        {
          std::unique_lock<std::mutex> lock(this->mutex);

          // if supplied a desired size at which to lock, use it.
          if (this->VAR(size) > VAR(reserved))
          {
            this->VAR(size) = VAR(reserved);
          }
          // set reserve so no further reservations occur.
          VAR(reserved) = Capacity + 1;


          while (VAR(written) < this->VAR(size)) {
            CV_WAIT(condvar, lock);
          }

          lock.unlock();
        }


        void flush() {
          std::unique_lock<std::mutex> lock(this->mutex);

          while (VAR(written) < this->VAR(size)) {
            CV_WAIT(condvar, lock);
          }

          lock.unlock();
        }


        /**
         * @brief Clears the buffer. (set the size to 0, leaving the Capacity and the memory allocation intact).  also blocks future writes
         * @note  const because the caller will have a const reference to the buffer.
         */
        void clear_and_block_writes()
        {
          std::unique_lock<std::mutex> lock(this->mutex);
          VAR(reserved) = Capacity + 1;
          this->VAR(size) = 0;
          VAR(written) = 0;

          lock.unlock();
          CV_NOTIFY_ALL(condvar);

        }


        /**
         * @brief Clears the buffer. (set the size to 0, leaving the Capacity and the memory allocation intact).  also unblocks the buffer for writes
         * @note  const because the caller will have a const reference to the buffer.
         */
        void clear_and_unblock_writes()
        {
          std::unique_lock<std::mutex> lock(this->mutex);
          VAR(reserved) = 0;
          this->VAR(size) = Capacity + 1;
          VAR(written) = 0;

          lock.unlock();
          CV_NOTIFY_ALL(condvar);
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
          size_t pos = reserve(count);
          _inserted = nullptr;

          // if fails, then already full
          if (pos == BaseType::FULL_FLAG)
          {
            // already read locked at this point

            return 0x0;  // full and not swapping.
          }
          else if (pos == BaseType::ALMOST_FULL_FLAG)
          { // filled to beyond capacity

            flush();
            return 0x2;  // full, swapping.
          }
          else
          { // valid position returned, so can write.

            // write
            std::memcpy(this->VAR(data_ptr) + pos, _data, count);
            complete_write(count); // all full buffers lock the read and unlock the writer

            // for DEBUGGING.
            _inserted = this->VAR(data_ptr) + pos;

            if ((pos + count) == Capacity)
            { // thread that JUST filled the buffer

              // wait for other threads to finish
              flush();

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
    template<size_t Capacity, size_t MetadataSize>
    class Buffer<bliss::concurrent::LockType::SPINLOCK, Capacity, MetadataSize> :
      public BufferBase<bliss::concurrent::LockType::SPINLOCK, Capacity, MetadataSize>
    {
        template<bliss::concurrent::LockType LT, size_t C, size_t MDS>
        friend std::ostream & operator<<(std::ostream &os, const Buffer<LT, C, MDS>& p);


      public:
        static constexpr bliss::concurrent::LockType LockType = bliss::concurrent::LockType::SPINLOCK;

      protected:
        using BufferType = Buffer<bliss::concurrent::LockType::SPINLOCK, Capacity, MetadataSize>;
        using BaseType = BufferBase<bliss::concurrent::LockType::SPINLOCK, Capacity, MetadataSize>;

        // instead of a single "size" variable and a "blocked" variable (to prevent further writes to the buffer)
        //   (which require concurrent updates to both)
        // we use Capacity+1 to indicate a blocked state: if reserved is set to Capacity+1, it is blocked..
        // also, because reservation and actual write may occur separately, we need a "reserved" and a "written"
        //   for separate incrementing.
        // having the "written" variable allows us to wait for all writes to complete.


        /// position of current start of free space.  for reservation.
        /// capacity + 1 indicates buffer is blocked/full.  other values indicate that buffer is unblocked.
        /// exchanges with size.
        VAR_T(size_t) reserved;

        /// amount of data written so far.  will not update beyond the FINAL size
        VAR_T(size_t) written;


        /// spinlock for locking access to the buffer.  available in both thread safe and unsafe versions so we on't need to extensively enable_if or inherit
        mutable INIT_ATOMIC_FLAG(spinlock);

      private:

        /// remove copy constructor and copy assignement operators.
        explicit DELETED_FUNC_DECL(Buffer(const BufferType& other));
        /// remove copy constructor and copy assignement operators.
        BufferType& DELETED_FUNC_DECL(operator=(const BufferType& other));


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
        explicit Buffer(BufferType && other) {
          while (other.spinlock.test_and_set(std::memory_order_acquire));

          this->VAR(size) = other.VAR(size); other.VAR(size) = 0;
          this->VAR(start_ptr) = other.VAR(start_ptr); other.VAR(start_ptr) = nullptr;
          this->VAR(data_ptr) = other.VAR(data_ptr); other.VAR(data_ptr) = nullptr;
          VAR(written) = other.VAR(written); other.VAR(written) = 0;
          VAR(reserved) = other.VAR(reserved); other.VAR(reserved) = Capacity + 1;

          other.spinlock.clear(std::memory_order_release);
        }



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

          if (this->VAR(start_ptr) != other.VAR(start_ptr))
          {

            while(spinlock.test_and_set(std::memory_order_acquire));
            while(other.spinlock.test_and_set(std::memory_order_acquire));

            // move the internal memory.
            this->VAR(start_ptr) = other.VAR(start_ptr);                 other.VAR(start_ptr) = nullptr;
            this->VAR(data_ptr) = other.VAR(data_ptr);                   other.VAR(data_ptr) = nullptr;
            this->VAR(size) = other.VAR(size);                           other.VAR(size) = 0;

            VAR(written) = (size_t)(other.VAR(written));                other.VAR(written) = 0;
            VAR(reserved) = (size_t)(other.VAR(reserved));              other.VAR(reserved) = Capacity + 1;

            other.spinlock.clear(std::memory_order_release);
            spinlock.clear(std::memory_order_release);
          }
          return *this;

        }
      protected:

        /**
         * @brief get the current approximate size of the Buffer.
         * @return    current size
         */
/*        size_t getApproximateSize() const
        {
          std::atomic_thread_fence(std::memory_order_acquire);
          return reserved;
        }
*/
        /**
         * @brief get the current written data size of the Buffer.
         * @note  written data may be scattered.
         * @return    current size
         */
/*        size_t getWrittenSize() const
        {
          std::atomic_thread_fence(std::memory_order_acquire);
          return written;
        }
*/
        /// Check if buffer has been blocked from further reservation but still waiting for all writes to finish. via spin lock
/*        bool is_flushing() const
        {
          while (spinlock.test_and_set(std::memory_order_acquire));
          bool res = (reserved > Capacity) && (written < this->size);
          spinlock.clear(std::memory_order_release);
          return res;
        }
*/

      public:


        /**
         * @brief   Checks if a buffer is empty.
         * @note    The return value is approximate due to threading
         *
         * @return  true if the buffer is empty, false otherwise.
         */
        bool isEmpty() const
        {
          while (spinlock.test_and_set(std::memory_order_acquire));
          bool res = (VAR(written) == 0);
          spinlock.clear(std::memory_order_release);
          return res;
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
        bool is_writing() const
        {
          while (spinlock.test_and_set(std::memory_order_acquire));
          bool res = VAR(written) < this->VAR(size);
          spinlock.clear(std::memory_order_release);
          return res;
        }

        /// check if the buffer is blocked from further writes and no threads need to write to it.  uses spinlock
        bool is_read_only() const
        {
          while (spinlock.test_and_set(std::memory_order_acquire));
          bool res = (VAR(reserved) > Capacity) && (VAR(written) >= this->VAR(size));
          spinlock.clear(std::memory_order_release);
          return res;
        }




        /// marked write as completed. spinlock version
        void complete_write(const int count)
        {
          while (spinlock.test_and_set(std::memory_order_acquire));
          VAR(written) += count;
          spinlock.clear(std::memory_order_release);
        }


        /**
         * @brief       reserves a position in the buffer at which to insert the count number of bytes.
         * @details     increment only if we have room (protected by mutex or spinlock).
         *                if full, return FULL_FLAG.
         *                if becoming over-capacity, mark as full and return ALMOST_FULL_FLAG to indicate it just became full.
         *                if becoming exactly full, mark as full and return the reservation
         *                else return the reservation.
         * @param       count of bytes to reserve
         * @return      the position in buffer where the reservation starts.
         */
        size_t reserve(const uint32_t count)
        {
          while (spinlock.test_and_set(std::memory_order_acquire));
          size_t curr = VAR(reserved);

          // this part reduces likelihood of accessing freed memory (a different thread tries to reserve to a freed region, i.e. )
          // does not prevent ABA problem (if the memory is reallocated.)
          if (curr >= Capacity)
          {
            curr = BaseType::FULL_FLAG;
          }
          else
          {
            if (curr + count > Capacity)  // filled to past capacity, so no valid position is returned.
            {
              VAR(reserved) = Capacity + 1 ;
              this->VAR(size) = curr;

              curr = BaseType::ALMOST_FULL_FLAG;

            }
            else if (curr + count == Capacity)  // filled to capacity, so valid position is returned.
            {
              VAR(reserved) = Capacity + 1 ;
              this->VAR(size) = Capacity;

            }  // else normal append
            else {
              VAR(reserved) += count;

            }
          }

          spinlock.clear(std::memory_order_release);

          return curr;
        }



        /**
         * @brief  prevents future writes.  block can be turned on while there are writes in progress.  once block is on, reserve will fail.
         *
         * @details purpose of this method is to swap the size and reserved, so that reserved is at max+1, and size is the smallest of all locking thread's pointers
         *         user can specify a size at which which the block occurs.
         */
        void block_writes()
        {
          while (spinlock.test_and_set(std::memory_order_acquire));

          // if supplied a desired size at which to lock, use it.

          size_t curr = VAR(reserved);

          // set reserve so no further reservations occur.
          VAR(reserved) = Capacity + 1;

          // choose smallest curr value amount all threads.
          if (this->VAR(size) > curr)
          {
            this->VAR(size) = curr;
          }

          spinlock.clear(std::memory_order_release);
        }



        /**
         * @brief  allow future writes.  block can be turned off while there are writes in progress.  once block is off, reserve may proceed
         *
         * @details purpose of this method is to swap the size and reserved, so that reserved is again below max, and size is the at max+1
         */
        void unblock_writes()
        {
          while (spinlock.test_and_set(std::memory_order_acquire));

          if (VAR(reserved) > Capacity)
            VAR(reserved) = this->VAR(size);

          // then reset size to max
          this->VAR(size) = Capacity + 1;

          spinlock.clear(std::memory_order_release);

        }



        /**
         * @brief   Sets the buffer as blocked and wait for all pending writes.
         */
        void block_and_flush()
        {
          while (spinlock.test_and_set(std::memory_order_acquire));

          // choose smallest curr value amount all threads.
          if (this->VAR(size) > VAR(reserved))
          {
            this->VAR(size) = VAR(reserved);
          }

          // set reserve so no further reservations occur.
          VAR(reserved) = Capacity + 1;

          while (VAR(written) < this->VAR(size)) {
            spinlock.clear(std::memory_order_release);
            _mm_pause();
            while (spinlock.test_and_set(std::memory_order_acquire));
          }
          
          spinlock.clear(std::memory_order_release);
        }


        void flush() {
          while (spinlock.test_and_set(std::memory_order_acquire));

          while (VAR(written) < this->VAR(size)) {
            spinlock.clear(std::memory_order_release);
            _mm_pause();
            while (spinlock.test_and_set(std::memory_order_acquire));
          }
          
          spinlock.clear(std::memory_order_release);
        }


        /**
         * @brief Clears the buffer. (set the size to 0, leaving the Capacity and the memory allocation intact).  also blocks future writes
         * @note  const because the caller will have a const reference to the buffer.
         */
        void clear_and_block_writes()
        {
          while (spinlock.test_and_set(std::memory_order_acquire)) ;
          VAR(reserved) = Capacity + 1;
          this->VAR(size) = 0;
          VAR(written) = 0;
          spinlock.clear(std::memory_order_release);
        }

        /**
         * @brief Clears the buffer. (set the size to 0, leaving the Capacity and the memory allocation intact).  also unblocks the buffer for writes
         * @note  const because the caller will have a const reference to the buffer.
         */
        void clear_and_unblock_writes()
        {
          while (spinlock.test_and_set(std::memory_order_acquire));
          VAR(reserved) = 0;
          this->VAR(size) = Capacity + 1;
          VAR(written) = 0;

          spinlock.clear(std::memory_order_release);

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
          size_t pos = reserve(count);
          _inserted = nullptr;

          // if fails, then already full
          if (pos == BaseType::FULL_FLAG)
          {
            // already read locked at this point

            return 0x0;  // full and not swapping.
          }
          else if (pos == BaseType::ALMOST_FULL_FLAG)
          { // filled to beyond capacity

            flush();  // wait for other writes to complete
            return 0x2;  // full, swapping.
          }
          else
          { // valid position returned, so can write.

            // write
            std::memcpy(this->VAR(data_ptr) + pos, _data, count);
            complete_write(count); // all full buffers lock the read and unlock the writer

            // for DEBUGGING.
            _inserted = this->VAR(data_ptr) + pos;

            if ((pos + count) == Capacity)
            { // thread that JUST filled the buffer

              // wait for other threads to finish
              flush();
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
    template<size_t Capacity, size_t MetadataSize>
    class Buffer<bliss::concurrent::LockType::NONE, Capacity, MetadataSize> :
      public BufferBase<bliss::concurrent::LockType::NONE, Capacity, MetadataSize>
    {

        template<bliss::concurrent::LockType LT, size_t C, size_t MDS>
        friend std::ostream & operator<<(std::ostream &os, const Buffer<LT, C, MDS>& p);


      public:
        static constexpr bliss::concurrent::LockType LockType = bliss::concurrent::LockType::NONE;


      protected:
        using BufferType = Buffer<bliss::concurrent::LockType::NONE, Capacity, MetadataSize>;
        using BaseType = BufferBase<bliss::concurrent::LockType::NONE, Capacity, MetadataSize>;
      
        // since not threaded, uses a single size member variable  and blocked flag. - reserved == written == size

        /// inidicate if the buffer is accepting further append
        mutable VAR_T(bool) blocked;


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
            : BaseType(std::move(other), l),
              blocked(other.VAR(blocked))
        {
          other.VAR(blocked) = Capacity + 1;
        };

        /// remove copy constructor and copy assignement operators.
        explicit DELETED_FUNC_DECL(Buffer(const BufferType& other));
        /// remove copy constructor and copy assignement operators.
        BufferType& DELETED_FUNC_DECL(operator=(const BufferType& other));


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

          if (this->VAR(start_ptr) != other.VAR(start_ptr))
          {

            std::unique_lock<std::mutex> myLock(this->mutex, std::defer_lock),
                otherLock(other.mutex, std::defer_lock);
            std::lock(myLock, otherLock);

            this->move_assign(std::move(other), myLock, otherLock);

            /// move the internal memory.
            VAR(blocked) = other.VAR(blocked);  other.VAR(blocked) = true;
          }
          return *this;
        }

//      protected:

        /**
         * @brief get the current approximate size of the Buffer.
         * @return    current size
         */
/*        size_t getApproximateSize() const
        {
          return this->size;
        }
*/        /**
         * @brief get the current written data size of the Buffer.
         * @note  written data may be scattered.
         * @return    current size
         */
/*        size_t getWrittenSize() const
        {
          return this->size;
        }
*/
        /// Check if buffer has been blocked from further reservation but still waiting for all writes to finish.
/*        bool is_flushing() const
        {
          return false;
        }
*/


      public:

        /**
         * @brief Checks if a buffer is empty.
         *
         * @return    true if the buffer is empty, false otherwise.
         */
        bool isEmpty() const
        {
          return this->VAR(size) == 0;
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
          return VAR(blocked);
        }

        void complete_write(const uint32_t count) {}


        /**
          * @brief       reserves a position in the buffer at which to insert count number of bytes.
          * @details     increment only if we have room
         *                if full, return FULL_FLAG.
         *                if becoming over-capacity, mark as full and return ALMOST_FULL_FLAG to indicate it just became full.
         *                if becoming exactly full, mark as full and return the reservation
         *                else return the reservation.
          * @param count
          * @return      pointer at which to insert.
          */
        size_t reserve(const uint32_t count)
        {
          size_t curr = this->VAR(size);

          if (VAR(blocked))
          {
            curr = BaseType::FULL_FLAG;
          }
          else
          {
            if (curr + count > Capacity)  // filled to past capacity, so no valid position is returned.
            {
              VAR(blocked) = true;
              curr = BaseType::ALMOST_FULL_FLAG;

            } else if (curr + count == Capacity)  // filled to capacity, so valid position is returned.
            {
              VAR(blocked) = true ;
              this->VAR(size) = Capacity;
            }  // else normal append
            else
              this->VAR(size) += count;

          }
          return curr;
        }


        /**
         * @brief  prevents future writes.  block can be turned on while there are writes in progress.  once block is on, reserve will fail.
         *
         * @details purpose of this method is to swap the size and reserved, so that reserved is at max+1, and size is the smallest of all locking thread's pointers
         *         user can specify a size at which which the block occurs.
         */
        void block_writes()
        {
          VAR(blocked) = true;
        }

        /**
         * @brief  allow future writes.  block can be turned off while there are writes in progress.  once block is off, reserve may proceed
         *
         * @details purpose of this method is to swap the size and reserved, so that reserved is again below max, and size is the at max+1
         */
        void unblock_writes()
        {
          VAR(blocked) = false;
        }

        /**
         * @brief   Sets the buffer as read only and wait for all pending writes.
         */
        void block_and_flush()
        {
          VAR(blocked) = true;
        }



        /**
         * @brief Clears the buffer. (set the size to 0, leaving the Capacity and the memory allocation intact).  also blocks future writes
         * @note  const because the caller will have a const reference to the buffer.
         */
        void clear_and_block_writes()
        {
          // blocked
          VAR(blocked) = true;
          this->VAR(size) = 0;
        }

        /**
         * @brief Clears the buffer. (set the size to 0, leaving the Capacity and the memory allocation intact).  also unblocks the buffer for writes
         * @note  const because the caller will have a const reference to the buffer.
         */
        void clear_and_unblock_writes()
        {
          // blocked
          VAR(blocked) = false;
          this->VAR(size) = 0;
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
          size_t pos = reserve(count);
          _inserted = nullptr;

          // if fails, then already full
          if (pos == BaseType::FULL_FLAG)
          {
            // already read locked.
            return 0x0;  // full and not swapping.
          }
          else if (pos == BaseType::ALMOST_FULL_FLAG)
          { // filled to beyond capacity

            // no other threads, no wait for other thread to finish writing.
            return 0x2;  // full, swapping.
          }
          else
          { // valid position returned, so can write.

            // write
            std::memcpy(this->VAR(data_ptr) + pos, _data, count);

            // for DEBUGGING.
            _inserted = this->VAR(data_ptr) + pos;

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
     * as we increment during reserve and reset to Cap + 1, we could have ABA type of threading behavior, but it is not an issue.
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
    template<size_t Capacity, size_t MetadataSize>
    class Buffer<bliss::concurrent::LockType::LOCKFREE, Capacity, MetadataSize> :
      public BufferBase<bliss::concurrent::LockType::LOCKFREE, Capacity, MetadataSize>
    {
        template<bliss::concurrent::LockType LT, size_t C, size_t MDS>
        friend std::ostream & operator<<(std::ostream &os, const Buffer<LT, C, MDS>& p);

      public:
        static constexpr bliss::concurrent::LockType LockType = bliss::concurrent::LockType::LOCKFREE;


      protected:
        using BufferType = Buffer<bliss::concurrent::LockType::LOCKFREE, Capacity, MetadataSize>;        
        using BaseType = BufferBase<bliss::concurrent::LockType::LOCKFREE, Capacity, MetadataSize>;

        // instead of a single "size" variable and a "blocked" variable (to prevent further writes to the buffer)
        //   (which require concurrent updates to both)
        // we use the sign bit to indicate "blocked".
        // also, because reservation and actual write may occur separately, we need a "reserved" and a "written"
        //   for separate incrementing.
        // having the "written" variable allows us to wait for all writes to complete.

        /// pointer to current head of reservation.
        /// capacity + 1 indicates buffer is blocked/full.  other values indicate that buffer is unblocked.
        /// exchanges with size.
        std::atomic<size_t> reserved;


        /// represent amount of data written.  will not update beyond the FINAL size
        std::atomic<size_t> written;

      private:

        /// remove copy constructor and copy assignement operators.
        explicit DELETED_FUNC_DECL(Buffer(const BufferType& other));
        /// remove copy constructor and copy assignement operators.
        BufferType& DELETED_FUNC_DECL(operator=(const BufferType& other));


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
        {
          std::atomic_thread_fence(std::memory_order_acquire);

          this->VAR(start_ptr) = other.VAR(start_ptr);
          this->VAR(data_ptr) = other.VAR(data_ptr);

    			reserved.store(other.reserved.exchange(Capacity+1, std::memory_order_relaxed), std::memory_order_relaxed);
		    	written.store(other.written.exchange(0, std::memory_order_relaxed), std::memory_order_relaxed);
		    	this->size.store(other.size.exchange(0, std::memory_order_relaxed), std::memory_order_relaxed);

          other.VAR(start_ptr) = nullptr;
          other.VAR(data_ptr) = nullptr;

          std::atomic_thread_fence(std::memory_order_release);

        }

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

          if (this->VAR(start_ptr) != other.VAR(start_ptr))
          {

            std::atomic_thread_fence(std::memory_order_acquire);

            this->VAR(start_ptr) = other.VAR(start_ptr);
            this->VAR(data_ptr) = other.VAR(data_ptr);

    	  		reserved.store(other.reserved.exchange(Capacity+1, std::memory_order_relaxed), std::memory_order_relaxed);
		      	written.store(other.written.exchange(0, std::memory_order_relaxed), std::memory_order_relaxed);
		      	this->size.store(other.size.exchange(0, std::memory_order_relaxed), std::memory_order_relaxed);

            other.VAR(start_ptr) = nullptr;
            other.VAR(data_ptr) = nullptr;

            std::atomic_thread_fence(std::memory_order_release);
          }
          return *this;
        }

//      protected:
//
//        /**
//         * @brief get the current approximate size of the Buffer.
//         * @return    current size
//         */
//        size_t getApproximateSize() const
//        {
//          return reserved.load(std::memory_order_relaxed);
//        }
//
//        /**
//         * @brief get the current written data size of the Buffer.
//         * @note  written data may be scattered.
//         * @return    current size
//         */
//        size_t getWrittenSize() const
//        {
//          return written.load(std::memory_order_relaxed);
//        }
//
//        /// Check if buffer has been blocked from further reservation but still waiting for all writes to finish.
//        bool is_flushing() const
//        {
//
//        	std::atomic_thread_fence(std::memory_order_acquire);
//        	auto s = this->VAR(size);
//          bool res = (reserved.load(std::memory_order_relaxed) > Capacity) && (s > written.load(std::memory_order_relaxed));
//          std::atomic_thread_fence(std::memory_order_release);
//          return res;
//        }
//
//
//
//      public:

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
          size_t v = written.load(std::memory_order_acquire);
          return this->getSize() > v;
        }

        /// check if the buffer is blocked from further writes and no threads need to write to it.
        bool is_read_only() const
        {
          // is read only - really only for debugging external to this class
        	std::atomic_thread_fence(std::memory_order_acquire);
        	bool res =  (reserved.load(std::memory_order_relaxed) > Capacity) &&
        	    (this->getSize() <= written.load(std::memory_order_relaxed));
          std::atomic_thread_fence(std::memory_order_release);
          return res;
        }


        /// marked write as completed.
        void complete_write(const int count)
        {
          written.fetch_add(count, std::memory_order_acq_rel);
        }

      protected:


        /**
         * @brief       reserves a position in the buffer at which to insert the count number of bytes.
         * @details     increment only if we have room (protected by mutex or spinlock).
         *                if full, return FULL_FLAG.
         *                if becoming over-capacity, mark as full and return ALMOST_FULL_FLAG to indicate it just became full.
         *                if becoming exactly full, mark as full and return the reservation
         *                else return the reservation.
         * @param       count of bytes to reserve
         * @return      the position in buffer where the reservation starts.
         */
        size_t reserve(const uint32_t count)
        {
            assert(count <= Capacity);

            //if (reserved.load(std::memory_order_acquire) > Capacity) return BaseType::FULL_FLAG;   // slows down the reserve.

            

            // fetch_add here - if too many calls we may wrap around.  rely on resetting to Capacity + 1.
            size_t curr = reserved.fetch_add(count, std::memory_order_acquire);
  
            size_t next = curr + count;
            //--- 4 cases:
            //  curr >= Capacity          full/disabled.
            //       > Capacity - count   just filled and no room for last reserve
            //       == Capacity - count  just filled and exactly room for last reserve
            //       < Capacity - count   normal

            if (curr >= Capacity) {
              // full or blocked.

              // reset reserved only when we JUST FILL the buffer, or when there have been too many attempts to reserve a full buffer.
              //  otherwise let the reserved value grow.  this reduces ABA problem very significantly.
              if (curr > BaseType::DISABLED) {
                reserved.compare_exchange_strong(next, Capacity + 1, std::memory_order_acq_rel);
              }

              curr = BaseType::FULL_FLAG;
            } else if (next > Capacity) {
              // just filled. curr position not a valid insertion point.

              // since only 1 thread reaching here (resv is at next > Capacity so no subsequent threads will get here),
              // just set size, no need to check to make sure we store minimum or maximum
              this->size.store(curr, std::memory_order_relaxed);  // new size is the prev reserved.  reserved will be set to Capacity + 1
              std::atomic_thread_fence(std::memory_order_release);
      
              curr = BaseType::ALMOST_FULL_FLAG;
            } else if (next == Capacity) {
              // just filled.  curr position is a valid insertion point.

              // reset reserved only when we JUST FILL the buffer, or when there have been too many attempts to reserve a full buffer.
              //  otherwise let the reserved value grow.  this reduces ABA problem very significantly.
              reserved.compare_exchange_strong(next, Capacity + 1, std::memory_order_acq_rel);

              // since only 1 thread reaching here, just set size, no need to check to make sure we store minimum.              
              this->size.store(Capacity, std::memory_order_relaxed);
              // curr is okay.
              std::atomic_thread_fence(std::memory_order_release);

            } 
            // else normal case. nothing to do.

            std::atomic_thread_fence(std::memory_order_release);
            return curr;
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
        void block_writes()
        {
          // this line below sequences changes to reserved between diff threads calling blocked writes or reserve or unblock writes.
          // size update depends on value of reserved so we would not have 2 updates to size by different threads.
          size_t curr = reserved.exchange(Capacity + 1, std::memory_order_acq_rel);
          // curr could be greater than capacity, or not.
          size_t s = this->getSize();

          if (curr > Capacity) {  // already disabled.
//            WARNINGF("block_writes called, but already blocked. reserved: %ld, size: %ld", curr, s);
            return;
          }

          // thre should only be 1 thread reaching here.  reaching here means  reserve does not change size.
          if (this->size.exchange(curr, std::memory_order_relaxed) <= Capacity) {
            ERRORF("block_writes called, but size already set. reserved: %ld, size: %ld", curr, s);
          }

//          // multiple threads:  smallest wins.
//          bool done = false;
//          while (s > curr && !done) {
//            done = this->size.compare_exchange_weak(s, curr, std::memory_order_acq_rel);
//          }
//          std::atomic_thread_fence(std::memory_order_release);

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
//          size_t curr = reserved.load(std::memory_order_relaxed);
//          size_t end = size.exchange(Capacity + 1, std::memory_order_relaxed);
//
//          bool stop = false;
//          while (curr > Capacity && !stop)
//          {
//            stop = reserved.compare_exchange_weak(curr, end, std::memory_order_relaxed);
//          }
          size_t s = this->size.exchange(Capacity + 1, std::memory_order_acq_rel);
          size_t curr = reserved.load(std::memory_order_relaxed);

          if (s > Capacity) {
            //WARNINGF("unblock_writes called, but another thread already beat this one to it.  reserved: %ld, size: %ld", curr, s);
            return;
          }

          if (this->reserved.exchange(s, std::memory_order_relaxed) <= Capacity) {
            ERRORF("unblock_writes called, but reserved already unblocked. reserved: %ld, size: %ld", curr, s);
          }

//          bool done = false;
//        	while (curr > s && !done) {
//            done = reserved.compare_exchange_weak(curr, s, std::memory_order_acq_rel);
//        	}
          std::atomic_thread_fence(std::memory_order_release);
        }


        /**
         * @brief   Sets the buffer as read only and wait for all pending writes.
         */
        void block_and_flush()
        {
          this->block_writes();
          this->flush();
        }

        void flush() {
          std::atomic_thread_fence(std::memory_order_acquire);
          while (written.load(std::memory_order_acquire) < this->getSize()) {
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
          this->size.store(0, std::memory_order_relaxed);
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
          this->size.store(Capacity + 1, std::memory_order_relaxed);
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
          size_t pos = reserve(count);
          _inserted = nullptr;

          // if fails, then already full
          if (pos == BaseType::FULL_FLAG)
          {
            // already read locked.

            return 0x0;  // full and not swapping.
          }
          else if (pos == BaseType::ALMOST_FULL_FLAG)
          { // filled to beyond capacity

           flush();  // wait for other writes to complete
            return 0x2;  // full, swapping.
          }
          else
          { // valid position returned, so can write.
            // for DEBUGGING.
            std::atomic_thread_fence(std::memory_order_acquire);
            _inserted = this->VAR(data_ptr) + pos;

            // write
            std::memcpy(_inserted, _data, count);
            std::atomic_thread_fence(std::memory_order_release);
            complete_write(count); // all full buffers lock the read and unlock the writer


            if ((pos + count) == Capacity)
            { // thread that JUST filled the buffer

              // wait for other threads to finish
              flush();
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
          //uint8_t* ptr = this->VAR(data_ptr);
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
    template<bliss::concurrent::LockType LT, size_t C, size_t MDS>
    std::ostream & operator<<(std::ostream &ost, const Buffer<LT, C, MDS>& buffer)
    {
      ost << "LockType=" << static_cast<int>(LT)
      << " BUFFER: this->VAR(data_ptr)/size="
      << static_cast<const void *>(buffer.start_ptr) << "/"
      << buffer.getSize() << "," << buffer.getSize() << "/" << buffer.getMetadataSize() << "+"
      << buffer.getCapacity() << " R? " << (buffer.is_read_only() ? "y" : "n")
      << " W? "
      << (buffer.is_writing() ? "y" : "n") << std::flush;

      return ost;
    }


  } /* namespace io */
} /* namespace bliss */

#endif /* BUFFER_HPP_ */
