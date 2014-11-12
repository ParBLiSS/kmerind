/**
 * @file		Buffer.hpp
 * @ingroup bliss::io
 * @author	tpan
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

#include "concurrent/concurrent.hpp"   // ThreadSafety boolean constants
#include "utils/logging.h"

#include <xmmintrin.h>  // _mm_pause, instead of usleep.

namespace bliss
{
  namespace io
  {

  // TODO: double check end_ptr initialization.  require that end_ptr be greater than byteCount at all times, but how does that work with read lock/unlock?

    /**
     * @brief Thread safe/unsafe Memory Buffer, which is a fixed size allocated block of memory where data can be appended to, and read from via pointer.
     * @details this class uses unique_ptr internally to manage the memory, and supports only MOVE semantics.
     *
     * Thread Safety is enforced using atomic variables - writing to offsets that are atomically computed.
     * Memory ordering uses acquire/consume and release, and avoids seq_cst.
     *
     * Only when necessary, mutex lock is used. (in move constructor, move assignment operator, and in append function)
     * See append function for information about why mutex lock is needed.
     *
     * thread-safe and thread-unsafe versions can be move constructed/assigned from each other.
     *
     * life cycle:  pool acquire():  buffer unblock().  -> allow rw.
     *              application:  buffer block() -> allow ro
     *              pool release():  buffer clear().
     *
     * @tparam  ThreadSafety  controls whether this class is thread safe or not
     */
    template<bliss::concurrent::ThreadSafety ThreadSafety>
    class Buffer
    {
      /*
       * Declare Buffer<!ThreadSafety> class as a friend class, so we can reference its member functions and variables
       * in move constructor and assignment operators that moves data between instances with different thread safety.
       *
       * alternative is to use inheritance, which would require virtual functions that are potentially expensive.
       */
      friend class Buffer<!ThreadSafety>;
      template<bliss::concurrent::ThreadSafety TS>
      friend std::ostream & operator<<(std::ostream &os, const Buffer<TS>& p);


      protected:

        /// maximum capacity.  unrealistic to use size_t - can't possibly allocate.

      	/// internal data storage
        mutable uint8_t* start_ptr; // const, does not change
        mutable uint8_t* max_ptr;   // const, does not change

        volatile typename std::conditional<ThreadSafety, std::atomic<uint8_t*>, uint8_t*>::type curr_ptr;
        mutable typename std::conditional<ThreadSafety, std::atomic<uint8_t*>, uint8_t*>::type end_ptr;   // changes ONLY when read_lock, or when full
                                                                                                          // store min of all attempts.

        /// Reference Count of threads performing append update.  to ensure that all updates are done before used for sending.
        /// using signbit to indicate flushing/full  - NO MORE.
        volatile typename std::conditional<ThreadSafety, std::atomic<uint8_t*>, uint8_t*>::type byteCount;

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
        Buffer(Buffer<ThreadSafety>&& other, const std::lock_guard<std::mutex> &l)
      	  : start_ptr(other.start_ptr), max_ptr(other.max_ptr),
      	    curr_ptr((uint8_t*)(other.curr_ptr)), end_ptr((uint8_t*)(other.end_ptr)),
      	    byteCount((uint8_t*)(other.byteCount)) {

          //DEBUG("BUFFER private move contructor 1 called");

          other.start_ptr = nullptr;
          other.end_ptr =   nullptr;
          other.max_ptr =   nullptr;
          other.byteCount = nullptr;
          other.curr_ptr = max_ptr + 1;

        };

        /**
         * @brief Move constructor with a mutex lock on the move source object.
         * @details  This constructor is private and only use as constructor delegation target.
         * the source object is locked before this function is called and data moved to the newly constructed object.
         * This version copies between 2 objects with the DIFFERENT thread safety.
         *
         * @param other   the source Buffer
         * @param l       the mutex lock on the source Buffer.
         */
        Buffer(Buffer<!ThreadSafety>&& other, const std::lock_guard<std::mutex> &l)
          : start_ptr(other.start_ptr), max_ptr(other.max_ptr),
            curr_ptr((uint8_t*)(other.curr_ptr)), end_ptr((uint8_t*)(other.end_ptr)),
            byteCount((uint8_t*)other.byteCount) {

          other.start_ptr = nullptr;
          other.end_ptr =   nullptr;
          other.max_ptr =   nullptr;
          other.byteCount = nullptr;
          other.curr_ptr = max_ptr + 1;
        }

      public:
        /**
         * @brief Normal constructor.  Allocate and initialize memory with size specified as parameter.
         * @param _capacity   The maximum capacity of the Buffer in bytes
         */
        explicit Buffer(const uint32_t _capacity)
        {
          if (_capacity == 0)
            throw std::invalid_argument("Buffer constructor parameter capacity is given as 0");

          std::lock_guard<std::mutex> lock(mutex);
          start_ptr = new uint8_t[_capacity];

          memset(start_ptr, 0, _capacity * sizeof(uint8_t));
          // buffer empty
          byteCount = start_ptr;
          max_ptr = start_ptr + _capacity;
          // buffer blocked.
          end_ptr = start_ptr;
          curr_ptr = max_ptr + 1;   // block from insertion.

          //printf("NEW BUFF:  start: %p, max_ptr %p, curr_ptr %p, end %p, byteCount %p\n", start_ptr, max_ptr, (uint8_t*)curr_ptr, (uint8_t*)end_ptr, (uint8_t*)byteCount);
        };

        /**
         * @brief reate a new Buffer using the memory specified along with allocated memory's size.
         * @param _data   The pointer/address of preallocated memory block
         * @param count   The size of preallocated memory block
         */
        Buffer(void* _data, const int32_t count)  {
          if (count == 0)
            throw std::invalid_argument("Buffer constructor parameter count is given as 0");
          if (_data == nullptr)
            throw std::invalid_argument("Buffer constructor parameter _data is given as nullptr");

          // buffer full
          std::lock_guard<std::mutex> lock(mutex);

          start_ptr = reinterpret_cast<uint8_t*>(_data);
          max_ptr = start_ptr + count;
          byteCount = start_ptr + count;
          end_ptr = start_ptr + count;
          // and blocked.
          curr_ptr = max_ptr + 1;
        }

        /**
         * @brief Destructor.  deallocate memory manually.
         */
        virtual ~Buffer() {

           //printf("deleting.  %p being cleared.\n", start_ptr);

          std::lock_guard<std::mutex> lock(mutex);

          max_ptr = nullptr;
          byteCount = nullptr;
          // blocked.
          end_ptr = nullptr;
          curr_ptr = max_ptr + 1;


          if (start_ptr != nullptr) {
            delete [] start_ptr;
            start_ptr = nullptr;
          }
        };

        /**
         * @brief Move constructs from a Buffer with the SAME ThreadSafety property
         * @details  internal data memory moved by std::unique_ptr semantics.
         * Thread Safe always, using approach from http://www.justsoftwaresolutions.co.uk/threading/thread-safe-copy-constructors.html
         *
         * Constructs a mutex lock and then delegates to another constructor.
         *
         *
         * @param other   Source object to move
         */
        explicit Buffer(Buffer<ThreadSafety>&& other) : Buffer<ThreadSafety>(std::move(other), std::lock_guard<std::mutex>(other.mutex) ) {};
        //explicit Buffer(Buffer<ThreadSafety>&& other) = delete;


        /**
         * @brief Move constructs from a Buffer with the DIFFERENT ThreadSafety property
         * @details internal data memory moved by std::unique_ptr semantics.
         * Thread Safe always, using approach from http://www.justsoftwaresolutions.co.uk/threading/thread-safe-copy-constructors.html
         *
         * Constructs a mutex lock and then delegates to another constructor.
         *
         * @param other   Source object to move
         */
        explicit Buffer(Buffer<!ThreadSafety>&& other) : Buffer<ThreadSafety>(std::move(other), std::lock_guard<std::mutex>(other.mutex) ) {};

        /**
         * @brief Move assignment operator, between Buffers of the SAME ThreadSafety property.
         * @details  Internal data memory moved by std::unique_ptr semantics.
         * The move is done in a thread safe way always.
         *
         * @param other     Source Buffer to move
         * @return          target Buffer reference
         */
        Buffer<ThreadSafety>& operator=(Buffer<ThreadSafety>&& other) {
          //DEBUG("BUFFER public move assignment op 1 called");

          if (this->start_ptr != other.start_ptr) {

            std::unique_lock<std::mutex> myLock(mutex, std::defer_lock),
                                         otherLock(other.mutex, std::defer_lock);
            std::lock(myLock, otherLock);

            if (start_ptr != nullptr) {
              delete [] start_ptr;
              start_ptr = nullptr;
            }

            /// move the internal memory.

            start_ptr = other.start_ptr;
            byteCount = (uint8_t*)(other.byteCount);
            end_ptr = (uint8_t*)(other.end_ptr);
            max_ptr = other.max_ptr;
            curr_ptr = (uint8_t*)(other.curr_ptr);

            other.end_ptr = nullptr;
            other.start_ptr = nullptr;
            other.max_ptr =   nullptr;
            other.byteCount = nullptr;
            other.curr_ptr = max_ptr + 1;
          }
          return *this;
        }
//        Buffer<ThreadSafety>& operator=(Buffer<ThreadSafety>&& other) = delete;

        /**
         * @brief Move assignment operator between Buffers of the DIFFERENT ThreadSafety property.
         * @details  Internal data memory moved by std::unique_ptr semantics.
         * The move is done in a thread safe way always.
         *
         * @param other     Source Buffer to move
         * @return          target Buffer reference
         */
        Buffer<ThreadSafety>& operator=(Buffer<!ThreadSafety>&& other) {
          //DEBUG("BUFFER public move assignment op 2 called");

          if (this->start_ptr != other.start_ptr) {

            std::unique_lock<std::mutex> myLock(mutex, std::defer_lock),
                                         otherLock(other.mutex, std::defer_lock);
            std::lock(myLock, otherLock);

            if (start_ptr != nullptr) {
              delete [] start_ptr;
              start_ptr = nullptr;
            }

            /// move the internal memory.
            start_ptr = other.start_ptr;
            byteCount = (uint8_t*)(other.byteCount);
            end_ptr = (uint8_t*)(other.end_ptr);
            max_ptr = other.max_ptr;
            curr_ptr = (uint8_t*)(other.curr_ptr);

            other.end_ptr = nullptr;
            other.start_ptr = nullptr;
            other.max_ptr =   nullptr;
            other.byteCount = nullptr;
            other.curr_ptr = max_ptr + 1;

          }
          return *this;
        }



        /// remove copy constructor and copy assignement operators.
        explicit Buffer(const Buffer<ThreadSafety>& other) = delete;
        /// remove copy constructor and copy assignement operators.
        Buffer<ThreadSafety>& operator=(const Buffer<ThreadSafety>& other) = delete;

        Buffer() = delete;


        /**
         * @brief get the current size of the Buffer.
         * @tparam TS Choose thread safe vs unsafe implementation. defaults to same as the parent class.
         * @return    current size
         */
        const int64_t getSize() const {
          std::atomic_thread_fence(std::memory_order_seq_cst);

          return (uint8_t*)end_ptr - start_ptr;

        }

      protected:

        /**
         * @brief get the current size of the Buffer.
         * @tparam TS Choose thread safe vs unsafe implementation. defaults to same as the parent class.
         * @return    current size
         */
        const int64_t getApproximateSize() const {
          std::atomic_thread_fence(std::memory_order_seq_cst);

          return (uint8_t*)curr_ptr - start_ptr;
        }

        const int64_t getWrittenSize() const {
          std::atomic_thread_fence(std::memory_order_seq_cst);

          return (uint8_t*)byteCount - start_ptr;
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
        template <typename T>
        operator T*() const {
          return reinterpret_cast<T*>(start_ptr);
        }


        /**
         * @brief get the capacity of the buffer.
         * @return    maximum capacity of buffer.
         */
        const int64_t getCapacity() const {
          return max_ptr - start_ptr;
        }


        /**
         * @brief Checks if a buffer is empty.
         * @note The return value is not precise - between getSize and return, other threads may have modified the size.
         * For "isEmpty" error is infrequent.
         *
         * @return    true if the buffer is empty, false otherwise.
         */
        const bool isEmpty() const {
          return (getWrittenSize() == 0);
        }


        // TODO:  make this more clear.
        // need 1. flag to indicate no more new writers
        //      2. wait for current writers to finish
        //      3. 3 phases:  writing, flushing, reading.  -> write_lock, wirte_unlock, flush_begin, read_lock, read_unlock
        //          write_lock (counter ++)
        //          wirte_unlock (counter --)
        //          flush (MSB = 1)
        //          read_lock ( wait for MSB==1 && counter == 0)
        //          read_unlock ( MSB = 0)
        // state checks:  is_writing (MSB==0|1 && counter > 0)
        //                is_flushing (MSB==1 && counter> 0)
        //                is_reading ( MSB==1 && counter == 0)

        inline bool is_writing() const {
          std::atomic_thread_fence(std::memory_order_seq_cst);
          return (uint8_t*)byteCount < (uint8_t*)end_ptr;
        }

        inline bool is_full() const {
          std::atomic_thread_fence(std::memory_order_seq_cst);
          return (uint8_t*)curr_ptr > max_ptr;
        }

        inline bool is_flushing() const {
        	return (is_full()) && (is_writing());
        }

        inline bool is_reading() const {
          return (is_full()) && (!is_writing());
        }


      protected:
        /// increment lock only if flush bit is not set.  else do nothing.  cannot increment in all cases then decrement if flush bit is set - with threading could cause race condition.
    	/// reserves a place for storing the thread's new data.
        ///  this is designed so only 1 pointer needs to be atomically incremented.  to determine if all threads wrote to the buffer, have a pointer = start + completed count so far
        ///  chase the current end pointer.
        /// if a buffer is full (curr_ptr > max_ptr == start + capacity)
        template<bliss::concurrent::ThreadSafety TS = ThreadSafety>
        inline typename std::enable_if<TS, uint8_t*>::type write_lock(const uint32_t count) {

          std::atomic_thread_fence(std::memory_order_seq_cst);
          uint8_t* ptr = curr_ptr.fetch_add(count);  // always increment

          if (ptr > max_ptr) {  // full.  previously saved to end_ptr
            return nullptr;
          } else {
            if (ptr + count > max_ptr) {  // just filled.  set end_ptr.  curr_ptr is already set.

            	// min reduction to get the smallest pointer as end pointer
                std::atomic_thread_fence(std::memory_order_seq_cst);
              uint8_t* end = (uint8_t*)end_ptr;
              bool stop = false;            // use CAS to set the minimum. same as end_ptr = (end > ptr ? ptr : end);
              while (end > ptr && !stop) {
                stop = end_ptr.compare_exchange_weak(end, ptr);
              }
            } else {
              // normal.
            }
          }
          return ptr;
        }

        template<bliss::concurrent::ThreadSafety TS = ThreadSafety>
        inline typename std::enable_if<!TS, uint8_t*>::type write_lock(const uint32_t count) {
          uint8_t* ptr = curr_ptr;
          curr_ptr += count;
          if (ptr > max_ptr) {
            return nullptr;  //if ptr == max_ptr, buffer must be FILLED in this cycle.  need to process it as if not full yet.
          } else {
            if (ptr + count > max_ptr) {

              if (end_ptr > ptr) end_ptr = ptr;
            } else {
              // normal write.
            }
          }
        	return ptr;
        }

        inline void write_unlock(const int count) {
          std::atomic_thread_fence(std::memory_order_seq_cst);

        	byteCount += count;
        }



      public:

        //===== read lock trumps write.
        /// read lock prevents future writes.  read lock can be turned on while there are writes in progress.  once read lock is on, write_lock will returns nullptr.
        // thread safe - multiple threads calling is the same as 1 thread calling
        // purpose of this method is to swap the end_ptr and curr_ptr, so that curr_ptr is at max+1, and end_ptr is the smallest of all locking thread's pointers
        template<bliss::concurrent::ThreadSafety TS = ThreadSafety>
        inline typename std::enable_if<TS, void>::type read_lock(uint8_t* _curr = nullptr) {  // not using xor since it toggles and write_unlock is not going to change the sign bit (so we don't need xor)
          // disable curr_ptr and
          uint8_t* end = (uint8_t*)end_ptr;
          uint8_t* curr = curr_ptr.exchange(max_ptr + 1);

          std::atomic_thread_fence(std::memory_order_seq_cst);  // make sure all writes are done.

          curr = (_curr ? _curr : curr);
          bool stop = false;


          while (end > curr && !stop) {
            stop = end_ptr.compare_exchange_weak(end, curr);  // get smallest, indicating the first one to reach read_lock
          }
          std::atomic_thread_fence(std::memory_order_seq_cst);  // make sure all writes are done.
        }
        template<bliss::concurrent::ThreadSafety TS = ThreadSafety>
        inline typename std::enable_if<!TS, void>::type read_lock(uint8_t* _curr = nullptr) {  // not using xor since it toggles and write_unlock is not going to change the sign bit (so we don't need xor)
          // disable curr_ptr and
          uint8_t* curr = (_curr ? _curr : curr_ptr);
          curr_ptr = max_ptr + 1;

          if (end_ptr > curr) {
            end_ptr = curr;
          }
        }

        /// read unlock allows future writes.  read unlock blocks until all pending writes are done before .
        template<bliss::concurrent::ThreadSafety TS = ThreadSafety>
        inline typename std::enable_if<TS, void>::type read_unlock() {
          // put the end_ptr into curr if it hasn't been done.
          uint8_t* curr = (uint8_t*)curr_ptr;
          uint8_t* end = end_ptr.exchange(max_ptr+1);
          std::atomic_thread_fence(std::memory_order_seq_cst);  // make sure all writes are done.

          bool stop = false;
          while (curr > max_ptr && ! stop) {
            stop = curr_ptr.compare_exchange_weak(curr, end);
          }
          std::atomic_thread_fence(std::memory_order_seq_cst);  // make sure all writes are done.

        }
        template<bliss::concurrent::ThreadSafety TS = ThreadSafety>
        inline typename std::enable_if<!TS, void>::type read_unlock() {
          // put the end_ptr into curr if it hasn't been done.
          if (curr_ptr > max_ptr)
            curr_ptr = end_ptr;
          // then reset end_ptr to max
          end_ptr = max_ptr + 1;
        }

      protected:
        void wait_for_writes() {
          while (is_writing()) _mm_pause();
        }

        // has to ensure flush bit is set else return
        void wait_for_flush() {
            std::atomic_thread_fence(std::memory_order_seq_cst);  // make sure all writes are done.

          if ((uint8_t*)curr_ptr <= max_ptr)  // check flush bit.
            throw std::logic_error("ERROR: wait for flush called but flush bit is not set");
          std::atomic_thread_fence(std::memory_order_seq_cst);  // make sure all writes are done.

          while (is_writing()) {
            _mm_pause();  // all threads attempting to lock read are waiting until writing is done.
            if ((uint8_t*)curr_ptr <= max_ptr)
              throw std::logic_error("ERROR: wait for flush called but flush bit is not set");
            std::atomic_thread_fence(std::memory_order_seq_cst);  // make sure all writes are done.
          }
        }

      public:

        // ONLY CALL IF BUFFER IS NOT FULL AND WE ARE FORCING THE FLUSH
        inline void flush_and_set_size() {
            //std::cout << "DEBUG: read unlocked buffer: " << *(this) << std::endl << std::flush;

          read_lock<ThreadSafety>();
          //std::cout << "DEBUG: read locked buffer: " << *(this) << std::endl << std::flush;

          std::atomic_thread_fence(std::memory_order_seq_cst);  // make sure all writes are done.

          // now everyone wait for the writes to finish during flush period.
          wait_for_flush();
         // std::cout << "DEBUG: flushed buffer: " << *(this) << std::endl << std::flush;

          std::atomic_thread_fence(std::memory_order_seq_cst);  // make sure all writes are done.
        }






        /**
         * @brief Clears the buffer. (set the size to 0, leaving the capacity and the memory allocation intact)
         * @note
         * const because the caller will have a const reference to the buffer.
         */
        void clear() {
          std::lock_guard<std::mutex> lock(mutex);
          // TODO: remove after finish debugging.
          memset(start_ptr, 0, getCapacity());
          // blocked
          curr_ptr = max_ptr +1;
          end_ptr = start_ptr;
          byteCount = start_ptr;
        }


        /**
         * @brief Append data to the buffer, THREAD SAFE.
         * @details The function updates the current occupied size of the Buffer
         *  and memcopies the supplied data into the internal memory block.
         *
         *  This is the THREAD SAFE version using MUTEX LOCK.
         *
         * Semantically, a "false" return value means there is not enough room for the new data, not that the buffer is flushing.
         *
         * method is const because the caller will have const reference to the Buffer.
         *
         * NOTE: we can't use memory ordering alone.  WE have to lock with a mutex or use CAS in a loop
         * Suppose count can be different for each call, and there are 3 threads calling with count1, count2, and count3 in that order.
         *  further suppose s+count1 > capacity, count1 + count2 < count1, count1
         * then it could happen that
         *    thread1 fetch_add s' to s+count1 > capacity,
         *        thread 2 fetch_add s' to s + count1+count2 > capacity,
         *    thread1 fetch_sub s' to s+count2 < capacity,
         *            thread3 fetch_add s' to s+count2 + count3 < capacity,
         *        thread2 fetch_sub s' to s+count3,
         *            thread3 writes to s' at s+count2.
         * In the end thread 3 has written the data between s+count2 and s+count2+count3,  while s has been set to s+count3.
         * To avoid this situation, this set of calls has to be locked with mutex (else we need to use CAS in a loop)
         *
         *
         * @tparam TS Choose thread safe vs unsafe implementation. defaults to same as the parent class.
         * @param[in] _data   pointer to data to be copied into the Buffer..  this SHOULD NOT be shared between threads
         * @param[in] count   number of bytes to be copied into the Buffer
         * @return            unsigned char, bit 0 indicated whether operation succeeded, bit 1 indicating whether buffer swap is needed.  single primitive type is better as it can be returned in 1 atomic op.
         */
        template<bliss::concurrent::ThreadSafety TS = ThreadSafety>
        typename std::enable_if<TS, unsigned int>::type append(const void* _data, const uint32_t count) {

          if (count == 0 || _data == nullptr)
        	  throw std::invalid_argument("_data is nullptr");

          //if (curr_ptr > max_ptr) printf("write_lock after full.  data is %p\n", start_ptr);


          uint8_t* ptr = write_lock<TS>(count);
          // issue here:
          //    the writerCount may be delayed relative to ptr computation.
          // because ptr is monotonically increasing, if we could couple writeCount increment to that then we
          // can guarantee that at swap time, we have the maximum writeCount we need.
          // however,  writerCount is not well coupled, and can be delayed, so wait_for_flush is not effective.
          // POSSIBLE SOLUTION:  use a second pointer/counter that increments by "count" not 1.  when that is equal to end_ptr, then we have completed.
          //    essentially, curr_ptr is the reserved count, and this second ptr is the actual written count.

          // After the F&A, there are 3 types of threads:  target write area is
          //  1. completely within buffer:  these threads should proceed to memcpy with the local ptr var
          //  2. completely outside buffer: these threads will not memcpy.  they all disabled buffer, and
          //     only 1 thread from 2) or 3) should swap in a new buffer atomically.
          //      - if curr buffer is not disabled (already swapped), then don't swap further.
          //  3. crossing buffer boundary: this thread will not memcpy.  disabled.  and it should retract the pointer advance.



          // if fails, then already flushing
          if (ptr == nullptr) {  // write lock returns nullptr if ptr points to outside of the range.
            //printf("flushing     shared pointer %p, count %d, start %p, max_ptr %p, size %d, writers %d\n", (uint8_t*)pointer, count, start_ptr, max_ptr, (int32_t)size, (CountType)writerCount);

        	//printf(" FULL.  end_ptr = %p, curr_ptr = %p. writerCount = %d.\n", end_ptr.load(), curr_ptr.load(), writerCount.load());

        	  // already read locked.

            return 0x0;  // full and not swapping.
          } else {  // ptr starts inside allocated buffer range.

            if ((ptr + count) > this->max_ptr) { // thread that filled the buffer


              //// prvent from being swapped out (prevent move constructor and assignment operator)
              std::unique_lock<std::mutex> lock(mutex);


              if (ptr > this->max_ptr) {
                fprintf(stdout, "FAIL: ptr is %p, larger than maxptr %p\n", ptr, this->max_ptr);
                std::cout << " buffer: " << *this << std::endl << std::flush;
              }

              if (this->start_ptr == nullptr)
                std::cout << "DEBUG FAIL null data ptr.  swap requested. buffer: " << *this << std::endl << std::flush;

              std::atomic_thread_fence(std::memory_order_seq_cst);  // make sure all writes are done.
              if ((uint8_t*)curr_ptr <= max_ptr) {
                fprintf(stdout, "FAIL before flush: flush bit should be set\n");
                std::cout << " buffer: " << *this << std::endl << std::flush;
              }


              wait_for_flush();

              if (! is_reading()) {
                fprintf(stdout, "FAIL after flush: all writes should be done and flush set. curr %p, max %p, byteCount %p\n", (uint8_t*)curr_ptr, max_ptr, (uint8_t*)byteCount);
                std::cout << " buffer: " << *this << std::endl << std::flush;
              }

              if ((getSize() / sizeof(int)) != (getCapacity() / sizeof(int))) {
                fprintf(stdout , "FAIL IN BUFFER:  NOT %lu elements. got %ld.", (getCapacity() / sizeof(int)), getSize() / sizeof(int));
                std::cout << "FAIL BUF: " << *this << std::endl << std::flush;

                return 0x4;
              }

              //printf("prep for swap done.  data is %p\n", start_ptr);


              lock.unlock();
              return 0x2;  // was not full, now full.

            } else { // can insert.  may make buffer full
            	//printf(" COPY.  end_ptr = %p, curr_ptr = %p. writerCount = %d.\n", end_ptr.load(), curr_ptr.load(), writerCount.load());
            	std::atomic_thread_fence(std::memory_order_seq_cst);  // unlock only after memcpy.

              if (ptr < start_ptr || ptr >= max_ptr) printf("ERROR: ptr %p is outside of the buffer range %p to %p\n", ptr, start_ptr, max_ptr);

              if ((uint8_t*)curr_ptr > max_ptr && ptr >= max_ptr) printf("ERROR: writing after full at %p.  data is %p\n", ptr, start_ptr);


              // write
              std::memcpy(ptr, _data, count);
              std::atomic_thread_fence(std::memory_order_seq_cst);  // unlock only after memcpy.
              write_unlock(count);  // all full buffers lock the read and unlock the writer
              std::atomic_thread_fence(std::memory_order_seq_cst);  // unlock only after memcpy.

              return 0x1;   // not full, successfully inserted.
            }

          }

        }

        /**
         * @brief Append data to the buffer, THREAD UNSAFE.
         * @details  The function updates the current occupied size of the Buffer
         * and memcopies the supplied data into the internal memory block.
         *
         *  This is the THREAD UNSAFE version.
         *
         * Semantically, a "false" return value means there is not enough room for the new data, not that the buffer is flushing.
         *
         * method is const because the caller will have const reference to the Buffer.
         *
         * @tparam TS Choose thread safe vs unsafe implementation. defaults to same as the parent class.
         * @param[in] _data   pointer to data to be copied into the Buffer
         * @param[in] count   number of bytes to be copied into the Buffer
         * @return            unsigned int, bit 0 indicated whether operation succeeded, bit 1 indicating if swap is needed.
         */
        template<bliss::concurrent::ThreadSafety TS = ThreadSafety>
        typename std::enable_if<!TS, unsigned int >::type append(const void* _data, const uint32_t count) {

          if (count == 0 || _data == nullptr)
            throw std::invalid_argument("_data has 0 count or is _nullptr");


          uint8_t* ptr = write_lock<TS>(count);


          if (ptr == nullptr) {

            return 0x0;
          }

          if ((ptr + count) > max_ptr) {  // append that filled the buffer

            wait_for_flush();

            return 0x2;
          } else {
            if (ptr < start_ptr || ptr >= max_ptr) printf("ERROR: ptr %p is outside of the buffer range %p to %p\n", ptr, start_ptr, max_ptr);

            std::memcpy(ptr, _data, count);
            write_unlock(count);

            return 0x1;
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
    template<bliss::concurrent::ThreadSafety ThreadSafety>
    std::ostream& operator<<(std::ostream& ost, const Buffer<ThreadSafety> & buffer)
    {
      ost << "THREAD "<< (ThreadSafety ? "SAFE" : "UNSAFE") << " BUFFER: data_ptr/end_ptr=" << static_cast <const void *>(buffer.start_ptr) << "/" << static_cast <const void *>(buffer.end_ptr)
        << " currptr/maxptr=" << static_cast <const void *>((uint8_t*)(buffer.curr_ptr)) << "/" << static_cast <const void *>(buffer.max_ptr)
        << " byteCount=" << static_cast <const void *>((uint8_t*)(buffer.byteCount))
        << " approx,size/cap=" << buffer.getApproximateSize() << "," << buffer.getSize() << "/" << buffer.getCapacity()
        << " R? " << (buffer.is_reading() ? "y" : "n") << " F? " << (buffer.is_flushing() ? "y" : "n") << " W? " << (buffer.is_writing() ? "y" : "n") << std::flush;


      return ost;
    }

  } /* namespace io */
} /* namespace bliss */

#endif /* BUFFER_HPP_ */
