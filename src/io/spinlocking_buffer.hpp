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

#include "xmmintrin.h"  // _mm_pause, instead of usleep.

namespace bliss
{
  namespace io
  {

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


      public:
      typedef int   CountType;

      protected:

        /// maximum capacity.  unrealistic to use size_t - can't possibly allocate.
        mutable uint32_t capacity;

        /// internal data storage
        mutable uint8_t* start_ptr; // const, does not change
        mutable uint8_t* max_ptr;   // const, does not change

        volatile typename std::conditional<ThreadSafety, std::atomic<uint8_t*>, uint8_t*>::type curr_ptr;
        mutable typename std::conditional<ThreadSafety, std::atomic<uint8_t*>, uint8_t*>::type end_ptr;   // changes ONLY when read_lock, or when full
                                                                                                          // store min of all attempts.

        static constexpr uint8_t* MAX_PTR_ADDR = reinterpret_cast<uint8_t*>(UINTPTR_MAX);


        /// Reference Count of threads performing append update.  to ensure that all updates are done before used for sending.
        /// using signbit to indicate flushing/full  - NO MORE.
        volatile typename std::conditional<ThreadSafety, std::atomic<CountType>, CountType>::type writerCount;    // track how many are appending.

        mutable std::atomic_flag spinlock = ATOMIC_FLAG_INIT;

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
      	  : capacity(other.capacity), start_ptr(other.start_ptr), max_ptr(other.max_ptr),
      	    curr_ptr((uint8_t*)(other.curr_ptr)), end_ptr((uint8_t*)(other.end_ptr)),  writerCount(other.writerCount)
      	    {

          //DEBUG("BUFFER private move contructor 1 called");

          other.capacity = 0;
          other.start_ptr = nullptr;
          other.end_ptr = nullptr;
          other.max_ptr = nullptr;
          other.curr_ptr = MAX_PTR_ADDR;
          other.writerCount = 0;
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
          : capacity(other.capacity), start_ptr(other.start_ptr), max_ptr(other.max_ptr),
            curr_ptr((uint8_t*)(other.curr_ptr)), end_ptr((uint8_t*)(other.end_ptr)), writerCount(other.writerCount) {

          other.capacity = 0;
          other.start_ptr = nullptr;
          other.end_ptr = nullptr;
          other.max_ptr = nullptr;
          other.curr_ptr = MAX_PTR_ADDR;
          other.writerCount = 0;
        }

      public:
        /**
         * @brief Normal constructor.  Allocate and initialize memory with size specified as parameter.
         * @param _capacity   The maximum capacity of the Buffer in bytes
         */
        explicit Buffer(const uint32_t _capacity) : capacity(_capacity), start_ptr(new uint8_t[_capacity]),
          writerCount(0)
        {
          if (_capacity == 0)
            throw std::invalid_argument("Buffer constructor parameter capacity is given as 0");

          memset(start_ptr, 0, _capacity * sizeof(uint8_t));
          end_ptr = start_ptr;
          max_ptr = start_ptr + _capacity;
          curr_ptr = MAX_PTR_ADDR;   // block from insertion.

          //printf("NEW BUFF:  start: %p, max_ptr %p, curr_ptr %p, size %d, writers %d\n", start_ptr, max_ptr, (uint8_t*)curr_ptr, (int32_t)size, (CountType)writerCount);
        };

        /**
         * @brief reate a new Buffer using the memory specified along with allocated memory's size.
         * @param _data   The pointer/address of preallocated memory block
         * @param count   The size of preallocated memory block
         */
        Buffer(void* _data, const int32_t count) : capacity(count), start_ptr(static_cast<uint8_t*>(_data)),
            writerCount(0) {
          if (count == 0)
            throw std::invalid_argument("Buffer constructor parameter count is given as 0");
          if (_data == nullptr)
            throw std::invalid_argument("Buffer constructor parameter _data is given as nullptr");

          end_ptr = start_ptr + count;
          max_ptr = start_ptr + count;
          curr_ptr = MAX_PTR_ADDR;
        }

        /**
         * @brief Destructor.  deallocate memory manually.
         */
        virtual ~Buffer() {
          std::lock_guard<std::mutex> lock(mutex);
          if (start_ptr != nullptr) {
            delete [] start_ptr;
            start_ptr = nullptr;
          }
          end_ptr = nullptr;
          max_ptr = nullptr;
          curr_ptr = MAX_PTR_ADDR;
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

            /// move the internal memory.

            capacity = other.capacity;
            start_ptr = other.start_ptr;
            end_ptr = (uint8_t*)(other.end_ptr);
            max_ptr = other.max_ptr;
            curr_ptr = (uint8_t*)(other.curr_ptr);
            writerCount = CountType(other.writerCount);

            other.capacity = 0;
            other.end_ptr = nullptr;
            other.start_ptr = other.max_ptr = nullptr;
            other.curr_ptr = MAX_PTR_ADDR;
            other.writerCount = 0;
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

            /// move the internal memory.
            capacity = other.capacity;
            start_ptr = other.start_ptr;
            end_ptr = (uint8_t*)(other.end_ptr);
            max_ptr = other.max_ptr;
            curr_ptr = (uint8_t*)(other.curr_ptr);
            writerCount = CountType(other.writerCount);

            other.capacity = 0;
            other.end_ptr = nullptr;
            other.start_ptr = other.max_ptr = nullptr;
            other.curr_ptr = MAX_PTR_ADDR;
            other.writerCount = 0;

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
          return end_ptr - start_ptr;

        }

      protected:

        /**
         * @brief get the current size of the Buffer.
         * @tparam TS Choose thread safe vs unsafe implementation. defaults to same as the parent class.
         * @return    current size
         */
        int64_t getApproximateSize() const {
          return curr_ptr - start_ptr;
        }

      public:

        // TODO:  make this more clear.
        // need 1. flag to indicate no more new writers
        //      2. wait for current writers to finish
        //      3. 3 phases:  writing, flushing, reading.  -> write_lock, wirte_unlock, flush_begin, read_lock, read_unlock
        //          write_lock (counter ++)
        //          wirte_unlock (counter --)
        //          flush_begin (MSB = 1)
        //          read_lock ( wait for MSB==1 && counter == 0)
        //          read_unlock ( MSB = 0)
        // state checks:  is_writing (MSB==0|1 && counter > 0)
        //                is_flushing (MSB==1 && counter> 0)
        //                is_reading ( MSB==1 && counter == 0)
        //   if unsigned type, can use highest bit to manage this, and reinterpret_cast to sign type for comparison.

        inline bool is_writing() const {
        	// not flushing (MSB == 0 | 1) and counter > 0
          //std::atomic_thread_fence(std::memory_order_seq_cst);
        	return writerCount > 0;
        }

        inline bool is_flushing() const {
        	// flushing (MSB == 1) and some writes (lower bits > 0) == same as value in range (lowest(), 0), exclusive of both ends;
          //std::atomic_thread_fence(std::memory_order_seq_cst);
          bool result = false;
          while (spinlock.test_and_set());
          result = (curr_ptr == MAX_PTR_ADDR) && (writerCount > 0);
          spinlock.clear();

        	return result;
        	// reinterpret_cast does not generate any op codes.
        }

        // TODO
        inline bool is_reading() const {
        	// flushing (MSB == 1) and no write.
          //std::atomic_thread_fence(std::memory_order_seq_cst);
          bool result = false;
          while (spinlock.test_and_set());
          result = (curr_ptr == MAX_PTR_ADDR) && (writerCount == 0);
          spinlock.clear();

          return result;
        }


        /// increment lock only if flush bit is not set.  else do nothing.  cannot increment in all cases then decrement if flush bit is set - with threading could cause race condition.
    	// increment in unsigned and 2's complement follow the same directionality when signbit is excluded.
    	// so this means that sign bit can change by another thread and the lower bits remain valid, even is reordered between threads.
    	// e.g. 100 + 1 = 101. 101 & 011 = 001.  001 + 1 = 010. 010 - 1 = 001, 001 - 1 = 000
    	// so we can atomically increment, then compare to lowest() (100), and undo or not.
    	// we rely on 2 things: 1. no 0x7FFFFFFF threads in single node, so no chance of signbit being changed from overflow.
        // 2. i and i+1 have the same lower bit patterns whether i is negative or positive.
        // return bool indicating success or failure. since we increment only if not flushing, we need to known success/fail
        // to avoid extra unlocks.
       // no chance of overflow with addition unless we have machines with 2B threads.  goal is to preserve sign bit during this op.
        template<bliss::concurrent::ThreadSafety TS = ThreadSafety>
        inline typename std::enable_if<TS, uint8_t*>::type write_lock(const uint32_t count) {


          while (spinlock.test_and_set());

          uint8_t* ptr = (uint8_t*)curr_ptr;
          if (ptr == MAX_PTR_ADDR) {
            ptr = nullptr;
          } else if (ptr > max_ptr) {
            curr_ptr += count;
          
            ptr = nullptr;  //if ptr == max_ptr, buffer must be FILLED in this cycle.  need to process it as if not full yet.
          } else {
            curr_ptr += count;
            ++writerCount; // important:  only increment if not full.
          }
          spinlock.clear();
          return ptr;
        }
        template<bliss::concurrent::ThreadSafety TS = ThreadSafety>
        inline typename std::enable_if<!TS, uint8_t*>::type write_lock(const uint32_t count) {
          uint8_t* ptr = (uint8_t*)curr_ptr;
          if (ptr == MAX_PTR_ADDR) {
            ptr = nullptr;
          } else if (ptr > max_ptr) {
            curr_ptr += count;

            ptr = nullptr;  //if ptr == max_ptr, buffer must be FILLED in this cycle.  need to process it as if not full yet.
          } else {
            curr_ptr += count;
            ++writerCount; // important:  only increment if not full.
          }
          return ptr;
        }




        /// cannot decrement past 0 ( with or without sign bit).  should be called with knowledge of whether try_lock_write succeeded or not.
        /// throws exception if there are no writers
        // goal is to preserve sign bit for the entire duration of the op.
        template<bliss::concurrent::ThreadSafety TS = ThreadSafety>
        inline typename std::enable_if<TS, void>::type write_unlock() {
          while (spinlock.test_and_set());
          if (writerCount <= 0) {
            spinlock.clear();
            throw std::logic_error("ERROR: write_unlock thread unsafe on 0 writers");
          }
          --writerCount;
          spinlock.clear();
        }

        template<bliss::concurrent::ThreadSafety TS = ThreadSafety>
        inline typename std::enable_if<!TS, void>::type write_unlock() {
          if (writerCount <= 0)
            throw std::logic_error("ERROR: write_unlock thread unsafe on 0 writers");
          --writerCount;
        }





        //===== read lock trumps write.
        /// read lock prevents future writes.  read lock can be turned on while there are writes in progress.  once read lock is on, lock_write will fail.
        // thread safe - multiplr threads calling is the same as 1 thread calling
        template<bliss::concurrent::ThreadSafety TS = ThreadSafety>
        inline typename std::enable_if<TS, void>::type read_lock(uint8_t* _curr = nullptr) {  // not using xor since it toggles and write_unlock is not going to change the sign bit (so we don't need xor)

          while (spinlock.test_and_set());

          uint8_t* curr = (_curr ? _curr : (uint8_t*)curr_ptr);
          std::atomic_thread_fence(std::memory_order_seq_cst);
          curr_ptr = MAX_PTR_ADDR;

          if (end_ptr > curr) {
            end_ptr = curr;
          }
          spinlock.clear();
        }
        template<bliss::concurrent::ThreadSafety TS = ThreadSafety>
        inline typename std::enable_if<!TS, void>::type read_lock(uint8_t* _curr = nullptr) {  // not using xor since it toggles and write_unlock is not going to change the sign bit (so we don't need xor)
          // disable curr_ptr and
          uint8_t* curr = (_curr ? _curr : (uint8_t*)curr_ptr);
          curr_ptr = MAX_PTR_ADDR;

          if (end_ptr > curr) {
            end_ptr = curr;
          }
        }


        /// read unlock allows future writes.  read unlock blocks until all pending writes are done before .
        template<bliss::concurrent::ThreadSafety TS = ThreadSafety>
        inline typename std::enable_if<TS, void>::type read_unlock() {
          while (spinlock.test_and_set());
          if (curr_ptr == MAX_PTR_ADDR) {
            curr_ptr = (uint8_t*)end_ptr;
            end_ptr = MAX_PTR_ADDR;
          }
          spinlock.clear();
        }
        template<bliss::concurrent::ThreadSafety TS = ThreadSafety>
        inline typename std::enable_if<!TS, void>::type read_unlock() {
          // put the end_ptr into curr if it hasn't been done.
          if (curr_ptr == MAX_PTR_ADDR) {
            curr_ptr = end_ptr;
            end_ptr = MAX_PTR_ADDR;
          }
        }

        void wait_for_writes() {
          while (is_writing()) _mm_pause();
        }

        // has to ensure flush bit is set else return
        void wait_for_flush() {
          if (curr_ptr != MAX_PTR_ADDR)  // check flush bit.
            throw std::logic_error("ERROR: wait for flush called but flush bit is not set");
          //std::atomic_thread_fence(std::memory_order_seq_cst);  // make sure all writes are done.

          while(writerCount > 0 ) {
            _mm_pause();  // all threads attempting to lock read are waiting until writing is done.
            if (curr_ptr != MAX_PTR_ADDR)
              throw std::logic_error("ERROR: wait for flush called but flush bit is not set");
            ///std::atomic_thread_fence(std::memory_order_seq_cst);  // make sure all writes are done.
          }
        }


        // ONLY CALL IF BUFFER IS NOT FULL AND WE ARE FORCING THE FLUSH
        inline void flush_and_set_size() {
//            std::cout << "DEBUG: read unlocked buffer: " << *(this) << std::endl << std::flush;

          read_lock<ThreadSafety>();
//          std::cout << "DEBUG: read locked buffer: " << *(this) << std::endl << std::flush;

//          std::atomic_thread_fence(std::memory_order_seq_cst);  // make sure all writes are done.

          // now everyone wait for the writes to finish during flush period.
          wait_for_flush();
//          std::cout << "DEBUG: flushed buffer: " << *(this) << std::endl << std::flush;

//          std::atomic_thread_fence(std::memory_order_seq_cst);  // make sure all writes are done.
        }




      public:


        /**
         * @brief Clears the buffer. (set the size to 0, leaving the capacity and the memory allocation intact)
         * @note
         * const because the caller will have a const reference to the buffer.
         */
        void clear() {
          while (spinlock.test_and_set());
          writerCount = 0;
          // TODO: remove after finish debugging.
          memset(start_ptr, 0, capacity * sizeof(uint8_t));
          curr_ptr = MAX_PTR_ADDR;
          end_ptr = start_ptr;
          spinlock.clear();
        }

        /**
         * @brief get the capacity of the buffer.
         * @return    maximum capacity of buffer.
         */
        const uint32_t & getCapacity() const {
          return capacity;
        }

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
         * @brief Checks if a buffer is empty.
         * @note The return value is not precise - between getSize and return, other threads may have modified the size.
         * For "isEmpty" error is infrequent.
         *
         * @return    true if the buffer is empty, false otherwise.
         */
//        const bool isEmpty() const {
//          return (uint8_t*)curr_ptr <= start_ptr;
//        }


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


          uint8_t* ptr = write_lock<TS>(count);

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

            // don't need to write unlock since it was not locked.  but should read lock
            read_lock<TS>();  // not okay to call multiple times.

            return 0x0;  // full and not swapping.
          } else {  // ptr starts inside allocated buffer range.

            if ((ptr + count) > this->max_ptr) { // thread that filled the buffer

            	//printf(" JUST FULL.  end_ptr = %p, curr_ptr = %p. writerCount = %d.\n", end_ptr.load(), curr_ptr.load(), writerCount.load());

            	// read lock
              read_lock<TS>(ptr);
              write_unlock<TS>();

              // now wait for all other writes to complete
              //std::atomic_thread_fence(std::memory_order_seq_cst);  // only update size after all writes are done.


              if (ptr > this->max_ptr) {
                fprintf(stdout, "FAIL: ptr is %p, larger than maxptr %p\n", ptr, this->max_ptr);
                std::cout << " buffer: " << *this << std::endl << std::flush;
              }

              if (this->start_ptr == nullptr)
                std::cout << "FAIL null data ptr.  swap requested. buffer: " << *this << std::endl << std::flush;

              if (curr_ptr != MAX_PTR_ADDR) {
                fprintf(stdout, "FAIL before flush: flush bit should be set\n");
                std::cout << " buffer: " << *this << std::endl << std::flush;
              }

              // TODO: other threads may still be writing.  the problem here is that wait_for_flush is only waiting for writerCount to go to 0,
              // which may not be same as all pending threads finished writing.
              // counter point is that read_lock_write_unlock should have at least set flush bit to true.  not seeing that in testing.

              wait_for_flush();
              if (! is_reading()) {
                fprintf(stdout, "FAIL after flush: all writes should be done and flush set.\n");
                std::cout << " buffer: " << *this << std::endl << std::flush;
              }

              atomic_thread_fence(std::memory_order_seq_cst);   // commits all changes up to this point to memory


              if ((getSize() / sizeof(int)) != (capacity / sizeof(int))) {
                fprintf(stdout , "FAIL IN BUFFER:  NOT %lu elements. got %ld.", (capacity / sizeof(int)), getSize() / sizeof(int));
                std::cout << " buffer: " << *this << std::endl << std::flush;
              }




              return 0x2;  // was not full, now full.

            } else { // can insert.  may make buffer full
            	//printf(" COPY.  end_ptr = %p, curr_ptr = %p. writerCount = %d.\n", end_ptr.load(), curr_ptr.load(), writerCount.load());


              // write
              std::memcpy(ptr, _data, count);
              //std::atomic_thread_fence(std::memory_order_seq_cst);  // unlock only after memcpy.
              write_unlock<TS>();  // all full buffers lock the read and unlock the writer
              return 0x1;   // not full, successfully inserted.
            }

          }





//            atomic_thread_fence(std::memory_order_seq_cst);   // commits all changes up to this point to memory

           // free(ldata);
          //return std::make_pair(appended, swap);
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


          //DEBUG("Thread Unsafe buffer append");
          uint8_t* ptr = write_lock<TS>(count);


          if (ptr == nullptr) {

            read_lock<TS>();
            return 0x0;
          }

          if ((ptr + count) > max_ptr) {  // append that filled the buffer
            read_lock<TS>(ptr);
            write_unlock<TS>();  // also takes away an writer.

            wait_for_flush();

            //size = ptr - start_ptr; // final
            return 0x2;
          } else {


            std::memcpy(ptr, _data, count);
            write_unlock<TS>();

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
        << " approx,size/cap=" << buffer.getApproximateSize() << "," << buffer.getSize() << "/" << buffer.capacity
        << " writerCount=" << std::hex << (typename Buffer<ThreadSafety>::CountType)(buffer.writerCount) << std::dec
        << " R? " << (buffer.is_reading() ? "y" : "n") << " F? " << (buffer.is_flushing() ? "y" : "n") << " W? " << (buffer.is_writing() ? "y" : "n") << std::flush;


      return ost;
    }


    template<bliss::concurrent::ThreadSafety ThreadSafety>
    constexpr uint8_t* Buffer<ThreadSafety>::MAX_PTR_ADDR;

  } /* namespace io */
} /* namespace bliss */

#endif /* BUFFER_HPP_ */
