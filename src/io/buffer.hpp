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
        /// type of internal data pointer.
        typedef std::unique_ptr<uint8_t[], std::default_delete<uint8_t[]> >  DataType;


      protected:

        /// maximum capacity.  unrealistic to use size_t - can't possibly allocate.
        mutable uint32_t capacity;

        /// internal data storage
        mutable DataType data;
        mutable uint8_t* max_pointer;


        /**
         * @brief current occupied size of the buffer
         * @details depending on thread safety, either an atomic or a raw data type.
         *
         * @note mutable so a const Buffer object can still call append and modify the internal. (conceptually, size and data are constant members of the Buffer)
         */
        mutable typename std::conditional<ThreadSafety, std::atomic<int32_t>, int32_t>::type size;

        /**
         * @brief indicating whether the buffer is accepting data or not.
         * @details depending on thread safety, either an atomic or a raw data type.
         *          this is useful when swapping buffers in bufferPool or in MessageBuffers classes.
         *
         * @note mutable so a const Buffer object can still call append and modify the internal. (conceptually, size and data are constant members of the Buffer)
         */
//        volatile typename std::conditional<ThreadSafety, std::atomic<bool>, bool>::type flushing;   // no more future append.  can still have pending appends


        volatile typename std::conditional<ThreadSafety, std::atomic<uint8_t*>, uint8_t*>::type pointer;

        /// Reference Count of threads performing append update.  to ensure that all updates are done before used for sending.
        /// using signbit to indicate flushing/full
        volatile typename std::conditional<ThreadSafety, std::atomic<int>, int>::type writerCount;    // track how many are appending.


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
      	  : capacity(other.capacity), data(std::move(other.data)), max_pointer(other.max_pointer),
      	    size(int32_t(other.size)),
      	   // flushing(bool(other.flushing)),
      	    pointer((uint8_t*)(other.pointer)), writerCount(0) {

          //DEBUG("BUFFER private move contructor 1 called");

          // other.flushing = true;
          other.capacity = 0;
          other.size = -1;
          other.data = nullptr;
          other.pointer = nullptr;
          other.max_pointer = nullptr;
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
          : capacity(other.capacity), data(std::move(other.data)), max_pointer(other.max_pointer),
            size(int32_t(other.size)),
           // flushing(bool(other.flushing)),
            pointer((uint8_t*)(other.pointer)), writerCount(0) {

          //other.flushing = true;
          other.capacity = 0;
          other.size = -1;
          other.data = nullptr;
          other.pointer = nullptr;
          other.max_pointer = nullptr;
          other.writerCount = 0;
        }

      public:
        /**
         * @brief Normal constructor.  Allocate and initialize memory with size specified as parameter.
         * @param _capacity   The maximum capacity of the Buffer in bytes
         */
        explicit Buffer(const uint32_t _capacity) : capacity(_capacity), data(new uint8_t[_capacity]), size(-1), // flushing(false),
        writerCount(0)
        {
          if (_capacity == 0)
            throw std::invalid_argument("Buffer constructor parameter capacity is given as 0");

          memset(data.get(), 0, _capacity * sizeof(uint8_t));
          pointer = data.get();
          max_pointer = data.get() + _capacity;

          //printf("NEW BUFF:  start: %p, max_pointer %p, pointer %p, size %d, writers %d\n", data.get(), max_pointer, (uint8_t*)pointer, (int32_t)size, (int)writerCount);
        };

        /**
         * @brief reate a new Buffer using the memory specified along with allocated memory's size.
         * @param _data   The pointer/address of preallocated memory block
         * @param count   The size of preallocated memory block
         */
        Buffer(void* _data, const int32_t count) : capacity(count), data(static_cast<uint8_t*>(_data)), size(count), //flushing(false),
            writerCount(0) {
          if (count == 0)
            throw std::invalid_argument("Buffer constructor parameter count is given as 0");
          if (_data == nullptr)
            throw std::invalid_argument("Buffer constructor parameter _data is given as nullptr");

          pointer = data.get();
          max_pointer = data.get() + count;
        }

        /**
         * @brief Destructor.  deallocation of memory is automatic.
         */
        virtual ~Buffer() {};

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

          if (this->data != other.data) {

            std::unique_lock<std::mutex> myLock(mutex, std::defer_lock),
                                         otherLock(other.mutex, std::defer_lock);
            std::lock(myLock, otherLock);

            /// move the internal memory.
//            flushing = bool(other.flushing);
//            other.flushing = true;

            capacity = other.capacity;
            size = int32_t(other.size);
            data = std::move(other.data);
            pointer = (uint8_t*)(other.pointer);
            max_pointer = other.max_pointer;
            writerCount = int(other.writerCount);

            other.capacity = 0;
            other.size = -1;
            other.data = nullptr;
            other.pointer = nullptr;
            other.max_pointer = nullptr;
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

          if (this->data != other.data) {

            std::unique_lock<std::mutex> myLock(mutex, std::defer_lock),
                                         otherLock(other.mutex, std::defer_lock);
            std::lock(myLock, otherLock);

            /// move the internal memory.
   //         flushing = bool(other.flushing);             other.flushing = true;

            capacity = other.capacity;
            size = int32_t(other.size);
            data = std::move(other.data);
            pointer = (uint8_t*)(other.pointer);
            max_pointer = other.max_pointer;
            writerCount = int(other.writerCount);

            other.capacity = 0;
            other.size = -1;
            other.data = nullptr;
            other.pointer = nullptr;
            other.max_pointer = nullptr;
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
        template<bliss::concurrent::ThreadSafety TS = ThreadSafety>
        typename std::enable_if<TS, const int32_t>::type getSize() const {
          //std::atomic_thread_fence(std::memory_order_seq_cst);
          // wait until updating is done.
          //while(!isReading()) { _mm_pause(); };  // after all updates, no threads touching this.  so size would have been updated by a thread if flushing, else unset.
          return size.load(std::memory_order_seq_cst);

        }

        /**
         * @brief get the current size of the Buffer.
         * @tparam TS Choose thread safe vs unsafe implementation. defaults to same as the parent class.
         * @return    current size
         */
        template<bliss::concurrent::ThreadSafety TS = ThreadSafety>
        typename std::enable_if<!TS, const int32_t>::type getSize() const {
          //while(!isReading()) { _mm_pause(); };  // after all updates, no threads touching this.  so size would have been updated by a thread if flushing, else unset.
          return size;

        }

      protected:

        /**
         * @brief get the current size of the Buffer.
         * @tparam TS Choose thread safe vs unsafe implementation. defaults to same as the parent class.
         * @return    current size
         */
        template<bliss::concurrent::ThreadSafety TS = ThreadSafety>
        typename std::enable_if<TS, const int32_t>::type getApproximateSize() const {
          //std::atomic_thread_fence(std::memory_order_seq_cst);

          return pointer.load(std::memory_order_seq_cst) - data.get();
        }

        /**
         * @brief get the current size of the Buffer.
         * @tparam TS Choose thread safe vs unsafe implementation. defaults to same as the parent class.
         * @return    current size
         */
        template<bliss::concurrent::ThreadSafety TS = ThreadSafety>
        typename std::enable_if<!TS, const int32_t>::type getApproximateSize() const {
          return pointer - data.get();
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
          std::atomic_thread_fence(std::memory_order_seq_cst);

        	return (writerCount & std::numeric_limits<int>::max()) > 0;
        }

        inline bool is_flushing() const {
        	// flushing (MSB == 1) and some writes (lower bits > 0) == same as value in range (lowest(), 0), exclusive of both ends;
          std::atomic_thread_fence(std::memory_order_seq_cst);

          int x = writerCount;
        	return (x > std::numeric_limits<int>::lowest()) && (x < 0);
        	// reinterpret_cast does not generate any op codes.
        }

        inline bool is_reading() const {
        	// flushing (MSB == 1) and no write.
          std::atomic_thread_fence(std::memory_order_seq_cst);
        	return writerCount == std::numeric_limits<int>::lowest();
        }

        /// increment lock only if flush bit is not set.  else do nothing
    	// increment in unsigned and 2's complement follow the same directionality when signbit is excluded.
    	// so this means that sign bit can change by another thread and the lower bits remain valid, even is reordered between threads.
    	// e.g. 100 + 1 = 101. 101 & 011 = 001.  001 + 1 = 010. 010 - 1 = 001, 001 - 1 = 000
    	// so we can atomically increment, then compare to lowest() (100), and undo or not.
    	// we rely on 2 things: 1. no 0x7FFFFFFF threads in single node, so no chance of signbit being changed from overflow.
        // 2. i and i+1 have the same lower bit patterns whether i is negative or positive.
        // return bool indicating success or failure. since we increment only if not flushing, we need to known success/fail
        // to avoid extra unlocks.
       // no chance of overflow with addition unless we have machines with 2B threads.  goal is to preserve sign bit during this op.
        inline bool write_lock() {
        	// atomic fetch_add as post increment returns old WriterCount value.
        	// check if old value indicates MSB is set (< 0).  if so, undo.
        	if ( 0 > writerCount++ ) {
        		writerCount--;
            std::atomic_thread_fence(std::memory_order_seq_cst);  // make sure all writes are done.

        		return false;
        	} else {
            std::atomic_thread_fence(std::memory_order_seq_cst);  // make sure all writes are done.

        		return true;
        	}
        }


        /// cannot decrement past 0 ( with or without sign bit).  should be called with knowledge of whether try_lock_write succeeded or not.
        /// throws exception if there are no writers
        // goal is to preserve sign bit for the entire duration of the op.
        template<bliss::concurrent::ThreadSafety TS = ThreadSafety>
        inline typename std::enable_if<TS, void>::type write_unlock() {
          int lcount = writerCount;

          if ((lcount & std::numeric_limits<int>::max()) == 0)
            throw std::logic_error("ERROR: write_unlock 1 on 0 writers");

          std::atomic_thread_fence(std::memory_order_seq_cst);  // make sure all writes are done.

          // now update
          while (!writerCount.compare_exchange_weak(lcount, lcount - 1)) {
            if ((lcount & std::numeric_limits<int>::max()) == 0)
              throw std::logic_error("ERROR: write_unlock 2 on 0 writers");
            std::atomic_thread_fence(std::memory_order_seq_cst);  // make sure all writes are done.
          }
        }

        template<bliss::concurrent::ThreadSafety TS = ThreadSafety>
        inline typename std::enable_if<!TS, void>::type write_unlock() {
          if ((writerCount & std::numeric_limits<int>::max()) == 0)
            throw std::logic_error("ERROR: write_unlock thread unsafe on 0 writers");
          --writerCount;
        }





        //===== read lock trumps write.
        /// read lock prevents future writes.  read lock can be turned on while there are writes in progress.  once read lock is on, lock_write will fail.
        // thread safe - multiplr threads calling is the same as 1 thread calling
        inline void read_lock() {  // not using xor since it toggles and write_unlock is not going to change the sign bit (so we don't need xor)
          writerCount |= std::numeric_limits<int>::lowest();
          std::atomic_thread_fence(std::memory_order_seq_cst);  // make sure all writes are done.
        }

        /// read unlock allows future writes.  read unlock blocks until all pending writes are done before .
        inline void read_unlock() {  // not using xor.
          writerCount &= std::numeric_limits<int>::max();
          std::atomic_thread_fence(std::memory_order_seq_cst);  // make sure all writes are done.

        }

        void wait_for_writes() {
          while (is_writing()) _mm_pause();
        }

        // has to ensure flush bit is set else return
        void wait_for_flush() {
          int lcount = writerCount;
          if ((lcount & std::numeric_limits<int>::lowest()) == 0)  // check flush bit.
            throw std::logic_error("ERROR: wait for flush called but flush bit is not set");
          std::atomic_thread_fence(std::memory_order_seq_cst);  // make sure all writes are done.

          while((lcount & std::numeric_limits<int>::max()) > 0 ) {
            _mm_pause();  // all threads attempting to lock read are waiting until writing is done.
            lcount = writerCount;
            if ((lcount & std::numeric_limits<int>::lowest()) == 0)
              throw std::logic_error("ERROR: wait for flush called but flush bit is not set");
            std::atomic_thread_fence(std::memory_order_seq_cst);  // make sure all writes are done.

          }
        }


        // ONLY CALL IF NOT FULL.
        template<bliss::concurrent::ThreadSafety TS = ThreadSafety>
        inline typename std::enable_if<TS, void>::type flush_and_set_size() {
          int old = writerCount.fetch_or(std::numeric_limits<int>::lowest());  // set the flush bit and get the old value
          std::atomic_thread_fence(std::memory_order_seq_cst);  // make sure all writes are done.

          // now everyone wait for the writes to finish during flush period.
          wait_for_flush();
          std::atomic_thread_fence(std::memory_order_seq_cst);  // make sure all writes are done.


          if ((old & std::numeric_limits<int>::lowest()) == 0) {// old  value did  not have flush bit set
            // so now set the size to the approximate size;

            // now update the size.  all writes are done, so if buffer became full during wait_for_flush, then size will be set.
            // so we can use compare exchange.  (this cannot happen before the size set when full.)
            int temp = -1;  // unset size is -1
            size.compare_exchange_strong(temp, getApproximateSize<TS>());   // only change if size is not yet set.  what if buffer size is supposed to be 0?
            std::atomic_thread_fence(std::memory_order_seq_cst);  // make sure all writes are done.

          }
        }
        template<bliss::concurrent::ThreadSafety TS = ThreadSafety>
        inline typename std::enable_if<!TS, void>::type flush_and_set_size() {
          int old = writerCount;
          writerCount |= std::numeric_limits<int>::lowest();  // set flush
          // wait for all writes to finish.
          wait_for_flush();

          if ((old & std::numeric_limits<int>::lowest()) == 0) {// old value did  not have flush bit set
            // if flushing was false before,  then update size.  else size was already updated previously.
            // now update the size
            if (size == -1)  // only update if size was not set.
              size = getApproximateSize<TS>();
          }
        }


      protected:

        template<bliss::concurrent::ThreadSafety TS = ThreadSafety>
        inline typename std::enable_if<TS, void>::type read_lock_write_unlock() {
          int lcount = writerCount;

          if ((lcount & std::numeric_limits<int>::max()) == 0)
            throw std::logic_error("ERROR: read_lock_write_unlock 1 on 0 writers");

          std::atomic_thread_fence(std::memory_order_seq_cst);  // make sure all writes are done.

          // now update
          while (!writerCount.compare_exchange_weak(lcount, ( (lcount - 1) | std::numeric_limits<int>::lowest() ) ) ) {
            if ((lcount & std::numeric_limits<int>::max()) == 0)
              throw std::logic_error("ERROR: read_lock_write_unlock 2 on 0 writers");
            std::atomic_thread_fence(std::memory_order_seq_cst);  // make sure all writes are done.
          }
        }
        template<bliss::concurrent::ThreadSafety TS = ThreadSafety>
        inline typename std::enable_if<!TS, void>::type read_lock_write_unlock() {
          if ((writerCount & std::numeric_limits<int>::max()) == 0)
            throw std::logic_error("ERROR: read_lock_write_unlock thread unsafe on 0 writers");
          --writerCount;
          writerCount |= std::numeric_limits<int>::lowest();
        }


      public:


        /**
         * @brief Clears the buffer. (set the size to 0, leaving the capacity and the memory allocation intact)
         * @note
         * const because the caller will have a const reference to the buffer.
         */
        void clear() {
          writerCount = 0;
          size = -1;
          // TODO: remove after finish debugging.
          memset(data.get(), 0, capacity * sizeof(uint8_t));
          pointer = data.get();
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
          return reinterpret_cast<T*>(data.get());
        }


        /**
         * @brief Checks if a buffer is empty.
         * @note The return value is not precise - between getSize and return, other threads may have modified the size.
         * For "isEmpty" error is infrequent.
         *
         * @return    true if the buffer is empty, false otherwise.
         */
        const bool isEmpty() const {
          return (uint8_t*)pointer <= data.get() || (size == 0 || size == -1);
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

//          volatile bool swap = false;
//          volatile bool appended = false;


          // if fails, then already flushing
          if (!write_lock()) {  // lock checks if buffer is being read/flushing, if no, then increment writer count.
            //printf("flushing     shared pointer %p, count %d, start %p, max_pointer %p, size %d, writers %d\n", (uint8_t*)pointer, count, data.get(), max_pointer, (int32_t)size, (int)writerCount);
            return 0x0;
          }


            // After the F&A, there are 3 types of threads:  target write area is
            //  1. completely within buffer:  these threads should proceed to memcpy with the local ptr var
            //  2. completely outside buffer: these threads will not memcpy.  they all disabled buffer, and
            //     only 1 thread from 2) or 3) should swap in a new buffer atomically.
            //      - if curr buffer is not disabled (already swapped), then don't swap further.
            //  3. crossing buffer boundary: this thread will not memcpy.  disabled.  and it should retract the pointer advance.


            // ONE thread at a time can get the position to copy in data.  If flushing, then the buffer is flushing.
            std::atomic_thread_fence(std::memory_order_seq_cst);    // compute ptr strictly after determining can lock write.
            uint8_t* ptr = pointer.fetch_add(count, std::memory_order_seq_cst);  // no memory ordering needed within mutex loc
            std::atomic_thread_fence(std::memory_order_seq_cst);    // branch and memcpy only after computing the ptr.

            if ((this->max_pointer - ptr) >= static_cast<std::ptrdiff_t>(count) ) {
              // has room to copy.


              std::memcpy(ptr, _data, count);

              std::atomic_thread_fence(std::memory_order_seq_cst);  // unlock only after memcpy.
              //appended = true;
              write_unlock<TS>();

              atomic_thread_fence(std::memory_order_seq_cst);   // commits all changes up to this point to memory

              return 0x1;
            } else {
              read_lock_write_unlock<TS>();  // all full buffers lock the read and unlock the writer
              std::atomic_thread_fence(std::memory_order_seq_cst);  // ensure these are completed.
//              if ((writerCount & std::numeric_limits<int>::lowest()) == 0) {
//                fprintf(stdout, "FAIL: flush bit should be set\n");
//              }


              if (ptr <= this->max_pointer) {  // thread that made the buffer flushing.
                // now wait for all other writes to complete
                std::atomic_thread_fence(std::memory_order_seq_cst);  // only update size after all writes are done.

                if ((writerCount & std::numeric_limits<int>::lowest()) == 0) {
                  fprintf(stdout, "FAIL before flush: flush bit should be set\n");
                }
                wait_for_flush();
                if ((writerCount & std::numeric_limits<int>::lowest()) == 0) {
                  fprintf(stdout, "FAIL after flush: flush bit should be set\n");
                }

//                std::atomic_thread_fence(std::memory_order_seq_cst);  // only update size after all writes are done.

                this->size.store(ptr - data.get(), std::memory_order_seq_cst);  // overwrite any previous size with this one's
                if ((writerCount & std::numeric_limits<int>::lowest()) == 0) {
                  fprintf(stdout, "FAIL after setsize: flush bit should be set\n");
                }

//                swap = true;
                atomic_thread_fence(std::memory_order_seq_cst);   // commits all changes up to this point to memory


                if ((size.load() / sizeof(int)) != (capacity / sizeof(int))) {
                  fprintf(stdout , "IN BUFFER:  NOT 2047 elements. got %ld.", size.load() / sizeof(int));
                  std::cout << " buffer: " << *this << std::endl << std::flush;
                }
                if ((writerCount & std::numeric_limits<int>::lowest()) == 0) {
                  fprintf(stdout, "FAIL before return: flush bit should be set\n");
                }

                return 0x2;
              } else {
                return 0x0;
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
          if (!write_lock()) return 0x0;

          if ((pointer + count) > max_pointer) {
            read_lock_write_unlock<TS>();  // also takes away an writer.

            if (pointer <= max_pointer)  {// original pointer is within buffer - put it back. however, this could cause multiple threads to reach this point.
              wait_for_flush();

              size = pointer - data.get(); // final
              return 0x2;
            } else
              return 0x0;
          } else {


            std::memcpy(pointer, _data, count);
            pointer += count;
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
      ost << "THREAD "<< (ThreadSafety ? "SAFE" : "UNSAFE") << " BUFFER: data_ptr=" << static_cast <const void *>(buffer.data.get())
        << " ptr/maxptr=" << static_cast <const void *>((uint8_t*)(buffer.pointer)) << "/" << static_cast <const void *>(buffer.max_pointer)
        << " approx,size/cap=" << buffer.getApproximateSize() << "," << int(buffer.size) << "/" << buffer.capacity
        << " writerCount=" << std::hex << int(buffer.writerCount) << std::dec << std::endl << std::flush;


      return ost;
    }




  } /* namespace io */
} /* namespace bliss */

#endif /* BUFFER_HPP_ */
