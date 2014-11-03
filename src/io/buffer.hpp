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
#include <unistd.h>  // usleep

#include <atomic>
#include <mutex>

#include <memory>     // unique_ptr
#include <utility>    // move, forward, swap, make_pair
#include <stdexcept>

#include "concurrent/concurrent.hpp"   // ThreadSafety boolean constants
#include "utils/logging.h"

#include "xmmintrin.h"

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
        typename std::conditional<ThreadSafety, std::atomic<int32_t>, int32_t>::type size;

        /**
         * @brief indicating whether the buffer is accepting data or not.
         * @details depending on thread safety, either an atomic or a raw data type.
         *          this is useful when swapping buffers in bufferPool or in MessageBuffers classes.
         *
         * @note mutable so a const Buffer object can still call append and modify the internal. (conceptually, size and data are constant members of the Buffer)
         */
        typename std::conditional<ThreadSafety, std::atomic<bool>, bool>::type full;   // no more future append.  can still have pending appends


        typename std::conditional<ThreadSafety, std::atomic<uint8_t*>, uint8_t*>::type pointer;

        /// Reference Count of threads performing append update.  to ensure that all updates are done before used for sending.
        typename std::conditional<ThreadSafety, std::atomic<int>, int>::type writerCount;    // track how many are appending.  writerCount == int::MAX means reading thread is active.


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
      	    //full(bool(other.full)),
      	    pointer((uint8_t*)(other.pointer)), writerCount(0) {

          //DEBUG("BUFFER private move contructor 1 called");

          //other.full = true;
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
            //full(bool(other.full)),
            pointer((uint8_t*)(other.pointer)), writerCount(0) {

          DEBUG("BUFFER private move contructor 2 called");


          //other.full = true;
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
        explicit Buffer(const uint32_t _capacity) : capacity(_capacity), data(new uint8_t[_capacity]), size(-1), //full(false),
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
        Buffer(void* _data, const int32_t count) : capacity(count), data(static_cast<uint8_t*>(_data)), size(count), //full(false),
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
            //full = bool(other.full);
            //other.full = true;

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
            //full = bool(other.full);             other.full = true;

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
        typename std::enable_if<TS, const int32_t>::type getFinalSize() const {
          //std::atomic_thread_fence(std::memory_order_seq_cst);
          // wait until updating is done.
          //while(!isReading()) { _mm_pause(); };  // after all updates, no threads touching this.  so size would have been updated by a thread if full, else unset.
          return size.load(std::memory_order_seq_cst);

        }

        /**
         * @brief get the current size of the Buffer.
         * @tparam TS Choose thread safe vs unsafe implementation. defaults to same as the parent class.
         * @return    current size
         */
        template<bliss::concurrent::ThreadSafety TS = ThreadSafety>
        typename std::enable_if<!TS, const int32_t>::type getFinalSize() const {
          //while(!isReading()) { _mm_pause(); };  // after all updates, no threads touching this.  so size would have been updated by a thread if full, else unset.
          return size;

        }

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


        /**
         * @brief get the status of whether buffer is accepting new appends
         * @tparam TS Choose thread safe vs unsafe implementation. defaults to same as the parent class.
         * @return    true if buffer is NOT accepting more data
         */
//        template<bliss::concurrent::ThreadSafety TS = ThreadSafety>
//        typename std::enable_if<TS, const bool>::type isFull() const {
//          return full.load(std::memory_order_seq_cst);
//        }
        /**
         * @brief get the status of whether buffer is accepting new appends
         * @tparam TS Choose thread safe vs unsafe implementation. defaults to same as the parent class.
         * @return    true if buffer is NOT accepting more data
         */
//        template<bliss::concurrent::ThreadSafety TS = ThreadSafety>
//        typename std::enable_if<!TS, const bool>::type isFull() const {
//          return full;
//        }
        /**
         * @brief block the buffer from accepting more data.
         * @tparam TS Choose thread safe vs unsafe implementation. defaults to same as the parent class.
         */
//        template<bliss::concurrent::ThreadSafety TS = ThreadSafety>
//        typename std::enable_if<TS, void>::type markFull() {
//          bool init = false;
//          if (full.compare_exchange_strong(init, true, std::memory_order_seq_cst)) {
//            //writerCount.fetch_sub(1, std::memory_order_seq_cst);   // transition from true to false, so add one writer
//            size.store(getApproximateSize<TS>());  // one thread writes here.  if blocking because full, size will be overwritten by the correct one.
//          } // else already false.
//        }
        /**
         * @brief markFull the buffer from accepting more data.
         * @tparam TS Choose thread safe vs unsafe implementation. defaults to same as the parent class.
         */
//        template<bliss::concurrent::ThreadSafety TS = ThreadSafety>
//        typename std::enable_if<!TS, void>::type markFull() {
//          if (!full) { // transition from true to false, so remove one writer
//            //--writerCount;
//            size = getApproximateSize<TS>();  // one thread writes here.  if blocking because full, size will be overwritten by the correct one.
//          }
//          full = true;
//        }

        /**
         * @brief unmarkFull the buffer to accept more data.
         * @tparam TS Choose thread safe vs unsafe implementation. defaults to same as the parent class.
         */
//        template<bliss::concurrent::ThreadSafety TS = ThreadSafety>
//        typename std::enable_if<TS, void>::type unmarkFull() {
//          bool init = true;
//          if (full.compare_exchange_strong(init, false, std::memory_order_seq_cst)) {
//            //writerCount.fetch_add(1, std::memory_order_seq_cst);   // transition from true to false, so add one writer
//          } // else already false.
//        }
        /**
         * @brief unmarkFull the buffer to accept more data.
         * @tparam TS Choose thread safe vs unsafe implementation. defaults to same as the parent class.
         */
//        template<bliss::concurrent::ThreadSafety TS = ThreadSafety>
//        typename std::enable_if<!TS, void>::type unmarkFull() {
//          if (full) { // transition from true to false, so add one writer
//            //++writerCount;
//          }
//          full = false;
//        }

        inline void lock_write() {
          ++writerCount;  // normal ++.
        }
        inline void unlock_write() {
          --writerCount;  // normal --.
        }

        inline bool isWriting() const {
          return writerCount > 0;
        }
        inline bool isQuiet() const {
          return writerCount == 0;
        }

        inline void lock_read() {  // only one thread can call this.
          full = true;
          while(isWriting()) _mm_pause();  // all threads attempting to lock read are waiting until writing is done.
//          std::atomic_thread_fence(std::memory_order_seq_cst);  // ensure approximate size is updated correctly
//          size = getApproximateSize<ThreadSafety>();
//          std::atomic_thread_fence(std::memory_order_seq_cst);  // ensure subsequent call to change size (by the
        }

        inline void unlock_read() {
          full = false;
          assert(writerCount == 0);
        }

        inline bool isReading() const {
          return bool(full);
        }


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
         * Semantically, a "false" return value means there is not enough room for the new data, not that the buffer is full.
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
         * @return            bool indicated whether operation succeeded, bool indicating whether buffer swap is needed.
         */
        template<bliss::concurrent::ThreadSafety TS = ThreadSafety>
        typename std::enable_if<TS, std::pair<bool, bool> >::type append(const void* _data, const uint32_t count) {


//


          if (count == 0) {
              //printf("append EMPTY...\n");
          	  return std::make_pair(false, false);
          }
          if (_data == nullptr)
        	  throw std::invalid_argument("_data is nullptr");

          assert(capacity != 1);
          if (capacity == 0)
            throw std::length_error("ThreadSafe Buffer capacity is 0");


          // local copy of data
//          void* ldata = malloc(count);
//          std::unique_lock<std::mutex> lock(mutex);   // block() and size need to be synchronized consistently across threads.
//          memcpy(ldata, _data, count);
//          lock.unlock();


          bool swap = false;
          bool appended = false;


          // ONE thread at a time can get the position to copy in data.  If full, then the buffer is full.


//          std::unique_lock<std::mutex> lock(mutex);
//
//          //DEBUG("Thread Safe buffer append");
//          if (isFull<TS>()) {
//            //printf("full     shared pointer %p, count %d, start %p, max_pointer %p, size %d, writers %d\n", (uint8_t*)pointer, count, data.get(), max_pointer, (int32_t)size, (int)writerCount);
//            return std::make_pair(false, false);
//          }
//
//            writerCount.fetch_add(1, std::memory_order_seq_cst);
//          uint8_t* ptr = pointer.fetch_add(count, std::memory_order_seq_cst);  // no memory ordering needed within mutex lock
//          if ((ptr + count) > max_pointer) { // FULL   // max_pointer may be swapped during buffer move.  TEST: disable that and use pointer to buffers only
//
//           lock.unlock();
//
//
//            // After the F&A, there are 3 types of threads:  target write area is
//        	  //  1. completely within buffer:	these threads should proceed to memcpy with the local ptr var
//        	  //  2. completely outside buffer:	these threads will not memcpy.  they all disabled buffer, and
//        	  //     only 1 thread from 2) or 3) should swap in a new buffer atomically.
//        	  //     	- if curr buffer is not disabled (already swapped), then don't swap further.
//        	  //  3. crossing buffer boundary: this thread will not memcpy.  disabled.  and it should retract the pointer advance.
//          	if (ptr <= max_pointer)  {// original pointer is within buffer - put it back. however, this could cause multiple threads to reach this point.
//          	  size.store(ptr - start, std::memory_order_seq_cst);  // overwrite the size with this one's
//          	  swap = true;
//
//            block<TS>();   // will also take away an writer.
//              //printf("FULL    SWAP ptr %p, shared pointer %p, count %d, start %p, max_pointer %p, size %d, ptrdiff %ld, swap %s, writers %d\n", ptr, (uint8_t*)pointer, count, data.get(), max_pointer, (uint32_t)size, (ptr - data.get()), (swap ? "Y" : "N"), (int)writerCount);
//          	} // else swap is false;
//          	//else {
//              //printf("FULL NO SWAP ptr %p, shared pointer %p, count %d, start %p, max_pointer %p, size %d, ptrdiff %ld, swap %s, writers %d\n", ptr, (uint8_t*)pointer, count, data.get(), max_pointer, (uint32_t)size, (ptr - data.get()), (swap ? "Y" : "N"), (int)writerCount);
//          	// }
//
//            	// at this point, buffer is full, and can be replaced.  this may be faster than we can copy the data.
//            				// to fix this - compute the pointer where data is to be copied before unlock.  For threads that
//            				// pass the capacity test, they continue to memcpy, knowing where the addresses are.  The failed ones will
//            				// cause a buffer swap, but the passed threads has the absolute address.
//            writerCount.fetch_sub(1, std::memory_order_seq_cst);
//          } else {
//            lock.unlock();
//
//            appended = true;
//            std::memcpy(ptr, _data, count);
//            writerCount.fetch_sub(1, std::memory_order_seq_cst);
//          }


// same number of missing and extra elements.
//            std::unique_lock<std::mutex> lock(mutex);

            //DEBUG("Thread Safe buffer append");
            if (isReading()) {  // already full and being read.
            //  free(ldata);
              //printf("full     shared pointer %p, count %d, start %p, max_pointer %p, size %d, writers %d\n", (uint8_t*)pointer, count, data.get(), max_pointer, (int32_t)size, (int)writerCount);
              return std::make_pair(false, false);
            }


            lock_write();
            std::atomic_thread_fence(std::memory_order_seq_cst);  // pointer++ only after lock_write.

             // if fails, then already full


            // After the F&A, there are 3 types of threads:  target write area is
            //  1. completely within buffer:  these threads should proceed to memcpy with the local ptr var
            //  2. completely outside buffer: these threads will not memcpy.  they all disabled buffer, and
            //     only 1 thread from 2) or 3) should swap in a new buffer atomically.
            //      - if curr buffer is not disabled (already swapped), then don't swap further.
            //  3. crossing buffer boundary: this thread will not memcpy.  disabled.  and it should retract the pointer advance.

            uint8_t* ptr = pointer.fetch_add(count, std::memory_order_seq_cst);  // no memory ordering needed within mutex loc
            std::atomic_thread_fence(std::memory_order_seq_cst);    // branch and memcpy only after computing the ptr.

            if ((this->max_pointer - ptr) >= static_cast<std::ptrdiff_t>(count) ) {
              // has room to copy.


              std::memcpy(ptr, _data, count);

              std::atomic_thread_fence(std::memory_order_seq_cst);  // unlock only after memcpy.
              unlock_write();
              appended = true;


            } else {

              unlock_write();
              std::atomic_thread_fence(std::memory_order_seq_cst);  // lock_read needs to happen after unlock_write.
              lock_read();

              if (ptr <= this->max_pointer) {
                // lock it for read.


                this->size.store(ptr - data.get(), std::memory_order_seq_cst);  // overwrite the size with this one's

                swap = true;
              }


            }

            atomic_thread_fence(std::memory_order_seq_cst);   // commits all changes up to this point to memory

           // free(ldata);
          return std::make_pair(appended, swap);
        }

        /**
         * @brief Append data to the buffer, THREAD UNSAFE.
         * @details  The function updates the current occupied size of the Buffer
         * and memcopies the supplied data into the internal memory block.
         *
         *  This is the THREAD UNSAFE version.
         *
         * Semantically, a "false" return value means there is not enough room for the new data, not that the buffer is full.
         *
         * method is const because the caller will have const reference to the Buffer.
         *
         * @tparam TS Choose thread safe vs unsafe implementation. defaults to same as the parent class.
         * @param[in] _data   pointer to data to be copied into the Buffer
         * @param[in] count   number of bytes to be copied into the Buffer
         * @return            bool indicated whether operation succeeded.
         */
        template<bliss::concurrent::ThreadSafety TS = ThreadSafety>
        typename std::enable_if<!TS, std::pair<bool, bool> >::type append(const void* _data, const uint32_t count) {

          //DEBUG("Thread Unsafe buffer append");
          if (isReading()) return std::make_pair(false, false);


            if (count == 0)
            	  return std::make_pair(true, false);
            if (_data == nullptr)
          	  throw std::invalid_argument("_data is nullptr");


          assert(capacity != 1);
          if (capacity == 0)
            throw std::length_error("Thread Unsafe Buffer capacity is 0");

          lock_write();
          bool swap = false;
          if ((pointer + count) > max_pointer) {

            unlock_write();
            lock_read();  // also takes away an writer.

            if (pointer <= max_pointer)  {// original pointer is within buffer - put it back. however, this could cause multiple threads to reach this point.

              size = pointer - data.get(); // final
              swap = true;
            }
            atomic_thread_fence(std::memory_order_seq_cst);   // commits all changes up to this point to memory

            return std::make_pair(false, swap);
          }


          std::memcpy(pointer, _data, count);
          unlock_write();

          pointer += count;
          --writerCount;
          atomic_thread_fence(std::memory_order_seq_cst);   // commits all changes up to this point to memory

          return std::make_pair(true, swap);
        }

//        /**
//         * @brief Append data to the buffer, THREAD SAFE.
//         * @details  The function updates the current occupied size of the Buffer
//         * and memcopies the supplied data into the internal memory block.
//         *
//         * This is the THREAD SAFE version using Compare And Swap operation in a loop, lookfree.
//         *  sync version is faster on 4 threads.  using compare_exchange_strong (_weak causes extra loops but is supposed to be faster)
//         *    using relaxed and release vs acquire and acq_rel - no major difference.
//         *
//         * Semantically, a "false" return value means there is not enough room for the new data, not that the buffer is full.
//         *
//         * method is const because the caller will have const reference to the Buffer.
//         *
//         * NOTE: we can't use memory ordering alone.  WE have to lock with a mutex or use CAS in a loop
//         * See Thread Safe "append" method documentation for rationale.
//         *
//         * @tparam TS Choose thread safe vs unsafe implementation. defaults to same as the parent class.
//         * @param[in] _data   pointer to data to be copied into the Buffer
//         * @param[in] count   number of bytes to be copied into the Buffer
//         * @return            bool indicated whether operation succeeded.
//         */
//        template<bliss::concurrent::ThreadSafety TS = ThreadSafety>
//        typename std::enable_if<TS, const bool>::type append_lockfree(const void* typed_data, const uint32_t count) const {
//
//            if (count == 0)
//            	  return true;
//            if (_data == nullptr)
//          	  throw std::invalid_argument("_data is nullptr");
//
//
//          assert(capacity != 1);
//          if (capacity == 0)
//            throw std::length_error("ThreadSafe LockFree Buffer capacity is 0");
//
//          // since this is not using a lock, each iteration needs to check when acquiring a position whether the
//          // buffer is full or not.  If full, then new calls here would not go through
//          if (isFull<TS>()) return false;
//
//          // if the new position does not exceed capacity, then we can let it through. (data only)
//
//
//          // try with compare and exchange weak.
//          uint32_t s = size.load(std::memory_order_seq_cst);
//          uint32_t ns = s + count;
//          if (ns > capacity) {
//            markFull<TS>();
//            return false;
//          }
//          while (!size.compare_exchange_strong(s, ns, std::memory_order_seq_cst, std::memory_order_seq_cst)) {
//            // choosing strong, slower but less spurious "false" return, and should allow the "return" statement to more reliably be reached.
//            ns = s + count;
//            if (ns > capacity) {
//              markFull<TS>();
//              return false;
//            }
//          }
//
//          std::memcpy(data.get() + s, typed_data, count);
//          return true;
//        }

    };

  } /* namespace io */
} /* namespace bliss */

#endif /* BUFFER_HPP_ */
