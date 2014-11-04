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
        mutable typename std::conditional<ThreadSafety, std::atomic<int32_t>, int32_t>::type size;

        /**
         * @brief indicating whether the buffer is accepting data or not.
         * @details depending on thread safety, either an atomic or a raw data type.
         *          this is useful when swapping buffers in bufferPool or in MessageBuffers classes.
         *
         * @note mutable so a const Buffer object can still call append and modify the internal. (conceptually, size and data are constant members of the Buffer)
         */
        volatile typename std::conditional<ThreadSafety, std::atomic<bool>, bool>::type full;   // no more future append.  can still have pending appends


        volatile typename std::conditional<ThreadSafety, std::atomic<uint8_t*>, uint8_t*>::type pointer;

        /// Reference Count of threads performing append update.  to ensure that all updates are done before used for sending.
        volatile typename std::conditional<ThreadSafety, std::atomic<int>, int>::type writerCount;    // track how many are appending.  writerCount == int::MAX means reading thread is active.


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
      	    full(bool(other.full)),
      	    pointer((uint8_t*)(other.pointer)), writerCount(0) {

          //DEBUG("BUFFER private move contructor 1 called");

          other.full = true;
          other.capacity = 0;
          other.size = 0;
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
            full(bool(other.full)),
            pointer((uint8_t*)(other.pointer)), writerCount(0) {

          other.full = true;
          other.capacity = 0;
          other.size = 0;
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
        explicit Buffer(const uint32_t _capacity) : capacity(_capacity), data(new uint8_t[_capacity]), size(0), full(false),
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
        Buffer(void* _data, const int32_t count) : capacity(count), data(static_cast<uint8_t*>(_data)), size(count), full(false),
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
            full = bool(other.full);
            other.full = true;

            capacity = other.capacity;
            size = int32_t(other.size);
            data = std::move(other.data);
            pointer = (uint8_t*)(other.pointer);
            max_pointer = other.max_pointer;
            writerCount = int(other.writerCount);

            other.capacity = 0;
            other.size = 0;
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
            full = bool(other.full);             other.full = true;

            capacity = other.capacity;
            size = int32_t(other.size);
            data = std::move(other.data);
            pointer = (uint8_t*)(other.pointer);
            max_pointer = other.max_pointer;
            writerCount = int(other.writerCount);

            other.capacity = 0;
            other.size = 0;
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


        // TODO:  make this more clear.
        // need 1. flag to indicate no more new writers
        //      2. wait for current writers to finish
        //      3. 3 phases:  writing, flushing, reading.  -> write_lock, wirte_unlock, flush_begin, read_lock, read_unlock
        //          write_lock (counter ++)
        //          wirte_unlock (counter --)
        //          flush_begin (writable = false)
        //          read_lock ( wait for writable = false && counter == 0)
        //          read_unlock ( writable = true)
        // state checks:  is_writing (writable == true && counter > 0)
        //                is_flushing (writable == false && counter> 0)
        //                is_reading ( writable == false && counter == 0)
        //   if unsigned type, can use highest bit to manage this.



        inline bool isReading() const {
          return bool(full);
        }

        //===== read lock trumps write.
        /// read lock prevents future writes.  read lock can be turned on while there are writes in progress.  once read lock is on, lock_write will fail.
        inline void lock_read2() {  // only one thread can call this.
          full = true;
        }

        /// read unlock allows future writes.  read unlock blocks until all pending writes are done before .
        inline void unlock_read() {
          full = false;
        }

        template<bliss::concurrent::ThreadSafety TS = ThreadSafety>
        inline typename std::enable_if<TS, void>::type force_lock_read() {
          bool temp = false;
          // wait for all writes to finish.
          wait_for_all_writes();

          std::atomic_thread_fence(std::memory_order_seq_cst);  // make sure all writes are done.
          if (full.compare_exchange_strong(temp, true)) {
            // if full was set to false before,  then update size.  else size was already updated previously.
            // now update the size
            size = getApproximateSize<TS>();
          }
        }

        template<bliss::concurrent::ThreadSafety TS = ThreadSafety>
        inline typename std::enable_if<!TS, void>::type force_lock_read() {
          // wait for all writes to finish.
          wait_for_all_writes();

          if (full == false) {
            full = true;
            // if full was set to false before,  then update size.  else size was already updated previously.
            // now update the size
            size = getApproximateSize<TS>();
          }
        }




        //------ protected?
        inline bool try_lock_write() {
          // no DCAS type of method so can't combine these 2 lines.
          if (isReading()) return false;
          ++writerCount;  // normal ++.
          return true;
        }
        //------ protected?
        /// cannot decrement past 0.  should be called with knowledge of whether try_lock_write succedded or not.
        inline void unlock_write() {
          // uses post increment (== writerCount.fetch_sub(1))
          if (writerCount-- <= 0) {  // atomic in decrement, post decrement so old value is returned, so compare to 0
            writerCount++;  // error, restore old value.
            throw std::logic_error("unmatched unlock_write present.  Please fix.");
          }

        }


        void wait_for_all_writes() {
          if (!isReading()) {  //require read lock to provide guarantee that we don't see new lock_writes occur after this method starts.
            throw std::logic_error("wait_for_all_writes should be called only after lock_read");
          }
          while(isWriting<ThreadSafety>()) {
            _mm_pause();  // all threads attempting to lock read are waiting until writing is done.
          }
//          std::atomic_thread_fence(std::memory_order_seq_cst);  // ensure approximate size is updated correctly
//          size = getApproximateSize<ThreadSafety>();
//          std::atomic_thread_fence(std::memory_order_seq_cst);  // ensure subsequent call to change size (by the
        }

        template<bliss::concurrent::ThreadSafety TS = ThreadSafety>
        inline typename std::enable_if<TS, bool>::type isWriting() const {
          return writerCount.load() > 0;
        }
        template<bliss::concurrent::ThreadSafety TS = ThreadSafety>
        inline typename std::enable_if<!TS, bool>::type isWriting() const {
          return writerCount > 0;
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
         * @return            unsigned char, bit 0 indicated whether operation succeeded, bit 1 indicating whether buffer swap is needed.  single primitive type is better as it can be returned in 1 atomic op.
         */
        template<bliss::concurrent::ThreadSafety TS = ThreadSafety>
        typename std::enable_if<TS, unsigned int>::type append(const void* _data, const uint32_t count) {

          if (count == 0 || _data == nullptr)
        	  throw std::invalid_argument("_data is nullptr");

//          volatile bool swap = false;
//          volatile bool appended = false;

          uint8_t* start = data.get();

          // if fails, then already full
          if (!try_lock_write()) {  // lock checks if buffer is being read/full, and then increment writer count.
            //printf("full     shared pointer %p, count %d, start %p, max_pointer %p, size %d, writers %d\n", (uint8_t*)pointer, count, data.get(), max_pointer, (int32_t)size, (int)writerCount);
            return 0x0;
          }


            // After the F&A, there are 3 types of threads:  target write area is
            //  1. completely within buffer:  these threads should proceed to memcpy with the local ptr var
            //  2. completely outside buffer: these threads will not memcpy.  they all disabled buffer, and
            //     only 1 thread from 2) or 3) should swap in a new buffer atomically.
            //      - if curr buffer is not disabled (already swapped), then don't swap further.
            //  3. crossing buffer boundary: this thread will not memcpy.  disabled.  and it should retract the pointer advance.


            // ONE thread at a time can get the position to copy in data.  If full, then the buffer is full.
            std::atomic_thread_fence(std::memory_order_seq_cst);    // compute ptr strictly after determining can lock write.
            uint8_t* ptr = pointer.fetch_add(count, std::memory_order_seq_cst);  // no memory ordering needed within mutex loc
            std::atomic_thread_fence(std::memory_order_seq_cst);    // branch and memcpy only after computing the ptr.

            if ((this->max_pointer - ptr) >= static_cast<std::ptrdiff_t>(count) ) {
              // has room to copy.


              std::memcpy(ptr, _data, count);

              std::atomic_thread_fence(std::memory_order_seq_cst);  // unlock only after memcpy.
              //appended = true;
              unlock_write();

              atomic_thread_fence(std::memory_order_seq_cst);   // commits all changes up to this point to memory

              return 0x1;
            } else {
              lock_read2();  // any full buffer can lock the read.
              unlock_write();  // then it unlock the write.
              std::atomic_thread_fence(std::memory_order_seq_cst);  // ensure these are completed.


              if (ptr <= this->max_pointer) {  // thread that made the buffer full.
                // now wait for all other writes to complete
                std::atomic_thread_fence(std::memory_order_seq_cst);  // only update size after all writes are done.

                wait_for_all_writes();

                std::atomic_thread_fence(std::memory_order_seq_cst);  // only update size after all writes are done.

                this->size.store(ptr - data.get(), std::memory_order_seq_cst);  // overwrite the size with this one's

//                swap = true;
                atomic_thread_fence(std::memory_order_seq_cst);   // commits all changes up to this point to memory


                if ((size.load() / sizeof(int)) != (capacity / sizeof(int))) fprintf(stderr, "IN BUFFER:  NOT 2047 elements. got %ld. other info: data %p, approxsize=%d, buffer read locked? %s  buffer write lock %s\n", size.load() / sizeof(int), data.get(), getApproximateSize<TS>(), isReading() ? "yes" : "no", isWriting() ? "yes": "no");

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
         * Semantically, a "false" return value means there is not enough room for the new data, not that the buffer is full.
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
          if (!try_lock_write()) return 0x0;

          if ((pointer + count) > max_pointer) {
            lock_read2();  // also takes away an writer.
            unlock_write();

            if (pointer <= max_pointer)  {// original pointer is within buffer - put it back. however, this could cause multiple threads to reach this point.
              wait_for_all_writes();

              size = pointer - data.get(); // final
              return 0x2;
            } else
              return 0x0;
          } else {


            std::memcpy(pointer, _data, count);
            pointer += count;
            unlock_write();

            return 0x1;
          }
        }



    };

  } /* namespace io */
} /* namespace bliss */

#endif /* BUFFER_HPP_ */
