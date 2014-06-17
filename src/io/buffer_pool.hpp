/**
 * @file		buffer_pool.hpp
 * @ingroup
 * @author	tpan
 * @brief   a pool of memory buffers.
 * @details provides a reusable pool of buffers.  Each buffer is a block of preallocated memory that can be copied into.
 *          this dynamically grows if capacity is not specified.
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef BUFFERPOOL_HPP_
#define BUFFERPOOL_HPP_


#include <mutex>
#include <queue>
#include <unordered_set>

#include "io/buffer.hpp"
#include "concurrent/concurrent.hpp"


namespace bliss
{
  namespace io
  {

    static constexpr size_t UNLIMITED = 0;

    template<bliss::concurrent::ThreadSafety TS>
    class BufferPool;


    /**
     * @class			bliss::io::BufferPool
     * @brief     thread-safe version of buffer pool
     * @details   this can only contain thread-safe version of buffers.
     *
     */
    template<>
    class BufferPool<bliss::concurrent::THREAD_SAFE>
    {
      public:
        typedef bliss::io::Buffer<bliss::concurrent::THREAD_SAFE>     BufferType;

      protected:

        size_t capacity;

        std::queue<BufferType> available;
        std::unordered_set<BufferType> inUse;

      public:
        BufferPool(const size_t _capacity = bliss::io::UNLIMITED) : capacity(_capacity) {
          if (_capacity != bliss::io::UNLIMITED) {
            available.
          }
        };
        virtual ~BufferPool() {};

        BufferType& getBuffer() {

        }

        BufferType&
    };



    /**
     * @class     bliss::io::BufferPool
     * @brief     non-thread-safe version of buffer pool
     * @details   this can only contain non-thread-safe version of buffers.
     *
     */
    template<>
    class BufferPool<bliss::concurrent::THREAD_UNSAFE>
    {
      protected:
        size_t capacity;

      public:
        BufferPool(const size_t capacity = bliss::io::UNLIMITED) : capacity(_capacity) {};
        virtual ~BufferPool() {};
    };


  } /* namespace io */
} /* namespace bliss */

#endif /* BUFFERPOOL_HPP_ */
