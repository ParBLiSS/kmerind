/**
 * @file		BufferPool.hpp
 * @ingroup
 * @author	tpan
 * @brief   a pool of memory buffers.
 * @details provides a reusable pool of buffers.  Each buffer is a block of preallocated memory that can be copied into.
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef BUFFERPOOL_HPP_
#define BUFFERPOOL_HPP_

#include <string.h>

#include <mutex>

namespace bliss
{
  namespace io
  {

    /**
     * @class			bliss::io::BufferPool
     * @brief
     * @details
     *
     */
    template<bool THREAD_SAFE = true>
    class BufferPool
    {
      public:
        BufferPool();
        virtual ~BufferPool();
    };

  } /* namespace io */
} /* namespace bliss */

#endif /* BUFFERPOOL_HPP_ */
