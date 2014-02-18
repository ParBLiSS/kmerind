/**
 * fastq_loader.cpp
 *
 *  Created on: Feb 18, 2014
 *      Author: tpan
 */

#include <sys/mman.h>
#include <cstdlib>  // for itoa.
#include <cerrno>

#include "io/fastq_loader.hpp"


namespace bliss
{
  namespace io
  {

    fastq_loader::~fastq_loader()
    {
      if (data != nullptr)
        munmap(data, range.length)
    }

    fastq_loader::fastq_loader(std::string const & filename,
                               file_loader::range_type const & range)
     : file_loader(filename)
    {
      // do something (store/compute) the range.
    }

    void fastq_loader::map() throw(io_exception) {
      data = (char*)mmap(nullptr, range.length, PROT_READ, MAP_PRIVATE, file_handle, range.offset );

      if (data == MAP_FAILED)
      {
        int myerr = errno;
        std::string msg = "ERROR in mmap: " + itoa(myerr) + ": " + strerror(myerr);
        throw io_exception(msg);
      }
    }

    void fastq_loader::unmap() throw(io_exception) {
      munmap(data, range.length);
    }




  } /* namespace io */
} /* namespace bliss */
