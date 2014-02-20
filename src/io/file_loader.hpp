/**
 * file_loader.hpp
 *
 *  Created on: Feb 18, 2014
 *      Author: tpan
 */

#ifndef FILE_LOADER_HPP_
#define FILE_LOADER_HPP_

#include <string>
#include <cstring>
#include <exception>

#include <unistd.h>  // sysconf

#include <sys/stat.h>
#include <fcntl.h> // for open

#include "iterators/range.hpp"


namespace bliss
{
  namespace io
  {

    class io_exception : public std::exception
    {
      protected:
        std::string message;

      public:
        io_exception(const char* _msg)
        {
          message = _msg;
        }

        io_exception(std::string const & _msg)
        {
          message = _msg;
        }


        virtual const char* what() const throw()
        {
          return message.c_str();
        }
    };

    /**
     *
     */
    class file_loader
    {
      public:

        typedef bliss::iterator::range<size_t> range_type;

        /**
         * opens the file and save in file handle
         */
        file_loader()
        {
          page_size = sysconf(_SC_PAGE_SIZE);
        };


        file_loader(std::string const & _filename, range_type const & _range)
          : filename(_filename), range(_range)
        {
          file_handle = open(filename.c_str(), O_RDONLY);
          page_size = sysconf(_SC_PAGE_SIZE);

          // now mmap the file
          map();
        };

        /**
         * closes the file.
         */
        virtual ~file_loader()
        {
          unmap();
          close(file_handle);
        };

        /**
         * return data pointed to by range.start.  not same as range.block_start.
         */
        char* getData()
        {
          return data;
        };

        range_type getRange() const
        {
          return range;
        };


      protected:
        std::string filename;
        range_type range;  // offset in file from where to read

        size_t page_size;
        int file_handle;// file handle

        char* data;
        char* aligned_data;  // memmapped data, page aligned.  strictly internal


        /**
         * TODO: test on remotely mounted file system.
         */

        /**
         *
         */
        void map() throw(io_exception) {
          range = range.align_to_page(page_size);

          aligned_data = (char*)mmap(nullptr,
                                     range.end - range.block_start,
                                     PROT_READ, MAP_PRIVATE,
                                     file_handle, range.block_start );

          if (aligned_data == MAP_FAILED)
          {
            std::stringstream ss;
            int myerr = errno;
            ss << "ERROR in mmap: "  << myerr << ": " << strerror(myerr);
            throw io_exception(ss.str());
          }

          data = aligned_data + (range.start - range.block_start);
        }

        void unmap() throw(io_exception) {
          if (aligned_data != nullptr || aligned_data != MAP_FAILED)
            munmap(aligned_data, range.end - range.block_start);
        }






    };

  } /* namespace functional */
} /* namespace bliss */
#endif /* FILE_LOADER_HPP_ */
