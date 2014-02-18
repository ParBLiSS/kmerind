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

        io_exception(std::string& _msg)
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

        typedef bliss::iterator::RangeType<size_t> range_type;

        /**
         * opens the file and save in file handle
         */
        file_loader(std::string const & filename)
          : range(), data(0)
        {
          file_handle = open(filename.c_str(), O_RDONLY);
          page_size = sysconf(_SC_PAGE_SIZE);
        };

        /**
         * closes the file.
         */
        virtual ~file_loader()
        {
          close(file_handle);
        };

        virtual void map() throw(io_exception) = 0;
        virtual void unmap() throw(io_exception) = 0;

        char* getData()
        {
          return data;
        };

        void setRange(bliss::iterator::range<size_t> const &_range) {
          range = _range;
        }

        bliss::iterator::range<size_t>& getRange() const
        {
          return range;
        };


      protected:
        int file_handle;// file handle

        range_type range;  // offset in file from where to read

        char* data;     // memmapped data.

        size_t page_size;
    };

  } /* namespace functional */
} /* namespace bliss */
#endif /* FILE_LOADER_HPP_ */
