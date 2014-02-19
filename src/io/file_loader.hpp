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

        typedef bliss::iterator::range<size_t> range_type;

        /**
         * opens the file and save in file handle
         */
        file_loader(std::string const & filename)
          : range(), data(std::nullptr)
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

        /**
         * return data pointed to by range.start.  not same as range.block_start.
         */
        char* getData()
        {
          return data + (range.start - range.block_start);
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

        char* data;
        char* aligned_data;  // memmapped data, page aligned.

        size_t page_size;

        /**
         * TODO: test on remotely mounted file system.
         */

        /**
         * TODO: change this so file_loader calls this.
         */
        char* map(file_loader::range_type & r) throw(io_exception) {
          assert(r.is_page_aligned(page_size), "ERROR: range is not page aligned.");

          char* output = (char*)mmap(nullptr, r.end - r.block_start, PROT_READ, MAP_PRIVATE, file_handle, r.block_start );

          if (output == MAP_FAILED)
          {
            std::stringstream ss;
            int myerr = errno;
            ss << "ERROR in mmap: "  << myerr << ": " << strerror(myerr);
            throw io_exception(ss.str());
          }
        }

        void unmap(char * input, file_loader::range_type & r) throw(io_exception) {
          assert(r.is_page_aligned(page_size), "ERROR: range is not page aligned.");

          munmap(input, r.end - r.block_start);
        }






    };

  } /* namespace functional */
} /* namespace bliss */
#endif /* FILE_LOADER_HPP_ */
