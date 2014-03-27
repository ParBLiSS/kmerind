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
#include <sys/mman.h>
#include <fcntl.h> // for open
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>

#include <type_traits>

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

        io_exception(const std::string &_msg)
        {
          message = _msg;
        }

        virtual const char* what() const throw ()
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
        file_loader(const std::string &_filename, bool preload = false)
            : filename(_filename), data(nullptr), aligned_data(nullptr),
              preloaded(preload)
        {
          //printf("FILE_LOADER called.  not loading.\n");

          file_handle = open(filename.c_str(), O_RDONLY);
          page_size = sysconf(_SC_PAGE_SIZE);
        }
        ;

        file_loader(const std::string &_filename, const range_type &_range,
                    bool preload = false)
            : filename(_filename), range(_range), data(nullptr),
              aligned_data(nullptr), preloaded(preload)
        {
          //printf("FILE_LOADER called with range.  loading..\n");

          file_handle = open(filename.c_str(), O_RDONLY);
          page_size = sysconf(_SC_PAGE_SIZE);

          // now mmap the file
          load();
        }
        ;

        /**
         * closes the file.
         */
        virtual ~file_loader()
        {
          //printf("~FILE_LOADER called.\n");

          unload();

          //printf("~FILE_LOADER unloaded.\n");

          close(file_handle);
        }
        ;

        /**
         * return data pointed to by range.start.  not same as range.block_start.
         */
        char* getData()
        {
          return data;
        }
        ;

        range_type getRange() const
        {
          return range;
        }
        ;

      protected:
        std::string filename;
        range_type range;  // offset in file from where to read

        size_t page_size;
        int file_handle;  // file handle

        char* data;
        char* aligned_data;  // memmapped data, page aligned.  strictly internal

        bool preloaded;

        /**
         * TODO: test on remotely mounted file system.
         */

        /**
         *
         */
        void load() throw (io_exception)
        {
          range = range.align_to_page(page_size);

          aligned_data = (char*)mmap(nullptr, range.end - range.block_start,
                                     PROT_READ,
                                     MAP_PRIVATE, file_handle,
                                     range.block_start);

          if (aligned_data == MAP_FAILED)
          {
            std::stringstream ss;
            int myerr = errno;
            ss << "ERROR in mmap: " << myerr << ": " << strerror(myerr);
            throw io_exception(ss.str());
          }

          if (preloaded)
          {
            // allocate space
            data = new char[range.end - range.start + 1];
            //copy data over
            memcpy(data, aligned_data + (range.start - range.block_start),
                   sizeof(char) * (range.end - range.start));

            data[range.end - range.start] = 0;

            // close the input
            munmap(aligned_data, range.end - range.block_start);

            aligned_data = data;

          }
          else
          {

            data = aligned_data + (range.start - range.block_start);

          }
        }

        void unload()
        {
          if (preloaded)
          {
            if (data != nullptr)
              delete [] data;
          }
          else
          {
            if (aligned_data != nullptr || aligned_data != MAP_FAILED)
              munmap(aligned_data, range.end - range.block_start);
          }
          aligned_data = nullptr;
          data = nullptr;

          printf("unloading complete.\n");
        }

    };

  } /* namespace functional */
} /* namespace bliss */
#endif /* FILE_LOADER_HPP_ */
