/**
 * @file		TextInspector.cpp
 * @ingroup
 * @author	Tony Pan <tpan7@gatech.edu>
 * @brief
 * @details
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */


#include <string>
#include <iostream>
#include <iterator>  // for ostream_iterator
#include <cstdlib>   // for strtoul

#include <fcntl.h>      // for open
#include <sys/stat.h>   // block size.
#include <unistd.h>     // sysconf
#include <sys/mman.h>   // mmap
#include <cstring>      // memcpy, strerror

#include <cassert>
#include "io/file_loader.hpp"
#include "partition/range.hpp"


int main(int argc, char** argv) {

  assert(argc > 1);

  std::string filename;
  filename.assign(argv[1]);



  /// get file size.
  struct stat filestat;
  int ret = stat(filename.c_str(), &filestat);
  if (ret < 0 ) {
    BL_ERROR( "ERROR in file open to get size." );
    exit(-1);
  }
  size_t file_size = static_cast<size_t>(filestat.st_size);



  /// open the file and get a handle.
  int file_handle = open(filename.c_str(), O_RDONLY);
  if (file_handle == -1)
  {
    int myerr = errno;
    BL_FATAL("ERROR in file open: ["  << filename << "] error " << myerr << ": " << strerror(myerr));
    USED_BY_LOGGER_ONLY(myerr);
  }

  size_t pos = 0;
  size_t end_pos = file_size;

  if (argc > 2) {
    pos = ::std::strtoul(argv[2], NULL, 0);
    pos = ::std::min(pos, file_size);
  }
  if (argc > 3) {
    end_pos = ::std::strtoul(argv[3], NULL, 0);
    end_pos = ::std::min(end_pos, file_size);
  }

  size_t map_size = end_pos - pos;

  size_t page_size = sysconf(_SC_PAGE_SIZE);
  size_t map_offset = pos % page_size;
  size_t map_pos = pos - map_offset;

  // mmap
  unsigned char* data = (unsigned char*)mmap64(nullptr, end_pos - map_pos,
                                       PROT_READ,
                                       MAP_PRIVATE, file_handle,
                                       map_pos);

  std::ostream_iterator<unsigned char> oit(std::cout);
  std::copy(data + map_offset, data + map_offset + map_size, oit);
  std::cout << std::endl;

  ///
  munmap(data, end_pos - map_pos);

}
