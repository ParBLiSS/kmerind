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


#include <fcntl.h>      // for open
#include <sys/stat.h>   // block size.
#include <unistd.h>     // sysconf
#include <sys/mman.h>   // mmap
#include <cstring>      // memcpy, strerror

#include <cassert>
#include "io/file_loader.hpp"
#include "partition/range.hpp"


int main(int argc, char** argv) {

  assert(argc > 0);

  std::string filename;
  filename.assign(argv[1]);



  /// get file size.
  struct stat filestat;
  int ret = stat(filename.c_str(), &filestat);
  if (ret < 0 ) {
    ERROR( "ERROR in file open to get size." );
    exit(-1);
  }
  size_t file_size = static_cast<size_t>(filestat.st_size);


  /// open the file and get a handle.
  int file_handle = open(filename.c_str(), O_RDONLY);
  if (file_handle == -1)
  {
    int myerr = errno;
    FATAL("ERROR in file open: ["  << filename << "] error " << myerr << ": " << strerror(myerr));
  }


  // mmap
  unsigned char* data = (unsigned char*)mmap(nullptr, file_size,
                                       PROT_READ,
                                       MAP_PRIVATE, file_handle,
                                       0);



  std::cout << "please enter start and end positions.  enter a negative number to exit." << std::endl;
  ///
  int pos = 0;
  int endpos = 0;
  do {
    pos = 0;
    endpos = 0;
    std::cin >> pos >> endpos;
    if (pos >= 0) {
      std::cout << "char at pos " << pos << " - " << endpos << std::endl;
      std::ostream_iterator<unsigned char> oit(std::cout);
      std::copy(data + pos, data + endpos, oit);
      std::cout << std::endl;
    }
  } while (pos >= 0);

  ///
  munmap(data, file_size);


}
