/*
 * Copyright 2015 Georgia Institute of Technology
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/**
 * @file		TextInspector.cpp
 * @ingroup
 * @author	Tony Pan <tpan7@gatech.edu>
 * @brief
 * @details
 *
 */


#include <string>
#include <iostream>
#include <iterator>  // for ostream_iterator
#include <cstdlib>   // for strtoul
#include <cstdio>    // for fopen, fseek, fread

#include <fcntl.h>      // for open
#include <sys/stat.h>   // block size.
#include <unistd.h>     // sysconf
#include <sys/mman.h>   // mmap
#include <cstring>      // memcpy, strerror
#include <limits>

#include <cassert>
#include "io/file_loader.hpp"
#include "partition/range.hpp"

#include "tclap/CmdLine.h"


int main(int argc, char** argv) {

  // Get the value parsed by each arg.
  std::string filename;
  size_t pos = 0;
  size_t end_pos = ::std::numeric_limits<size_t>::max();
  bool map = false;

  // Wrap everything in a try block.  Do this every time,
  // because exceptions will be thrown for problems.
  try {

    // Define the command line object, and insert a message
    // that describes the program. The "Command description message"
    // is printed last in the help text. The second argument is the
    // delimiter (usually space) and the last one is the version number.
    // The CmdLine object parses the argv array based on the Arg objects
    // that it contains.
    TCLAP::CmdLine cmd("Extracts text within an offset range from a file", ' ', "0.1");

    // Define a value argument and add it to the command line.
    // A value arg defines a flag and a type of value that it expects,
    // such as "-n Bishop".
    TCLAP::ValueArg<std::string> fileArg("f", "file", "Text file path", true, "", "string", cmd);

    TCLAP::ValueArg<size_t> beginArg("b", "begin", "Start of range of file to extract", false, 0UL, "unsigned long int", cmd);
    TCLAP::ValueArg<size_t> endArg("e", "end", "End of range of file to extract", false, ::std::numeric_limits<size_t>::max(), "unsigned long int", cmd);


    // Define a switch and add it to the command line.
    // A switch arg is a boolean argument and only defines a flag that
    // indicates true or false.  In this example the SwitchArg adds itself
    // to the CmdLine object as part of the constructor.  This eliminates
    // the need to call the cmd.add() method.  All args have support in
    // their constructors to add themselves directly to the CmdLine object.
    // It doesn't matter which idiom you choose, they accomplish the same thing.
    TCLAP::SwitchArg mmapSwitch("m", "mmap","Use MMAP for file reading.  (else fread)", cmd, false);

    // Parse the argv array.
    cmd.parse( argc, argv );

    // Get the value parsed by each arg.
    filename = fileArg.getValue();
    pos = beginArg.getValue();
    end_pos = endArg.getValue();
    map = mmapSwitch.getValue();

    // Do what you intend.


  } catch (TCLAP::ArgException &e)  // catch any exceptions
  {
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    exit(-1);
  }


  // if no file, will die.

  /// get file size.
  struct stat filestat;
  int ret = stat(filename.c_str(), &filestat);
  if (ret < 0 ) {
    BL_FATAL( "ERROR in file open to get size." );
  }
  size_t file_size = static_cast<size_t>(filestat.st_size);

  pos = ::std::min(pos, file_size);
  end_pos = ::std::min(end_pos, file_size);

  size_t size = end_pos - pos;

  if (map) {
    /// open the file and get a handle.
    int file_handle = open(filename.c_str(), O_RDONLY);
    if (file_handle == -1)
    {
      int myerr = errno;
      BL_FATAL("ERROR in file open: ["  << filename << "] error " << myerr << ": " << strerror(myerr));
      USED_BY_LOGGER_ONLY(myerr);
    }

    size_t page_size = sysconf(_SC_PAGE_SIZE);
    size_t map_offset = pos % page_size;
    size_t map_pos = pos - map_offset;

    // mmap
    unsigned char* data = (unsigned char*)mmap64(nullptr, end_pos - map_pos,
                                         PROT_READ,
                                         MAP_PRIVATE, file_handle,
                                         map_pos);

    std::ostream_iterator<unsigned char> oit(std::cout);
    std::copy(data + map_offset, data + map_offset + size, oit);
    std::cout << std::endl;

    ///
    munmap(data, end_pos - map_pos);
  } else {
    FILE * f = fopen(filename.c_str(), "r");
    fseek(f, pos, SEEK_SET );
    unsigned char * data = new unsigned char[size + 1];

    (void) fread(data, sizeof(unsigned char), size, f);

    std::ostream_iterator<unsigned char> oit(std::cout);
    std::copy(data, data + size, oit);
    std::cout << std::endl;

    fclose(f);
  }
}
