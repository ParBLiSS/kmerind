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
 * file_loader_test.cpp
 *   No separate MPI test.  FileLoader has no logical differences between MPI and non-MPI versions.  Only difference is MPI_comm is stored for convenience.
 *   note that this may change in the future if MPI_IO is used.
 *
 *  Created on: Feb 18, 2014
 *      Author: Tony Pan <tpan7@gatech.edu>
 */


#include "io/test/file_loader_mpi_test_fixture.cpp"

// TODO: need to test mmap with multibyte elements.

typedef ::testing::Types<
    FileLoader<unsigned char, 0, bliss::io::BaseFileParser, false, false>,
    FileLoader<char,          0, bliss::io::BaseFileParser, false, false>,
    FileLoader<int16_t,       0, bliss::io::BaseFileParser, false, false>,
    FileLoader<uint16_t,      0, bliss::io::BaseFileParser, false, false>,
    FileLoader<int32_t,       0, bliss::io::BaseFileParser, false, false>,
    FileLoader<uint32_t,      0, bliss::io::BaseFileParser, false, false>,
    FileLoader<int64_t,       0, bliss::io::BaseFileParser, false, false>,
    FileLoader<uint64_t,      0, bliss::io::BaseFileParser, false, false>,
    FileLoader<unsigned char, 0, bliss::io::BaseFileParser, false, false>,
    FileLoader<unsigned char, 0, bliss::io::BaseFileParser, false, true>,
    FileLoader<unsigned char, 0, bliss::io::BaseFileParser, true, false>,
    FileLoader<unsigned char, 0, bliss::io::BaseFileParser, true, true>,
    FileLoader<unsigned char, 1, bliss::io::BaseFileParser, false, false>,
    FileLoader<unsigned char, 2, bliss::io::BaseFileParser, false, false>,
    FileLoader<unsigned char, 3, bliss::io::BaseFileParser, false, false>,
    FileLoader<unsigned char, 4, bliss::io::BaseFileParser, false, false>
> FileLoaderTestDefaultPartitionerTypes;
INSTANTIATE_TYPED_TEST_CASE_P(Bliss_Default, FileLoaderTest, FileLoaderTestDefaultPartitionerTypes);

