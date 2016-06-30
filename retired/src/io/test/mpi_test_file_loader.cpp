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
    FileLoader<unsigned char, 0, bliss::io::BaseFileParser, false, false, bliss::partition::CyclicPartitioner<bliss::partition::range<size_t> >, bliss::partition::BlockPartitioner<bliss::partition::range<size_t> > >,
    FileLoader<unsigned char, 0, bliss::io::BaseFileParser, false, false, bliss::partition::BlockPartitioner<bliss::partition::range<size_t> >, bliss::partition::BlockPartitioner<bliss::partition::range<size_t> > >,
    FileLoader<unsigned char, 0, bliss::io::BaseFileParser, false, false, bliss::partition::DemandDrivenPartitioner<bliss::partition::range<size_t> >, bliss::partition::CyclicPartitioner<bliss::partition::range<size_t> > >,
    FileLoader<unsigned char, 0, bliss::io::BaseFileParser, false, false, bliss::partition::CyclicPartitioner<bliss::partition::range<size_t> >, bliss::partition::CyclicPartitioner<bliss::partition::range<size_t> > >,
    FileLoader<unsigned char, 0, bliss::io::BaseFileParser, false, false, bliss::partition::BlockPartitioner<bliss::partition::range<size_t> >, bliss::partition::CyclicPartitioner<bliss::partition::range<size_t> > >,
    FileLoader<unsigned char, 0, bliss::io::BaseFileParser, false, false, bliss::partition::DemandDrivenPartitioner<bliss::partition::range<size_t> >, bliss::partition::DemandDrivenPartitioner<bliss::partition::range<size_t> > >,
    FileLoader<unsigned char, 0, bliss::io::BaseFileParser, false, false, bliss::partition::CyclicPartitioner<bliss::partition::range<size_t> >, bliss::partition::DemandDrivenPartitioner<bliss::partition::range<size_t> > >,
    FileLoader<unsigned char, 0, bliss::io::BaseFileParser, false, false, bliss::partition::BlockPartitioner<bliss::partition::range<size_t> >, bliss::partition::DemandDrivenPartitioner<bliss::partition::range<size_t> > >,
    FileLoader<unsigned char, 4, bliss::io::BaseFileParser, false, false, bliss::partition::CyclicPartitioner<bliss::partition::range<size_t> >, bliss::partition::BlockPartitioner<bliss::partition::range<size_t> > >,
    FileLoader<unsigned char, 4, bliss::io::BaseFileParser, false, false, bliss::partition::BlockPartitioner<bliss::partition::range<size_t> >, bliss::partition::BlockPartitioner<bliss::partition::range<size_t> > >,
    FileLoader<unsigned char, 4, bliss::io::BaseFileParser, false, false, bliss::partition::DemandDrivenPartitioner<bliss::partition::range<size_t> >, bliss::partition::CyclicPartitioner<bliss::partition::range<size_t> > >,
    FileLoader<unsigned char, 4, bliss::io::BaseFileParser, false, false, bliss::partition::CyclicPartitioner<bliss::partition::range<size_t> >, bliss::partition::CyclicPartitioner<bliss::partition::range<size_t> > >,
    FileLoader<unsigned char, 4, bliss::io::BaseFileParser, false, false, bliss::partition::BlockPartitioner<bliss::partition::range<size_t> >, bliss::partition::CyclicPartitioner<bliss::partition::range<size_t> > >,
    FileLoader<unsigned char, 4, bliss::io::BaseFileParser, false, false, bliss::partition::DemandDrivenPartitioner<bliss::partition::range<size_t> >, bliss::partition::DemandDrivenPartitioner<bliss::partition::range<size_t> > >,
    FileLoader<unsigned char, 4, bliss::io::BaseFileParser, false, false, bliss::partition::CyclicPartitioner<bliss::partition::range<size_t> >, bliss::partition::DemandDrivenPartitioner<bliss::partition::range<size_t> > >,
    FileLoader<unsigned char, 4, bliss::io::BaseFileParser, false, false, bliss::partition::BlockPartitioner<bliss::partition::range<size_t> >, bliss::partition::DemandDrivenPartitioner<bliss::partition::range<size_t> > >

> FileLoaderTestTypes;
INSTANTIATE_TYPED_TEST_CASE_P(Bliss, FileLoaderTest, FileLoaderTestTypes);

