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


#include "bliss-config.hpp"    // for location of data.

#if defined(USE_MPI)
#include "mpi.h"
#endif

// include google test
#include <gtest/gtest.h>
#include <cstdint> // for uint64_t, etc.
#include <string>

#include "io/file_loader.hpp"
#include "io/file.hpp"

#include "io/test/file_loader_test_fixtures.hpp"


using namespace bliss::io;




// indicate this is a typed test
TYPED_TEST_CASE_P(FileLoadProcedureTest);


// normal test cases
TYPED_TEST_P(FileLoadProcedureTest, open)
{

#ifdef USE_MPI
	int rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0) {
#endif

  typedef TypeParam FileLoaderType;
  constexpr size_t overlap = FileLoaderType::get_overlap_size();
  typedef typename FileLoaderType::RangeType RangeType;

  // get fileName

  FileLoaderType loader(this->fileName, 1, 0, 1, 64 * 1024, 64 * 1024);   // NOTE: larger block is MUCH faster.  2 * threads * 2048 bytes take a LONG time because of mmap each time...
  RangeType fr = loader.getFileRange();

  auto block = loader.getNextL1Block();
  RangeType r = block.getRange();

  // start at file starting point
  ASSERT_EQ(fr.start, r.start);
  ASSERT_EQ(fr.end, r.end);

  bool same = true;

  ::bliss::io::MemData data = ::bliss::io::memmap::load_file(this->fileName, overlap);

  ASSERT_TRUE(data.mem_range.size() == r.size());

  same = std::equal(block.begin(), block.end(), data.data);

  ::bliss::io::unload_data(data);

    ASSERT_TRUE(same);

#ifdef USE_MPI
	}
#endif
}


#ifdef USE_MPI
TYPED_TEST_P(FileLoadProcedureTest, open_mpi)
{
  typedef TypeParam                          FileLoaderType;
  constexpr size_t overlap = FileLoaderType::get_overlap_size();
  typedef typename FileLoaderType::RangeType RangeType;

  // get this->fileName

  FileLoaderType loader(this->fileName, MPI_COMM_WORLD, 1, 64 * 1024, 64 * 1024 );

  auto block = loader.getNextL1Block();
  RangeType r = block.getRange();
  bool same = true;

  ::bliss::io::MemData data = ::bliss::io::parallel::memmap::load_file(this->fileName, MPI_COMM_WORLD,  overlap);

  if (data.valid_range.size() != r.size()) {
	  std::cout << "file open: " << data.mem_range << ::std::endl;
	  std::cout << "file open: " << data.valid_range << ::std::endl;
	  std::cout << "file loader: " << r << ::std::endl;
  }

  ASSERT_TRUE(data.valid_range.size() == r.size());

  same = std::equal(block.begin(), block.end(), data.data);

  ::bliss::io::unload_data(data);

  ASSERT_TRUE(same);


  MPI_Barrier(MPI_COMM_WORLD);

}
#endif




// now register the test cases
REGISTER_TYPED_TEST_CASE_P(FileLoadProcedureTest,
#ifdef USE_MPI
		open_mpi,
#endif
 open
);


// TODO: need to test mmap with multibyte elements.

typedef ::testing::Types<
    FileLoader<unsigned char, 0, bliss::io::BaseFileParser, false, false, bliss::partition::BlockPartitioner<bliss::partition::range<size_t> >, bliss::partition::BlockPartitioner<bliss::partition::range<size_t> > >,
    FileLoader<unsigned char, 0, bliss::io::BaseFileParser, false,  true, bliss::partition::BlockPartitioner<bliss::partition::range<size_t> >, bliss::partition::BlockPartitioner<bliss::partition::range<size_t> > >,
    FileLoader<unsigned char, 0, bliss::io::BaseFileParser,  true, false, bliss::partition::BlockPartitioner<bliss::partition::range<size_t> >, bliss::partition::BlockPartitioner<bliss::partition::range<size_t> > >,
    FileLoader<unsigned char, 0, bliss::io::BaseFileParser,  true,  true, bliss::partition::BlockPartitioner<bliss::partition::range<size_t> >, bliss::partition::BlockPartitioner<bliss::partition::range<size_t> > >,
    FileLoader<unsigned char, 1, bliss::io::BaseFileParser, false,  true, bliss::partition::BlockPartitioner<bliss::partition::range<size_t> >, bliss::partition::BlockPartitioner<bliss::partition::range<size_t> > >,
    FileLoader<unsigned char, 2, bliss::io::BaseFileParser, false,  true, bliss::partition::BlockPartitioner<bliss::partition::range<size_t> >, bliss::partition::BlockPartitioner<bliss::partition::range<size_t> > >,
    FileLoader<unsigned char, 3, bliss::io::BaseFileParser, false,  true, bliss::partition::BlockPartitioner<bliss::partition::range<size_t> >, bliss::partition::BlockPartitioner<bliss::partition::range<size_t> > >,
    FileLoader<unsigned char, 4, bliss::io::BaseFileParser, false,  true, bliss::partition::BlockPartitioner<bliss::partition::range<size_t> >, bliss::partition::BlockPartitioner<bliss::partition::range<size_t> > >
> FileLoadProcedureTestTypes;

//typedef ::testing::Types< FileLoader<unsigned char, 0, bliss::io::BaseFileParser, false, false> > FileLoadProcedureTestTypes;

INSTANTIATE_TYPED_TEST_CASE_P(Bliss, FileLoadProcedureTest, FileLoadProcedureTestTypes);


int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);



#if defined(USE_MPI)
  MPI_Init(&argc, &argv);
#endif

  result = RUN_ALL_TESTS();

#if defined(USE_MPI)
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
#endif

  return result;
}



