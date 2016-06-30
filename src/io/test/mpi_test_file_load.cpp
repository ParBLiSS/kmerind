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
#include "mxx/env.hpp"


TODO: test file load direct, serially. - just with base file types.  no parsers needed
verify that file is read correctly

TODO: test file load direct, MPI.  test with 3 parsers.  test with different overlaps for basefileparser and fastaparser
first verify that overlaps are correct.
next verify that there is complete coverage
finally verify that file is read correctly.

NEXT: test file load and parse to sequence  (fasta and fastq only.)  pick POSIX for this test.  k = 31.
FINALLY:  test k-mer parse (fasta and fastq only).  pick POSIX for this test.  k = 31



using namespace bliss::io;


typedef BaseFileParser<unsigned char *> ParserType;


class FileLoadProcedureTest : public FileLoaderTest
{
protected:
	~FileLoadProcedureTest() {};

	template <typename file_loader>
	void open_seq(file_loader & fobj) {

		  using RangeType = typename file_loader::range_type;

		  // get fileName

			bliss::io::file_data fdata = fobj.read_file();
			ASSERT_EQ(fdata.getRange().size(), fobj.get_file_size());
			ASSERT_EQ(fdata.getRange().end, fobj.get_file_size());

			ASSERT_EQ(fdata.in_mem_range_bytes.size(), fobj.get_file_size());
			ASSERT_EQ(fdata.in_mem_range_bytes.end, fobj.get_file_size());

			ASSERT_EQ(fdata.parent_range_bytes.size(), fobj.get_file_size());
			ASSERT_EQ(fdata.parent_range_bytes.end, fobj.get_file_size());

		  bool same = true;



		  same = std::equal(block.begin(), block.end(), fdata.begin());

		    ASSERT_TRUE(same);
	}

	template <typename file_loader>
	void open_mpi(file_loader & fobj, mxx::comm const & comm) {
		  using RangeType = ::bliss::partition::range<size_t>;

		  // get this->fileName

		  FileLoaderType loader(this->fileName, comm, 1, 64 * 1024, 64 * 1024 );

		  auto block = loader.getNextL1Block();
		  RangeType r = block.getRange();
		  bool same = true;

			::bliss::io::file_data fdata = fobj.read_file();

		  if (fdata.valid_range_bytes.size() != r.size()) {
			  std::cout << "file filename: " << fobj.get_filename() << std::endl;
			  std::cout << "file open in_mem: " << fdata.in_mem_range_bytes << ::std::endl;
			  std::cout << "file open valid: " << fdata.valid_range_bytes << ::std::endl;
			  std::cout << "file open parent: " << fdata.parent_range_bytes << ::std::endl;

			  std::cout << "this filename: " << this->fileName << ::std::endl;
			  std::cout << "fileloader filename: " << loader.getFilename() << ::std::endl;
			  std::cout << "fileloader block: " << r << ::std::endl;

		  }

		  ASSERT_TRUE(fdata.valid_range_bytes.size() == r.size());

		  same = std::equal(block.begin(), block.end(), fdata.begin());

		  ASSERT_TRUE(same);

	}
};


// indicate this is a typed test
TYPED_TEST_CASE_P(FileLoadProcedureTest);


// normal test cases
TYPED_TEST_P(FileLoadProcedureTest, open_mmap)
{
#ifdef USE_MPI
	::mxx::comm comm;
	if (comm.rank() == 0) {
#endif

		bliss::io::mmap_file fobj(this->fileName);

		this->open_seq(fobj);

#ifdef USE_MPI
	}
#endif
}

// normal test cases
TYPED_TEST_P(FileLoadProcedureTest, open_stdio)
{
#ifdef USE_MPI
	::mxx::comm comm;
	if (comm.rank() == 0) {
#endif

		bliss::io::stdio_file fobj(this->fileName);

		this->open_seq(fobj);

#ifdef USE_MPI
	}
#endif
}

// normal test cases
TYPED_TEST_P(FileLoadProcedureTest, open_posix)
{
#ifdef USE_MPI
  ::mxx::comm comm;
  if (comm.rank() == 0) {
#endif

    bliss::io::posix_file fobj(this->fileName);

    this->open_seq(fobj);

#ifdef USE_MPI
  }
#endif
}



#ifdef USE_MPI
TYPED_TEST_P(FileLoadProcedureTest, open_mmap_mpi)
{
//  ::mxx::comm comm;
    constexpr size_t overlap = TypeParam::get_overlap_size();

  ::bliss::io::parallel::partitioned_file<::bliss::io::mmap_file, ParserType> fobj(this->fileName, overlap, MPI_COMM_WORLD);

  this->open_mpi(fobj, MPI_COMM_WORLD);

  MPI_Barrier(MPI_COMM_WORLD);
//  comm.barrier();
}

TYPED_TEST_P(FileLoadProcedureTest, open_stdio_mpi)
{
  ::mxx::comm comm;
    constexpr size_t overlap = TypeParam::get_overlap_size();

  ::bliss::io::parallel::partitioned_file<::bliss::io::stdio_file, ParserType> fobj(this->fileName, overlap, comm);

  this->open_mpi(fobj, comm);

  comm.barrier();
}


TYPED_TEST_P(FileLoadProcedureTest, open_posix_mpi)
{
  ::mxx::comm comm;
    constexpr size_t overlap = TypeParam::get_overlap_size();

  ::bliss::io::parallel::partitioned_file<::bliss::io::posix_file, ParserType> fobj(this->fileName, overlap, comm);

  this->open_mpi(fobj, comm);

  comm.barrier();
}

TYPED_TEST_P(FileLoadProcedureTest, open_mmap_shared_mpi)
{
//	::mxx::comm comm;
	  constexpr size_t overlap = TypeParam::get_overlap_size();

	::bliss::io::parallel::partitioned_file<::bliss::io::mmap_file, ParserType, ::bliss::io::parallel::base_shared_fd_file> fobj(this->fileName, overlap, MPI_COMM_WORLD);

	this->open_mpi(fobj, MPI_COMM_WORLD);

	MPI_Barrier(MPI_COMM_WORLD);
//	comm.barrier();
}


TYPED_TEST_P(FileLoadProcedureTest, open_posix_shared_mpi)
{
  ::mxx::comm comm;
    constexpr size_t overlap = TypeParam::get_overlap_size();

  ::bliss::io::parallel::partitioned_file<::bliss::io::posix_file, ParserType, ::bliss::io::parallel::base_shared_fd_file> fobj(this->fileName, overlap, comm);

  this->open_mpi(fobj, comm);

  comm.barrier();
}


TYPED_TEST_P(FileLoadProcedureTest, open_mpiio_mpi)
{
	::mxx::comm comm;
	  constexpr size_t overlap = TypeParam::get_overlap_size();

	::bliss::io::parallel::mpiio_file<ParserType> fobj(this->fileName, overlap, comm);

	this->open_mpi(fobj, comm);

	comm.barrier();
}
#endif




// now register the test cases
REGISTER_TYPED_TEST_CASE_P(FileLoadProcedureTest,
#ifdef USE_MPI
		open_mmap_mpi,
		open_stdio_mpi,
    open_posix_mpi,
    open_mmap_shared_mpi,
    open_posix_shared_mpi,
		open_mpiio_mpi,
#endif
 open_mmap,
 open_stdio,
 open_posix
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
  ::mxx::env e(argc, argv);
  ::mxx::comm comm;

#endif

  result = RUN_ALL_TESTS();

#if defined(USE_MPI)
  comm.barrier();
#endif

  return result;
}



