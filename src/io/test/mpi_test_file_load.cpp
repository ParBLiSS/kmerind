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
#include "mxx/env.hpp"
#endif

// include google test
#include <gtest/gtest.h>
#include <cstdint> // for uint64_t, etc.
#include <string>
#include <utility>  // for pair
#include <type_traits>  // for integral_constant


#include "io/test/file_load_test_fixtures.hpp"

#include "io/file_loader.hpp"
#include "io/fastq_loader.hpp"
#include "io/fasta_loader.hpp"
#include "io/file.hpp"



//TODO: test file load direct, serially.
// - just with base file types.  no parsers needed
//verify that file is read correctly
//
//TODO: test file load direct, MPI.
//test with 3 parsers.  test with different overlaps for basefileparser and fastaparser
//first verify that overlaps are correct (0 for fastq parser.)
//next verify that there is complete coverage
//finally verify that file is read correctly.
// for fastq parser, verify that first char is @
//
//NEXT: test file load and parse to sequence  (fasta and fastq only.)  pick POSIX for this test.  k = 31.
//FINALLY:  test k-mer parse (fasta and fastq only).  pick POSIX for this test.  k = 31



template <typename file_loader>
class FileSequentialLoadTest : public FileLoaderTest
{

  protected:
    using ValueType = typename FileLoaderTest::ValueType;
    using InputIterType = typename FileLoaderTest::InputIterType;


	~FileSequentialLoadTest() {};

	virtual void open(file_loader & fobj) {

		// get fileName
		bliss::io::file_data fdata = fobj.read_file();

		ASSERT_TRUE(fobj.size() > 0);

		ASSERT_EQ(fdata.getRange().size(), fobj.size());
		ASSERT_EQ(fdata.getRange().end, fobj.size());

		ASSERT_EQ(fdata.in_mem_range_bytes.size(), fobj.size());
		ASSERT_EQ(fdata.in_mem_range_bytes.end, fobj.size());

		ASSERT_EQ(fdata.parent_range_bytes.size(), fobj.size());
		ASSERT_EQ(fdata.parent_range_bytes.end, fobj.size());

		bool same = true;

		ValueType * data = new ValueType[fobj.size()];
		this->readFilePOSIX(this->fileName, 0, fobj.size(), data);

		same = equal(data, fdata.begin(), fobj.size(), true);

		delete [] data;
		ASSERT_TRUE(same);

	}
};



// indicate this is a typed test
TYPED_TEST_CASE_P(FileSequentialLoadTest);


// normal test cases
TYPED_TEST_P(FileSequentialLoadTest, read)
{
#ifdef USE_MPI
	  ::mxx::comm comm;
	if (comm.rank() == 0) {
#endif

		TypeParam fobj(this->fileName);

		this->open(fobj);

#ifdef USE_MPI
	}
#endif
}



// now register the test cases
REGISTER_TYPED_TEST_CASE_P(FileSequentialLoadTest,
 read
);



typedef ::testing::Types<
		bliss::io::mmap_file,
		bliss::io::stdio_file,
		bliss::io::posix_file
> FileSequentialLoadTestTypes;

//typedef ::testing::Types< FileLoader<unsigned char, 0, bliss::io::BaseFileParser, false, false> > FileSequentialLoadTestTypes;

INSTANTIATE_TYPED_TEST_CASE_P(Bliss, FileSequentialLoadTest, FileSequentialLoadTestTypes);


#ifdef USE_MPI

template <typename file_loader>
class FileMPILoadTest : public FileLoaderTest
{
protected:
	~FileMPILoadTest() {};

    virtual void SetUp()
    {
      this->fileName.assign(PROJ_SRC_DIR);
      this->fileName.append("/test/data/test.medium.fasta");  // does not matter which file to use.

      // get file size
      struct stat filestat;
      stat(this->fileName.c_str(), &filestat);
      size_t fileSize = static_cast<size_t>(filestat.st_size);

      ASSERT_EQ(1092580UL, fileSize);
    }

	virtual void open(typename std::tuple_element<0, file_loader>::type & fobj, size_t const & overlap, mxx::comm const & comm) {

		::bliss::io::file_data fdata = fobj.read_file();

		ASSERT_TRUE(fobj.size() > 0);

		// parent range should match.
		ASSERT_EQ(fdata.parent_range_bytes.size(), fobj.size());
		ASSERT_EQ(fdata.parent_range_bytes.end, fobj.size());

		// get the ranges to makes sure they are in memory.
		ASSERT_TRUE(fdata.in_mem_range_bytes.start >= fdata.parent_range_bytes.start);
		ASSERT_TRUE(fdata.in_mem_range_bytes.end <= fdata.parent_range_bytes.end);

		ASSERT_TRUE(fdata.in_mem_range_bytes.start <= fdata.valid_range_bytes.start);
		ASSERT_TRUE(fdata.in_mem_range_bytes.end >= fdata.valid_range_bytes.end);


		// get all the sizes to make sure total is same as file size.
		size_t region_size = fdata.valid_range_bytes.size() - (comm.rank() == (comm.size() - 1) ? 0 : overlap);
		region_size = ::mxx::allreduce(region_size, comm);
		ASSERT_EQ(region_size, fobj.size());

		// make sure the aggregate range is same as parent range.
		std::vector<size_t> begins = mxx::allgather(fdata.valid_range_bytes.start, comm);
		std::vector<size_t> ends = mxx::allgather(fdata.valid_range_bytes.end, comm);

		ASSERT_EQ(begins.front(), 0UL);
		ASSERT_EQ(ends.back(), fobj.size());

		for (int i = 1; i < comm.size(); ++i) {
			ASSERT_TRUE(ends[i-1] == (begins[i] + overlap) );
		}

		// make sure the files read are the same.
		bool same = true;

		ValueType * data = new ValueType[fdata.valid_range_bytes.size()];
		this->readFilePOSIX(this->fileName,
				fdata.valid_range_bytes.start,
				fdata.valid_range_bytes.size(), data);

		same = equal(data, fdata.begin(), fdata.valid_range_bytes.size(), true);

		delete [] data;
		ASSERT_TRUE(same);

	}
};


// indicate this is a typed test
TYPED_TEST_CASE_P(FileMPILoadTest);

TYPED_TEST_P(FileMPILoadTest, read)
{
	  ::mxx::comm comm;

	  constexpr size_t overlap = std::tuple_element<1, TypeParam>::type::value;

	  typename std::tuple_element<0, TypeParam>::type fobj(this->fileName, overlap, comm);

	  this->open(fobj, overlap, comm);

	  comm.barrier();
}


// now register the test cases
REGISTER_TYPED_TEST_CASE_P(FileMPILoadTest,
 read
);


typedef ::testing::Types<
	std::pair<::bliss::io::parallel::partitioned_file<::bliss::io::mmap_file,  ::bliss::io::BaseFileParser<typename ::bliss::io::file_data::const_iterator> >, std::integral_constant<size_t, 0> >,
	std::pair<::bliss::io::parallel::partitioned_file<::bliss::io::stdio_file, ::bliss::io::BaseFileParser<typename ::bliss::io::file_data::const_iterator> >, std::integral_constant<size_t, 0> >,
	std::pair<::bliss::io::parallel::partitioned_file<::bliss::io::posix_file, ::bliss::io::BaseFileParser<typename ::bliss::io::file_data::const_iterator> >, std::integral_constant<size_t, 0> >,
	std::pair<::bliss::io::parallel::partitioned_file<::bliss::io::mmap_file,  ::bliss::io::BaseFileParser<typename ::bliss::io::file_data::const_iterator> , ::bliss::io::parallel::base_shared_fd_file>, std::integral_constant<size_t, 0> >,
	std::pair<::bliss::io::parallel::partitioned_file<::bliss::io::posix_file, ::bliss::io::BaseFileParser<typename ::bliss::io::file_data::const_iterator> , ::bliss::io::parallel::base_shared_fd_file>, std::integral_constant<size_t, 0> >,
	std::pair<::bliss::io::parallel::mpiio_file<::bliss::io::BaseFileParser<typename ::bliss::io::file_data::const_iterator> >, std::integral_constant<size_t, 0> >,
	std::pair<::bliss::io::parallel::partitioned_file<::bliss::io::mmap_file,  ::bliss::io::FASTAParser<typename ::bliss::io::file_data::const_iterator> >, std::integral_constant<size_t, 0> >,
	std::pair<::bliss::io::parallel::partitioned_file<::bliss::io::stdio_file, ::bliss::io::FASTAParser<typename ::bliss::io::file_data::const_iterator> >, std::integral_constant<size_t, 0> >,
	std::pair<::bliss::io::parallel::partitioned_file<::bliss::io::posix_file, ::bliss::io::FASTAParser<typename ::bliss::io::file_data::const_iterator> >, std::integral_constant<size_t, 0> >,
	std::pair<::bliss::io::parallel::partitioned_file<::bliss::io::mmap_file,  ::bliss::io::FASTAParser<typename ::bliss::io::file_data::const_iterator> , ::bliss::io::parallel::base_shared_fd_file>,  std::integral_constant<size_t, 0> >,
	std::pair<::bliss::io::parallel::partitioned_file<::bliss::io::posix_file, ::bliss::io::FASTAParser<typename ::bliss::io::file_data::const_iterator> , ::bliss::io::parallel::base_shared_fd_file>,  std::integral_constant<size_t, 0> >,
	std::pair<::bliss::io::parallel::mpiio_file<::bliss::io::FASTAParser<typename ::bliss::io::file_data::const_iterator> >, std::integral_constant<size_t, 0> >
> FileMPILoadTestTypes;

INSTANTIATE_TYPED_TEST_CASE_P(Bliss, FileMPILoadTest, FileMPILoadTestTypes);



template <typename file_loader>
class FASTQMPILoadTest : public FileLoaderTest
{
protected:
	~FASTQMPILoadTest() {};

	virtual void open(typename std::tuple_element<0, file_loader>::type & fobj, size_t const & overlap, mxx::comm const & comm) {

		::bliss::io::file_data fdata = fobj.read_file();

		ASSERT_TRUE(fobj.size() > 0);

		// parent range should match.
		ASSERT_EQ(fdata.parent_range_bytes.size(), fobj.size());
		ASSERT_EQ(fdata.parent_range_bytes.end, fobj.size());

		// get the ranges to makes sure they are in memory.
		ASSERT_TRUE(fdata.in_mem_range_bytes.start >= fdata.parent_range_bytes.start);
		ASSERT_TRUE(fdata.in_mem_range_bytes.end <= fdata.parent_range_bytes.end);

		ASSERT_TRUE(fdata.in_mem_range_bytes.start <= fdata.valid_range_bytes.start);
		ASSERT_TRUE(fdata.in_mem_range_bytes.end >= fdata.valid_range_bytes.end);


		// get all the sizes to make sure total is same as file size.
		size_t region_size = fdata.valid_range_bytes.size();
		region_size = ::mxx::allreduce(region_size, comm);
		ASSERT_EQ(region_size, fobj.size());

		// make sure the aggregate range is same as parent range.
		std::vector<size_t> begins = mxx::allgather(fdata.valid_range_bytes.start, comm);
		std::vector<size_t> ends = mxx::allgather(fdata.valid_range_bytes.end, comm);

		ASSERT_EQ(begins.front(), 0UL);
		ASSERT_EQ(ends.back(), fobj.size());

		for (int i = 1; i < comm.size(); ++i) {
			ASSERT_TRUE(ends[i-1] == begins[i] );
		}

		// make sure the files read are the same.

		if (fdata.valid_range_bytes.size() > 0) {
			bool same = true;

			ValueType * data = new ValueType[fdata.valid_range_bytes.size()];
			this->readFilePOSIX(this->fileName,
					fdata.valid_range_bytes.start,
					fdata.valid_range_bytes.size(), data);

			same = equal(data, fdata.begin(), fdata.valid_range_bytes.size(), true);
			ASSERT_TRUE(same);

			ASSERT_EQ(data[0], '@');
			ASSERT_EQ(*(fdata.begin()), '@');

			delete [] data;
		}
	}
};

// indicate this is a typed test
TYPED_TEST_CASE_P(FASTQMPILoadTest);

TYPED_TEST_P(FASTQMPILoadTest, read)
{
	  ::mxx::comm comm;

	  constexpr size_t overlap = std::tuple_element<1, TypeParam>::type::value;
  typename std::tuple_element<0, TypeParam>::type fobj(this->fileName, overlap, comm);

  this->open(fobj, overlap, comm);

  comm.barrier();
}


// now register the test cases
REGISTER_TYPED_TEST_CASE_P(FASTQMPILoadTest,
 read
);


typedef ::testing::Types<
	std::pair<::bliss::io::parallel::partitioned_file<::bliss::io::mmap_file,  ::bliss::io::FASTQParser<typename ::bliss::io::file_data::const_iterator> >, std::integral_constant<size_t, 0> >,
	std::pair<::bliss::io::parallel::partitioned_file<::bliss::io::stdio_file, ::bliss::io::FASTQParser<typename ::bliss::io::file_data::const_iterator> >, std::integral_constant<size_t, 0> >,
	std::pair<::bliss::io::parallel::partitioned_file<::bliss::io::posix_file, ::bliss::io::FASTQParser<typename ::bliss::io::file_data::const_iterator> >, std::integral_constant<size_t, 0> >,
	std::pair<::bliss::io::parallel::partitioned_file<::bliss::io::mmap_file,  ::bliss::io::FASTQParser<typename ::bliss::io::file_data::const_iterator> , ::bliss::io::parallel::base_shared_fd_file>, std::integral_constant<size_t, 0> >,
	std::pair<::bliss::io::parallel::partitioned_file<::bliss::io::posix_file, ::bliss::io::FASTQParser<typename ::bliss::io::file_data::const_iterator> , ::bliss::io::parallel::base_shared_fd_file>, std::integral_constant<size_t, 0> >,
	std::pair<::bliss::io::parallel::mpiio_file<::bliss::io::FASTQParser<typename ::bliss::io::file_data::const_iterator> >, std::integral_constant<size_t, 0> >
> FASTQMPILoadTestTypes;

INSTANTIATE_TYPED_TEST_CASE_P(Bliss, FASTQMPILoadTest, FASTQMPILoadTestTypes);



#endif



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



