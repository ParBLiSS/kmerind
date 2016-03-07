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
 * fastaloader_test.cpp
 *
 *  Created on: Oct 30, 2014
 *      Author: cjain
 */

#include "bliss-config.hpp"    // for location of data.

#if defined(USE_MPI)
#include "mpi.h"
#endif

// include google test
#include <gtest/gtest.h>
#include <cstdint> // for uint64_t, etc.
#include <string>
#include <limits>

#include "common/sequence.hpp"
#include "io/sequence_iterator.hpp"

#include "io/file.hpp"
#include "io/fasta_loader.hpp"
#include "io/file_loader.hpp"
#include <mxx/reduction.hpp>

#include "io/test/file_loader_test_fixtures.hpp"

using namespace bliss::io;

static constexpr size_t block_size = 32768;

// for these tests, we are counting sequences, not kmers, so no overlap.
typedef FASTALoader<unsigned char, 0> FileLoaderType;
typedef FASTAParser<unsigned char *> ParserType;


class FASTAParseProcedureTest : public FileParserTest<FileLoaderType >
{
	// can't compare the sequences directly since offsets may be different.
protected:

	template <typename file_loader>
	void parse_seq(file_loader & fobj ) {
		  typedef typename FileLoaderType::RangeType RangeType;

		  // get fileName

		  FileLoaderType loader(this->fileName, 1, 0, 1, 64 * 1024, 64 * 1024);   // NOTE: larger block is MUCH faster.  2 * threads * 2048 bytes take a LONG time because of mmap each time...
		  RangeType fr = loader.getFileRange();

		  auto block = loader.getNextL1Block();
		  RangeType r = block.getRange();

		  // start at file starting point
		  ASSERT_EQ(fr.start, r.start);
		  ASSERT_EQ(fr.end, r.end);


			bliss::io::file_data fdata = fobj.read_file();

		  ASSERT_TRUE(fdata.in_mem_range_bytes.size() == r.size());

		  bool same = std::equal(block.begin(), block.end(), fdata.begin());

		    ASSERT_TRUE(same);


			this->elemCount = 0;

			// from FileLoader type, get the block iter type and range type
			using BlockIterType = unsigned char *;
			using SeqIterType = ::bliss::io::SequencesIterator<BlockIterType, bliss::io::FASTAParser >;

			ParserType l1parser = loader.getSeqParser();


			size_t offset = l1parser.init_parser(fdata.data.data(), loader.getFileRange(), fdata.in_mem_range_bytes, fdata.valid_range_bytes);

			  //==  and wrap the chunk inside an iterator that emits Reads.
			  SeqIterType seqs_start(l1parser, fdata.data.data() + (fdata.valid_range_bytes.start - fdata.in_mem_range_bytes.start),
					  fdata.data.data() + (fdata.valid_range_bytes.end - fdata.in_mem_range_bytes.start), offset);
			  SeqIterType seqs_end(fdata.data.data() + (fdata.valid_range_bytes.end - fdata.in_mem_range_bytes.start));

			  //== loop over the reads
			  for (; seqs_start != seqs_end; ++seqs_start)
			  {
				if (fdata.valid_range_bytes.start <= (*seqs_start).id.get_pos()) ++(this->elemCount);
		  //            std::cout << *seqs_start << ", ";
		  //            std::cout << std::distance((*seqs_start).seq_begin, (*seqs_start).seq_end) << std::endl;
				//      std::cout << std::distance(l1.begin(), (*seqs_start).qual_begin) << ", " << std::distance(l1.begin(), (*seqs_start).qual_end) << std::endl;
			  }



	}


	template <typename file_loader>
	void parse_mpi(file_loader & fobj, MPI_Comm comm) {
		  typedef typename FileLoaderType::RangeType RangeType;

		  // get this->fileName

		  FileLoaderType loader(this->fileName, comm, 1, 64 * 1024, 64 * 1024 );

		  auto block = loader.getNextL1Block();
		  RangeType r = block.getRange();
		  bool same = true;

			::bliss::io::file_data fdata = fobj.read_file();

	  if (fdata.valid_range_bytes.size() != r.size()) {
		  std::cout << "file open in mem: " << fdata.in_mem_range_bytes << ::std::endl;
		  std::cout << "file open valid: " << fdata.valid_range_bytes << ::std::endl;
		  std::cout << "file loader: " << r << ::std::endl;
	  }


	//	  std::cout << " orig mem " << fdata.in_mem_range_bytes << " valid " << fdata.valid_range_bytes << std::endl;
		  ASSERT_TRUE(fdata.valid_range_bytes.size() == r.size());

		  same = std::equal(block.begin(), block.end(), fdata.begin());

		  ASSERT_TRUE(same);

		  this->elemCount = 0;

		// from FileLoader type, get the block iter type and range type
		using BlockIterType = unsigned char *;
		using SeqIterType = ::bliss::io::SequencesIterator<BlockIterType, bliss::io::FASTAParser >;

		ParserType l1parser = loader.getSeqParser();

		size_t offset = l1parser.init_parser(fdata.data.data(), loader.getFileRange(), fdata.in_mem_range_bytes, fdata.valid_range_bytes, comm);

		  //==  and wrap the chunk inside an iterator that emits Reads.
		  SeqIterType seqs_start(l1parser, fdata.data.data() + (fdata.valid_range_bytes.start - fdata.in_mem_range_bytes.start),
				  fdata.data.data() + (fdata.valid_range_bytes.end - fdata.in_mem_range_bytes.start), offset);
		  SeqIterType seqs_end(fdata.data.data() + (fdata.valid_range_bytes.end - fdata.in_mem_range_bytes.start));

		  //== loop over the reads
		  for (; seqs_start != seqs_end; ++seqs_start)
		  {
			if (fdata.valid_range_bytes.start <= (*seqs_start).id.get_pos()) ++(this->elemCount);
	  //            std::cout << *seqs_start << ", ";
	  //            std::cout << std::distance((*seqs_start).seq_begin, (*seqs_start).seq_end) << std::endl;
			//      std::cout << std::distance(l1.begin(), (*seqs_start).qual_begin) << ", " << std::distance(l1.begin(), (*seqs_start).qual_end) << std::endl;
		  }
	//	  int rank;
	//	  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	//	  std::cout << "rank " << rank << " elem count " << this->elemCount << std::endl;

	//	  if (this->elemCount == 0) {
	//		  std::cout << " mem " << fdata.in_mem_range_bytes << " valid " << fdata.valid_range_bytes << std::endl;
	//	  }


	}

};


TEST_P(FASTAParseProcedureTest, parse_mmap)
{
#ifdef USE_MPI
	int rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0) {
#endif


		bliss::io::mmap_file fobj(this->fileName);

		this->parse_seq(fobj);

#ifdef USE_MPI
	}
#endif

}

TEST_P(FASTAParseProcedureTest, parse_stdio)
{
#ifdef USE_MPI
	int rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0) {
#endif


		bliss::io::stdio_file fobj(this->fileName);

		this->parse_seq(fobj);

#ifdef USE_MPI
	}
#endif

}



#ifdef USE_MPI
TEST_P(FASTAParseProcedureTest, parse_mmap_mpi)
{
	  constexpr size_t overlap = FileLoaderType::get_overlap_size();

	::bliss::io::parallel::partitioned_file<::bliss::io::mmap_file, ParserType> fobj(this->fileName, overlap, MPI_COMM_WORLD);

	this->parse_mpi(fobj, MPI_COMM_WORLD);

	  MPI_Barrier(MPI_COMM_WORLD);
}
#endif



#ifdef USE_MPI
TEST_P(FASTAParseProcedureTest, parse_stdio_mpi)
{
	  constexpr size_t overlap = FileLoaderType::get_overlap_size();

	::bliss::io::parallel::partitioned_file<::bliss::io::stdio_file, ParserType> fobj(this->fileName, overlap, MPI_COMM_WORLD);

	this->parse_mpi(fobj, MPI_COMM_WORLD);

	  MPI_Barrier(MPI_COMM_WORLD);
}
#endif



#ifdef USE_MPI
TEST_P(FASTAParseProcedureTest, parse_mpiio_mpi)
{
	  constexpr size_t overlap = FileLoaderType::get_overlap_size();

	::bliss::io::parallel::mpiio_file<ParserType> fobj(this->fileName, overlap, MPI_COMM_WORLD);

	this->parse_mpi(fobj, MPI_COMM_WORLD);

	  MPI_Barrier(MPI_COMM_WORLD);
}
#endif





INSTANTIATE_TEST_CASE_P(Bliss, FASTAParseProcedureTest, ::testing::Values(
    TestFileInfo(243, 13790, std::string("/test/data/natural.fasta")),
    TestFileInfo(250, 14625, std::string("/test/data/natural.withN.fasta")),
    TestFileInfo(4, 512, std::string("/test/data/test.fasta")),
    TestFileInfo(5000, 1092580, std::string("/test/data/test.medium.fasta")),
    TestFileInfo(6, 940, std::string("/test/data/test2.fasta"))
));

//INSTANTIATE_TEST_CASE_P(Bliss, FASTAParseProcedureTest, ::testing::Values(
//		TestFileInfo(250, 14625, std::string("/test/data/natural.fasta"))
//		));

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

