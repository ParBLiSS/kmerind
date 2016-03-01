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
 * fastqloader_test.cpp
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
#include <limits>

#include "common/sequence.hpp"
#include "io/sequence_iterator.hpp"

#include "io/file.hpp"
#include "io/fastq_loader.hpp"
#include "io/file_loader.hpp"
#include <mxx/reduction.hpp>

#include "io/test/file_loader_test_fixtures.hpp"

using namespace bliss::io;


typedef FASTQLoader<unsigned char> FileLoaderType;
typedef FASTQParser<unsigned char *> ParserType;


class FASTQParseProcedureTest : public FileParserTest<FileLoaderType >
{};


TEST_P(FASTQParseProcedureTest, parse)
{
#ifdef USE_MPI
	int rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0) {
#endif


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


		  ::bliss::io::MemData data = ::bliss::io::memmap::load_file<ParserType >(this->fileName, overlap);

		  ASSERT_TRUE(data.data[data.valid_range.start - data.mem_range.start] == '@');

			this->elemCount = 0;

			// from FileLoader type, get the block iter type and range type
			using BlockIterType = unsigned char *;
			using SeqIterType = ::bliss::io::SequencesIterator<BlockIterType, bliss::io::FASTQParser >;

			ParserType l1parser = loader.getSeqParser();


			size_t offset = l1parser.init_parser(data.data, loader.getFileRange(), data.mem_range, data.valid_range);

			  //==  and wrap the chunk inside an iterator that emits Reads.
			  SeqIterType seqs_start(l1parser, data.data + (data.valid_range.start - data.mem_range.start),
					  data.data + (data.valid_range.end - data.mem_range.start), offset);
			  SeqIterType seqs_end(data.data + (data.valid_range.end - data.mem_range.start));

			  //== loop over the reads
			  for (; seqs_start != seqs_end; ++seqs_start)
			  {
				if (data.valid_range.start <= (*seqs_start).id.get_pos()) ++(this->elemCount);
		  //            std::cout << *seqs_start << ", ";
		  //            std::cout << std::distance((*seqs_start).seq_begin, (*seqs_start).seq_end) << std::endl;
				//      std::cout << std::distance(l1.begin(), (*seqs_start).qual_begin) << ", " << std::distance(l1.begin(), (*seqs_start).qual_end) << std::endl;
			  }


		  ::bliss::io::unload_data(data);



#ifdef USE_MPI
	}
#endif

}




#ifdef USE_MPI
TEST_P(FASTQParseProcedureTest, parse_mpi)
{
	  constexpr size_t overlap = FileLoaderType::get_overlap_size();

	  // get this->fileName

	  FileLoaderType loader(this->fileName, MPI_COMM_WORLD, 1, 64 * 1024, 64 * 1024 );

	  auto block = loader.getNextL1Block();

	  ::bliss::io::MemData data = ::bliss::io::parallel::memmap::load_file<ParserType >(this->fileName, MPI_COMM_WORLD,  overlap);

	  if (data.valid_range.size() > 0) {
		  if (data.data[data.valid_range.start - data.mem_range.start] != '@') {
			  std::cout << "data at pos " << (data.valid_range.start - data.mem_range.start) << " for "
					  << data.valid_range << " in mem " << data.mem_range << " is "
					  << data.data[data.valid_range.start - data.mem_range.start] << std::endl;
		  }

		  ASSERT_TRUE(data.data[data.valid_range.start - data.mem_range.start] == '@');
	  }

//	  std::cout << " orig mem " << data.mem_range << " valid " << data.valid_range << std::endl;

	  this->elemCount = 0;

	// from FileLoader type, get the block iter type and range type
	using BlockIterType = unsigned char *;
	using SeqIterType = ::bliss::io::SequencesIterator<BlockIterType, bliss::io::FASTQParser >;

	ParserType l1parser = loader.getSeqParser();


	size_t offset = l1parser.init_parser(data.data, loader.getFileRange(), data.mem_range, data.valid_range, MPI_COMM_WORLD);

	  //==  and wrap the chunk inside an iterator that emits Reads.
	  SeqIterType seqs_start(l1parser, data.data + (data.valid_range.start - data.mem_range.start),
			  data.data + (data.valid_range.end - data.mem_range.start), offset);
	  SeqIterType seqs_end(data.data + (data.valid_range.end - data.mem_range.start));

	  //== loop over the reads
	  for (; seqs_start != seqs_end; ++seqs_start)
	  {
		if (data.valid_range.start <= (*seqs_start).id.get_pos()) ++(this->elemCount);
  //            std::cout << *seqs_start << ", ";
  //            std::cout << std::distance((*seqs_start).seq_begin, (*seqs_start).seq_end) << std::endl;
		//      std::cout << std::distance(l1.begin(), (*seqs_start).qual_begin) << ", " << std::distance(l1.begin(), (*seqs_start).qual_end) << std::endl;
	  }
//	  int rank;
//	  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//	  std::cout << "rank " << rank << " elem count " << this->elemCount << std::endl;

//	  if (this->elemCount == 0) {
//		  std::cout << " mem " << data.mem_range << " valid " << data.valid_range << std::endl;
//	  }

	  ::bliss::io::unload_data(data);

	  MPI_Barrier(MPI_COMM_WORLD);
}
#endif




INSTANTIATE_TEST_CASE_P(Bliss, FASTQParseProcedureTest, ::testing::Values(
    TestFileInfo(243, 27580, std::string("/test/data/natural.fastq")),
    TestFileInfo(250, 29250, std::string("/test/data/natural.withN.fastq")),
    TestFileInfo(140, 18761, std::string("/test/data/test.medium.fastq")),
    TestFileInfo(7, 939, std::string("/test/data/test.small.fastq")),
    TestFileInfo(254562, 34111308, std::string("/test/data/test.fastq"))
));

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



