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

#include "io/fastq_loader.hpp"
#include "io/file_loader.hpp"
#include <mxx/reduction.hpp>

#include "io/test/file_loader_test_fixtures.hpp"

using namespace bliss::io;


typedef FASTQLoader<unsigned char> FASTQLoaderType;


class FASTQParserTest : public FileParserTest<FASTQLoaderType >
{};



TEST_P(FASTQParserTest, parse)
{
	  this->elemCount = 0;

#ifdef USE_MPI
	int rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0) {
#endif

  // from FileLoader type, get the block iter type and range type
  using BlockIterType = typename FASTQLoaderType::L1BlockType::iterator;
  using SeqIterType = ::bliss::io::SequencesIterator<BlockIterType, bliss::io::FASTQParser >;

  FASTQLoaderType loader(this->fileName, 1, 0, 1, 2048, 2048 );

  auto l1 = loader.getNextL1Block();

  while (l1.getRange().size() > 0) {

    //==  and wrap the chunk inside an iterator that emits Reads.
    SeqIterType seqs_start(bliss::io::FASTQParser<BlockIterType>(), l1.begin(), l1.end(), l1.getRange().start);
    SeqIterType seqs_end(l1.end());

    //== loop over the reads
    for (; seqs_start != seqs_end; ++seqs_start, ++(this->elemCount))
    {
      //      std::cout << "Rank " << rank << "/" << nprocs << " seq: " << (*seqs_start).id.id << ", ";
      //      std::cout << std::distance(l1.begin(), (*seqs_start).seq_begin) << ", " << std::distance(l1.begin(), (*seqs_start).seq_end) << ", ";
      //      std::cout << std::distance(l1.begin(), (*seqs_start).qual_begin) << ", " << std::distance(l1.begin(), (*seqs_start).qual_end) << std::endl;
    }

    l1 = loader.getNextL1Block();
  }

#ifdef USE_MPI
	}
#endif

}


#ifdef USE_OPENMP
TEST_P(FASTQParserTest, parse_omp)
{
	  this->elemCount = 0;
#ifdef USE_MPI
	int rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0) {
#endif

  using BlockIterType = typename FASTQLoaderType::L1BlockType::iterator;
  using SeqIterType = ::bliss::io::SequencesIterator<BlockIterType, bliss::io::FASTQParser >;


  int nthreads = 4;

  FASTQLoaderType loader(this->fileName, 1, 0, nthreads, 2048, 2048 * nthreads * 2 );

  auto l1 = loader.getNextL1Block();
  size_t localKmerCount = 0;

  while (l1.getRange().size() > 0) {

    localKmerCount = 0;
#pragma omp parallel num_threads(nthreads) shared(loader) reduction(+:localKmerCount)
   {

      int tid = omp_get_thread_num();

      auto l2 = loader.getNextL2Block(tid);

      while (l2.getRange().size() > 0) {

        // from FileLoader type, get the block iter type and range type


        //==  and wrap the chunk inside an iterator that emits Reads.
        SeqIterType seqs_start(bliss::io::FASTQParser<BlockIterType>(), l2.begin(), l2.end(), l2.getRange().start);
        SeqIterType seqs_end(l2.end());

        //== loop over the reads
        for (; seqs_start != seqs_end; ++seqs_start, ++localKmerCount)
        {
          //      std::cout << "Rank " << rank << "/" << nprocs << " seq: " << (*seqs_start).id.id << ", ";
          //      std::cout << std::distance(l1.begin(), (*seqs_start).seq_begin) << ", " << std::distance(l1.begin(), (*seqs_start).seq_end) << ", ";
          //      std::cout << std::distance(l1.begin(), (*seqs_start).qual_begin) << ", " << std::distance(l1.begin(), (*seqs_start).qual_end) << std::endl;
        }

        l2 = loader.getNextL2Block(tid);
      }
    }  // end omp parallel region.

   this->elemCount += localKmerCount;

   l1 = loader.getNextL1Block();
  }  // end l1 while.

#ifdef USE_MPI
	}
#endif
}
#endif



#ifdef USE_MPI
TEST_P(FASTQParserTest, parse_mpi)
{
  this->elemCount = 0;

  // from FileLoader type, get the block iter type and range type
  using BlockIterType = typename FASTQLoaderType::L1BlockType::iterator;
  using SeqIterType = ::bliss::io::SequencesIterator<BlockIterType, bliss::io::FASTQParser >;

  int nthreads = 1;

  FASTQLoaderType loader(this->fileName, MPI_COMM_WORLD, nthreads, 2048, 2048 * nthreads * 2);

  int rank, nprocs;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  auto l1 = loader.getNextL1Block();

  while (l1.getRange().size() > 0) {


    //==  and wrap the chunk inside an iterator that emits Reads.
    SeqIterType seqs_start(bliss::io::FASTQParser<BlockIterType>(), l1.begin(), l1.end(), l1.getRange().start);
    SeqIterType seqs_end(l1.end());

    //== loop over the reads
    for (; seqs_start != seqs_end; ++seqs_start, ++this->elemCount)
    {
//
//
//            std::cout << "Rank " << rank << "/" << nprocs << " seq: " << (*seqs_start).id.get_id() << ", offset " << (*seqs_start).seq_offset << " size " << (*seqs_start).record_size << " ";
//            std::cout << std::distance(l1.begin(), (*seqs_start).seq_begin) << ", " << std::distance(l1.begin(), (*seqs_start).seq_end) << ", ";
//            std::cout << std::distance(l1.begin(), (*seqs_start).qual_begin) << ", " << std::distance(l1.begin(), (*seqs_start).qual_end) << std::endl;
    }

    l1 = loader.getNextL1Block();
  }
}
#endif


#ifdef USE_MPI
#ifdef USE_OPENMP
TEST_P(FASTQParserTest, parse_mpi_omp)
{
  this->elemCount = 0;

  using BlockIterType = typename FASTQLoaderType::L1BlockType::iterator;
  using SeqIterType = ::bliss::io::SequencesIterator<BlockIterType, bliss::io::FASTQParser >;


  int nthreads = 4;

  FASTQLoaderType loader(this->fileName, MPI_COMM_WORLD, nthreads, 2048, 2048 * nthreads * 2);

  auto l1 = loader.getNextL1Block();
  size_t localKmerCount = 0;

  while (l1.getRange().size() > 0) {

    localKmerCount = 0;

    //std::cout << "L1 block " << l1.getRange() << std::endl;

#pragma omp parallel num_threads(nthreads) shared(loader) reduction(+:localKmerCount)
   {

      int tid = omp_get_thread_num();

      auto l2 = loader.getNextL2Block(tid);

      while (l2.getRange().size() > 0) {

        // from FileLoader type, get the block iter type and range type
	//	std::cout << "thread " << tid << " L1 block " << l1.getRange() << " L2 block " << l2.getRange() << std::endl;


        //==  and wrap the chunk inside an iterator that emits Reads.
        SeqIterType seqs_start(bliss::io::FASTQParser<BlockIterType>(), l2.begin(), l2.end(), l2.getRange().start);
        SeqIterType seqs_end(l2.end());

        //== loop over the reads
        for (; seqs_start != seqs_end; ++seqs_start, ++localKmerCount)
        {
          //      std::cout << "Rank " << rank << "/" << nprocs << " seq: " << (*seqs_start).id.id << ", ";
          //      std::cout << std::distance(l1.begin(), (*seqs_start).seq_begin) << ", " << std::distance(l1.begin(), (*seqs_start).seq_end) << ", ";
          //      std::cout << std::distance(l1.begin(), (*seqs_start).qual_begin) << ", " << std::distance(l1.begin(), (*seqs_start).qual_end) << std::endl;
        }

        l2 = loader.getNextL2Block(tid);
      }
    }  // end omp parallel region.

   this->elemCount += localKmerCount;

   l1 = loader.getNextL1Block();
  }  // end l1 while.

}
#endif
#endif


INSTANTIATE_TEST_CASE_P(Bliss, FASTQParserTest, ::testing::Values(
    TestFileInfo(243, 27580, std::string("/test/data/natural.fastq")),
    TestFileInfo(250, 29250, std::string("/test/data/natural.withN.fastq")),
    TestFileInfo(7, 939, std::string("/test/data/test.debruijn.small.fastq")),
    TestFileInfo(1, 134, std::string("/test/data/test.debruijn.tiny.fastq")),
    TestFileInfo(254562, 34111308, std::string("/test/data/test.fastq")),
    TestFileInfo(140, 18761, std::string("/test/data/test.medium.fastq")),
    TestFileInfo(7, 939, std::string("/test/data/test.small.fastq")),
    TestFileInfo(1, 33194, std::string("/test/data/test.unitiq1.fastq")),
    TestFileInfo(1, 7144, std::string("/test/data/test.unitiq1.short2.fastq")),
    TestFileInfo(1, 14824, std::string("/test/data/test.unitiq1.short.fastq")),
    TestFileInfo(1, 29296, std::string("/test/data/test.unitiq2.fastq")),
    TestFileInfo(2, 62490, std::string("/test/data/test.unitiqs.fastq"))
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



