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

#include "io/fasta_loader.hpp"
#include "io/file_loader.hpp"
#include <mxx/reduction.hpp>

#include "io/test/file_loader_test_fixtures.hpp"

using namespace bliss::io;

static constexpr size_t block_size = 32768;

// for these tests, we are counting sequences, not kmers, so no overlap.
typedef FASTALoader<unsigned char, 0> FASTALoaderType;


class FASTAParserTest : public FileParserTest<FASTALoaderType >
{};



TEST_P(FASTAParserTest, parse)
{
	this->elemCount = 0;


#ifdef USE_MPI
	int rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0) {
#endif

  // from FileLoader type, get the block iter type and range type
  using BlockIterType = typename FASTALoaderType::L1BlockType::iterator;
  using SeqIterType = ::bliss::io::SequencesIterator<BlockIterType, bliss::io::FASTAParser >;

  FASTALoaderType loader(this->fileName, 1, 0, 1, block_size, block_size );

  auto l1 = loader.getNextL1Block();

  bliss::io::FASTAParser<BlockIterType> l1parser = loader.getSeqParser();

  while (l1.getRange().size() > 0) {

    l1parser.init_for_iterator(l1.begin(), loader.getFileRange(), l1.getRange(), l1.getRange());

    //==  and wrap the chunk inside an iterator that emits Reads.
    SeqIterType seqs_start(l1parser, l1.begin(), l1.end(), l1.getRange().start);
    SeqIterType seqs_end(l1.end());

    //== loop over the reads
    for (; seqs_start != seqs_end; ++seqs_start)
    {
    	if (l1.getRange().start <= (*seqs_start).id.pos_in_file) ++(this->elemCount);
//            std::cout << *seqs_start << ", ";
//            std::cout << std::distance((*seqs_start).seqBegin, (*seqs_start).seqEnd) << std::endl;
      //      std::cout << std::distance(l1.begin(), (*seqs_start).qualBegin) << ", " << std::distance(l1.begin(), (*seqs_start).qualEnd) << std::endl;
    }

    l1 = loader.getNextL1Block();
  }
#ifdef USE_MPI
	}
#endif
}


#ifdef USE_OPENMP
TEST_P(FASTAParserTest, parse_omp)
{
  this->elemCount = 0;

#ifdef USE_MPI
	int rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0) {
#endif

  using BlockIterType = typename FASTALoaderType::L1BlockType::iterator;
  using SeqIterType = ::bliss::io::SequencesIterator<BlockIterType, bliss::io::FASTAParser >;


  int nthreads = 4;

  FASTALoaderType loader(this->fileName, 1, 0, nthreads, block_size, block_size * nthreads * 2 );

  auto l1 = loader.getNextL1Block();
  size_t localKmerCount = 0;

  bliss::io::FASTAParser<BlockIterType> l1parser = loader.getSeqParser();

  while (l1.getRange().size() > 0) {

    l1parser.init_for_iterator(l1.begin(), loader.getFileRange(), l1.getRange(), l1.getRange());

    localKmerCount = 0;
#pragma omp parallel num_threads(nthreads) shared(loader, l1parser) reduction(+:localKmerCount)
   {
    	  bliss::io::FASTAParser<typename FASTALoaderType::L2BlockType::iterator> l2parser = l1parser;

      int tid = omp_get_thread_num();

      auto l2 = loader.getNextL2Block(tid);


      while (l2.getRange().size() > 0) {

        // from FileLoader type, get the block iter type and range type


        //==  and wrap the chunk inside an iterator that emits Reads.
        SeqIterType seqs_start(l2parser, l2.begin(), l2.end(), l2.getRange().start);
        SeqIterType seqs_end(l2.end());

        //printf("tid %d first: pos in file %lu, l2 range start: %lu\n", tid, (*seqs_start).id.pos_in_file, l2.getRange().start);

        //== loop over the reads
        int i = 0;
//        size_t last_pos, last_seq_start, last_seq_end;
        for (; seqs_start != seqs_end; ++seqs_start, ++i)
        {
        	if (l2.getRange().start <= (*seqs_start).id.pos_in_file) ++localKmerCount;
        	// else printf("tid %d excluded:  pos in file %lu < l2 range start %lu, id = %d\n", tid, (*seqs_start).id.pos_in_file, l2.getRange().start, i);
          //      std::cout << "Rank " << rank << "/" << nprocs << " seq: " << (*seqs_start).id.id << ", ";
          //      std::cout << std::distance(l1.begin(), (*seqs_start).seqBegin) << ", " << std::distance(l1.begin(), (*seqs_start).seqEnd) << ", ";
          //      std::cout << std::distance(l1.begin(), (*seqs_start).qualBegin) << ", " << std::distance(l1.begin(), (*seqs_start).qualEnd) << std::endl;

//        	last_pos = (*seqs_start).id.pos_in_file;
//        	last_seq_start = last_pos + (*seqs_start).seq_begin_offset;
//        	last_seq_end = last_pos + (*seqs_start).length;
        }
        //printf("tid %d last: pos in file [%lu, %lu, %lu), l2 range start: %lu, i=%d\n", tid, last_pos, last_seq_start, last_seq_end, l2.getRange().start, i);

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
TEST_P(FASTAParserTest, parse_mpi)
{
  // from FileLoader type, get the block iter type and range type
  using BlockIterType = typename FASTALoaderType::L1BlockType::iterator;
  using SeqIterType = ::bliss::io::SequencesIterator<BlockIterType, bliss::io::FASTAParser >;

  FASTALoaderType loader(this->fileName, MPI_COMM_WORLD, 1, block_size, block_size );

  auto l1 = loader.getNextL1Block();
  this->elemCount = 0;

  bliss::io::FASTAParser<BlockIterType> l1parser = loader.getSeqParser();


  while (l1.getRange().size() > 0) {
    l1parser.init_for_iterator(l1.begin(), loader.getFileRange(), l1.getRange(), l1.getRange(), MPI_COMM_WORLD);


    //==  and wrap the chunk inside an iterator that emits Reads.
    SeqIterType seqs_start(l1parser, l1.begin(), l1.end(), l1.getRange().start);
    SeqIterType seqs_end(l1.end());

    //== loop over the reads
    for (; seqs_start != seqs_end; ++seqs_start)
    {
    	if (l1.getRange().start <= (*seqs_start).id.pos_in_file) ++(this->elemCount);

      //      std::cout << "Rank " << rank << "/" << nprocs << " seq: " << (*seqs_start).id.id << ", ";
      //      std::cout << std::distance(l1.begin(), (*seqs_start).seqBegin) << ", " << std::distance(l1.begin(), (*seqs_start).seqEnd) << ", ";
      //      std::cout << std::distance(l1.begin(), (*seqs_start).qualBegin) << ", " << std::distance(l1.begin(), (*seqs_start).qualEnd) << std::endl;
    }

    l1 = loader.getNextL1Block();
  }

}
#endif


#ifdef USE_MPI
#ifdef USE_OPENMP
TEST_P(FASTAParserTest, parse_mpi_omp)
{

  using BlockIterType = typename FASTALoaderType::L1BlockType::iterator;
  using SeqIterType = ::bliss::io::SequencesIterator<BlockIterType, bliss::io::FASTAParser >;


  int nthreads = 4;

  FASTALoaderType loader(this->fileName, MPI_COMM_WORLD, nthreads, block_size, block_size * nthreads * 2);

  auto l1 = loader.getNextL1Block();
  this->elemCount = 0;
  size_t localKmerCount = 0;


  bliss::io::FASTAParser<BlockIterType> l1parser = loader.getSeqParser();

  while (l1.getRange().size() > 0) {
    l1parser.init_for_iterator(l1.begin(), loader.getFileRange(), l1.getRange(), l1.getRange(), MPI_COMM_WORLD);


    localKmerCount = 0;
#pragma omp parallel num_threads(nthreads) shared(loader, l1parser, l1) reduction(+:localKmerCount)
   {
    	bliss::io::FASTAParser<typename FASTALoaderType::L2BlockType::iterator> l2parser = l1parser;
      int tid = omp_get_thread_num();

      auto l2 = loader.getNextL2Block(tid);

      while (l2.getRange().size() > 0) {

        // from FileLoader type, get the block iter type and range type


        //==  and wrap the chunk inside an iterator that emits Reads.
        SeqIterType seqs_start(l2parser, l2.begin(), l2.end(), l2.getRange().start);
        SeqIterType seqs_end(l2.end());

        //== loop over the reads
        for (; seqs_start != seqs_end; ++seqs_start)
        {
        	if (l2.getRange().start <= (*seqs_start).id.pos_in_file) ++localKmerCount;
          //      std::cout << "Rank " << rank << "/" << nprocs << " seq: " << (*seqs_start).id.id << ", ";
          //      std::cout << std::distance(l1.begin(), (*seqs_start).seqBegin) << ", " << std::distance(l1.begin(), (*seqs_start).seqEnd) << ", ";
          //      std::cout << std::distance(l1.begin(), (*seqs_start).qualBegin) << ", " << std::distance(l1.begin(), (*seqs_start).qualEnd) << std::endl;
        }

// /       printf("tid %d l1 range [%lu, %lu), l2 range [%lu, %lu) local kmer count: %lu\n", tid, l1.getRange().start, l1.getRange().end, l2.getRange().start, l2.getRange().end, localKmerCount);
        l2 = loader.getNextL2Block(tid);
      }

    }  // end omp parallel region.

   this->elemCount += localKmerCount;

   l1 = loader.getNextL1Block();
  }  // end l1 while.

}
#endif
#endif


INSTANTIATE_TEST_CASE_P(Bliss, FASTAParserTest, ::testing::Values(
    TestFileInfo(250, 14625, std::string("/test/data/natural.fasta")),
    TestFileInfo(4, 512, std::string("/test/data/test.fasta")),
    TestFileInfo(5000, 1092580, std::string("/test/data/test.medium.fasta")),
    TestFileInfo(6, 940, std::string("/test/data/test2.fasta"))
));

//INSTANTIATE_TEST_CASE_P(Bliss, FASTAParserTest, ::testing::Values(
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

