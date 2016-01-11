/**
 * fastaIterator_test.cpp
 *
 *  Created on: Nov 14, 2014
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
#include <iostream>

#include "io/fastq_loader.hpp"
#include "io/file_loader.hpp"
#include "common/kmer.hpp"
#include "common/alphabets.hpp"
#include "io/mxx_support.hpp"
#include "partition/partitioner.hpp"
#include "partition/range.hpp"
#include "index/kmer_index.hpp"

#include "io/test/file_loader_test_fixtures.hpp"

using namespace bliss::io;


static constexpr size_t kmer_size = 35;
typedef FASTQLoader<unsigned char> FASTQLoaderType;

class FASTQIteratorTest : public KmerReaderTest<FASTQLoaderType > {};

// Single thread, single process
TEST_P(FASTQIteratorTest, read)
{
	  this->elemCount = 0;

#ifdef USE_MPI
	int rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0) {
#endif

  typedef typename FASTQLoaderType::InputIteratorType InputIterType;
  typedef typename std::iterator_traits<InputIterType>::value_type ValueType;

  // from FileLoader type, get the block iter type and range type
  using BlockIterType = typename FASTQLoaderType::L1BlockType::iterator;

  // define sequence parser
  using SeqIterType = ::bliss::io::SequencesIterator<BlockIterType, bliss::io::FASTQParser >;

  //Define Kmer parser
  typedef bliss::common::Kmer<kmer_size, bliss::common::DNA5, uint32_t> KmerType;
  typedef std::pair<KmerType, bliss::common::ShortSequenceKmerId> TupleType;
  typedef std::vector<TupleType>   ResultVecType;
  bliss::index::kmer::KmerPositionTupleParser<TupleType > kmer_parser;

  //== process the chunk of data
  FASTQLoaderType loader(this->fileName, 1, 0, 1, 2048, 2048 );

  auto l1 = loader.getNextL1Block();

  ValueType* gold = new ValueType[kmer_size+1];
  gold[kmer_size] = 0;

  ResultVecType result;
  bool same = true;
  bool local_same;

  while (l1.getRange().size() > 0) {

    //==  and wrap the chunk inside an iterator that emits Reads.
    SeqIterType seqs_start(bliss::io::FASTQParser<BlockIterType>(), l1.begin(), l1.end(), l1.getRange().start);
    SeqIterType seqs_end(l1.end());

    //== loop over the reads

    for (; seqs_start != seqs_end; ++seqs_start)
    {
    	auto seq = *seqs_start;
      result.clear();

//      printf("sequence record: id %lu, offset %lu, local offset %lu, length %lu\n", seq.id.pos_in_file, seq.seq_begin_offset, seq.seq_offset, seq.record_size);

      ::fsc::back_emplace_iterator<ResultVecType > emplace_iter(result);
      emplace_iter = kmer_parser(seq, emplace_iter);

      // compare the results.
      for (size_t i = 0; i < result.size(); ++i) {
        this->readFilePOSIX(this->fileName, result[i].second.get_pos(), kmer_size, gold);

        std::string KmertoString = bliss::utils::KmerUtils::toASCIIString(result[i].first);

        local_same = equal(KmertoString.begin(), gold, kmer_size);
        same &= local_same;

        if (!local_same) {
        	BL_ERROR("sequence record: " << seq);
          BL_ERRORF("i %lu id: pos %lu, id %lu, file %d\n", i, result[i].second.get_pos(), result[i].second.get_id(), result[i].second.get_file_id());
          BL_ERRORF("i %lu pos %lu gold: [%s]\npos %lu test: [%s]\n", i, result[i].second.get_pos(), gold, result[i].second.get_pos(), KmertoString.c_str());
        }
      }

      this->elemCount += result.size();
    }

    l1 = loader.getNextL1Block();
  }
  EXPECT_TRUE(same);

  delete [] gold;
#ifdef USE_MPI
	}
#endif
}

#ifdef USE_MPI
// Single thread, multiple processes
TEST_P(FASTQIteratorTest, read_mpi)
{
  typedef typename FASTQLoaderType::InputIteratorType InputIterType;
  typedef typename std::iterator_traits<InputIterType>::value_type ValueType;

  // from FileLoader type, get the block iter type and range type
  using BlockIterType = typename FASTQLoaderType::L1BlockType::iterator;

  // define sequence parser
  using SeqIterType = ::bliss::io::SequencesIterator<BlockIterType, bliss::io::FASTQParser >;

  //Define Kmer parser
  typedef bliss::common::Kmer<kmer_size, bliss::common::DNA5, uint32_t> KmerType;
  typedef std::pair<KmerType, bliss::common::ShortSequenceKmerId> TupleType;
  typedef std::vector<TupleType>   ResultVecType;

  bliss::index::kmer::KmerPositionTupleParser<TupleType > kmer_parser;

  //== process the chunk of data
  FASTQLoaderType loader(this->fileName, MPI_COMM_WORLD, 1, 2048, 2048 );

  auto l1 = loader.getNextL1Block();

  ValueType* gold = new ValueType[kmer_size+1];
  gold[kmer_size] = 0;

  ResultVecType result;

  bool same = true;
  this->elemCount = 0;

  while (l1.getRange().size() > 0) {


    //==  and wrap the chunk inside an iterator that emits Reads.
    SeqIterType seqs_start(bliss::io::FASTQParser<BlockIterType>(), l1.begin(), l1.end(), l1.getRange().start);
    SeqIterType seqs_end(l1.end());

    //== loop over the reads

    for (; seqs_start != seqs_end; ++seqs_start)
    {
      result.clear();

      ::fsc::back_emplace_iterator<ResultVecType > emplace_iter(result);
      emplace_iter = kmer_parser(*seqs_start, emplace_iter);

      // compare the results.
      for (size_t i = 0; i < result.size(); ++i) {
        this->readFilePOSIX(this->fileName, result[i].second.get_pos(), kmer_size, gold);

        std::string KmertoString = bliss::utils::KmerUtils::toASCIIString(result[i].first);

        same &= equal(KmertoString.begin(), gold, kmer_size);
      }

      this->elemCount += result.size();
    }

    l1 = loader.getNextL1Block();
  }

  EXPECT_TRUE(same);

  delete [] gold;

}
#endif


#ifdef USE_OPENMP
// Single process, multiple threads
TEST_P(FASTQIteratorTest, read_omp)
{
  this->elemCount = 0;

#ifdef USE_MPI
	int rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0) {
#endif

  typedef typename FASTQLoaderType::InputIteratorType InputIterType;
  typedef typename std::iterator_traits<InputIterType>::value_type ValueType;

  // from FileLoader type, get the block iter type and range type
  using BlockIterType = typename FASTQLoaderType::L2BlockType::iterator;

  // define sequence parser
  using SeqIterType = ::bliss::io::SequencesIterator<BlockIterType, bliss::io::FASTQParser >;

  //Define Kmer parser
  typedef bliss::common::Kmer<kmer_size, bliss::common::DNA5, uint32_t> KmerType;
  typedef std::pair<KmerType, bliss::common::ShortSequenceKmerId> TupleType;
  typedef std::vector<TupleType>   ResultVecType;



  int nthreads = 4;

  //== process the chunk of data
  FASTQLoaderType loader(this->fileName, 1, 0, nthreads, 2048, 2048 * nthreads * 2 );

  bool same = true;
  size_t localKmerCount = 0;
  std::string filename = this->fileName;

  auto l1 = loader.getNextL1Block();

  while (l1.getRange().size() > 0) {

    localKmerCount = 0;


#pragma omp parallel num_threads(nthreads) shared(filename, same, loader) reduction(+:localKmerCount)
    {
      int tid = omp_get_thread_num();
      bliss::index::kmer::KmerPositionTupleParser<TupleType > kmer_parser;

      bool localcomp = true;

      //Get L2 block for this thread
      auto l2 = loader.getNextL2Block(tid);

      ValueType* gold = new ValueType[kmer_size+1];
      gold[kmer_size] = 0;

      ResultVecType result;
      ::fsc::back_emplace_iterator<ResultVecType > emplace_iter(result);

      while (l2.getRange().size() > 0) {

        //==  and wrap the chunk inside an iterator that emits Reads.
        SeqIterType seqs_start(bliss::io::FASTQParser<BlockIterType>(), l2.begin(), l2.end(), l2.getRange().start);
        SeqIterType seqs_end(l2.end());

        //== loop over the reads

        for (; seqs_start != seqs_end; ++seqs_start)
        {
          result.clear();

          // generate kmers and save them
          kmer_parser(*seqs_start, emplace_iter);

          // compare the results.
          for (size_t i = 0; i < result.size(); ++i) {
            this->readFilePOSIX(filename, result[i].second.get_pos(), kmer_size, gold);

            std::string KmertoString = bliss::utils::KmerUtils::toASCIIString(result[i].first);

            localcomp &= equal(KmertoString.begin(), gold, kmer_size);
          }  // end kmers for

          localKmerCount += result.size();
        }  // end sequences for

        l2 = loader.getNextL2Block(tid);
      }  // end L2 while

      delete [] gold;

#pragma omp critical
      same &= localcomp;


    }  // end omp parallel

    this->elemCount += localKmerCount;

    l1 = loader.getNextL1Block();
  }// end L1 while


  EXPECT_TRUE(same);
#ifdef USE_MPI
	}
#endif
}
#endif




// Single process, multiple threads
#ifdef USE_MPI
#ifdef USE_OPENMP
TEST_P(FASTQIteratorTest, read_omp_mpi)
{
  typedef typename FASTQLoaderType::InputIteratorType InputIterType;
  typedef typename std::iterator_traits<InputIterType>::value_type ValueType;

  // from FileLoader type, get the block iter type and range type
  using BlockIterType = typename FASTQLoaderType::L2BlockType::iterator;

  // define sequence parser
  using SeqIterType = ::bliss::io::SequencesIterator<BlockIterType, bliss::io::FASTQParser >;

  //Define Kmer parser
  typedef bliss::common::Kmer<kmer_size, bliss::common::DNA5, uint32_t> KmerType;
  typedef std::pair<KmerType, bliss::common::ShortSequenceKmerId> TupleType;
  typedef std::vector<TupleType>   ResultVecType;


  int nthreads = 4;

  //== process the chunk of data
  FASTQLoaderType loader(this->fileName, MPI_COMM_WORLD, nthreads, 2048, 2048 * nthreads * 2);

  bool same = true;
  this->elemCount = 0;
  size_t localKmerCount = 0;
  std::string filename = this->fileName;

  auto l1 = loader.getNextL1Block();

  while (l1.getRange().size() > 0) {

    localKmerCount = 0;


#pragma omp parallel num_threads(nthreads) shared(filename, same, loader) reduction(+:localKmerCount)
    {
      int tid = omp_get_thread_num();
      bliss::index::kmer::KmerPositionTupleParser<TupleType > kmer_parser;

      bool localcomp = true;

      //Get L2 block for this thread
      auto l2 = loader.getNextL2Block(tid);

      ValueType* gold = new ValueType[kmer_size+1];
      gold[kmer_size] = 0;

      ResultVecType result;
      ::fsc::back_emplace_iterator<ResultVecType > emplace_iter(result);

      while (l2.getRange().size() > 0) {

        //==  and wrap the chunk inside an iterator that emits Reads.
        SeqIterType seqs_start(bliss::io::FASTQParser<BlockIterType>(), l2.begin(), l2.end(), l2.getRange().start);
        SeqIterType seqs_end(l2.end());

        //== loop over the reads

        for (; seqs_start != seqs_end; ++seqs_start)
        {
          result.clear();

          // generate kmers and save them
          kmer_parser(*seqs_start, emplace_iter);

          // compare the results.
          for (size_t i = 0; i < result.size(); ++i) {
            this->readFilePOSIX(filename, result[i].second.get_pos(), kmer_size, gold);

            std::string KmertoString = bliss::utils::KmerUtils::toASCIIString(result[i].first);

            localcomp &= equal(KmertoString.begin(), gold, kmer_size);
          }  // end kmers for

          localKmerCount += result.size();
        }  // end sequences for

        l2 = loader.getNextL2Block(tid);
      }  // end L2 while

      delete [] gold;

#pragma omp critical
      same &= localcomp;


    }  // end omp parallel

    this->elemCount += localKmerCount;

    l1 = loader.getNextL1Block();
  }// end L1 while


  EXPECT_TRUE(same);
}
#endif
#endif

// first number is number of kmers (duplicate or not).  second number is file size.
INSTANTIATE_TEST_CASE_P(Bliss, FASTQIteratorTest, ::testing::Values(
    TestFileInfo(500, 29250, std::string("/test/data/natural.fastq")),
    TestFileInfo(3640, 18761, std::string("/test/data/test.medium.fastq")),
    TestFileInfo(182, 939, std::string("/test/data/test.small.fastq")),
    TestFileInfo(6618612, 34111308, std::string("/test/data/test.fastq"))
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
