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

#include "io/fasta_loader.hpp"
#include "io/file_loader.hpp"
#include "common/kmer.hpp"
#include "common/alphabets.hpp"
//#include "io/fasta_iterator.hpp"
#include "io/mxx_support.hpp"
#include "partition/partitioner.hpp"
#include "partition/range.hpp"
#include "index/kmer_index.hpp"

#include "io/test/file_loader_test_fixtures.hpp"

using namespace bliss::io;


// NOTE:  there are a few problems with the old fasta iterator and loader
//        1. fasta_loader is not a subclass of fileloader, so it does not use the same configurePartition calls.
//        2. fasta_loader is not using overlap information
//        3. fasta_iterator is processing then filtering. so it's wasting time on headers
//        4. fasta_iterator is not handling the EOL characters.

static constexpr size_t block_size = 32768;

static constexpr size_t kmer_size = 35;
typedef FASTALoader<unsigned char, kmer_size - 1> FASTALoaderType;

class FASTAIteratorTest : public KmerReaderTest<FASTALoaderType > {};

// Single thread, single process
TEST_P(FASTAIteratorTest, read)
{
	  this->elemCount = 0;

#ifdef USE_MPI
	int rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0) {
#endif


  typedef typename FASTALoaderType::InputIteratorType InputIterType;
  typedef typename std::iterator_traits<InputIterType>::value_type ValueType;

  // from FileLoader type, get the block iter type and range type
  using BlockIterType = typename FASTALoaderType::L1BlockType::iterator;

  // define sequence parser
  using SeqIterType = ::bliss::io::SequencesIterator<BlockIterType, bliss::io::FASTAParser >;

  //Define Kmer parser
  typedef bliss::common::Kmer<kmer_size, bliss::common::DNA5, uint32_t> KmerType;
  typedef std::pair<KmerType, bliss::common::LongSequenceKmerId> TupleType;
  typedef std::vector<TupleType>   ResultVecType;
  bliss::index::kmer::KmerPositionTupleParser<TupleType > kmer_parser;

  //== process the chunk of data
  FASTALoaderType loader(this->fileName, 1, 0, 1, block_size, block_size);

  auto l1 = loader.getNextL1Block();
  auto l1parser = loader.getSeqParser();

  ValueType* gold = new ValueType[kmer_size+1];
  gold[kmer_size] = 0;

  ResultVecType result;
  bool same = true;
  bool local_same = true;

  while (l1.getRange().size() > 0) {
    l1parser.init_for_iterator(l1.begin(), loader.getFileRange(), l1.getRange(), l1.getRange());

    //==  and wrap the chunk inside an iterator that emits Reads.
    SeqIterType seqs_start(l1parser, l1.begin(), l1.end(), l1.getRange().start);
    SeqIterType seqs_end(l1.end());

    //== loop over the reads
    for (; seqs_start != seqs_end; ++seqs_start)
    {
    	auto seq = *seqs_start;

      result.clear();

      ::fsc::back_emplace_iterator<ResultVecType > emplace_iter(result);
      emplace_iter = kmer_parser(seq, emplace_iter);

      // compare the results.
//  	printf("sequence record: id %lu, offset %lu, local offset %lu, length %lu\n", seq.id.pos_in_file, seq.seq_begin_offset, seq.local_offset, seq.length);
      for (size_t i = 0; i < result.size(); ++i) {
        this->readFilePOSIX(this->fileName, result[i].second.get_pos(), kmer_size, gold);

        std::string KmertoString = bliss::utils::KmerUtils::toASCIIString(result[i].first);

        local_same = equal(KmertoString.begin(), gold, kmer_size);
        same &= local_same;

        if (!local_same) {
          printf("sequence record: id %lu, offset %lu, local offset %lu, length %lu\n", seq.id.pos_in_file, seq.seq_begin_offset, seq.local_offset, seq.length);
          printf("i %lu id: pos %lu, id %lu, file %d\n", i, result[i].second.get_pos(), result[i].second.get_id(), result[i].second.get_file_id());
          printf("i %lu pos %lu gold: [%s]\npos %lu test: [%s]\n", i, result[i].second.get_pos(), gold, result[i].second.get_pos(), KmertoString.c_str());
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
TEST_P(FASTAIteratorTest, read_mpi)
{
  typedef typename FASTALoaderType::InputIteratorType InputIterType;
  typedef typename std::iterator_traits<InputIterType>::value_type ValueType;

  // from FileLoader type, get the block iter type and range type
  using BlockIterType = typename FASTALoaderType::L1BlockType::iterator;

  // define sequence parser
  using SeqIterType = ::bliss::io::SequencesIterator<BlockIterType, bliss::io::FASTAParser >;

  //Define Kmer parser
  typedef bliss::common::Kmer<kmer_size, bliss::common::DNA5, uint32_t> KmerType;
  typedef std::pair<KmerType, bliss::common::LongSequenceKmerId> TupleType;
  typedef std::vector<TupleType>   ResultVecType;

  bliss::index::kmer::KmerPositionTupleParser<TupleType > kmer_parser;

  //== process the chunk of data
  FASTALoaderType loader(this->fileName, MPI_COMM_WORLD, 1, block_size, block_size);

  auto l1 = loader.getNextL1Block();
  auto l1parser = loader.getSeqParser();

  ValueType* gold = new ValueType[kmer_size+1];
  gold[kmer_size] = 0;

  ResultVecType result;

  bool same = true;
  this->elemCount = 0;

  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  while (l1.getRange().size() > 0) {
    l1parser.init_for_iterator(l1.begin(), loader.getFileRange(), l1.getRange(), l1.getRange(), MPI_COMM_WORLD);

//    printf("block range [%lu, %lu)\n", l1.getRange().start, l1.getRange().end);

    //==  and wrap the chunk inside an iterator that emits Reads.
    SeqIterType seqs_start(l1parser, l1.begin(), l1.end(), l1.getRange().start);
    SeqIterType seqs_end(l1.end());

    //== loop over the reads
    //printf("rank %d l1 block [%lu, %lu)\n", rank, l1.getRange().start, l1.getRange().end);

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

//      printf("result size: %lu\n", result.size());
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
TEST_P(FASTAIteratorTest, read_omp)
{
	  this->elemCount = 0;

#ifdef USE_MPI
	int rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0) {
#endif
  typedef typename FASTALoaderType::InputIteratorType InputIterType;
  typedef typename std::iterator_traits<InputIterType>::value_type ValueType;

  // from FileLoader type, get the block iter type and range type
  using BlockIterType = typename FASTALoaderType::L2BlockType::iterator;

  // define sequence parser
  using SeqIterType = ::bliss::io::SequencesIterator<BlockIterType, bliss::io::FASTAParser >;

  //Define Kmer parser
  typedef bliss::common::Kmer<kmer_size, bliss::common::DNA5, uint32_t> KmerType;
  typedef std::pair<KmerType, bliss::common::LongSequenceKmerId> TupleType;
  typedef std::vector<TupleType>   ResultVecType;



  int nthreads = 4;

  //== process the chunk of data
  FASTALoaderType loader(this->fileName, 1, 0, nthreads, block_size, block_size * nthreads * 2 );

  bool same = true;
  bool local_same[nthreads];
  std::string filename = this->fileName;

  auto l1 = loader.getNextL1Block();
  bliss::io::FASTAParser<BlockIterType> l1parser = loader.getSeqParser();

  while (l1.getRange().size() > 0) {

    l1parser.init_for_iterator(l1.begin(), loader.getFileRange(), l1.getRange(), l1.getRange());


    size_t localKmerCount = 0;


#pragma omp parallel num_threads(nthreads) shared(filename, local_same, loader, l1parser) reduction(+:localKmerCount)
    {
      int tid = omp_get_thread_num();

      bliss::index::kmer::KmerPositionTupleParser<TupleType > kmer_parser;

      bliss::io::FASTAParser<typename FASTALoaderType::L2BlockType::iterator> l2parser = l1parser;

      bool threadcomp = true;
      local_same[tid] = true;

      //Get L2 block for this thread
      auto l2 = loader.getNextL2Block(tid);

      ValueType* gold = new ValueType[kmer_size+1];
      gold[kmer_size] = 0;

      ResultVecType result;
      ::fsc::back_emplace_iterator<ResultVecType > emplace_iter(result);

      while (l2.getRange().size() > 0) {

        //==  and wrap the chunk inside an iterator that emits Reads.
        SeqIterType seqs_start(l2parser, l2.begin(), l2.end(), l2.getRange().start);
        SeqIterType seqs_end(l2.end());

        //== loop over the reads

        for (; seqs_start != seqs_end; ++seqs_start)
        {
          auto seq = *seqs_start;
          result.clear();

          // generate kmers and save them
          emplace_iter = kmer_parser(seq, emplace_iter);

          // compare the results.
          for (size_t i = 0; i < result.size(); ++i) {
            this->readFilePOSIX(filename, result[i].second.get_pos(), kmer_size, gold);

            std::string KmertoString = bliss::utils::KmerUtils::toASCIIString(result[i].first);

            threadcomp = equal(KmertoString.begin(), gold, kmer_size);

            if (!threadcomp) {
              BL_DEBUGF("tid %d l2 block: start, %lu end %lu\n", tid, l2.getRange().start, l2.getRange().end);
              BL_DEBUGF("tid %d sequence record: id %lu, offset %lu, local offset %lu, length %lu\n", tid, seq.id.pos_in_file, seq.seq_begin_offset, seq.local_offset, seq.length);
              BL_DEBUGF("tid %d i %lu id: pos %lu, id %lu, file %d\n", tid, i, result[i].second.get_pos(), result[i].second.get_id(), result[i].second.get_file_id());
              BL_DEBUGF("tid %d\ti %lu\tpos %lu gold: [%s]\n\t\tpos %lu test: [%s]\n", tid,  i, result[i].second.get_pos(), gold, result[i].second.get_pos(), KmertoString.c_str());
            }

            local_same[tid] &= threadcomp;
          }  // end kmers for

          localKmerCount += result.size();
        }  // end sequences for

        l2 = loader.getNextL2Block(tid);
      }  // end L2 while

      delete [] gold;

    }  // end omp parallel

    for (int i = 0; i < nthreads; ++i) {
      same &= local_same[i];
    }

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
TEST_P(FASTAIteratorTest, read_omp_mpi)
{
  typedef typename FASTALoaderType::InputIteratorType InputIterType;
  typedef typename std::iterator_traits<InputIterType>::value_type ValueType;

  // from FileLoader type, get the block iter type and range type
  using BlockIterType = typename FASTALoaderType::L2BlockType::iterator;

  // define sequence parser
  using SeqIterType = ::bliss::io::SequencesIterator<BlockIterType, bliss::io::FASTAParser >;

  //Define Kmer parser
  typedef bliss::common::Kmer<kmer_size, bliss::common::DNA5, uint32_t> KmerType;
  typedef std::pair<KmerType, bliss::common::LongSequenceKmerId> TupleType;
  typedef std::vector<TupleType>   ResultVecType;



  int nthreads = 4;

  //== process the chunk of data
  FASTALoaderType loader(this->fileName, MPI_COMM_WORLD, nthreads, block_size, block_size * 2 * nthreads );

  bool same = true;
  this->elemCount = 0;
  std::string filename = this->fileName;

  auto l1 = loader.getNextL1Block();
  bliss::io::FASTAParser<BlockIterType> l1parser = loader.getSeqParser();

  while (l1.getRange().size() > 0) {

    l1parser.init_for_iterator(l1.begin(), loader.getFileRange(), l1.getRange(), l1.getRange(), MPI_COMM_WORLD);
    bliss::io::FASTAParser<typename FASTALoaderType::L2BlockType::iterator> l2parser = l1parser;

    size_t localKmerCount = 0;


#pragma omp parallel num_threads(nthreads) shared(filename, same, loader, l2parser) reduction(+:localKmerCount)
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
        SeqIterType seqs_start(l2parser, l2.begin(), l2.end(), l2.getRange().start);
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

INSTANTIATE_TEST_CASE_P(Bliss, FASTAIteratorTest, ::testing::Values(
    TestFileInfo(500, 14625, std::string("/test/data/natural.fasta")),
    TestFileInfo(246, 940, std::string("/test/data/test2.fasta")),
    TestFileInfo(335000, 1092580, std::string("/test/data/test.medium.fasta")),
    TestFileInfo(64, 512, std::string("/test/data/test.fasta"))
));
//INSTANTIATE_TEST_CASE_P(Bliss, FASTAIteratorTest, ::testing::Values(
//    TestFileInfo(246, 940, std::string("/test/data/test2.fasta")),
//    TestFileInfo(64, 512, std::string("/test/data/test.fasta"))
//    ));

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
