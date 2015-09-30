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

#ifdef USE_OPENMP
#include "omp.h"
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
#include "io/fasta_iterator.hpp"
#include "io/mxx_support.hpp"
#include "partition/partitioner.hpp"
#include "partition/range.hpp"

using namespace bliss::io;


// TODO:  there are a few problems with fasta iterator and loader.
//        1. fasta_loader is not a subclass of fileloader, so it does not use the same configurePartition calls.
//        2. fasta_loader is not use overlap information
//        3. fasta_iterator is processing then filtering. so it's wasting time on headers
//        4. fasta_iterator is not handling the EOL characters.

#if 0

template<typename Iter1, typename Iter2>
bool equal(const Iter1 &i1, const Iter2 &i2, size_t len) {
  Iter1 ii1(i1);
  Iter2 ii2(i2);
  for (size_t i = 0; i < len; ++i, ++ii1, ++ii2) {
    if (*ii1 != *ii2) return false;
  }
  return true;
}

struct TestFileInfo {
    size_t totalKmerCount;
    size_t fileSize;
    std::string filename;

    TestFileInfo(size_t const & _kmer_count, size_t const & _file_size, std::string const & _file_name) : totalKmerCount(_kmer_count), fileSize(_file_size), filename(_file_name) {};
};

class FASTAIteratorTest2 : public ::testing::TestWithParam<TestFileInfo>
{
  protected:
    virtual void SetUp()
    {
      TestFileInfo const & p = GetParam();

      kmerCount = 0;

      fileName.assign(PROJ_SRC_DIR);
      fileName.append(p.filename);

      // get file size
      struct stat filestat;
      stat(fileName.c_str(), &filestat);
      size_t fileSize = static_cast<size_t>(filestat.st_size);

      ASSERT_EQ(p.fileSize, fileSize);
      //ASSERT_EQ(34526831, fileSize);
      //ASSERT_EQ(114, fileSize);
    }

    virtual void TearDown() {
      TestFileInfo const & p = GetParam();

      size_t totalKmerCount = mxx::allreduce(kmerCount, std::plus<size_t>());
      EXPECT_EQ(p.totalKmerCount, totalKmerCount);
    }

    static void readFilePOSIX(const std::string &fileName, const size_t& offset,
                              const size_t& length, unsigned char* result)
    {
      FILE *fp = fopen(fileName.c_str(), "r");
      fseek(fp, offset * sizeof(unsigned char), SEEK_SET);
      unsigned char c;
      unsigned char* curr = result;
      int read = 0;
      size_t total_read = 0;
      for (size_t i = 0; i < length;) {
        read = fread_unlocked(&c, sizeof(unsigned char), 1, fp);
        total_read += read;
        if (read > 0) {  // only if read something.
			if ((c >= 64 && c <= 90) || (c >= 97 && c <=122)) {
			  *curr = c;
			  ++curr;
			  ++i;
			}
        }
      }
      fclose(fp);
      EXPECT_GT(0UL, total_read);

    }


    std::string fileName;
    size_t kmerCount;
    static const int kmerSize = 35;

};

// Single thread, single process
TEST_P(FASTAIteratorTest2, read)
{

	this->kmerCount = 0;

#ifdef USE_MPI
	int rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0) {
#endif


  typedef FileLoader<unsigned char, this->kmerSize, bliss::io::BaseFileParser, false, false> FILELoaderType;
  typedef typename FILELoaderType::InputIteratorType InputIterType;
  typedef typename std::iterator_traits<InputIterType>::value_type ValueType;

  //Define Kmer type
  typedef bliss::common::Kmer<this->kmerSize, bliss::common::DNA5, uint32_t> KmerType;

  //Define type for FASTALoader 
  typedef FASTALoader<typename FILELoaderType::InputIteratorType, KmerType> FASTALoaderType;
  typename FASTALoaderType::vectorType vectorReturned;

  typedef typename FILELoaderType::RangeType RangeType;
  typename FILELoaderType::L1BlockType d;


  int nprocs = 1;
  int rank = 0;

  //Load the file in FILE Loader
  FILELoaderType loader(this->fileName, nprocs, rank);

  d = loader.getNextL1Block();
  RangeType r = d.getRange();

  //Get the vector of fasta headers from FASTA Loader
  FASTALoaderType obj(MPI_COMM_WORLD);
  obj.countSequenceStarts(d.begin(), loader.getFileRange() , r, vectorReturned);
  EXPECT_LT(0UL , vectorReturned.size());


  //Begin iteration
  FASTAParser<typename FILELoaderType::L1BlockType::iterator, KmerType> FASTAIterObj(d.begin(), d.end(), r, loader.getFileRange(), vectorReturned);
  ValueType* gold = new ValueType[this->kmerSize];

  //For crosschecking the contents
  for(auto iter=FASTAIterObj.begin(); iter != FASTAIterObj.end(); iter++)
  {
    std::string KmertoString = bliss::utils::KmerUtils::toASCIIString(FASTAIterObj.getKmer(iter));
    //std::cout << FASTAIterObj.getOffset(iter) << ", Kmer# " << totalKmerCount+1 << " :  " << KmertoString << std::endl;
    auto offset = FASTAIterObj.getOffset(iter);

    //Fetch Kmer using simple file reading 
    this->readFilePOSIX(this->fileName, offset, this->kmerSize, gold);

    //Compare the contents, by constructing a new kmer and equating
    std::string directFromFile((char *)gold, this->kmerSize);
    //std::cout << "Should match with " << " :  " << directFromFile << std::endl; 
    this->kmerCount++;


    bool comp = (KmertoString == directFromFile);

    if (!comp) {
      printf("not same:\n\tkmer:\t%s\n\tfile:\t%s\n", KmertoString.c_str(), directFromFile.c_str());
    }
    EXPECT_EQ(true, comp);
  }

  delete [] gold;
#ifdef USE_MPI
	}
#endif
}

#ifdef USE_MPI
// Single thread, multiple processes
TEST_P(FASTAIteratorTest2, read_mpi)
{
  typedef FileLoader<unsigned char, this->kmerSize, bliss::io::BaseFileParser, false, false> FILELoaderType;
  typedef typename FILELoaderType::InputIteratorType InputIterType;
  typedef typename std::iterator_traits<InputIterType>::value_type ValueType;

  //Define Kmer type
  typedef bliss::common::Kmer<this->kmerSize, bliss::common::DNA5, uint32_t> KmerType;

  typedef FASTALoader<typename FILELoaderType::InputIteratorType, KmerType> FASTALoaderType;
  typename FASTALoaderType::vectorType vectorReturned;

  typedef typename FILELoaderType::RangeType RangeType;
  typename FILELoaderType::L1BlockType d;

  //Set the rank and number of processors using MPI
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  //Load the file in FILE Loader
  FILELoaderType loader(this->fileName, MPI_COMM_WORLD);

  d = loader.getNextL1Block();
  RangeType r = d.getRange();

  //Get the vector of fasta headers from FASTA Loader
  FASTALoaderType obj(MPI_COMM_WORLD);
  obj.countSequenceStarts(d.begin(), loader.getFileRange() , r, vectorReturned);
  EXPECT_LT(0UL, vectorReturned.size());

  //Begin iteration
  FASTAParser<typename FILELoaderType::L1BlockType::iterator, KmerType> FASTAIterObj(d.begin(), d.end(), r, loader.getFileRange(), vectorReturned);
  ValueType* gold;
  gold = new ValueType[this->kmerSize];

  //For crosschecking the contents
  for(auto iter=FASTAIterObj.begin(); iter != FASTAIterObj.end(); iter++)
  {
    std::string KmertoString = bliss::utils::KmerUtils::toASCIIString(FASTAIterObj.getKmer(iter));
    //std::cout << FASTAIterObj.getOffset(iter) << ", RANK#,Kmer# " << rank << "," << totalKmerCount+1 << " :  " << KmertoString << std::endl; 
    auto offset = FASTAIterObj.getOffset(iter);

    //Fetch Kmer using simple file reading 
    this->readFilePOSIX(this->fileName, offset, this->kmerSize, gold);

    //Compare the contents, by constructing a new kmer and equating
    std::string directFromFile((char *)gold, this->kmerSize);
    //std::cout << "Should match with " << " :  " << directFromFile << std::endl; 
    this->kmerCount++;

    bool comp = (KmertoString == directFromFile);
    EXPECT_EQ(true, comp);
  }

  delete [] gold;
  MPI_Barrier(MPI_COMM_WORLD);
}
#endif

#ifdef USE_OPENMP
// Single process, multiple threads
TEST_P(FASTAIteratorTest2, read_omp_demand)
{
	  this->kmerCount = 0;

#ifdef USE_MPI
	int rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0) {
#endif

  //Define L2 partitioner type to be block partition
  //TODO: Make L2BlockPartitioner threadsafe
  typedef FileLoader<unsigned char, this->kmerSize, bliss::io::BaseFileParser, false, false> FILELoaderType;
  typedef typename FILELoaderType::InputIteratorType InputIterType;
  typedef typename std::iterator_traits<InputIterType>::value_type ValueType;

  //Define Kmer type
  typedef bliss::common::Kmer<this->kmerSize, bliss::common::DNA5, uint32_t> KmerType;

  typedef FASTALoader<typename FILELoaderType::InputIteratorType, KmerType> FASTALoaderType;
  typename FASTALoaderType::vectorType vectorReturned;

  typedef typename FILELoaderType::RangeType RangeType;
  typename FILELoaderType::L1BlockType d;

  int rank = 0;
  int nprocs = 1;
  int nthreads = 4;

  //Load the file in FILE Loader
  FILELoaderType loader(this->fileName, nprocs, rank, nthreads, 1024);

  d = loader.getNextL1Block();
  RangeType r = d.getRange();

  //Get the vector of fasta headers from FASTA Loader
  FASTALoaderType obj(MPI_COMM_WORLD);
  obj.countSequenceStarts(d.begin(), loader.getFileRange() , r, vectorReturned);

  EXPECT_LT(0UL, vectorReturned.size());

  size_t localKmerCount = 0;
  bool localcomp[4];
  bool comp;

#pragma omp parallel num_threads(4) shared(comp, localcomp, loader) reduction(+:localKmerCount)
  {
    int tid = omp_get_thread_num();

    localcomp[tid] = true;

    //Get L2 block for this thread
    typename FILELoaderType::L2BlockType d2 = loader.getNextL2Block(tid);

    RangeType r2 = d2.getRange();
    ValueType* gold = new ValueType[this->kmerSize];
    while(r2.end - r2.start > 0)
    {
      //std::cout << "Thread# " << omp_get_thread_num() << ", range: " << r2.start << "," << r2.end << std::endl;

      //Begin iteration
      FASTAParser<typename FILELoaderType::L1BlockType::iterator, KmerType> FASTAIterObj(d2.begin(), d2.end(), r2, loader.getFileRange(), vectorReturned);

      //For crosschecking the contents
      for(auto iter=FASTAIterObj.begin(); iter != FASTAIterObj.end(); iter++)
      {
        std::string KmertoString = bliss::utils::KmerUtils::toASCIIString(FASTAIterObj.getKmer(iter));
        //std::cout << FASTAIterObj.getOffset(iter) << ", Thread#,Kmer# " << omp_get_thread_num() << "," << totalKmerCount+1 << " :  " << KmertoString << std::endl; 
        auto offset = FASTAIterObj.getOffset(iter);

        //Fetch Kmer using simple file reading 
        this->readFilePOSIX(this->fileName, offset, this->kmerSize, gold);

        //Compare the contents, by constructing a new kmer and equating
        std::string directFromFile((char *)gold, this->kmerSize);
        //std::cout << "Should match with " << " :  " << directFromFile<< std::endl; 

        ++localKmerCount;

        localcomp[tid] &= (KmertoString == directFromFile);
      }

      //Get next L2 block
      d2 = loader.getNextL2Block(tid);
      r2 = d2.getRange();
    }
    delete [] gold;

#pragma omp critical
    comp &= localcomp[tid];


  } // end omp pragma

  this->kmerCount = localKmerCount;

  //Total kmers discovered should match with #kmers found in sequential test
  EXPECT_EQ(true, comp);

#ifdef USE_MPI
	}
#endif
}
#endif

// Single process, multiple threads
#ifdef USE_OPENMP
TEST_P(FASTAIteratorTest2, read_omp_block)
{
	  this->kmerCount = 0;

#ifdef USE_MPI
	int rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0) {
#endif

  //Define L2 partitioner type to be block partition
  typedef bliss::partition::BlockPartitioner<bliss::partition::range<size_t> > L2BlockPartitionType;
  typedef FileLoader<unsigned char, this->kmerSize, bliss::io::BaseFileParser, false, false, L2BlockPartitionType> FILELoaderType;
  typedef typename FILELoaderType::InputIteratorType InputIterType;
  typedef typename std::iterator_traits<InputIterType>::value_type ValueType;

  //Define Kmer type
  typedef bliss::common::Kmer<this->kmerSize, bliss::common::DNA5, uint32_t> KmerType;

  typedef FASTALoader<typename FILELoaderType::InputIteratorType, KmerType> FASTALoaderType;
  typename FASTALoaderType::vectorType vectorReturned;

  typedef typename FILELoaderType::RangeType RangeType;
  typename FILELoaderType::L1BlockType d;

  int rank = 0;
  int nprocs = 1;
  int nthreads = 4;

  //Load the file in FILE Loader
  FILELoaderType loader(this->fileName, nprocs, rank, nthreads, 1024);

  d = loader.getNextL1Block();
  RangeType r = d.getRange();

  //Get the vector of fasta headers from FASTA Loader
  FASTALoaderType obj(MPI_COMM_WORLD);
  obj.countSequenceStarts(d.begin(), loader.getFileRange() , r, vectorReturned);

  EXPECT_LT(0UL , vectorReturned.size());

  size_t localKmerCount = 0;
  bool localcomp[4];
  bool comp;

#pragma omp parallel num_threads(4) shared(comp, localcomp, loader) reduction(+:localKmerCount)
  {
    int tid = omp_get_thread_num();


    localcomp[tid] = true;

    //Get L2 block for this thread
    typename FILELoaderType::L2BlockType d2 = loader.getNextL2Block(tid);

    RangeType r2 = d2.getRange();
    ValueType* gold = new ValueType[this->kmerSize];
    while(r2.end - r2.start > 0)
    {
      //std::cout << "Thread# " << omp_get_thread_num() << ", range: " << r2.start << "," << r2.end << std::endl;

      //Begin iteration
      FASTAParser<typename FILELoaderType::L1BlockType::iterator, KmerType> FASTAIterObj(d2.begin(), d2.end(), r2, loader.getFileRange(), vectorReturned);

      //For crosschecking the contents
      for(auto iter=FASTAIterObj.begin(); iter != FASTAIterObj.end(); iter++)
      {
        std::string KmertoString = bliss::utils::KmerUtils::toASCIIString(FASTAIterObj.getKmer(iter));
        //std::cout << FASTAIterObj.getOffset(iter) << ", Thread#,Kmer# " << omp_get_thread_num() << "," << totalKmerCount+1 << " :  " << KmertoString << std::endl;
        auto offset = FASTAIterObj.getOffset(iter);

        //Fetch Kmer using simple file reading
        this->readFilePOSIX(this->fileName, offset, this->kmerSize, gold);

        //Compare the contents, by constructing a new kmer and equating
        std::string directFromFile((char *)gold, this->kmerSize);
        //std::cout << "Should match with " << " :  " << directFromFile<< std::endl;

        ++localKmerCount;

        localcomp[tid] &= (KmertoString == directFromFile);
      }

      //Get next L2 block
      d2 = loader.getNextL2Block(tid);
      r2 = d2.getRange();
    }
    delete [] gold;

#pragma omp critical
    comp &= localcomp[tid];


  } // end omp pragma

  this->kmerCount = localKmerCount;

  //Total kmers discovered should match with #kmers found in sequential test
  EXPECT_EQ(true, comp);
#ifdef USE_MPI
	}
#endif
}
#endif



// Single process, multiple threads
#ifdef USE_MPI
#ifdef USE_OPENMP
TEST_P(FASTAIteratorTest2, read_omp_mpi)
{
  //Define L2 partitioner type to be block partition
  typedef FileLoader<unsigned char, this->kmerSize, bliss::io::BaseFileParser, false, false> FILELoaderType;
  typedef typename FILELoaderType::InputIteratorType InputIterType;
  typedef typename std::iterator_traits<InputIterType>::value_type ValueType;

  //Define Kmer type
  typedef bliss::common::Kmer<this->kmerSize, bliss::common::DNA5, uint32_t> KmerType;

  typedef FASTALoader<typename FILELoaderType::InputIteratorType, KmerType> FASTALoaderType;
  typename FASTALoaderType::vectorType vectorReturned;

  typedef typename FILELoaderType::RangeType RangeType;
  typename FILELoaderType::L1BlockType d;

  int rank = 0;
  int nprocs = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  int nthreads = 4;

  //Load the file in FILE Loader
  FILELoaderType loader(this->fileName, MPI_COMM_WORLD, nthreads, 1024);

  d = loader.getNextL1Block();
  RangeType r = d.getRange();

  //Get the vector of fasta headers from FASTA Loader
  FASTALoaderType obj(MPI_COMM_WORLD);
  obj.countSequenceStarts(d.begin(), loader.getFileRange() , r, vectorReturned);

  EXPECT_LT(0UL, vectorReturned.size());

  size_t localKmerCount = 0;
  bool localcomp[4];
  bool comp;

#pragma omp parallel num_threads(4) shared(comp, localcomp, loader) reduction(+:localKmerCount)
  {
    int tid = omp_get_thread_num();


    localcomp[tid] = true;

    //Get L2 block for this thread
    typename FILELoaderType::L2BlockType d2 = loader.getNextL2Block(tid);

    RangeType r2 = d2.getRange();
    ValueType* gold = new ValueType[this->kmerSize];
    while(r2.end - r2.start > 0)
    {
      //std::cout << "Thread# " << omp_get_thread_num() << ", range: " << r2.start << "," << r2.end << std::endl;

      //Begin iteration
      FASTAParser<typename FILELoaderType::L1BlockType::iterator, KmerType> FASTAIterObj(d2.begin(), d2.end(), r2, loader.getFileRange(), vectorReturned);

      //For crosschecking the contents
      for(auto iter=FASTAIterObj.begin(); iter != FASTAIterObj.end(); iter++)
      {
        std::string KmertoString = bliss::utils::KmerUtils::toASCIIString(FASTAIterObj.getKmer(iter));
        //std::cout << FASTAIterObj.getOffset(iter) << ", Thread#,Kmer# " << omp_get_thread_num() << "," << totalKmerCount+1 << " :  " << KmertoString << std::endl;
        auto offset = FASTAIterObj.getOffset(iter);

        //Fetch Kmer using simple file reading
        this->readFilePOSIX(this->fileName, offset, this->kmerSize, gold);

        //Compare the contents, by constructing a new kmer and equating
        std::string directFromFile((char *)gold, this->kmerSize);
        //std::cout << "Should match with " << " :  " << directFromFile<< std::endl;

        ++localKmerCount;

        localcomp[tid] &= (KmertoString == directFromFile);
      }

      //Get next L2 block
      d2 = loader.getNextL2Block(tid);
      r2 = d2.getRange();
    }
    delete [] gold;

#pragma omp critical
    comp &= localcomp[tid];


  } // end omp pragma

  this->kmerCount = localKmerCount;

  //Total kmers discovered should match with #kmers found in sequential test
  EXPECT_EQ(true, comp);

}
#endif
#endif

INSTANTIATE_TEST_CASE_P(Bliss, FASTAIteratorTest2, ::testing::Values(
    TestFileInfo(500, 14625, std::string("/test/data/natural.fasta")),
    TestFileInfo(246, 940, std::string("/test/data/test2.fasta")),
    TestFileInfo(335000, 1092580, std::string("/test/data/test.medium.fasta")),
    TestFileInfo(64, 512, std::string("/test/data/test.fasta"))
));
#endif

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
