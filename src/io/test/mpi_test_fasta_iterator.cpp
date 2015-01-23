/**
 * fastaIterator_test.cpp
 *
 *  Created on: Nov 14, 2014
 *      Author: cjain
 */


#include "config.hpp"    // for location of data.

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
#include "common/Kmer.hpp"
#include "common/alphabets.hpp"
#include "io/fasta_iterator.hpp"

using namespace bliss::io;


template<typename Iter1, typename Iter2>
bool equal(const Iter1 &i1, const Iter2 &i2, size_t len) {
  Iter1 ii1(i1);
  Iter2 ii2(i2);
  for (size_t i = 0; i < len; ++i, ++ii1, ++ii2) {
    if (*ii1 != *ii2) return false;
  }
  return true;
}

template <typename T>
class FASTAIteratorTest : public ::testing::Test
{
  protected:
    virtual void SetUp()
    {
      fileName.assign(PROJ_SRC_DIR);
      fileName.append("/test/data/test.medium.fasta");

      // get file size
      struct stat filestat;
      stat(fileName.c_str(), &filestat);
      fileSize = static_cast<size_t>(filestat.st_size);

      //ASSERT_EQ(34526831, fileSize);
      //ASSERT_EQ(114, fileSize);
    }

    static void readFilePOSIX(const std::string &fileName, const size_t& offset,
                              const size_t& length, T* result)
    {
      FILE *fp = fopen(fileName.c_str(), "r");
      fseek(fp, offset * sizeof(T), SEEK_SET);
      size_t read = fread_unlocked(result, sizeof(T), length, fp);
      fclose(fp);

      ASSERT_GT(read, 0);
    }

    std::string fileName;
    size_t fileSize;      // in units of T.
    //This variable to add a check that the count of Kmers is equal 
    //irrespective of noprocs and nothreads.
    static int totalKmerCount_seq;

};

// indicate this is a typed test
TYPED_TEST_CASE_P(FASTAIteratorTest);

// Single thread, single process
TYPED_TEST_P(FASTAIteratorTest, BufferingChunks)
{
  typedef FileLoader<TypeParam, false, false> FILELoaderType;
  typedef typename FILELoaderType::InputIteratorType InputIterType;
  typedef typename std::iterator_traits<InputIterType>::value_type ValueType;

  //Define Kmer type
  const int len = 35;
  typedef bliss::Kmer<len, DNA5, uint32_t> KmerType;

  //Define type for FASTALoader 
  typedef FASTALoader<typename FILELoaderType::InputIteratorType, KmerType> FASTALoaderType;
  typename FASTALoaderType::vectorType vectorReturned;

  typedef typename FILELoaderType::RangeType RangeType;
  typename FILELoaderType::L1BlockType d;

  int rank = 0;
  int nprocs = 1;


  //Load the file in FILE Loader
  FILELoaderType loader(nprocs, rank, this->fileName);

  d = loader.getNextL1Block();
  RangeType r = d.getRange();

  //Get the vector of fasta headers from FASTA Loader
  FASTALoaderType obj;
  obj.countSequenceStarts(d.begin(), loader.getFileRange() , r, vectorReturned);
  ASSERT_LT(0 , vectorReturned.size());

#ifdef USE_MPI
  //This test is supposed to be run on a single core only
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if(rank==0)
  {
#endif

  //Begin iteration
  FASTAParser<typename FILELoaderType::L1BlockType::iterator, KmerType> FASTAIterObj(d.begin(), d.end(), r, loader.getFileRange(), vectorReturned);
  ValueType* gold;
  gold = new ValueType[len];

  int totalKmerCount = 0;
  //For crosschecking the contents
  for(auto iter=FASTAIterObj.begin(); iter != FASTAIterObj.end(); iter++)
  {
    std::string KmertoString = FASTAIterObj.getKmer(iter).toAlphabets();
    //std::cout << FASTAIterObj.getOffset(iter) << ", Kmer# " << totalKmerCount+1 << " :  " << KmertoString << std::endl; 
    auto offset = FASTAIterObj.getOffset(iter);

    //Fetch Kmer using simple file reading 
    FASTAIteratorTest<TypeParam>::readFilePOSIX(this->fileName, offset, len, gold);

    //Compare the contents, by constructing a new kmer and equating
    std::string directFromFile((char *)gold, len);
    //std::cout << "Should match with " << " :  " << directFromFile << std::endl; 
    totalKmerCount++;

    bool comp = (KmertoString == directFromFile);
    ASSERT_EQ(true, comp);
  }
  FASTAIteratorTest<TypeParam>::totalKmerCount_seq = totalKmerCount; 
  delete [] gold;

#ifdef USE_MPI
  }
#endif
}

#ifdef USE_MPI
// Single thread, multiple processes
TYPED_TEST_P(FASTAIteratorTest, BufferingChunks_MultiProc)
{
  typedef FileLoader<TypeParam, false, false> FILELoaderType;
  typedef typename FILELoaderType::InputIteratorType InputIterType;
  typedef typename std::iterator_traits<InputIterType>::value_type ValueType;

  //Define Kmer type
  const int len = 35;
  typedef bliss::Kmer<len, DNA5, uint32_t> KmerType;

  typedef FASTALoader<typename FILELoaderType::InputIteratorType, KmerType> FASTALoaderType;
  typename FASTALoaderType::vectorType vectorReturned;

  typedef typename FILELoaderType::RangeType RangeType;
  typename FILELoaderType::L1BlockType d;


  //Set the rank and number of processors using MPI
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  //Load the file in FILE Loader
  FILELoaderType loader(MPI_COMM_WORLD, this->fileName);

  d = loader.getNextL1Block();
  RangeType r = d.getRange();

  //Get the vector of fasta headers from FASTA Loader
  FASTALoaderType obj;
  obj.countSequenceStarts(d.begin(), loader.getFileRange() , r, vectorReturned);
  ASSERT_LT(0 , vectorReturned.size());

  //Begin iteration
  FASTAParser<typename FILELoaderType::L1BlockType::iterator, KmerType> FASTAIterObj(d.begin(), d.end(), r, loader.getFileRange(), vectorReturned);
  ValueType* gold;
  gold = new ValueType[len];

  int totalKmerCount = 0;
  //For crosschecking the contents
  for(auto iter=FASTAIterObj.begin(); iter != FASTAIterObj.end(); iter++)
  {
    std::string KmertoString = FASTAIterObj.getKmer(iter).toAlphabets();
    //std::cout << FASTAIterObj.getOffset(iter) << ", RANK#,Kmer# " << rank << "," << totalKmerCount+1 << " :  " << KmertoString << std::endl; 
    auto offset = FASTAIterObj.getOffset(iter);

    //Fetch Kmer using simple file reading 
    FASTAIteratorTest<TypeParam>::readFilePOSIX(this->fileName, offset, len, gold);

    //Compare the contents, by constructing a new kmer and equating
    std::string directFromFile((char *)gold, len);
    //std::cout << "Should match with " << " :  " << directFromFile << std::endl; 
    totalKmerCount++;

    bool comp = (KmertoString == directFromFile);
    ASSERT_EQ(true, comp);
  }

  //Reduce sum of Kmer count on root processor
  int totalKmerCountSUM;
  MPI_Reduce(&totalKmerCount, &totalKmerCountSUM, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

  //Total kmers discovered should match with #kmers found in sequential test
  if(rank == 0)
    ASSERT_EQ(totalKmerCountSUM, FASTAIteratorTest<TypeParam>::totalKmerCount_seq);

  delete [] gold;
}
#endif

// now register the test cases
//REGISTER_TYPED_TEST_CASE_P(FASTAIteratorTest, FASTAIteration, BufferingChunks);
#ifdef USE_MPI
REGISTER_TYPED_TEST_CASE_P(FASTAIteratorTest, BufferingChunks, BufferingChunks_MultiProc);
#endif
#ifndef USE_MPI
REGISTER_TYPED_TEST_CASE_P(FASTAIteratorTest, BufferingChunks);
#endif


typedef ::testing::Types<unsigned char> FASTAIteratorTestTypes;
INSTANTIATE_TYPED_TEST_CASE_P(Bliss, FASTAIteratorTest, FASTAIteratorTestTypes);

//Define static member of FASTAIteratorTest
template <typename T>int FASTAIteratorTest<T>::totalKmerCount_seq;

int main(int argc, char* argv[]) 
{
    int result = 0;

    ::testing::InitGoogleTest(&argc, argv);

    //int i = 0;
    //char hostname[256];
    //gethostname(hostname, sizeof(hostname));
    //printf("PID %d on %s ready for attach\n", getpid(), hostname);
    //fflush(stdout);
    //while (0 == i)
      //sleep(5);

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
