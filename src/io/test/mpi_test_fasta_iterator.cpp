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
      fileName.append("/test/data/test.fasta");

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

};

// indicate this is a typed test
TYPED_TEST_CASE_P(FASTAIteratorTest);

/*
// normal test cases
TYPED_TEST_P(FASTAIteratorTest, FASTAIteration)
{
  typedef FileLoader<TypeParam, false, false> FILELoaderType;
  typedef FASTALoader<typename FILELoaderType::InputIteratorType> FASTALoaderType;
  typedef typename FILELoaderType::RangeType RangeType;

  // get this->fileName
  int rank = 3;
  int nprocs = 7;

  FILELoaderType loader(nprocs, rank, this->fileName);

  RangeType r;
  typename FILELoaderType::L1BlockType d;
  d = loader.getNextL1Block();
  r = d.getRange();

  typename FASTALoaderType::vectorType vectorReturned;

  FASTALoaderType obj;
  obj.countSequenceStarts(d.begin(), loader.getFileRange() , r, vectorReturned);
  ASSERT_LT(0 , vectorReturned.size());

  //Begin iteration
  typedef bliss::Kmer<35, 2, uint32_t> KmerType;
  FASTAParser<typename FILELoaderType::L1BlockType::iterator, KmerType> FASTAIterObj(d.begin(), d.end(), loader.getFileRange(), vectorReturned);
}
*/

// normal test cases
TYPED_TEST_P(FASTAIteratorTest, BufferingChunks)
{
  typedef FileLoader<TypeParam, false, false> FILELoaderType;
  typedef typename FILELoaderType::InputIteratorType InputIterType;
  typedef typename std::iterator_traits<InputIterType>::value_type ValueType;

  typedef FASTALoader<typename FILELoaderType::InputIteratorType> FASTALoaderType;
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

  //Begin iteration
  const int len = 35;
  typedef bliss::Kmer<len, 2, uint32_t> KmerType;
  FASTAParser<typename FILELoaderType::L1BlockType::iterator, KmerType> FASTAIterObj(d.begin(), d.end(), loader.getFileRange(), vectorReturned);
  ValueType* gold;
  gold = new ValueType[len];

  std::cout << "Test started" << std::endl;

  auto startIter = FASTAIterObj.begin();
  auto startoffset = FASTAIterObj.getOffset(startIter);

  printf("BL : %lu \n", startoffset);

  auto endIter = FASTAIterObj.end();

  auto iter=startIter;
  int i = 1;
  for(auto iter=startIter; iter != endIter; iter++)
  {
    std::cout << FASTAIterObj.getOffset(iter) << ", Kmer# " << i++ << " :  " << FASTAIterObj.getKmer(iter).toString() << std::endl; 
    /*
    auto offset = FASTAIterObj.getOffset(iter);
    auto thisKmer = FASTAIterObj.getKmer(iter);

    //Fetch Kmer using simple file reading 
    FASTAIteratorTest<TypeParam>::readFilePOSIX(this->fileName, offset, len, gold);

    //Compare the contents, by constructing a new kmer and equating
    KmerType otherKmer;  
    otherKmer.fillFromChars(gold);

    bool comp = (thisKmer == otherKmer);
    ASSERT_EQ(true, comp);
    std::cout << "Check successful" << std::endl;
    */
  }
  delete [] gold;
}

// now register the test cases
//REGISTER_TYPED_TEST_CASE_P(FASTAIteratorTest, FASTAIteration, BufferingChunks);
REGISTER_TYPED_TEST_CASE_P(FASTAIteratorTest, BufferingChunks);


typedef ::testing::Types<unsigned char> FASTAIteratorTestTypes;
INSTANTIATE_TYPED_TEST_CASE_P(Bliss, FASTAIteratorTest, FASTAIteratorTestTypes);

int main(int argc, char* argv[]) {
    int result = 0;

    ::testing::InitGoogleTest(&argc, argv);


#if defined(USE_MPI)
    MPI_Init(&argc, &argv);
#endif
    result = RUN_ALL_TESTS();
#if defined(USE_MPI)
    MPI_Finalize();
#endif
    return result;
}
