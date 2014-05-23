/**
 * fastqloader_test.cpp
 *
 *  Created on: Feb 18, 2014
 *      Author: tpan
 */


#include "config.hpp"    // for location of data.

#if defined(USE_MPI)
#include "mpi.h"
#endif

// include google test
#include <gtest/gtest.h>
#include <cstdint> // for uint64_t, etc.
#include <string>

#include "io/fastq_loader.hpp"

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
class FASTQLoaderTest : public ::testing::Test
{
  protected:
    virtual void SetUp()
    {
      fileName.assign(PROJ_SRC_DIR);
      fileName.append("/test/data/test.fastq");

      // get file size
      struct stat filestat;
      stat(fileName.c_str(), &filestat);
      fileSize = static_cast<size_t>(filestat.st_size);

      ASSERT_EQ(34111308, fileSize);

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
TYPED_TEST_CASE_P(FASTQLoaderTest);




// normal test cases
TYPED_TEST_P(FASTQLoaderTest, OpenWithRange)
{
  typedef FASTQLoader<TypeParam, false, false> FASTQLoaderType;
  typedef typename FASTQLoaderType::RangeType RangeType;

  // get this->fileName
  int rank = 3;
  int nprocs = 7;

  FASTQLoaderType loader(this->fileName, nprocs, rank);

  RangeType r = loader.getNextPartitionRange(rank);

  loader.load(r);

//  std::cout << "range: " << r2 << std::endl;
  ASSERT_EQ('@', loader.getData().begin()[0]);
  ASSERT_EQ('@', loader.getData().end()[0]);
//  std::cout << " characters = '" << loader.getData().begin()[0]  << "'" << std::endl;
//  std::cout << " characters = '" << loader.getData().end()[0]  << "'" << std::endl;
  loader.unload();
}


TYPED_TEST_P(FASTQLoaderTest, OpenConsecutiveRanges)
{
  typedef FASTQLoader<TypeParam, false, false> FASTQLoaderType;
  typedef typename FASTQLoaderType::RangeType RangeType;

  // get this->fileName
  int nprocs = 7;
  typename RangeType::ValueType lastEnd = 0;

  FASTQLoaderType loader(this->fileName, nprocs);
  RangeType r;

  for (int rank = 0; rank < nprocs; ++rank) {

    r = loader.getNextPartitionRange(rank);
    loader.load(r);

    ASSERT_EQ(lastEnd, r.start);
    lastEnd = r.end;

    //  std::cout << "range: " << r2 << std::endl;
      ASSERT_EQ('@', loader.getData().begin()[0]);
      ASSERT_EQ('@', loader.getData().end()[0]);
    //  std::cout << " characters = '" << loader.getData().begin()[0]  << "'" << std::endl;
    //  std::cout << " characters = '" << loader.getData().end()[0]  << "'" << std::endl;
      loader.unload();
  }
}


// now register the test cases
REGISTER_TYPED_TEST_CASE_P(FASTQLoaderTest, OpenWithRange, OpenConsecutiveRanges);





template <typename Loader>
class FASTQLoaderBufferTest : public ::testing::Test
{
  protected:
    typedef typename Loader::InputIteratorType InputIterType;
    typedef typename std::iterator_traits<InputIterType>::value_type ValueType;


    virtual void SetUp()
    {
      fileName.assign(PROJ_SRC_DIR);
      fileName.append("/test/data/test.fastq");

      // get file size
      struct stat filestat;
      stat(fileName.c_str(), &filestat);
      fileSize = static_cast<size_t>(filestat.st_size);

      ASSERT_EQ(34111308, fileSize);

    }

    static void readFilePOSIX(const std::string &fileName, const size_t& offset,
                              const size_t& length, ValueType* result)
    {
      FILE *fp = fopen(fileName.c_str(), "r");
      fseek(fp, offset * sizeof(ValueType), SEEK_SET);
      size_t read = fread_unlocked(result, sizeof(ValueType), length, fp);
      fclose(fp);

      ASSERT_GT(read, 0);
    }

    std::string fileName;
    size_t fileSize;      // in units of T.

};

// indicate this is a typed test
TYPED_TEST_CASE_P(FASTQLoaderBufferTest);


// normal test cases
TYPED_TEST_P(FASTQLoaderBufferTest, BufferingChunks)
{
  typedef TypeParam FASTQLoaderType;
  typedef typename FASTQLoaderType::InputIteratorType InputIterType;
  typedef typename std::iterator_traits<InputIterType>::value_type ValueType;

  typedef typename FASTQLoaderType::RangeType RangeType;
  typedef typename FASTQLoaderType::DataBlockType blockType;

  // get this->fileName
  int rank = 3;
  int nprocs = 7;

  int nThreads = 2;


  FASTQLoaderType loader(this->fileName, nprocs, rank, nthreads, 2048);

  RangeType r = loader.getNextPartitionRange(rank);
  loader.load(r);

  typename RangeType::ValueType lastEnd = r.start;

  RangeType r2;

  TypeParam* gold;

  RangeType r2 = loader.getNextChunkRange(0);
  blockType data = loader.getChunk(0, r2);
  size_t len = r2.size();

  while (len  > 0) {


    ASSERT_EQ(lastEnd, r2.start);
    lastEnd = r2.end;


    gold = new TypeParam[len + 1];
    FASTQLoaderTest<TypeParam>::readFilePOSIX(this->fileName, r2.start, len, gold);

    ASSERT_TRUE(equal(gold, data.begin(), len));
    delete [] gold;

    r2 = loader.getNextChunkRange(0 ^ 1);
    data = loader.getChunk(0 ^ 1, r2);
    len = r2.size();

  }


  loader.unload();
}


// now register the test cases
REGISTER_TYPED_TEST_CASE_P(FASTQLoaderBufferTest, BufferingChunks);








typedef ::testing::Types<unsigned char> FASTQLoaderTestTypes;
INSTANTIATE_TYPED_TEST_CASE_P(Bliss, FASTQLoaderTest, FASTQLoaderTestTypes);


typedef ::testing::Types<bliss::io::FASTQLoader<unsigned char, false, false>,
    bliss::io::FASTQLoader<unsigned char, false, true>,
    bliss::io::FASTQLoader<unsigned char, true, false>,
    bliss::io::FASTQLoader<unsigned char, true, true> >
    FASTQLoaderBufferTestTypes;
INSTANTIATE_TYPED_TEST_CASE_P(Bliss, FASTQLoaderBufferTest, FASTQLoaderBufferTestTypes);



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
