/**
 * fastqloader_test.cpp
 *
 *  Created on: Feb 18, 2014
 *      Author: tpan
 */

#include "io/fastq_partition_helper.hpp"

#include "config.hpp"    // for location of data.

#if defined(USE_MPI)
#include "mpi.h"
#endif

// include google test
#include <gtest/gtest.h>
#include <cstdint> // for uint64_t, etc.
#include <string>

#include "io/file_loader.hpp"

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
class FASTQPartitionHelperTest : public ::testing::Test
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
TYPED_TEST_CASE_P(FASTQPartitionHelperTest);




// normal test cases
TYPED_TEST_P(FASTQPartitionHelperTest, AdjustRange)
{
  typedef file_loader<TypeParam, false, false> FileLoaderType;
  typedef typename FileLoaderType::RangeType RangeType;

  // get this->fileName
  int rank = 3;
  int nprocs = 7;

  FileLoaderType loader(this->fileName);

  RangeType r = loader.getFileRange().block_partition(nprocs, rank);
  loader.setRange(r);
  loader.adjustRange(bliss::io::FASTQPartitionHelper<typename FileLoaderType::IteratorType, RangeType>());
  RangeType r2 = loader.getRange();
  loader.load();

//  std::cout << "range: " << r2 << std::endl;
  ASSERT_EQ('@', loader.getData().begin()[0]);
  ASSERT_EQ('@', loader.getData().end()[0]);
//  std::cout << " characters = '" << loader.getData().begin()[0]  << "'" << std::endl;
//  std::cout << " characters = '" << loader.getData().end()[0]  << "'" << std::endl;
  loader.unload();
}


TYPED_TEST_P(FASTQPartitionHelperTest, AdjustConsecutiveRanges)
{
  typedef file_loader<TypeParam, false, false> FileLoaderType;
  typedef typename FileLoaderType::RangeType RangeType;

  // get this->fileName
  int nprocs = 7;
  typename RangeType::ValueType lastEnd = 0;

  FileLoaderType loader(this->fileName);
  RangeType r, r2;
  bliss::io::FASTQPartitionHelper<typename FileLoaderType::IteratorType, RangeType> ip;

  for (int rank = 0; rank < nprocs; ++rank) {

    r = loader.getFileRange().block_partition(nprocs, rank);
    loader.setRange(r);
    loader.adjustRange(ip);
    r2 = loader.getRange();

    ASSERT_EQ(lastEnd, r2.start);
    lastEnd = r2.end;
  }
}

// normal test cases
TYPED_TEST_P(FASTQPartitionHelperTest, BufferChunks)
{
  typedef file_loader<TypeParam, true, false> FileLoaderType;
  typedef typename FileLoaderType::RangeType RangeType;

  // get this->fileName
  int rank = 3;
  int nprocs = 7;

  bliss::io::FASTQPartitionHelper<typename FileLoaderType::InputIteratorType, RangeType> ip;
  bliss::io::FASTQPartitionHelper<typename FileLoaderType::IteratorType,      RangeType> ip2;

  FileLoaderType loader(this->fileName);

  RangeType r = loader.getFileRange().block_partition(nprocs, rank);
  loader.setRange(r);
  loader.adjustRange(ip);
  loader.load();

  r = loader.getRange();
  typename RangeType::ValueType lastEnd = r.start;

  RangeType r2;

  size_t len;
  TypeParam* gold;


  typename FileLoaderType::DataBlockType data = loader.getNextChunkAtomic(ip2, 2048);

  r2 = data.getRange();
  len = r2.end - r2.start;

  while (len  > 0) {


    ASSERT_EQ(lastEnd, r2.start);
    lastEnd = r2.end;


    gold = new TypeParam[len + 1];
    FASTQPartitionHelperTest<TypeParam>::readFilePOSIX(this->fileName, r2.start, len, gold);

    ASSERT_TRUE(equal(gold, data.begin(), len));
    delete [] gold;

    data= loader.getNextChunkAtomic(ip2, 2048);
    r2 = data.getRange();
    len = r2.end - r2.start;

  }


  loader.unload();
}

// normal test cases
TYPED_TEST_P(FASTQPartitionHelperTest, UnbufferChunks)
{
  typedef file_loader<TypeParam, false, false> FileLoaderType;
  typedef typename FileLoaderType::RangeType RangeType;

  // get this->fileName
  int rank = 3;
  int nprocs = 7;

  bliss::io::FASTQPartitionHelper<typename FileLoaderType::InputIteratorType, RangeType> ip;
  bliss::io::FASTQPartitionHelper<typename FileLoaderType::IteratorType,      RangeType> ip2;

  FileLoaderType loader(this->fileName);

  RangeType r = loader.getFileRange().block_partition(nprocs, rank);
  loader.setRange(r);
  loader.adjustRange(ip);
  loader.load();

  r = loader.getRange();
  typename RangeType::ValueType lastEnd = r.start;

  RangeType r2;

  size_t len;
  TypeParam* gold;

  typename FileLoaderType::DataBlockType data = loader.getNextChunkAtomic(ip2, 2048);

  r2 = data.getRange();
  len = r2.end - r2.start;

  while (len  > 0) {


    ASSERT_EQ(lastEnd, r2.start);
    lastEnd = r2.end;


    gold = new TypeParam[len + 1];
    FASTQPartitionHelperTest<TypeParam>::readFilePOSIX(this->fileName, r2.start, len, gold);

    ASSERT_TRUE(equal(gold, data.begin(), len));
    delete [] gold;

    data= loader.getNextChunkAtomic(ip2, 2048);
    r2 = data.getRange();
    len = r2.end - r2.start;

  }

  loader.unload();
}

// normal test cases
TYPED_TEST_P(FASTQPartitionHelperTest, BufferChunksWithPreload)
{
  typedef file_loader<TypeParam, true, true> FileLoaderType;
  typedef typename FileLoaderType::RangeType RangeType;

  // get this->fileName
  int rank = 3;
  int nprocs = 7;

  bliss::io::FASTQPartitionHelper<typename FileLoaderType::InputIteratorType, RangeType> ip;
  bliss::io::FASTQPartitionHelper<typename FileLoaderType::IteratorType,      RangeType> ip2;

  FileLoaderType loader(this->fileName);

  RangeType r = loader.getFileRange().block_partition(nprocs, rank);
  loader.setRange(r);
  loader.adjustRange(ip);
  loader.load();

  r = loader.getRange();
  typename RangeType::ValueType lastEnd = r.start;

  RangeType r2;

  size_t len;
  TypeParam* gold;

  typename FileLoaderType::DataBlockType data = loader.getNextChunkAtomic(ip2, 2048);

  r2 = data.getRange();
  len = r2.end - r2.start;

  while (len  > 0) {


    ASSERT_EQ(lastEnd, r2.start);
    lastEnd = r2.end;


    gold = new TypeParam[len + 1];
    FASTQPartitionHelperTest<TypeParam>::readFilePOSIX(this->fileName, r2.start, len, gold);

    ASSERT_TRUE(equal(gold, data.begin(), len));
    delete [] gold;

    data= loader.getNextChunkAtomic(ip2, 2048);
    r2 = data.getRange();
    len = r2.end - r2.start;

  }

  loader.unload();
}

// normal test cases
TYPED_TEST_P(FASTQPartitionHelperTest, UnbufferChunksWithPreload)
{
  typedef file_loader<TypeParam, false, true> FileLoaderType;
  typedef typename FileLoaderType::RangeType RangeType;

  // get this->fileName
  int rank = 3;
  int nprocs = 7;

  bliss::io::FASTQPartitionHelper<typename FileLoaderType::InputIteratorType, RangeType> ip;
  bliss::io::FASTQPartitionHelper<typename FileLoaderType::IteratorType,      RangeType> ip2;

  FileLoaderType loader(this->fileName);

  RangeType r = loader.getFileRange().block_partition(nprocs, rank);
  loader.setRange(r);
  loader.adjustRange(ip);
  loader.load();

  r = loader.getRange();
  typename RangeType::ValueType lastEnd = r.start;

  RangeType r2;

  size_t len;
  TypeParam* gold;

  typename FileLoaderType::DataBlockType data = loader.getNextChunkAtomic(ip2, 2048);

  r2 = data.getRange();
  len = r2.end - r2.start;

  while (len  > 0) {


    ASSERT_EQ(lastEnd, r2.start);
    lastEnd = r2.end;


    gold = new TypeParam[len + 1];
    FASTQPartitionHelperTest<TypeParam>::readFilePOSIX(this->fileName, r2.start, len, gold);

    ASSERT_TRUE(equal(gold, data.begin(), len));
    delete [] gold;

    data= loader.getNextChunkAtomic(ip2, 2048);
    r2 = data.getRange();
    len = r2.end - r2.start;

  }

  loader.unload();
}


// now register the test cases
REGISTER_TYPED_TEST_CASE_P(FASTQPartitionHelperTest, AdjustRange, AdjustConsecutiveRanges,
                           BufferChunks, UnbufferChunks, BufferChunksWithPreload, UnbufferChunksWithPreload);


//template <typename T>
//class FileLoaderDeathTest : public ::testing::Test
//{
//  protected:
//    virtual void SetUp()
//    {
//    }
//
//};
//
//// indicate this is a typed test
//TYPED_TEST_CASE_P(FileLoaderDeathTest);
//
//
//// TODO negative test cases
//TYPED_TEST_P(FileLoaderDeathTest, NoFilename)
//{
//  typedef file_loader<TypeParam> FileLoaderType;
//
//  std::string fileName;
//
//  std::string err_regex = ".*file_loader.hpp.* Assertion .* failed.*";
//  EXPECT_EXIT(FileLoaderType(fileName), ::testing::KilledBySignal(SIGABRT), err_regex);
//
//}
//
//TYPED_TEST_P(FileLoaderDeathTest, BadFilename)
//{
//  typedef file_loader<TypeParam> FileLoaderType;
//
//  std::string fileName;
//  fileName.assign("bad");
//
//  std::string err_regex = ".*file_loader.hpp.* Exception .* ERROR .* error 2: No such file or directory";
//  EXPECT_THROW(FileLoaderType(fileName), bliss::io::io_exception);
//}
//
//
//// now register the test cases
//REGISTER_TYPED_TEST_CASE_P(FileLoaderDeathTest, NoFilename, BadFilename);




typedef ::testing::Types<unsigned char> FASTQPartitionHelperTestTypes;
INSTANTIATE_TYPED_TEST_CASE_P(Bliss, FASTQPartitionHelperTest, FASTQPartitionHelperTestTypes);
//INSTANTIATE_TYPED_TEST_CASE_P(Bliss, FileLoaderDeathTest, FASTQPartitionHelperTestTypes);


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
