/**
 * file_loader_test.cpp

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

#include "io/file_loader.hpp"

using namespace bliss::io;


template<typename Iter1, typename Iter2>
bool equal(const Iter1 &i1, const Iter2 &i2, size_t len, bool print = false) {
  Iter1 ii1(i1);
  Iter2 ii2(i2);
  for (size_t i = 0; i < len; ++i, ++ii1, ++ii2) {
    if (print && *ii1 != *ii2) std::cout << i << " ";

    if (*ii1 != *ii2) return false;
  }
  return true;
}


template <typename T>
class FileLoaderTest : public ::testing::Test
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

      //printf("filename %s, offset %lu, length %lu\n", fileName.c_str(), offset, length);

      ASSERT_GT(read, 0);
    }

    std::string fileName;
    size_t fileSize;      // in units of T.

};

// indicate this is a typed test
TYPED_TEST_CASE_P(FileLoaderTest);


// normal test cases
TYPED_TEST_P(FileLoaderTest, OpenWithFullRange)
{
  typedef FileLoader<TypeParam, false, false> FileLoaderType;
  typedef typename FileLoaderType::RangeType RangeType;

  // get fileName

  FileLoaderType loader(this->fileName);
  RangeType r = loader.getFileRange();

  size_t len = r.end - r.start;
  TypeParam* gold = new TypeParam[len + 1];
  FileLoaderTest<TypeParam>::readFilePOSIX(this->fileName, r.start, len, gold);

  loader.load(r);
  int comp = memcmp(gold, loader.getData().begin(), len * sizeof(TypeParam));
  ASSERT_EQ(0, comp);
  delete [] gold;
  loader.unload();
}

// normal test cases
TYPED_TEST_P(FileLoaderTest, PreloadWithFullRange)
{
  typedef FileLoader<TypeParam, false, true> FileLoaderType;
  typedef typename FileLoaderType::RangeType RangeType;

  // get fileName

  FileLoaderType loader(this->fileName);
  RangeType r = loader.getFileRange();

  size_t len = r.end - r.start;
  TypeParam* gold = new TypeParam[len + 1];
  FileLoaderTest<TypeParam>::readFilePOSIX(this->fileName, r.start, len, gold);

  loader.load(r);
  ASSERT_TRUE(equal(gold, loader.getData().begin(), len));
  delete [] gold;
  loader.unload();

}


TYPED_TEST_P(FileLoaderTest, OpenWithRange)
{
  typedef FileLoader<TypeParam, false, false> FileLoaderType;
  typedef typename FileLoaderType::RangeType RangeType;

  // get this->fileName
  int rank = 3;
  int nprocs = 7;

  FileLoaderType loader(this->fileName, nprocs, rank);

  RangeType r = loader.getNextPartitionRange(rank);


  size_t len = r.end - r.start;
  TypeParam* gold = new TypeParam[len + 1];
  FileLoaderTest<TypeParam>::readFilePOSIX(this->fileName, r.start, len, gold);

  loader.load(r);
  int comp = memcmp(gold, loader.getData().begin(), len * sizeof(TypeParam));
  ASSERT_EQ(0, comp);
  delete [] gold;
  loader.unload();
}

TYPED_TEST_P(FileLoaderTest, OpenWithAlignedRange)
{
  typedef FileLoader<TypeParam, false, false> FileLoaderType;
  typedef typename FileLoaderType::RangeType RangeType;

  // get this->fileName
  int rank = 3;
  int nprocs = 7;


  FileLoaderType loader(this->fileName, nprocs, rank);

  RangeType r = loader.getNextPartitionRange(rank);


  RangeType ra = r.align_to_page(sysconf(_SC_PAGE_SIZE));

  loader.load(ra);

  size_t len = ra.end - ra.start;
  TypeParam* gold = new TypeParam[len + 1];
  FileLoaderTest<TypeParam>::readFilePOSIX(this->fileName, ra.start, len, gold);
  int comp = memcmp(gold, loader.getData().begin(), len * sizeof(TypeParam));
  ASSERT_EQ(0, comp);

  len = r.end - r.start;
  FileLoaderTest<TypeParam>::readFilePOSIX(this->fileName, r.start, len, gold);
  comp = memcmp(gold, loader.getData().begin(), len * sizeof(TypeParam));
  ASSERT_EQ(0, comp);


  delete [] gold;
  loader.unload();
}


// normal test cases
TYPED_TEST_P(FileLoaderTest, PreloadWithRange)
{
  typedef FileLoader<TypeParam, false, true> FileLoaderType;
  typedef typename FileLoaderType::RangeType RangeType;

  // get this->fileName
  int rank = 3;
  int nprocs = 7;

  FileLoaderType loader(this->fileName, nprocs, rank);

  RangeType r = loader.getNextPartitionRange(rank);

  size_t len = r.end - r.start;
  TypeParam* gold = new TypeParam[len + 1];
  FileLoaderTest<TypeParam>::readFilePOSIX(this->fileName, r.start, len, gold);

  loader.load(r);
  ASSERT_TRUE(equal(gold, loader.getData().begin(), len));
  delete [] gold;
  loader.unload();
}


TYPED_TEST_P(FileLoaderTest, OpenConsecutiveRanges)
{
  typedef FileLoader<TypeParam, false, false> FileLoaderType;
  typedef typename FileLoaderType::RangeType RangeType;

  // get this->fileName
  int nprocs = 7;

  FileLoaderType loader(this->fileName, nprocs);
  RangeType r;


  for (int rank = 0; rank < nprocs; ++rank) {

    loader.resetPartitionRange();

    r = loader.getNextPartitionRange(rank);
    //std::cout << "partition for rank " << rank << ": " << r << std::endl;

    size_t len = r.end - r.start;
    TypeParam* gold = new TypeParam[len + 1];
    FileLoaderTest<TypeParam>::readFilePOSIX(this->fileName, r.start, len, gold);

    loader.load(r);
    int comp = memcmp(gold, loader.getData().begin(), len * sizeof(TypeParam));
    ASSERT_EQ(0, comp);
    delete [] gold;
    loader.unload();
  }
}

// now register the test cases
REGISTER_TYPED_TEST_CASE_P(FileLoaderTest, OpenWithFullRange, PreloadWithFullRange,
                           OpenWithRange, OpenWithAlignedRange, PreloadWithRange, OpenConsecutiveRanges);


template <typename Loader>
class FileLoaderBufferTest : public ::testing::Test
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
TYPED_TEST_CASE_P(FileLoaderBufferTest);


// normal test cases
TYPED_TEST_P(FileLoaderBufferTest, BufferingChunks)
{
  typedef TypeParam FileLoaderType;
  typedef typename FileLoaderType::InputIteratorType InputIterType;
  typedef typename std::iterator_traits<InputIterType>::value_type ValueType;

  typedef typename FileLoaderType::RangeType RangeType;
  typedef typename FileLoaderType::DataBlockType blockType;

  // get this->fileName
  int rank = 3;
  int nprocs = 7;

  int nThreads = 2;

  FileLoaderType loader(this->fileName, nprocs, rank, nThreads, 2048);

  RangeType r = loader.getNextPartitionRange(rank);


  loader.load(r);
  typename RangeType::ValueType lastEnd = r.start;

  ValueType* gold;

  RangeType r2 = loader.getNextChunkRange(0);
  blockType data = loader.getChunk(0, r2);
  size_t len = r2.size();

  while (len  > 0) {

    ASSERT_EQ(lastEnd, r2.start);
    lastEnd = r2.end;

    gold = new ValueType[len + 1];
    FileLoaderBufferTest<TypeParam>::readFilePOSIX(this->fileName, r2.start, len, gold);

    ASSERT_TRUE(equal(gold, data.begin(), len));
    delete [] gold;

    r2 = loader.getNextChunkRange(0 ^ 1);
    data = loader.getChunk(0 ^ 1, r2);
    len = r2.size();
  }

  loader.unload();
}



// now register the test cases
REGISTER_TYPED_TEST_CASE_P(FileLoaderBufferTest, BufferingChunks);


template <typename T>
class FileLoaderDeathTest : public ::testing::Test
{
  protected:
    virtual void SetUp()
    {
    }

};

// indicate this is a typed test
TYPED_TEST_CASE_P(FileLoaderDeathTest);


// TODO negative test cases
TYPED_TEST_P(FileLoaderDeathTest, NoFilename)
{
  typedef FileLoader<TypeParam> FileLoaderType;

  std::string fileName;

  std::string err_regex = ".*file_loader.hpp.* Assertion .* failed.*";
  EXPECT_EXIT(new FileLoaderType(fileName), ::testing::KilledBySignal(SIGABRT), err_regex);

}

TYPED_TEST_P(FileLoaderDeathTest, BadFilename)
{
  typedef FileLoader<TypeParam> FileLoaderType;

  std::string fileName;
  fileName.assign("bad");

  std::string err_regex = ".*file_loader.hpp.* Exception .* ERROR .* error 2: No such file or directory";
  ASSERT_THROW(new FileLoaderType(fileName), bliss::io::IOException);
}

TYPED_TEST_P(FileLoaderDeathTest, BadNProcs)
{
  typedef FileLoader<TypeParam> FileLoaderType;

  std::string fileName;
  fileName.append("/test/data/test.fastq");

//  std::string err_regex = ".*file_loader.hpp.* Exception .* ERROR .* error 2: No such file or directory";
  std::string err_regex = ".*file_loader.hpp.* Assertion .* failed.*";
  ASSERT_EXIT(new FileLoaderType(fileName, 0), ::testing::KilledBySignal(SIGABRT), err_regex);
}

TYPED_TEST_P(FileLoaderDeathTest, BadRank)
{
  typedef FileLoader<TypeParam> FileLoaderType;

  std::string fileName;
  fileName.append("/test/data/test.fastq");

  std::string err_regex = ".*file_loader.hpp.* Assertion .* failed.*";
  ASSERT_EXIT(new FileLoaderType(fileName, 1, -1), ::testing::KilledBySignal(SIGABRT), err_regex);
}

TYPED_TEST_P(FileLoaderDeathTest, BadNThreads)
{
  typedef FileLoader<TypeParam> FileLoaderType;

  std::string fileName;
  fileName.append("/test/data/test.fastq");

  std::string err_regex = ".*file_loader.hpp.* Assertion .* failed.*";
  ASSERT_EXIT(new FileLoaderType(fileName, 1, 0, 0), ::testing::KilledBySignal(SIGABRT), err_regex);
}

TYPED_TEST_P(FileLoaderDeathTest, BadChunkSize)
{
  typedef FileLoader<TypeParam> FileLoaderType;

  std::string fileName;
  fileName.append("/test/data/test.fastq");

  std::string err_regex = ".*file_loader.hpp.* Assertion .* failed.*";
  ASSERT_EXIT(new FileLoaderType(fileName, 1, 0, 0, 0), ::testing::KilledBySignal(SIGABRT), err_regex);
}



// now register the test cases
REGISTER_TYPED_TEST_CASE_P(FileLoaderDeathTest, NoFilename, BadFilename, BadNProcs, BadRank, BadNThreads, BadChunkSize);




typedef ::testing::Types<unsigned char, char, int16_t, uint16_t, int32_t, uint32_t, int64_t, uint64_t>
    FileLoaderTestTypes;
INSTANTIATE_TYPED_TEST_CASE_P(Bliss, FileLoaderTest, FileLoaderTestTypes);
INSTANTIATE_TYPED_TEST_CASE_P(Bliss, FileLoaderDeathTest, FileLoaderTestTypes);


typedef ::testing::Types<bliss::io::FileLoader<unsigned char, false, false>,
    bliss::io::FileLoader<unsigned char, false, true>,
    bliss::io::FileLoader<unsigned char, true, false>,
    bliss::io::FileLoader<unsigned char, true, true> >
    FileLoaderBufferTestTypes;
INSTANTIATE_TYPED_TEST_CASE_P(Bliss, FileLoaderBufferTest, FileLoaderBufferTestTypes);

