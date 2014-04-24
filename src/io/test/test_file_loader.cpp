/**
 * file_loader_test.cpp
 *
 *  Created on: Feb 18, 2014
 *      Author: tpan
 */

// include google test
#include <gtest/gtest.h>
#include <cstdint> // for uint64_t, etc.
#include <string>

#include "config.hpp"    // for location of data.
#include "io/file_loader.hpp"

using namespace bliss::io;

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
      nElem = fileSize / sizeof(T);
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
    size_t nElem;

};

// indicate this is a typed test
TYPED_TEST_CASE_P(FileLoaderTest);


// normal test cases
TYPED_TEST_P(FileLoaderTest, OpenWithFullRange)
{
  typedef file_loader<TypeParam> FileLoaderType;
  typedef typename FileLoaderType::RangeType RangeType;

  // get fileName
  int rank = 0;
  int nprocs = 1;

  RangeType r =
      RangeType::block_partition(nprocs, rank, 0, this->nElem);

  FileLoaderType loader(this->fileName, r);

  size_t len = r.end - r.start;
  TypeParam* gold = new TypeParam[len + 1];
  FileLoaderTest<TypeParam>::readFilePOSIX(this->fileName, r.start, len, gold);

  int comp = memcmp(gold, loader.begin(), len * sizeof(TypeParam));
  ASSERT_EQ(0, comp);
  delete [] gold;
}

// normal test cases
TYPED_TEST_P(FileLoaderTest, PreloadWithFullRange)
{
  typedef file_loader<TypeParam> FileLoaderType;
  typedef typename FileLoaderType::RangeType RangeType;

  // get fileName
  int rank = 0;
  int nprocs = 1;

  RangeType r =
      RangeType::block_partition(nprocs, rank, 0, this->nElem);

  FileLoaderType loader(this->fileName, r, 0.1f);

  size_t len = r.end - r.start;
  TypeParam* gold = new TypeParam[len + 1];
  FileLoaderTest<TypeParam>::readFilePOSIX(this->fileName, r.start, len, gold);

  int comp = memcmp(gold, loader.begin(), len * sizeof(TypeParam));
  ASSERT_EQ(0, comp);
  delete [] gold;
}

// normal test cases
TYPED_TEST_P(FileLoaderTest, AttemptPreloadWithFullRange)
{
  typedef file_loader<TypeParam> FileLoaderType;
  typedef typename FileLoaderType::RangeType RangeType;

  // get fileName
  int rank = 0;
  int nprocs = 1;

  RangeType r =
      RangeType::block_partition(nprocs, rank, 0, this->nElem);

  FileLoaderType loader(this->fileName, r, 0.0001f);

  size_t len = r.end - r.start;
  TypeParam* gold = new TypeParam[len + 1];
  FileLoaderTest<TypeParam>::readFilePOSIX(this->fileName, r.start, len, gold);

  int comp = memcmp(gold, loader.begin(), len * sizeof(TypeParam));
  ASSERT_EQ(0, comp);
  delete [] gold;
}


TYPED_TEST_P(FileLoaderTest, OpenWithRange)
{
  typedef file_loader<TypeParam> FileLoaderType;
  typedef typename FileLoaderType::RangeType RangeType;

  // get this->fileName
  int rank = 3;
  int nprocs = 7;

  RangeType r =
      RangeType::block_partition(nprocs, rank, 0, this->nElem);

  FileLoaderType loader(this->fileName, r);

  size_t len = r.end - r.start;
  TypeParam* gold = new TypeParam[len + 1];
  FileLoaderTest<TypeParam>::readFilePOSIX(this->fileName, r.start, len, gold);

  int comp = memcmp(gold, loader.begin(), len * sizeof(TypeParam));
  ASSERT_EQ(0, comp);
  delete [] gold;
}

TYPED_TEST_P(FileLoaderTest, OpenWithAlignedRange)
{
  typedef file_loader<TypeParam> FileLoaderType;
  typedef typename FileLoaderType::RangeType RangeType;

  // get this->fileName
  int rank = 3;
  int nprocs = 7;

  RangeType r =
      RangeType::block_partition(nprocs, rank, 0, this->nElem);

  RangeType ra = r.align_to_page(
      sysconf(_SC_PAGE_SIZE));

  FileLoaderType loader(this->fileName, ra);

  RangeType r2 = loader.getRange();
  size_t len = r2.end - r2.start;

  TypeParam* gold = new TypeParam[len + 1];
  FileLoaderTest<TypeParam>::readFilePOSIX(this->fileName, r2.start, len, gold);

  int comp = memcmp(gold, loader.begin(), len * sizeof(TypeParam));
  ASSERT_EQ(0, comp);
  delete [] gold;
}


// normal test cases
TYPED_TEST_P(FileLoaderTest, PreloadWithRange)
{
  typedef file_loader<TypeParam> FileLoaderType;
  typedef typename FileLoaderType::RangeType RangeType;

  // get this->fileName
  int rank = 3;
  int nprocs = 7;

  RangeType r =
      RangeType::block_partition(nprocs, rank, 0, this->nElem);

  FileLoaderType loader(this->fileName, r, 0.1f);

  RangeType r2 = loader.getRange();
  size_t len = r2.end - r2.start;

  TypeParam* gold = new TypeParam[len + 1];
  FileLoaderTest<TypeParam>::readFilePOSIX(this->fileName, r2.start, len, gold);

  int comp = memcmp(gold, loader.begin(), len * sizeof(TypeParam));
  ASSERT_EQ(0, comp);
  delete [] gold;
}
// normal test cases
TYPED_TEST_P(FileLoaderTest, AttemptPreloadWithRange)
{
  typedef file_loader<TypeParam> FileLoaderType;
  typedef typename FileLoaderType::RangeType RangeType;

  // get this->fileName
  int rank = 3;
  int nprocs = 7;

  RangeType r =
      RangeType::block_partition(nprocs, rank, 0, this->nElem);

  FileLoaderType loader(this->fileName, r, 0.0001f);

  RangeType r2 = loader.getRange();
  size_t len = r2.end - r2.start;

  TypeParam* gold = new TypeParam[len + 1];
  FileLoaderTest<TypeParam>::readFilePOSIX(this->fileName, r2.start, len, gold);


  int comp = memcmp(gold, loader.begin(), len * sizeof(TypeParam));
  ASSERT_EQ(0, comp);
  delete [] gold;
}


// now register the test cases
REGISTER_TYPED_TEST_CASE_P(FileLoaderTest, OpenWithFullRange, PreloadWithFullRange, AttemptPreloadWithFullRange,
                           OpenWithRange, OpenWithAlignedRange, PreloadWithRange, AttemptPreloadWithRange);


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
  typedef file_loader<TypeParam> FileLoaderType;
  typedef typename FileLoaderType::RangeType RangeType;

  RangeType r;
  r.start = 0;
  r.end = 1;

  std::string fn;

  std::string err_regex = ".*file_loader.hpp.* Assertion .* failed.*";
  EXPECT_EXIT(FileLoaderType(fn, r), ::testing::KilledBySignal(SIGABRT), err_regex);

}

TYPED_TEST_P(FileLoaderDeathTest, BadFilename)
{
  typedef file_loader<TypeParam> FileLoaderType;
  typedef typename FileLoaderType::RangeType RangeType;

  RangeType r;
  r.start = 0;
  r.end = 1;
  std::string fn;
  fn.assign("bad");

  std::string err_regex = ".*file_loader.hpp.* Exception .* ERROR .* error 2: No such file or directory";
  EXPECT_THROW(FileLoaderType(fn, r), bliss::io::io_exception);
}


TYPED_TEST_P(FileLoaderDeathTest, EmptyRange)
{
  typedef file_loader<TypeParam> FileLoaderType;
  typedef typename FileLoaderType::RangeType RangeType;

  std::string fn;
  fn.assign(PROJ_SRC_DIR);
  fn.append("/test/data/empty.fastq");

  RangeType r;
  r.start = 0;
  r.end = 0;

  std::string err_regex = ".*file_loader.hpp.* Assertion .* failed.*";
  EXPECT_EXIT(FileLoaderType(fn, r), ::testing::KilledBySignal(SIGABRT), err_regex);
}

// now register the test cases
REGISTER_TYPED_TEST_CASE_P(FileLoaderDeathTest, NoFilename, BadFilename, EmptyRange);




typedef ::testing::Types<unsigned char, char, int16_t, uint16_t, int32_t, uint32_t,
    int64_t, uint64_t> FileLoaderTestTypes;
INSTANTIATE_TYPED_TEST_CASE_P(Bliss, FileLoaderTest, FileLoaderTestTypes);
INSTANTIATE_TYPED_TEST_CASE_P(Bliss, FileLoaderDeathTest, FileLoaderTestTypes);
