/*
 * Copyright 2015 Georgia Institute of Technology
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/**
 * @file    file_loader_test_fixtures.hpp
 * @ingroup
 * @author  tpan
 * @brief
 * @details
 *

 */
#ifndef SRC_IO_TEST_FILE_LOADER_TEST_FIXTURES_HPP_
#define SRC_IO_TEST_FILE_LOADER_TEST_FIXTURES_HPP_



#include "bliss-config.hpp"    // for location of data.

// include google test
#include <gtest/gtest.h>
#include <cstdint> // for uint64_t, etc.
#include <string>

#ifdef USE_MPI
#include "mxx/env.hpp"
#include "io/mxx_support.hpp"
#include <mxx/collective.hpp>
#include <mxx/big_collective.hpp>
#endif


struct TestFileInfo {
    size_t seqCount;
    size_t kmerCount;
    size_t fileSize;
    std::string filename;

    TestFileInfo(size_t const & _elem_count, size_t const& _kmer_count, size_t const & _file_size, std::string const & _file_name) :
      seqCount(_elem_count), kmerCount(_kmer_count), fileSize(_file_size), filename(_file_name) {};
};


template<typename Iter1, typename Iter2>
bool equal(const Iter1 &i1, const Iter2 &i2, size_t len, bool print = false) {
  Iter1 ii1(i1);
  Iter2 ii2(i2);

  for (size_t i = 0; i < len; ++i, ++ii1, ++ii2) {
    if (*ii1 != *ii2) {
      if (print) std::cout << i << ": \"" << *ii1 << "\"!=\"" << *ii2 << "\"" << std::endl;
      return false;
    }
  }
  return true;
}



/**
 * @brief test for FileLoader, with various template parameters controlling the behavior (caching, overlapping, etc)
 */
class FileLoadTypeParamTest : public ::testing::Test
{
  protected:
    using ValueType = unsigned char;
    using InputIterType = ValueType*;

    virtual ~FileLoadTypeParamTest() {};

    virtual void readFilePOSIX(const std::string &fileName, const size_t offset,
                              const size_t length, ValueType* result)
    {
//    	std::cout << "open " << fileName << " offset " << offset << " length " << length << " address " << static_cast<void*>(result) << std::endl;

      FILE *fp = fopen(fileName.c_str(), "r");

      int res = fseek(fp, offset * sizeof(ValueType), SEEK_SET);

      if (res == -1) {
        std::stringstream ss;
        int myerr = errno;
        ss << "ERROR in file seek: " << myerr << ": " << strerror(myerr);
        throw ::std::logic_error(ss.str());
      }

      size_t read = fread_unlocked(result, sizeof(ValueType), length, fp);
      fclose(fp);

      ASSERT_GT(read, 0UL);
    }

    std::string fileName;
};

/**
 * @brief test fixture for file parsers, which segments file into sequences.
 */
class FileLoadValueParamTest : public ::testing::TestWithParam<TestFileInfo>
{
  protected:
    using ValueType = unsigned char;
    using InputIterType = ValueType*;

    virtual ~FileLoadValueParamTest() {};

    virtual void SetUp()
    {
      TestFileInfo const & p = GetParam();

      fileName.assign(PROJ_SRC_DIR);
      fileName.append(p.filename);

      // get file size
      struct stat filestat;
      int result = stat(fileName.c_str(), &filestat);

      if (result == -1) {
        std::stringstream ss;
        int myerr = errno;
        ss << "ERROR in file stat: " << myerr << ": " << strerror(myerr);
        throw ::std::logic_error(ss.str());
      }


      size_t fileSize = static_cast<size_t>(filestat.st_size);

      if (p.fileSize != fileSize) printf("filename : %s\n", fileName.c_str());

      ASSERT_EQ(p.fileSize, fileSize);
    }


    virtual void readFilePOSIX(const std::string &fileName, const size_t offset,
                              const size_t length, ValueType* result)
    {
      FILE *fp = fopen(fileName.c_str(), "r");

      int res = fseek(fp, offset * sizeof(ValueType), SEEK_SET);

      if (res == -1) {
        std::stringstream ss;
        int myerr = errno;
        ss << "ERROR in file seek: " << myerr << ": " << strerror(myerr);
        throw ::std::logic_error(ss.str());
      }

      size_t read = fread_unlocked(result, sizeof(ValueType), length, fp);
      fclose(fp);

      ASSERT_GT(read, 0UL);
    }

    std::string fileName;

};


/**
 * @brief test fixture for file parsers, which segments file into sequences.
 */
class FileParserTest : public FileLoadValueParamTest
{
  protected:
    using ValueType = typename FileLoadValueParamTest::ValueType;
    using InputIterType = typename FileLoadValueParamTest::InputIterType;

    virtual ~FileParserTest() {};

    virtual void SetUp()
    {
    	FileLoadValueParamTest::SetUp();

      seqCount = 0;
      kmerCount = 0;

#if defined(USE_MPI)
      MPI_Barrier(MPI_COMM_WORLD);
#endif

    }

    virtual void TearDown() {

      if (this->seqCount < std::numeric_limits<size_t>::max() ) {

        TestFileInfo const & p = this->GetParam();

#ifdef USE_MPI
        MPI_Barrier(MPI_COMM_WORLD);
        size_t totalseqCount = mxx::allreduce(seqCount);
#else
        size_t totalseqCount = seqCount;
#endif
        EXPECT_EQ(p.seqCount, totalseqCount);
      }

      if (this->kmerCount < std::numeric_limits<size_t>::max() ) {

        TestFileInfo const & p = this->GetParam();

#ifdef USE_MPI
        MPI_Barrier(MPI_COMM_WORLD);
        size_t totalKmerCount = mxx::allreduce(kmerCount);
//        std::cout << "rank " << mxx::comm().rank() << " local kmer count = " << kmerCount << " total " << totalKmerCount << std::endl;

#else
        size_t totalKmerCount = kmerCount;
#endif
        EXPECT_EQ(p.kmerCount, totalKmerCount);
      }
    }

    size_t seqCount;
    size_t kmerCount;
};

/**
 * @brief Kmer Reader from file.   posix reader skips non sequence characters such as EOL.
 */
class KmerReaderTest : public FileParserTest
{
  protected:
    using ValueType = typename FileParserTest::ValueType;
    using InputIterType = typename FileParserTest::InputIterType;

    virtual ~KmerReaderTest() {};

    /**
     * @param length_hint   size of range in the file to read.  actual returned results may be smaller
     *
     */
    virtual void readFilePOSIX(const std::string &fileName, const size_t offset,
                               const size_t length_hint, ValueType* result)
    {
      // length is the target reading region.
      // actual output contains only ASCII alphabetic values.

      if (length_hint <= 0) return;

      FILE *fp = fopen(fileName.c_str(), "r");


      int res = fseek(fp, offset * sizeof(ValueType), SEEK_SET);

      if (res == -1) {
        std::stringstream ss;
        int myerr = errno;
        ss << "ERROR in file seek: " << myerr << ": " << strerror(myerr);
        throw ::std::logic_error(ss.str());
      }


      ValueType c;
      ValueType* curr = result;
      int read = 0;
      size_t total_read = 0;
      for (size_t i = 0; i < length_hint; ++i) {
        // each time we read, the position moves forward 1.
        read = fread_unlocked(&c, sizeof(ValueType), 1, fp);
        if (read > 0) {  // only if read something.
          total_read += read;
          if ((c >= 64 && c <= 90) || (c >= 97 && c <=122)) {
            // but only if we read a alphabetic character.
            *curr = std::toupper(c);
            ++curr;
          }
        } else {
          break;
        }
      }
      fclose(fp);
      EXPECT_GT(total_read, 0UL);

    }

};


/**
 * @brief test for FileLoader, with various template parameters controlling the behavior (caching, overlapping, etc)
 */
class FileLoaderTest : public FileLoadTypeParamTest
{
  protected:
    using ValueType = typename FileLoadTypeParamTest::ValueType;
    using InputIterType = typename FileLoadTypeParamTest::InputIterType;

    virtual ~FileLoaderTest() {};

    virtual void SetUp()
    {
      this->fileName.assign(PROJ_SRC_DIR);
      this->fileName.append("/test/data/test.fastq");  // does not matter which file to use.

      // get file size
      struct stat filestat;

      int result = stat(this->fileName.c_str(), &filestat);

      if (result == -1) {
        std::stringstream ss;
        int myerr = errno;
        ss << "ERROR in file stat: " << myerr << ": " << strerror(myerr);
        throw ::std::logic_error(ss.str());
      }

      size_t fileSize = static_cast<size_t>(filestat.st_size);

      ASSERT_EQ(34111308UL, fileSize);

    }

};




#endif /* SRC_IO_TEST_FILE_LOADER_TEST_FIXTURES_HPP_ */
