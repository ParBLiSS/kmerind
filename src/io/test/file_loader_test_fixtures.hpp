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

#include "io/file_loader.hpp"

#ifdef USE_MPI
#include "mpi.h"
#include "io/mxx_support.hpp"
#include <mxx/collective.hpp>
#include <mxx/big_collective.hpp>
#endif

struct TestFileInfo {
    size_t elemCount;
    size_t fileSize;
    std::string filename;

    TestFileInfo(size_t const & _elem_count, size_t const & _file_size, std::string const & _file_name) :
      elemCount(_elem_count), fileSize(_file_size), filename(_file_name) {};
};


template<typename Iter1, typename Iter2>
bool equal(const Iter1 &i1, const Iter2 &i2, size_t len, bool print = false) {
  Iter1 ii1(i1);
  Iter2 ii2(i2);

  for (size_t i = 0; i < len; ++i, ++ii1, ++ii2) {
    if (*ii1 != *ii2) {
      if (print) std::cout << i << " " << *ii1 << "!=" << *ii2 << std::endl;
      return false;
    }
  }
  return true;
}



/**
 * @brief test for FileLoader, with various template parameters controlling the behavior (caching, overlapping, etc)
 */
template <typename Loader>
class FileLoadTypeParamTest : public ::testing::Test
{
  protected:
    typedef typename Loader::InputIteratorType                        InputIterType;
    typedef typename std::iterator_traits<InputIterType>::value_type  ValueType;

    virtual ~FileLoadTypeParamTest() {};

    virtual void readFilePOSIX(const std::string &fileName, const size_t offset,
                              const size_t length, ValueType* result)
    {
      FILE *fp = fopen(fileName.c_str(), "r");
      fseek(fp, offset * sizeof(ValueType), SEEK_SET);
      size_t read = fread_unlocked(result, sizeof(ValueType), length, fp);
      fclose(fp);

      ASSERT_GT(read, 0UL);
    }

    std::string fileName;
};

/**
 * @brief test fixture for file parsers, which segments file into sequences.
 */
template <typename Loader>
class FileLoadValueParamTest : public ::testing::TestWithParam<TestFileInfo>
{
  protected:
    typedef typename Loader::InputIteratorType                        InputIterType;
    typedef typename std::iterator_traits<InputIterType>::value_type  ValueType;

    virtual ~FileLoadValueParamTest() {};

    virtual void SetUp()
    {
      TestFileInfo const & p = GetParam();

      fileName.assign(PROJ_SRC_DIR);
      fileName.append(p.filename);

      // get file size
      struct stat filestat;
      stat(fileName.c_str(), &filestat);
      size_t fileSize = static_cast<size_t>(filestat.st_size);

      if (p.fileSize != fileSize) printf("filename : %s\n", fileName.c_str());

      ASSERT_EQ(p.fileSize, fileSize);
    }


    virtual void readFilePOSIX(const std::string &fileName, const size_t offset,
                              const size_t length, ValueType* result)
    {
      FILE *fp = fopen(fileName.c_str(), "r");
      fseek(fp, offset * sizeof(ValueType), SEEK_SET);
      size_t read = fread_unlocked(result, sizeof(ValueType), length, fp);
      fclose(fp);

      ASSERT_GT(read, 0UL);
    }

    std::string fileName;

};


/**
 * @brief test fixture for file parsers, which segments file into sequences.
 */
template <typename Loader>
class FileParserTest : public FileLoadValueParamTest<Loader>
{
  protected:
    typedef typename Loader::InputIteratorType                        InputIterType;
    typedef typename std::iterator_traits<InputIterType>::value_type  ValueType;

    virtual ~FileParserTest() {};

    virtual void SetUp()
    {
    	FileLoadValueParamTest<Loader>::SetUp();

      elemCount = 0;

#if defined(USE_MPI)
      MPI_Barrier(MPI_COMM_WORLD);
#endif

    }

    virtual void TearDown() {

      if (this->elemCount < std::numeric_limits<size_t>::max() ) {

        TestFileInfo const & p = this->GetParam();

#ifdef USE_MPI
        MPI_Barrier(MPI_COMM_WORLD);
        size_t totalElemCount = mxx::allreduce(elemCount);
#else
        size_t totalElemCount = elemCount;
#endif
        EXPECT_EQ(p.elemCount, totalElemCount);
      }
    }

    size_t elemCount;
};

/**
 * @brief Kmer Reader from file.   posix reader skips non sequence characters such as EOL.
 */
template <typename Loader>
class KmerReaderTest : public FileParserTest<Loader>
{
  protected:
    typedef typename Loader::InputIteratorType                        InputIterType;
    typedef typename std::iterator_traits<InputIterType>::value_type  ValueType;

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
      fseek(fp, offset * sizeof(ValueType), SEEK_SET);
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

    static constexpr size_t get_kmer_size() { return Loader::get_overlap_size(); }

};


/**
 * @brief test for FileLoader, with various template parameters controlling the behavior (caching, overlapping, etc)
 */
template <typename Loader>
class FileLoaderTest : public FileLoadTypeParamTest<Loader>
{
  protected:
    typedef typename Loader::InputIteratorType                        InputIterType;
    typedef typename std::iterator_traits<InputIterType>::value_type  ValueType;

    virtual ~FileLoaderTest() {};

    virtual void SetUp()
    {
      this->fileName.assign(PROJ_SRC_DIR);
      this->fileName.append("/test/data/test.fastq");  // does not matter which file to use.

      // get file size
      struct stat filestat;
      stat(this->fileName.c_str(), &filestat);
      size_t fileSize = static_cast<size_t>(filestat.st_size);

      ASSERT_EQ(34111308UL, fileSize);

    }

};




#endif /* SRC_IO_TEST_FILE_LOADER_TEST_FIXTURES_HPP_ */
