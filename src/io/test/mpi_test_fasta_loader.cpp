/**
 * fastaloader_test.cpp
 *
 *  Created on: Oct 30, 2014
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

#include "io/fasta_loader.hpp"
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
class FASTALoaderTest : public ::testing::Test
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
TYPED_TEST_CASE_P(FASTALoaderTest);


// normal test cases
TYPED_TEST_P(FASTALoaderTest, OpenWithRange)
{
  typedef FileLoader<TypeParam, false, false> FILELoaderType;
  typedef FASTALoader<typename FILELoaderType::InputIteratorType> FASTALoaderType;
  typedef typename FILELoaderType::RangeType RangeType;

  // get this->fileName
  int rank = 3;
  int nprocs = 7;

  FILELoaderType loader(nprocs, rank, this->fileName);
  FASTALoaderType fastaPreparse();

  RangeType r;
  typename FILELoaderType::L1BlockType d;
  d = loader.getNextL1Block();
  r = d.getRange();

  typename FASTALoaderType::vectorType vectorReturned;

  FASTALoaderType obj;
  obj.countSequenceStarts(d.begin(), loader.getFileRange() , r, vectorReturned);
  ASSERT_LT(0 , vectorReturned.size());

}

// now register the test cases
REGISTER_TYPED_TEST_CASE_P(FASTALoaderTest, OpenWithRange);


typedef ::testing::Types<unsigned char> FASTALoaderTestTypes;
INSTANTIATE_TYPED_TEST_CASE_P(Bliss, FASTALoaderTest, FASTALoaderTestTypes);


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
