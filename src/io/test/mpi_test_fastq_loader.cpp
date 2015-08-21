/**
 * fastqloader_test.cpp
 *
 *  Created on: Feb 18, 2014
 *      Author: Tony Pan <tpan7@gatech.edu>
 */


#include "bliss-config.hpp"    // for location of data.

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
      //ASSERT_EQ(29250, fileSize);
#if defined(USE_MPI)
  MPI_Barrier(MPI_COMM_WORLD);
#endif

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
TYPED_TEST_P(FASTQLoaderTest, MPIOpenWithRange)
{
  typedef FASTQLoader<TypeParam, true, false, false> FASTQLoaderType;

  // get this->fileName
  int rank = 0;
  int nprocs = 1;

#if defined(USE_MPI)
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

#endif

  FASTQLoaderType loader(nprocs, rank, this->fileName );

  auto l1 = loader.getNextL1Block();
  l1.getRange();


//  printf("Rank %d/%d file range: %lu, %lu, range [%lu, %lu)\n", rank, nprocs, loader.getFileRange().start, loader.getFileRange().end, r.start, r.end);


//  INFO( "range: " << r2 );
  ASSERT_EQ('@', loader.getCurrentL1Block().begin()[0]);
  if (rank == nprocs - 1)
    ASSERT_EQ('\0', loader.getCurrentL1Block().end()[0]);
  else
    ASSERT_EQ('@', loader.getCurrentL1Block().end()[0]);
//  INFO( " characters = '" << loader.getData().begin()[0]  << "'" );
//  INFO( " characters = '" << loader.getData().end()[0]  << "'" );

#if defined(USE_MPI)
  MPI_Barrier(MPI_COMM_WORLD);
#endif
}

// normal test cases
TYPED_TEST_P(FASTQLoaderTest, OpenWithRange)
{
  typedef FASTQLoader<TypeParam, true, false, false> FASTQLoaderType;

  // get this->fileName
  int rank = 3;
  int nprocs = 7;



  FASTQLoaderType loader(nprocs, rank, this->fileName );

  loader.getNextL1Block();



//  INFO( "range: " << r2 );
  ASSERT_EQ('@', loader.getCurrentL1Block().begin()[0]);
  ASSERT_EQ('@', loader.getCurrentL1Block().end()[0]);
//  INFO( " characters = '" << loader.getData().begin()[0]  << "'" );
//  INFO( " characters = '" << loader.getData().end()[0]  << "'" );
#if defined(USE_MPI)
  MPI_Barrier(MPI_COMM_WORLD);
#endif

}


TYPED_TEST_P(FASTQLoaderTest, OpenConsecutiveRanges)
{
  typedef FASTQLoader<TypeParam, true, false, false> FASTQLoaderType;
  typedef typename FASTQLoaderType::RangeType RangeType;

  // get this->fileName
  int nprocs = 7;
  typename RangeType::ValueType lastEnd = 0;

  RangeType r;
  typename FASTQLoaderType::L1BlockType d;

  for (int rank = 0; rank < nprocs; ++rank) {
    {
      FASTQLoaderType loader(nprocs, rank, this->fileName);

      d = loader.getNextL1Block();
      r = d.getRange();

      ASSERT_EQ(lastEnd, r.start);
      lastEnd = r.end;

      //  INFO( "range: " << r2 );
        ASSERT_EQ('@', d.begin()[0]);
        if (r.end == loader.getFileRange().end)
          ASSERT_EQ(0, d.end()[0]);
        else
          ASSERT_EQ('@', d.end()[0]);
      //  INFO( " characters = '" << loader.getData().begin()[0]  << "'" );
      //  INFO( " characters = '" << loader.getData().end()[0]  << "'" );
    }
  }
#if defined(USE_MPI)
  MPI_Barrier(MPI_COMM_WORLD);
#endif

}


// now register the test cases
REGISTER_TYPED_TEST_CASE_P(FASTQLoaderTest, MPIOpenWithRange, OpenWithRange, OpenConsecutiveRanges);





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
#if defined(USE_MPI)
  MPI_Barrier(MPI_COMM_WORLD);
#endif

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
  typedef typename FASTQLoaderType::L2BlockType blockType;

  // get this->fileName
  int rank = 3;
  int nprocs = 7;

  int nThreads = 2;


  FASTQLoaderType loader(nprocs, rank, this->fileName, nThreads, 2048);

  typename FASTQLoaderType::L1BlockType d = loader.getNextL1Block();

  typename RangeType::ValueType lastEnd = d.getRange().start;


  ValueType* gold;
  int i = 0;

  blockType data = loader.getNextL2Block(0);
  RangeType r2 = data.getRange();
  size_t len = r2.size();

  while (len  > 0) {


    ASSERT_EQ(lastEnd, r2.start);
    lastEnd = r2.end;


    gold = new ValueType[len + 1];
    FASTQLoaderBufferTest<TypeParam>::readFilePOSIX(this->fileName, r2.start, len, gold);

    ASSERT_TRUE(equal(gold, data.begin(), len));
    delete [] gold;

    i ^= 1;
    data = loader.getNextL2Block(i);
    r2 = data.getRange();
    len = r2.size();

  }
#if defined(USE_MPI)
  MPI_Barrier(MPI_COMM_WORLD);
#endif

}


// now register the test cases
REGISTER_TYPED_TEST_CASE_P(FASTQLoaderBufferTest, BufferingChunks);








typedef ::testing::Types<unsigned char> FASTQLoaderTestTypes;
INSTANTIATE_TYPED_TEST_CASE_P(Bliss, FASTQLoaderTest, FASTQLoaderTestTypes);


typedef ::testing::Types<
    bliss::io::FASTQLoader<unsigned char, true, false, false>,
    bliss::io::FASTQLoader<unsigned char, true, false, true>,
    bliss::io::FASTQLoader<unsigned char, true, true, false>,
    bliss::io::FASTQLoader<unsigned char, true, true, true>,
    bliss::io::FASTQLoader<unsigned char, false, false, false>,
    bliss::io::FASTQLoader<unsigned char, false, false, true>,
    bliss::io::FASTQLoader<unsigned char, false, true, false>,
    bliss::io::FASTQLoader<unsigned char, false, true, true> >
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
