/**
 * fastaloader_test.cpp
 *
 *  Created on: Oct 30, 2014
 *      Author: cjain
 */


#include "bliss-config.hpp"    // for location of data.

#if defined(USE_MPI)
#include "mpi.h"
#endif

// include google test
#include <gtest/gtest.h>
#include <cstdint> // for uint64_t, etc.
#include <string>

#include "io/fasta_loader.hpp"
#include "io/file_loader.hpp"
#include "common/kmer.hpp"
#include "common/alphabets.hpp"

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
      fileName.append("/test/data/test2.fasta");

      // get file size
      struct stat filestat;
      stat(fileName.c_str(), &filestat);
      fileSize = static_cast<size_t>(filestat.st_size);

//      ASSERT_EQ(34526831, fileSize);
      //ASSERT_EQ(512, fileSize);
//      ASSERT_EQ(940, fileSize);

//      ASSERT_EQ(14625, fileSize);
      MPI_Barrier(MPI_COMM_WORLD);

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

//
//// normal test cases
//TYPED_TEST_P(FASTALoaderTest, OpenWithRange)
//{
//  typedef FileLoader<TypeParam, false, false> FILELoaderType;
//
//  //Define Kmer type
//  const int len = 35;
//  typedef typename bliss::common::Kmer<len, bliss::common::DNA5, uint32_t> KmerType;
//
//  typedef FASTALoader<typename FILELoaderType::InputIteratorType, KmerType> FASTALoaderType;
//  typedef typename FILELoaderType::RangeType RangeType;
//
//  int rank = 3;
//  int nprocs = 7;
//
//  FILELoaderType loader(nprocs, rank, this->fileName);
//
//  RangeType r;
//  typename FILELoaderType::L1BlockType d;
//  d = loader.getNextL1Block();
//  r = d.getRange();
//
//  typename FASTALoaderType::vectorType vectorReturned;
//
//  FASTALoaderType obj(MPI_COMM_WORLD);
//  obj.countSequenceStarts(d.begin(), loader.getFileRange() , r, vectorReturned);
//
//
//
//  for (auto e : vectorReturned) {
//    std::cout << "Rank " << rank << "/" << nprocs << " element: " << e.first << ", " << e.second << std::endl;
//  }
//
//
//  std::cout << "Rank "  << rank << "/" << nprocs << " file range: " << loader.getFileRange() << ", range " << r << ", vector returned size " << vectorReturned.size() << std::endl;
//
//  ASSERT_LT(0 , vectorReturned.size());
//
//  MPI_Barrier(MPI_COMM_WORLD);
//
//}
//
//
//// normal test cases
//TYPED_TEST_P(FASTALoaderTest, OpenWithRange2)
//{
//  typedef FileLoader<TypeParam, false, false> FILELoaderType;
//
//  //Define Kmer type
//  const int len = 35;
//  typedef typename bliss::common::Kmer<len, bliss::common::DNA5, uint32_t> KmerType;
//
//  typedef FASTALoader<typename FILELoaderType::InputIteratorType, KmerType> FASTALoaderType;
//  typedef typename FILELoaderType::RangeType RangeType;
//
//  int rank = 3;
//  int nprocs = 7;
//
//  FILELoaderType loader(nprocs, rank, this->fileName);
//
//  RangeType r;
//  typename FILELoaderType::L1BlockType d;
//  d = loader.getNextL1Block();
//  r = d.getRange();
//
//  typename FASTALoaderType::vectorType vectorReturned;
//
//  FASTALoaderType obj(MPI_COMM_WORLD);
//  obj.countSequenceStarts2(d.begin(), loader.getFileRange() , r, vectorReturned);
//
//
//
//  for (auto e : vectorReturned) {
//    std::cout << "Rank " << rank << "/" << nprocs << " element: " << e.first << ", " << e.second << std::endl;
//  }
//
//
//  std::cout << "Rank "  << rank << "/" << nprocs << " file range: " << loader.getFileRange() << ", range " << r << ", vector returned size " << vectorReturned.size() << std::endl;
//
//  ASSERT_LT(0 , vectorReturned.size());
//
//  MPI_Barrier(MPI_COMM_WORLD);
//
//}
//


// normal test cases
TYPED_TEST_P(FASTALoaderTest, MPIOpenWithRange)
{
  typedef FileLoader<TypeParam, false, false> FILELoaderType;

  //Define Kmer type
  const int len = 35;
  typedef typename bliss::common::Kmer<len, bliss::common::DNA5, uint32_t> KmerType;

  typedef FASTALoader<typename FILELoaderType::InputIteratorType, KmerType> FASTALoaderType;
  typedef typename FILELoaderType::RangeType RangeType;

  int rank = 0;
  int nprocs = 1;

#if defined(USE_MPI)
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

#endif

  FILELoaderType loader(nprocs, rank, this->fileName);

  RangeType r;
  typename FILELoaderType::L1BlockType d;
  d = loader.getNextL1Block();
  r = d.getRange();

  typename FASTALoaderType::vectorType vectorReturned;

  FASTALoaderType obj(MPI_COMM_WORLD);
  obj.countSequenceStarts(d.begin(), loader.getFileRange() , r, vectorReturned);



  for (auto e : vectorReturned) {
    std::cout << "Rank " << rank << "/" << nprocs << " element: " << e.first << ", " << e.second << std::endl;
  }


  std::cout << "Rank "  << rank << "/" << nprocs << " file range: " << loader.getFileRange() << ", range " << r << ", vector returned size " << vectorReturned.size() << std::endl;

  ASSERT_LT(0 , vectorReturned.size());

  MPI_Barrier(MPI_COMM_WORLD);

}



// normal test cases
TYPED_TEST_P(FASTALoaderTest, MPIOpenWithRange2)
{
  typedef FileLoader<TypeParam, false, false> FILELoaderType;

  //Define Kmer type
  const int len = 35;
  typedef typename bliss::common::Kmer<len, bliss::common::DNA5, uint32_t> KmerType;

  typedef FASTALoader<typename FILELoaderType::InputIteratorType, KmerType> FASTALoaderType;
  typedef typename FILELoaderType::RangeType RangeType;

  int rank = 0;
  int nprocs = 1;

#if defined(USE_MPI)
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

#endif

  FILELoaderType loader(nprocs, rank, this->fileName);

  RangeType r;
  typename FILELoaderType::L1BlockType d = loader.getNextL1Block();
  r = d.getRange();

  typename FASTALoaderType::vectorType vectorReturned;

  FASTALoaderType obj(MPI_COMM_WORLD);
  obj.countSequenceStarts2(d.begin(), loader.getFileRange() , r, vectorReturned);



  for (auto e : vectorReturned) {
    std::cout << "Rank " << rank << "/" << nprocs << " element: " << e.first << ", " << e.second << std::endl;
  }


  std::cout << "Rank "  << rank << "/" << nprocs << " file range: " << loader.getFileRange() << ", range " << r << ", vector returned size " << vectorReturned.size() << std::endl;

  ASSERT_LT(0 , vectorReturned.size());

  MPI_Barrier(MPI_COMM_WORLD);
}



// now register the test cases
REGISTER_TYPED_TEST_CASE_P(FASTALoaderTest, MPIOpenWithRange, MPIOpenWithRange2);


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
