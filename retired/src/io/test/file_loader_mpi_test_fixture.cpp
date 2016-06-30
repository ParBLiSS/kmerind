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
 * file_loader_test.cpp
 *   No separate MPI test.  FileLoader has no logical differences between MPI and non-MPI versions.  Only difference is MPI_comm is stored for convenience.
 *   note that this may change in the future if MPI_IO is used.
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

#include "io/file_loader.hpp"

#include "io/test/file_loader_test_fixtures.hpp"

using namespace bliss::io;


// indicate this is a typed test
TYPED_TEST_CASE_P(FileLoaderTest);


// normal test cases
TYPED_TEST_P(FileLoaderTest, open)
{

#ifdef USE_MPI
	int rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0) {
#endif

  typedef TypeParam FileLoaderType;
  typedef typename FileLoaderType::RangeType RangeType;
  using BaseValueType = typename FileLoaderTest<TypeParam>::ValueType;

  // get fileName

  FileLoaderType loader(this->fileName, 1, 0, 1, 64 * 1024, 64 * 1024);   // NOTE: larger block is MUCH faster.  2 * threads * 2048 bytes take a LONG time because of mmap each time...
  RangeType fr = loader.getFileRange();

  auto block = loader.getNextL1Block();
  RangeType r = block.getRange();

  // start at file starting point
  ASSERT_EQ(fr.start, r.start);

//  auto range_end = r.end;
  bool same = true;

  while (r.size() > 0) {
    BaseValueType* gold = new BaseValueType[r.size() + 1];
    FileLoaderTest<TypeParam>::readFilePOSIX(this->fileName, r.start, r.size(), gold);
    same &= std::equal(gold, gold + r.size(), block.begin());
    delete [] gold;

    block = loader.getNextL1Block();
    r = block.getRange();
//
//
//    ASSERT_EQ(r.start + FileLoaderType::get_overlap_size() + loader.get_L1_stride(), range_end);
//
//    range_end = r.end;
  }
  ASSERT_TRUE(same);

#ifdef USE_MPI
	}
#endif
}

#ifdef USE_MPI
TYPED_TEST_P(FileLoaderTest, open_mpi)
{
  typedef TypeParam                          FileLoaderType;
  typedef typename FileLoaderType::RangeType RangeType;
  using BaseValueType = typename FileLoaderTest<TypeParam>::ValueType;

  // get this->fileName

  FileLoaderType loader(this->fileName, MPI_COMM_WORLD, 1, 64 * 1024, 64 * 1024 );

  auto block = loader.getNextL1Block();
  RangeType r = block.getRange();
  bool same = true;

  while (r.size() > 0) {

    BaseValueType* gold = new BaseValueType[r.size() + 1];
    FileLoaderTest<TypeParam>::readFilePOSIX(this->fileName, r.start, r.size(), gold);
    same &= std::equal(gold, gold + r.size(), block.begin());
    delete [] gold;

    block = loader.getNextL1Block();
    r = block.getRange();
  }

  ASSERT_TRUE(same);

  MPI_Barrier(MPI_COMM_WORLD);

}
#endif


#ifdef USE_OPENMP
TYPED_TEST_P(FileLoaderTest, open_omp)
{
#ifdef USE_MPI
	int rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0) {
#endif


  typedef TypeParam                          FileLoaderType;
  typedef typename FileLoaderType::RangeType RangeType;
  using BaseValueType = typename FileLoaderTest<TypeParam>::ValueType;

  constexpr int nthreads = 4;


  FileLoaderType loader(this->fileName, 1, 0, nthreads, 64 * 1024, 1024 * 1024);

  auto l1block = loader.getNextL1Block();
  RangeType r1 = l1block.getRange();

  bool local_same[nthreads];
  bool same = true;

  while (r1.size() > 0) {

#pragma omp parallel num_threads(nthreads) shared(loader, local_same)
   {
    int tid = omp_get_thread_num();
//	   printf("thread %d start ", tid);

    auto l2block = loader.getNextL2Block(tid);
    RangeType r2 = l2block.getRange();

    local_same[tid] = true;

    while (r2.size() > 0) {
      BaseValueType* gold = new BaseValueType[r2.size() + 1];
      FileLoaderTest<TypeParam>::readFilePOSIX(this->fileName, r2.start, r2.size(), gold);
      local_same[tid] &= std::equal(gold, gold + r2.size(), l2block.begin());

      delete [] gold;

      l2block = loader.getNextL2Block(tid);
      r2 = l2block.getRange();
    }
//	   printf("thread %d end ", tid);

   } // end omp parallel

   l1block = loader.getNextL1Block();
    r1 = l1block.getRange();
  }

  for (int i = 0; i < nthreads; ++i) {
	  same &= local_same[i];
  }

  ASSERT_TRUE(same);

#ifdef USE_MPI
	}
#endif
}
#endif

#ifdef USE_MPI
#ifdef USE_OPENMP
TYPED_TEST_P(FileLoaderTest, open_mpi_omp)
{
  typedef TypeParam                          FileLoaderType;
  typedef typename FileLoaderType::RangeType RangeType;
  using BaseValueType = typename FileLoaderTest<TypeParam>::ValueType;

  constexpr int nthreads = 4;

  FileLoaderType loader(this->fileName, MPI_COMM_WORLD, nthreads, 64 * 1024, 1024 * 1024);

  auto l1block = loader.getNextL1Block();
  RangeType r1 = l1block.getRange();

  bool local_same[nthreads];
  bool same = true;

  while (r1.size() > 0) {

#pragma omp parallel num_threads(nthreads) shared(loader, local_same)
   {

    int tid = omp_get_thread_num();
    local_same[tid] = true;

    auto l2block = loader.getNextL2Block(tid);
    RangeType r2 = l2block.getRange();

    while (r2.size() > 0) {
      BaseValueType* gold = new BaseValueType[r2.size() + 1];
      FileLoaderTest<TypeParam>::readFilePOSIX(this->fileName, r2.start, r2.size(), gold);
      local_same[tid] &= std::equal(gold, gold + r2.size(), l2block.begin());
      delete [] gold;

      l2block = loader.getNextL2Block(tid);
      r2 = l2block.getRange();
    }

   } // end omp parallel

   for (int i = 0; i < nthreads; ++i) {
	   same &= local_same[i];
   }

   l1block = loader.getNextL1Block();
    r1 = l1block.getRange();
  }

  ASSERT_TRUE(same);
  MPI_Barrier(MPI_COMM_WORLD);
}
#endif
#endif


// now register the test cases
REGISTER_TYPED_TEST_CASE_P(FileLoaderTest,
#ifdef USE_OPENMP
		open_omp,
#endif
#ifdef USE_MPI
		open_mpi,
#ifdef USE_OPENMP
		open_mpi_omp,
#endif
#endif
		open);


// TODO: need to test mmap with multibyte elements.



int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

#if defined(USE_MPI)
  MPI_Init(&argc, &argv);
#endif

  result = RUN_ALL_TESTS();

#if defined(USE_MPI)
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
#endif

  return result;
}


//
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
//  typedef TestableFileLoader<TypeParam> FileLoaderType;
//
//  std::string fileName;
//
////  std::string err_regex = ".*file_loader.hpp.* Assertion .* failed.*";
////  ASSERT_EXIT(new FileLoaderType(1, 0, fileName), ::testing::KilledBySignal(SIGABRT), err_regex);
//  std::string err_regex = "ERROR: Filename Length is less than 1";
//  ASSERT_THROW(new FileLoaderType(1, 0, fileName), bliss::io::IOException);
//
//}
//
//TYPED_TEST_P(FileLoaderDeathTest, BadFilename)
//{
//  typedef TestableFileLoader<TypeParam> FileLoaderType;
//
//  std::string fileName;
//  fileName.assign("bad");
//
//  std::string err_regex = ".*file_loader.hpp.* Exception .* ERROR .* error 2: No such file or directory";
//  ASSERT_THROW(new FileLoaderType(1, 0, fileName), bliss::io::IOException);
//}
//
//TYPED_TEST_P(FileLoaderDeathTest, BadNProcs)
//{
//  typedef TestableFileLoader<TypeParam> FileLoaderType;
//
//  std::string fileName;
//  fileName.assign(PROJ_SRC_DIR);
//  fileName.append("/test/data/test.fastq");
//
////  std::string err_regex = ".*file_loader.hpp.* Exception .* ERROR .* error 2: No such file or directory";
//  std::string err_regex = ".*file_loader.hpp.* Assertion .* failed.*";
//  ASSERT_EXIT(new FileLoaderType(0 , 0, fileName), ::testing::KilledBySignal(SIGABRT), err_regex);
//}
//
//TYPED_TEST_P(FileLoaderDeathTest, BadRank)
//{
//  typedef TestableFileLoader<TypeParam> FileLoaderType;
//
//  std::string fileName;
//  fileName.assign(PROJ_SRC_DIR);
//  fileName.append("/test/data/test.fastq");
//
//  std::string err_regex = ".*file_loader.hpp.* Assertion .* failed.*";
//  ASSERT_EXIT(new FileLoaderType(1, -1, fileName), ::testing::KilledBySignal(SIGABRT), err_regex);
//}
//
//TYPED_TEST_P(FileLoaderDeathTest, BadNThreads)
//{
//  typedef TestableFileLoader<TypeParam> FileLoaderType;
//
//  std::string fileName;
//  fileName.assign(PROJ_SRC_DIR);
//  fileName.append("/test/data/test.fastq");
//
//  std::string err_regex = ".*file_loader.hpp.* Assertion .* failed.*";
//  ASSERT_EXIT(new FileLoaderType(1, 0, fileName, 0), ::testing::KilledBySignal(SIGABRT), err_regex);
//}
//
//TYPED_TEST_P(FileLoaderDeathTest, BadChunkSize)
//{
//  typedef TestableFileLoader<TypeParam> FileLoaderType;
//
//  std::string fileName;
//  fileName.assign(PROJ_SRC_DIR);
//  fileName.append("/test/data/test.fastq");
//
//  std::string err_regex = ".*file_loader.hpp.* Assertion .* failed.*";
//  ASSERT_EXIT(new FileLoaderType(1, 0, fileName, 0, 0), ::testing::KilledBySignal(SIGABRT), err_regex);
//}
//
//
//
//// now register the test cases
//REGISTER_TYPED_TEST_CASE_P(FileLoaderDeathTest, NoFilename, BadFilename, BadNProcs, BadRank, BadNThreads, BadChunkSize);



//INSTANTIATE_TYPED_TEST_CASE_P(Bliss, FileLoaderDeathTest, FileLoaderTestTypes);


