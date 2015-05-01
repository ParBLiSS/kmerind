/**
 * @file		quicktest.cpp
 * @ingroup
 * @author	Tony Pan <tpan7@gatech.edu>
 * @brief
 * @details
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */


#include <vector>
#include <limits>
#include <cstdlib>
#include <unistd.h>

#include "sys/types.h"
#include "sys/sysinfo.h"

#include "common/kmer.hpp"
#include "common/alphabets.hpp"
#include "common/alphabet_traits.hpp"

#include "retired/kmer_index_element.hpp"
#include "io/sequence_iterator.hpp"
#include "common/sequence.hpp"
#include "io/fastq_loader.hpp"


using namespace std;
using namespace bliss::index;
using namespace bliss::io;


typedef bliss::common::Kmer<21, bliss::common::DNA, uint64_t> KmerType;


typedef KmerIndexElement<KmerType > KmerIndexType1;
typedef KmerIndexElementWithId<KmerType, bliss::io::FASTQ::SequenceId > KmerIndexType2;
typedef KmerIndexElementWithIdAndQuality<KmerType, bliss::io::FASTQ::SequenceId, float > KmerIndexType3;
typedef KmerIndexElementWithIdAndQuality<KmerType, bliss::io::FASTQ::SequenceId, double > KmerIndexType4;


void checkMemUsed(long long &phyMemUsed, long long &swapUsed, bool print) {
  //from http://stackoverflow.com/questions/63166/how-to-determine-cpu-and-memory-consumption-from-inside-a-process
  struct sysinfo memInfo;

  sysinfo (&memInfo);
  phyMemUsed = (memInfo.totalram - memInfo.freeram) * memInfo.mem_unit;
  swapUsed = (memInfo.totalswap - memInfo.freeswap) * memInfo.mem_unit;
  if (print)
    INFOF("physical mem used %lld, swap used %lld\n", phyMemUsed, swapUsed);
}

void memUsedvsBaseline(const long long &phyMemUsed, const long long &swapUsed, long long &phyMemUsed2, long long &swapUsed2) {
  checkMemUsed(phyMemUsed2, swapUsed2, false);

  INFOF("physical mem used new %lld, swap used new %lld\n", phyMemUsed2 - phyMemUsed, swapUsed2 - swapUsed);
}


int main(int argc, char** argv) {
  srand(1);

  long long phyMemUsed, swapUsed;
  long long phyMemUsed2, swapUsed2;



  int size = 1000000;

  INFOF("\n");
  checkMemUsed(phyMemUsed, swapUsed, true);
  {
    INFOF("KmerIndexElement copy.  element size %lu\n", sizeof(KmerIndexType1));
    std::vector<KmerIndexType1 > test;
    INFOF("\tcreated vector\n");

    memUsedvsBaseline(phyMemUsed, swapUsed, phyMemUsed2, swapUsed2);

    for (int i = 0; i < size; ++i) {
      KmerIndexType1 kmer;
//      kmer.kmer = rand() % std::numeric_limits<uint64_t>::max();
      kmer.kmer.nextFromChar(rand() % (std::numeric_limits<uint8_t>::max() + 1));
      test.push_back(kmer);
    }
    INFOF("\tpopulated vector %lu\n", test.size());
    memUsedvsBaseline(phyMemUsed, swapUsed, phyMemUsed2, swapUsed2);


    test.clear();
    INFOF("\tcleared vector %lu\n", test.size());
    memUsedvsBaseline(phyMemUsed, swapUsed, phyMemUsed2, swapUsed2);

  }
  INFOF("\tdeallocated vector\n");
  memUsedvsBaseline(phyMemUsed, swapUsed, phyMemUsed2, swapUsed2);
  sleep(2);

  INFOF("\n");
  checkMemUsed(phyMemUsed, swapUsed, true);
  {
    INFOF("KmerIndexElement with move.  element size %lu\n", sizeof(KmerIndexType1));
    std::vector<KmerIndexType1 > test;
    INFOF("\tcreated vector\n");

    memUsedvsBaseline(phyMemUsed, swapUsed, phyMemUsed2, swapUsed2);

    for (int i = 0; i < size; ++i) {
      KmerIndexType1 kmer;
//      kmer.kmer = rand() % std::numeric_limits<uint64_t>::max();
      kmer.kmer.nextFromChar(rand() % (std::numeric_limits<uint8_t>::max() + 1));
      test.push_back(std::move(kmer));
    }
    INFOF("\tpopulated vector %lu\n", test.size());
    memUsedvsBaseline(phyMemUsed, swapUsed, phyMemUsed2, swapUsed2);


    test.clear();
    INFOF("\tcleared vector %lu\n", test.size());
    memUsedvsBaseline(phyMemUsed, swapUsed, phyMemUsed2, swapUsed2);

  }
  INFOF("\tdeallocated vector\n");
  memUsedvsBaseline(phyMemUsed, swapUsed, phyMemUsed2, swapUsed2);
  sleep(2);


  INFOF("\n");
  checkMemUsed(phyMemUsed, swapUsed, true);
  {
    INFOF("Shared KmerIndexElement with move.  element size %lu\n", sizeof(KmerIndexType1));
    std::vector<KmerIndexType1 > test;
    INFOF("\tcreated vector\n");

    memUsedvsBaseline(phyMemUsed, swapUsed, phyMemUsed2, swapUsed2);

    KmerIndexType1 kmer;
    for (int i = 0; i < size; ++i) {
//      kmer.kmer = rand() % std::numeric_limits<uint64_t>::max();
      kmer.kmer.nextFromChar(rand() % (std::numeric_limits<uint8_t>::max() + 1));
      test.push_back(std::move(kmer));
    }
    INFOF("\tpopulated vector %lu\n", test.size());
    memUsedvsBaseline(phyMemUsed, swapUsed, phyMemUsed2, swapUsed2);


    test.clear();
    INFOF("\tcleared vector %lu\n", test.size());
    memUsedvsBaseline(phyMemUsed, swapUsed, phyMemUsed2, swapUsed2);

  }
  INFOF("\tdeallocated vector\n");
  memUsedvsBaseline(phyMemUsed, swapUsed, phyMemUsed2, swapUsed2);
  sleep(2);





  INFOF("\n");
  checkMemUsed(phyMemUsed, swapUsed, true);
  {
    INFOF("KmerIndexElementWithId.  element size %lu\n", sizeof(KmerIndexType2));
    std::vector<KmerIndexType2 > test;
    INFOF("\tcreated vector\n");

    memUsedvsBaseline(phyMemUsed, swapUsed, phyMemUsed2, swapUsed2);

    for (int i = 0; i < size; ++i) {
      KmerIndexType2 kmer;
//      kmer.kmer = rand() % std::numeric_limits<uint64_t>::max();
      kmer.kmer.nextFromChar(rand() % (std::numeric_limits<uint8_t>::max() + 1));
      kmer.id.file_pos = rand() % std::numeric_limits<uint64_t>::max();

      test.push_back(std::move(kmer));
    }
    INFOF("\tpopulated vector %lu\n", test.size());
    memUsedvsBaseline(phyMemUsed, swapUsed, phyMemUsed2, swapUsed2);


    test.clear();
    INFOF("\tcleared vector %lu\n", test.size());
    memUsedvsBaseline(phyMemUsed, swapUsed, phyMemUsed2, swapUsed2);

  }
  INFOF("\tdeallocated vector\n");
  memUsedvsBaseline(phyMemUsed, swapUsed, phyMemUsed2, swapUsed2);
  sleep(2);


  INFOF("\n");
  checkMemUsed(phyMemUsed, swapUsed, true);
  {
    INFOF("KmerIndexElementWithIdAndQuality.  element size %lu\n", sizeof(KmerIndexType3));
    std::vector<KmerIndexType3 > test;
    INFOF("\tcreated vector\n");

    memUsedvsBaseline(phyMemUsed, swapUsed, phyMemUsed2, swapUsed2);

    for (int i = 0; i < size; ++i) {
      KmerIndexType3 kmer;
//      kmer.kmer = rand() % std::numeric_limits<uint64_t>::max();
      kmer.kmer.nextFromChar(rand() % (std::numeric_limits<uint8_t>::max() + 1));
      kmer.id.file_pos = rand() % std::numeric_limits<uint64_t>::max();
      kmer.qual = float(rand()) / float(RAND_MAX);

      test.push_back(std::move(kmer));
    }
    INFOF("\tpopulated vector %lu\n", test.size());
    memUsedvsBaseline(phyMemUsed, swapUsed, phyMemUsed2, swapUsed2);


    test.clear();
    INFOF("\tcleared vector %lu\n", test.size());
    memUsedvsBaseline(phyMemUsed, swapUsed, phyMemUsed2, swapUsed2);

  }
  INFOF("\tdeallocated vector\n");
  memUsedvsBaseline(phyMemUsed, swapUsed, phyMemUsed2, swapUsed2);
  sleep(2);


  INFOF("\n");
  checkMemUsed(phyMemUsed, swapUsed, true);
  {
    INFOF("KmerIndexElementWithIdAndQuality.  element size %lu\n", sizeof(KmerIndexType4));
    std::vector<KmerIndexType4 > test;
    INFOF("\tcreated vector\n");

    memUsedvsBaseline(phyMemUsed, swapUsed, phyMemUsed2, swapUsed2);

    for (int i = 0; i < size; ++i) {
      KmerIndexType4 kmer;
//      kmer.kmer = rand() % std::numeric_limits<uint64_t>::max();
      kmer.kmer.nextFromChar(rand() % (std::numeric_limits<uint8_t>::max() + 1));
      kmer.id.file_pos = rand() % std::numeric_limits<uint64_t>::max();
      kmer.qual = float(rand()) / float(RAND_MAX);

      test.push_back(std::move(kmer));
    }
    INFOF("\tpopulated vector %lu\n", test.size());
    memUsedvsBaseline(phyMemUsed, swapUsed, phyMemUsed2, swapUsed2);


    test.clear();
    INFOF("\tcleared vector %lu\n", test.size());
    memUsedvsBaseline(phyMemUsed, swapUsed, phyMemUsed2, swapUsed2);

  }
  INFOF("\tdeallocated vector\n");
  memUsedvsBaseline(phyMemUsed, swapUsed, phyMemUsed2, swapUsed2);
  sleep(2);


  return 0;
}
