/**
 * @file		quicktest.cpp
 * @ingroup
 * @author	tpan
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


#include "index/kmer_index_element.hpp"
#include "io/fastq_iterator.hpp"


using namespace std;
using namespace bliss::index;
using namespace bliss::io;

typedef KmerIndexElement<KmerSize<21>, uint64_t > KmerIndexType1;
typedef KmerIndexElementWithId<KmerSize<21>, uint64_t, FASTQSequenceId > KmerIndexType2;
typedef KmerIndexElementWithIdAndQuality<KmerSize<21>, uint64_t, FASTQSequenceId, float > KmerIndexType3;
typedef KmerIndexElementWithIdAndQuality<KmerSize<21>, uint64_t, FASTQSequenceId, double > KmerIndexType4;


void checkMemUsed(long long &phyMemUsed, long long &swapUsed, bool print) {
  //from http://stackoverflow.com/questions/63166/how-to-determine-cpu-and-memory-consumption-from-inside-a-process
  struct sysinfo memInfo;

  sysinfo (&memInfo);
  phyMemUsed = (memInfo.totalram - memInfo.freeram) * memInfo.mem_unit;
  swapUsed = (memInfo.totalswap - memInfo.freeswap) * memInfo.mem_unit;
  if (print)
    printf("physical mem used %lld, swap used %lld\n", phyMemUsed, swapUsed);
}

void memUsedvsBaseline(const long long &phyMemUsed, const long long &swapUsed, long long &phyMemUsed2, long long &swapUsed2) {
  checkMemUsed(phyMemUsed2, swapUsed2, false);

  printf("physical mem used new %lld, swap used new %lld\n", phyMemUsed2 - phyMemUsed, swapUsed2 - swapUsed);
}


int main(int argc, char** argv) {
  srand(1);

  long long phyMemUsed, swapUsed;
  long long phyMemUsed2, swapUsed2;



  int size = 1000000;

  printf("\n");
  checkMemUsed(phyMemUsed, swapUsed, true);
  {
    printf("KmerIndexElement copy.  element size %lu\n", sizeof(KmerIndexType1));
    std::vector<KmerIndexType1 > test;
    printf("\tcreated vector\n");

    memUsedvsBaseline(phyMemUsed, swapUsed, phyMemUsed2, swapUsed2);

    for (int i = 0; i < size; ++i) {
      KmerIndexType1 kmer;
      kmer.kmer = rand() % std::numeric_limits<uint64_t>::max();
      test.push_back(kmer);
    }
    printf("\tpopulated vector %lu\n", test.size());
    memUsedvsBaseline(phyMemUsed, swapUsed, phyMemUsed2, swapUsed2);


    test.clear();
    printf("\tcleared vector %lu\n", test.size());
    memUsedvsBaseline(phyMemUsed, swapUsed, phyMemUsed2, swapUsed2);

  }
  printf("\tdeallocated vector\n");
  memUsedvsBaseline(phyMemUsed, swapUsed, phyMemUsed2, swapUsed2);
  sleep(2);

  printf("\n");
  checkMemUsed(phyMemUsed, swapUsed, true);
  {
    printf("KmerIndexElement with move.  element size %lu\n", sizeof(KmerIndexType1));
    std::vector<KmerIndexType1 > test;
    printf("\tcreated vector\n");

    memUsedvsBaseline(phyMemUsed, swapUsed, phyMemUsed2, swapUsed2);

    for (int i = 0; i < size; ++i) {
      KmerIndexType1 kmer;
      kmer.kmer = rand() % std::numeric_limits<uint64_t>::max();
      test.push_back(std::move(kmer));
    }
    printf("\tpopulated vector %lu\n", test.size());
    memUsedvsBaseline(phyMemUsed, swapUsed, phyMemUsed2, swapUsed2);


    test.clear();
    printf("\tcleared vector %lu\n", test.size());
    memUsedvsBaseline(phyMemUsed, swapUsed, phyMemUsed2, swapUsed2);

  }
  printf("\tdeallocated vector\n");
  memUsedvsBaseline(phyMemUsed, swapUsed, phyMemUsed2, swapUsed2);
  sleep(2);


  printf("\n");
  checkMemUsed(phyMemUsed, swapUsed, true);
  {
    printf("Shared KmerIndexElement with move.  element size %lu\n", sizeof(KmerIndexType1));
    std::vector<KmerIndexType1 > test;
    printf("\tcreated vector\n");

    memUsedvsBaseline(phyMemUsed, swapUsed, phyMemUsed2, swapUsed2);

    KmerIndexType1 kmer;
    for (int i = 0; i < size; ++i) {
      kmer.kmer = rand() % std::numeric_limits<uint64_t>::max();
      test.push_back(std::move(kmer));
    }
    printf("\tpopulated vector %lu\n", test.size());
    memUsedvsBaseline(phyMemUsed, swapUsed, phyMemUsed2, swapUsed2);


    test.clear();
    printf("\tcleared vector %lu\n", test.size());
    memUsedvsBaseline(phyMemUsed, swapUsed, phyMemUsed2, swapUsed2);

  }
  printf("\tdeallocated vector\n");
  memUsedvsBaseline(phyMemUsed, swapUsed, phyMemUsed2, swapUsed2);
  sleep(2);





  printf("\n");
  checkMemUsed(phyMemUsed, swapUsed, true);
  {
    printf("KmerIndexElementWithId.  element size %lu\n", sizeof(KmerIndexType2));
    std::vector<KmerIndexType2 > test;
    printf("\tcreated vector\n");

    memUsedvsBaseline(phyMemUsed, swapUsed, phyMemUsed2, swapUsed2);

    for (int i = 0; i < size; ++i) {
      KmerIndexType2 kmer;
      kmer.kmer = rand() % std::numeric_limits<uint64_t>::max();
      kmer.id.composite = rand() % std::numeric_limits<uint64_t>::max();

      test.push_back(std::move(kmer));
    }
    printf("\tpopulated vector %lu\n", test.size());
    memUsedvsBaseline(phyMemUsed, swapUsed, phyMemUsed2, swapUsed2);


    test.clear();
    printf("\tcleared vector %lu\n", test.size());
    memUsedvsBaseline(phyMemUsed, swapUsed, phyMemUsed2, swapUsed2);

  }
  printf("\tdeallocated vector\n");
  memUsedvsBaseline(phyMemUsed, swapUsed, phyMemUsed2, swapUsed2);
  sleep(2);


  printf("\n");
  checkMemUsed(phyMemUsed, swapUsed, true);
  {
    printf("KmerIndexElementWithIdAndQuality.  element size %lu\n", sizeof(KmerIndexType3));
    std::vector<KmerIndexType3 > test;
    printf("\tcreated vector\n");

    memUsedvsBaseline(phyMemUsed, swapUsed, phyMemUsed2, swapUsed2);

    for (int i = 0; i < size; ++i) {
      KmerIndexType3 kmer;
      kmer.kmer = rand() % std::numeric_limits<uint64_t>::max();
      kmer.id.composite = rand() % std::numeric_limits<uint64_t>::max();
      kmer.qual = float(rand()) / float(RAND_MAX);

      test.push_back(std::move(kmer));
    }
    printf("\tpopulated vector %lu\n", test.size());
    memUsedvsBaseline(phyMemUsed, swapUsed, phyMemUsed2, swapUsed2);


    test.clear();
    printf("\tcleared vector %lu\n", test.size());
    memUsedvsBaseline(phyMemUsed, swapUsed, phyMemUsed2, swapUsed2);

  }
  printf("\tdeallocated vector\n");
  memUsedvsBaseline(phyMemUsed, swapUsed, phyMemUsed2, swapUsed2);
  sleep(2);


  printf("\n");
  checkMemUsed(phyMemUsed, swapUsed, true);
  {
    printf("KmerIndexElementWithIdAndQuality.  element size %lu\n", sizeof(KmerIndexType4));
    std::vector<KmerIndexType4 > test;
    printf("\tcreated vector\n");

    memUsedvsBaseline(phyMemUsed, swapUsed, phyMemUsed2, swapUsed2);

    for (int i = 0; i < size; ++i) {
      KmerIndexType4 kmer;
      kmer.kmer = rand() % std::numeric_limits<uint64_t>::max();
      kmer.id.composite = rand() % std::numeric_limits<uint64_t>::max();
      kmer.qual = float(rand()) / float(RAND_MAX);

      test.push_back(std::move(kmer));
    }
    printf("\tpopulated vector %lu\n", test.size());
    memUsedvsBaseline(phyMemUsed, swapUsed, phyMemUsed2, swapUsed2);


    test.clear();
    printf("\tcleared vector %lu\n", test.size());
    memUsedvsBaseline(phyMemUsed, swapUsed, phyMemUsed2, swapUsed2);

  }
  printf("\tdeallocated vector\n");
  memUsedvsBaseline(phyMemUsed, swapUsed, phyMemUsed2, swapUsed2);
  sleep(2);


  return 0;
}
