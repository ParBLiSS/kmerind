/**
 *
 * test reading the file
 *
 * 1,2,3,4 MPI procs, each reading a different part.
 *
 * step 1.  evaluate the speed of mmap, fopen, ifstream.
 */

#include <string>     // for std::string
#include <sys/stat.h>  // for stat
#include <sys/mman.h>  // for mmap
#include <cassert>
#include <cstdio>   // for fopen
#include <fcntl.h>  // for open
#include <iostream>
#include <fstream>
#include <chrono>
#include <unistd.h>


#include "mpi.h"

#include "utils/logging.h"
#include "config.hpp"

template<typename T1, typename T2>
T1 computeOffset(const T1& total, const T1& blocksize, const T2& np, const T2& pid)
{
  assert(total > 0);
  assert(np > 0);
  assert(pid >= 0 && pid < np);

  if (np == 1)
    return 0;

  // spread the number of blocks first.
  T1 nblock = total / blocksize;

  T1 div = nblock / static_cast<T1>(np);
  T1 rem = nblock % static_cast<T1>(np);
  if (static_cast<T1>(pid) < rem)
    return static_cast<T1>(pid) * (div + 1) * blocksize;
  else
    return (static_cast<T1>(pid) * div + rem) * blocksize;
}

template<typename T1, typename T2>
T1 computeLength(const T1& total, const T1& blocksize, const T2& np, const T2& pid)
{
  assert(total > 0);
  assert(np > 0);
  assert(pid >= 0 && pid < np);

  if (np == 1)
    return total;

  T1 nblock = total / blocksize;

  T1 div = nblock / static_cast<T1>(np);
  T1 rem = nblock % static_cast<T1>(np);
  if (static_cast<T1>(pid) < rem)
    return (div + 1) * blocksize;
  else if (pid == (np - 1))
    return total - ((static_cast<T1>(pid) * div + rem) * blocksize);
  else
    return div * blocksize;
}

void readFileStream(std::string const& filename, const uint64_t& offset, const uint64_t& length, char* result)
{
  std::ifstream ifs(filename);
  ifs.seekg(offset, ifs.beg);

  ifs.read(result, length);
  if (ifs)
    std::cout << "read " << length << " at " << offset << std::endl;
  else
    std::cout << "error: only " << ifs.gcount() << " could be read\n";

  ifs.close();
}



void readFileC(std::string const& filename, const uint64_t& offset, const uint64_t& length, char* result)
{
  FILE *fp = fopen(filename.c_str(),"r");
  fseek(fp, offset, SEEK_SET);
  fread(result, sizeof(char), length, fp);
  fclose(fp);
  std::cout << "read " << length << " at " << offset << std::endl;

}

int readFileMMap(std::string const& filename, const uint64_t& offset, const uint64_t& length, char*& result)
{
  int fp = open(filename.c_str(), O_RDONLY);

  // change offset to align by page size.

  // allow read only. any modifications are not written back to file.  read ahead.  do not swap
//  char* temp = (char*)mmap(nullptr, p_length, PROT_READ, MAP_PRIVATE | MAP_POPULATE | MAP_LOCKED, fp, p_offset );
  result = (char*)mmap(nullptr, length, PROT_READ, MAP_PRIVATE, fp, offset );

  if (result == MAP_FAILED)
  {
    int myerr = errno;
    fprintf(stderr, "ERROR in mmap: %d: %s\n", myerr, strerror(myerr));
    close(fp);
    return -1;
  }

  return fp;
}

void closeFileMMap(int fp, char* temp, uint64_t length)
{

  munmap(temp, length);

  close(fp);
  std::cout << "all characters read successfully.\n";

}


int main(int argc, char* argv[]) {
  LOG_INIT();


  std::string filename("/home/tpan/src/bliss/test/test.fastq");

  int rank = 0, nprocs = 1;
#ifdef USE_MPI
  // initialize MPI
  MPI_Init(&argc, &argv);

  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == 0)
    std::cout << "USE_MPI is set" << std::endl;
#endif

  // first thread gets the file size.
  uint64_t file_size = 0;
  if (rank == 0)
  {
    struct stat filestat;
    stat(filename.c_str(), &filestat);
    file_size = static_cast<uint64_t>(filestat.st_size);
  }

#ifdef USE_MPI
  // broadcast to all
  MPI_Bcast(&file_size, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
#endif

  fprintf(stderr, "file size is %ld\n", file_size);


  // then compute the offsets for current node
  uint64_t page_size = sysconf(_SC_PAGE_SIZE);
  uint64_t offset = computeOffset(file_size, page_size, nprocs, rank);
  uint64_t length = computeLength(file_size, page_size, nprocs, rank);

  // now open the file and begin reading
  std::chrono::high_resolution_clock::time_point t1, t2;
  std::chrono::duration<double> time_span1;
  std::chrono::duration<double> time_span2;
  std::chrono::duration<double> time_span3;
  std::chrono::duration<double> time_span;



  char* buffer1 = new char[length];
  memset(buffer1, 0, length * sizeof(char));
  t1 = std::chrono::high_resolution_clock::now();
  readFileStream(filename, offset, length, buffer1);
  t2 = std::chrono::high_resolution_clock::now();
  time_span1 = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
  INFO("read from file stream " << "elapsed time: " << time_span1.count() << "s.");

  char* buffer2 = new char[length];
  memset(buffer2, 0, length * sizeof(char));
  t1 = std::chrono::high_resolution_clock::now();
  readFileC(filename, offset, length, buffer2);
  t2 = std::chrono::high_resolution_clock::now();
  time_span2 = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
  INFO("read from C file " << "elapsed time: " << time_span2.count() << "s.");

//  char* buffer3 = new char[length];
//  memset(buffer3, 0, length * sizeof(char));
  char* buffer3;
  t1 = std::chrono::high_resolution_clock::now();
  int fp = readFileMMap(filename, offset, length, buffer3);
  t2 = std::chrono::high_resolution_clock::now();
  time_span3 = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
  INFO("read from MMap " << "elapsed time: " << time_span3.count() << "s.");


  // check to see if these are identical
  t1 = std::chrono::high_resolution_clock::now();
  int areSame = memcmp(buffer1, buffer2, length * sizeof(char));
  if (areSame != 0)
    fprintf(stderr, "buffer 1 and buffer 2 are different\n");
  t2 = std::chrono::high_resolution_clock::now();
  time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1) + time_span2;
  INFO("compared 1 2 " << "elapsed time: " << time_span.count() << "s.");

  t1 = std::chrono::high_resolution_clock::now();
  areSame = memcmp(buffer1, buffer3, length * sizeof(char));
  if (areSame != 0)
    fprintf(stderr, "buffer 1 and buffer 3 are different\n");
  t2 = std::chrono::high_resolution_clock::now();
  time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1) + time_span3;
  INFO("compared 1 3 " << "elapsed time: " << time_span.count() << "s.");



  delete [] buffer1;
  delete [] buffer2;
  //delete [] buffer3;
  closeFileMMap(fp, buffer3, length);


#ifdef USE_MPI
  MPI_Finalize();

#endif

  return 0;


}
