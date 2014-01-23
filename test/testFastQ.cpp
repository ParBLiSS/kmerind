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

#include "mpi.h"

#include "utils/logging.h"
#include "config.hpp"

template<typename T1, typename T2>
T1 computeOffset(T1 total, T2 np, T2 pid)
{
  assert(total > 0);
  assert(np > 0);
  assert(pid >= 0 && pid < np);

  if (np == 1)
    return 0;

  T1 div = total / static_cast<T1>(np);
  T2 rem = total % static_cast<T1>(np);
  if (pid < rem)
    return pid * (div + 1);
  else
    return pid * div + rem;
}

template<typename T1, typename T2>
T1 computeLength(T1 total, T2 np, T2 pid)
{
  assert(total > 0);
  assert(np > 0);
  assert(pid >= 0 && pid < np);

  if (np == 1)
    return total;

  T1 div = total / static_cast<T1>(np);
  T2 rem = total % static_cast<T1>(np);
  if (pid < rem)
    return div + 1;
  else
    return div;
}

void readFileStream(std::string const& filename, uint64_t offset, uint64_t length, char* result)
{
  std::ifstream ifs(filename);
  ifs.seekg(offset, ifs.beg);

  ifs.read(result, length);
  if (ifs)
    std::cout << "all characters read successfully.\n";
  else
    std::cout << "error: only " << ifs.gcount() << " could be read\n";

  ifs.close();
}



void readFileC(std::string const& filename, uint64_t offset, uint64_t length, char* result)
{
  FILE *fp = fopen(filename.c_str(),"r");
  fseek(fp, offset, SEEK_SET);
  fread(result, sizeof(char), length, fp);
  fclose(fp);
  std::cout << "all characters read successfully.\n";

}

void readFileMMap(std::string const& filename, uint64_t offset, uint64_t length, char* result)
{
  int fp = open(filename.c_str(), O_RDONLY);

  // change offset to align by page size.
  long page_size = sysconf(_SC_PAGE_SIZE);
  uint64_t p_offset = (offset / page_size) * page_size;
  uint64_t fp_offset = offset - p_offset;
  uint64_t p_length = length + fp_offset;

  // allow read only. any modifications are not written back to file.  read ahead.  do not swap
  char* temp = (char*)mmap(nullptr, p_length, PROT_READ, MAP_PRIVATE | MAP_POPULATE | MAP_LOCKED, fp, p_offset );

  memcpy(result, temp + fp_offset, length );

  munmap(temp, p_length);

  close(fp);
  std::cout << "all characters read successfully.\n";

}


int main(int argc, char* argv[]) {
  LOG_INIT();


  std::string filename("/home/tpan/src/pbil/test/test.fastq");

  int rank = 0, nprocs = 1;
#ifdef USE_MPI
  std::cout << "USE_MPI is set" << std::endl;

  // initialize MPI
  MPI_Init(&argc, &argv);

  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  // first thread gets the file size.
  uint64_t file_size;
  if (rank == 0)
  {
    struct stat filestat;
    stat(filename.c_str(), &filestat);
    file_size = static_cast<uint64_t>(filestat.st_size);

    // broadcast to all
#ifdef USE_MPI
    MPI_Bcast(&file_size, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
#endif
  }

  // then compute the offsets for current node
  uint64_t offset = computeOffset(file_size, nprocs, rank);
  uint64_t length = computeLength(file_size, nprocs, rank);

  // now open the file and begin reading
  char* buffer1 = new char[length];
  memset(buffer1, 0, length * sizeof(char));
  readFileStream(filename, offset, length, buffer1);

  char* buffer2 = new char[length];
  memset(buffer2, 0, length * sizeof(char));
  readFileC(filename, offset, length, buffer2);

  char* buffer3 = new char[length];
  memset(buffer3, 0, length * sizeof(char));
  readFileMMap(filename, offset, length, buffer3);


  delete [] buffer1;
  delete [] buffer2;
  delete [] buffer3;

#ifdef USE_MPI
  MPI_Finalize();

#endif

  return 0;


}
