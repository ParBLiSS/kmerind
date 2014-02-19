/**
 *
 * test reading the file
 *
 * 1,2,3,4 MPI procs, each reading a different part.
 *
 * step 1.  evaluate the speed of mmap, fopen, ifstream. - DONE. mmap is faster to get to an array representation.
 *    for streaming, there may not be much difference
 *
 * step 2.  verify reads are correct, including using overlaps - DONE
 * step 3.  parse fastq.
 *          issue here is if we have long sequence, with quality scores, some nodes may only have quality score and
 *          some may have only sequence.  for shorter reads, we also get that same situation at the boundaries.
 *
 *          best to do a 2 pass approach.  first pass mark the start of each read and corresponding quality score start
 *          then if # reads >> nprocs, split by reads else split by length
 *
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

#include "iterators/range.hpp"






void readFileStream(std::string const& filename, const uint64_t& offset, const uint64_t& length, char* result)
{
  std::ifstream ifs(filename);
  ifs.seekg(offset, ifs.beg);

  std::streamsize read = ifs.readsome(result, length);
  if (ifs)
    std::cout << "read " << read << " of "<< length << " at " << offset << std::endl;
  else
    std::cout << "error: only " << ifs.gcount() << " could be read\n";

  ifs.close();
}



void readFileC(std::string const& filename, const uint64_t& offset, const uint64_t& length, char* result)
{
  FILE *fp = fopen(filename.c_str(),"r");
  fseek(fp, offset, SEEK_SET);
  size_t read = fread_unlocked(result, sizeof(char), length, fp);
  fclose(fp);
  std::cout << "read " << read << " of "<< length << " at " << offset << std::endl;

}



/**
 * parse FASTQ actually
 * boundary cases:
 * 1. partition contains part of a sequence but more than 1 line.
 * 2. partition contains part of a line. e.g. only sequence or only quality score.
 *
 * if long sequence, have each node read a few segmented portions.
 *
 */



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

  if (rank == nprocs - 1)
    fprintf(stderr, "file size is %ld\n", file_size);

  // first generate an approximate partition.
  bliss::iterator::range r = bliss::iterator::range::block_partition(file_size, nprocs, rank);



  // now open the file and begin reading
  std::chrono::high_resolution_clock::time_point t1, t2;
  std::chrono::duration<double> time_span1;
  std::chrono::duration<double> time_span2;
  std::chrono::duration<double> time_span3;
  std::chrono::duration<double> time_span;


  size_t len = r.end - r.start;
  char* buffer1 = new char[len];
  memset(buffer1, 0, len * sizeof(char));
  t1 = std::chrono::high_resolution_clock::now();
  readFileStream(filename, r.start, len, buffer1);
  t2 = std::chrono::high_resolution_clock::now();
  time_span1 = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
  INFO("read from file stream " << "elapsed time: " << time_span1.count() << "s.");

  char* buffer2 = new char[len];
  memset(buffer2, 0, len * sizeof(char));
  t1 = std::chrono::high_resolution_clock::now();
  readFileC(filename, r.start, len, buffer2);
  t2 = std::chrono::high_resolution_clock::now();
  time_span2 = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
  INFO("read from C file " << "elapsed time: " << time_span2.count() << "s.");

//  char* buffer3 = new char[length];
//  memset(buffer3, 0, length * sizeof(char));
  t1 = std::chrono::high_resolution_clock::now();
  // now open the file
  bliss::io::fastq_loader loader(filename, r, file_size);
  t2 = std::chrono::high_resolution_clock::now();
  time_span3 = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
  INFO("read from MMap " << "elapsed time: " << time_span3.count() << "s.");

  // check to see if these are identical between the 3 methods - SAME.
  t1 = std::chrono::high_resolution_clock::now();
  int areSame = memcmp(buffer1, buffer2, len * sizeof(char));
  if (areSame != 0)
    fprintf(stderr, "buffer 1 and buffer 2 are different\n");
  t2 = std::chrono::high_resolution_clock::now();
  time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1) + time_span2;
  INFO("compared 1 2 " << "elapsed time: " << time_span.count() << "s.");

  uint64_t internaloffset = (r.start - rangeAligned.offset);
  t1 = std::chrono::high_resolution_clock::now();
  areSame = memcmp(buffer1, buffer3 + internaloffset, len * sizeof(char));
  if (areSame != 0)
    fprintf(stderr, "buffer 1 and buffer 3 are different\n");
  t2 = std::chrono::high_resolution_clock::now();
  time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1) + time_span3;
  INFO("compared 1 3 " << "elapsed time: " << time_span.count() << "s.");


  // write them out, concat on file system, compare to original - SAME.
  std::stringstream ss;
  std::ofstream ofs;
  ss << "result.1." << rank << ".txt";
  ofs.open(ss.str().c_str());
  ss.str(std::string());
  if (rank < (nprocs - 1))
    ofs.write(buffer1, len - overlap);
  else
    ofs.write(buffer1, len);
  ofs.close();

  ss << "result.2." << rank << ".txt";
  ofs.open(ss.str().c_str());
  ss.str(std::string());
  if (rank < (nprocs - 1))
    ofs.write(buffer2, len - overlap);
  else
    ofs.write(buffer2, len);
  ofs.close();

  ss << "result.3." << rank << ".txt";
  ofs.open(ss.str().c_str());
  ss.str(std::string());
  if (rank < (nprocs - 1))
    ofs.write(buffer3 + internaloffset, len - overlap);
  else
    ofs.write(buffer3 + internaloffset, len);
  ofs.close();

  delete [] buffer1;
  delete [] buffer2;

  RangeType<uint64_t> readAlignedRange;
  readAlignedRange = adjustRange(buffer3+internaloffset, range, file_size, rank, nprocs, readAlignedRange);
  printf("rank %d range read aligned: %ld %ld\n", rank, readAlignedRange.offset, readAlignedRange.length);
  closeFileMMap(fp, buffer3, rangeAligned.length);

  buffer1 = new char[readAlignedRange.length];
  memset(buffer1, 0, readAlignedRange.length * sizeof(char));
  t1 = std::chrono::high_resolution_clock::now();
  readFileStream(filename, readAlignedRange.offset, readAlignedRange.length, buffer1);
  t2 = std::chrono::high_resolution_clock::now();
  time_span1 = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
  INFO("read from file stream read aligned " << "elapsed time: " << time_span1.count() << "s.");

  rangeAligned = alignRange(readAlignedRange, rangeAligned, file_size, nprocs);
  printf("rank %d range read aligned block aligned: %ld %ld\n", rank, rangeAligned.offset, rangeAligned.length);

  t1 = std::chrono::high_resolution_clock::now();
  fp = readFileMMap(filename, rangeAligned.offset, rangeAligned.length, buffer3);
  t2 = std::chrono::high_resolution_clock::now();
  time_span3 = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
  INFO("read from MMap read aligned " << "elapsed time: " << time_span3.count() << "s.");

  internaloffset = (readAlignedRange.offset - rangeAligned.offset);
  t1 = std::chrono::high_resolution_clock::now();
  areSame = memcmp(buffer1, buffer3 + internaloffset, readAlignedRange.length * sizeof(char));
  if (areSame != 0)
    fprintf(stderr, "buffer 1 and buffer 3 are different\n");
  t2 = std::chrono::high_resolution_clock::now();
  time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1) + time_span3;
  INFO("compared 1 3 read aligned" << "elapsed time: " << time_span.count() << "s.");

  ss << "result.5." << rank << ".txt";
  ofs.open(ss.str().c_str());
  ss.str(std::string());
  if (rank < (nprocs - 1))
    ofs.write(buffer3 + internaloffset, readAlignedRange.length - overlap);
  else
    ofs.write(buffer3 + internaloffset, readAlignedRange.length);
  ofs.close();


  delete [] buffer1;
  //delete [] buffer3;
  closeFileMMap(fp, buffer3, readAlignedRange.length);


#ifdef USE_MPI
  MPI_Finalize();

#endif


  return 0;


}
