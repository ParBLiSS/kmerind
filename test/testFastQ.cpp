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
#include <cassert>
#include <iostream>
#include <fstream>
#include <chrono>
#include <unistd.h>

#include "mpi.h"

#include "utils/logging.h"
#include "config.hpp"

#include "iterators/range.hpp"
#include "io/file_loader.hpp"
#include "io/fastq_loader.hpp"






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

  /////////////// now try to open the file


  // first generate an approximate partition.
  bliss::io::file_loader::range_type r = bliss::io::file_loader::range_type::block_partition(file_size, nprocs, rank);
  std::cout << rank << " equipart: " << r << std::endl;
  std::cout << rank << " test block aligned: " << r.align_to_page(sysconf(_SC_PAGE_SIZE)) << std::endl;

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

  // check to see if these are identical between the 3 methods - SAME.
  t1 = std::chrono::high_resolution_clock::now();
  int areSame = memcmp(buffer1, buffer2, len * sizeof(char));
  if (areSame != 0)
    fprintf(stderr, "buffer 1 and buffer 2 are different\n");
  t2 = std::chrono::high_resolution_clock::now();
  time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1) + time_span2;
  INFO("compared 1 2 " << "elapsed time: " << time_span.count() << "s.");

  // write them out, concat on file system, compare to original - SAME.
  std::stringstream ss;
  std::ofstream ofs;
  ss << "result.1." << rank << ".txt";
  ofs.open(ss.str().c_str());
  ss.str(std::string());
  if (rank < (nprocs - 1))
    ofs.write(buffer1, len - r.overlap);
  else
    ofs.write(buffer1, len);
  ofs.close();

  ss << "result.2." << rank << ".txt";
  ofs.open(ss.str().c_str());
  ss.str(std::string());
  if (rank < (nprocs - 1))
    ofs.write(buffer2, len - r.overlap);
  else
    ofs.write(buffer2, len);
  ofs.close();

  delete [] buffer1;
  delete [] buffer2;

  {
    t1 = std::chrono::high_resolution_clock::now();
    // now open the file
    bliss::io::fastq_loader loader(filename, r, file_size);
    t2 = std::chrono::high_resolution_clock::now();
    time_span3 = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    INFO("MMap " << "elapsed time: " << time_span3.count() << "s.");

    std::cout << rank << " record adjusted " << loader.getRange() << std::endl;

    len = loader.getRange().end - loader.getRange().start;
    buffer1 = new char[len];
    memset(buffer1, 0, len * sizeof(char));
    t1 = std::chrono::high_resolution_clock::now();
    readFileStream(filename, loader.getRange().start, len, buffer1);
    t2 = std::chrono::high_resolution_clock::now();
    time_span1 = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    INFO("read from file stream " << "elapsed time: " << time_span1.count() << "s.");


    t1 = std::chrono::high_resolution_clock::now();
    areSame = memcmp(buffer1, loader.getData(), len * sizeof(char));
    if (areSame != 0)
      fprintf(stderr, "buffer 1 and buffer 3 are different\n");
    t2 = std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1) + time_span3;
    INFO("compared 1 3 " << "elapsed time: " << time_span.count() << "s.");

//// get eols
//    t1 = std::chrono::high_resolution_clock::now();
//    loader.test();
//    t2 = std::chrono::high_resolution_clock::now();
//    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1) + time_span3;
//    INFO("test " << "elapsed time: " << time_span.count() << "s.");
//
//    // get eols
//        t1 = std::chrono::high_resolution_clock::now();
//        loader.test2();
//        t2 = std::chrono::high_resolution_clock::now();
//        time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1) + time_span3;
//        INFO("test2 " << "elapsed time: " << time_span.count() << "s.");


//    // test parsing the sequences.
    t1 = std::chrono::high_resolution_clock::now();
    loader.get_sequence_positions();
    t2 = std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1) + time_span3;
    INFO("get reads " << "elapsed time: " << time_span.count() << "s.");


    //    // test parsing the sequences.
        t1 = std::chrono::high_resolution_clock::now();
        loader.assign_sequence_ids();
        t2 = std::chrono::high_resolution_clock::now();
        time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1) + time_span3;
        INFO("assign Ids " << "elapsed time: " << time_span.count() << "s.");



    ss << "result.3." << rank << ".txt";
    ofs.open(ss.str().c_str());
    ss.str(std::string());
    if (rank < (nprocs - 1))
      ofs.write(loader.getData(), len - r.overlap);
    else
      ofs.write(loader.getData(), len);
    ofs.close();

    delete [] buffer1;

  }



#ifdef USE_MPI
  MPI_Finalize();

#endif


  return 0;


}
