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

template<typename T1>
struct RangeType {
    T1 offset;
    T1 length;
    T1 overlap;
};

/**
 * compute the offset and length of the range for rank pid.
 * takes into account the block size (e.g. page size) to force alignment
 *
 * uses overlap.  kernel handles bringing whole page in for the last (not full) page.
 */
template<typename T1, typename T2>
RangeType<T1> computeRange(const T1& total, const T1& blocksize, const T1& overlap, const T2& np, const T2& pid)
{
  assert(total > 0);
  assert(blocksize > 0);
  assert(overlap >= 0);
  assert(np > 0);
  assert(pid >= 0 && pid < np);

  RangeType<T1> output;
  output.overlap = overlap;

  if (np == 1)
  {
    output.offset = 0;
    output.length = total;
    return output;
  }

  // spread the number of blocks first.
  T1 nblock = total / blocksize;

  T1 div = nblock / static_cast<T1>(np);
  T1 rem = nblock % static_cast<T1>(np);
  if (static_cast<T1>(pid) < rem)
  {
    output.offset = static_cast<T1>(pid) * (div + 1) * blocksize;
    output.length = (div + 1) * blocksize + overlap;
  }
  else
  {
    output.offset = (static_cast<T1>(pid) * div + rem) * blocksize;
    if (pid == (np - 1))
      output.length = total - ((static_cast<T1>(pid) * div + rem) * blocksize);
    else
      output.length = div * blocksize + overlap;
  }
  return output;
}


/**
 * adjust the range based on the first part of the partitions read from disk
 */
template<typename T1, typename T2>
RangeType<T1> adjustRange(RangeType<T1> const original, const T1& blocksize, RangeType<T1> const & left, RangeType<T1> const & right)
{
}

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
 * TODO: test on remotely mounted file system.
 */
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

}

/**
 *  for coordinated distributed quality score reads to map back to original.
 *    short sequence - partition and parse at read boundaries.  in FASTQ format
 *    long sequence - 2 mmapped regions?  should be FASTA or BED formats.   (1 sequence at 150B bases).
 *    older sequencers may have longer reads in FASTA.
 *
 * 1 read dataset contains 1 or more fastq files.
 *  number of files few.
 *  treat as virtual file?  with virtual file size offset to filename mapping?  assumption is that distributed file system
 *  can service a small number of files to disjoint sets of clients better than 1 file to all clients.  but this may introduce
 *  issues where different nodes may reach the same file at different times.
 *
 *  === assume files are relatively large.  simplify things by processing 1 file at a time.
 *
 *
 * 1 fastq/fasta file contains 1 or more sequences.
 *  6B sequences at 50 bases, or smaller
 *
 *  === assume files are relatively large.  simplify things by processing 1 file at a time.
 *
 * ids:
 *  file id.          (fid to filename map)
 *  read id in file   (read id to seq/qual offset map)
 *  position in read
 *
 * need to scan to compute read ids.
 *
 *
 * 1. determine read start positions
 * scan fastq from partial file to get the read start positions.
 * since parsing may start anywhere, can't rely on line 2 and 4's length equality to help.
 * for the same reason, can't rely on having \n@ or \n+ on lines 1 and 3.
 * DO rely on no newline char within sequence or quality.
 *
 * standard pattern is @x*
 *                     l*
 *                     +x*
 *                     q*
 *                     @x*
 *                     ...
 *
 *    where x is any char except \n.  does not occur at position 1.
 *          l is sequence alphabet so no + or @
 *          q is quality char including @, +, l, t (all others).
 *
 * for any 2 lines, possible sequences are
 *
 *  @x*     lines 1 and 2. no ambiguity
 *  l*
 *
 *  l*      lines 2 and 3, no ambiguity
 *  +x*
 *
 *  +x*     subcases:   +x* \n @q*  - ambiguous: 3-4 or 4-1.  resolve with next line (should be @x*) or with previous line (l*)
 *  q*                  +x* \n +q*  - lines 3 and 4, no ambiguity
 *                      +x* \n lq*  - lines 3 and 4, no ambiguity
 *                      +x* \n tq*  - lines 3 and 4, no ambiguity
 *
 *  q*      subcases:   @q* \n @x*  - lines 4 and 1, no ambiguity
 *  @x*                 +q* \n @x*  - ambiguous: 3-4 or 4-1.  resolve with next line (should be l*) or with previous line (+x*)
 *                      lq* \n @x*  - lines 4 and 1, no ambiguity
 *                      tq* \n @x*  - lines 4 and 1, no ambiguity
 *
 *  ambiguity is whether a pair of lines can be interpreted to be at more than 1 positions
 *
 *  use @ and + as anchors.
 *
 *  Boundary Case: \n@ or \n+ split by boundary.  this is resolvable by either overlap read range (not preferred)
 *  or by storing offsets after the 2nd or 4th \n character (at the end of range)
 *  Boundary case: a partition contains no line break (unknown)
 *  Boundary case: a partition contains only 1 line break (ambiguous)
 *  Boundary case: a partition contains only 2 line breaks and it's ambiguous.
 *
 *  boundary cases are unlikely, but can occur.
 *
 *  === for now, check and report instead of solving it.
 *
 *  one way to deal with these boundary cases is to communicate with neighbors.  (issue is if we need to go a few hops away).
 *  another way is to gather, then compute/scatter - 6B reads so 24B characters may be too much for the headnode directly
 *
 *  ===  for reads, expected lengths are maybe up to 1K in length, and files contain >> 1 seq
 *
 *
 *  hybrid.  if any node has ambiguity and has a small total read count, then use the gather method  (large sequences assumption), and have close to 2p reads.
 *      else (all have more than 2 line breaks) if any ambiguous, then enter into communication phase
 *  alternatively, do log(p) ordering - comm with neighbor, see if ambiguity resolved.  if not, comm with neighbor 2^i hops away.  repeat.
 *
 *
 *  issue:  can we fit the raw string into memory?  can we fit the result into memory?
 *    result has to be  able to fit.  N/P uint64_t.
 *    if result fits, the raw string has to fit.
 *
 *   implementation: assume file is read completely into memory.
 *
 *
 *  2 pass algorithm: since we need read ids.
 */
template<typename T1>
void scanFastQ(char const* raw, RangeType<T1> const & range,
               std::vector<uint64_t>& seqOffsets, std::vector<uint64_t>& qualOffsets,
               int rank, int nprocs) {

  char* curr = raw;

  // need to look at 2 or 3 chars.  read 4 lines because of the line 2-3 combo below needs offset to next line 1.
  char first[4];
  uint64_t offsets[4] = {0, 0, 0, 0};

  // scan through to get the first At or Plus
  bool newlineChar = false;
  int currLineId = -1;        // current line id in the buffer
  int i = 0;
  while (i < length && currLineId < 4)
  {
    // encountered a newline.  mark newline found, increment currLineId.
    if (*curr == '\n' && !newlineChar)
    {
      newlineChar = true;  // toggle on
    }
    else if (*curr != '\n' && newlineChar) // first char
    {
      ++currLineId;
      first[currLineId] = *curr;
      offsets[currLineId] = i + offset;
      newlineChar = false;  // toggle off
    }
//    else  // other characters in the line - don't care.

    ++i;
    ++curr;
  }

  uint64_t offsetsToPrevNode[2] = {0, 0};  // sequence, quality
  uint64_t startOffset = 0;   // output the current line id.
  int firstLineId = 0;

  ////// determine the position within a read record based on the first char of the first 3 lines.
  if (first[0] == '@')
  {
    if (first[1] != '@')  // lines 1,2
    {
      startOffset = offsets[0];
      firstLineId = 1;
    }
    else  // lines 4,1
    {
      offsetsToPrevNode[1] = offsets[0];
      startOffset = offsets[1];
      firstLineId = 4;
    }
  }
  else if (first[0] == '+')
  {
    if (first[1] == '@') // ambiguous
    {
      if (first[2] != '@')  // lines 4, 1, 2
      {
        offsetsToPrevNode[1] = offsets[0];
        startOffset = offsets[1];
        firstLineId = 4;
      }
      else  // lines 3, 4, 1
      {
        offsetsToPrevNode[1] = offsets[1];
        startOffset = offsets[2];
        firstLineId = 3;
      }
    }
    else  // lines 3, 4 (+, ^@)
    {
      offsetsToPrevNode[1] = offsets[1];
      startOffset = offsets[2];
      firstLineId = 3;
    }
  }
  else if (first[1] == '+')  // lines 2, 3;
  {
    offsetsToPrevNode[0] = offsets[0];
    offsetsToPrevNode[1] = offsets[2];
    startOffset = offsets[3];
    firstLineId = 2;
  }
  else if (first[1] == '@')  // lines 4,1
  {
    offsetsToPrevNode[1] = offsets[0];
    startOffset = offsets[1];
    firstLineId = 4;
  }

  ////// now move around data fragment  (do this to avoid rereading from file system)
  // if firstLineId = 2 or 3, move prev node's data to here.
  if (firstLineId == 2 || firstLineId == 3)
  {
    if (rank > 0)
    {
      // receive async (post this first)
    }

    if (rank < nprocs - 1)
    {
      // send async
    }

    // sync
    MPI_Barrier(MPI_WORLD_COMM);
  }
  else // if firstLineId = 4 or 1, move this node's data to prev node.
  {
    if (rank < nprocs - 1)
    {
      // receive async (post this first)
    }
    else
    {
      //send
    }

    // sync
    MPI_Barrier(MPI_WORLD_COMM);
  }


  seqOffsets.reserve(length);
  qualOffsets.reserve(length);



  while (i < length) {
    if (*curr == '\n')
    {
      ++i;
      ++curr;
      if (*curr == '@')
      {
        // tentatively line 1

      }
    }
    ++i;
    ++curr;
  }

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

  fprintf(stderr, "file size is %ld\n", file_size);


  // then compute the offsets for current node
  uint64_t page_size = sysconf(_SC_PAGE_SIZE);
  uint64_t overlap = 20;
  RangeType<uint64_t> range = computeRange(file_size, page_size, overlap, nprocs, rank);
  uint64_t offset = range.offset;
  uint64_t length = range.length;

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


  // check to see if these are identical between the 3 methods - SAME.
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


  // write them out, concat on file system, compare to original - SAME.
  std::stringstream ss;
  std::ofstream ofs;
  ss << "result.1." << rank << ".txt";
  ofs.open(ss.str().c_str());
  ss.str(std::string());
  if (rank < (nprocs - 1))
    ofs.write(buffer1, length - overlap);
  else
    ofs.write(buffer1, length);
  ofs.close();

  ss << "result.2." << rank << ".txt";
  ofs.open(ss.str().c_str());
  ss.str(std::string());
  if (rank < (nprocs - 1))
    ofs.write(buffer1, length - overlap);
  else
    ofs.write(buffer1, length);
  ofs.close();

  ss << "result.3." << rank << ".txt";
  ofs.open(ss.str().c_str());
  ss.str(std::string());
  if (rank < (nprocs - 1))
    ofs.write(buffer1, length - overlap);
  else
    ofs.write(buffer1, length);
  ofs.close();


  delete [] buffer1;
  delete [] buffer2;
  //delete [] buffer3;
  closeFileMMap(fp, buffer3, length);


#ifdef USE_MPI
  MPI_Finalize();

#endif

  return 0;


}
