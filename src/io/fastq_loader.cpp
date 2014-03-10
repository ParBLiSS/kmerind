/**
 * fastq_loader.cpp
 *
 *  Created on: Feb 18, 2014
 *      Author: tpan
 */

#include <sys/mman.h>
#include <sstream>
#include <cerrno>

#include "mpi.h"

#include "io/fastq_loader.hpp"
#include "iterators/filter_iterator.hpp"
#include "config.hpp"

namespace bliss
{
  namespace io
  {

    fastq_loader::~fastq_loader()
    {
    }  //  superclass destructor handles everything.

    fastq_loader::fastq_loader(const std::string &_filename,
                               const file_loader::range_type &_range,
                               const size_t &total, bool preload)
        : file_loader(_filename, preload)
    {
      filename = _filename;

      // adjust the range
      // internally opens the file twice with _range and shifted_range, search for record start/end
      // and construct new range and then open the file again.
      range = align_to_sequence(_filename, _range, total);

      // use the superclass constructor to map, with adjusted range.
      load();

    }

    fastq_loader::iterator fastq_loader::begin()
    {
      bliss::iterator::fastq_parser<char*> parser(data, range.start);

      return iterator(parser, data, data + (range.end - range.start));

    }

    fastq_loader::iterator fastq_loader::end()
    {
      bliss::iterator::fastq_parser<char*> parser(data, range.start);

      return iterator(parser, data + (range.end - range.start));

    }

    /**
     *  for coordinated distributed quality score reads to map back to original.
     *    short sequence - partition and parse at read boundaries.  in FASTQ format
     *    long sequence - partition in block aligned form and read.  in FASTA.   (up to 150B bases).
     *    older sequencers may have longer reads in FASTA.
     *
     * 1 read dataset contains 1 or more fastq files.
     *  number of files few.
     *
     *  === assume files are relatively large.  simplify things by processing 1 file at a time.
     *     *
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
     *  issue:  can we fit the raw string into memory?  can we fit the result into memory?
     *    result has to be  able to fit.  N/P uint64_t.
     *    if result fits, the raw string has to fit.
     *
     *   implementation: assume file is read completely into memory.
     *
     *
     *  2 pass algorithm: since we need read ids.
     *
     */

    /**
     * search for first occurence of @ from an arbitrary starting point.
     */
    size_t fastq_loader::find_sequence_start(
        char const* _data, const file_loader::range_type &range)
            throw (io_exception)
    {
      // need to look at 2 or 3 chars.  read 4 lines because of the line 2-3 combo below needs offset to next line 1.
      char first[4];
      size_t offsets[4] =
      { 0, 0, 0, 0 };   // units:  offset from beginning of file.

      // scan through to get the first At or Plus
      bool newlineChar = false;
      int currLineId = -1;        // current line id in the buffer
      size_t i = range.start;
      if (i == 0) // beginning of file, treat specially, since there is no preceeding \n for the @
      {
        // no preceding \n.  populate directly
        ++currLineId;
        first[currLineId] = *_data;
        offsets[currLineId] = i;

        ++i;
        ++_data;
      }

      while (i < range.end && currLineId < 4)
      {
        // encountered a newline.  mark newline found, increment currLineId.
        if (*_data == '\n' && !newlineChar)
        {
          newlineChar = true;  // toggle on
        }
        else if (*_data != '\n' && newlineChar) // first char
        {
          ++currLineId;
          first[currLineId] = *_data;
          offsets[currLineId] = i;
          newlineChar = false;  // toggle off
        }
        //    else  // other characters in the line - don't care.

        ++i;
        ++_data;
      }

      ////// determine the position within a read record based on the first char of the first 3 lines.
      //     and adjust the starting positions and lengths
      // always shift the offset to the right (don't want to try to read to the end to get an end offset.
      size_t new_pos = range.end;

      if (first[0] == '@')
      {
        if (first[1] != '@')  // lines 1,2
        {
          new_pos = offsets[0];
        }
        else  // lines 4,1
        {
          new_pos = offsets[1];
        }
      }
      else if (first[0] == '+')
      {
        if (first[1] == '@') // ambiguous
        {
          if (first[2] != '@')  // lines 4, 1, 2
          {
            new_pos = offsets[1];
          }
          else  // lines 3, 4, 1
          {
            new_pos = offsets[2];
          }
        }
        else  // lines 3, 4 (+, ^@)
        {
          new_pos = offsets[2];
        }
      }
      else if (first[1] == '+')  // lines 2, 3;
      {
        new_pos = offsets[3];
      }
      else if (first[1] == '@')  // lines 4,1
      {
        new_pos = offsets[1];
      }
      else
      {
        std::stringstream ss;
        ss << "WARNING in file processing: file segment " << range.start
           << " - " << range.end << " does not contain valid FASTQ markers.";
        throw io_exception(ss.str());
      }
      return new_pos;
    }

    /**
     * adjust the range's start and end
     */
    file_loader::range_type fastq_loader::align_to_sequence(
        const std::string &_filename, const file_loader::range_type &input,
        const size_t &total) throw (io_exception)
    {

      // get the new starting offset.
      size_t newStart = input.start;
      if (input.start > 0)
      {  // not first block
         // open the file for search of first record.
        file_loader f(_filename, input);

        // search for the new start
        newStart = find_sequence_start(f.getData(), f.getRange());

        // clean up when exiting the block, including close the file handle.
      }

      size_t newEnd = input.end;
      if (input.end < total)
      {   // not the last block

        // shift the range over to the next block.  keep all else same.
        file_loader::range_type next(input);
        next.start = input.end;
        next.end = input.end + (input.end - input.start);
        if (next.end > total)
          next.end = total;

        // then map the next range.
        file_loader f(_filename, next);

        // search for the new end
        newEnd = find_sequence_start(f.getData(), f.getRange());

        // clean up when if statement exits
      }

      // construct the new output range to use.
      file_loader::range_type output(newStart, newEnd, input.overlap,
                                     input.step);
      return output;

    }

    void fastq_loader::get_sequence_positions() throw (io_exception)
    {

      seqPositions.reserve((range.end - range.start) / 8); // at minimal, have @\nA\n+\n9\n for a record

      iterator fqiter = this->begin();
      iterator fqend = this->end();

      // walk through the iterator and
      for (; fqiter != fqend; ++fqiter)
      {
        seqPositions.push_back(*fqiter);
      }

      int max = seqPositions.size();
      uint32_t maxSeqID = 0;

      for (int i = max; i > 0; --i)
      {
        bliss::iterator::fastq_sequence<char*> result = seqPositions[max - i];

//        std::cout << "seq " << (max - i) << ": name=[" << std::string(result.name, result.name_end) << "], ";
//        std::cout << "id=[" << std::hex << result.id << std::endl;
        //        std::cout << "seq=[" << std::string(result.seq, result.seq_end) << "], ";
//        std::cout << "qual=[" << std::string(result.qual, result.qual_end) << "]" << std::endl;
        //printf("seq %d id = %lx, pos %d, file %d, seq msb %d, seq %d\n", (max-i), result.id, result.ids.pos, result.ids.file_id, result.ids.seq_msb, result.ids.seq_id);

        maxSeqID =
            (maxSeqID > result.ids.seq_id ? maxSeqID : result.ids.seq_id);
      }
      printf("maxId: %d\n", maxSeqID);

    }

    void fastq_loader::assign_sequence_ids() throw (io_exception)
    {
      // mpi based.

      int local_seq_count = seqPositions.size();
      seqIdStart = 0;
      int rank = 0, nprocs = 1;

#if defined(USE_MPI)

      MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);

#endif

//      printf("rank %d of %d has %d sequences\n", rank, nprocs, local_seq_count);

#if defined(USE_MPI)
      MPI_Exscan(&local_seq_count, &seqIdStart, 1, MPI_INT, MPI_SUM,
                 MPI_COMM_WORLD);

#endif

//      printf("rank %d of %d has prefix %ld\n", rank, nprocs, seqIdStart);

      for (int i = 0; i < local_seq_count; ++i)
      {
        seqPositions[i].id = i + seqIdStart;
      }

    }

    void fastq_loader::test()
    {
      char * curr = data;

      char** positions = new char*[(range.end - range.start) / 8];
      //std::vector<char*> positions;
//      positions.reserve((range.end - range.start) / 8);
      int i = 0;
      for (; curr < data + (range.end - range.start); ++curr, ++i)
      {
        if (*curr == '@')
          positions[i] = curr;
      }

      bool info = true;
      for (int j = 0; j < i; ++j)
      {
        info |= (positions[j] != 0);
      }

      delete[] positions;
    }

    void fastq_loader::test2()
    {
      char * curr = data;

      char** positions = (char**)malloc(
          (range.end - range.start) / 8 * sizeof(char*));
      //std::vector<char*> positions;
//      positions.reserve((range.end - range.start) / 8);
      int i = 0;
      for (; curr < data + (range.end - range.start); ++curr, ++i)
      {
        if (*curr == '@')
          positions[i] = (curr);
      }

      bool info = true;
      for (int j = 0; j < i; ++j)
      {
        info |= (positions[j] != 0);
      }

      free(positions);
    }

  } /* namespace io */
} /* namespace bliss */
