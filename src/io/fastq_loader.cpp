/**
 * fastq_loader.cpp
 *
 *  Created on: Feb 18, 2014
 *      Author: tpan
 */

#include <sys/mman.h>
#include <sstream>
#include <cerrno>

#include "io/fastq_loader.hpp"


namespace bliss
{
  namespace io
  {

    fastq_loader::~fastq_loader()
    {
      if (data != nullptr)
        unmap(data, range);
    }

    fastq_loader::fastq_loader(std::string const & filename,
                               file_loader::range_type const & _range,
                               size_t const &total)
     : file_loader(filename)
    {
      // get adjustment first.
      range = adjust_to_record_boundaries(_range, total);

      // open the file
      data = map(range);
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
     *
     *  input_range:  offset is the global offset in file.  length is the length of the raw array.
     *  read through raw for first 4 lines, and adjust the input range.
     */

    /**
     * search for first occurence of @ from an arbitrary starting point.
     */
    size_t search_for_record_boundary(char const* curr, file_loader::range_type & range) throw(io_exception) {
       // need to look at 2 or 3 chars.  read 4 lines because of the line 2-3 combo below needs offset to next line 1.
       char first[4];
       uint64_t offsets[4] = {0, 0, 0, 0};   // units:  offset from beginning of file.

       // scan through to get the first At or Plus
       bool newlineChar = false;
       int currLineId = -1;        // current line id in the buffer
       int i = range.start;
       if (i == 0)  // beginning of file, treat specially, since there is no preceeding \n for the @
       {
         // no preceding \n.  populate directly
         ++currLineId;
         first[currLineId] = *curr;
         offsets[currLineId] = i;

         ++i;
         ++curr;
       }

       while (i < range.end && currLineId < 4)
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
           offsets[currLineId] = i;
           newlineChar = false;  // toggle off
         }
         //    else  // other characters in the line - don't care.

         ++i;
         ++curr;
       }


       ////// determine the position within a read record based on the first char of the first 3 lines.
       //     and adjust the starting positions and lengths
       // always shift the offset to the right (don't want to try to read to the end to get an end offset.
       size_t new_start = range.end;

       if (first[0] == '@')
       {
         if (first[1] != '@')  // lines 1,2
         {
           new_start = offsets[0];
         }
         else  // lines 4,1
         {
           new_start = offsets[1];
         }
       }
       else if (first[0] == '+')
       {
         if (first[1] == '@') // ambiguous
         {
           if (first[2] != '@')  // lines 4, 1, 2
           {
             new_start = offsets[1];
           }
           else  // lines 3, 4, 1
           {
             new_start = offsets[2];
           }
         }
         else  // lines 3, 4 (+, ^@)
         {
           new_start= offsets[2];
         }
       }
       else if (first[1] == '+')  // lines 2, 3;
       {
         new_start = offsets[3];
       }
       else if (first[1] == '@')  // lines 4,1
       {
         new_start = offsets[1];
       }
       else
       {
         std::stringstream ss;
         ss <<  "WARNING in file processing: file segment " << range.start << " - " << range.end << " does not contain valid FASTQ markers.";
         throw io_exception(ss.str());
       }
       return new_start;
    }

    /**
     * adjust the range's start and end
     */
    file_loader::range_type fastq_loader::adjust_to_record_boundaries(file_loader::range_type const & input,
                                                                             size_t const & total) throw(io_exception) {

      // init output range to same as input range
      file_loader::range_type aligned = input.align_to_page(page_size);

      // map the curr block using the page aligned range.
      char* curr = map(aligned);
      // search for the new start
      size_t newStart = search_for_record_boundary(curr + (aligned.start - aligned.block_start), aligned);
      // clean up.
      unmap(curr, aligned);

      // get the new starting offset.

      size_t newEnd = input.end;
      if (input.end < total)  {
        // not the last block

        // shift the range over to the next block.  keep all else same.
        file_loader::range_type next_range(input);
        next_range.start = input.end;
        next_range.end = input.end + (input.end - input.start);
        if (next_range.end > total) next_range.end = total;

        // then align the range to page boundary (front only)
        aligned = next_range.align_to_page(page_size);

        // map the next block using the page aligned range.
        curr = map(aligned);

        // search for the new end
        newEnd = search_for_record_boundary(curr + (aligned.start - aligned.block_start), aligned);

        // clean up.
        unmap(curr, aligned);
      }

      // construct the new output range to use.
      file_loader::range_type output(input);
      output.start = newStart;
      output.end = newEnd;
      // and page align it.
      return output.align_to_page(page_size);

    }



  } /* namespace io */
} /* namespace bliss */
