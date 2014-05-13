/**
 * fastq_partition_helper.hpp
 *
 *  Created on: Feb 18, 2014
 *      Author: tpan
 */

#ifndef FASTQ_PARTITION_HELPER_HPP_
#define FASTQ_PARTITION_HELPER_HPP_

#include <sstream>

#include "config.hpp"

#include "common/base_types.hpp"
#include "io/file_loader.hpp"

namespace bliss
{
  namespace io
  {


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
     * for the same reason, can't rely on having \n@ or \n+ on lines 1 or 3.
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
     *  1 pass algorithm: we can use offset as ids.
     *
     * this class is now a functor.  The file_loader class takes this class as a parameter, for the purpose of repartitioning.
     *
     */
    template<typename Iterator, typename Range>
    struct FASTQPartitionHelper
    {
//
//        IteratorType begin(bool copying = false) {
//          return IteratorType(ParserType(data, range.start, copying), data, data + (range.end - range.start));
//        }
//        IteratorType end(bool copying = false) {
//          return IteratorType(ParserType(data, range.start, copying), data + (range.end - range.start));
//
//        }
        typedef typename Range::ValueType                                 SizeType;
        typedef typename std::iterator_traits<Iterator>::value_type       ValueType;


//        /**
//         * search for first occurence of @ from an arbitrary starting point.
//         *
//         * @param _data     start of iterator.
//         * @param target    this indicates the start and end position in the file.  NOTE that start position of 0 is treated differently.
//         *                  any other target is at the beginning of the mpi_proc range, do not count the first line.
//         * @return          position of the start of next read sequence (@).  if there is no complete sequence, return "end"
//         */
//        OffsetType operator()(const Iterator _data, const Range &target)
//                                                       throw (bliss::io::io_exception) {
//          Range t = target & parent; // intersection to bound target to between parent's ends.
//          if (target.start == target.end) return target.start;
//
//          // need to look at 2 or 3 chars.  read 4 lines because of the line 2-3 combo below needs offset to next line 1.
//          ValueType first[4];
//          OffsetType offsets[4] =
//          { target.end, target.end, target.end, target.end };   // units:  offset from beginning of file.
//
//          const Iterator data = _data;
//          // scan through to get the first At or Plus
//          bool newlineChar = false;
//          int currLineId = -1;        // current line id in the buffer
//          OffsetType i = target.start;
//
//          // TODO:
//          static_assert(false, "TODO:  review this logic to make sure it works for both the multithread and mpi parallel reading");
//
//          if (target.start == parent.start) // beginning of parent range, treat specially, since there is no preceeding \n for the @
//          {                                 // all others risk losing the first @ (will be caught by the previous)
//                                            // however, if parent.start is in the middle of a line, can't tell.
//                                            // approach is to place more weight on the other lines than the first line.
//            // no preceding \n.  populate directly
//            ++currLineId;
//            first[currLineId] = *data;
//            offsets[currLineId] = i;
//
//            ++i;
//            ++data;
//          }
//
//
//          ValueType c;
//          bool isEOL;
//          while (i < target.end && currLineId < 4)
//          {
//            c = *data;
//            isEOL = c == '\n';
//
//            // encountered a newline.  mark newline found, increment currLineId.
//            if (isEOL && !newlineChar)
//            {
//              newlineChar = true;  // toggle on
//            }
//            else if (newlineChar && !isEOL) // first char
//            {
//              ++currLineId;
//              first[currLineId] = c;
//              offsets[currLineId] = i;
//              newlineChar = false;  // toggle off
//            }
//            //    else  // other characters in the line - don't care.
//
//            ++i;
//            ++data;
//          }
//
//          ////// determine the position within a read record based on the first char of the first 3 lines.
//          //     and adjust the starting positions and lengths
//          // always shift the offset to the right (don't want to try to read to the end to get an end offset.
//          OffsetType new_pos = target.end;
//
//          if (first[0] == '@')
//          {
//            if (first[1] != '@')  // lines 1,2
//            {
//              new_pos = offsets[0];
//            }
//            else  // lines 4,1
//            {
//              new_pos = offsets[1];
//            }
//          }
//          else if (first[0] == '+')
//          {
//            if (first[1] == '@') // ambiguous
//            {
//              if (first[2] != '@')  // lines 4, 1, 2
//              {
//                new_pos = offsets[1];
//              }
//              else  // lines 3, 4, 1
//              {
//                new_pos = offsets[2];
//              }
//            }
//            else  // lines 3, 4 (+, ^@)
//            {
//              new_pos = offsets[2];
//            }
//          }
//          else if (first[1] == '+')  // lines 2, 3;
//          {
//            new_pos = offsets[3];
//          }
//          else if (first[1] == '@')  // lines 4,1
//          {
//            new_pos = offsets[1];
//          }
//          else
//          {
//            // is it an error not to find a fastq marker?
//            std::stringstream ss;
//            ss << "WARNING in file processing: file segment " << target.start
//               << " - " << target.end << " does not contain valid FASTQ markers.";
//            throw bliss::io::io_exception(ss.str());
//          }
//          return new_pos;
//
//        }

          /**
           * search for first occurence of @ from an arbitrary starting point.  Specifically, we are looking for the pairing of @ and + separated by 1 line
           *    e.g. .*\n[@][^\n]*\n[ATCGN]+\n[+].*\n.*
           *    note that this combination can occur at 2 places only.  at beginning of read/sequence, or if read name contains @ and the partition happens to start there.
           *
           *    we can also look for  + followed by @ in 2 lines.
           *    e.g. [ATCGN]*\n[+][^\n]*\n[^\n]+\n[@].*\n.*
           *
           *    using just 1 of these patterns would require us to search through up to 7 lines.
           *      pattern 1 works better if partition starts on lines 3 or 4, and pattern 2 works better if partition starts on 1 or 2.
           *      either case need to look at 5 lines.
           *
           * standard pattern is @x*
           *                     l*
           *                     +x*
           *                     q*
           *                     @x*
           *                     ....
           *                  x could be @, +, or other char except "\n".
           *                  l is alphabet, so no @, +, or "\n"
           *                  q could be @, +, or other char within some ascii range., no "\n"
           *
           *                  block start on line 1:
           *                    @ at start, or @ at some unfortunate position in the middle of name.  + is unambiguous.
           *                        if first block (target.start=parent.start) we would have @ at start.
           *                        if partition is right at beginning of line 1, we can afford to skip to next line 1.
           *                            also takes care of partition in the middle of line 1.
           *                  block start on line 2:
           *                    can't have @ anywhere in this line.  not possible
           *                  block start on line 3:
           *                    @ at some unfortunate position in the middle of line 3.  line 5 = line 1 can't have + at start.  not possible
           *                  block start on line 4:
           *                    quality line unfortunately starts with @ or contains @.  line 6 == line 2 can't have + at start.  not possible.
           *
           *  Algorithm:  read 4 lines, not including the part before the first "\n". treat first partition as if it was preceded with a "\n"
           *      look through the saved results, and look for @..+ or +..@.
           *
           *
           * @param _data     start of iterator.
           * @param target    this indicates the start and end position in the file.  NOTE that start position of 0 is treated differently.
           *                  any other target is at the beginning of the mpi_proc range, do not count the first line.
           * @return          position of the start of next read sequence (@).  if there is no complete sequence, return "end"
           */
          const SizeType operator()(const Iterator &_data, const Range &parent, const Range &target) const
                                                         throw (bliss::io::io_exception) {

            Range t = target & parent; // intersection to bound target to between parent's ends.
            if (t.start == t.end) return t.start;

            SizeType i = t.start;
            const Iterator data = _data;
            bool wasEOL = false;

            ////// remove leading newlines
            while ((i < t.end) && (*data == '\n')) {
              wasEOL = true;
              ++data;
              ++i;
            }

            //// if no more, end.
            if (i == t.end)  // this part only contained \n
              return t.end;

            ///// if at beginning of parent partition, treat as if previous line was EOL
            if (t.start == parent.start) { // beginning of parent range, treat specially, since there is no preceeding \n for the @
                                                // all other partitions will lose the part before the first "\n@" (will be caught by the previous partition)
              wasEOL = true;
            }



            ///// now read 4 lines
            // already has the first line first char..
            ValueType first[4];
            SizeType offsets[4] =
            { t.end, t.end, t.end, t.end };   // units:  offset from beginning of file.

            int currLineId = 0;        // current line id in the buffer

            // read the rest of the lines.
            ValueType c;
            while ((i < t.end) && currLineId < 4)
            {
              c = *data;

              // encountered a newline.  mark newline found.
              if (c == '\n')
              {
                wasEOL = true;  // toggle on
              }
              else
              {
                if (wasEOL) // first char
                {
                  first[currLineId] = c;
                  offsets[currLineId] = i;
                  wasEOL = false;  // toggle off
                  ++currLineId;
                }
              }
              //    else  // other characters in the line - don't care.

              ++data;
              ++i;
            }

            //////////// determine the position of a read by looking for @...+ or +...@
            // at this point, first[0] is pointing to first char after first newline, or in the case of first block, the first char.
            // everything in "first" are right after newline.
            // ambiguity from quality line?  no.
            if (first[0] == '@' && first[2] == '+')  return offsets[0];
            if (first[1] == '@' && first[3] == '+')  return offsets[1];
            if (first[2] == '@' && first[0] == '+')  return offsets[2];
            if (first[3] == '@' && first[1] == '+')  return offsets[3];


            // is it an error not to find a fastq marker?
            std::stringstream ss;
            ss << "ERROR in file processing: file segment \n" << "\t\t"
               << t << "\n\t\t(original " << target << ")\n"
               << "\t\tdoes not contain valid FASTQ markers.";
            throw bliss::io::io_exception(ss.str());

          }

    };

  } /* namespace io */
} /* namespace bliss */
#endif /* FASTQ_PARTITION_HELPER_HPP_ */
