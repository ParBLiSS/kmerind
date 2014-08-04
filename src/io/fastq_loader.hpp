/**
 * fastq_partitioner.hpp
 *
 *  Created on: Feb 18, 2014
 *      Author: tpan
 */

#ifndef FASTQ_PARTITIONER_HPP_
#define FASTQ_PARTITIONER_HPP_

#include <sstream>
#include <iterator>  // for ostream_iterator
#include <algorithm> // for copy.

#include "config.hpp"

#include "common/base_types.hpp"
#include "io/file_loader.hpp"

namespace bliss
{
  namespace io
  {

    ///// subclass of FileLoader.  using static polymorphism via CRTP.


    /**
     *  for coordinated distributed quality score reads to map back to original.
     *    short sequence - partition and parse at read boundaries.  in FASTQ format
     *    long sequence - partition in block aligned form and read.  in FASTA.   (up to 150B bases).
     *    older sequencers may have longer reads in FASTA.
     *
     *
     *  === assume files are relatively large.  simplify things by processing 1 file at a time.
     *
     * ids:
     *  file id.          (fid to filename map)
     *  read id in file   (read id to seq/qual offset map)
     *  position in read
     *
     *  ===  for reads, expected lengths are maybe up to 1K in length, and files contain >> 1 seq
     *
     *  issue:  can we fit the raw string into memory?  can we fit the result into memory?
     *    result has to be  able to fit.  N/P uint64_t.
     *    if result fits, the raw string has to fit.
     *
     *  1 pass algorithm: we can use offset as ids.
     *
     */
    template<typename T, bool Preloading = false, bool Buffering = true,
        typename L1PartitionerT = bliss::partition::BlockPartitioner<bliss::partition::range<size_t> >,
        typename L2PartitionerT = bliss::partition::DemandDrivenPartitioner<bliss::partition::range<size_t> > >
    class FASTQLoader : public FileLoader<T, Preloading, Buffering, L1PartitionerT, L2PartitionerT,
                                          FASTQLoader<T, Preloading, Buffering, L1PartitionerT, L2PartitionerT> >
    {
      protected:
        typedef FileLoader<T, Preloading, Buffering, L1PartitionerT, L2PartitionerT,
            FASTQLoader<T, Preloading, Buffering, L1PartitionerT, L2PartitionerT> >    SuperType;

      public:
        /// exposing types from super
        typedef typename SuperType::InputIteratorType                   InputIteratorType;


        typedef typename SuperType::RangeType                            RangeType;

      protected:
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
         *  first block in the parent range is treated differently.
         *
         * @param _data       start of iterator.
         * @param partRange   the "full" range to which the target range belongs.  used to determine if the target is a "first" block and last block.  has to start and end with valid delimiters (@)
         * @param loadedRange the portion of the "full" range that's loaded in memory (via mmap or in DataBlock).  does NOT have to start and end with @
         * @param searchRange the range in which to search for a record.  the start and end position in the file.  does NOT have to start and end with @
         * @return            position of the start of next read sequence (@).  if there is no complete sequence, and is not at end of partition, throw exception.
         */
        template<typename Iterator>
        const size_t findStart(const Iterator &_data, const RangeType &partRange, const RangeType &loadedRange, const RangeType &searchRange) const
        throw (bliss::io::IOException) {

          typedef typename std::iterator_traits<Iterator>::value_type  ValueType;

          assert(partRange.contains(loadedRange));
          RangeType t = RangeType::intersect(searchRange, loadedRange); // intersection to bound target to between parent's ends.
          if (t.start == t.end) return t.start;

          size_t i = t.start;
          Iterator data(_data);
          bool wasEOL;

          ///// if at beginning of parent partition, treat as if previous line was EOL
          if (t.start == partRange.start) { // beginning of parent range, treat specially, since there is no preceeding \n for the @
            // all other partitions will lose the part before the first "\n@" (will be caught by the previous partition)
            // if this is called to partition, then partRange is the whole thing, with first part start with @ (implicit \n prior).
            // if this is called to chunk, then partRange is the loaded partition, already aligned to @, so previous is \n.
            wasEOL = true;
          } else {
            std::advance(data, t.start - loadedRange.start);
            wasEOL = false;
          }


          ////// remove leading newlines
          while ((i < t.end) && (*data == '\n')) {
            wasEOL = true;
            ++data;
            ++i;
          }

          //// if no more, end.
          if (i == t.end)  // this part only contained \n
            return t.end;





          ///// now read 4 lines
          // already has the first line first char..
          ValueType first[4] =
          { 0, 0, 0, 0 };
          size_t offsets[4] =
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

//          printf("chars: %c %c %c %c, %lu - %lu, lines: %d\n", first[0], first[1], first[2], first[3], target.start, target.end, currLineId);

          //////////// determine the position of a read by looking for @...+ or +...@
          // at this point, first[0] is pointing to first char after first newline, or in the case of first block, the first char.
          // everything in "first" are right after newline.
          // ambiguity from quality line?  no.
          if (first[0] == '@' && first[2] == '+')  return offsets[0];
          if (first[1] == '@' && first[3] == '+')  return offsets[1];
          if (first[2] == '@' && first[0] == '+')  return offsets[2];
          if (first[3] == '@' && first[1] == '+')  return offsets[3];


          /// finally, if nothing was found, because we are at the end of the partRange, then return the end.
          if (t.end == partRange.end) {
//            printf("\nchars: %c %c %c %c, %lu - %lu, lines: %d\n", first[0], first[1], first[2], first[3], t.start, t.end, currLineId);
            return t.end;
          }

          // is it an error not to find a fastq marker?
          std::stringstream ss;
          ss << "WARNING in file processing: file segment \n" << "\t\t"
              << t << "\n\t\t(original " << searchRange << ")\n"
              << "\t\tdoes not contain valid FASTQ markers.\n String:";
          std::ostream_iterator<typename std::iterator_traits<Iterator>::value_type> oit(ss);
          Iterator s(_data);
          std::advance(s, (t.start - loadedRange.start));
          Iterator e(_data);
          std::advance(e, (t.end - loadedRange.start));

          std::copy(s, e, oit);

//          printf("\nchars: %c %c %c %c, %lu - %lu, lines: %d\n", first[0], first[1], first[2], first[3], t.start, t.end, currLineId);

          throw bliss::io::IOException(ss.str());

        }



      public:
#if defined(USE_MPI)
        explicit FASTQLoader(const std::string &_filename, const MPI_Comm& _comm, const int _nThreads = 1, const size_t _chunkSize = 1) throw (IOException)
                             : SuperType(_filename, _comm, _nThreads, _chunkSize)
        {}
#endif

        explicit FASTQLoader(const std::string &_filename,
                            const int _nProcs = 1, const int _rank = 0,
                            const int _nThreads = 1, const size_t _chunkSize = 1) throw (IOException)
                            : SuperType(_filename, _nProcs, _rank, _nThreads, _chunkSize)
        {};

        virtual ~FASTQLoader() {};

        /// defining move constructor will disallow automatic copy constructor.
        /// move constructor and move assign operator
        FASTQLoader(FASTQLoader<T, Preloading, Buffering, L1PartitionerT, L2PartitionerT>&& other) : SuperType(std::forward(other)) {};

        FASTQLoader<T, Preloading, Buffering, L1PartitionerT, L2PartitionerT>& operator=(FASTQLoader<T, Preloading, Buffering, L1PartitionerT, L2PartitionerT>&& other) {
          if (this != &other) {
            this->SuperType::operator =(other);
          }
          return *this;
        }


        /**
         * search at the start and end of the block partition for some matching condition,
         * and set as new partition positions.
         */
        RangeType getNextL1BlockRangeImpl(const size_t pid) {

          // get basic range
          RangeType hint = this->L1Partitioner.getNext(pid);

          // modify it.
          hint.intersect(this->fileRange);
          size_t length = hint.size();

          RangeType next = hint + length;
          next.intersect(this->fileRange);

          // concatenate current and next ranges
          RangeType loadRange = RangeType::merge(hint, next);
          typename RangeType::ValueType block_start = RangeType::align_to_page(loadRange, this->pageSize);

          // map the content
          typename SuperType::InputIteratorType mappedData = this->map(loadRange);
          typename SuperType::InputIteratorType searchData = mappedData + (loadRange.start - block_start);

          // NOTE: assume that the range is large enough that we would not run into trouble with partition search not finding a starting position.

          // for output
          RangeType output(hint);

          // search for new start using finder
          try {
            output.start = findStart(searchData, this->fileRange, loadRange, hint);
            output.end = findStart(searchData, this->fileRange, loadRange, next);


          } catch (IOException& ex) {
            // did not find the end, so set e to next.end.
            // TODO: need to handle this scenario better - should keep search until end.
            WARNING(ex.what());

            printf("curr range: partition hint %lu-%lu, next %lu-%lu, file_range %lu-%lu\n",
                   hint.start, hint.end, next.start, next.end, this->fileRange.start, this->fileRange.end);
            printf("got an exception search for partition:  %s \n", ex.what());
            output.start = hint.end;
            output.end = hint.end;

            // either start or end are not found so return an empty range.
          }

          // clean up and unmap
          this->unmap(mappedData, loadRange);

          // readjust
          return output;
        }


        RangeType getNextL2BlockRangeImpl(const size_t tid) {
          assert(this->loaded);
          RangeType srcRange = this->L1Block.getRange();

          RangeType hint = this->L2Partitioner.getNext(tid);    // chunkPartitioner has this as atomic if needed
          size_t length = hint.size();

          // chunk size is already set to at least 2x record size - set during FileLoader::load()


          RangeType output(hint);

          if (length == 0) {
            //printf("WARNING: search range for next chunk is empty.  should not be here often\n");
          } else {
            /// search for start position.
            RangeType next = hint + length;
            next.intersect(srcRange);


            try {
              // search for start.
              output.start = findStart(this->L1Block.begin(), srcRange, srcRange, hint);
              output.end = findStart(this->L1Block.begin(), srcRange, srcRange, next);
            } catch (IOException& ex) {
              // did not find the end, so set e to next.end.
              // TODO: need to handle this scenario better - should keep search until end.
              WARNING(ex.what());

              printf("curr range: chunk %lu, hint %lu-%lu, next %lu-%lu, srcData range %lu-%lu, mmap_range %lu-%lu\n",
                     this->L2BlockSize, hint.start, hint.end, next.start, next.end, srcRange.start, srcRange.end, this->L1Block.getRange().start, this->L1Block.getRange().end);
              printf("got an exception search:  %s \n", ex.what());

              // either start or end are not found, so return nothing.
              output.start = hint.end;
              output.end = hint.end;
            }

          }

          //DEBUG("range " << range << " chunk " << cs << " chunkPos " << chunkPos);

          return output;
        }


        size_t getRecordSizeImpl(int iterations = 3) {

          //// TODO: if a different partitioner is used, seqSize may be incorrect.
          ////   seqSize is a property of the partitioner applied to the data.

          assert(this->loaded);

          size_t s, e;
          size_t ss = 0;
          RangeType parent(this->L1Block.getRange());
          RangeType r(this->L1Block.getRange());

          s = findStart(this->L1Block.begin(), this->fileRange, parent, r);

          for (int i = 0; i < iterations; ++i) {
            r.start = s + 1;   // advance by 1, in order to search for next entry.
            r.intersect(parent);
            e = findStart(this->L1Block.begin(), this->fileRange, parent, r);
            ss = std::max(ss, (e - s));
            s = e;
          }

          // DEBUG("Sequence size is " << ss);
          return ss;
        }


    };

  } /* namespace io */
} /* namespace bliss */
#endif /* FASTQ_PARTITIONER_HPP_ */
