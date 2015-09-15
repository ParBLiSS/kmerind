/**
 * @file    fasta_loader.hpp
 * @ingroup bliss::io
 * @author  cjain
 * @brief   contains a FileLoader subclass to perform distributed and concurrent FASTA file access.
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */

#ifndef FASTA_PARTITIONER_HPP_
#define FASTA_PARTITIONER_HPP_

#include <sstream>
#include <iterator>  // for ostream_iterator
#include <algorithm> // for copy.

#include "bliss-config.hpp"

#include "common/sequence.hpp"
#include "common/base_types.hpp"
#include "io/file_loader.hpp"
#include "mxx/datatypes.hpp"

#include "io/mxx_support.hpp"
#include "iterators/filter_iterator.hpp"


// include MPI
#if defined(USE_MPI)
#include "mpi.h"
#endif

namespace bliss
{
  namespace io
  {
    /// class to indicate FASTA format
    struct FASTA {

    };





    /**
     * @class bliss::io::FASTQParser
     * @brief Functoid encapsulating increment functionality traversing data in an iterator record by record, by parsing and searching for record boundaries.
     * @details   The purpose of this class is to abstract the increment operations in
     *    an iterator that transforms a block of data into a collection of FASTQ records. Using this class allows us to reuse
     *    an sequence parsing iterator, but potentially for different sequence record formats.
     *
     *    This class assumes that the iterator at start is pointing to the beginning of a FASTQ record already.
     *    Increment is done by walking through the data and recording the positions of the start and end of each
     *    of the 4 lines in a FASTQ record.  Because we are always searching to the end of a record (by reading 4 lines always),
     *    successive calls to increment always are aligned correctly with record boundaries
     *
     * @note  NOT THREAD SAFE.
     *    Also, this parser needs to store global information for entire file, so is only compatible with block partitioner for L1.
     *    L1 partitioner's result needs to be stored.
     *
     * @tparam Iterator   The underlying iterator to be traversed to generate a Sequence
     */
    template <typename Iterator>
    class FASTAParser2 : public bliss::io::BaseFileParser<Iterator >
    {

        // let another BaseFileParser with some other iterator type be a friend so we can convert.
        template <typename Iterator2> friend class FASTAParser2;

      protected:
//        struct NotEOL {
//          bool operator()(typename std::iterator_traits<Iterator>::value_type const & x) {
//            return (x != bliss::io::BaseFileParser<Iterator>::eol && x != bliss::io::BaseFileParser<Iterator>::cr );
//          }
//        } filter;

        /// internal filter iterator to remove eol and cr characters.
//        using BaseCharIterator = typename ::bliss::iterator::filter_iterator<NotEOL, Iterator>;
        using BaseCharIterator = Iterator;

      public:
        /// sequence
        using SequenceType = typename ::bliss::common::Sequence<BaseCharIterator>;
        using SequenceIdType = typename SequenceType::IdType;


      protected:
        /**
         * @typedef RangeType
         * @brief   range object types
         */
        using RangeType = bliss::partition::range<size_t>;

        /**
         * meanings are record start, sequence start, sequence end, and record id.
         */
        using SequenceOffsetsType = ::std::tuple<typename RangeType::ValueType, typename RangeType::ValueType, typename RangeType::ValueType, size_t >;
        using SequenceVecType = std::vector<SequenceOffsetsType >;
        /**
         * @brief internal storage marking each sequence.
         */
        SequenceVecType sequences;


      public:

        /// default constructor.
        FASTAParser2() {};

        /// default destructor
        ~FASTAParser2() {};



        /// converting constructor.
        template <typename Iterator2>
        FASTAParser2(FASTAParser2<Iterator2> const & other) : sequences(other.sequences.begin(), other.sequences.end()) { }
        /// converting assignment operator that can transform the base iterator type.
        template <typename Iterator2>
        FASTAParser2<Iterator>& operator=(FASTAParser2<Iterator2> const & other) {
          sequences.assign(other.sequences.begin(), other.sequences.end());
          return *this;
        }



//        /**
//         * @brief  given a block, find the starting point that best aligns with the first sequence object (here, just just the actual start)
//         * @note   not overridden since here we just look for the beginning of the assigned block, same as the BaseFileParser logic.
//         */
//        const std::size_t findStart(const Iterator &_data, const RangeType &parentRange, const RangeType &inMemRange, const RangeType &searchRange) const
//        throw (bliss::io::IOException) {
//          return bliss::io::BaseFileParser<Iterator>::findStart(_data, parentRange, inMemRange, searchRange);
//        }

        // TODO: consolidate this.

#ifdef USE_MPI
        /**
         * @brief   Keep track of all '>' and the following EOL character to save the positions of the fasta record headers for each fasta record in the local block
         * @attention   This function should be called by all the MPI ranks in MPI communicator
         * @note        This function should be called by only 1 thread in a process.
         *              This function assumes that the searchRange are NOT OVERLAPPING between processes.  Else we could have messy duplicates.
         *
         * @tparam      Iterator      type of iterator for data to traverse.  raw pointer if no L2Buffering, or vector's iterator if buffering.
         * @param[in]   _data         start of iterator.
         * @param[in]   parentRange   the "full" range to which the inMemRange belongs.  used to determine if the target is a "first" block and last block.  has to start and end with valid delimiters (@)
         * @param[in]   inMemRange    the portion of the "full" range that's loaded in memory (e.g. in DataBlock).  does NOT have to start and end with @
         * @param[in]   searchRange   the range in which to search for a record.  the start and end position in the file.  does NOT have to start and end with @. should be inside inMemRange
         * @param[out]  sequences  vector of positions of the fasta sequence headers.
         *              Each pair in the vector represents the position of '>' and '\n' in the fasta record header
         */
        void init_for_iterator(const Iterator &_data, const RangeType &parentRange, const RangeType &inMemRange, const RangeType &searchRange, MPI_Comm comm)
          throw (bliss::io::IOException) {

          // now populate sequences

          if (!::std::is_same<typename std::iterator_traits<Iterator>::iterator_category, std::random_access_iterator_tag>::value)
            WARNING("WARNING: countSequenceStarts input iterator is not random access iterator.  lower.");

          using TT = typename std::iterator_traits<Iterator>::value_type;

          // make sure searchRange is inside parent Range.
          RangeType r = RangeType::intersect(searchRange, parentRange);


          int p = 1, myRank = 0;
          MPI_Comm_size(comm, &p);
          MPI_Comm_rank(comm, &myRank);

          // construct the search Range so they do not overlap.  This is necessary so that the sequences tuples are computed correctly.
          if (myRank == (p-1))
            mxx::left_shift(r.start, comm);
          else 
            r.end = mxx::left_shift(r.start, comm);

          // clear the starting position storage.
          sequences.clear();

          TT prev_char = '\n';  // default to EOL char.

          //============== OUTLINE
          // search for "\n>"
          // if internal, fine.
          // if at boundary,
          // if first char of proc i is ">", and rank is 0, can assume this is a starting point.  else, if prev char is "\n", then this is a start.


          DEBUGF("Rank %d search range [%lu, %lu)\n", myRank, r.start, r.end);


          // ===== get the previous character from the next smaller, non-empty rank.
          // do by excluding the empty procs
          MPI_Comm in_comm;
          mxx2::split_communicator_by_function(comm, [&r](){ return (r.size() == 0) ? MPI_UNDEFINED : 1; }, in_comm);

          // if no data, then in_comm is null, and does not participate further.
          if (in_comm == MPI_COMM_NULL) return;

          //============== GET LAST CHAR FROM PREVIOUS PROCESSOR, to check if first char is at start of line.
            // for the nonempty ones, shift right by 1.
           int in_p = 1, in_rank = 0;
           MPI_Comm_size(in_comm, &in_p);
           MPI_Comm_rank(in_comm, &in_rank);


            // get the current last character.
            TT last_char = 0;  // default for empty.  note that the real last char from the last non-empty proc is shifted to the right place.
            {
              Iterator last = _data;
              if (r.size() > 0) {
                std::advance(last, r.size() - 1);
                last_char = *last;
              }
            }
            // right shift using the non-empty procs.
            prev_char = mxx::right_shift(last_char, in_comm);

            if (in_rank == 0) prev_char = '\n';
            //=====DONE===== GET LAST CHAR FROM PREVIOUS PROCESSOR, to check if first char is at start of line.


            //============== GET THE POSITION OF START OF EACH LINE.
            // local storage
            using LL = std::pair< typename RangeType::ValueType, uint8_t>;
            std::vector< LL > line_starts;

            // initialize the iterators
            Iterator it = _data;
            Iterator end = it;
            std::advance(end, r.size());
            auto i = searchRange.start;

            // ==== find all the beginning of line positions - exclude consecutive eols.

            // now determine if the first character is a line start
            if (prev_char == '\n') {
              // previous char is eol, so add the position here if it's not another eol, and encode it for header vs not header
              line_starts.emplace_back(i, ((*it == ';') || (*it == '>')) ? 1 : 0);
            }

            // now search for occurrences of '\n', since we've already handled the case where the first char is '>'
            Iterator it2 = it;
            ++it;
            ++i;
            while (it != end) {
              if (*it2 == '\n'){
                // previous char is eol, so add the position here if it's not another eol, and encode it for header vs not
                line_starts.emplace_back(i, ((*it == ';') || (*it == '>')) ? 1 : 0);
              }
              ++it;
              ++it2;
              ++i;
            }
            if (in_rank == (in_p - 1)) {
              // add a eof entry.
              line_starts.emplace_back(parentRange.end, 1);
            }
            //=====DONE===== GET THE POSITION OF START OF EACH LINE.


            //============== KEEP JUST THE FIRST POSITION OF CONSECUTIVE HEADER OR SEQUENCE LINES  (== unique)
            auto new_end = line_starts.end();
            // filter out duplicates locally - i.e. consecutive lines with > or ; (true), and consecutive lines with standard alphabet (false).
            new_end = mxx2::unique_contiguous(line_starts, in_comm, [](LL const & x, LL const & y) {
              return x.second == y.second;
            });
            
            line_starts.erase(new_end, line_starts.end());

            //=====DONE===== KEEP JUST THE FIRST POSITION OF CONSECUTIVE HEADER OR SEQUENCE LINES

//            for (auto x : line_starts)
//              DEBUGF("R %d line starts %lu, new record? %d\n", myRank, x.first, x.second);


            //============== CONVERT FROM LINE START POSITION TO SEQUENCE OBJECTS
            int in_rank2 = in_rank;
            int in_p2 = in_p;

            // split comm again for non-empty line start vectors
            MPI_Comm in_comm2;
            mxx2::split_communicator_by_function(in_comm, [&line_starts](){ return line_starts.size() == 0 ? MPI_UNDEFINED : 1; }, in_comm2);

            // get first from next proc (
            if (in_comm2 != MPI_COMM_NULL) {
              MPI_Comm_rank(in_comm2, &in_rank2);
              MPI_Comm_size(in_comm2, &in_p2);

              // do a shift
              LL temp = line_starts.front();
              LL next = mxx::left_shift(temp, in_comm2);

              // insert into array.
              if (in_rank2 < (in_p2 - 1)) {
                line_starts.emplace_back(next);
              }

              MPI_Comm_free(&in_comm2);
            }

            // get the second from next proc (so we can construct elements with header start, seq start, and seq end.)
            // each proc should have at least 2, except for the last proc may have 1.
            // split comm again for non-empty line start vectors
            mxx2::split_communicator_by_function(in_comm, [&line_starts](){ return line_starts.size() <= 1 ? MPI_UNDEFINED : 1; }, in_comm2);

            // get first from next proc (

            if (in_comm2 != MPI_COMM_NULL) {
              MPI_Comm_rank(in_comm2, &in_rank2);
              MPI_Comm_size(in_comm2, &in_p2);

              // do a shift
              LL temp = line_starts[1];
              LL next = mxx::left_shift(temp, in_comm2);

              // insert into array.
              if (in_rank2 < (in_p2 - 1)) {
                line_starts.emplace_back(next);
              }

              MPI_Comm_free(&in_comm2);
            }

            // now only work with procs with at least 1 complete entry (header start, seq start, seq end
            mxx2::split_communicator_by_function(in_comm, [&line_starts](){ return line_starts.size() < ((line_starts.size() > 0 && line_starts.front().second == 0) ? 4 : 3) ? MPI_UNDEFINED : 1; }, in_comm2);


            if (in_comm2 != MPI_COMM_NULL) {


              // insert first entry into sequences
              typename RangeType::ValueType pos, pos2, pos3;
              int k = 0;

              if (line_starts[k].second == 0 ) // not starting with a header.  move to next
                ++k;

              pos = line_starts[k].first;  // first start.
              ++k;

              // insert all other entries into sequences.   should be globally unique at this point.
              for ( ; (k + 1) < line_starts.size() ; k += 2) {

                pos2 = line_starts[k].first;
                pos3 = line_starts[k + 1].first;

                // insert if there is a next line or not.
                sequences.emplace_back(pos, pos2, pos3, k/2);

                pos = pos3;
              }

              // get prefix sum of entries.
              size_t entries = sequences.size();
              entries = mxx2::scan::exscan(entries, std::plus<size_t>(), in_comm2);
              int in_rank2 = 0;
              MPI_Comm_rank(in_comm2, &in_rank2);
              MPI_Comm_free(&in_comm2);
              if (in_rank2 > 0) {
                  for (size_t ii = 0; ii < sequences.size(); ++ii) {
                    std::get<3>(sequences[ii]) += entries;
                  }
              }



            }
            //=====DONE===== CONVERT FROM LINE START POSITION TO SEQUENCE OBJECTS


//            for (auto x : sequences)
//              DEBUGF("R %d before spreading [%lu, %lu, %lu), id %lu\n", myRank, std::get<0>(x), std::get<1>(x), std::get<2>(x), std::get<3>(x));


            using SL = SequenceOffsetsType;
            //============== NOW SPREAD THE LAST ENTRY FORWARD (right).  empty range (r) does not participate

            //== first define segments.  empty is 0, and non-empty is 1 (start of segment)
            uint8_t sstart = sequences.size() == 0 ? 0 : 1;
            // convert to unique segments
            size_t segf = mxx2::segment::to_unique_segment_id_from_start<size_t>(sstart, 0, in_comm);

            // get the last entry.
            SL headerPos = SL();
            if (sstart) headerPos = sequences.back();

            //== do seg scan (to fill in empty procs).

            // now do inclusive scan to spread the values forward
            headerPos = mxx2::segmented_scan::scan(headerPos, segf, [](SL & x, SL & y){ return (std::get<0>(x) > std::get<0>(y) ? x : y); }, in_comm);

            // then shift by 1
            headerPos = mxx::right_shift(headerPos, in_comm);

            // seg scan to fill empty then shift is NOT the same as seg exscan


            //== and merge the results with next.  If there is no sequences content, then need to insert.
            if (in_rank > 0) {
              // no special case for overlap.  downstream processing should be aware of the overlap.
              if ((sequences.size() == 0) || (std::get<0>(sequences.front()) > r.start)) {
                sequences.insert(sequences.begin(), headerPos);
              }

            }
            //=====DONE===== NOW SPREAD THE LAST ENTRY FORWARD.  empty range (r) does not participate

            if (in_comm != MPI_COMM_NULL) MPI_Comm_free(&in_comm);


//            for (auto x : sequences)
//              DEBUGF("R %d header range [%lu, %lu, %lu), id %lu\n", myRank, std::get<0>(x), std::get<1>(x), std::get<2>(x), std::get<3>(x));
        }
#endif

        /**
         * @brief   Keep track of all '>' and the following EOL character to save the positions of the fasta record headers for each fasta record in the local block
         * @attention   This function should be called by all the MPI ranks in MPI communicator
         * @note        This function should be called by only 1 thread in a process.
         *              This function assumes that the searchRange are NOT OVERLAPPING between processes.  Else we could have messy duplicates.
         *
         * @tparam      Iterator      type of iterator for data to traverse.  raw pointer if no L2Buffering, or vector's iterator if buffering.
         * @param[in]   _data         start of iterator.
         * @param[in]   parentRange   the "full" range to which the inMemRange belongs.  used to determine if the target is a "first" block and last block.  has to start and end with valid delimiters (@)
         * @param[in]   inMemRange    the portion of the "full" range that's loaded in memory (e.g. in DataBlock).  does NOT have to start and end with @
         * @param[in]   searchRange   the range in which to search for a record.  the start and end position in the file.  does NOT have to start and end with @. should be inside inMemRange
         * @param[out]  sequences  vector of positions of the fasta sequence headers.
         *              Each pair in the vector represents the position of '>' and '\n' in the fasta record header
         */
        void init_for_iterator(const Iterator &_data, const RangeType &parentRange, const RangeType &inMemRange, const RangeType &searchRange)
          throw (bliss::io::IOException) {

          // now populate sequences

          if (!::std::is_same<typename std::iterator_traits<Iterator>::iterator_category, std::random_access_iterator_tag>::value)
            WARNING("WARNING: countSequenceStarts input iterator is not random access iterator.  lower.");

          using TT = typename std::iterator_traits<Iterator>::value_type;

          // make sure searchRange is inside parent Range.
          RangeType r = RangeType::intersect(searchRange, parentRange);

          // clear the starting position storage.
          sequences.clear();

          TT prev_char = '\n';  // default to EOL char.

          //============== OUTLINE
          // search for "\n>"
          // if internal, fine.
          // if at boundary,
          // if first char of proc i is ">", and rank is 0, can assume this is a starting point.  else, if prev char is "\n", then this is a start.

            //=====DONE===== GET LAST CHAR FROM PREVIOUS PROCESSOR, to check if first char is at start of line.


            //============== GET THE POSITION OF START OF EACH LINE.
            // local storage
            using LL = std::pair< typename RangeType::ValueType, uint8_t>;
            std::vector< LL > line_starts;

            // initialize the iterators
            Iterator it = _data;
            Iterator end = it;
            std::advance(end, r.size());
            auto i = searchRange.start;

            // ==== find all the beginning of line positions - exclude consecutive eols.

            // now determine if the first character is a line start
            if (prev_char == '\n') {
              // previous char is eol, so add the position here if it's not another eol, and encode it for header vs not header
              line_starts.emplace_back(i, ((*it == ';') || (*it == '>')) ? 1 : 0);
            }

            // now search for occurrences of '\n', since we've already handled the case where the first char is '>'
            Iterator it2 = it;
            ++it;
            ++i;
            while (it != end) {
              if (*it2 == '\n'){
                // previous char is eol, so add the position here if it's not another eol, and encode it for header vs not
                line_starts.emplace_back(i, ((*it == ';') || (*it == '>')) ? 1 : 0);
              }
              ++it;
              ++it2;
              ++i;
            }
            // add a very last line to mark end of file.
            line_starts.emplace_back(parentRange.end, 1);
            //=====DONE===== GET THE POSITION OF START OF EACH LINE.


            //============== KEEP JUST THE FIRST POSITION OF CONSECUTIVE HEADER OR SEQUENCE LINES  (== unique)
            auto new_end = line_starts.end();
            // filter out duplicates locally - i.e. consecutive lines with > or ; (true), and consecutive lines with standard alphabet (false).
            new_end = std::unique(line_starts.begin(), line_starts.end(), [](LL const & x, LL const & y) {
              return x.second == y.second;
            });
            line_starts.erase(new_end, line_starts.end());

            //=====DONE===== KEEP JUST THE FIRST POSITION OF CONSECUTIVE HEADER OR SEQUENCE LINES



            //============== CONVERT FROM LINE START POSITION TO SEQUENCE OBJECTS
              // insert first entry into sequences
              typename RangeType::ValueType pos, pos2, pos3;
              size_t k = 0;

              if (line_starts[k].second == 0 ) // not starting with a header.  move to next
                ++k;

              pos = line_starts[k].first;  // first start.
              ++k;

              // insert all other entries into sequences
              for ( ; (k + 1) < line_starts.size() ; k += 2) {

                pos2 = line_starts[k].first;
                pos3 = line_starts[k + 1].first;

                // insert if there is a next line or not.
                sequences.emplace_back(pos, pos2, pos3, k/2);

                pos = pos3;
              }

            //=====DONE===== CONVERT FROM LINE START POSITION TO RANGE FOR HEADERS


//            for (auto x : sequences)
//              DEBUGF("R 0 header range [%lu, %lu, %lu), id %lu\n", std::get<0>(x), std::get<1>(x), std::get<2>(x), std::get<3>(x));


        }
        /**
         * @brief searches the thread shared vector of sequence starts for where the offset belongs.  called as first increment.
         * @details  this method is almost the same as increment, but does a binary search whereas increment does not.
         *
         * @param iter
         * @param end
         * @param offset   offset from the file's start.
         * @param output
         * @return
         */
        size_t increment(Iterator &iter, const Iterator &end, size_t &offset, SequenceType &output) throw (bliss::io::IOException) {

          if (sequences.size() == 0) throw std::logic_error("calling FASTAParser increment without first initializingfor iterator.");

          // end.  return.
          if (iter == end) {
        	  WARNINGF("FASTA LOADER searching increment iter == end\n");
              output = SequenceType(SequenceIdType(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max()),
  	  	  	  	    0, 0, 0,
  	  	  	  	    end, end);
              return 0;
          }

          size_t input_dist = std::distance(iter, end);

          RangeType input(offset, offset + input_dist);

          // search for position using offset and the sequence's third element + std::upper_bound.
          // then we get first position where seq_end > offset, i.e. offset is inside the header or in the sequence.
          auto seqveciter = std::upper_bound(sequences.begin(), sequences.end(),
        		  offset, [](size_t const & off, SequenceOffsetsType const & seq){ return off < std::get<2>(seq); });
          size_t seq_offset = std::distance(sequences.begin(), seqveciter);

          if (seqveciter == sequences.end()) {
        	  WARNINGF("FASTA LOADER searching increment seqveciter == sequences.end()\n");
            // did not find a matching one.
            offset += std::distance(iter, end);
            iter = end;
            output = SequenceType(SequenceIdType(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max()),
	  	  	  	    0, 0, 0,
	  	  	  	    end, end);
            return seq_offset;
          }

          auto seq = *seqveciter;


          // the offset tuple from sequences, position 1 and 2 mark the begin and end of the dna sequence
//          if (std::get<1>(seq) > std::get<2>(seq)) {
//            DEBUGF("seq: offset %lu start %lu seq %lu end %lu, id %lu\n", offset, std::get<0>(seq), std::get<1>(seq), std::get<2>(seq), std::get<3>(seq));
//          }
          RangeType found(std::get<1>(seq), std::get<2>(seq));

          // intersect them to identify the valid region.
          RangeType valid = RangeType::intersect(input, found);

          if (valid.size() == 0) {
        	  DEBUGF("FASTA LOADER search increment valid.size() == 0.  seq: [%lu, %lu, %lu).  block range [%lu, %lu)\n",
        			  std::get<0>(seq), std::get<1>(seq), std::get<2>(seq),
					  input.start, input.end);
            // no valid.
            iter = end;
            offset += input_dist;
            output = SequenceType(SequenceIdType(std::get<0>(seq), std::get<3>(seq)),
	  	  	  	    std::get<2>(seq) - std::get<0>(seq),
                  std::get<1>(seq) - std::get<0>(seq),
                    valid.start - std::get<0>(seq),
                    end, end);
            return seq_offset;
          } // else there is valid

          // compute the output sequencetype = start to end of valid
          std::advance(iter, valid.start - offset);
          Iterator start = iter;
          std::advance(iter, valid.size());

          output = SequenceType(SequenceIdType(std::get<0>(seq), std::get<3>(seq)),
                                std::get<2>(seq) - std::get<0>(seq),
                                std::get<1>(seq) - std::get<0>(seq),
                                valid.start - std::get<0>(seq),
								start, iter);

//          printf("Sequence Id: seq_offset=%lu\n", seq_offset);

          // compute the new offset = end of valid.
          offset = valid.end;

          return ++seq_offset;
        }



        /**
         * @brief increments the iterator to the beginning of the next record, while saving the current record in memory and update the record id.
         * @details   since the sequences are pre computed, we just need to adjust iter and offset appropriately.
         * @note      ASSUMPTION IS THAT init_for_iterator HAS ALREADY BEEN CALLED.  ELSE THE SAME SEQUENCE MAY BE ACCESSED BY DIFFERENT THREADS.
         *
         * @param[in/out] iter          source iterator, pointing to data to traverse
         * @param[in]     end           end of the source iterator - not to go past.
         * @param[in/out] offset        starting position in units of source character types (e.g. char) from the beginning of the file.
         * @param[in/out] output        updated sequence type, values updated.
         * @throw IOException           if parse failed, throw exception
         */
        size_t increment(Iterator &iter, const Iterator &end, size_t &offset, size_t & seq_offset, SequenceType& output) throw (bliss::io::IOException)
        {

          if (sequences.size() == 0) throw std::logic_error("calling FASTAParser increment without first initializingfor iterator.");

          // end.  return.
          if (iter == end) {
        	  WARNINGF("FASTA LOADER iter == end\n");
              output = SequenceType(SequenceIdType(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max()),
  	  	  	  	    0, 0, 0,
					end, end);
              return seq_offset;
          }

          // use the seq_offset to get the offset group and convert to sequence type.

          // do a linear scan search to find the correct sequence
          // again, search until the seq_end > offset, i.e. offset is in header or seq.
          while ((seq_offset < sequences.size()) && (offset >= std::get<2>(sequences[seq_offset]))) {
        	  ++seq_offset;
          }

          size_t input_dist = std::distance(iter, end);

          // no more.  return.
          if (seq_offset >= sequences.size()) {
            if (offset > std::get<2>(sequences.back()))
              WARNINGF("FASTA LOADER did not find sequence for offset %lu.  last one ended at %lu\n", offset, std::get<2>(sequences.back()));
            offset += input_dist;
            iter = end;
            output = SequenceType(SequenceIdType(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max()),
	  	  	  	    0, 0, 0,
	  	  	  	    end, end);
            return seq_offset;
          }


          auto seq = sequences[seq_offset];

          // compute the valid range.  get input, get the selected, then intersect them.  this is in case iter-end is only part of a sequence.
          RangeType input(offset, offset + input_dist);

//          if (std::get<1>(seq) > std::get<2>(seq)) {
//            DEBUGF("seq: offset %lu start %lu seq %lu end %lu, id %lu\n", offset, std::get<0>(seq), std::get<1>(seq), std::get<2>(seq), std::get<3>(seq));
//          }
          RangeType found(std::get<1>(seq), std::get<2>(seq));

          RangeType valid = RangeType::intersect(input, found);

          // or potentially completely just header.
          if (valid.size() == 0) {
            // no valid.
        	  DEBUGF("FASTA LOADER valid.size() == 0.  seq: [%lu, %lu, %lu).  block range [%lu, %lu)\n",
        			  std::get<0>(seq), std::get<1>(seq), std::get<2>(seq),
					  input.start, input.end);

            offset += input_dist;
            iter = end;
            output = SequenceType(SequenceIdType(std::get<0>(seq), std::get<3>(seq)),
	  	  	  	    std::get<2>(seq) - std::get<0>(seq),
                  std::get<1>(seq) - std::get<0>(seq),
                    valid.start - std::get<0>(seq),
					end, end);
            return seq_offset;
          } // else there is valid

          // now convert to sequence type
          std::advance(iter, valid.start - offset);
          Iterator start = iter;
          std::advance(iter, valid.size());

          output = SequenceType(SequenceIdType(std::get<0>(seq), std::get<3>(seq)),
        		  	  	  	    std::get<2>(seq) - std::get<0>(seq),
                            std::get<1>(seq) - std::get<0>(seq),
                                valid.start - std::get<0>(seq),
								start, iter);

//          printf("Sequence Id: seq_offset=%lu\n", seq_offset);

          offset = valid.end;

          return ++seq_offset;
        }


    };



    /**
     * @class FASTQLoader
     * @brief FileLoader subclass specialized for the FASTQ file format, uses CRTP to enforce interface consistency.
     * @details   FASTQLoader understands the FASTQ file format, and enforces that the
     *            L1 and L2 partition boundaries occur at FASTQ sequence record boundaries.
     *
     *            FASTQ files allow storage of per-base quality score and is typically used
     *              to store a collection of short sequences.
     *            Long sequences (such as complete genome) can be found in FASTA formatted files
     *              these CAN be read using standard FileLoader, provided an appropriate overlap
     *              is specified (e.g. for kmer reading, overlap should be k-1).
     *
     *          for reads, expected lengths are maybe up to 1K in length, and files contain >> 1 seq
     *
     *          Majority of the interface is inherited from FileLoader.
     *          This class modifies the behaviors of GetNextL1Block and getNextL2Block,
     *            as described in FileLoader's "Advanced Usage".  Specifically, the range is modified
     *            for each block by calling "findStart" to align the start and end to FASTA sequence
     *            record boundaries.
     *
     * @note    Processes 1 file at a time.  For paired end reads, a separate subclass is appropriate.
     *          Sequences are identified by the following set of attributes:
     *
     *             file id.          (via filename to id map)
     *             read id in file   (e.g. using the offset at whcih read occurs in the file as id)
     *             position in read  (useful for kmer)
     *
     *          using the offset as id allows us to avoid d assignment via prefix sum.
     *
     * @tparam  T                 type of each element read from file
     * @tparam  Overlap           overlap between blocks at L1 or L2.  should be more than 0
     * @tparam  L1Buffering       bool indicating if L1 partition blocks should be buffered.  default to false
     * @tparam  L2Buffering       bool indicating if L2 partition blocks should be buffered.  default to true (to avoid contention between threads)
     * @tparam  L1PartitionerT    Type of the Level 1 Partitioner to generate the range of the file to load from disk
     * @tparam  L2PartitionerT    L2 partitioner, default to DemandDrivenPartitioner
     */
    template<typename T, size_t Overlap,
        bool L2Buffering = false,
        bool L1Buffering = true,
        typename L2PartitionerT = bliss::partition::DemandDrivenPartitioner<bliss::partition::range<size_t> > >
    using FASTALoader2 = FileLoader<T, Overlap, FASTAParser2, L2Buffering, L1Buffering, L2PartitionerT,
        bliss::partition::BlockPartitioner<bliss::partition::range<size_t> > >;







    /**
     * @class FASTALoader
     * @brief FileLoader subclass specialized for the FASTA file format
     */
    template<typename Iterator, typename KmerType>
    class FASTALoader
    {
      protected:

        /// constant representing the EOL character
        static const typename std::iterator_traits<Iterator>::value_type eol = '\n';

        /// alias for this type, for use inside this class.
        using type = FASTALoader<Iterator, KmerType>;

      public:
        //==== exposing types from base class

        /// Type of range object.
        typedef bliss::partition::range<size_t>     RangeType;

        //Type of offsets and corresponding location vector
        typedef std::size_t  offSetType;
        typedef std::vector <std::tuple<offSetType, offSetType, offSetType> >   vectorType;

        /// MPI communicator used by fasta loader.
        MPI_Comm comm;

        /**
         * @brief   Constructor of this class
         */
        FASTALoader(const MPI_Comm& _comm)
          : comm(_comm)
        {}

      protected:
//        /**
//         * @brief  search for first EOL character in a iterator.
//         * @param[in/out] iter    iterator to advance
//         * @param[in]     end     position to stop the traversal
//         * @param[in/out] offset  the global offset within the file.
//         * @return        True if EOL char is found, and false if not found.
//         */
//        inline bool findEOL(Iterator& iter, const size_t end, size_t &offset) const{
//          while (*iter != type::eol) {
//            if (offset == end) {
//              return false;
//            }
//            ++iter;
//            ++offset;
//          }
//          return true;
//        }
//
//#ifdef USE_MPI
//        /**
//         * @brief MPI point to point communication to message single element (Rank 0 -> 1, 1 -> 2, ....-> end).
//         * @param[in]     iter      iterator to advance
//         * @param[in]     sendInfo  element to send to next pid
//         * @param[in]     datatype  MPI Data type of the element to message
//         * @param[out]    recvInfo  element received
//         */
//        template<typename type>
//        void MPICommForward(type *sendInfo, type *recvInfo, MPI_Datatype datatype) const
//        {
//            int noProcs, myRank;
//            MPI_Comm_size(comm, &noProcs);
//            MPI_Comm_rank(comm, &myRank);
//            MPI_Request request_send, request_recv;
//
//            //2.1 To fetch and store the information if we have to resolve previous block's broken header problem
//            if (myRank < noProcs -1 )
//            {
//              const int rank_dest = myRank + 1;
//              const int msg_tag_send = 1000 + rank_dest;
//              MPI_Isend(const_cast<type *>(sendInfo), 1, datatype, rank_dest, msg_tag_send, comm, &request_send);
//            }
//            if (myRank > 0)
//            {
//              const int rank_src = myRank - 1;
//              const int msg_tag_recv = 1000 + myRank;
//              MPI_Irecv(recvInfo, 1, datatype, rank_src, msg_tag_recv, comm, &request_recv);
//            }
//            MPI_Status stat;
//            if (myRank < noProcs - 1) MPI_Wait(&request_send, &stat);
//            if (myRank > 0) MPI_Wait(&request_recv, &stat);
//        }
//#endif
//
//
//#ifdef USE_MPI
//        /**
//         * @brief MPI point to point communication to message single element (Rank 0 <- 1, 1 <- 2, ....<- end).
//         * @tparam      type      datatype of the element
//         * @param[in]   iter      iterator to advance
//         * @param[in]   sendInfo  element to send to next pid
//         * @param[in]   datatype  MPI Data type of the element to message
//         * @param[out]  recvInfo  element received
//         */
//        template<typename type>
//        void MPICommBackward(const type *sendInfo, type *recvInfo, MPI_Datatype datatype) const
//        {
//            int noProcs, myRank;
//            MPI_Comm_size(comm, &noProcs);
//            MPI_Comm_rank(comm, &myRank);
//            MPI_Request request_send, request_recv;
//
//            //2.1 To fetch and store the information if we have to resolve previous block's broken header problem
//            if (myRank > 0)
//            {
//              const int rank_dest = myRank - 1;
//              const int msg_tag_send = 1000 + rank_dest;
//              MPI_Isend(const_cast<type *>(sendInfo), 1, datatype, rank_dest, msg_tag_send, comm, &request_send);
//            }
//            if (myRank < noProcs - 1)
//            {
//              const int rank_src = myRank + 1;
//              const int msg_tag_recv = 1000 + myRank;
//              MPI_Irecv(recvInfo, 1, datatype, rank_src, msg_tag_recv, comm, &request_recv);
//            }
//            MPI_Status stat;
//            if (myRank > 0) MPI_Wait(&request_send, &stat);
//            if (myRank < noProcs - 1) MPI_Wait(&request_recv, &stat);
//        }
//#endif

      public:
//        /**
//         * @brief   Keep track of all '>' and the following EOL character to save the positions of the fasta record headers for each fasta record in the local block
//         * @attention   This function should be called by all the MPI ranks in MPI communicator
//         * @tparam      Iterator      type of iterator for data to traverse.  raw pointer if no L2Buffering, or vector's iterator if buffering.
//         * @param[in]   _data         start of iterator.
//         * @param[in]   parentRange   the "full" range to which the inMemRange belongs.  used to determine if the target is a "first" block and last block.  has to start and end with valid delimiters (@)
//         * @param[in]   inMemRange    the portion of the "full" range that's loaded in memory (e.g. in DataBlock).  does NOT have to start and end with @
//         * @param[in]   searchRange   the range in which to search for a record.  the start and end position in the file.  does NOT have to start and end with @. should be inside inMemRange
//         * @param[out]  localStartLocStore  vector of positions of the fasta sequence headers.
//         *              Each pair in the vector represents the position of '>' and '\n' in the fasta record header
//         */
//        void countSequenceStarts(const Iterator _data, const RangeType &parentRange, const RangeType &searchRange, vectorType &localStartLocStore) const
//          throw (bliss::io::IOException) {
//
//            RangeType t = searchRange;
//
//            //== set up the iterator for later
//            // initialize the counter
//            std::size_t i = t.start;
//
//#ifdef USE_MPI
//            int noProcs, myRank;
//            MPI_Comm_size(comm, &noProcs);
//            MPI_Comm_rank(comm, &myRank);
//
//
//#endif
//
//
//            // set iterator for the data
//            Iterator data(_data);
//            //== at this point, data points to start of the search range.
//
//            //indicates whether the local block needs MPI communication to know first sequence start
//            bool changeFirstStartLocation;
//
//            //indicates if we encounter '>' without the following '\n' in our search range
//            //False by default
//            bool lastBrokenHeader = false;
//
//            //=== look for all the starting fasta record starting points in the search space
//            //Consider the first position seperately.
//            //Push the beginning start position (might be changed after mpi communication)
//            if ((*data == '>') || (*data == ';')) {
//              std::size_t storeI = i;
//              changeFirstStartLocation = false;
//              bool EOLfound = findEOL(data, t.end, i);
//              localStartLocStore.push_back(std::make_pair(storeI,i));
//              //Case if first character is '>' and there is no EOL following it in the search space
//              if (!EOLfound)
//                lastBrokenHeader = true;
//
//            }
//            else {
//              changeFirstStartLocation = true;
//              //To avoid its contribution in the MPI parallel maximum computation, using 0 as temporary value
//              localStartLocStore.push_back(std::make_pair(0, 0));
//            }
//
//
//            //Above code will parse through the '>' character so that redundant check is not made below
//
//            while (i < t.end) {
//              if ((*data == '>') || (*data == ';'))
//              {
//                std::size_t storeI = i;
//                bool EOLfound = findEOL(data, t.end, i);
//                //Push the position i
//                localStartLocStore.push_back(std::make_pair(storeI,i));
//
//                if (!EOLfound)
//                  lastBrokenHeader = true;
//              }
//              ++data;
//              ++i;
//            }
//
//#ifdef USE_MPI
//            ///1. Each block who doesn't see the start '>' of the first fasta record needs to fetch it from previous blocks
//            //MPI Collective operation can help here
//            offSetType localMaximum = localStartLocStore.back().first;
//            offSetType newFirstStartLocation;
//
//            mxx::datatype<offSetType> dt;
//            MPI_Exscan(&localMaximum, &newFirstStartLocation, 1, dt.type(), MPI_MAX, comm);
//
//            //Update the first starting location (first element in the pair)
//            if (changeFirstStartLocation == true)
//              localStartLocStore.front().first = newFirstStartLocation;
//
//
//            ///2. Resolve the blocks which reported last broken header, meaning couldn't parse '\n' after the last '>' symbol
//            //2.1 To fetch and store the information if we have to resolve previous block's broken header problem
//            bool resolveLastBrokenHeader = false;
//            MPICommForward<bool>(&lastBrokenHeader, &resolveLastBrokenHeader, MPI_BYTE);
//
//            //2.2 Find the position of EOL if resolveLastBrokenHeader is TRUE
//            offSetType iCpy = 0;
//            if(resolveLastBrokenHeader)
//            {
//              Iterator dataCpy(_data);
//              iCpy = t.start;
//              findEOL(dataCpy, t.end, iCpy);
//            }
//
//            //2.3 Send the EOL index back to the neighbour
//            offSetType nextBlockEOLIndex;
//            MPICommBackward<offSetType>(&iCpy, &nextBlockEOLIndex, dt.type());
//
//            //2.4 Update the local vector
//            if (lastBrokenHeader == true)
//              localStartLocStore.back().second = nextBlockEOLIndex;
//
//            ///3. Each block should know the end '\n' of the header of the last fasta record in the previous block
//            localMaximum = localStartLocStore.back().second;
//
//            offSetType newFirstHeaderEndLocation;
//
//            MPI_Exscan(&localMaximum, &newFirstHeaderEndLocation, 1, dt.type(), MPI_MAX, comm);
//
//            //Update the first starting location (first element in the pair)
//            if (changeFirstStartLocation == true)
//              localStartLocStore.front().second = newFirstHeaderEndLocation;
//
//            //Due to overlap in the search range, its crucial to know if there is nearby header in the next block
//            offSetType extendedEnd = std::min(parentRange.end, t.end + KmerType::size -1);
//            while (i < extendedEnd) {
//              if ((*data == '>') || (*data == ';'))
//              {
//                std::size_t storeI = i;
//                findEOL(data, extendedEnd , i);
//                //Push the position i
//                localStartLocStore.push_back(std::make_pair(storeI,i));
//
//              }
//              ++data;
//              ++i;
//            }
//
//
//#endif
//          }



        /**
         * @brief   Keep track of all '>' and the following EOL character to save the positions of the fasta record headers for each fasta record in the local block
         * @attention   This function should be called by all the MPI ranks in MPI communicator
         * @tparam      Iterator      type of iterator for data to traverse.  raw pointer if no L2Buffering, or vector's iterator if buffering.
         * @param[in]   _data         start of iterator.
         * @param[in]   parentRange   the "full" range to which the inMemRange belongs.  used to determine if the target is a "first" block and last block.  has to start and end with valid delimiters (@)
         * @param[in]   inMemRange    the portion of the "full" range that's loaded in memory (e.g. in DataBlock).  does NOT have to start and end with @
         * @param[in]   searchRange   the range in which to search for a record.  the start and end position in the file.  does NOT have to start and end with @. should be inside inMemRange
         * @param[out]  localStartLocStore  vector of positions of the fasta sequence headers.
         *              Each pair in the vector represents the position of '>' and '\n' in the fasta record header
         */
        void countSequenceStarts(const Iterator _data, const RangeType &parentRange, const RangeType &searchRange, vectorType &localStartLocStore) const
          throw (bliss::io::IOException) {

          if (!::std::is_same<typename std::iterator_traits<Iterator>::iterator_category, std::random_access_iterator_tag>::value)
            WARNING("WARNING: countSequenceStarts input iterator is not random access iterator.  lower.");

          using TT = typename std::iterator_traits<Iterator>::value_type;

          // make sure searchRange is inside parent Range.
          RangeType r = RangeType::intersect(searchRange, parentRange);

          int p = 1, myRank = 0;
          int in_p = 1, in_rank = 0;

          // clear the starting position storage.
          localStartLocStore.clear();
          TT prev_char = '\n';  // default to EOL char.


          //============== OUTLINE
          // search for "\n>"
          // if internal, fine.
          // if at boundary,
          // if first char of proc i is ">", and rank is 0, can assume this is a starting point.  else, if prev char is "\n", then this is a start.

#ifdef USE_MPI
          MPI_Comm_size(comm, &p);
          MPI_Comm_rank(comm, &myRank);

          // ===== get the previous character from the next smaller, non-empty rank.
          // do by excluding the empty procs
          MPI_Comm in_comm;
          mxx2::split_communicator_by_function(comm, [&r](){ return (r.size() == 0) ? MPI_UNDEFINED : 1; }, in_comm);

          // if no data, then in_comm is null, and does not participate further.
          if (in_comm == MPI_COMM_NULL) return;


           //============== GET LAST CHAR FROM PREVIOUS PROCESSOR, to check if first char is at start of line.

            // for the nonempty ones, shift right by 1.
           MPI_Comm_size(in_comm, &in_p);
           MPI_Comm_rank(in_comm, &in_rank);


            // get the current last character.
            TT last_char = 0;  // default for empty.  note that the real last char from the last non-empty proc is shifted to the right place.
            {
              Iterator last = _data;
              if (r.size() > 0) {
                std::advance(last, r.size() - 1);
                last_char = *last;
              }
            }
            // right shift using the non-empty procs.
            prev_char = mxx::right_shift(last_char, in_comm);

            if (in_rank == 0) prev_char = '\n';
#endif
            //=====DONE===== GET LAST CHAR FROM PREVIOUS PROCESSOR, to check if first char is at start of line.


            //============== GET THE POSITION OF START OF EACH LINE.
            // local storage
            using LL = std::pair< typename RangeType::ValueType, uint8_t>;
            std::vector< LL > line_starts;

            // initialize the iterators
            Iterator it = _data;
            Iterator end = it;
            std::advance(end, r.size());
            auto i = searchRange.start;

            // ==== find all the beginning of line positions - exclude consecutive eols.

            // now determine if the first character is a line start
            if (prev_char == '\n') {
              // previous char is eol, so add the position here if it's not another eol, and encode it for header vs not header
              line_starts.emplace_back(i, ((*it == ';') || (*it == '>')) ? 1 : 0);
            }

            // now search for occurrences of '\n', since we've already handled the case where the first char is '>'
            Iterator it2 = it;
            ++it;
            ++i;
            while (it != end) {
              if (*it2 == '\n'){
                // previous char is eol, so add the position here if it's not another eol, and encode it for header vs not
                line_starts.emplace_back(i, ((*it == ';') || (*it == '>')) ? 1 : 0);
              }
              ++it;
              ++it2;
              ++i;
            }
#ifdef USE_MPI
            if (in_rank == (in_p - 1)) {
              // add a eof entry.
              line_starts.emplace_back(parentRange.end, 1);
            }
#endif
            //=====DONE===== GET THE POSITION OF START OF EACH LINE.


            //============== KEEP JUST THE FIRST POSITION OF CONSECUTIVE HEADER OR SEQUENCE LINES  (== unique)
            auto new_end = line_starts.end();
#ifdef USE_MPI
            // filter out duplicates locally - i.e. consecutive lines with > or ; (true), and consecutive lines with standard alphabet (false).
            new_end = mxx2::unique_contiguous(line_starts, in_comm, [](LL const & x, LL const & y) {
              return x.second == y.second;
            });

#else
            // filter out duplicates locally - i.e. consecutive lines with > or ; (true), and consecutive lines with standard alphabet (false).
            new_end = std::unique(line_starts.begin(), line_starts.end(), [](LL const & x, LL const & y) {
              return x.second == y.second;
            });
#endif
            line_starts.erase(new_end, line_starts.end());

            //=====DONE===== KEEP JUST THE FIRST POSITION OF CONSECUTIVE HEADER OR SEQUENCE LINES



            //============== CONVERT FROM LINE START POSITION TO RANGE FOR HEADERS
            int in_rank2 = in_rank;
            int in_p2 = in_p;

#ifdef USE_MPI
            // split comm again for non-empty line start vectors
            MPI_Comm in_comm2;
            mxx2::split_communicator_by_function(in_comm, [&line_starts](){ return line_starts.size() == 0 ? MPI_UNDEFINED : 1; }, in_comm2);

            // get first from next proc (
            if (in_comm2 != MPI_COMM_NULL) {
              MPI_Comm_rank(in_comm2, &in_rank2);
              MPI_Comm_size(in_comm2, &in_p2);

              // do a shift
              LL temp = line_starts.front();
              LL next = mxx::left_shift(temp, in_comm2);

              // insert into array.
              if (in_rank2 < (in_p2 - 1)) {
                line_starts.emplace_back(next);
              }

              MPI_Comm_free(&in_comm2);
            }

            // get the second from next proc (so we can construct elements with header start, seq start, and seq end.)
            // each proc should have at least 2, except for the last proc may have 1.
            // split comm again for non-empty line start vectors
            mxx2::split_communicator_by_function(in_comm, [&line_starts](){ return line_starts.size() <= 1 ? MPI_UNDEFINED : 1; }, in_comm2);

            // get first from next proc (

            if (in_comm2 != MPI_COMM_NULL) {
              MPI_Comm_rank(in_comm2, &in_rank2);
              MPI_Comm_size(in_comm2, &in_p2);

              // do a shift
              LL temp = line_starts[1];
              LL next = mxx::left_shift(temp, in_comm2);

              // insert into array.
              if (in_rank2 < (in_p2 - 1)) {
                line_starts.emplace_back(next);
              }

              MPI_Comm_free(&in_comm2);
            }

            // now only work with procs with at least 1 complete entry (header start, seq start, seq end
            mxx2::split_communicator_by_function(in_comm, [&line_starts](){ return line_starts.size() < ((line_starts.size() > 0 && line_starts.front().second == 0) ? 4 : 3) ? MPI_UNDEFINED : 1; }, in_comm2);


            if (in_comm2 != MPI_COMM_NULL) {
              MPI_Comm_rank(in_comm2, &in_rank2);
              MPI_Comm_size(in_comm2, &in_p2);
              MPI_Comm_free(&in_comm2);
#endif

              // insert first entry into localStartLocStore
              typename RangeType::ValueType pos, pos2, pos3;
              int k = 0;

              if (line_starts[k].second == 0 ) // not starting with a header.  move to next
                ++k;

              pos = line_starts[k].first;  // first start.
              ++k;

              // insert all other entries into localStartLocStore
              for ( ; (k + 1) < line_starts.size() ; k += 2) {

                pos2 = line_starts[k].first;
                pos3 = line_starts[k + 1].first;

                // insert if there is a next line or not.
                localStartLocStore.emplace_back(pos, pos2, pos3);

                pos = pos3;
              }

#ifdef USE_MPI
            }
#endif
            //=====DONE===== CONVERT FROM LINE START POSITION TO RANGE FOR HEADERS


//            for (auto x : localStartLocStore)
//              DEBUGF("R %d header range [%lu, %lu)\n", myRank, x.first, x.second);


#ifdef USE_MPI

            using SL = std::tuple<typename RangeType::ValueType, typename RangeType::ValueType, typename RangeType::ValueType>;
            //============== NOW SPREAD THE LAST ENTRY FORWARD.  empty range (r) does not participate

            //== first define segments.  empty is 0, and non-empty is 1 (start of segment)
            uint8_t sstart = localStartLocStore.size() == 0 ? 0 : 1;
            // convert to unique segments
            size_t segf = mxx2::segment::to_unique_segment_id_from_start<size_t>(sstart, 0, in_comm);

            // get the last entry.
            SL headerPos = SL();
            if (sstart) headerPos = localStartLocStore.back();

            //== do scan then shift forward by 1.  this is not the same as exscan when we're doing segmented scan.

            // now do inclusive scan to spread the values forward
            headerPos = mxx2::segmented_scan::scan(headerPos, segf, [](SL & x, SL & y){ return (std::get<0>(x) > std::get<0>(y) ? x : y); }, in_comm);

            // then shift by 1
            headerPos = mxx::right_shift(headerPos, in_comm);

            //== and merge the results with next.  If there is no localStartLocStore content, then need to insert.
            if (in_rank > 0) {
              if ((localStartLocStore.size() == 0) || (std::get<0>(localStartLocStore.front()) > (r.start + KmerType::size))) {
                localStartLocStore.insert(localStartLocStore.begin(), headerPos);
              }

            }
            //=====DONE===== NOW SPREAD THE LAST ENTRY FORWARD.  empty range (r) does not participate

            if (in_comm != MPI_COMM_NULL) MPI_Comm_free(&in_comm);
#endif
        }
    };
  }
}

#endif /* FASTA_PARTITIONER_HPP_ */
