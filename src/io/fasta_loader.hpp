/*
 * Copyright 2015 Georgia Institute of Technology
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/**
 * @file    fasta_loader.hpp
 * @ingroup bliss::io
 * @author  cjain
 * @brief   contains a FileLoader subclass to perform distributed and concurrent FASTA file access.
 *
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

#include "io/mxx_support.hpp"
#include <mxx/datatypes.hpp>
#include <mxx/sort.hpp> // for mxx::unique

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
    class FASTAParser : public bliss::io::BaseFileParser<Iterator >
    {

        // let another BaseFileParser with some other iterator type be a friend so we can convert.
        template <typename Iterator2>
        friend class FASTAParser;

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
        size_t seq_offset;

      public:

        /// default constructor.
        FASTAParser() : seq_offset(::std::numeric_limits<size_t>::max()) {};

        /// default destructor
        ~FASTAParser() {};



        /// converting constructor.
        template <typename Iterator2>
        FASTAParser(FASTAParser<Iterator2> const & other) : sequences(other.sequences), seq_offset(other.seq_offset) { }
        /// converting assignment operator that can transform the base iterator type.
        template <typename Iterator2>
        FASTAParser<Iterator>& operator=(FASTAParser<Iterator2> const & other) {
          sequences.assign(other.sequences.begin(), other.sequences.end());
          seq_offset = other.seq_offset;
          return *this;
        }

        using bliss::io::BaseFileParser<Iterator>::find_overlap_end;
        using bliss::io::BaseFileParser<Iterator>::reset;
        using bliss::io::BaseFileParser<Iterator>::init_parser;
        using bliss::io::BaseFileParser<Iterator>::find_first_record;
        using bliss::io::BaseFileParser<Iterator>::get_next_record;
        using bliss::io::BaseFileParser<Iterator>::get_record_size;
        using bliss::io::BaseFileParser<Iterator>::should_parse;

// use base class definition
//        /**
//         * @brief  given a block, find the starting point that best aligns with the first sequence object (here, just just the actual start)
//         * @note   not overridden since here we just look for the beginning of the assigned block, same as the BaseFileParser logic.
//         */
//        virtual std::size_t find_first_record(const Iterator &_data, const RangeType &parentRange, const RangeType &inMemRange, const RangeType &searchRange)
//        {
//          init_for_iterator(_data, parentRange, inMemRange, searchRange);
//          return bliss::io::BaseFileParser<Iterator>::find_first_record(_data, parentRange, inMemRange, searchRange);
//        }

        virtual typename RangeType::ValueType find_overlap_end(const Iterator &_data,
        		const RangeType &parentRange, const RangeType &inMemRange,
				typename RangeType::ValueType end, size_t overlap ) {
          // if overlap is 0, then no need to compute.
          if (overlap == 0) return end;

          RangeType olr(end, inMemRange.end);
          Iterator curr(_data);
          ::std::advance(curr, end - inMemRange.start);

          // then advance Overlap characters, excluding eol characters.
          auto pos = olr.start;
          for (size_t i = 0; (pos < olr.end) && (i < overlap); ++curr, ++pos) {
            // TODO: probably should not skip eol characters.  this is searching the tail end.
            if ((*curr != bliss::io::BaseFileParser<typename ::std::iterator_traits<Iterator>::pointer>::cr) &&
                (*curr != bliss::io::BaseFileParser<typename ::std::iterator_traits<Iterator>::pointer>::eol)) ++i;
          }

          // new end.
          return pos;
        }


        virtual void reset() {
          seq_offset = ::std::numeric_limits<size_t>::max();
          sequences.clear();
        }


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
        virtual size_t init_parser(const Iterator &_data, const RangeType &parentRange, const RangeType &inMemRange, const RangeType &searchRange, const mxx::comm& comm)
        {

//          std::cout << "rank " << comm.rank()<< "   parent " << parentRange << std::endl;
//          std::cout << "rank " << comm.rank()<< "   inmem " << inMemRange << std::endl;
//          std::cout << "rank " << comm.rank()<< "   search " << searchRange << std::endl;

          //== range checking
          if(!parentRange.contains(inMemRange)) {
            ::std::cout << "parent: " << parentRange << " , in mem: " << inMemRange << ::std::endl;
            throw std::invalid_argument("ERROR: Parent Range does not contain inMemRange");
          }

          // now populate sequences

          if (!::std::is_same<typename std::iterator_traits<Iterator>::iterator_category, std::random_access_iterator_tag>::value)
            BL_WARNING("WARNING: countSequenceStarts input iterator is not random access iterator.  lower.");

          if (sequences.size() > 0) {
            BL_WARNING("WARNING: fasta_parser init called without reset first.  previous records are cleared.");
            sequences.clear();
          }

          using TT = typename std::iterator_traits<Iterator>::value_type;

          // make sure searchRange is inside parent Range.
          RangeType r = RangeType::intersect(searchRange, parentRange);

          // construct the search Range so they do not overlap.  This is necessary so that the sequences tuples are computed correctly.
          if (comm.rank() == (comm.size()-1))
            (void) mxx::left_shift(r.start, comm);
          else
            r.end = mxx::left_shift(r.start, comm);

          // clear the starting position storage.

          TT prev_char = '\n';  // default to EOL char.

          //============== OUTLINE
          // search for "\n>"
          // if internal, fine.
          // if at boundary,
          // if first char of proc i is ">", and rank is 0, can assume this is a starting point.  else, if prev char is "\n", then this is a start.


          //BL_WARNINGF("Rank %d search range [%lu, %lu)\n", comm.rank(), r.start, r.end);


          // ===== get the previous character from the next smaller, non-empty rank.
          // do by excluding the empty procs
          mxx::comm in_comm = comm.split(r.size() == 0);

          // if no data, then in_comm is null, and does not participate further.
          if (r.size() == 0) {
            //::std::cout << r;
            return r.start;
          }  // else we are using in_comm != MPI_COMM_NULL

          //============== GET LAST CHAR FROM PREVIOUS PROCESSOR, to check if first char is at start of line.
            // for the nonempty ones, shift right by 1.

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

            if (in_comm.rank() == 0) prev_char = '\n';
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
            if (in_comm.rank() == (in_comm.size() - 1)) {
              // add a eof entry.
              line_starts.emplace_back(parentRange.end, 1);
            }
            //=====DONE===== GET THE POSITION OF START OF EACH LINE.


            //============== KEEP JUST THE FIRST POSITION OF CONSECUTIVE HEADER OR SEQUENCE LINES  (== unique)
            auto new_end = line_starts.end();
            // filter out duplicates locally - i.e. consecutive lines with > or ; (true), and consecutive lines with standard alphabet (false).
            auto second_eq = [](LL const & x, LL const & y) { return x.second == y.second; };
            new_end = mxx::unique(line_starts.begin(), line_starts.end(), second_eq, in_comm);

            line_starts.erase(new_end, line_starts.end());

            //=====DONE===== KEEP JUST THE FIRST POSITION OF CONSECUTIVE HEADER OR SEQUENCE LINES

//            for (auto x : line_starts)
//              BL_DEBUGF("R %d line starts %lu, new record? %d\n", myRank, x.first, x.second);


            //============== CONVERT FROM LINE START POSITION TO SEQUENCE OBJECTS

            // split comm again for non-empty line start vectors
            in_comm.with_subset(line_starts.size() > 0, [&](const mxx::comm& subcomm) {
              // do a shift
              LL temp = line_starts.front();
              LL next = mxx::left_shift(temp, subcomm);

              // insert into array if not last proc
              if (subcomm.rank() < (subcomm.size() - 1)) {
                line_starts.emplace_back(next);
              }
            });

            // get the second from next proc (so we can construct elements with header start, seq start, and seq end.)
            // each proc should have at least 2, except for the last proc may have 1.
            // split comm again for non-empty line start vectors

            // get first from next proc (

            in_comm.with_subset(line_starts.size() > 1, [&](const mxx::comm& subcomm) {
              // do a shift
              LL temp = line_starts[1];
              LL next = mxx::left_shift(temp, subcomm);

              // insert into array.
              if (subcomm.rank() < (subcomm.size() - 1)) {
                line_starts.emplace_back(next);
              }
            });


            //mxx2::split_communicator_by_function(in_comm, [&line_starts](){ return line_starts.size() < ((line_starts.size() > 0 && line_starts.front().second == 0) ? 4 : 3) ? MPI_UNDEFINED : 1; }, in_comm2);


            //if (in_comm2 != MPI_COMM_NULL) {

            // now only work with procs with at least 1 complete entry (header start, seq start, seq end
            bool has_complete = line_starts.size() >= ((line_starts.size() > 0 && line_starts.front().second == 0) ? 4 : 3);
            in_comm.with_subset(has_complete, [&](const mxx::comm& subcomm) {

              // insert first entry into sequences
              typename RangeType::ValueType pos, pos2, pos3;
              size_t k = 0;

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
              entries = mxx::exscan(entries, std::plus<size_t>(), subcomm);
              if (subcomm.rank() > 0) {
                  for (size_t ii = 0; ii < sequences.size(); ++ii) {
                    std::get<3>(sequences[ii]) += entries;
                  }
              }
            });
            //=====DONE===== CONVERT FROM LINE START POSITION TO SEQUENCE OBJECTS


//            for (auto x : sequences)
//              BL_DEBUGF("R %d before spreading [%lu, %lu, %lu), id %lu\n", myRank, std::get<0>(x), std::get<1>(x), std::get<2>(x), std::get<3>(x));


            using SL = SequenceOffsetsType;
            //============== NOW SPREAD THE LAST ENTRY FORWARD (right).  empty range (r) does not participate

            // get the last entry.
            //SL headerPos = SL();
            //if (sequences.size() > 0) headerPos = sequences.back();
            std::pair<int, SL> lastentry = (sequences.size() > 0) ? std::make_pair(in_comm.rank(), sequences.back()) : std::make_pair(-1, SL());

            //== do seg scan (to fill in empty procs).
            // now do inclusive scan to spread the values forward
            auto maxpair = [](const std::pair<int, SL>& x, const std::pair<int, SL>& y) { return x.first >= y.first ? x : y; };   // when used in this context with this function, scan does same as segscan.
            lastentry = mxx::scan(lastentry, maxpair, in_comm);


            //== first define segments.  empty is 0, and non-empty is 1 (start of segment)
            //uint8_t sstart = sequences.size() == 0 ? 0 : 1;
            // convert to unique segments
            //size_t segf = mxx2::segment::to_unique_segment_id_from_start<size_t>(sstart, 0, in_comm);


            //headerPos = mxx2::segmented_scan::scan(headerPos, segf, [](SL & x, SL & y){ return (std::get<0>(x) > std::get<0>(y) ? x : y); }, in_comm);

            // then shift by 1
            SL headerPos = mxx::right_shift(lastentry.second, in_comm);

            // seg scan to fill empty then shift is NOT the same as seg exscan


            //== and merge the results with next.  If there is no sequences content, then need to insert.
            if (in_comm.rank() > 0) {
              // no special case for overlap.  downstream processing should be aware of the overlap.
              if ((sequences.size() == 0) || (std::get<0>(sequences.front()) > r.start)) {
                sequences.insert(sequences.begin(), headerPos);
              }

            }
            //=====DONE===== NOW SPREAD THE LAST ENTRY FORWARD.  empty range (r) does not participate
            size_t total_sequences = mxx::allreduce(sequences.size(), comm);
            if (total_sequences == 0) throw std::logic_error("ERROR: no sequences found distributed.  wrong input file type?");

            //if (in_comm != MPI_COMM_NULL) MPI_Comm_free(&in_comm);


//            for (auto x : sequences)
//              BL_DEBUGF("R %d header range [%lu, %lu, %lu), id %lu\n", myRank, std::get<0>(x), std::get<1>(x), std::get<2>(x), std::get<3>(x));

            // no op, other than to use intersect to make sure that the start is valid and within parent, and in memry ranges.
            return RangeType::intersect(searchRange, inMemRange).start;
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
        virtual size_t init_parser(const Iterator &_data, const RangeType &parentRange, const RangeType &inMemRange, const RangeType &searchRange)
        {

          //== range checking
          if(!parentRange.contains(inMemRange)) throw std::invalid_argument("ERROR: Parent Range does not contain inMemRange");

          // now populate sequences

          if (!::std::is_same<typename std::iterator_traits<Iterator>::iterator_category, std::random_access_iterator_tag>::value)
            BL_WARNING("WARNING: countSequenceStarts input iterator is not random access iterator.  lower.");

          if (sequences.size() > 0) {
            BL_WARNING("WARNING: fasta_parser init called without reset first.  previous records are cleared.");
            // clear the starting position storage.
            sequences.clear();
          }

          using TT = typename std::iterator_traits<Iterator>::value_type;

          // make sure searchRange is inside parent Range.
          RangeType r = RangeType::intersect(searchRange, parentRange);

          if (r.size() == 0) return searchRange.start;

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

              if (sequences.size() == 0) {
                BL_WARNING("WARNING: no sequences found serially. either incorrect file type, or a part of file that does not contain sequence start was read.");
              }

//            for (auto x : sequences)
//              BL_DEBUGF("R 0 header range [%lu, %lu, %lu), id %lu\n", std::get<0>(x), std::get<1>(x), std::get<2>(x), std::get<3>(x));

              // no op, other than to use intersect to make sure that the start is valid and within parent, and in memry ranges.
              return RangeType::intersect(searchRange, inMemRange).start;

        }



        /**
         * @brief increments the iterator to the beginning of the next record, while saving the current record in memory and update the record id.
         * @details   since the sequences are pre computed, we just need to adjust iter and offset appropriately.
         * @note      ASSUMPTION IS THAT init_for_iterator HAS ALREADY BEEN CALLED.  ELSE THE SAME SEQUENCE MAY BE ACCESSED BY DIFFERENT THREADS.
         *
         * @param[in/out] iter          source iterator, pointing to data to traverse
         * @param[in]     end           end of the source iterator - not to go past.
         * @param[in/out] offset        starting position in units of source character types (e.g. char) from the beginning of the file.
         * @return                      sequence object
         */
        SequenceType get_next_record(Iterator &iter, const Iterator &end, size_t &offset)
        {

          if (sequences.size() == 0) throw std::logic_error("calling FASTAParser increment without first initializing for iterator.");

          // end.  return.
          if (iter == end) {
        	  BL_WARNINGF("FASTA LOADER get next record iter == end\n");
            return SequenceType(SequenceIdType(offset),
  	  	  	  	    0, // length
  	  	  	  	    0, // offset in record
  	  	  	  	    iter, iter);
          }

          if (seq_offset == ::std::numeric_limits<size_t>::max()) {  // first search.  do binary
            // search for position using offset and the sequence's third element + std::upper_bound.
            // then we get first position where seq_end > offset, i.e. offset is inside the header or in the sequence.
            auto seqveciter = std::upper_bound(sequences.begin(), sequences.end(),
                offset, [](size_t const & off, SequenceOffsetsType const & seq){
                  return off < ::std::get<2>(seq);
                });
            seq_offset = std::distance(sequences.begin(), seqveciter);

          } else {  // subsequent search.  do linear.
            // do a linear scan search to find the correct sequence
            // again, search until the seq_end > offset, i.e. offset is in header or seq.
            while ((seq_offset < sequences.size()) &&
                (offset >= std::get<2>(sequences[seq_offset]))) {
              ++seq_offset;
            }
          }

          // check if we found one.
          size_t input_dist = std::distance(iter, end);

          // no more.  return.
          if (seq_offset >= sequences.size()) {
            BL_WARNINGF("FASTA LOADER did not find sequence for offset %lu.  last one ended at %lu\n", offset, std::get<2>(sequences.back()));
            offset += input_dist;
            iter = end;

            return SequenceType(SequenceIdType(offset),
                    0, // length
                    0, // offset in record
                    iter, iter);
          }

          // use the seq_offset to get the offset group and convert to sequence type.
          auto seq = sequences[seq_offset];

          // compute the valid range.  get input, get the selected, then intersect them.  this is in case iter-end is only part of a sequence.
          RangeType input(offset, offset + input_dist);

          // the offset tuple from sequences, position 1 and 2 mark the begin and end of the dna sequence
//          if (std::get<1>(seq) > std::get<2>(seq)) {
//            BL_DEBUGF("seq: offset %lu start %lu seq %lu end %lu, id %lu\n", offset, std::get<0>(seq), std::get<1>(seq), std::get<2>(seq), std::get<3>(seq));
//          }
          RangeType found(std::get<1>(seq), std::get<2>(seq));

          // intersect them to identify the valid region.
          RangeType valid = RangeType::intersect(input, found);

          // or potentially completely just header.
          if (valid.size() == 0) {
            // no valid.
        	  BL_DEBUGF("FASTA LOADER valid.size() == 0.  seq: [%lu, %lu, %lu).  block range [%lu, %lu)\n",
        			  std::get<0>(seq), std::get<1>(seq), std::get<2>(seq),
					  input.start, input.end);

            // no valid.
            offset += input_dist;
            iter = end;
            return SequenceType(SequenceIdType(std::get<0>(seq), std::get<3>(seq)),
	  	  	  	      std::get<2>(seq) - std::get<0>(seq), // record length
                    offset - std::get<0>(seq),  // offset in record
                    iter, iter);
          } // else there is valid range of sequence data

          // now convert to sequence type
          std::advance(iter, valid.start - offset);
          Iterator start = iter;
          std::advance(iter, valid.size());
          offset = valid.end;
	
					++seq_offset;

          return SequenceType(SequenceIdType(std::get<0>(seq), std::get<3>(seq)),
                              std::get<2>(seq) - std::get<0>(seq),
                                valid.start - std::get<0>(seq),
                                start, iter);

        }

        /**
         * @brief   get the average record size in the supplied range
         * @return  return the records size and the internal data size
         */
        virtual ::std::pair<size_t, size_t> get_record_size(const Iterator &_data, const RangeType &parentRange, const RangeType &inMemRange, const RangeType &searchRange, size_t const count) {
          if (sequences.size() == 0) throw std::logic_error("calling FASTAParser get_record_size without first initializing for iterator.");

          if (count == 0) throw std::invalid_argument("ERROR: called FASTAParser get_record_size with count == 0");

          // now initialize the search.  no need to call the first increment separately since we know that we have to start from seq id 0.
          size_t seq_data_len = 0;
          size_t record_size = 0;

          size_t max = ::std::min(sequences.size(), count);

          auto seq = sequences[0];

          size_t i = 0;
          for (; i < max; ++i) {
            seq = sequences[i];
            record_size += ::std::get<2>(seq) - ::std::get<0>(seq);
            seq_data_len += ::std::get<2>(seq) - ::std::get<1>(seq);
          }

          seq_data_len /= i;
          record_size /= i;

          if (seq_data_len == 0) throw std::logic_error("ERROR: estimated sequence data size is 0");
          if (record_size == 0) throw std::logic_error("ERROR: estimated record size is 0");

          return std::make_pair(record_size, seq_data_len);

        }
#ifdef USE_MPI
       virtual ::std::pair<size_t, size_t> get_record_size(const Iterator &_data, const RangeType &parentRange, const RangeType &inMemRange, const RangeType &validRange, mxx::comm const & comm, size_t const count) {
         if (sequences.size() == 0) throw std::logic_error("calling FASTAParser get_record_size without first initializing for iterator.");

         if (count == 0) throw std::invalid_argument("ERROR: called FASTAParser get_record_size with count == 0");

         // now initialize the search.  no need to call the first increment separately since we know that we have to start from seq id 0.
         size_t seq_data_len = 0;
         size_t record_size = 0;


         size_t max = ::std::min(sequences.size(), (count + comm.size() - 1) / comm.size());

         auto seq = sequences[0];

         size_t i = 0;
         for (; i < max; ++i) {
           seq = sequences[i];
           record_size += ::std::get<2>(seq) - ::std::get<0>(seq);
           seq_data_len += ::std::get<2>(seq) - ::std::get<1>(seq);
         }

         i = mxx::allreduce(i, comm);  // total count

         // aggregate.
         seq_data_len = mxx::allreduce(seq_data_len, comm) / i;
         record_size = mxx::allreduce(record_size, comm) / i;


         if (seq_data_len == 0) throw std::logic_error("ERROR: estimated sequence data size is 0");
         if (record_size == 0) throw std::logic_error("ERROR: estimated record size is 0");

         return std::make_pair(record_size, seq_data_len);
       }


       virtual bool should_parse(const RangeType & range, const mxx::comm & comm) {
         return mxx::any_of(range.size() > 0, comm);
       }

#endif
    };


  }
}

#endif /* FASTA_PARTITIONER_HPP_ */
