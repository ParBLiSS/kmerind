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

#include "common/base_types.hpp"
#include "io/file_loader.hpp"
#include "mxx/datatypes.hpp"

#include "mxx/datatypes.hpp"

#include "io/mxx_support.hpp"


// include MPI
#if defined(USE_MPI)
#include "mpi.h"
#endif

namespace bliss
{
  namespace io
  {
    //dummy struct to indicate FASTA format
    struct FASTA {
        using SequenceId = size_t;
    };


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
        typedef std::vector <std::pair<offSetType, offSetType>>   vectorType;

        /// MPI communicator used by fasta loader.
        MPI_Comm comm;

        /**
         * @brief   Constructor of this class
         */
        FASTALoader(const MPI_Comm& _comm)
          : comm(_comm)
        {}

      protected:
        /**
         * @brief  search for first EOL character in a iterator.
         * @param[in/out] iter    iterator to advance
         * @param[in]     end     position to stop the traversal
         * @param[in/out] offset  the global offset within the file.
         * @return        True if EOL char is found, and false if not found.
         */
        inline bool findEOL(Iterator& iter, const size_t end, size_t &offset) const{
          while (*iter != type::eol) {
            if (offset == end) { 
              return false;
            }
            ++iter;
            ++offset;
          }
          return true;
        }

#ifdef USE_MPI
        /**
         * @brief MPI point to point communication to message single element (Rank 0 -> 1, 1 -> 2, ....-> end).
         * @param[in]     iter      iterator to advance
         * @param[in]     sendInfo  element to send to next pid
         * @param[in]     datatype  MPI Data type of the element to message
         * @param[out]    recvInfo  element received 
         */
        template<typename type>
        void MPICommForward(type *sendInfo, type *recvInfo, MPI_Datatype datatype) const
        {
            int noProcs, myRank;
            MPI_Comm_size(comm, &noProcs); 
            MPI_Comm_rank(comm, &myRank);
            MPI_Request request_send, request_recv; 

            //2.1 To fetch and store the information if we have to resolve previous block's broken header problem
            if (myRank < noProcs -1 )
            {
              const int rank_dest = myRank + 1;
              const int msg_tag_send = 1000 + rank_dest;
              MPI_Isend(const_cast<type *>(sendInfo), 1, datatype, rank_dest, msg_tag_send, comm, &request_send);
            }
            if (myRank > 0)
            {
              const int rank_src = myRank - 1;
              const int msg_tag_recv = 1000 + myRank;
              MPI_Irecv(recvInfo, 1, datatype, rank_src, msg_tag_recv, comm, &request_recv);
            }
            MPI_Status stat;   
            if (myRank < noProcs - 1) MPI_Wait(&request_send, &stat);
            if (myRank > 0) MPI_Wait(&request_recv, &stat);
        }
#endif


#ifdef USE_MPI
        /**
         * @brief MPI point to point communication to message single element (Rank 0 <- 1, 1 <- 2, ....<- end).
         * @tparam      type      datatype of the element 
         * @param[in]   iter      iterator to advance
         * @param[in]   sendInfo  element to send to next pid
         * @param[in]   datatype  MPI Data type of the element to message
         * @param[out]  recvInfo  element received 
         */
        template<typename type>
        void MPICommBackward(const type *sendInfo, type *recvInfo, MPI_Datatype datatype) const
        {
            int noProcs, myRank;
            MPI_Comm_size(comm, &noProcs); 
            MPI_Comm_rank(comm, &myRank);
            MPI_Request request_send, request_recv; 

            //2.1 To fetch and store the information if we have to resolve previous block's broken header problem
            if (myRank > 0)
            {
              const int rank_dest = myRank - 1;
              const int msg_tag_send = 1000 + rank_dest;
              MPI_Isend(const_cast<type *>(sendInfo), 1, datatype, rank_dest, msg_tag_send, comm, &request_send);
            }
            if (myRank < noProcs - 1)
            {
              const int rank_src = myRank + 1;
              const int msg_tag_recv = 1000 + myRank;
              MPI_Irecv(recvInfo, 1, datatype, rank_src, msg_tag_recv, comm, &request_recv);
            }
            MPI_Status stat;   
            if (myRank > 0) MPI_Wait(&request_send, &stat);
            if (myRank < noProcs - 1) MPI_Wait(&request_recv, &stat);
        }
#endif

      public:
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

            RangeType t = searchRange;

            //== set up the iterator for later
            // initialize the counter
            std::size_t i = t.start;

#ifdef USE_MPI
            int noProcs, myRank;
            MPI_Comm_size(comm, &noProcs);
            MPI_Comm_rank(comm, &myRank);


#endif


            // set iterator for the data
            Iterator data(_data);
            //== at this point, data points to start of the search range.

            //indicates whether the local block needs MPI communication to know first sequence start
            bool changeFirstStartLocation;
            
            //indicates if we encounter '>' without the following '\n' in our search range
            //False by default
            bool lastBrokenHeader = false;

            //=== look for all the starting fasta record starting points in the search space
            //Consider the first position seperately.
            //Push the beginning start position (might be changed after mpi communication)
            if ((*data == '>') || (*data == ';')) {
              std::size_t storeI = i;
              changeFirstStartLocation = false;
              bool EOLfound = findEOL(data, t.end, i);
              localStartLocStore.push_back(std::make_pair(storeI,i));
              //Case if first character is '>' and there is no EOL following it in the search space
              if (!EOLfound)
                lastBrokenHeader = true;

            }
            else {
              changeFirstStartLocation = true;
              //To avoid its contribution in the MPI parallel maximum computation, using 0 as temporary value
              localStartLocStore.push_back(std::make_pair(0, 0));   
            }


            //Above code will parse through the '>' character so that redundant check is not made below

            while (i < t.end) {
              if ((*data == '>') || (*data == ';'))
              {
                std::size_t storeI = i;
                bool EOLfound = findEOL(data, t.end, i);
                //Push the position i 
                localStartLocStore.push_back(std::make_pair(storeI,i));

                if (!EOLfound)
                  lastBrokenHeader = true;
              }
              ++data;
              ++i;
            }

#ifdef USE_MPI
            ///1. Each block who doesn't see the start '>' of the first fasta record needs to fetch it from previous blocks
            //MPI Collective operation can help here
            offSetType localMaximum = localStartLocStore.back().first;
            offSetType newFirstStartLocation;

            mxx::datatype<offSetType> dt;
            MPI_Exscan(&localMaximum, &newFirstStartLocation, 1, dt.type(), MPI_MAX, comm);

            //Update the first starting location (first element in the pair)
            if (changeFirstStartLocation == true)
              localStartLocStore.front().first = newFirstStartLocation;


            ///2. Resolve the blocks which reported last broken header, meaning couldn't parse '\n' after the last '>' symbol
            //2.1 To fetch and store the information if we have to resolve previous block's broken header problem
            bool resolveLastBrokenHeader = false;
            MPICommForward<bool>(&lastBrokenHeader, &resolveLastBrokenHeader, MPI_BYTE);
            
            //2.2 Find the position of EOL if resolveLastBrokenHeader is TRUE
            offSetType iCpy = 0;
            if(resolveLastBrokenHeader)
            {
              Iterator dataCpy(_data);
              iCpy = t.start;
              findEOL(dataCpy, t.end, iCpy);
            }

            //2.3 Send the EOL index back to the neighbour        
            offSetType nextBlockEOLIndex;
            MPICommBackward<offSetType>(&iCpy, &nextBlockEOLIndex, dt.type());

            //2.4 Update the local vector
            if (lastBrokenHeader == true)
              localStartLocStore.back().second = nextBlockEOLIndex;

            ///3. Each block should know the end '\n' of the header of the last fasta record in the previous block
            localMaximum = localStartLocStore.back().second;

            offSetType newFirstHeaderEndLocation;

            MPI_Exscan(&localMaximum, &newFirstHeaderEndLocation, 1, dt.type(), MPI_MAX, comm);

            //Update the first starting location (first element in the pair)
            if (changeFirstStartLocation == true)
              localStartLocStore.front().second = newFirstHeaderEndLocation;

            //Due to overlap in the search range, its crucial to know if there is nearby header in the next block
            offSetType extendedEnd = std::min(parentRange.end, t.end + KmerType::size -1);
            while (i < extendedEnd) {
              if ((*data == '>') || (*data == ';'))
              {
                std::size_t storeI = i;
                findEOL(data, extendedEnd , i);
                //Push the position i 
                localStartLocStore.push_back(std::make_pair(storeI,i));

              }
              ++data;
              ++i;
            }


#endif
          }



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
        void countSequenceStarts2(const Iterator _data, const RangeType &parentRange, const RangeType &searchRange, vectorType &localStartLocStore) const
          throw (bliss::io::IOException) {

          if (!::std::is_same<typename std::iterator_traits<Iterator>::iterator_category, std::random_access_iterator_tag>::value)
            WARNING("WARNING: countSequenceStarts input iterator is not random access iterator.  lower.");

          using TT = typename std::iterator_traits<Iterator>::value_type;

          // make sure searchRange is inside parent Range.
          RangeType r = RangeType::intersect(searchRange, parentRange);

          int p = 1, myRank = 0;

#ifdef USE_MPI
          MPI_Comm_size(comm, &p);
          MPI_Comm_rank(comm, &myRank);

          // ===== get the previous character from the next smaller, non-empty rank.
          // do by excluding the empty procs
          MPI_Comm in_comm;
          mxx2::split_communicator_by_function(comm, [&r](){ return (r.size() == 0) ? MPI_UNDEFINED : 1; }, in_comm);
#endif

          // clear the starting position storage.
          localStartLocStore.clear();

          // if no data, then in_comm is null, and does not participate further.
          if (in_comm == MPI_COMM_NULL) return;


          //============== OUTLINE
          // search for "\n>"
          // if internal, fine.
          // if at boundary,
          // if first char of proc i is ">", and rank is 0, can assume this is a starting point.  else, if prev char is "\n", then this is a start.


           //============== GET LAST CHAR FROM PREVIOUS PROCESSOR, to check if first char is at start of line.

           TT prev_char = '\n';  // default to EOL char.
           int seg = 1;
#ifdef USE_MPI
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
            //=====DONE===== GET THE POSITION OF START OF EACH LINE.


            //============== KEEP JUST THE FIRST POSITION OF CONSECUTIVE HEADER OR SEQUENCE LINES  (== unique)
            auto new_end = line_starts.end();
#if defined(USE_MPI)
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

            // get first from next proc.
            if (in_comm2 != MPI_COMM_NULL) {
              // do a shift
              LL temp = line_starts.front();
              LL next = mxx::left_shift(temp, in_comm2);

              MPI_Comm_rank(in_comm2, &in_rank2);
              MPI_Comm_size(in_comm2, &in_p2);
              MPI_Comm_free(&in_comm2);
#endif

              // insert first entry into localStartLocStore
              typename RangeType::ValueType pos, pos2;
              auto lit = line_starts.begin();

              if (lit->second == 0 ) // not starting with a header.  move to next
                ++lit;

              // insert all other entries into localStartLocStore
              for (;lit != line_starts.end();) {

                if (lit->second == 1)
                  pos = lit->first;
                else
                  ERRORF("ERROR: consecutive non-header lines.");

                ++lit;
                if (lit != line_starts.end()) {  // has next line
                  temp = *lit;
                  ++lit;
                } else {
                  if (in_rank2 < (in_p2 - 1))
                    temp = next;
                  else {
                    ERRORF("ERROR: ended with header start without end.");
                    break;
                  }

                }

                if (temp.second == 0)
                  pos2 = temp.first;
                else
                  ERRORF("ERROR: consecutive header lines.");


                // insert if there is a next line or not.
                localStartLocStore.emplace_back(pos, pos2);

              }

#ifdef USE_MPI
            }
#endif
            //=====DONE===== CONVERT FROM LINE START POSITION TO RANGE FOR HEADERS


            for (auto x : localStartLocStore)
              DEBUGF("R %d header range [%lu, %lu)\n", myRank, x.first, x.second);


#ifdef USE_MPI

            using SL = std::pair<typename RangeType::ValueType, typename RangeType::ValueType>;
            //============== NOW SPREAD THE LAST ENTRY FORWARD.  empty range (r) does not participate

            // first define segments.  empty is 0, and non-empty is 1 (start of segment)
            uint8_t sstart = localStartLocStore.size() == 0 ? 0 : 1;
            // convert to unique segments
            size_t segf = mxx2::segment::to_unique_segment_id_from_start<size_t>(sstart, 0, in_comm);

            // get the last entry.
            SL headerPos = SL();
            if (sstart) headerPos = localStartLocStore.back();

            // now do inclusive scan to spread the values forward
            headerPos = mxx2::segmented_scan::scan(headerPos, segf, [](SL & x, SL & y){ return (x.first > y.first ? x : y); }, in_comm);

            // then shift by 1
            headerPos = mxx::right_shift(headerPos, in_comm);

            // and merge the results with next

            if ((in_rank > 0) && (localStartLocStore.front().first > (r.start + KmerType::size))) {
              localStartLocStore.insert(localStartLocStore.begin(), headerPos);
            }
            //=====DONE===== NOW SPREAD THE LAST ENTRY FORWARD.  empty range (r) does not participate

#endif

            // finally, for compatibility, shift the second element back by 1
            std::for_each(localStartLocStore.begin(), localStartLocStore.end(), [](SL & x) {
              x.second --;
            });

#ifdef USE_MPI
            if (in_comm != MPI_COMM_NULL) MPI_Comm_free(&in_comm);
#endif
        }
    };
  }
}

#endif /* FASTA_PARTITIONER_HPP_ */
