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

#include "config.hpp"

#include "common/base_types.hpp"
#include "io/file_loader.hpp"
#include "boost/mpi/datatype.hpp"

// include MPI
#if defined(USE_MPI)
#include "mpi.h"
#endif

namespace bliss
{
  namespace io
  {
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
            MPI_Comm_size(MPI_COMM_WORLD, &noProcs); 
            MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
            MPI_Request request_send, request_recv; 

            //2.1 To fetch and store the information if we have to resolve previous block's broken header problem
            if (myRank < noProcs -1 )
            {
              const int rank_dest = myRank + 1;
              const int msg_tag_send = 1000 + rank_dest;
              MPI_Isend(const_cast<type *>(sendInfo), 1, datatype, rank_dest, msg_tag_send, MPI_COMM_WORLD, &request_send);
            }
            if (myRank > 0)
            {
              const int rank_src = myRank - 1;
              const int msg_tag_recv = 1000 + myRank;
              MPI_Irecv(recvInfo, 1, datatype, rank_src, msg_tag_recv, MPI_COMM_WORLD, &request_recv);
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
            MPI_Comm_size(MPI_COMM_WORLD, &noProcs); 
            MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
            MPI_Request request_send, request_recv; 

            //2.1 To fetch and store the information if we have to resolve previous block's broken header problem
            if (myRank > 0)
            {
              const int rank_dest = myRank - 1;
              const int msg_tag_send = 1000 + rank_dest;
              MPI_Isend(const_cast<type *>(sendInfo), 1, datatype, rank_dest, msg_tag_send, MPI_COMM_WORLD, &request_send);
            }
            if (myRank < noProcs - 1)
            {
              const int rank_src = myRank + 1;
              const int msg_tag_recv = 1000 + myRank;
              MPI_Irecv(recvInfo, 1, datatype, rank_src, msg_tag_recv, MPI_COMM_WORLD, &request_recv);
            }
            MPI_Status stat;   
            if (myRank > 0) MPI_Wait(&request_send, &stat);
            if (myRank < noProcs - 1) MPI_Wait(&request_recv, &stat);
        }
#endif

      public:
        /**
         * @brief   Keep track of all '>' and the following EOL character to save the positions of the fasta record headers for each fasta record in the local block
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
            if (*data == '>') {
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

            //Above code will parse through the '>' character so that renundant check is not made below

            while (i < t.end) {
              if(*data == '>') 
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

            MPI_Exscan(&localMaximum, &newFirstStartLocation, 1, boost::mpi::get_mpi_datatype(newFirstStartLocation), MPI_MAX, MPI_COMM_WORLD);

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
            MPICommBackward<offSetType>(&iCpy, &nextBlockEOLIndex, boost::mpi::get_mpi_datatype(nextBlockEOLIndex));

            //2.4 Update the local vector
            if (lastBrokenHeader == true)
              localStartLocStore.back().second = nextBlockEOLIndex;

            ///3. Each block should know the end '\n' of the header of the last fasta record in the previous block
            localMaximum = localStartLocStore.back().second;

            offSetType newFirstHeaderEndLocation;

            MPI_Exscan(&localMaximum, &newFirstHeaderEndLocation, 1, boost::mpi::get_mpi_datatype(newFirstHeaderEndLocation), MPI_MAX, MPI_COMM_WORLD);

            //Update the first starting location (first element in the pair)
            if (changeFirstStartLocation == true)
              localStartLocStore.front().second = newFirstHeaderEndLocation;

            //Due to overlap in the search range, its crucial to know if there is nearby header in the next block
            offSetType extendedEnd = std::min(parentRange.end, t.end + KmerType::getKmerSize() -1);
            while (i < extendedEnd) {
              if(*data == '>') 
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
    };
  }
}

#endif /* FASTA_PARTITIONER_HPP_ */
