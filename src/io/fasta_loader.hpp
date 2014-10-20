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

// include MPI
#include <mpi.h>

namespace bliss
{
  namespace io
  {


    /**
     * @class FASTALoader
     * @brief FileLoader subclass specialized for the FASTA file format
     *
     *
     * @tparam  T                 type of each element read from file
     * @tparam  L1Buffering       bool indicating if L1 partition blocks should be buffered.  default to false
     * @tparam  L2Buffering       bool indicating if L2 partition blocks should be buffered.  default to true (to avoid contention between threads)
     * @tparam  L1PartitionerT    Type of the Level 1 Partitioner to generate the range of the file to load from disk
     * @tparam  L2PartitionerT    L2 partitioner, default to DemandDrivenPartitioner
     */
    template<typename T,
      bool L2Buffering = true,
      bool L1Buffering = false,
      typename L2PartitionerT = bliss::partition::DemandDrivenPartitioner<bliss::partition::range<size_t> >,
      typename L1PartitionerT = bliss::partition::BlockPartitioner<bliss::partition::range<size_t> > >
        class FASTALoader : public FileLoader<T, L2Buffering, L1Buffering, L2PartitionerT, L1PartitionerT, void>
    {
      protected:
        /// base class type (FileLoader with the specific template parameters)
        typedef FileLoader<T, L2Buffering, L1Buffering, L2PartitionerT, L1PartitionerT,
                FASTQLoader<T, L2Buffering, L1Buffering, L2PartitionerT, L1PartitionerT> >    SuperType;

        friend class FileLoader<T, L2Buffering, L1Buffering, L2PartitionerT, L1PartitionerT,
               FASTQLoader<T, L2Buffering, L1Buffering, L2PartitionerT, L1PartitionerT> >;

        /// constant representing the EOL character
        static const typename std::iterator_traits<Iterator>::value_type eol = '\n';

      public:
        //==== exposing types from base class
        /// Type of iterator for traversing the input data.
        typedef typename SuperType::InputIteratorType                   InputIteratorType;

        /// Type of range object.
        typedef typename SuperType::RangeType                            RangeType;

      protected:
        /**
         * @brief  search for first EOL character in a iterator.
         * @param[in/out] iter    iterator to advance
         * @param[in]     end     position to stop the traversal
         * @param[in/out] offset  the global offset within the file.
         * @return        position of EOL where the EOL char is found, or 0 if not found.
         */
        template<typename Iterator>
          inline Iterator& findEOL(Iterator& iter, const size_t end, size_t &offset) {
            while (*iter != type::eol) {
            if (offset == end) { 
              return 0;
            }
            ++iter;
            ++offset;
          }
            return offset;
          }

        /**
         * @brief MPI point to point communication to message single element (Rank 0 -> 1, 1 -> 2, ....-> end).
         * @param[in]     iter      iterator to advance
         * @param[in]     sendInfo  element to send to next pid
         * @param[in]     datatype  MPI Data type of the element to message
         * @param[out]    recvInfo  element received 
         */
        template<typename type>
        void MPICommForward(type *sendInfo, type *recvInfo, MPI_Datatype datatype)
        {
            int noProcs, myRank;
            MPI_Comm_size(MPI_COMM_WORLD, &noProcs); 
            MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
            MPI_Request request_send, request_recv; 

            //2.1 To fetch and store the information if we have to resolve previous block's broken header problem
            bool resolveLastBrokenHeader = false;
            if (myRank < noProcs -1 )
            {
              const int rank_dest = myRank + 1;
              const int msg_tag_send = 1000 + rank_dest;
              MPI_ISend(sendInfo, 1, datatype, rank_dest, msg_tag_send, MPI_COMM_WORLD, &request_send);
            }
            if (myRank > 0)
            {
              const int rank_src = myRank - 1;
              const int msg_tag_recv = 1000 + rank_src;
              MPI_IRecv(recvInfo, 1, datatype, rank_src, msg_tag_recv, MPI_COMM_WORLD, &request_recv);
            }
            MPI_Status stat;   
            MPI_Wait(&request_recv, &stat);
            MPI_Wait(&request_send, &stat);
        }

        /**
         * @brief MPI point to point communication to message single element (Rank 0 <- 1, 1 <- 2, ....<- end).
         * @tparam      type      datatype of the element 
         * @param[in]   iter      iterator to advance
         * @param[in]   sendInfo  element to send to next pid
         * @param[in]   datatype  MPI Data type of the element to message
         * @param[out]  recvInfo  element received 
         */
        template<typename type>
        void MPICommBackward(type *sendInfo, type *recvInfo, MPI_Datatype datatype)
        {
            int noProcs, myRank;
            MPI_Comm_size(MPI_COMM_WORLD, &noProcs); 
            MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
            MPI_Request request_send, request_recv; 

            //2.1 To fetch and store the information if we have to resolve previous block's broken header problem
            bool resolveLastBrokenHeader = false;
            if (myRank > 0)
            {
              const int rank_dest = myRank - 1;
              const int msg_tag_send = 1000 + rank_dest;
              MPI_ISend(sendInfo, 1, datatype, rank_dest, msg_tag_send, MPI_COMM_WORLD, &request_send);
            }
            if (myRank < noProcs - 1)
            {
              const int rank_src = myRank + 1;
              const int msg_tag_recv = 1000 + rank_src;
              MPI_IRecv(recvInfo, 1, datatype, rank_src, msg_tag_recv, MPI_COMM_WORLD, &request_recv);
            }
            MPI_Status stat;   
            MPI_Wait(&request_recv, &stat);
            MPI_Wait(&request_send, &stat);
        }



        /**
         * @brief   Keep track of all '>' and the following EOL character to save the positions of the fasta record headers for each fasta record in the local block
         * @tparam  Iterator      type of iterator for data to traverse.  raw pointer if no L2Buffering, or vector's iterator if buffering.
         * @param   _data         start of iterator.
         * @param   parentRange   the "full" range to which the inMemRange belongs.  used to determine if the target is a "first" block and last block.  has to start and end with valid delimiters (@)
         * @param   inMemRange    the portion of the "full" range that's loaded in memory (e.g. in DataBlock).  does NOT have to start and end with @
         * @param   searchRange   the range in which to search for a record.  the start and end position in the file.  does NOT have to start and end with @. should be inside inMemRange
         * @return vector of positions of the fasta sequence headers.
         */

        template<typename Iterator>
          const std::vector <pair<std::size_t, std::size_t>> countSequenceStarts(const Iterator _data, const RangeType &parentRange, const RangeType &inMemRange, const RangeType &searchRange) const
          throw (bliss::io::IOException) {

            typedef typename std::iterator_traits<Iterator>::value_type  ValueType;

            //== range checking
            assert(parentRange.contains(inMemRange));
            RangeType t = RangeType::intersect(searchRange, inMemRange); // intersection to bound target to between parent's ends.

            //== set up the iterator for later
            // initialize the counter
            std::size_t i = t.start;

            // set iterator for the data
            Iterator data(_data);
            //== at this point, data points to start of the search range.

            //Local vector to store all the fasta record starting locations within the local search space
            //Each pair in the vector represents the position of '>' and '\n' in the fasta record header
            std::vector <pair<std::size_t, std::size_t>> localStartLocStore;

            //indicates whether the local block needs MPI communication to know first sequence start
            bool changeFirstStartLocation;
            
            //indicates if we encounter '>' without the following '\n' in our search range
            //False by default
            bool lastBrokenHeader = false;

            //=== look for all the starting fasta record starting points in the search space
            //Consider the first position seperately.
            //Push the beginning start position (might be changed after mpi communication)
            if (*data == '>') {
              storeI = i;
              changeFirstStartLocation = false;
              std::size_t endIndex = findEOL<Iterator>(data, t.end, i);
              localStartLocStore.push_back(make_pair(storeI,endIndex));
              //Case if first character is '>' and there is no EOL following it in the search space
              if (!endIndex)
                lastBrokenHeader = true;

            }
            else {
              changeFirstStartLocation = true;
              //To avoid its contribution in the MPI parallel maximum computation, using 0 as temporary value
              localStartLocStore.push_back(make_pair(0, 0));   
            }

            //Above code will parse through the '>' character so that renundant check is not made below

            while (i < t.end) {
              if(*data == '>') 
              {
                storeI = i;
                std::size_t endIndex = findEOL<Iterator>(data, t.end, i);
                //Push the position i 
                localStartLocStore.push_back(make_pair(storeI,endIndex));
                if (!endIndex)
                  lastBrokenHeader = true;
              }
              ++data;
              ++i;
            }
            

#ifdef USE_MPI
            ///1. Each block who doesn't see the start '>' of the first fasta record needs to fetch it from previous blocks
            //MPI Collective operation can help here
            int localMaximum = localStartLocStore.back().first;
            int newFirstStartLocation;

            MPI_Allreduce(&localMaximum, &newFirstStartLocation, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

            //Update the first starting location (first element in the pair)
            if (changeFirstStartLocation == true)
              myvector.front().first = newFirstStartLocation;


            ///2. Resolve the blocks which reported last broken header, meaning couldn't parse '\n' after the last '>' symbol
            //2.1 To fetch and store the information if we have to resolve previous block's broken header problem
            bool resolveLastBrokenHeader = false;
            MPICommForward<bool>(&lastBrokenHeader, &resolveLastBrokenHeader, MPI_BYTE);
            
            //2.2 Find the position of EOL if resolveLastBrokenHeader is TRUE
            const std::size_t EOLIndex;
            if(resolveLastBrokenHeader)
            {
              Iterator dataCpy(_data);
              std::size_t iCpy = t.start;
              EOLIndex = findEOL<Iterator>(dataCpy, t.end, iCpy);
            }

            //2.3 Send the EOL index back to the neighbour        
            const std::size_t nextBlockEOLIndex;
            MPICommBackward<std::size_t>(&EOLIndex, &nextBlockEOLIndex, MPI_INT);

            //2.4 Update the local vector
            if (lastBrokenHeader == true)
              myvector.back().second = EOLIndex;

            ///3. Each block should know the end '\n' of the header of the last fasta record in the previous block
            localMaximum = localStartLocStore.back().second;
            int newFirstHeaderEndLocation;

            MPI_Allreduce(&localStartLocStore, &newFirstHeaderEndLocation, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

            //Update the first starting location (first element in the pair)
            if (changeFirstStartLocation == true)
              myvector.front().second = newFirstHeaderEndLocation;
#endif
            //return the vector 
            return localStartLocStore;
          }
    }
  }
}

#endif /* FASTA_PARTITIONER_HPP_ */
