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
 * @file    file_loader.hpp
 * @ingroup io
 * @author  Tony Pan <tpan7@gatech.edu>
 * @brief   contains a generic FileLoader class to perform distributed and concurrent file access.
 *
 */

#ifndef FILE_LOADER_HPP_
#define FILE_LOADER_HPP_

#include "bliss-config.hpp"

#if defined(USE_MPI)
#include "mpi.h"
#endif

#if defined(USE_OPENMP)
#include "omp.h"
#endif


#include <string>
#include <cstring>      // memcpy, strerror
#include <exception>    // ioexception
#include <sstream>      // stringstream
#include <memory>

#include <unistd.h>     // sysconf
#include <sys/stat.h>   // block size.
#include <sys/mman.h>   // mmap
#include <fcntl.h>      // for open

#include <sys/sysinfo.h>  // for meminfo
#include <type_traits>

#include "partition/range.hpp"
#include "partition/partitioner.hpp"
//#include "io/data_block.hpp"
#include "io/io_exception.hpp"
#include "utils/logging.h"
#include "common/sequence.hpp"
#include <mxx/comm.hpp> // for mxx::comm


namespace bliss
{
namespace io
{
  struct BaseFile {

  };

  // TODO: simplify  possibly by separating L1 and L2 partitioning.


  /**
   * @brief  Partitions input data into sequences, e.g. blocks or reads.
   * @details  Iterator template parameter refers to the iterator type that is currently being traversed.
   * 		 this applies to the sequence object, and the following methods:
   * 		 get_next_record, and find_first_record.
   * 		 
   * 		 init method can potentially use a different iterator type. this is useful when
   * 		   state information is extracted from perhaps a containing block (e.g. L1Block),
   * 		   then used with a sub block (e.g. L2Block).  
   * 		 	
   * @note   there are 2 (non-interchangeable) scenarios when we may need to change the iterator type.
   *    1. L1 at process level, L2 at thread level.  SeqParser state is calculated at L1.
   *    2. L1 at rank 0 of a sub communicator (processes), L2 is broadcast from rank0 of the rest of processes.  SeqParser state is calculated at L1.
   *
   *    to handle the first case: approaches are
   *    1a. to allow a parser to operate on different iterator types when invoking find_first_record, increment.
   *          additional template parameter for these functions.  SequenceType will need to change.
   *          API in SequenceIterator, BaseFileParser, FASTQParser, FASTAParser all need to change
   *    1b. convert the seq parser from L1 to L2, but copy the state of L1.
   *          can be shared between threads.  additional api for convert only
   *
   *    to handle the second case.
   *    2a. convert the seq parser from L1 to L2, communicate (and copy) the state of L1.
   *          additional api for convert and broadcast. minimal api change. requires L1 to be block (or cyclic with synchronization).
   *    2b. create seq parser for L2, and initialize using all processes.
   *          no api change, but works only when L1 and L2 are block partitioned (or cyclic with synchronization).
   *          potentially less communication, and Parser does not need to have a broadcast semantic.
   *
   *    1a is not general, potentially complicated, and a source of potential errors.
   *    1b and 2a share common apis.
   *    2b requires application to specifically call init at L2.  Also, already initialized L1 during data load.
   *
   *    choose 2b and 1b.  Sequence Iterator does not call init_for_iterator, mostly for flexibility and compatibility with case 1.
   *         
   * @tparam Iterator   The underlying iterator to be traversed to generate a Sequence
   *
   */
  template <typename Iterator = unsigned char* >
   class BaseFileParser {

      // let another BaseFileParser with some other iterator type be a friend so we can convert.
      template <typename Iterator2>
      friend class BaseFileParser;

     public:
//      using SequenceIdType = SeqIdType;
      using SequenceType = typename ::bliss::common::Sequence<Iterator, ::bliss::common::SequenceId>;
      using SequenceIdType = typename SequenceType::IdType;

      /// static constant for end of line.  note this is same for unicode as well
      static constexpr unsigned char eol = '\n';
      /// static constant for carriage return.  note this is same for unicode as well
      static constexpr unsigned char cr = '\r';

     protected:
       /**
        * @typedef RangeType
        * @brief   range object types
        */
      using RangeType = bliss::partition::range<size_t>;

       /**
        * @brief  search for first non-EOL character in a iterator, returns the stopping position as an iterator.  also records offset.
        * @details       iter can point to the previous EOL, or a nonEOL character.
        * @param[in/out] iter    iterator to advance
        * @param[in]     end     position to stop the traversal
        * @param[in/out] offset  the global offset of the sequence record within the file.  used as record id.
        * @return        iterator at the new position, where the Non EOL char is found, or end.
        */
       template <typename IT = Iterator,
           typename = typename ::std::enable_if<(::std::is_same<typename ::std::iterator_traits<IT>::value_type, char>::value ||
                                                ::std::is_same<typename ::std::iterator_traits<IT>::value_type, unsigned char>::value)
                                               >::type >
       inline IT findNonEOL(IT& iter, const IT& end, size_t &offset) const {
         while ((iter != end) && ((*iter == eol) || (*iter == cr))) {
           ++iter;
           ++offset;
         }
         return iter;
       }

       /**
        * @brief  search for first EOL character in a iterator, returns the stopping position as an iterator.  also records offset.
        * @details       iter can point to a nonEOL char, or an EOL char.
        * @param[in/out] iter    iterator to advance
        * @param[in]     end     position to stop the traversal
        * @param[in/out] offset  the global offset of the sequence record within the file.  used as record id.
        * @return        iterator at the new position, where the EOL char is found, or end.
        */
       template <typename IT = Iterator,
           typename = typename ::std::enable_if<(::std::is_same<typename ::std::iterator_traits<IT>::value_type, char>::value ||
                                                ::std::is_same<typename ::std::iterator_traits<IT>::value_type, unsigned char>::value)
                                               >::type >
       inline IT findEOL(IT& iter, const IT& end, size_t &offset) const {
         while ((iter != end) && ((*iter != eol) && (*iter != cr) ) ) {
           ++iter;
           ++offset;
         }
         return iter;
       }


       /**
        * @brief constructs an IOException object with the relevant debug data.
        * @param errType     string indicating the source of the error.  user choice
        * @param start       iterator pointing to beginning of the data in question
        * @param end         iterator pointing to end of the data in question
        * @param startOffset offset for the beginning of the data in question
        * @param endOffset   offset for the end of the data in question
        */
       void handleError(const std::string& errType, const Iterator &start, const Iterator &end, const size_t& startOffset, const size_t& endOffset) throw (bliss::io::IOException) {
         std::stringstream ss;
//#ifdef USE_MPI
//         int rank = 0;
//         MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//         ss << "rank: " << rank << " ";
//#endif

         ss << "ERROR: " << errType << " in " << startOffset << " to " << endOffset << std::endl;
         ss << "  offending string is \"";
         std::ostream_iterator<typename std::iterator_traits<Iterator>::value_type> oit(ss);
         std::copy(start, end, oit);
         ss << "\".";
         throw ::std::logic_error(ss.str());
       }

       /**
        * @brief print out an Warning, for example when there is malformed partial data.
        * @param errType     string indicating the source of the error.  user choice
        * @param start       iterator pointing to beginning of the data in question
        * @param end         iterator pointing to end of the data in question
        * @param startOffset offset for the beginning of the data in question
        * @param endOffset   offset for the end of the data in question
        */
       void handleWarning(const std::string& errType, const Iterator &start, const Iterator &end, const size_t& startOffset, const size_t& endOffset) {
         BL_WARNING("WARNING: " << errType << " in " << startOffset << " to " << endOffset);
       }


     public:

       /// default constructor.
       BaseFileParser() {};

       /// default destructor
       virtual ~BaseFileParser() {};


       /// converting constructor.
       template <typename Iterator2>
       BaseFileParser(BaseFileParser<Iterator2> const & other) {}
       /// converting assignment operator that can transform the base iterator type.
       template <typename Iterator2>
       BaseFileParser<Iterator>& operator=(BaseFileParser<Iterator2> const & other) {
         return *this;
       }


       /**
        * @brief given a block/range, find the starting point of the first sequence object (here, just the actual start)
        * @note  virtual function, for convenience.
        *         used by getNextL1BlockRange to find exact range for an L1 Partition.
        *
        *         Also can be used to do other initializations.
        *
        * @param _data
        * @param parentRange
        * @param inMemRange
        * @param searchRange
        * @return
        */
       virtual std::size_t find_first_record(const Iterator &_data, const RangeType &parentRange, const RangeType &inMemRange, const RangeType &searchRange)
       {
         //== range checking
         if(!parentRange.contains(inMemRange)) {
           ::std::cout << "parent: " << parentRange << " , in mem: " << inMemRange << ::std::endl;
           throw std::invalid_argument("ERROR: Parent Range does not contain inMemRange");
         }

         // no op, other than to use intersect to make sure that the start is valid and within parent, and in memry ranges.
         return RangeType::intersect(searchRange, inMemRange).start;
       }
#ifdef USE_MPI
       virtual std::size_t find_first_record(const Iterator &_data, const RangeType &parentRange, const RangeType &inMemRange, const RangeType &searchRange, const mxx::comm& comm) {
         return find_first_record(_data, parentRange, inMemRange, searchRange);
       }
#endif

       virtual typename RangeType::ValueType find_overlap_end(const Iterator &_data,
    		   const RangeType &parentRange, const RangeType &inMemRange, typename RangeType::ValueType end, size_t overlap ) {
         return std::max(inMemRange.start, std::min(end + overlap, inMemRange.end));
       }


       /// reset any parser internal state so the parser can be reused on a different data range.  overridden in subclass.
       virtual void reset() {}

       virtual bool should_parse(const RangeType & range) {
         return range.size() > 0;
       }
#ifdef USE_MPI
       virtual bool should_parse(const RangeType & range, const mxx::comm & comm) {
         return range.size() > 0;
       }
#endif

       /// initializes parser for a particular data range.  overridden in subclass.  returns start of first record.  _data points to first element in mem.
       virtual std::size_t init_parser(const Iterator &_data, const RangeType &parentRange, const RangeType &inMemRange, const RangeType &searchRange) {
         return this->find_first_record(_data, parentRange, inMemRange, searchRange);
       }
#ifdef USE_MPI
       /// initializes the parser.  only useful for FASTA parser for now.  Assumes searchRange do NOT overlap between processes.   _data points to first element in mem.
       virtual std::size_t init_parser(const Iterator &_data, const RangeType &parentRange, const RangeType &inMemRange, const RangeType &searchRange, const mxx::comm& comm)
       {
         return this->find_first_record(_data, parentRange, inMemRange, searchRange);
       };
#endif



       /**
        * @brief increment to get next sequence object
        * @note   not virtual because SequenceType is different for each subclass, so the signatures are different.  We rely on compiler to enforce the correct one.
        *         this is in general not an issue as we do not need polymorphic Sequence Parsers - they are specific to file types.
        *         this is probably better for performance anyways.
        *
        *         Used by FileLoader's getRecordSize to compute an average record size.
        *
        * @param iter           start of a sequence object.  this is the beginning of a record, not just DNA sequence
        * @param end            end of a sequence object.  this is the end of a record, not just DNA sequence
        * @param offset         offset in the file for the start of the record.
        * @param seq_offset     index in the local (shared) vector of sequence breaks, if there is one (e.g. FASTA).  used by FASTA to lookup the nearest complete sequence record info.
        * @param output         updated sequence object.
        * @return               next seq id offset.
        */
       SequenceType get_next_record(Iterator & iter, const Iterator & end, size_t & offset)
       {
         Iterator orig_iter = iter;
         size_t orig_offset = offset;

         size_t dist = std::distance(iter, end);

         offset += dist;
         iter = end;

         return SequenceType(SequenceIdType(orig_offset), dist, 0, orig_iter, end);
       }


       /**
        * @brief   get the average record size in the supplied range
        * @return  return the records size and the internal data size
        */
       virtual ::std::pair<size_t, size_t> get_record_size(const Iterator &_data, const RangeType &parentRange, const RangeType &inMemRange, const RangeType &validRange, size_t const count) {
         size_t dist = validRange.size();
         return ::std::pair<size_t, size_t>(dist,dist);
       }
#ifdef USE_MPI
       virtual ::std::pair<size_t, size_t> get_record_size(const Iterator &_data, const RangeType &parentRange, const RangeType &inMemRange, const RangeType &validRange, mxx::comm const & comm, size_t const count) {
         size_t dist = validRange.size();
         return ::std::pair<size_t, size_t>(dist,dist);
       }
#endif

   };
   /// template class' static variable definition (declared and initialized in class)
  template <typename Iterator>
   constexpr unsigned char bliss::io::BaseFileParser<Iterator >::eol;
  template <typename Iterator>
   constexpr unsigned char bliss::io::BaseFileParser<Iterator>::cr;


  } /* namespace functional */
} /* namespace bliss */
#endif /* FILE_LOADER_HPP_ */
