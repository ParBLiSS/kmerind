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
 * @file    kmer_file_helper.hpp
 * @ingroup io
 * @author  Tony Pan <tpan7@gatech.edu>
 * @brief convenience helper functions for reading (and conditional reading) of blocks and parse into sequences and kmers.
 * @details 3 primary file readers are used in convenience functions.
 *      MPI-IO
 *      POSIX
 *      MMAP
 *
 *    several convenience functions are provided for different ways of reading the files.
 *      1. parser all reads for all kmers.
 *      2. split reads by N and parse segments into kmers.
 *      3. discard reads with N and parse remainders into kmers
 *      4. conditionally parse, allow different operations to occur based on some sequence predicate?
 *      5. conditionally emit k-mers based on some predicate, e.g. count index content.
 *      6. non-kmer producing operations, such as first and last valid k-mer.
 *
 */
#ifndef KMER_FILE_HELPER_HPP_
#define KMER_FILE_HELPER_HPP_

#include "bliss-config.hpp"

#if defined(USE_MPI)
#include "mpi.h"
#endif

//#if defined(USE_OPENMP)
//#include "omp.h"
//#endif

#include <unistd.h>     // sysconf
#include <sys/stat.h>   // block size.
#include <tuple>        // tuple and utility functions
#include <utility>      // pair and utility functions.
#include <type_traits>
#include <cctype>       // tolower.

#include "io/file.hpp"
#include "io/fastq_loader.hpp"
#include "io/fasta_loader.hpp"
//#include "io/fasta_iterator.hpp"

#include "utils/logging.h"
#include "utils/file_utils.hpp"
#include "common/kmer.hpp"
#include "common/base_types.hpp"
#include "common/sequence.hpp"
#include "utils/kmer_utils.hpp"
#include "index/kmer_hash.hpp"
#include "common/kmer_transform.hpp"
#include "containers/fsc_container_utils.hpp"

#include "io/kmer_parser.hpp"
#include "io/mxx_support.hpp"

#include "io/sequence_iterator.hpp"

#include "utils/benchmark_utils.hpp"
#include "utils/file_utils.hpp"

#include <fstream> // debug only
#include <iostream>  // debug only
#include <sstream>  // debug only

namespace bliss
{
namespace io
{




/**
 * @tparam MapType    container type
 * @tparam KmerParser   functor to generate kmer (tuple) from input.  specified here so we specialize for different index.  note KmerParser needs to be supplied with a data type.
 */
struct KmerFileHelper {


  /**
   * @brief  generate kmers or kmer tuples for 1 block of raw data.
   * @note   requires that SeqParser be passed in and operates on the Block's Iterators.
   *          Mostly, this is because we need to broadcast the state of SeqParser to procs on the same node (L1 seq info) and recreate on child procs a new SeqParser (for L2)
   * @tparam SeqParser    parser type for extracting sequences.  supports FASTQ and FASTA.  template template parameter, param is iterator
   * @tparam KmerParser           parser type for generating Kmer.  supports kmer, kmer+pos, kmer+count, kmer+pos/qual.
   * @tparam BlockType    input partition type, supports in memory (vector) vs memmapped.
   * @param partition
   * @param result        output vector.  should be pre allocated.
   */
  template <typename KmerParser, template <typename> class SeqParser, template <typename,  template <typename> class> class SeqIterType, typename BlockType>
  static std::pair<size_t, size_t> read_block(BlockType const & partition,
      SeqParser<typename BlockType::iterator> const &seq_parser,
      std::vector<typename KmerParser::value_type>& result) {

    // from FileLoader type, get the block iter type and range type
    using CharIterType = typename BlockType::const_iterator;

    //using SeqIterType = ::bliss::io::SequencesIterator<CharIterType, SeqParser >;

    //== sequence parser type
    KmerParser kmer_parser;
    ::bliss::utils::file::NotEOL not_eol;

    //== process the chunk of data

    //==  and wrap the chunk inside an iterator that emits Reads.
    SeqIterType<CharIterType, SeqParser> seqs_start(seq_parser, partition.cbegin(), partition.in_mem_cend(), partition.getRange().start);
    SeqIterType<CharIterType, SeqParser> seqs_end(partition.in_mem_cend());

    ::fsc::back_emplace_iterator<std::vector<typename KmerParser::value_type> > emplace_iter(result);

    size_t before = result.size();
    size_t seqs = 0;

    //== loop over the reads
    for (; seqs_start != seqs_end; ++seqs_start)
    {
      auto seq = *seqs_start;
      if (seq.seq_size() == 0) continue;
      //      std::cout << "** seq: " << (*seqs_start).id.id << ", ";
      //      ostream_iterator<typename std::iterator_traits<typename SeqType::IteratorType>::value_type> osi(std::cout);
      //      std::copy((*seqs_start).seq_begin, (*seqs_start).seq_end, osi);
      //      std::cout << std::endl;

      size_t start_offset = seq.seq_global_offset();

      // if seq data starts outside of valid, then skip
      if (start_offset >= partition.valid_range_bytes.end) {
        continue;
      }

      // check if last.  if yes, and seqParser is a FASTAParser, then inspect and change if needed
      if (::std::is_same<SeqParser<CharIterType>, ::bliss::io::FASTAParser<CharIterType> >::value) {
        // if seq data ends in overlap region, then go at most k-1 characters from end of valid range.
        if ((start_offset + seq.seq_size()) >= partition.valid_range_bytes.end) {
          // scan for k-1 characters, from the valid range end.
          auto endd = seq.seq_begin + (partition.valid_range_bytes.end - start_offset);
          size_t steps = KmerParser::kmer_type::size - 1;
          size_t count = 0;

          // iterate and find the k-1 chars in overlap, starting from valid end.  should be less than current seq end.
          while ((endd != seq.seq_end) && (count < steps)) {
            if (not_eol(*endd)) {
              ++count;
            }

            ++endd;
          }


          seq.seq_end = endd;
        }
      }

      emplace_iter = kmer_parser(seq, emplace_iter);
      if ((seq.seq_offset == seq.seq_begin_offset) ||
          (start_offset > partition.valid_range_bytes.start)) ++seqs;

      //      std::cout << "Last: pos - kmer " << result.back() << std::endl;

    }

    return std::make_pair(seqs, result.size() - before);
  }


  /**
   * @brief initialize the sequence parser, estimate capacity and reserver, and then call read_block to parse the actual data.
   */
  template <typename KmerParser, template <typename> class SeqParser, template <typename,  template <typename> class> class SeqIterType, typename BlockType>
  static std::pair<size_t, size_t> parse_file_data(const BlockType & partition,
                         std::vector<typename KmerParser::value_type>& result) {
      std::pair<size_t, size_t> read = {0, 0};

     constexpr int kmer_size = KmerParser::kmer_type::size;

      BL_BENCH_INIT(file);
      {
        // not reusing the SeqParser in loader.  instead, reinitializing one.
        BL_BENCH_START(file);
        SeqParser<typename BlockType::const_iterator> seq_parser;
        seq_parser.init_parser(partition.in_mem_cbegin(), partition.parent_range_bytes, partition.in_mem_range_bytes, partition.getRange());
        BL_BENCH_END(file, "mark_seqs", partition.getRange().size());

        //== reserve
        BL_BENCH_START(file);
        // modifying the local index directly here causes a thread safety issue, since callback thread is already running.
        // index reserve internally sends a message to itself.

        // call after getting first L1Block to ensure that file is loaded.  (rank 0 reads and broadcast)
        size_t record_size = 0;
        size_t seq_len = 0;
        std::tie(record_size, seq_len) = seq_parser.get_record_size(partition.cbegin(), partition.parent_range_bytes, partition.getRange(), partition.getRange(), 10);
        size_t est_size = (record_size == 0) ? 0 : (partition.getRange().size() + record_size - 1) / record_size;  // number of records
        est_size *= (seq_len < kmer_size) ? 0 : (seq_len - kmer_size + 1) ;  // number of kmers in a record
        result.reserve(result.size() + est_size);
        BL_BENCH_END(file, "reserve", est_size);

        BL_BENCH_START(file);
        //=== copy into array
        if (partition.getRange().size() > 0) {
          read = read_block<KmerParser, SeqParser, SeqIterType>(partition, seq_parser, result);
        }
        BL_BENCH_END(file, "read", read.second);
        // std::cout << "Last: pos - kmer " << result.back() << std::endl;
      }

      BL_BENCH_REPORT_NAMED(file, "index:read_file_data");
      return read;

  }

  template <typename FileType>
  static ::bliss::io::file_data open_file(const std::string & filename, const size_t overlap) {
        // file extension determines SeqParserType
        std::string extension = ::bliss::utils::file::get_file_extension(filename);
        std::transform(extension.begin(), extension.end(), extension.begin(), ::tolower);
        if ((extension.compare("fastq") != 0) && (extension.compare("fasta") != 0)) {
          throw std::invalid_argument("input filename extension is not supported.");
        }

        FileType fobj(filename, overlap);
        return fobj.read_file();
  }

  /**
   * @brief read a file's content and generate kmers, place in a vector as return result.
   * @note  static so can be used without instantiating a internal map.
   * @tparam SeqParser    parser type for extracting sequences.  supports FASTQ and FASTA.   template template parameter, param is iterator
   * @tparam KmerParser   parser type for generating Kmer.  supports kmer, kmer+pos, kmer+count, kmer+pos/qual.
   */
  template <typename FileType, typename KmerParser, template <typename> class SeqParser, template <typename,  template <typename> class> class SeqIterType>
  static std::pair<size_t, size_t> read_file(const std::string & filename,
                         std::vector<typename KmerParser::value_type>& result) {

      std::pair<size_t, size_t> read = {0, 0};

      constexpr int kmer_size = KmerParser::kmer_type::size;

      BL_BENCH_INIT(file);
      {  // ensure that fileloader is closed at the end.

        BL_BENCH_START(file);
        ::bliss::io::file_data partition = open_file<FileType>(filename, kmer_size - 1);
        BL_BENCH_END(file, "open", partition.getRange().size());

        //std::cout << "rank " << _comm.rank() << " mpiio " << partition.getRange() << " in mem " << partition.in_mem_range_bytes << std::endl;


        // not reusing the SeqParser in loader.  instead, reinitializing one.
        BL_BENCH_START(file);
        read = parse_file_data<KmerParser, SeqParser, SeqIterType>(partition, result);
        BL_BENCH_END(file, "read", read.second);
        // std::cout << "Last: pos - kmer " << result.back() << std::endl;
      }


      BL_BENCH_REPORT_NAMED(file, "io:read_file");
      return read;
  }


  /**
   * @brief read a file's content and generate kmers, place in a vector as return result.
   * @note  static so can be used wihtout instantiating a internal map.
   * @tparam SeqParser    parser type for extracting sequences.  supports FASTQ and FASTA.   template template parameter, param is iterator
   * @tparam KmerParser   parser type for generating Kmer.  supports kmer, kmer+pos, kmer+count, kmer+pos/qual.
   */
  template <typename KmerParser, template <typename> class SeqParser, template <typename,  template <typename> class> class SeqIterType>
  static std::pair<size_t, size_t> read_file_mpiio(const std::string & filename,
                         std::vector<typename KmerParser::value_type>& result) {

      return read_file<::bliss::io::parallel::mpiio_file<SeqParser >,
          KmerParser, SeqParser, SeqIterType>(filename, result);
  }


  /**
   * @brief read a file's content and generate kmers, place in a vector as return result.
   * @note  static so can be used wihtout instantiating a internal map.
   * @tparam SeqParser    parser type for extracting sequences.  supports FASTQ and FASTA.   template template parameter, param is iterator
   * @tparam KmerParser   parser type for generating Kmer.  supports kmer, kmer+pos, kmer+count, kmer+pos/qual.
   */
  template <typename KmerParser, template <typename> class SeqParser, template <typename,  template <typename> class> class SeqIterType>
  static std::pair<size_t, size_t> read_file_mmap(const std::string & filename,
                        std::vector<typename KmerParser::value_type>& result) {

      // partitioned file with mmap or posix is only slightly faster than mpiio and may result in more jitter when congested.
      return read_file<::bliss::io::parallel::partitioned_file<::bliss::io::mmap_file, SeqParser >,
          KmerParser, SeqParser, SeqIterType>(filename, result);

  }


  /**
   * @brief read a file's content and generate kmers, place in a vector as return result.
   * @note  static so can be used wihtout instantiating a internal map.
   * @tparam SeqParser    parser type for extracting sequences.  supports FASTQ and FASTA.   template template parameter, param is iterator
   * @tparam KmerParser   parser type for generating Kmer.  supports kmer, kmer+pos, kmer+count, kmer+pos/qual.
   */
  template <typename KmerParser, template <typename> class SeqParser, template <typename,  template <typename> class> class SeqIterType>
  static std::pair<size_t, size_t> read_file_posix(const std::string & filename,
                         std::vector<typename KmerParser::value_type>& result) {



      // partitioned file with mmap or posix do not seem to be much faster than mpiio and may result in more jitter when congested.
      return read_file<::bliss::io::parallel::partitioned_file<::bliss::io::posix_file, SeqParser >,
          KmerParser, SeqParser, SeqIterType>(filename, result);

  }


#if defined(USE_MPI)

  /**
   * @brief initialize the sequence parser, estimate capacity and reserver, and then call read_block to parse the actual data.
   */
  template <typename KmerParser, template <typename> class SeqParser, template <typename,  template <typename> class> class SeqIterType, typename BlockType>
  static  ::std::pair<size_t, size_t> parse_file_data(const BlockType & partition,
                         std::vector<typename KmerParser::value_type>& result, const mxx::comm & _comm) {
      ::std::pair<size_t, size_t> read = {0,0};

     constexpr int kmer_size = KmerParser::kmer_type::size;

      BL_BENCH_INIT(file);
      {
        // not reusing the SeqParser in loader.  instead, reinitializing one.
        BL_BENCH_START(file);
        SeqParser<typename BlockType::const_iterator> seq_parser;
        seq_parser.init_parser(partition.in_mem_cbegin(), partition.parent_range_bytes, partition.in_mem_range_bytes, partition.getRange(), _comm);
        BL_BENCH_END(file, "mark_seqs", partition.getRange().size());

        //== reserve
        BL_BENCH_START(file);
        // modifying the local index directly here causes a thread safety issue, since callback thread is already running.
        // index reserve internally sends a message to itself.

        // call after getting first L1Block to ensure that file is loaded.  (rank 0 reads and broadcast)
        size_t record_size = 0;
        size_t seq_len = 0;
        std::tie(record_size, seq_len) = seq_parser.get_record_size(partition.cbegin(), partition.parent_range_bytes, partition.getRange(), partition.getRange(), _comm, 10);
        size_t est_size = (record_size == 0) ? 0 : (partition.getRange().size() + record_size - 1) / record_size;  // number of records
        est_size *= (seq_len < kmer_size) ? 0 : (seq_len - kmer_size + 1) ;  // number of kmers in a record
        result.reserve(result.size() + est_size);
        BL_BENCH_END(file, "reserve", est_size);

        BL_BENCH_START(file);
        //=== copy into array
        if (partition.getRange().size() > 0) {
          read = read_block<KmerParser, SeqParser, SeqIterType>(partition, seq_parser, result);
        }
        BL_BENCH_END(file, "read", read.first);
        // std::cout << "Last: pos - kmer " << result.back() << std::endl;
      }

      BL_BENCH_REPORT_MPI_NAMED(file, "index:read_file_data", _comm);
      return read;

  }

  template <typename FileType>
  static ::bliss::io::file_data open_file(const std::string & filename, const size_t overlap, const mxx::comm & _comm) {
        // file extension determines SeqParserType
        std::string extension = ::bliss::utils::file::get_file_extension(filename);
        std::transform(extension.begin(), extension.end(), extension.begin(), ::tolower);
        if ((extension.compare("fastq") != 0) && (extension.compare("fasta") != 0)) {
          throw std::invalid_argument("input filename extension is not supported.");
        }

        FileType fobj(filename, overlap, _comm);
        return fobj.read_file();
  }

  /**
   * @brief read a file's content and generate kmers, place in a vector as return result.
   * @note  static so can be used without instantiating a internal map.
   * @tparam SeqParser    parser type for extracting sequences.  supports FASTQ and FASTA.   template template parameter, param is iterator
   * @tparam KmerParser   parser type for generating Kmer.  supports kmer, kmer+pos, kmer+count, kmer+pos/qual.
   */
  template <typename FileType, typename KmerParser, template <typename> class SeqParser, template <typename,  template <typename> class> class SeqIterType>
  static  ::std::pair<size_t, size_t> read_file(const std::string & filename,
                         std::vector<typename KmerParser::value_type>& result,
                         const mxx::comm & _comm) {

      ::std::pair<size_t, size_t> read = {0, 0};

      constexpr int kmer_size = KmerParser::kmer_type::size;

      BL_BENCH_INIT(file);
      {  // ensure that fileloader is closed at the end.

        BL_BENCH_START(file);
        ::bliss::io::file_data partition = open_file<FileType>(filename, kmer_size - 1, _comm);
        BL_BENCH_END(file, "open", partition.getRange().size());

        //std::cout << "rank " << _comm.rank() << " mpiio " << partition.getRange() << " in mem " << partition.in_mem_range_bytes << std::endl;


        // not reusing the SeqParser in loader.  instead, reinitializing one.
        BL_BENCH_START(file);
        read = parse_file_data<KmerParser, SeqParser, SeqIterType>(partition, result, _comm);
        BL_BENCH_END(file, "read", read.first);
        // std::cout << "Last: pos - kmer " << result.back() << std::endl;
      }


      BL_BENCH_REPORT_MPI_NAMED(file, "io:read_file", _comm);
      return read;
  }


  /**
   * @brief read a file's content and generate kmers, place in a vector as return result.
   * @note  static so can be used wihtout instantiating a internal map.
   * @tparam SeqParser    parser type for extracting sequences.  supports FASTQ and FASTA.   template template parameter, param is iterator
   * @tparam KmerParser   parser type for generating Kmer.  supports kmer, kmer+pos, kmer+count, kmer+pos/qual.
   */
  template <typename KmerParser, template <typename> class SeqParser, template <typename,  template <typename> class> class SeqIterType>
  static  ::std::pair<size_t, size_t> read_file_mpiio(const std::string & filename,
                         std::vector<typename KmerParser::value_type>& result,
                         const mxx::comm & _comm) {

      return read_file<::bliss::io::parallel::mpiio_file<SeqParser >,
          KmerParser, SeqParser, SeqIterType>(filename, result, _comm);
  }


  /**
   * @brief read a file's content and generate kmers, place in a vector as return result.
   * @note  static so can be used wihtout instantiating a internal map.
   * @tparam SeqParser    parser type for extracting sequences.  supports FASTQ and FASTA.   template template parameter, param is iterator
   * @tparam KmerParser   parser type for generating Kmer.  supports kmer, kmer+pos, kmer+count, kmer+pos/qual.
   */
  template <typename KmerParser, template <typename> class SeqParser, template <typename,  template <typename> class> class SeqIterType>
  static  ::std::pair<size_t, size_t> read_file_mmap(const std::string & filename,
                        std::vector<typename KmerParser::value_type>& result,
                        const mxx::comm & _comm) {

      // partitioned file with mmap or posix is only slightly faster than mpiio and may result in more jitter when congested.
      return read_file<::bliss::io::parallel::partitioned_file<::bliss::io::mmap_file, SeqParser >,
          KmerParser, SeqParser, SeqIterType>(filename, result, _comm);

  }


  /**
   * @brief read a file's content and generate kmers, place in a vector as return result.
   * @note  static so can be used wihtout instantiating a internal map.
   * @tparam SeqParser    parser type for extracting sequences.  supports FASTQ and FASTA.   template template parameter, param is iterator
   * @tparam KmerParser   parser type for generating Kmer.  supports kmer, kmer+pos, kmer+count, kmer+pos/qual.
   */
  template <typename KmerParser, template <typename> class SeqParser, template <typename,  template <typename> class> class SeqIterType>
  static ::std::pair<size_t, size_t> read_file_posix(const std::string & filename,
                         std::vector<typename KmerParser::value_type>& result,
                         const mxx::comm & _comm) {



      // partitioned file with mmap or posix do not seem to be much faster than mpiio and may result in more jitter when congested.
      return read_file<::bliss::io::parallel::partitioned_file<::bliss::io::posix_file, SeqParser >,
          KmerParser, SeqParser, SeqIterType>(filename, result, _comm);

  }
#endif




};


} /* namespace io */
} /* namespace bliss */

#endif /* KMER_FILE_HELPER_HPP_ */
